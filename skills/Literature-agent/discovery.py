#!/usr/bin/env python3
"""
Literature Discovery & Citation Network Module
================================================
基于 Semantic Scholar API 的文献发现与引用网络分析模块。
支持从种子文献向外扩展（引用/被引/相关推荐/关键词搜索），构建真实引用网络，
并进行 PageRank、桥梁文献、新兴聚类、盲区检测等分析。

API 限制（免费版）：
- 无需 API Key 即可使用，但有 IP 速率限制（约 100 requests/5min）
- 建议设置 rate_limit_delay >= 1.0 秒
- 如需要更高配额，可在 https://www.semanticscholar.org/product/api 申请 Key
"""

import os
import sys
import time
import json
import hashlib
import argparse
from typing import List, Dict, Optional, Set, Tuple, Any
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from pathlib import Path

import requests
import networkx as nx
import numpy as np

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    plt = None
    sns = None


class SemanticScholarRateLimitCircuitOpen(RuntimeError):
    """Semantic Scholar returned 429 and this run should stop SS requests."""


# ==================== 配置 ====================
DISCOVERY_DEFAULTS = {
    "rate_limit_delay": 1.0,
    "relevance_threshold": 0.25,
    "max_papers": 500,
    "max_depth": 1,
    "breadth_limit": 20,
}

# 光酶/光氧化还原/自由基领域关键词库（用于相关性评估）
DOMAIN_KEYWORDS = [
    # 核心概念
    'photoenzym', 'photobiocatal', 'enzyme photocatalysis',
    'photoredox', 'photocatalysis', 'visible light', 'photoinduced',
    'radical', 'radicals', 'minisci', 'giese',
    'energy transfer', 'triplet', 'hat', 'hydrogen atom transfer',
    'eda complex', 'electron donor acceptor', 'lmct',
    # 反应类型
    'decarboxylation', 'decarboxylative', 'decarboxyl',
    'arylation', 'alkylation', 'cycloaddition', 'halogenation',
    'hydrofunctionalization', 'cross-coupling', 'c-h functionalization',
    # 酶
    'p450', 'cytochrome', 'ered', 'ene-reductase',
    'old yellow enzyme', 'oye', 'fap', 'fatty acid photodecarboxylase',
    'photodecarboxylase', 'flavin', 'fad', 'fmn',
    'lipase', 'esterase', 'hydrolase', 'calb', 'candida antarctica',
    'biocatalysis', 'enzymatic', 'enzyme', 'cofactor', 'active site',
    # 选择性
    'stereoselective', 'enantioselective', 'regioselective',
    'site-selective', 'chemoselective', 'diastereoselective',
    # 机制研究
    'epr', 'esr', 'transient absorption', 'stern-volmer',
    'quantum yield', 'radical clock', 'isotope labeling',
    'dft', 'density functional', 'md simulation', 'molecular dynamics',
    'crystal structure', 'control experiment',
    # 催化剂/光敏剂
    'ru(bpy)', 'ir(ppy)', 'organic photocatalyst', 'quinolinone',
    # 底物
    'carboxylic acid', 'alkene', 'arene', 'heteroarene', 'pyridine',
    'sulfonamide', 'amine', 'carbonyl', 'imine', 'alkyl halide',
    'drug-like', 'natural product', 'amino acid',
]


@dataclass
class DiscoveredPaper:
    """发现文献的轻量数据结构，可与现有 PaperNode 互转"""
    doi: str
    title: str
    authors: List[str] = field(default_factory=list)
    year: int = 0
    abstract: Optional[str] = None
    citation_count: int = 0
    reference_count: int = 0
    fields_of_study: List[str] = field(default_factory=list)
    publication_types: List[str] = field(default_factory=list)
    venue: str = ""
    # 发现元数据
    relevance_score: float = 0.0
    discovery_depth: int = 0
    source_seed: str = ""
    discovery_path: List[str] = field(default_factory=list)

    def to_paper_node_dict(self) -> Dict[str, Any]:
        """转换为与现有 PaperNode 兼容的 dict"""
        pid = self.doi.replace("/", "_").replace(".", "_") if self.doi else hashlib.md5(self.title.encode()).hexdigest()[:10]
        return {
            "id": pid,
            "title": self.title,
            "authors": ", ".join(self.authors[:3]) + (" et al." if len(self.authors) > 3 else ""),
            "year": self.year,
            "doi": self.doi,
            "journal": self.venue or "Unknown",
            "abstract": self.abstract or "",
            "url": f"https://doi.org/{self.doi}" if self.doi else "",
            "file_attachments": "",
            "link_attachments": "",
            "tags": {
                "relevance_score": f"{self.relevance_score:.2f}",
                "discovery_depth": str(self.discovery_depth),
                "citation_count": str(self.citation_count),
                "source_seed": self.source_seed,
                "is_seed": str(self.discovery_depth == 0),
            }
        }

    def to_markdown_card(self) -> str:
        """生成 Markdown 文献卡片"""
        lines = [
            f"### {self.title}",
            "",
            f"- **DOI**: [{self.doi}](https://doi.org/{self.doi})" if self.doi else "",
            f"- **Authors**: {', '.join(self.authors[:5])}{' et al.' if len(self.authors) > 5 else ''}",
            f"- **Year**: {self.year}" if self.year else "",
            f"- **Venue**: {self.venue}" if self.venue else "",
            f"- **Citations**: {self.citation_count}",
            f"- **Relevance**: {self.relevance_score:.2f}",
            f"- **Discovery depth**: {self.discovery_depth} (from {self.source_seed[:30]}...)" if self.discovery_depth > 0 else "- **Seed paper**",
            "",
        ]
        if self.abstract:
            lines.append(f"> {self.abstract[:400]}..." if len(self.abstract) > 400 else f"> {self.abstract}")
            lines.append("")
        return "\n".join(filter(None, lines))


class SemanticScholarClient:
    """Semantic Scholar Graph API 客户端（带缓存和限速）"""

    BASE_URL = "https://api.semanticscholar.org/graph/v1"

    def __init__(
        self,
        api_key: Optional[str] = None,
        rate_limit_delay: float = 1.0,
        max_retries: int = 5,
        backoff_factor: float = 2.0,
        max_backoff: float = 90.0,
        stop_on_rate_limit: bool = True,
    ):
        self.api_key = api_key
        self.delay = rate_limit_delay
        self.max_retries = max(0, int(max_retries))
        self.backoff_factor = max(1.0, float(backoff_factor))
        self.max_backoff = max(1.0, float(max_backoff))
        self.stop_on_rate_limit = bool(stop_on_rate_limit)
        self.rate_limit_stopped = False
        self.rate_limit_reason = ""
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": "LiteratureDiscoveryAgent/1.0"
        })
        if api_key:
            self.session.headers.update({"x-api-key": api_key})
        # 简单内存缓存
        self._cache: Dict[str, Any] = {}

    def _retry_wait_seconds(self, resp, attempt: int) -> float:
        retry_after = resp.headers.get("Retry-After") if resp is not None else None
        if retry_after:
            try:
                return min(self.max_backoff, max(self.delay, float(retry_after)))
            except ValueError:
                pass
        return min(self.max_backoff, max(self.delay, self.delay * (self.backoff_factor ** attempt)))

    def _request_with_retries(self, method: str, url: str, **kwargs):
        if self.rate_limit_stopped:
            raise SemanticScholarRateLimitCircuitOpen(
                self.rate_limit_reason or "Semantic Scholar rate limit circuit is open."
            )
        last_resp = None
        for attempt in range(self.max_retries + 1):
            if attempt == 0:
                time.sleep(self.delay)
            resp = self.session.request(method, url, timeout=30, **kwargs)
            last_resp = resp
            if resp.status_code != 429:
                return resp
            if self.stop_on_rate_limit:
                self.rate_limit_stopped = True
                self.rate_limit_reason = (
                    "Semantic Scholar returned 429. Stopping remaining Semantic Scholar "
                    "requests for this run; continuing with papers already collected."
                )
                print(f"      [SS] Rate limited (429). {self.rate_limit_reason}")
                return resp
            if attempt >= self.max_retries:
                break
            wait = self._retry_wait_seconds(resp, attempt)
            print(f"      [SS] Rate limited (429), retry {attempt + 1}/{self.max_retries} after {wait:.1f}s...")
            time.sleep(wait)
        if last_resp is not None and last_resp.status_code == 429:
            self.rate_limit_stopped = True
            self.rate_limit_reason = (
                "Semantic Scholar still returned 429 after retries. "
                "Stopping remaining Semantic Scholar requests for this run."
            )
        return last_resp

    def _get(self, endpoint: str, params: Optional[Dict] = None) -> Dict:
        cache_key = f"{endpoint}:{json.dumps(params or {}, sort_keys=True)}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        url = f"{self.BASE_URL}/{endpoint}"

        try:
            resp = self._request_with_retries("GET", url, params=params)
            resp.raise_for_status()
            data = resp.json()
            self._cache[cache_key] = data
            return data
        except requests.exceptions.RequestException as e:
            print(f"      [SS] API request failed: {e}")
            raise

    def search_papers(self, query: str, fields: Optional[List[str]] = None,
                     limit: int = 100, offset: int = 0) -> Dict:
        default_fields = [
            "title", "authors", "year", "abstract", "doi",
            "citationCount", "referenceCount", "fieldsOfStudy",
            "publicationTypes", "venue", "publicationVenue", "journal"
        ]
        params = {
            "query": query,
            "fields": ",".join(fields or default_fields),
            "limit": limit,
            "offset": offset,
        }
        return self._get("paper/search", params)

    def get_paper(self, doi: str, fields: Optional[List[str]] = None) -> Dict:
        default_fields = [
            "title", "authors", "year", "abstract", "doi",
            "citationCount", "referenceCount", "fieldsOfStudy",
            "publicationTypes", "venue", "publicationVenue", "journal"
        ]
        return self._get(f"paper/DOI:{doi}", {"fields": ",".join(fields or default_fields)})

    def get_citations(self, doi: str, fields: Optional[List[str]] = None,
                     limit: int = 100, offset: int = 0) -> Dict:
        # NOTE: /citations endpoint does NOT support 'doi' field; use 'externalIds' instead
        default_fields = [
            "title", "authors", "year", "abstract", "externalIds",
            "citationCount", "fieldsOfStudy", "venue", "publicationVenue", "journal"
        ]
        params = {"fields": ",".join(fields or default_fields), "limit": limit, "offset": offset}
        return self._get(f"paper/DOI:{doi}/citations", params)

    def get_references(self, doi: str, fields: Optional[List[str]] = None,
                      limit: int = 100, offset: int = 0) -> Dict:
        # NOTE: /references endpoint does NOT support 'doi' field; use 'externalIds' instead
        default_fields = [
            "title", "authors", "year", "abstract", "externalIds",
            "citationCount", "fieldsOfStudy", "venue", "publicationVenue", "journal"
        ]
        params = {"fields": ",".join(fields or default_fields), "limit": limit, "offset": offset}
        return self._get(f"paper/DOI:{doi}/references", params)

    def get_related(self, doi: str, fields: Optional[List[str]] = None,
                   limit: int = 100) -> Dict:
        # NOTE: /related endpoint does NOT support 'doi' field; use 'externalIds' instead
        default_fields = [
            "title", "authors", "year", "abstract", "externalIds",
            "citationCount", "fieldsOfStudy", "venue", "publicationVenue", "journal"
        ]
        return self._get(f"paper/DOI:{doi}/related",
                        {"fields": ",".join(fields or default_fields), "limit": limit})

    def get_papers_batch(self, ids: List[str], fields: Optional[List[str]] = None) -> Dict:
        """批量查询论文详情（POST /paper/batch）。ids 可以是 DOI 或 CorpusId 等。"""
        default_fields = [
            "title", "authors", "year", "abstract", "externalIds",
            "citationCount", "referenceCount", "fieldsOfStudy",
            "publicationTypes", "venue", "publicationVenue", "journal"
        ]
        payload = {
            "ids": ids,
            "fields": ",".join(fields or default_fields)
        }
        url = f"{self.BASE_URL}/paper/batch"
        try:
            resp = self._request_with_retries("POST", url, json=payload)
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.RequestException as e:
            print(f"      [SS] Batch API request failed: {e}")
            raise


def _extract_doi(paper) -> str:
    """兼容 PaperNode 和 dict 的 DOI 提取"""
    if hasattr(paper, 'doi'):
        return (paper.doi or "").strip()
    if isinstance(paper, dict):
        return (paper.get('doi') or "").strip()
    return ""


def _extract_title(paper) -> str:
    if hasattr(paper, 'title'):
        return paper.title or ""
    if isinstance(paper, dict):
        return paper.get('title', "")
    return ""


def _extract_abstract(paper) -> str:
    if hasattr(paper, 'abstract'):
        return paper.abstract or ""
    if isinstance(paper, dict):
        return paper.get('abstract', "")
    return ""


def _extract_year(paper) -> int:
    if hasattr(paper, 'year'):
        return int(paper.year) if paper.year else 0
    if isinstance(paper, dict):
        y = paper.get('year', 0)
        return int(y) if y else 0
    return 0


class LiteratureDiscovery:
    """文献发现引擎：从种子集向外扩展"""

    def __init__(self, client: SemanticScholarClient,
                 relevance_threshold: float = DISCOVERY_DEFAULTS["relevance_threshold"],
                 max_papers: int = DISCOVERY_DEFAULTS["max_papers"]):
        self.client = client
        self.relevance_threshold = relevance_threshold
        self.max_papers = max_papers

        self.discovered: Dict[str, DiscoveredPaper] = {}
        self.seed_dois: Set[str] = set()
        self.citation_edges: List[Tuple[str, str]] = []
        self._visited: Set[str] = set()

    @staticmethod
    def _paper_dict_to_dataclass(paper_dict: Dict, depth: int = 0,
                                 source: str = "",
                                 path: Optional[List[str]] = None,
                                 relevance: float = 0.0) -> Optional[DiscoveredPaper]:
        # DOI may be nested under externalIds (from citations/references endpoints)
        doi = paper_dict.get("doi")
        if not doi:
            ext = paper_dict.get("externalIds") or {}
            if isinstance(ext, dict):
                doi = ext.get("DOI")
        if not doi:
            return None

        authors = []
        for a in (paper_dict.get("authors") or []):
            if isinstance(a, dict):
                name = a.get("name", "")
            else:
                name = str(a)
            if name:
                authors.append(name)

        venue = paper_dict.get("venue") or ""
        publication_venue = paper_dict.get("publicationVenue") or {}
        if not venue and isinstance(publication_venue, dict):
            venue = publication_venue.get("name") or publication_venue.get("displayName") or ""
        journal_meta = paper_dict.get("journal") or {}
        if not venue and isinstance(journal_meta, dict):
            venue = journal_meta.get("name") or journal_meta.get("displayName") or ""

        return DiscoveredPaper(
            doi=doi,
            title=paper_dict.get("title") or "Unknown",
            authors=authors,
            year=paper_dict.get("year") or 0,
            abstract=paper_dict.get("abstract"),
            citation_count=paper_dict.get("citationCount") or 0,
            reference_count=paper_dict.get("referenceCount") or 0,
            fields_of_study=paper_dict.get("fieldsOfStudy") or [],
            publication_types=paper_dict.get("publicationTypes") or [],
            venue=venue,
            relevance_score=relevance,
            discovery_depth=depth,
            source_seed=source,
            discovery_path=(path or []) + [doi],
        )

    @staticmethod
    def _assess_relevance(paper_dict: Dict, seed_text_profile: str = "") -> float:
        """评估文献与种子集的相关性（基于领域关键词 + 学科过滤）"""
        text = ((paper_dict.get("title") or "") + " " +
                (paper_dict.get("abstract") or "")).lower()

        matches = sum(1 for term in DOMAIN_KEYWORDS if term in text)
        score = min(1.0, matches / 4.0)  # 匹配4个关键词即接近满分

        # 奖励有详细 abstract 的
        abstract = paper_dict.get("abstract") or ""
        if len(abstract) > 300:
            score = min(1.0, score + 0.08)
        elif len(abstract) > 100:
            score = min(1.0, score + 0.03)
        elif not abstract:
            score *= 0.7  # 无摘要降权

        # 学科过滤：如果领域信息明确偏离化学/生物，降权
        fields = [f.lower() for f in (paper_dict.get("fieldsOfStudy") or [])]
        if fields:
            has_relevant_field = any(
                kw in f for f in fields for kw in
                ['chemistry', 'biochem', 'organic', 'chemical',
                 'biology', 'biotechnology', 'materials', 'medicine',
                 'pharmacology', 'physics']
            )
            if not has_relevant_field:
                score *= 0.4

        # 如果和种子文本直接有词汇重叠，小幅加分
        if seed_text_profile:
            seed_words = set(seed_text_profile.lower().split())
            text_words = set(text.split())
            overlap = len(seed_words & text_words)
            if overlap > 20:
                score = min(1.0, score + 0.05)

        return score

    def expand(self, seed_papers: List[Any], mode: str = "both",
               max_depth: int = 1, breadth_limit: int = 20) -> List[DiscoveredPaper]:
        """
        从种子文献向外扩展。

        Args:
            seed_papers: List[PaperNode] 或 List[dict]，需含 doi
            mode: 'citations' | 'references' | 'both' | 'related'
            max_depth: 扩展深度（1=直接邻居，2=邻居的邻居）
            breadth_limit: 每篇文献最多取多少邻居
        """
        queue: List[Tuple[Any, int]] = []
        seed_text_profile = ""

        for p in seed_papers:
            doi = _extract_doi(p)
            # 跳过无效 DOI
            if not doi or doi.lower() in ('nan', 'none'):
                continue
            self.seed_dois.add(doi)
            self._visited.add(doi)
            # 将种子也注册进 discovered
            dp = self._paper_dict_to_dataclass({
                "doi": doi,
                "title": _extract_title(p),
                "abstract": _extract_abstract(p),
                "authors": [],
                "year": _extract_year(p),
                "citationCount": 0,
                "referenceCount": 0,
                "fieldsOfStudy": [],
                "publicationTypes": [],
            }, depth=0)
            if dp:
                self.discovered[doi] = dp

            seed_text_profile += f" {_extract_title(p)} {_extract_abstract(p)}"
            queue.append((p, 0))

        print(f"\n[Discovery] 种子: {len(self.seed_dois)} 篇 | 模式: {mode} | 深度: {max_depth} | 每节点上限: {breadth_limit}")

        while queue and len(self.discovered) < self.max_papers:
            if getattr(self.client, "rate_limit_stopped", False):
                print("[Discovery] Semantic Scholar rate-limit circuit is open; stopping expansion.")
                break
            current, depth = queue.pop(0)
            if depth >= max_depth:
                continue

            doi = _extract_doi(current)
            if not doi:
                continue

            neighbors_raw: List[Dict] = []

            if mode in ("citations", "both"):
                try:
                    resp = self.client.get_citations(doi, limit=breadth_limit)
                    for item in (resp.get("data") or []):
                        paper = item.get("citingPaper", {})
                        p_doi = paper.get("doi") or (paper.get("externalIds") or {}).get("DOI")
                        if p_doi:
                            paper["doi"] = p_doi  # normalize
                            neighbors_raw.append(paper)
                            self.citation_edges.append((p_doi, doi))
                except requests.exceptions.HTTPError as e:
                    if e.response is not None and e.response.status_code == 404:
                        pass  # SS 没收录这篇，静默跳过
                    else:
                        print(f"      [!] citations fetch failed for {doi[:30]}: {e}")
                except Exception as e:
                    print(f"      [!] citations fetch failed for {doi[:30]}: {e}")

            if mode in ("references", "both"):
                if getattr(self.client, "rate_limit_stopped", False):
                    break
                try:
                    resp = self.client.get_references(doi, limit=breadth_limit)
                    for item in (resp.get("data") or []):
                        paper = item.get("citedPaper", {})
                        p_doi = paper.get("doi") or (paper.get("externalIds") or {}).get("DOI")
                        if p_doi:
                            paper["doi"] = p_doi  # normalize
                            neighbors_raw.append(paper)
                            self.citation_edges.append((doi, p_doi))
                except requests.exceptions.HTTPError as e:
                    if e.response is not None and e.response.status_code == 404:
                        pass  # SS 没收录这篇，静默跳过
                    else:
                        print(f"      [!] references fetch failed for {doi[:30]}: {e}")
                except Exception as e:
                    print(f"      [!] references fetch failed for {doi[:30]}: {e}")

            if mode == "related":
                if getattr(self.client, "rate_limit_stopped", False):
                    break
                try:
                    resp = self.client.get_related(doi, limit=breadth_limit)
                    for paper in (resp.get("data") or []):
                        p_doi = paper.get("doi") or (paper.get("externalIds") or {}).get("DOI")
                        if p_doi:
                            paper["doi"] = p_doi  # normalize
                            neighbors_raw.append(paper)
                except requests.exceptions.HTTPError as e:
                    if e.response is not None and e.response.status_code == 404:
                        pass
                    else:
                        print(f"      [!] related fetch failed for {doi[:30]}: {e}")
                except Exception as e:
                    print(f"      [!] related fetch failed for {doi[:30]}: {e}")

            for paper in neighbors_raw:
                p_doi = paper.get("doi")
                if not p_doi or p_doi in self._visited:
                    continue

                relevance = self._assess_relevance(paper, seed_text_profile)
                if relevance < self.relevance_threshold:
                    continue

                self._visited.add(p_doi)

                # 构建发现路径
                source_dp = self.discovered.get(doi)
                path = source_dp.discovery_path if source_dp else []

                dp = self._paper_dict_to_dataclass(
                    paper, depth=depth + 1, source=doi, path=path, relevance=relevance
                )
                if dp:
                    self.discovered[p_doi] = dp
                    queue.append((paper, depth + 1))

        discovered_count = len(self.discovered) - len(self.seed_dois)
        print(f"[Discovery] 完成。总计 {len(self.discovered)} 篇（种子 {len(self.seed_dois)} + 新发现 {discovered_count}）")
        return list(self.discovered.values())

    def search(self, query: str, limit: int = 50) -> List[DiscoveredPaper]:
        """关键词搜索并加入发现集"""
        try:
            resp = self.client.search_papers(query, limit=limit)
            new_papers = []
            for paper in (resp.get("data") or []):
                p_doi = paper.get("doi") or (paper.get("externalIds") or {}).get("DOI")
                if not p_doi or p_doi in self._visited:
                    continue
                relevance = self._assess_relevance(paper, "")
                if relevance >= self.relevance_threshold:
                    self._visited.add(p_doi)
                    paper["doi"] = p_doi  # normalize
                    dp = self._paper_dict_to_dataclass(paper, depth=0, source="search: " + query, relevance=relevance)
                    if dp:
                        self.discovered[p_doi] = dp
                        new_papers.append(dp)
            print(f"[Discovery] 搜索 '{query}' → {len(new_papers)} 篇相关")
            return new_papers
        except SemanticScholarRateLimitCircuitOpen as e:
            print(f"[Discovery] Search skipped after Semantic Scholar 429: {e}")
            return []
        except Exception as e:
            print(f"[Discovery] 搜索失败: {e}")
            return []


class CitationNetworkAnalyzer:
    """引用网络分析器：基于真实 paper→paper 有向图"""

    def __init__(self):
        self.G = nx.DiGraph()

    def build(self, discovery: LiteratureDiscovery):
        """从 LiteratureDiscovery 结果构建网络"""
        for doi, paper in discovery.discovered.items():
            self.G.add_node(
                doi,
                title=paper.title,
                year=paper.year,
                citation_count=paper.citation_count,
                venue=paper.venue,
                is_seed=(doi in discovery.seed_dois),
                relevance=paper.relevance_score,
                depth=paper.discovery_depth,
            )

        for citing, cited in discovery.citation_edges:
            if citing in self.G and cited in self.G:
                self.G.add_edge(citing, cited)

        print(f"[Network] 网络: {self.G.number_of_nodes()} 节点, {self.G.number_of_edges()} 边")

    def enrich_seed_citations(self, discovery: LiteratureDiscovery, client: SemanticScholarClient):
        """用 batch API 补充种子文献的真实引用数和期刊信息"""
        seed_dois = list(discovery.seed_dois)
        if not seed_dois:
            return
        print(f"[Network] 正在为 {len(seed_dois)} 篇种子补充引用数...")
        # 每次批量查询最多 100 篇
        batch_size = 100
        enriched = 0
        for i in range(0, len(seed_dois), batch_size):
            batch = seed_dois[i:i+batch_size]
            ids = [f"DOI:{d}" for d in batch]
            try:
                data = client.get_papers_batch(ids)
                for paper in (data or []):
                    ext = paper.get("externalIds") or {}
                    doi = ext.get("DOI")
                    if not doi:
                        continue
                    if doi in self.G:
                        self.G.nodes[doi]["citation_count"] = paper.get("citationCount") or 0
                        venue = LiteratureDiscovery._paper_dict_to_dataclass(paper)
                        self.G.nodes[doi]["venue"] = venue.venue if venue else (paper.get("venue") or "")
                        # 同时更新 discovery 对象中的种子信息
                        if doi in discovery.discovered:
                            discovery.discovered[doi].citation_count = paper.get("citationCount") or 0
                            discovery.discovered[doi].venue = venue.venue if venue else (paper.get("venue") or "")
                        enriched += 1
            except Exception as e:
                print(f"      [!] Batch enrichment failed: {e}")
        print(f"[Network] 已补充 {enriched}/{len(seed_dois)} 篇种子的引用数")

    def pagerank(self, top_n: int = 20, exclude_seeds: bool = False) -> List[Tuple[str, float]]:
        """PageRank 找核心文献"""
        if self.G.number_of_nodes() == 0:
            return []
        pr = nx.pagerank(self.G)
        items = [(n, s) for n, s in pr.items()
                 if not exclude_seeds or not self.G.nodes[n].get('is_seed', False)]
        items.sort(key=lambda x: x[1], reverse=True)
        return items[:top_n]

    def betweenness_bridges(self, top_n: int = 20) -> List[Tuple[str, float]]:
        """Betweenness Centrality 找桥梁文献（连接不同子领域的文献）"""
        if self.G.number_of_nodes() < 3:
            return []
        try:
            bc = nx.betweenness_centrality(self.G)
            # 优先非种子节点
            candidates = [(n, s) for n, s in bc.items()
                         if not self.G.nodes[n].get('is_seed', False)]
            candidates.sort(key=lambda x: x[1], reverse=True)
            return candidates[:top_n]
        except Exception as e:
            print(f"[Network] Betweenness 计算失败: {e}")
            return []

    def emerging_clusters(self, year_threshold: int = 2022, min_size: int = 3):
        """识别近年新兴聚类（弱连通分量）"""
        recent = [n for n, attr in self.G.nodes(data=True)
                 if attr.get('year', 0) >= year_threshold
                 and not attr.get('is_seed', False)]
        if len(recent) < min_size:
            return []

        subG = self.G.subgraph(recent).to_undirected()
        clusters = list(nx.connected_components(subG))

        scored = []
        for cluster in clusters:
            if len(cluster) < min_size:
                continue
            nodes = [self.G.nodes[n] for n in cluster]
            scored.append({
                'size': len(cluster),
                'avg_citations': np.mean([n.get('citation_count', 0) for n in nodes]),
                'avg_relevance': np.mean([n.get('relevance', 0) for n in nodes]),
                'papers': [{'doi': n,
                           'title': self.G.nodes[n].get('title', ''),
                           'year': self.G.nodes[n].get('year', 0),
                           'citations': self.G.nodes[n].get('citation_count', 0)}
                          for n in cluster]
            })

        scored.sort(key=lambda x: (x['size'], x['avg_citations']), reverse=True)
        return scored

    def detect_blindspots(self, seed_dois: Set[str]) -> Dict[str, List[Dict]]:
        """发现种子集的盲区：高被引缺失 + 新兴热点"""
        blindspots = {
            'missing_landmarks': [],
            'emerging_hotspots': [],
        }

        for n, attr in self.G.nodes(data=True):
            if attr.get('is_seed'):
                continue
            cc = attr.get('citation_count', 0)
            year = attr.get('year', 0)

            if cc >= 150:
                blindspots['missing_landmarks'].append({
                    'doi': n,
                    'title': attr.get('title', ''),
                    'citations': cc,
                    'year': year,
                    'reason': '引用网络中高被引枢纽（ landmarks ），但不在你的收藏中',
                })

            if year >= 2022 and cc >= 15:
                blindspots['emerging_hotspots'].append({
                    'doi': n,
                    'title': attr.get('title', ''),
                    'year': year,
                    'citations': cc,
                    'reason': '近年发表且引用增长迅速，可能是新兴热点',
                })

        blindspots['missing_landmarks'].sort(key=lambda x: x['citations'], reverse=True)
        blindspots['emerging_hotspots'].sort(key=lambda x: x['citations'], reverse=True)
        return blindspots

    def network_stats(self) -> Dict:
        """网络统计摘要"""
        n = self.G.number_of_nodes()
        m = self.G.number_of_edges()
        seeds = sum(1 for _, attr in self.G.nodes(data=True) if attr.get('is_seed'))

        stats = {
            'total_nodes': n,
            'total_edges': m,
            'seed_nodes': seeds,
            'discovered_nodes': n - seeds,
            'density': nx.density(self.G) if n > 1 else 0,
        }

        try:
            wcc = list(nx.weakly_connected_components(self.G))
            stats['weakly_connected_components'] = len(wcc)
            stats['largest_component_size'] = len(max(wcc, key=len)) if wcc else 0
        except Exception:
            pass

        return stats

    def export_gexf(self, path: str):
        """导出为 Gephi 可用的 GEXF 格式"""
        nx.write_gexf(self.G, path)
        print(f"[Network] GEXF 导出: {path}")

    def export_edgelist(self, path: str):
        """导出边列表 CSV"""
        import csv
        with open(path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['citing_doi', 'cited_doi', 'citing_title', 'cited_title'])
            for u, v in self.G.edges():
                writer.writerow([
                    u, v,
                    self.G.nodes[u].get('title', ''),
                    self.G.nodes[v].get('title', '')
                ])
        print(f"[Network] 边列表导出: {path}")

    def export_html(self, path: str):
        """导出可在线浏览的引文网络 HTML（无需 pyvis）。"""
        nodes = []
        for node, attr in self.G.nodes(data=True):
            nodes.append({
                "id": node,
                "title": attr.get("title", ""),
                "year": attr.get("year", 0),
                "citations": attr.get("citation_count", 0) or 0,
                "venue": attr.get("venue", ""),
                "is_seed": bool(attr.get("is_seed")),
                "relevance": attr.get("relevance", 0) or 0,
                "depth": attr.get("depth", 0),
            })
        edges = [{"source": u, "target": v} for u, v in self.G.edges()]
        data = json.dumps({"nodes": nodes, "edges": edges}, ensure_ascii=False)
        html = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Citation Network</title>
<style>
body{{margin:0;font-family:Arial,"Microsoft YaHei",sans-serif;background:#f7f7f4;color:#222}}
header{{padding:14px 18px;background:#fff;border-bottom:1px solid #ddd;display:flex;gap:10px;align-items:center;flex-wrap:wrap}}
input{{padding:8px 10px;border:1px solid #ccc;border-radius:6px;min-width:280px}}
#wrap{{display:grid;grid-template-columns:1fr 360px;height:calc(100vh - 62px)}}
canvas{{width:100%;height:100%;background:#fbfbf8}}
aside{{border-left:1px solid #ddd;background:#fff;padding:14px;overflow:auto}}
.muted{{color:#667;font-size:13px}} .chip{{display:inline-block;background:#edf1ed;border-radius:999px;padding:3px 7px;margin:2px;font-size:12px}}
a{{color:#356a88}}
</style></head>
<body><header><strong>Citation Network</strong><input id="q" placeholder="搜索标题/DOI/期刊"><span id="stats" class="muted"></span></header>
<div id="wrap"><canvas id="c"></canvas><aside id="info"><div class="muted">点击节点查看文献。</div></aside></div>
<script>
const DATA={data};
let nodes=DATA.nodes, edges=DATA.edges, pos=new Map(), selected=null;
const c=document.getElementById('c'), ctx=c.getContext('2d'), info=document.getElementById('info');
function fit(){{const r=c.getBoundingClientRect(); c.width=r.width*devicePixelRatio; c.height=r.height*devicePixelRatio; ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);}}
function layout(){{const r=c.getBoundingClientRect(); pos.clear(); const sorted=[...nodes].sort((a,b)=>(b.citations||0)-(a.citations||0)); sorted.forEach((n,i)=>{{const ring=n.is_seed?0.38:0.25+0.32*((n.depth||1)%3)/3; const ang=2*Math.PI*i/Math.max(1,sorted.length); const rad=Math.min(r.width,r.height)*ring; pos.set(n.id,[r.width/2+Math.cos(ang)*rad,r.height/2+Math.sin(ang)*rad]);}})}}
function draw(){{fit();layout();const r=c.getBoundingClientRect();ctx.clearRect(0,0,r.width,r.height);ctx.strokeStyle='rgba(40,50,50,.16)';edges.forEach(e=>{{const a=pos.get(e.source),b=pos.get(e.target);if(!a||!b)return;ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();}});nodes.forEach(n=>{{const p=pos.get(n.id),s=Math.max(3,Math.min(17,Math.log((n.citations||0)+1)*2.2));ctx.beginPath();ctx.fillStyle=n.is_seed?'#356a88':'#48735c';ctx.arc(p[0],p[1],s,0,Math.PI*2);ctx.fill();if(selected&&selected.id===n.id){{ctx.strokeStyle='#9a6b2f';ctx.lineWidth=3;ctx.stroke();ctx.lineWidth=1;}}}});document.getElementById('stats').textContent=`${{nodes.length}} nodes · ${{edges.length}} edges`;}}
function show(n){{selected=n;info.innerHTML=`<h2>${{esc(n.title||n.id)}}</h2><div class="muted">${{n.year||''}} · ${{esc(n.venue||'')}}</div><p><a href="https://doi.org/${{esc(n.id)}}" target="_blank">${{esc(n.id)}}</a></p><div><span class="chip">${{n.is_seed?'seed':'discovered'}}</span><span class="chip">${{n.citations||0}} citations</span><span class="chip">depth ${{n.depth||0}}</span><span class="chip">rel ${{Number(n.relevance||0).toFixed(2)}}</span></div>`;draw();}}
function esc(s){{return String(s??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));}}
c.onclick=e=>{{const r=c.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;let best=null,bd=20;nodes.forEach(n=>{{const p=pos.get(n.id);if(!p)return;const d=Math.hypot(p[0]-x,p[1]-y);if(d<bd){{bd=d;best=n;}}}});if(best)show(best);}};
document.getElementById('q').oninput=e=>{{const q=e.target.value.toLowerCase();nodes=DATA.nodes.filter(n=>[n.title,n.id,n.venue].join(' ').toLowerCase().includes(q));const keep=new Set(nodes.map(n=>n.id));edges=DATA.edges.filter(e=>keep.has(e.source)&&keep.has(e.target));draw();}};
window.onresize=draw;draw();
</script></body></html>"""
        with open(path, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"[Network] HTML 导出: {path}")


class DiscoveryReportGenerator:
    """为 Discovery 结果生成独立报告"""

    def __init__(self, discovery: LiteratureDiscovery, analyzer: CitationNetworkAnalyzer,
                 output_dir: Path):
        self.discovery = discovery
        self.analyzer = analyzer
        self.output_dir = output_dir

    def generate(self):
        self._write_discovery_summary()
        self._write_network_analysis()
        self._write_blindspot_report()
        self._write_discovered_papers_cards()
        self._plot_citation_distribution()

    def _write_discovery_summary(self):
        stats = self.analyzer.network_stats()
        lines = [
            "# 文献发现与网络扩展报告",
            "",
            f"> 生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}" if pd else f"> 生成时间: {time.strftime('%Y-%m-%d %H:%M')}",
            "",
            "## 一、扩展概览",
            "",
            f"- **种子文献**: {stats['seed_nodes']} 篇（来自你的 Zotero CSV）",
            f"- **新发现文献**: {stats['discovered_nodes']} 篇（通过 Semantic Scholar API 扩展）",
            f"- **引用网络规模**: {stats['total_nodes']} 节点, {stats['total_edges']} 条有向边",
            f"- **网络密度**: {stats['density']:.4f}",
            f"- **弱连通分量**: {stats.get('weakly_connected_components', 'N/A')}",
            f"- **最大连通分量**: {stats.get('largest_component_size', 'N/A')} 节点",
            "",
            "### 发现深度分布",
            "",
        ]

        depth_counts = Counter(p.discovery_depth for p in self.discovery.discovered.values())
        for d in sorted(depth_counts):
            label = "种子" if d == 0 else f"第 {d} 层扩展"
            lines.append(f"- **{label}**: {depth_counts[d]} 篇")

        lines.append("")
        lines.append("---")
        lines.append("")

        path = self.output_dir / "discovery_summary.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"[Report] 发现摘要: {path}")

    def _write_network_analysis(self):
        lines = ["# 引用网络拓扑分析", ""]

        # 解释
        lines.append("> **PageRank (PR Score)**：引用网络中的『稳态影响力』。可以理解为：如果一只随机游走的学术兔子在网络中跳来跳去，它最终停留在某篇文献上的概率。PR 越高，说明该文献被越多重要文献引用，是领域的枢纽/综述/里程碑。")
        lines.append("> ")
        lines.append("> **Betweenness Centrality (BC Score)**：桥梁价值。如果移除该文献，网络中多少最短路径会断裂？BC 高的文献通常是跨子领域的『连接器』，把光化学和酶工程、或合成和方法学连在一起。")
        lines.append("")

        # PageRank
        lines.append("## 二、网络核心（PageRank Top 20）")
        lines.append("")
        lines.append("| Rank | DOI | Title | Year | Venue | Citations | PR Score | Seed? |")
        lines.append("|------|-----|-------|------|-------|-----------|----------|-------|")
        for i, (doi, score) in enumerate(self.analyzer.pagerank(top_n=20), 1):
            attr = self.analyzer.G.nodes[doi]
            title = attr.get('title', '')[:55] + "..." if len(attr.get('title', '')) > 55 else attr.get('title', '')
            venue = (attr.get('venue', '') or '')[:20]
            cc = attr.get('citation_count', 0) or 0
            is_seed = "✅" if attr.get('is_seed') else ""
            lines.append(f"| {i} | [{doi[:20]}...](https://doi.org/{doi}) | {title} | {attr.get('year', '')} | {venue} | {cc} | {score:.4f} | {is_seed} |")
        lines.append("")

        # Bridges
        lines.append("## 三、桥梁文献（Betweenness Centrality Top 15）")
        lines.append("")
        lines.append("| Rank | DOI | Title | Year | Venue | Citations | BC Score |")
        lines.append("|------|-----|-------|------|-------|-----------|----------|")
        for i, (doi, score) in enumerate(self.analyzer.betweenness_bridges(top_n=15), 1):
            attr = self.analyzer.G.nodes[doi]
            title = attr.get('title', '')[:55] + "..." if len(attr.get('title', '')) > 55 else attr.get('title', '')
            venue = (attr.get('venue', '') or '')[:20]
            cc = attr.get('citation_count', 0) or 0
            lines.append(f"| {i} | [{doi[:20]}...](https://doi.org/{doi}) | {title} | {attr.get('year', '')} | {venue} | {cc} | {score:.4f} |")
        lines.append("")

        # Emerging clusters
        lines.append("## 四、新兴聚类（2022 年以来）")
        lines.append("")
        clusters = self.analyzer.emerging_clusters(year_threshold=2022, min_size=3)
        if not clusters:
            lines.append("> 未发现足够大的新兴聚类（聚类大小 < 3）")
        for i, c in enumerate(clusters[:5], 1):
            lines.append(f"### 聚类 {i}（{c['size']} 篇, 均引 {c['avg_citations']:.1f}）")
            lines.append("")
            for p in c['papers'][:8]:
                lines.append(f"- [{p['title'][:70]}...](https://doi.org/{p['doi']}) ({p['year']}, {p['citations']} citations)")
            if len(c['papers']) > 8:
                lines.append(f"- ... 以及另外 {len(c['papers']) - 8} 篇")
            lines.append("")

        path = self.output_dir / "network_analysis.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"[Report] 网络分析: {path}")

    def _write_blindspot_report(self):
        blindspots = self.analyzer.detect_blindspots(self.discovery.seed_dois)
        lines = ["# 文献盲区与机会分析", "", "> 基于数据驱动的盲区检测，替代静态规则检查", ""]

        lines.append("## 缺失的里程碑（Missing Landmarks）")
        lines.append("")
        lines.append("这些高被引文献在你的引用网络中处于枢纽位置，但你尚未收藏。它们通常是领域奠基性工作。")
        lines.append("")
        if not blindspots['missing_landmarks']:
            lines.append("> 未发现明显缺失的里程碑文献（阈值: citations >= 150）")
        else:
            lines.append("| Title | DOI | Year | Venue | Citations |")
            lines.append("|-------|-----|------|-------|-----------|")
            for item in blindspots['missing_landmarks'][:15]:
                title = item['title'][:50] + "..." if len(item['title']) > 50 else item['title']
                # 尝试从网络节点获取 venue
                venue = ""
                if item['doi'] in self.analyzer.G:
                    venue = (self.analyzer.G.nodes[item['doi']].get('venue') or '')[:20]
                lines.append(f"| {title} | [{item['doi'][:20]}...](https://doi.org/{item['doi']}) | {item['year']} | {venue} | {item['citations']} |")
            lines.append("")

        lines.append("## 新兴热点盲区（Emerging Hotspots）")
        lines.append("")
        lines.append("这些近年发表的文献已经获得显著引用，可能是快速发展的子方向，但你尚未关注。")
        lines.append("")
        if not blindspots['emerging_hotspots']:
            lines.append("> 未发现新兴热点盲区（阈值: 2022年后发表，citations >= 15）")
        else:
            lines.append("| Title | DOI | Year | Venue | Citations |")
            lines.append("|-------|-----|------|-------|-----------|")
            for item in blindspots['emerging_hotspots'][:15]:
                title = item['title'][:50] + "..." if len(item['title']) > 50 else item['title']
                venue = ""
                if item['doi'] in self.analyzer.G:
                    venue = (self.analyzer.G.nodes[item['doi']].get('venue') or '')[:20]
                lines.append(f"| {title} | [{item['doi'][:20]}...](https://doi.org/{item['doi']}) | {item['year']} | {venue} | {item['citations']} |")
            lines.append("")

        path = self.output_dir / "blindspot_analysis.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"[Report] 盲区分析: {path}")

    def _write_discovered_papers_cards(self):
        """为所有新发现的文献生成卡片"""
        new_papers = [p for p in self.discovery.discovered.values() if p.doi not in self.discovery.seed_dois]
        if not new_papers:
            return

        # 按引用数排序
        new_papers.sort(key=lambda p: p.citation_count, reverse=True)

        lines = [f"# 新发现文献卡片（共 {len(new_papers)} 篇）", ""]
        lines.append("按被引次数降序排列。")
        lines.append("")

        for p in new_papers:
            lines.append(p.to_markdown_card())
            lines.append("")

        path = self.output_dir / "discovered_papers.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"[Report] 文献卡片: {path}")

    def _plot_citation_distribution(self):
        """种子 vs 发现文献的引用分布对比"""
        if plt is None:
            return

        seed_cc = [p.citation_count for p in self.discovery.discovered.values() if p.doi in self.discovery.seed_dois]
        disc_cc = [p.citation_count for p in self.discovery.discovered.values() if p.doi not in self.discovery.seed_dois]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # 直方图
        axes[0].hist([seed_cc, disc_cc], bins=30, label=['Seed', 'Discovered'],
                    color=['#5B7C99', '#C9B8A8'], alpha=0.8)
        axes[0].set_xlabel('Citation Count')
        axes[0].set_ylabel('Number of Papers')
        axes[0].set_title('Citation Count Distribution')
        axes[0].legend()
        axes[0].set_yscale('log')

        # 年份分布
        seed_years = [p.year for p in self.discovery.discovered.values() if p.doi in self.discovery.seed_dois and p.year > 1990]
        disc_years = [p.year for p in self.discovery.discovered.values() if p.doi not in self.discovery.seed_dois and p.year > 1990]
        all_years = list(range(min(seed_years + disc_years), max(seed_years + disc_years) + 1))
        seed_counts = Counter(seed_years)
        disc_counts = Counter(disc_years)

        axes[1].bar([y - 0.2 for y in all_years], [seed_counts.get(y, 0) for y in all_years],
                   width=0.4, label='Seed', color='#5B7C99', alpha=0.8)
        axes[1].bar([y + 0.2 for y in all_years], [disc_counts.get(y, 0) for y in all_years],
                   width=0.4, label='Discovered', color='#C9B8A8', alpha=0.8)
        axes[1].set_xlabel('Year')
        axes[1].set_ylabel('Number of Papers')
        axes[1].set_title('Publication Year Distribution')
        axes[1].legend()

        plt.tight_layout()
        path = self.output_dir / "citation_distribution.png"
        plt.savefig(path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[Plot] 引用分布图: {path}")

    def _plot_discovery_sunburst(self):
        """发现来源的桑基/流向简化图：展示种子→扩展的流向"""
        if plt is None:
            return

        # 统计每层发现的来源分布
        source_counter = Counter()
        for p in self.discovery.discovered.values():
            if p.doi in self.discovery.seed_dois:
                continue
            source = p.source_seed[:20] + "..." if len(p.source_seed) > 20 else p.source_seed
            source_counter[source] += 1

        if not source_counter:
            return

        fig, ax = plt.subplots(figsize=(10, 8))
        items = source_counter.most_common(15)
        labels = [item[0] for item in items]
        values = [item[1] for item in items]

        colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
        ax.barh(range(len(labels)), values, color=colors)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel('Discovered Papers Count')
        ax.set_title('Top 15 Seed Sources for Discovered Papers')
        plt.tight_layout()
        path = self.output_dir / "discovery_sources.png"
        plt.savefig(path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[Plot] 发现来源图: {path}")


def discovery_cli():
    parser = argparse.ArgumentParser(description="Literature Discovery via Semantic Scholar")
    parser.add_argument("--input", "-i", required=True, help="种子文献 CSV 或 JSON（与主脚本格式兼容）")
    parser.add_argument("--output", "-o", default="./discovery_output", help="输出目录")
    parser.add_argument("--ss-api-key", help="Semantic Scholar API Key（可选）")
    parser.add_argument("--mode", choices=["citations", "references", "both", "related"],
                       default="both", help="扩展模式")
    parser.add_argument("--max-depth", type=int, default=1, help="扩展深度")
    parser.add_argument("--breadth-limit", type=int, default=20, help="每节点邻居上限")
    parser.add_argument("--relevance-threshold", type=float, default=0.25, help="相关性阈值")
    parser.add_argument("--max-papers", type=int, default=500, help="总文献上限")
    parser.add_argument("--rate-limit", type=float, default=1.0, help="API request interval in seconds")
    parser.add_argument("--ss-max-retries", type=int, default=0, help="Semantic Scholar 429 retry count; default 0 avoids long waits")
    parser.add_argument("--ss-continue-after-429", action="store_true", help="Continue later Semantic Scholar requests after 429")
    parser.add_argument("--search-limit", type=int, default=25, help="每条关键词搜索最多取多少篇")
    parser.add_argument("--max-search-queries", type=int, default=4, help="Max supplemental keyword searches; 0 disables them")
    parser.add_argument("--search", action="append", help="额外关键词搜索（可多次）")
    parser.add_argument("--export-gexf", action="store_true", help="导出 Gephi GEXF")
    parser.add_argument("--export-edgelist", action="store_true", help="导出边列表 CSV")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 读取种子
    import pandas as pd
    df = pd.read_csv(args.input)
    seed_papers = []
    for _, row in df.iterrows():
        doi = str(row.get('DOI', '')).strip()
        if doi:
            seed_papers.append({
                'doi': doi,
                'title': str(row.get('Title', '')),
                'abstract': str(row.get('Abstract Note', '')),
                'year': int(row.get('Publication Year', 0)) if pd.notna(row.get('Publication Year')) else 0,
            })

    print(f"[CLI] 加载种子 {len(seed_papers)} 篇")

    # 初始化
    client = SemanticScholarClient(
        api_key=args.ss_api_key,
        rate_limit_delay=args.rate_limit,
        max_retries=args.ss_max_retries,
        stop_on_rate_limit=not args.ss_continue_after_429,
    )
    discovery = LiteratureDiscovery(client,
                                    relevance_threshold=args.relevance_threshold,
                                    max_papers=args.max_papers)

    # 扩展
    discovery.expand(seed_papers, mode=args.mode, max_depth=args.max_depth,
                    breadth_limit=args.breadth_limit)

    # 额外搜索
    if args.search and args.max_search_queries != 0:
        search_queries = args.search[:args.max_search_queries] if args.max_search_queries > 0 else args.search
        if args.max_search_queries > 0 and len(args.search) > len(search_queries):
            print(f"[Discovery] Skipping {len(args.search) - len(search_queries)} extra keyword searches due to --max-search-queries={args.max_search_queries}")
        for q in search_queries:
            if getattr(client, "rate_limit_stopped", False):
                print("[Discovery] Semantic Scholar rate-limit circuit is open; skipping remaining keyword searches.")
                break
            discovery.search(q, limit=args.search_limit)
    elif args.search and args.max_search_queries == 0:
        print("[Discovery] Supplemental Semantic Scholar keyword searches disabled by --max-search-queries=0")

    # 网络分析
    analyzer = CitationNetworkAnalyzer()
    analyzer.build(discovery)

    # 报告
    reporter = DiscoveryReportGenerator(discovery, analyzer, output_dir)
    reporter.generate()

    # 导出
    if args.export_gexf:
        analyzer.export_gexf(str(output_dir / "citation_network.gexf"))
    if args.export_edgelist:
        analyzer.export_edgelist(str(output_dir / "citation_edges.csv"))

    print("\n[CLI] Discovery 完成。输出:")
    for f in sorted(output_dir.iterdir()):
        print(f"  - {f.name}")


if __name__ == "__main__":
    discovery_cli()
