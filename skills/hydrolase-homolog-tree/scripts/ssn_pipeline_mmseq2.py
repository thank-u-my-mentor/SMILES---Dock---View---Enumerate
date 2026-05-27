#!/usr/bin/env python3
# python ssn_pipeline_mmseq2_fixed.py all.fa --core-fasta core.fa --output alpha_beta_ssn_v25 --max-rep 1500 --cdhit-identity 0.35 --cdhit-coverage 0.75 --max-recon 4000

import json
import csv
import os
import subprocess
import shutil
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set
from collections import defaultdict
from dataclasses import dataclass
import argparse
import networkx as nx
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


class Config:
    BASE_DIR = Path(".")
    OUTPUT_DIR = BASE_DIR / "ssn_analysis"
    
    # === 两个永久缓存 ===
    KINGDOM_CACHE_FILE = BASE_DIR / "kingdom_global_cache.json"      # organism → kingdom
    UNIPROT_ORG_CACHE_FILE = BASE_DIR / "uniprot_to_organism.json"  # ID → organism

    CLUSTALO = shutil.which("clustalo") or "/mnt/l/WSL/softwares/conda_envs/md/bin/clustalo"
    CDHIT = shutil.which("cd-hit") or shutil.which("cdhit")
    SSN_IDENTITY_THRESHOLD = 0.25
    LENGTH_MIN = 200
    LENGTH_MAX = 600
    MSA_THREADS = 8
    TOP_K_EDGES = 30
    KIMI_API_KEY = os.getenv("KIMI_API_KEY")
    KIMI_MODEL = "moonshot-v1-128k"


@dataclass
class SequenceRecord:
    id: str
    uniprot_id: str
    species_code: str
    sequence: str
    kingdom: str = "Unknown"
    organism: str = "Unknown"
    phylum_class: str = "Unknown"
    length: int = 0
    cluster_id: int = -1
    is_representative: bool = False
    is_core: bool = False
    community: int = -1

    def __post_init__(self):
        self.length = len(self.sequence)


@dataclass
class SSNEdge:
    source: str
    target: str
    identity: float


class CHRSSNPipelineCDHIT:
    @staticmethod
    def sanitize_id(raw_id: str) -> str:
        safe = re.sub(r'[^a-zA-Z0-9_.:-]', '_', raw_id)
        return '_' + safe if safe and safe[0].isdigit() else (safe or '_node')

    def __init__(self, input_fasta: Path, output_dir: Optional[Path] = None,
                 max_rep: int = 2500, cdhit_identity: float = 0.35,
                 cdhit_coverage: float = 0.75, max_recon: int = 4000,
                 core_fasta: Optional[Path] = None):
        
        self.input_fasta = Path(input_fasta)
        self.config = Config()
        self.output_dir = output_dir or self.config.OUTPUT_DIR
        self.output_dir.mkdir(exist_ok=True, parents=True)

        # === 两个永久缓存 ===
        self.kingdom_cache_file = Config.KINGDOM_CACHE_FILE
        self.kingdom_cache = {}
        self.load_global_kingdom_cache()

        self.uniprot_org_cache_file = Config.UNIPROT_ORG_CACHE_FILE
        self.uniprot_org_cache = {}
        self.load_uniprot_org_cache()

        self.max_rep = max_rep
        self.cdhit_identity = cdhit_identity
        self.cdhit_coverage = cdhit_coverage
        self.max_recon = max_recon

        self.sequences: List[SequenceRecord] = []
        self.representative_seqs: List[SequenceRecord] = []
        self.edges: List[SSNEdge] = []
        self.nodes: Dict[str, dict] = {}
        
        self.core_fasta = Path(core_fasta) if core_fasta else None
        self.core_ids: Set[str] = set()

    # ====================== 缓存 ======================
    def load_global_kingdom_cache(self):
        if self.kingdom_cache_file.exists():
            try:
                with open(self.kingdom_cache_file, 'r', encoding='utf-8') as f:
                    raw = json.load(f)
                self.kingdom_cache = {k: v for k, v in raw.items() if k.startswith("org:")}
                print(f" ✓ 加载 organism→kingdom 缓存 ({len(self.kingdom_cache)} 条)")
            except Exception as e:
                print(f" ⚠ 加载 kingdom 缓存失败: {e}")
                self.kingdom_cache = {}

    def load_uniprot_org_cache(self):
        if self.uniprot_org_cache_file.exists():
            try:
                with open(self.uniprot_org_cache_file, 'r', encoding='utf-8') as f:
                    self.uniprot_org_cache = json.load(f)
                print(f" ✓ 加载 UniProt ID→organism 缓存 ({len(self.uniprot_org_cache)} 条)")
            except Exception as e:
                print(f" ⚠ 加载 ID 缓存失败: {e}")
                self.uniprot_org_cache = {}

    def save_caches(self):
        try:
            with open(self.kingdom_cache_file, 'w', encoding='utf-8') as f:
                json.dump(self.kingdom_cache, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f" ⚠ 保存 kingdom 缓存失败: {e}")
        try:
            with open(self.uniprot_org_cache_file, 'w', encoding='utf-8') as f:
                json.dump(self.uniprot_org_cache, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f" ⚠ 保存 ID→organism 缓存失败: {e}")

    # ====================== 工具方法 ======================
    def _clean_uniprot_id(self, raw: str) -> str:
        uid = str(raw).strip()
        uid = re.sub(r'^UniRef\d+_', '', uid)
        uid = re.sub(r'\.\d+$', '', uid)
        match = re.search(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', uid)
        return match.group(1) if match else uid

    def _standardize_organism_name(self, name: str) -> str:
        if not name or name == "Unknown":
            return "Unknown"
        name = re.sub(r'\s+var\.\s*|\s+subsp\.\s*|\s+cf\.\s*|\s+cv\.\s*|\s+f\.\s*', ' ', name, flags=re.IGNORECASE)
        name = re.sub(r'[^\w\s]', '', name)
        name = re.sub(r'\s+', ' ', name).strip().lower()
        return name
    def _parse_lineage(self, lineage) -> tuple:
        """超级稳健版：支持 UniProtKB / UniParc 各种格式"""
        if not lineage:
            return "Unknown", "Unknown"

        # 统一转为干净的字符串列表
        lineage_parts = []
        for item in lineage:
            if isinstance(item, dict):
                # UniProtKB / UniParc 常见字段
                name = (item.get("scientificName") or 
                       item.get("name") or 
                       item.get("taxonName") or 
                       str(item))
                lineage_parts.append(str(name).strip())
            else:
                lineage_parts.append(str(item).strip())

        if not lineage_parts:
            return "Unknown", "Unknown"

        full_lineage = " ".join(lineage_parts).lower()
        superkingdom = lineage_parts[0]
        phylum = lineage_parts[1] if len(lineage_parts) > 1 else "Unknown"

        kingdom = "Unknown"

        # ==================== 1. Kingdom 识别（优先级最高） ====================
        if "bacteria" in full_lineage:
            kingdom = "Bacteria"
        elif "archaea" in full_lineage:
            kingdom = "Archaea"
        elif "fungi" in full_lineage:
            kingdom = "Fungi"
        elif any(x in full_lineage for x in ["viridiplantae", "plantae", "embryophyta", "tracheophyta"]):
            kingdom = "Plant"
        elif any(x in full_lineage for x in ["metazoa", "animalia", "chordata", "arthropoda", "mammalia"]):
            kingdom = "Animal"
        elif "eukaryota" in full_lineage or "eukaryote" in full_lineage:
            kingdom = "Eukaryota"
        elif any(x in full_lineage for x in ["virus", "viruses", "viral"]):
            kingdom = "Virus"

        # ==================== 2. Phylum 智能提取 ====================
        # 更全面的关键词（你可以继续补充）
        phylum_map = {
            "actinobacteria": "Actinobacteria",
            "proteobacteria": "Proteobacteria",
            "firmicutes": "Firmicutes",
            "bacteroidetes": "Bacteroidetes",
            "cyanobacteria": "Cyanobacteria",
            "ascomycota": "Ascomycota",
            "basidiomycota": "Basidiomycota",
            "arthropoda": "Arthropoda",
            "chordata": "Chordata",
            "magnoliophyta": "Magnoliophyta",
            "euryarchaeota": "Euryarchaeota",
        }

        for part in lineage_parts:
            lower_part = part.lower()
            for key, value in phylum_map.items():
                if key in lower_part:
                    phylum = value
                    break
            if phylum != "Unknown":
                break

        # 如果还没找到，就取第二个作为 phylum（大多数情况下是门）
        if phylum == "Unknown" and len(lineage_parts) >= 2:
            phylum = lineage_parts[1]

        return kingdom, phylum
    
    def _fetch_full_info(self, uid: str):
            """最终干净版：避免重复打印 + 可靠反哺"""
            clean_uid = self._clean_uniprot_id(uid)
            if not clean_uid or len(clean_uid) < 4:
                return uid, {"organism": "Unknown", "kingdom": "Unknown", "phylum_class": "Unknown"}
    
            # ====================== UniProtKB ======================
            try:
                r = requests.get(f"https://rest.uniprot.org/uniprotkb/{clean_uid}.json",
                                 timeout=12, headers={'User-Agent': 'CHR-SSN-Pipeline-V25'})
                if r.status_code == 200:
                    result = r.json()
                    
                    org_data = result.get("organism", {})
                    org_name = (org_data.get("scientificName") or 
                               org_data.get("fullName") or 
                               str(org_data).strip() or "Unknown")
    
                    lineage = result.get("lineage") or []
                    if not lineage and isinstance(org_data, dict):
                        lineage = org_data.get("lineage") or []
                    if not lineage:
                        lineage = result.get("taxonLineage") or []
    
                    kingdom, phylum = self._parse_lineage(lineage)
    
                    # ==================== 反哺全局缓存 ====================
                    if org_name != "Unknown" and kingdom != "Unknown":
                        std_name = self._standardize_organism_name(org_name)
                        cache_key = f"org:{std_name}"
                        self.kingdom_cache[cache_key] = {
                            "kingdom": kingdom,
                            "phylum_class": phylum,
                            "literature": ""
                        }
                        print(f"    ♻ 反哺全局缓存: {org_name} → {kingdom}")
                    # ============================================================
    
                    self.uniprot_org_cache[clean_uid] = org_name
    
                    # 只打印一次关键信息
                    if org_name != "Unknown" or kingdom != "Unknown":
                        print(f"    ✓ {clean_uid} → {org_name} | {kingdom}")
                    return uid, {"organism": org_name, "kingdom": kingdom, "phylum_class": phylum}
            except:
                pass
    
            # ====================== UniParc ======================
            print(f"    🔄 {clean_uid} UniProtKB 无结果，尝试 UniParc...")
            try:
                r = requests.get(f"https://rest.uniprot.org/uniparc/search?query=uniprotkb:{clean_uid}&format=json&size=1",
                                 timeout=15, headers={'User-Agent': 'CHR-SSN-Pipeline-V25'})
                if r.status_code == 200:
                    results = r.json().get("results", [])
                    if results:
                        entry = results[0]
                        org_name = entry.get("organism", {}).get("scientificName") or "Unknown"
    
                        lineage = entry.get("taxonomicLineage", []) or []
                        if not lineage and isinstance(entry.get("organismTaxon"), dict):
                            lineage = entry["organismTaxon"].get("lineage", [])
    
                        kingdom, phylum = self._parse_lineage(lineage)
    
                        if org_name != "Unknown" and kingdom != "Unknown":
                            std_name = self._standardize_organism_name(org_name)
                            cache_key = f"org:{std_name}"
                            self.kingdom_cache[cache_key] = {
                                "kingdom": kingdom,
                                "phylum_class": phylum,
                                "literature": ""
                            }
                            print(f"    ♻ UniParc 反哺全局缓存: {org_name} → {kingdom}")
    
                        self.uniprot_org_cache[clean_uid] = org_name
                        print(f"    ✓ UniParc 找回: {clean_uid} → {org_name} | {kingdom}")
                        return uid, {"organism": org_name, "kingdom": kingdom, "phylum_class": phylum}
            except Exception as e:
                print(f"    ⚠ UniParc 失败: {e}")
    
            # 最终失败
            self.uniprot_org_cache[clean_uid] = "Unknown"
            print(f"    ✗ {clean_uid} → Unknown")
            return uid, {"organism": "Unknown", "kingdom": "Unknown", "phylum_class": "Unknown"}

    
    def _get_local_dict(self, organism_name: str) -> dict:
        """优先查询全局 organism → kingdom 缓存"""
        if not organism_name or organism_name == "Unknown":
            return {"kingdom": "Unknown", "phylum_class": "Unknown", "literature": ""}

        std_name = self._standardize_organism_name(organism_name)
        cache_key = f"org:{std_name}"
        
        if cache_key in self.kingdom_cache:
            cached = self.kingdom_cache[cache_key]
            if isinstance(cached, dict):
                return cached
        
        return {"kingdom": "Unknown", "phylum_class": "Unknown", "literature": ""}

    def _get_kimi_kingdom(self, organism_name: str) -> dict:
        """调用 Kimi 前增加严格防护"""
        print(f"    [DEBUG] _get_kimi_kingdom called with: '{organism_name}'")
        if not organism_name or organism_name in ("Unknown", "", "{}", "None", None):
            return {"kingdom": "Unknown", "phylum_class": "Unknown"}

        local = self._get_local_dict(organism_name)
        if local["kingdom"] != "Unknown":
            return local

        # ... 后面保持你原来的 Kimi 调用代码 ...
        for attempt in range(3):
            try:
                prompt = f"""你是一个专业的生物分类专家。
物种名称: {organism_name}
请严格返回JSON（不要任何其他文字）：
{{
  "kingdom": "Plant" | "Animal" | "Fungi" | "Bacteria" | "Archaea" | "Protist" | "Unknown",
  "phylum_class": "具体门/纲名称"
}}"""
                response = requests.post(
                    "https://api.moonshot.cn/v1/chat/completions",
                    headers={"Authorization": f"Bearer {Config.KIMI_API_KEY}", "Content-Type": "application/json"},
                    json={"model": Config.KIMI_MODEL, "messages": [{"role": "user", "content": prompt}], "temperature": 1, "max_tokens": 200},
                    timeout=30
                )
                if response.status_code == 429:
                    time.sleep(2 ** attempt + 2)
                    continue
                if response.status_code != 200:
                    return {"kingdom": "Unknown", "phylum_class": "Unknown"}

                content = response.json()["choices"][0]["message"]["content"].strip()
                if content.startswith("```"):
                    content = content.split("```")[1].replace("json", "").strip()
                info = json.loads(content)

                result = {
                    "kingdom": info.get("kingdom", "Unknown"),
                    "phylum_class": info.get("phylum_class", "Unknown")
                }
                std_name = self._standardize_organism_name(organism_name)
                self.kingdom_cache[f"org:{std_name}"] = result
                self.save_caches()
                print(f" ✓ Kimi 补全: {organism_name} → {result['kingdom']}")
                return result
            except Exception as e:
                print(f" ⚠ Kimi 失败 {attempt+1}/3: {e}")
                if attempt < 2:
                    time.sleep(2)
        return {"kingdom": "Unknown", "phylum_class": "Unknown"} 

    def get_full_info(self, uniprot_id: str) -> dict:
        clean_id = self._clean_uniprot_id(uniprot_id)
        organism = self.uniprot_org_cache.get(clean_id, "Unknown")
        if organism == "Unknown":
            return {"organism": "Unknown", "kingdom": "Unknown", "phylum_class": "Unknown"}
        local = self._get_local_dict(organism)
        return {
            "organism": organism,
            "kingdom": local.get("kingdom", "Unknown"),
            "phylum_class": local.get("phylum_class", "Unknown")
        }

    def reconnaissance_kingdoms(self, unique_uniprot_ids: List[str]):
        print("\n🔍 Reconnaissance Phase：UniProt ID → organism（永久缓存）...")
        to_query = []
        for uid in unique_uniprot_ids:
            if not uid: continue
            clean_uid = self._clean_uniprot_id(uid)
            if clean_uid in self.uniprot_org_cache and self.uniprot_org_cache[clean_uid] != "Unknown":
                continue
            to_query.append(uid)
        to_query = to_query[:self.max_recon]
        if not to_query:
            print("  ✓ 所有 UniProt ID 已有 organism 缓存")
            return
        print(f"  需要查询 {len(to_query)} 个新 UniProt ID...")
        pbar = tqdm(total=len(to_query), desc="UniProt 查询") if TQDM_AVAILABLE else None
        with ThreadPoolExecutor(max_workers=2) as executor:
            future_to_uid = {executor.submit(self._fetch_full_info, uid): uid for uid in to_query}
            for future in as_completed(future_to_uid):
                raw_uid, info = future.result()
                clean_uid = self._clean_uniprot_id(raw_uid)
                organism = info.get("organism", "Unknown")
                self.uniprot_org_cache[clean_uid] = organism
                if pbar:
                    pbar.update(1)
                if organism != "Unknown":
                    print(f"    ✓ {clean_uid} → {organism}")
                time.sleep(0.12)
        if pbar: pbar.close()
        self.save_caches()
        print(f"  ✓ ID→organism 缓存更新完成（当前共 {len(self.uniprot_org_cache)} 条）")

    # ====================== parse_fasta ======================
    def _load_core_ids(self) -> Set[str]:
        if not self.core_fasta or not self.core_fasta.exists():
            print("  ℹ 未提供核心序列文件，所有序列平等采样")
            return set()
        core_ids = set()
        with open(self.core_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    seq_id = line[1:].split()[0]
                    core_ids.add(seq_id)
        print(f"  ✓ 加载核心序列 {len(core_ids)} 条")
        return core_ids

    def _create_record(self, header: str, sequence: str) -> SequenceRecord:
        raw_id = header.split()[0] if header else ""
        uniprot_id = self._clean_uniprot_id(raw_id)
        species_code = "UNKNOWN"
        if '_' in raw_id and not raw_id.startswith('UniRef'):
            parts = raw_id.rsplit('_', 1)
            if len(parts) == 2 and len(parts[0]) >= 4:
                species_code = parts[1]
        is_core = raw_id in self.core_ids
        if is_core:
            print(f"    ⭐ 核心序列标记: {raw_id}")
        return SequenceRecord(id=header, uniprot_id=uniprot_id, species_code=species_code, sequence=sequence, is_core=is_core)

    def parse_fasta(self) -> List[SequenceRecord]:
        print(f"\n📥 解析 FASTA: {self.input_fasta}")
        self.core_ids = self._load_core_ids()
        records = []
        uniprot_ids = set()
        current_id = None
        current_seq = []
        with open(self.input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                    if current_id and current_seq:
                        rec = self._create_record(current_id, ''.join(current_seq))
                        uniprot_ids.add(rec.uniprot_id)
                        if (Config.LENGTH_MIN <= rec.length <= Config.LENGTH_MAX) or rec.is_core:
                            records.append(rec)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id and current_seq:
            rec = self._create_record(current_id, ''.join(current_seq))
            uniprot_ids.add(rec.uniprot_id)
            if (Config.LENGTH_MIN <= rec.length <= Config.LENGTH_MAX) or rec.is_core:
                records.append(rec)

        self.reconnaissance_kingdoms(list(uniprot_ids))

        for rec in records:
            info = self.get_full_info(rec.uniprot_id)
            rec.organism = info["organism"]
            rec.kingdom = info["kingdom"]
            rec.phylum_class = info["phylum_class"]

        print("\n🤖 Kimi + 全局字典 补全 Kingdom...")
        organism_to_records = defaultdict(list)
        for rec in records:
            if rec.organism != "Unknown" and rec.kingdom in ("Unknown", "Eukaryota"):
                organism_to_records[rec.organism].append(rec)

        to_call = [org for org in organism_to_records if self._get_local_dict(org)["kingdom"] == "Unknown"]
        print(f" 需要调用 Kimi 的 organism: {len(to_call)} 个")

        if to_call:
            kimi_count = 0
            with ThreadPoolExecutor(max_workers=2) as executor:
                future_to_org = {executor.submit(self._get_kimi_kingdom, org): org for org in to_call}
                for future in tqdm(as_completed(future_to_org), total=len(to_call), desc="Kimi 补全"):
                    org = future_to_org[future]
                    try:
                        info = future.result()
                        if info["kingdom"] != "Unknown":
                            for rec in organism_to_records[org]:
                                rec.kingdom = info["kingdom"]
                                rec.phylum_class = info["phylum_class"]
                            kimi_count += 1
                    except Exception as e:
                        print(f" ⚠ Kimi 异常: {org} | {e}")
                    time.sleep(1.0)
            print(f" ✓ Kimi 补全完成！成功 {kimi_count} 个")

        records = self.filter_valid_sequences(records)
        print(f"  🔍 当前 Unknown 数量: {sum(1 for r in records if r.kingdom == 'Unknown')}")
        return records
 
        
    def filter_valid_sequences(self, records: List[SequenceRecord]) -> List[SequenceRecord]:
        print(f"  ✓ 所有序列全部保留（共 {len(records)} 条）")
        return records

    # ====================== 聚类策略 ======================
    def cluster_sequences(self) -> List[SequenceRecord]:
        print(f"\n🔄 聚类降维 (target identity={self.cdhit_identity})...")
        
        # Tier 1: MMseqs2
        reps = self._cluster_mmseqs2()
        if reps:
            return self._finalize_representatives(reps, tool="MMseqs2")
        
        # Tier 2: CD-HIT
        if self.cdhit_identity >= 0.50 and self.config.CDHIT:
            reps = self._cluster_cdhit()
            if reps:
                return self._finalize_representatives(reps, tool="CD-HIT")
        
        # Tier 3: 长度分箱保底
        print("  → 回退到长度分箱 + 核心保留...")
        reps = self._cluster_by_length_diversity()
        return self._finalize_representatives(reps, tool="Length-Diversity")

    def _cluster_mmseqs2(self) -> Optional[List[SequenceRecord]]:
        mmseqs = shutil.which("mmseqs")
        if not mmseqs:
            print("  ℹ MMseqs2 未安装，跳过")
            return None
        
        effective_identity = max(self.cdhit_identity, 0.30)
        filtered_fasta = self.output_dir / "filtered_for_cluster.fasta"
        with open(filtered_fasta, 'w') as f:
            for rec in self.sequences:
                f.write(f">{rec.id}\n{rec.sequence}\n")
        
        output_base = self.output_dir / "mmseqs_cluster"
        tmp_dir = self.output_dir / "mmseqs_tmp"
        clu_tsv = self.output_dir / "mmseqs_cluster_cluster.tsv"
        
        cmd = [
            mmseqs, "easy-linclust", str(filtered_fasta), str(output_base), str(tmp_dir),
            "--min-seq-id", str(effective_identity),
            "--cov-mode", "1", "--threads", "4", "--remove-tmp-files", "1"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode != 0:
                print(f"  ⚠ MMseqs2 失败: {result.stderr[:500]}")
                return None
        except Exception as e:
            print(f"  ⚠ MMseqs2 异常: {e}")
            return None
        
        if not clu_tsv.exists():
            alt = self.output_dir / "mmseqs_cluster.tsv"
            if alt.exists():
                clu_tsv = alt
            else:
                return None
        
        clusters = defaultdict(list)
        with open(clu_tsv, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    clusters[parts[0]].append(parts[1])
        
        print(f"  ✓ MMseqs2 解析: {len(clusters)} 个 cluster")
        return self._select_representatives_from_clusters(clusters, "MMseqs2")

    def _cluster_cdhit(self) -> Optional[List[SequenceRecord]]:
        if not self.config.CDHIT:
            return None
        # ... (如果你之前有完整 _cluster_cdhit 可以保留，这里用简化版)
        print("  CD-HIT 暂未启用，回退 MMseqs2")
        return None

    def _cluster_by_length_diversity(self) -> List[SequenceRecord]:
        core_seqs = [r for r in self.sequences if r.is_core]
        non_core = [r for r in self.sequences if not r.is_core]
        
        bins = defaultdict(list)
        for rec in non_core:
            bins[rec.length // 50].append(rec)
        
        representatives = []
        for rec in core_seqs:
            rec.is_representative = True
            rec.cluster_id = "CORE"
            representatives.append(rec)
        
        remaining = self.max_rep - len(representatives)
        count = 0
        for bin_key in sorted(bins.keys()):
            if count >= remaining:
                break
            bin_seqs = sorted(bins[bin_key], key=lambda x: x.length, reverse=True)
            for rec in bin_seqs[:2]:
                if count >= remaining:
                    break
                rec.is_representative = True
                rec.cluster_id = f"LEN_{bin_key}"
                representatives.append(rec)
                count += 1
        print(f"  ✓ 长度分箱采样: {len(representatives)} 条")
        return representatives
    def _select_representatives_from_clusters(self, clusters: Dict[str, List[str]], tool: str) -> List[SequenceRecord]:
        """动态采样：每个 cluster 最多选 2 条 + 核心强制保留"""
        id_to_record = {rec.id: rec for rec in self.sequences}
        representatives = []
        core_kept = 0

        # 核心序列全部保留
        for rec in self.sequences:
            if rec.is_core:
                rec.is_representative = True
                rec.cluster_id = "CORE"
                representatives.append(rec)
                core_kept += 1

        remaining_quota = self.max_rep - len(representatives)
        if remaining_quota <= 0:
            return representatives[:self.max_rep]

        print(f"  📊 动态采样: {len(clusters)} 个 cluster，剩余配额 {remaining_quota}，每个 cluster 最多选 2 条...")

        sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)

        for rep_id, members in sorted_clusters:
            cluster_recs = [id_to_record[mid] for mid in members if mid in id_to_record]
            if not cluster_recs:
                continue
            
            cluster_recs.sort(key=lambda x: x.length, reverse=True)
            
            to_take = min(2, remaining_quota, len(cluster_recs))   # ← 改成 2
            
            for rec in cluster_recs[:to_take]:
                if len(representatives) >= self.max_rep:
                    break
                if not rec.is_representative:
                    rec.is_representative = True
                    rec.cluster_id = f"{tool}_{rep_id[:20]}"
                    representatives.append(rec)
                    remaining_quota -= 1

            if len(representatives) >= self.max_rep:
                break

        print(f"  ✓ {tool} 最终选择: {len(representatives)} 条代表序列（核心 {core_kept} 条）")
        return representatives
        
    def _finalize_representatives(self, representatives: List[SequenceRecord], tool: str) -> List[SequenceRecord]:
        """后处理：同步信息 + 二次 Kimi 补全（增加防护）"""
        print(f" 📌 {tool} 后处理：同步 Kingdom/Organism/Phylum...")
        
        for rep in representatives:
            info = self.get_full_info(rep.uniprot_id)
            rep.organism = info.get("organism", "Unknown")
            rep.phylum_class = info.get("phylum_class", "Unknown")
            
            # 关键防护：只有 organism 有效时才调用 Kimi
            if (rep.kingdom == "Unknown" and 
                rep.organism not in ("Unknown", "", "{}", "None", None)):
                kimi_info = self._get_kimi_kingdom(rep.organism)
                rep.kingdom = kimi_info["kingdom"]
                rep.phylum_class = kimi_info["phylum_class"]
        
        # 统计
        kingdom_stats = defaultdict(int)
        for rep in representatives:
            kingdom_stats[rep.kingdom] += 1
        print(f"  📊 代表序列 Kingdom 分布: {dict(kingdom_stats)}")
        
        core_in_reps = sum(1 for r in representatives if r.is_core)
        print(f"  ⭐ 核心序列保留: {core_in_reps}/{len(self.core_ids)} 条")
        return representatives

    def final_kimi_rescue_for_unknown(self):
        """最后一层 Kimi 补救（修复 400 错误）"""
        unknown_reps = [rec for rec in self.representative_seqs if rec.kingdom == "Unknown"]
        if not unknown_reps:
            print("  ✓ 无 Unknown 序列需要补救")
            return
    
        # === 防御：如果 API Key 是空/假值，直接跳过，避免无意义请求 ===
        if not Config.KIMI_API_KEY or str(Config.KIMI_API_KEY).strip().lower() in ("none", "", "your-api-key"):
            print("  ⚠ Config.KIMI_API_KEY 未设置，跳过 Kimi 补救。请检查环境变量或重新生成 Key。")
            return
    
        print(f"\n🔥 最后一层 Kimi 补救：对 {len(unknown_reps)} 个 Unknown 序列进行强力补救...")
    
        rescued = 0
        # === 关键修复：URL 末尾绝对不能有空格 ===
        url = "https://api.moonshot.cn/v1/chat/completions"
    
        with tqdm(unknown_reps[:200], desc="Kimi 补救", unit="id") as pbar:
            for rec in pbar:
                if not rec.uniprot_id or rec.uniprot_id in ("Unknown", "", None):
                    continue
    
                try:
                    prompt = f"""你是一个专业的生物信息学专家。
    UniProt ID: {rec.uniprot_id}
    请查询该蛋白的物种科学名称和生物分类 Kingdom。
    请严格返回以下 JSON 格式，不要任何其他文字：
    {{"organism": "物种科学名", "kingdom": "Bacteria|Fungi|Plant|Animal|Archaea|Protist|Unknown", "phylum_class": "门/纲名称"}}"""
    
                    response = requests.post(
                        url,
                        headers={
                            "Authorization": f"Bearer {Config.KIMI_API_KEY}",
                            "Content-Type": "application/json"
                        },
                        json={
                            "model": Config.KIMI_MODEL,
                            "messages": [{"role": "user", "content": prompt}],
                            "temperature": 0.0,
                            "max_tokens": 300
                        },
                        timeout=30
                    )
    
                    # === 调试：一旦遇到 400，立刻打印原始返回体，方便定位 ===
                    if response.status_code == 400:
                        err_msg = response.text[:300].replace('\n', ' ')
                        print(f"\n    ⚠ 400 错误详情: {err_msg}")
                        pbar.set_postfix_str("400 Bad Request")
                        continue
                    if response.status_code == 401:
                        print(f"\n    ⚠ 401 Unauthorized: API Key 已被吊销或无效，请重新生成。")
                        pbar.set_postfix_str("401 Auth Error")
                        return  # Key 无效，继续跑也没意义，直接退出
                    if response.status_code == 429:
                        pbar.set_postfix_str("429 Rate Limit")
                        time.sleep(5)
                        continue
                    if response.status_code != 200:
                        pbar.set_postfix_str(f"HTTP {response.status_code}")
                        continue
    
                    content = response.json()["choices"][0]["message"]["content"].strip()
                    if content.startswith("```"):
                        content = content.split("```")[1].replace("json", "").strip()
    
                    info = json.loads(content)
    
                    new_kingdom = str(info.get("kingdom", "")).strip()
                    new_organism = str(info.get("organism", "")).strip()
    
                    if new_kingdom and new_kingdom.lower() not in ["unknown", ""]:
                        rec.organism = new_organism or rec.organism
                        rec.kingdom = new_kingdom
                        rec.phylum_class = str(info.get("phylum_class", "")).strip()
                        rescued += 1
    
                        clean_id = self._clean_uniprot_id(rec.uniprot_id)
                        self.uniprot_org_cache[clean_id] = rec.organism
                        std_name = self._standardize_organism_name(rec.organism)
                        self.kingdom_cache[f"org:{std_name}"] = {
                            "kingdom": rec.kingdom,
                            "phylum_class": rec.phylum_class,
                            "literature": ""
                        }
    
                        pbar.set_postfix_str(f"成功 {rescued}")
                        print(f"    ✅ 补救: {rec.uniprot_id} → {rec.organism} | {rec.kingdom}")
    
                    time.sleep(1.1)
    
                except Exception as e:
                    pbar.set_postfix_str(f"异常: {str(e)[:30]}")
        for rec in self.representative_seqs:
            if rec.id in self.nodes:
                self.nodes[rec.id]['kingdom'] = rec.kingdom
                self.nodes[rec.id]['organism'] = rec.organism
                self.nodes[rec.id]['phylum_class'] = rec.phylum_class
        self.save_caches()
        print(f"  🔥 Kimi 补救完成！成功补救 {rescued} 个 Unknown 序列")
    def run_msa(self, rep_fasta: Path, aln_file: Path) -> bool:
        if not Path(self.config.CLUSTALO).exists():
            return False
        cmd = [self.config.CLUSTALO, "-i", str(rep_fasta), "-o", str(aln_file),
               "--outfmt", "fasta", "--force", "--threads", str(self.config.MSA_THREADS)]
        result = subprocess.run(cmd, capture_output=True)
        return result.returncode == 0

    def _parse_alignment(self, aln_file: Path) -> Dict[str, str]:
        sequences = {}
        current_id = None
        current_seq = []
        with open(aln_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq)
        return sequences

    def _calculate_aligned_identity(self, seq1: str, seq2: str) -> float:
        pairs = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
        if not pairs:
            return 0.0
        matches = sum(1 for a, b in pairs if a == b)
        return matches / len(pairs)

    def build_similarity_network(self, aln_file: Path) -> List[SSNEdge]:
        """基于 MMseqs2 all-vs-all 构建 SSN，彻底替代 MSA + pairwise identity"""
        print(f"\n🕸 构建序列相似性网络 (MMseqs2 all-vs-all)...")
        
        # aln_file 参数在此新版中实际代表 rep_fasta（保持函数签名兼容）
        rep_fasta = aln_file
        if not rep_fasta.exists():
            raise FileNotFoundError(f"代表序列文件不存在: {rep_fasta}")
        
        result_m8 = self.output_dir / "mmseqs_avA.m8"
        tmp_dir = self.output_dir / "mmseqs_tmp_avA"
        
        mmseqs = shutil.which("mmseqs")
        if not mmseqs:
            raise RuntimeError("MMseqs2 未安装，无法构建 SSN")
        
        cmd = [
            mmseqs, "easy-search",
            str(rep_fasta), str(rep_fasta), str(result_m8), str(tmp_dir),
            "--min-seq-id", str(self.config.SSN_IDENTITY_THRESHOLD),
            "-c", "0.7",               # alignment coverage ≥70%
            "--cov-mode", "1",
            "--format-output", "query,target,pident,alnlen,evalue,bits",
            "--threads", "8",
            "--remove-tmp-files", "1"
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        except subprocess.CalledProcessError as e:
            print(f"  ⚠ MMseqs2 失败: {e.stderr[:500]}")
            return []
        
        if not result_m8.exists():
            return []
        
        # 1. 读取所有唯一边（去重 + 保留最高 pident）
        edge_by_key: Dict[tuple, SSNEdge] = {}
        with open(result_m8, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                q, t, pident_str = parts[0], parts[1], parts[2]
                if q == t:
                    continue
                key = tuple(sorted([q, t]))
                pident = float(pident_str) / 100.0
                if key not in edge_by_key or pident > edge_by_key[key].identity:
                    edge_by_key[key] = SSNEdge(q, t, pident)
        
        # 2. 构建节点邻接表
        node_edges: Dict[str, List[SSNEdge]] = defaultdict(list)
        for e in edge_by_key.values():
            node_edges[e.source].append(e)
            node_edges[e.target].append(e)
        
        # 3. 每个节点只保留 TOP_K_EDGES 条最强边（EFI-EST 标准做法）
        kept_keys: Set[tuple] = set()
        for node, e_list in node_edges.items():
            e_list.sort(key=lambda x: x.identity, reverse=True)
            for e in e_list[:self.config.TOP_K_EDGES]:
                kept_keys.add(tuple(sorted([e.source, e.target])))
        
        final_edges = [edge_by_key[k] for k in kept_keys]
        print(f"  ✅ 共生成 {len(final_edges)} 条边（MMseqs2 all-vs-all，阈值 {self.config.SSN_IDENTITY_THRESHOLD}）")
        return final_edges

    def build_nodes(self):
        degree_count = defaultdict(int)
        for edge in self.edges:
            degree_count[edge.source] += 1
            degree_count[edge.target] += 1
        for rec in self.representative_seqs:
            self.nodes[rec.id] = {
                'id': rec.id,
                'uniprot_id': rec.uniprot_id,
                'species_code': rec.species_code,
                'kingdom': rec.kingdom,
                'organism': rec.organism,
                'phylum_class': rec.phylum_class,
                'length': rec.length,
                'cluster_id': rec.cluster_id,
                'degree': degree_count.get(rec.id, 0),
                'is_core': rec.is_core,
            }

    def detect_communities(self):
        print("\n🔍 Louvain 社区检测...")
        G = nx.Graph()
        for rec in self.representative_seqs:
            G.add_node(rec.id)
        for e in self.edges:
            G.add_edge(e.source, e.target)
        communities = nx.community.louvain_communities(G)
        community_map = {}
        for idx, comm in enumerate(communities):
            for node in comm:
                community_map[node] = idx
        for rec in self.representative_seqs:
            rec.community = community_map.get(rec.id, -1)
        print(f"  ✓ 检测到 {len(communities)} 个社区")

    def export_cytoscape_files(self):
        self.detect_communities()
        print(f"\n💾 导出 Cytoscape 文件...")
        with open(self.output_dir / "ssn_edges.csv", 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["source", "target", "identity", "interaction"])
            for e in self.edges:
                writer.writerow([e.source, e.target, f"{e.identity:.4f}", "pp"])
        headers = ["id", "uniprot_id", "species_code", "kingdom", "organism", "phylum_class",
                   "length", "cluster_id", "degree", "community", "is_core"]
        with open(self.output_dir / "ssn_nodes.csv", 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            for n in self.nodes.values():
                row = {**n, "community": next((r.community for r in self.representative_seqs if r.id == n['id']), -1)}
                writer.writerow(row)
        self._export_xgmml()
        print(f"  ✓ 导出完成")

    def _export_xgmml(self):
        file_path = self.output_dir / "ssn_network.xgmml"
        colors = {'Plant': '#228B22', 'Animal': '#4169E1', 'Fungi': '#FF8C00',
                  'Bacteria': '#DC143C', 'Archaea': '#9932CC', 'Protist': '#FFD700',
                  'Unknown': '#808080', 'Metagenome': '#9E9E9E'}
        with open(file_path, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n<graph label="CHR SSN" xmlns="http://www.cs.rpi.edu/XGMML">\n')
            for node in self.nodes.values():
                base_color = colors.get(node.get('kingdom', 'Unknown'), '#808080')
                is_core = node.get('is_core', False)
                border_width = 3.5 if is_core else 1.5
                node_size = 36 if is_core else 32
                safe_id = self.sanitize_id(node['id'])
                label = node["uniprot_id"]
                if is_core:
                    label = f"★{label}"
                f.write(f'  <node id="{safe_id}" label="{label}">\n')
                f.write(f'    <att name="kingdom" value="{node.get("kingdom", "Unknown")}"/>\n')
                f.write(f'    <att name="organism" value="{node.get("organism", "Unknown")}"/>\n')
                f.write(f'    <att name="is_core" value="{str(is_core).lower()}"/>\n')
                f.write(f'    <graphics fill="{base_color}" width="{border_width}" size="{node_size}" outline="#000000" labelSize="11" />\n')
                f.write('  </node>\n')
            for edge in self.edges:
                f.write(f'  <edge source="{self.sanitize_id(edge.source)}" target="{self.sanitize_id(edge.target)}">\n')
                f.write(f'    <att name="identity" value="{edge.identity:.4f}"/>\n')
                f.write('  </edge>\n')
            f.write('</graph>\n')
        print(f"  ✓ XGMML 已生成: {file_path}")

    def run_pipeline(self):
        print(f"\n{'#'*80}\n# CHR SSN Pipeline V25 - 双永久缓存版\n{'#'*80}\n")
        self.sequences = self.parse_fasta()
        self.representative_seqs = self.cluster_sequences()

        rep_fasta = self.output_dir / "representatives.fasta"
        with open(rep_fasta, 'w') as f:
            for rec in self.representative_seqs:
                f.write(f">{rec.id}\n{rec.sequence}\n")

        aln_file = self.output_dir / "representatives_aln.fasta"
        msa_ok = self.run_msa(rep_fasta, aln_file)
        if not msa_ok:
            print("  ⚠ ClustalO MSA 失败（SSN 不再依赖 MSA，继续运行）")

        self.edges = self.build_similarity_network(rep_fasta)
        
        # 先完成所有 taxonomy 补救，再构建 nodes，确保导出数据最新
        self.final_kimi_rescue_for_unknown()
        self.build_nodes()
        
        self.export_cytoscape_files()
        print(f"\n{'='*80}\n✅ V25 双缓存版运行完成！\n{'='*80}")


def main():
    parser = argparse.ArgumentParser(description='SSN Pipeline V25 - 双永久缓存版')
    parser.add_argument('input', help='输入 FASTA 文件')
    parser.add_argument('--core-fasta', default=None)
    parser.add_argument('--output', default=None)
    parser.add_argument('--max-rep', type=int, default=2500)
    parser.add_argument('--cdhit-identity', type=float, default=0.35)
    parser.add_argument('--cdhit-coverage', type=float, default=0.75)
    parser.add_argument('--max-recon', type=int, default=4000)

    args = parser.parse_args()

    pipeline = CHRSSNPipelineCDHIT(
        input_fasta=Path(args.input),
        output_dir=Path(args.output) if args.output else None,
        max_rep=args.max_rep,
        cdhit_identity=args.cdhit_identity,
        cdhit_coverage=args.cdhit_coverage,
        max_recon=args.max_recon,
        core_fasta=Path(args.core_fasta) if args.core_fasta else None
    )
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()
