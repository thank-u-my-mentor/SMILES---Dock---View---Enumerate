#!/usr/bin/env python3
"""
Zotero Photoenzymatic Radical Knowledge Agent
==============================================
从Zotero导出的CSV出发，构建光酶/自由基文献的知识脉络。

Usage:
    1. Zotero中选中你的20-30篇文献 → 右键 "导出条目" → 格式选CSV → 保存为 zotero_library.csv
    2. python zotero_photoenzyme_agent.py --input zotero_library.csv --output ./knowledge_base/
"""

from __future__ import annotations

import os
import sys
import json
import argparse
import re
import time
import html as html_lib
import sqlite3
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Set, Tuple
from collections import defaultdict, Counter
import hashlib

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from matplotlib.patches import Rectangle
import requests

try:
    from pyvis.network import Network
except ImportError:
    Network = None

# Discovery 模块（可选导入，失败不阻塞原有流程）
try:
    from discovery import (
        SemanticScholarClient, LiteratureDiscovery, CitationNetworkAnalyzer,
        DiscoveryReportGenerator, DiscoveredPaper
    )
    _DISCOVERY_AVAILABLE = True
except ImportError as _e:
    _DISCOVERY_AVAILABLE = False
    _DISCOVERY_IMPORT_ERROR = str(_e)


# ==================== 配置区 ====================
LLM_CONFIG = {
    "api_key": os.getenv("KIMI_API_KEY") or os.getenv("MOONSHOT_API_KEY") or "",
    "base_url": os.getenv("MOONSHOT_BASE_URL") or "https://api.moonshot.cn/v1",
    "model": "kimi-k2.6",  
    "temperature": 1,
    "max_tokens": 3600,
    "connect_timeout": 15,
    "read_timeout": 360,
    "json_mode": True,
    "enable_thinking": False
}

# 光酶自由基领域的结构化提取维度
EXTRACTION_SCHEMA = {
    "enzyme_category": "酶催化类型: [天然辅因子依赖型 | 协同光催化/酶催化 | 人工光酶(合成辅因子) | 无辅因子光酶]",
    "enzyme_family": "具体酶家族: [P450 | ERED/OYE | FAP(脂肪酸光脱羧酶) | Aldolase | 脂肪酶/酯酶 | 其他]",
    "radical_mechanism": "自由基产生机制: [光氧化还原(PET) | 能量转移(EnT) | 电子转移(ET) | 直接激发辅因子 | 底物直接激发]",
    "reaction_type": "反应类型: [C-H芳基化 | C-H烷基化 | 脱羧偶联 | 卤化 | 环加成 | Giese加成 | 氢原子转移(HAT) | 其他]",
    "bond_formation": "形成的键: [C-C | C-N | C-O | C-S | C-Se | C-B | C-Si]",
    "substrate_class": "底物类别: [羧酸 | 烯烃 | 芳烃 | 亚胺 | 羰基化合物 | 卤代物 | 磺酰胺]",
    "stereochemistry": "立体化学: [对映选择性 | 非对映选择性 | 无立体控制 | 动力学拆分]",
    "has_mechanistic_study": "是否含机制研究: [是 | 否]",
    "mechanistic_methods": "机制研究手段: [晶体结构 | EPR/ESR | 瞬态吸收光谱 | DFT计算 | 分子动力学 | 同位素标记 | 无]",
    "light_source": "光源: [可见光(蓝光/绿光) | UV | 太阳光 | LED | 未明确]",
    "cofactor": "辅因子/光敏剂: [FAD/FMN | NAD(P)H | Ru(bpy)3 | Ir(ppy)3 | 有机染料 | 无]",
    "key_innovation": "核心创新点(一句话)",
    "limitation": "局限性/挑战(一句话)",
    "future_direction": "未来方向暗示(一句话)"
}

SCHEMA_PROFILES = {
    "photoenzyme": EXTRACTION_SCHEMA,
    "photo-radical": {
        "topic_domain": "主题领域: [photoenzymatic radical | photoenzymatic catalysis | photoredox radical reaction | radical chemical reaction | enzymatic catalysis | organic photochemistry | other]",
        "paper_domain": "文献域判断: [photoenzymatic/enzymatic | enzymatic non-photo | organic photoredox/radical no-enzyme | organic synthesis no-enzyme | unrelated/blacklist]",
        "photoenzymatic_status": "是否已接入酶体系: [已接入photoenzymatic | 仅光化学/自由基反应 | 仅酶催化背景 | 可设想接入酶体系 | 未明确]",
        "enzyme_family": "涉及的酶家族: [P450 | ERED/OYE | FAP | flavoprotein | nonheme iron enzyme | lipase/esterase | alpha/beta-hydrolase | engineered enzyme | no enzyme | other | unclear]",
        "radical_or_excited_state_mechanism": "自由基/激发态机制: [PET/photoredox | EnT | HAT | EDA complex | LMCT | SET | triplet diradical | radical-polar crossover | direct substrate excitation | cofactor excitation | unclear]",
        "radical_source": "自由基来源: [carboxylic acid/decarboxylation | alkyl halide | sulfonamide | pyridinium salt | diazo/azo precursor | peroxide | amine | alcohol/aldehyde | alkene | aromatic/heteroaromatic substrate | enzyme cofactor | unclear]",
        "reaction_type": "反应类型: [C-H functionalization | arylation | alkylation | cycloaddition | rearrangement/migration | Minisci reaction | Giese addition | decarboxylative coupling | halogenation | hydrofunctionalization | radical-polar crossover | other]",
        "bond_formation": "形成或断裂的关键键: [C-C | C-N | C-O | C-S | C-X | C-H | C-B | C-Si | skeletal rearrangement | unclear]",
        "substrate_class": "底物类别: [alkene | pyridine/heteroarene | arene | carboxylic acid | sulfonamide | amine | carbonyl compound | imine | alkyl halide | drug-like molecule | amino acid derivative | lipid/ester | other]",
        "product_or_scaffold": "产物类别或目标骨架(短语)",
        "selectivity": "选择性: [enantioselective | diastereoselective | regioselective | chemoselective | site-selective | stereocontrol unclear | no selectivity focus]",
        "mechanistic_evidence": "机制证据: [EPR/ESR | transient absorption | Stern-Volmer | quantum yield | radical clock | isotope labeling | DFT | crystal structure | MD | control experiments | none mentioned]",
        "light_source": "光源: [blue LED | green LED | visible light | UV | solar | no light | unclear]",
        "catalyst_or_photosensitizer": "催化剂/光敏剂/辅因子: [flavin/FAD/FMN | NAD(P)H | Ru | Ir | organic photocatalyst | quinolinone | palladium | enzyme active site | no catalyst | unclear]",
        "enzyme_integration_opportunity": "人工判断辅助: 该反应是否值得考虑接入酶/蛋白口袋: [high | medium | low | unclear]",
        "lipase_or_hydrolase_relevance": "与脂肪酶或alpha/beta水解酶的关系: [direct | pocket/substrate-transfer analogy | stereocontrol analogy | hydrophobic-substrate analogy | no obvious relevance | unclear]",
        "manual_review_priority": "人工复核优先级: [high | medium | low]",
        "metadata_keywords": "来自API/数据库的关键词、学科、概念或期刊主题，多个用 | 分隔",
        "is_blacklisted": "是否应从最终注释集中剔除: [true | false]",
        "blacklist_reason": "若is_blacklisted=true，说明原因；尤其是无DOI、非化学/非酶/不可注释、Other_enzyme无法澄清时",
        "key_innovation": "核心创新点(中文一句话)",
        "limitation": "局限性/挑战(中文一句话)",
        "future_direction": "未来方向或可迁移启发(中文一句话)",
        "confidence": "提取置信度: 0.0-1.0"
    },
    "flavin-photoenzyme": {
        "search_track": "任务分轨: [A_used_flavin_photoenzyme_new_to_nature | B_characterized_flavin_enzyme_not_photoenzymatic | C_photoenzyme_non_flavin | D_organic_photochemistry_no_enzyme | E_background_or_unrelated]",
        "paper_domain": "文献域判断: [flavin photoenzymatic application | characterized flavin enzyme non-photo | photoenzymatic non-flavin | organic photoredox/radical no-enzyme | unrelated/blacklist]",
        "photoenzymatic_application_status": "是否已用于photoenzymatic/new-to-nature反应: [applied to new-to-nature reaction | native/known enzymology only | candidate not yet photoenzymatic | no enzyme | unclear]",
        "new_to_nature_reaction": "是否报道new-to-nature酶反应: [true | false | unclear]",
        "flavin_cofactor": "黄素辅因子类型: [FAD | FMN | flavin unspecified | no flavin | unclear]",
        "enzyme_family": "具体flavin酶家族: [ERED/OYE | ene-reductase | old yellow enzyme | fatty acid photodecarboxylase/FAP | flavoenzyme oxidase | flavin monooxygenase | flavin reductase/dehydrogenase | flavin-dependent halogenase | photolyase/cryptochrome | nitroreductase | berberine bridge enzyme | no enzyme | other | unclear]",
        "enzyme_name_or_target": "文中出现的具体酶名、突变体、UniProt/PDB或候选靶点；没有则填unclear",
        "reaction_type": "反应类型: [enantioselective radical cyclization | radical hydroalkylation | dehalogenation | decarboxylation | C-C bond formation | C-N bond formation | C-H functionalization | reduction/oxidation native | halogenation | other | unclear]",
        "bond_formation": "形成或断裂的关键键: [C-C | C-N | C-O | C-S | C-X | C-H | redox only | unclear]",
        "substrate_class": "底物类别: [alkene | imine | enoate/enone | carboxylic acid/fatty acid | alkyl halide | heteroarene | amino acid derivative | drug-like molecule | native substrate | other | unclear]",
        "selectivity": "选择性: [enantioselective | diastereoselective | regioselective | chemoselective | site-selective | native selectivity only | unclear]",
        "mechanistic_evidence": "机制/表征证据: [crystal structure | PDB/structure | steady-state kinetics | substrate scope | EPR/ESR | transient absorption | Stern-Volmer | radical clock | isotope labeling | mutagenesis | none mentioned]",
        "light_source": "光源: [blue LED | visible light | UV | no light | unclear]",
        "characterized_enzyme_evidence": "酶表征强度: [purified enzyme kinetics | activity assay only | structure only | genome annotation only | no enzyme evidence | unclear]",
        "candidate_for_photoenzyme_repurposing": "是否值得作为尚未被photoenzymatic应用的候选flavin酶: [high | medium | low | no]",
        "structure_similarity_hint": "结构相似或可迁移依据: [same fold/cofactor pocket | homolog of applied enzyme | similar substrate pocket | PDB available | only sequence annotation | unclear]",
        "manual_review_priority": "人工复核优先级: [high | medium | low]",
        "metadata_keywords": "来自API/数据库的关键词、学科、概念或期刊主题，多个用 | 分隔",
        "is_blacklisted": "是否应从最终注释集中剔除: [true | false]",
        "blacklist_reason": "若is_blacklisted=true，说明原因；尤其是无flavin、无酶、纯基因组注释、无正规酶学表征时",
        "key_innovation": "核心创新点(中文一句话)",
        "limitation": "局限性/挑战(中文一句话)",
        "future_direction": "下一步可追踪方向，特别是是否应找同源flavin酶或结构相似酶(中文一句话)",
        "confidence": "提取置信度: 0.0-1.0"
    },
    "lipase": {
        "topic_domain": "主题领域: [光酶催化 | 光氧化还原化学 | 脂肪酶/酯酶催化 | alpha/beta水解酶 | 自由基化学 | 其他]",
        "enzyme_category": "酶催化类型: [天然酶 | 工程化酶 | 协同光催化/酶催化 | 无酶小分子光催化 | 不涉及酶 | 未明确]",
        "enzyme_family": "具体酶家族: [脂肪酶/酯酶 | alpha/beta水解酶 | P450 | ERED/OYE | FAP | 非血红素铁酶 | Aldolase | 不涉及酶 | 其他]",
        "radical_mechanism": "自由基/激发态机制: [光氧化还原(PET) | 能量转移(EnT) | 电子转移(ET) | HAT | EDA复合物 | 辅因子直接激发 | 底物直接激发 | 非自由基 | 未明确]",
        "reaction_type": "反应类型: [C-H官能团化 | C-H芳基化 | C-H烷基化 | 环加成 | 重排 | 脱羧偶联 | 卤化 | Giese加成 | 水解/酯交换 | 对映选择性转化 | 其他]",
        "bond_formation": "形成或断裂的关键键: [C-C | C-N | C-O | C-S | C-H | C-X | ester C-O cleavage | amide C-N cleavage | 未明确]",
        "substrate_class": "底物类别: [羧酸 | 烯烃 | 芳烃/杂芳烃 | 亚胺 | 羰基化合物 | 卤代物 | 磺酰胺 | 酯/脂质 | 酰胺 | 药物片段 | 未明确]",
        "product_class": "产物类别或目标骨架(短语)",
        "stereochemistry": "选择性信息: [对映选择性 | 非对映选择性 | 区域选择性 | 化学选择性 | 无立体控制 | 未明确]",
        "has_mechanistic_study": "是否含机制研究: [是 | 否 | 未明确]",
        "mechanistic_methods": "机制研究手段: [晶体结构 | EPR/ESR | 瞬态吸收光谱 | DFT计算 | 分子动力学 | 同位素标记 | 对照实验 | 无]",
        "light_source": "光源: [蓝光LED | 绿光LED | 可见光 | UV | 太阳光 | 不涉及光 | 未明确]",
        "cofactor_or_photosensitizer": "辅因子/光敏剂/催化剂: [FAD/FMN | NAD(P)H | Ru(bpy)3 | Ir(ppy)3 | 有机染料 | 金属催化剂 | 酶活性位点残基 | 无 | 未明确]",
        "lipase_or_hydrolase_relevance": "与脂肪酶或alpha/beta水解酶研究的关联: [直接相关 | 可借鉴机制 | 可借鉴底物/口袋 | 仅背景相关 | 无明显相关]",
        "key_innovation": "核心创新点(中文一句话)",
        "limitation": "局限性/挑战(中文一句话)",
        "future_direction": "未来方向暗示(中文一句话)",
        "confidence": "提取置信度: 0.0-1.0"
    }
}

# 颜色配置（莫兰迪色系，符合你的审美偏好）
MORANDI_COLORS = {
    'C-H_arylation': '#8B9DAF', 'C-H_alkylation': '#A5B5C5',
    'decarboxylative_coupling': '#B8C5D0', 'halogenation': '#C9B8A8',
    'cycloaddition': '#A8B5A0', 'Giese_addition': '#B5A8B0',
    'HAT': '#A0A8B5', 'other': '#D0D0D0',
    'P450': '#E8D5C4', 'ERED': '#D4C4B0', 'FAP': '#C4B8A8',
    'Aldolase': '#B8A8C4', 'Lipase': '#A8C4B8', 'Other_enzyme': '#C4C4C4',
    'photoredox': '#F0E6D3', 'energy_transfer': '#E6E0F0',
    'direct_excitation': '#E0F0E6', 'ET': '#F0E0E6'
}

EMPTY_TAG_VALUES = {"", "N/A", "NA", "n/a", "None", "none", "未明确", "未知", "无"}
NON_CONCEPT_VALUES = EMPTY_TAG_VALUES | {"unclear"}

FIELD_ALIASES = {
    "enzyme_category": ["enzyme_category", "photoenzymatic_status"],
    "radical_mechanism": ["radical_or_excited_state_mechanism", "radical_mechanism"],
    "stereochemistry": ["selectivity", "stereochemistry"],
    "cofactor": ["catalyst_or_photosensitizer", "cofactor_or_photosensitizer", "cofactor"],
    "mechanistic_methods": ["mechanistic_evidence", "mechanistic_methods"],
}

OTHER_ENZYME_VALUES = {"Other_enzyme", "other", "unclear", "unknown", "其他", "未明确"}
NO_ENZYME_VALUES = {"no enzyme", "不涉及酶", "无酶小分子光催化"}
FALSE_ASE_WORDS = {
    "base", "case", "cease", "disease", "ease", "phase", "phrase", "please",
    "release", "showcase", "purchase", "increase", "decrease", "database",
    "interface", "surface", "chase", "erase", "grease",
}
KNOWN_ENZYME_FAMILIES = {
    "P450", "ERED/OYE", "FAP", "flavoprotein", "nonheme iron enzyme",
    "lipase/esterase", "alpha/beta-hydrolase", "Aldolase", "transaminase",
    "dehydrogenase", "monooxygenase", "peroxidase", "nitroreductase",
    "halogenase", "flavin-dependent halogenase", "photolyase/cryptochrome",
    "berberine bridge enzyme", "amine oxidase", "acyl-CoA dehydrogenase",
    "D-amino acid oxidase", "glucose oxidase", "sulfide:quinone oxidoreductase",
    "luciferase", "flavin reductase", "synthase", "kinase/phosphatase", "protease", "lyase",
    "hydrolase", "oxidase", "reductase", "engineered enzyme",
}

CLASSIFICATION_FIELDS = {
    "topic_domain",
    "paper_domain",
    "search_track",
    "photoenzymatic_status",
    "photoenzymatic_application_status",
    "new_to_nature_reaction",
    "enzyme_category",
    "enzyme_family",
    "flavin_cofactor",
    "characterized_enzyme_evidence",
    "candidate_for_photoenzyme_repurposing",
    "structure_similarity_hint",
    "radical_or_excited_state_mechanism",
    "radical_mechanism",
    "radical_source",
    "reaction_type",
    "bond_formation",
    "substrate_class",
    "selectivity",
    "stereochemistry",
    "mechanistic_evidence",
    "mechanistic_methods",
    "light_source",
    "catalyst_or_photosensitizer",
    "cofactor_or_photosensitizer",
    "cofactor",
    "enzyme_integration_opportunity",
    "lipase_or_hydrolase_relevance",
    "manual_review_priority",
    "metadata_keywords",
}


@dataclass
class ResearchCriteria:
    """User-defined screening anchors used before expensive LLM extraction."""

    field: str = ""
    enzyme_feature: str = ""
    reaction: str = ""
    goal: str = ""

    def enabled(self) -> bool:
        return any([self.field.strip(), self.enzyme_feature.strip(), self.reaction.strip(), self.goal.strip()])

    def hard_fields(self) -> List[Tuple[str, str]]:
        return [
            ("field", self.field.strip()),
            ("enzyme_feature", self.enzyme_feature.strip()),
            ("reaction", self.reaction.strip()),
        ]

    def to_prompt(self) -> str:
        if not self.enabled():
            return ""
        lines = ["Structured screening criteria for this run:"]
        if self.field.strip():
            lines.append(f"- Field/domain: {self.field.strip()}")
        if self.enzyme_feature.strip():
            lines.append(f"- Enzyme feature/cofactor: {self.enzyme_feature.strip()}")
        if self.reaction.strip():
            lines.append(f"- Reaction criterion: {self.reaction.strip()}")
        if self.goal.strip():
            lines.append(f"- Final goal: {self.goal.strip()}")
        lines.append(
            "- Treat field, enzyme feature, and reaction as hard screening anchors before spending LLM calls; "
            "use the final goal to rank evidence-backed candidates for later homolog/tree work."
        )
        return "\n".join(lines)


CRITERIA_STOPWORDS = {
    "about", "across", "after", "again", "against", "also", "among", "and",
    "are", "because", "been", "being", "between", "biology", "both", "but",
    "can", "chemical", "chemistry", "cofactor", "cofactors", "common", "could",
    "domain", "enzyme", "enzymes", "feature", "field", "find", "from", "goal",
    "have", "include", "includes", "including", "into", "literature", "more", "noun", "paper", "papers", "rare",
    "reaction", "reactions", "reported", "review", "search", "synthesis", "that",
    "the", "their", "these", "this", "through", "toward", "under", "use", "used",
    "using", "with", "within", "absent", "biological", "common", "metabolism", "useful",
}

CRITERIA_SYNONYM_GROUPS = {
    "photoenzymatic": [
        "photoenzymatic",
        "photoenzyme",
        "photoenzymology",
        "photobiocatalysis",
        "photobiocatalytic",
        "enzyme photocatalysis",
        "enzymatic photocatalysis",
        "light-driven enzyme",
        "light driven enzyme",
        "light-driven enzymatic",
        "photoinduced enzymatic",
        "photoenzymatic catalysis",
    ],
    "flavin": [
        "flavin",
        "flavoenzyme",
        "flavin-dependent",
        "flavin dependent",
        "FAD",
        "FMN",
        "riboflavin",
        "old yellow enzyme",
        "OYE",
        "ene-reductase",
        "ene reductase",
        "ERED",
        "fatty acid photodecarboxylase",
        "photodecarboxylase",
        "FAP",
    ],
    "new_to_nature": [
        "new-to-nature",
        "new to nature",
        "non-natural",
        "nonnatural",
        "unnatural",
        "abiological",
        "non-native",
        "non native",
        "noncanonical",
        "non-canonical",
        "promiscuous radical",
        "promiscuous reaction",
        "repurposed enzyme",
        "engineered new reaction",
        "synthetically useful reaction",
    ],
    "characterized_enzyme": [
        "purified enzyme",
        "activity assay",
        "enzyme assay",
        "steady-state kinetics",
        "kinetic parameter",
        "kcat",
        "km",
        "substrate scope",
        "crystal structure",
        "x-ray structure",
        "PDB",
        "mutagenesis",
    ],
}


def normalize_match_text(value: str) -> str:
    text = str(value or "").lower()
    for ch in ["\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2212"]:
        text = text.replace(ch, "-")
    text = text.replace("_", " ").replace("/", " ")
    return re.sub(r"\s+", " ", text).strip()


def dedupe_preserve_order(values: List[str]) -> List[str]:
    seen = set()
    result = []
    for value in values:
        value = str(value or "").strip()
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def expand_criteria_terms(raw_value: str, kind: str = "") -> List[str]:
    raw = str(raw_value or "").strip()
    if not raw:
        return []
    low = normalize_match_text(raw)
    terms: Set[str] = set()

    def add_group(group_name: str):
        for alias in CRITERIA_SYNONYM_GROUPS.get(group_name, []):
            terms.add(alias)

    if any(k in low for k in ["photoenzym", "photoenzyme", "photobiocatal", "photocatalytic enzyme"]):
        add_group("photoenzymatic")
    if any(k in raw for k in ["光酶", "光生物催化", "光酶催化"]):
        add_group("photoenzymatic")

    if any(k in low for k in ["flavin", "flavo", "fad", "fmn", "riboflavin", "yellow enzyme", "ene-reductase", "oye"]):
        add_group("flavin")
    if any(k in raw for k in ["黄素", "黄素依赖", "黄素辅酶"]):
        add_group("flavin")

    if any(k in low for k in ["new-to-nature", "new to nature", "non-natural", "nonnatural", "unnatural", "abiological", "non-native"]):
        add_group("new_to_nature")
    if any(k in raw for k in ["非天然", "罕见", "生物体内罕见", "化工领域常见", "新反应"]):
        add_group("new_to_nature")

    if any(k in low for k in ["characterized", "kinetic", "substrate scope", "crystal structure", "pdb", "activity assay"]):
        add_group("characterized_enzyme")

    for phrase in re.split(r"[;,\n，；。]+", raw):
        norm = normalize_match_text(phrase)
        if 3 <= len(norm) <= 70 and norm not in CRITERIA_STOPWORDS:
            terms.add(norm)

    for token in re.findall(r"[a-z0-9][a-z0-9+\-()]{1,}", low):
        clean = token.strip("-()")
        if not clean or clean in CRITERIA_STOPWORDS:
            continue
        if len(clean) >= 4 or clean in {"fad", "fmn", "oye", "fap"}:
            terms.add(clean)

    # Preserve a few useful two-word phrases from arbitrary noun fields.
    words = [w for w in re.findall(r"[a-z0-9][a-z0-9+\-()]{2,}", low) if w not in CRITERIA_STOPWORDS]
    for i in range(len(words) - 1):
        phrase = f"{words[i]} {words[i + 1]}"
        if len(phrase) <= 50:
            terms.add(phrase)

    return sorted(terms, key=lambda x: (-len(x), x))


def term_hits_in_text(terms: List[str], text: str) -> List[str]:
    norm_text = normalize_match_text(text)
    hits = []
    for term in terms:
        norm = normalize_match_text(term)
        if not norm:
            continue
        variants = {norm, norm.replace("-", " "), norm.replace(" ", "-")}
        found = False
        for variant in variants:
            if not variant:
                continue
            if len(variant) <= 3:
                pattern = rf"(?<![a-z0-9]){re.escape(variant)}(?![a-z0-9])"
                found = re.search(pattern, norm_text) is not None
            else:
                found = variant in norm_text
            if found:
                break
        if found:
            hits.append(term)
    return dedupe_preserve_order(hits)


def criteria_context_for_paper(paper: "PaperNode") -> str:
    tags = paper.tags or {}
    fields = [
        paper.title,
        paper.abstract,
        paper.journal,
        tags.get("metadata_keywords", ""),
        tags.get("enzyme_family", ""),
        tags.get("enzyme_name_or_target", ""),
        tags.get("flavin_cofactor", ""),
        tags.get("paper_domain", ""),
        tags.get("search_track", ""),
        tags.get("reaction_type", ""),
        tags.get("bond_formation", ""),
        tags.get("substrate_class", ""),
        tags.get("characterized_enzyme_evidence", ""),
        tags.get("structure_similarity_hint", ""),
    ]
    return "\n".join(str(f) for f in fields if f)


def apply_research_criteria_to_paper(paper: "PaperNode", criteria: ResearchCriteria) -> str:
    if not criteria.enabled():
        return ""
    context = criteria_context_for_paper(paper)
    hard_results = {}
    passed = 0
    required = 0
    for label, raw_value in criteria.hard_fields():
        if not raw_value:
            continue
        required += 1
        terms = expand_criteria_terms(raw_value, label)
        hits = term_hits_in_text(terms, context)
        hard_results[label] = hits
        if hits:
            passed += 1
        paper.tags[f"criteria_{label}_hits"] = " | ".join(hits[:12])

    missing = [label for label, raw_value in criteria.hard_fields() if raw_value and not hard_results.get(label)]
    feature_hit = bool(hard_results.get("enzyme_feature"))
    field_hit = bool(hard_results.get("field"))
    reaction_hit = bool(hard_results.get("reaction"))
    search_track = tag_get(paper.tags, "search_track", default="")
    candidate_rank = tag_get(paper.tags, "candidate_for_photoenzyme_repurposing", default="")
    evidence = tag_get(paper.tags, "characterized_enzyme_evidence", default="")

    if required == 0:
        status = "criteria_pass"
    elif passed == required:
        status = "criteria_pass"
    elif (
        feature_hit
        and (
            field_hit
            or reaction_hit
            or search_track == "B_characterized_flavin_enzyme_not_photoenzymatic"
            or candidate_rank in {"high", "medium"}
            or evidence in {"purified enzyme kinetics", "activity assay only", "structure only"}
        )
    ):
        status = "criteria_borderline"
    elif required >= 3 and passed >= 2:
        status = "criteria_borderline"
    else:
        status = "criteria_fail"

    paper.tags["criteria_status"] = status
    paper.tags["criteria_score"] = f"{passed}/{required}"
    paper.tags["criteria_missing"] = " | ".join(missing)
    if criteria.goal.strip():
        paper.tags["criteria_goal"] = criteria.goal.strip()
    if criteria.field.strip():
        if field_hit:
            paper.tags["field_match_status"] = "field_hit"
        elif feature_hit:
            paper.tags["field_match_status"] = "enzyme_only_field_missing"
        else:
            paper.tags["field_match_status"] = "field_missing"
    elif criteria.enzyme_feature.strip():
        paper.tags["field_match_status"] = "enzyme_scope_no_field"
    else:
        paper.tags["field_match_status"] = "no_field_constraint"
    if status == "criteria_pass":
        paper.tags["manual_review_priority"] = "high"
    elif status == "criteria_borderline" and tag_get(paper.tags, "manual_review_priority", default="low") == "low":
        paper.tags["manual_review_priority"] = "medium"
    return status


def apply_research_criteria(papers: List["PaperNode"], criteria: ResearchCriteria) -> Counter:
    counts = Counter()
    if not criteria.enabled():
        return counts
    for paper in papers:
        counts[apply_research_criteria_to_paper(paper, criteria)] += 1
    print(
        "[Criteria] "
        f"pass={counts.get('criteria_pass', 0)}, "
        f"borderline={counts.get('criteria_borderline', 0)}, "
        f"fail={counts.get('criteria_fail', 0)}"
    )
    return counts


def filter_papers_by_criteria(papers: List["PaperNode"], output_dir: Path) -> List["PaperNode"]:
    failed = [p for p in papers if p.tags.get("criteria_status") == "criteria_fail"]
    kept = [p for p in papers if p.tags.get("criteria_status") != "criteria_fail"]
    if failed:
        rows = []
        for paper in failed:
            rows.append({
                "title": paper.title,
                "doi": normalize_doi_for_filter(paper.doi),
                "year": paper.year,
                "journal": paper.journal,
                "criteria_score": paper.tags.get("criteria_score", ""),
                "criteria_missing": paper.tags.get("criteria_missing", ""),
                "field_hits": paper.tags.get("criteria_field_hits", ""),
                "enzyme_feature_hits": paper.tags.get("criteria_enzyme_feature_hits", ""),
                "reaction_hits": paper.tags.get("criteria_reaction_hits", ""),
                "paper_domain": tag_get(paper.tags, "paper_domain", default=""),
                "enzyme_family": tag_get(paper.tags, "enzyme_family", default=""),
                "reaction_type": tag_get(paper.tags, "reaction_type", default=""),
            })
        path = output_dir / "criteria_filtered_papers.csv"
        pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")
        print(f"[Criteria] Filtered out {len(failed)} criteria-fail papers for this run: {path}")
    return kept


def read_research_criteria_file(path_value: str) -> Dict[str, str]:
    if not path_value:
        return {}
    path = Path(path_value)
    text = path.read_text(encoding="utf-8-sig").strip()
    if not text:
        return {}
    if path.suffix.lower() == ".json":
        data = json.loads(text)
        return {str(k): str(v) for k, v in data.items() if v is not None}

    key_map = {
        "field": "field",
        "domain": "field",
        "领域": "field",
        "enzyme": "enzyme_feature",
        "enzyme_feature": "enzyme_feature",
        "cofactor": "enzyme_feature",
        "酶": "enzyme_feature",
        "辅因子": "enzyme_feature",
        "reaction": "reaction",
        "反应": "reaction",
        "goal": "goal",
        "objective": "goal",
        "target": "goal",
        "目标": "goal",
    }
    values: Dict[str, List[str]] = defaultdict(list)
    unkeyed = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        match = re.match(r"^([^:：]{1,40})[:：]\s*(.+)$", stripped)
        if match:
            key = match.group(1).strip().lower()
            value = match.group(2).strip()
            mapped = key_map.get(key)
            if mapped:
                values[mapped].append(value)
            else:
                unkeyed.append(stripped)
        else:
            unkeyed.append(stripped)
    if unkeyed:
        values["goal"].append(" ".join(unkeyed))
    return {k: " ".join(v).strip() for k, v in values.items()}


def build_research_criteria(args) -> ResearchCriteria:
    file_values = read_research_criteria_file(getattr(args, "criteria_file", "") or "")

    def merged(name: str, arg_name: str) -> str:
        chunks = []
        if file_values.get(name):
            chunks.append(file_values[name])
        for value in getattr(args, arg_name, None) or []:
            if str(value).strip():
                chunks.append(str(value).strip())
        return " ".join(chunks).strip()

    return ResearchCriteria(
        field=merged("field", "field"),
        enzyme_feature=merged("enzyme_feature", "enzyme_feature"),
        reaction=merged("reaction", "reaction"),
        goal=merged("goal", "goal"),
    )


def query_anchor_from_text(raw_value: str, preferred_group: str = "") -> str:
    raw = str(raw_value or "").strip()
    low = normalize_match_text(raw)
    if not raw:
        return ""
    if preferred_group == "field" and any(k in low for k in ["photoenzym", "photoenzyme", "photobiocatal"]):
        return '"photoenzymatic"'
    if preferred_group == "enzyme" and any(k in low for k in ["flavin", "fad", "fmn", "yellow enzyme", "flavo"]):
        return '"flavin-dependent"'
    if preferred_group == "reaction" and any(k in low for k in ["new-to-nature", "new to nature", "non-natural", "unnatural", "abiological"]):
        return '"new-to-nature"'
    tokens = [t for t in re.findall(r"[a-z0-9][a-z0-9+\-]{2,}", low) if t not in CRITERIA_STOPWORDS]
    if not tokens:
        return ""
    phrase = " ".join(tokens[:3])
    return f'"{phrase}"' if " " in phrase else phrase


def default_search_queries_for_criteria(criteria: ResearchCriteria) -> List[str]:
    if not criteria.enabled():
        return []
    field = query_anchor_from_text(criteria.field, "field")
    enzyme = query_anchor_from_text(criteria.enzyme_feature, "enzyme")
    reaction = query_anchor_from_text(criteria.reaction, "reaction")
    queries = []
    if field and enzyme and reaction:
        queries.append(f"{field} {enzyme} {reaction}")
    if field and enzyme:
        queries.append(f"{field} {enzyme} enzyme")
    if enzyme and reaction:
        queries.append(f"{enzyme} enzyme {reaction}")
    if field and reaction:
        queries.append(f"{field} {reaction}")
    if enzyme:
        queries.append(f"{enzyme} enzyme substrate scope kinetics structure")
    return dedupe_preserve_order(queries)


def is_empty_tag(value) -> bool:
    """True for placeholders that should be backfilled or hidden in summaries."""
    if value is None:
        return True
    if isinstance(value, (list, tuple, set)):
        return not any(not is_empty_tag(v) for v in value)
    text = str(value).strip()
    return text in EMPTY_TAG_VALUES or text.lower().startswith("n/a")


def is_non_concept_tag(value) -> bool:
    """True for values that may be displayed but should not form graph concepts."""
    if value is None:
        return True
    if isinstance(value, (list, tuple, set)):
        return not any(not is_non_concept_tag(v) for v in value)
    text = str(value).strip()
    return is_empty_tag(text) or text in NON_CONCEPT_VALUES


def tag_get(tags: Dict, *fields: str, default: str = "N/A") -> str:
    """Read a tag using aliases across schema profiles."""
    for field in fields:
        for candidate in FIELD_ALIASES.get(field, [field]):
            value = tags.get(candidate)
            if not is_empty_tag(value):
                return str(value).strip()
    return default


def merge_with_fallback(primary: Dict[str, str], fallback: Dict[str, str], schema: Dict[str, str]) -> Dict[str, str]:
    """Fill empty LLM fields with deterministic title/abstract keyword tags."""
    merged = dict(primary or {})
    for key in schema.keys():
        if is_empty_tag(merged.get(key)) and not is_empty_tag(fallback.get(key)):
            merged[key] = fallback[key]
        elif is_empty_tag(merged.get(key)) and key in CLASSIFICATION_FIELDS:
            merged[key] = "unclear"
        elif key not in merged:
            merged[key] = "unclear" if key in CLASSIFICATION_FIELDS else "N/A"

    alias_pairs = [
        ("radical_mechanism", "radical_or_excited_state_mechanism"),
        ("stereochemistry", "selectivity"),
        ("cofactor", "catalyst_or_photosensitizer"),
        ("mechanistic_methods", "mechanistic_evidence"),
        ("enzyme_category", "photoenzymatic_status"),
    ]
    for legacy, modern in alias_pairs:
        if legacy not in merged and modern in merged:
            merged[legacy] = merged[modern]
        if modern not in merged and legacy in merged:
            merged[modern] = merged[legacy]
    return merged


def parse_float(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def parse_int(value, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def normalize_enzyme_family(value: str, title: str = "", abstract: str = "") -> str:
    """Clean noisy enzyme-family labels, especially false *ase matches."""
    raw = str(value or "").strip()
    if is_empty_tag(raw):
        return "unclear"
    low = raw.lower()
    if low in {v.lower() for v in NO_ENZYME_VALUES}:
        return "no enzyme"
    if low in {v.lower() for v in OTHER_ENZYME_VALUES}:
        return "Other_enzyme"

    for known in KNOWN_ENZYME_FAMILIES:
        if low == known.lower():
            return known

    text = f"{title} {abstract}".lower()
    if low in FALSE_ASE_WORDS:
        return "no enzyme" if not any(k in text for k in ["enzyme", "enzymatic", "biocatal", "protein"]) else "Other_enzyme"

    if low.endswith("ase"):
        enzyme_context = any(k in text for k in [
            "enzyme", "enzymatic", "biocatal", "protein", "active site",
            "catalyzes", "catalyses", "mutant", "variant",
        ])
        return raw if enzyme_context else "no enzyme"

    if "enzyme" in low:
        return raw
    return raw


def is_other_enzyme_family(value: str) -> bool:
    normalized = normalize_enzyme_family(value)
    return normalized in OTHER_ENZYME_VALUES or normalized == "Other_enzyme"


def should_blacklist_paper(paper: PaperNode) -> Tuple[bool, str]:
    """Blacklist papers that are not useful for annotation after extraction."""
    if str(paper.tags.get("is_blacklisted", "")).lower() == "true":
        return True, paper.tags.get("blacklist_reason", "LLM判定不适合纳入")

    family = normalize_enzyme_family(tag_get(paper.tags, "enzyme_family"), paper.title, paper.abstract)
    paper_domain = tag_get(paper.tags, "paper_domain", default="")
    if paper_domain == "unrelated/blacklist":
        return True, "paper_domain=unrelated/blacklist"
    if paper_domain == "organic synthesis no-enzyme":
        return True, "纯有机合成且缺少光酶/自由基迁移价值"
    if family == "no enzyme":
        status = tag_get(paper.tags, "photoenzymatic_status", default="")
        opportunity = tag_get(paper.tags, "enzyme_integration_opportunity", default="")
        radical_mechanism = tag_get(paper.tags, "radical_mechanism", default="")
        reaction_type = tag_get(paper.tags, "reaction_type", default="")
        has_transfer_value = (
            paper_domain == "organic photoredox/radical no-enzyme"
            or status in {"已接入photoenzymatic", "可设想接入酶体系", "仅光化学/自由基反应"}
            and opportunity in {"high", "medium"}
            and radical_mechanism not in {"", "unclear", "N/A"}
            and reaction_type not in {"", "other", "unclear", "N/A"}
        )
        if not has_transfer_value:
            return True, "no enzyme且缺少明确酶口袋迁移价值"

    if paper.tags.get("is_discovered") == "true" and not normalize_doi_for_filter(paper.doi):
        return True, "discovered paper missing DOI"
    return False, ""


def normalize_doi_for_filter(doi: str) -> str:
    doi = (doi or "").strip()
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi, flags=re.I)
    doi = doi.replace("doi:", "").replace("DOI:", "").strip()
    if doi.lower() in {"", "nan", "none", "null"}:
        return ""
    return doi


def normalize_title_for_match(title: str) -> str:
    return re.sub(r"\W+", " ", (title or "").lower()).strip()


def is_broad_subject_label(value: str) -> bool:
    """Broad fields from metadata APIs are not journal titles."""
    label = re.sub(r"\s+", " ", (value or "").strip().lower())
    return label in {
        "medicine", "chemistry", "biology", "computer science", "materials science",
        "physics", "engineering", "psychology", "sociology", "political science",
        "economics", "mathematics", "geology", "geography", "history", "art",
        "business", "philosophy", "environmental science", "agricultural and food sciences",
        "education", "law",
    }


def should_replace_journal(value: str) -> bool:
    label = (value or "").strip()
    return not label or label.lower() in {"unknown", "nan", "none", "null"} or is_broad_subject_label(label)


def build_zotero_index_from_papers(papers: List[PaperNode]) -> Dict[str, Set[str]]:
    index = {"dois": set(), "titles": set(), "journals": set()}
    for paper in papers:
        doi = normalize_doi_for_filter(paper.doi).lower()
        if doi:
            index["dois"].add(doi)
        title = normalize_title_for_match(paper.title)
        if title:
            index["titles"].add(title)
    return index


def build_zotero_index_from_sqlite(sqlite_path: str) -> Dict[str, Set[str]]:
    index = {"dois": set(), "titles": set(), "journals": set()}
    if not sqlite_path:
        return index
    path = Path(sqlite_path)
    if not path.exists():
        print(f"[Zotero] SQLite不存在: {sqlite_path}")
        return index
    try:
        conn = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
        cur = conn.cursor()
        cur.execute("""
            SELECT itemDataValues.value
            FROM itemDataValues
            JOIN itemData ON itemData.valueID = itemDataValues.valueID
            JOIN fields ON fields.fieldID = itemData.fieldID
            WHERE lower(fields.fieldName) = 'doi'
        """)
        for (doi,) in cur.fetchall():
            norm = normalize_doi_for_filter(str(doi)).lower()
            if norm:
                index["dois"].add(norm)
        cur.execute("""
            SELECT itemDataValues.value
            FROM itemDataValues
            JOIN itemData ON itemData.valueID = itemDataValues.valueID
            JOIN fields ON fields.fieldID = itemData.fieldID
            WHERE lower(fields.fieldName) = 'title'
        """)
        for (title,) in cur.fetchall():
            norm = normalize_title_for_match(str(title))
            if norm:
                index["titles"].add(norm)
        cur.execute("""
            SELECT itemDataValues.value
            FROM itemDataValues
            JOIN itemData ON itemData.valueID = itemDataValues.valueID
            JOIN fields ON fields.fieldID = itemData.fieldID
            WHERE lower(fields.fieldName) IN (
                'publicationtitle', 'journaltitle', 'proceedingstitle', 'conferencename'
            )
        """)
        for (journal,) in cur.fetchall():
            label = re.sub(r"\s+", " ", str(journal or "").strip())
            if label:
                index["journals"].add(label)
        conn.close()
        print(f"[Zotero] SQLite索引: DOI {len(index['dois'])}, title {len(index['titles'])}")
    except Exception as e:
        print(f"[Zotero] SQLite索引失败: {e}")
    return index


def merge_zotero_indices(*indices: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    merged = {"dois": set(), "titles": set(), "journals": set()}
    for index in indices:
        merged["dois"].update(index.get("dois", set()))
        merged["titles"].update(index.get("titles", set()))
        merged["journals"].update(index.get("journals", set()))
    return merged


def annotate_zotero_status(papers: List[PaperNode], seed_index: Dict[str, Set[str]], library_index: Optional[Dict[str, Set[str]]] = None):
    library_index = library_index or seed_index
    for paper in papers:
        doi = normalize_doi_for_filter(paper.doi).lower()
        title = normalize_title_for_match(paper.title)
        in_seed = (doi and doi in seed_index.get("dois", set())) or (title and title in seed_index.get("titles", set()))
        in_library = (doi and doi in library_index.get("dois", set())) or (title and title in library_index.get("titles", set()))
        if in_seed:
            status = "zotero_seed"
        elif in_library:
            status = "already_in_zotero"
        else:
            status = "new_discovery"
        paper.tags["zotero_status"] = status


def rank_paper_for_quota(paper: PaperNode) -> Tuple[float, int, int]:
    """Prefer relevant/cited/recent papers when limiting ambiguous buckets."""
    return (
        parse_float(paper.tags.get("relevance_score"), 0.0),
        parse_int(paper.tags.get("citation_count"), 0),
        paper.year or 0,
    )


def paper_memory_key(paper: PaperNode) -> str:
    doi = normalize_doi_for_filter(paper.doi).lower()
    if doi:
        return f"doi:{doi}"
    title = re.sub(r"\s+", " ", (paper.title or "").lower()).strip()
    return f"title:{hashlib.md5(title.encode('utf-8')).hexdigest()}" if title else ""


def should_run_llm_for_paper(paper: PaperNode, base_tags: Dict[str, str], scope: str = "all") -> bool:
    if scope == "all":
        return True
    criteria_status = paper.tags.get("criteria_status", "")
    if scope == "criteria":
        return criteria_status in {"criteria_pass", "criteria_borderline"}
    if scope == "criteria-ambiguous":
        if criteria_status == "criteria_fail":
            return False
        if criteria_status in {"criteria_pass", "criteria_borderline"}:
            family = tag_get(base_tags, "enzyme_family")
            domain = tag_get(base_tags, "paper_domain")
            return (
                is_other_enzyme_family(family)
                or family in {"Other_enzyme", "unclear", "no enzyme"}
                or domain in {"", "unclear", "unrelated/blacklist"}
                or is_empty_tag(paper.abstract)
                or tag_get(base_tags, "flavin_cofactor", default="") in {"", "unclear"}
                or tag_get(base_tags, "new_to_nature_reaction", default="") in {"", "unclear"}
                or tag_get(base_tags, "characterized_enzyme_evidence", default="") in {"", "unclear"}
            )
        return False
    family = tag_get(base_tags, "enzyme_family")
    domain = tag_get(base_tags, "paper_domain")
    if scope == "enzyme-family":
        flavin = tag_get(base_tags, "flavin_cofactor", default="")
        text = criteria_context_for_paper(paper).lower()
        flavin_signal = flavin in {"FAD", "FMN", "flavin unspecified"} or any(
            k in text for k in ["flavin", "fad", "fmn", "flavoenzyme", "old yellow enzyme", "ene-reductase", "ered", "oye"]
        )
        return (
            is_other_enzyme_family(family)
            or family in {"Other_enzyme", "unclear", "no enzyme"}
            or domain in {"", "unclear", "unrelated/blacklist"}
            or flavin_signal
        )
    if scope == "other-only":
        return is_other_enzyme_family(family) or family in {"Other_enzyme", "unclear"} or domain in {"", "unclear"}
    if scope == "ambiguous":
        return (
            is_other_enzyme_family(family)
            or family in {"Other_enzyme", "unclear", "no enzyme"}
            or domain in {"", "unclear", "unrelated/blacklist"}
            or is_empty_tag(paper.abstract)
        )
    return True


class LiteratureMemory:
    """Persistent blacklist/review memory, modeled after the SSN JSON caches."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.records: Dict[str, Dict] = {}
        self.load()

    def load(self):
        if not self.path.exists():
            self.records = {}
            return
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                raw = json.load(f)
            self.records = raw.get("records", raw if isinstance(raw, dict) else {})
            print(f"[Memory] 加载黑名单/注释记忆 {len(self.records)} 条: {self.path}")
        except Exception as e:
            print(f"[Memory] 加载失败: {e}")
            self.records = {}

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "version": 1,
            "updated_at": pd.Timestamp.now().isoformat(),
            "records": self.records,
        }
        with open(self.path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        print(f"[Memory] 已保存黑名单/注释记忆 {len(self.records)} 条: {self.path}")

    def apply(self, paper: PaperNode) -> bool:
        key = paper_memory_key(paper)
        if not key or key not in self.records:
            return False
        rec = self.records[key]
        paper.tags.update({
            "memory_hit": "true",
            "memory_status": rec.get("status", ""),
            "paper_domain": rec.get("paper_domain", ""),
            "enzyme_family": rec.get("enzyme_family", paper.tags.get("enzyme_family", "")),
        })
        if rec.get("status") in {"deferred", "review_deferred"}:
            paper.tags["review_status"] = "deferred_from_memory"
            paper.tags["review_reason"] = rec.get("reason", "")
            return False
        if rec.get("status") == "blacklisted":
            reason = rec.get("reason", "命中永久黑名单")
            if "Other_enzyme超过" in reason or "Other_enzyme overflow" in reason:
                paper.tags["review_status"] = "deferred_from_legacy_blacklist"
                paper.tags["review_reason"] = reason
                paper.tags["is_blacklisted"] = "false"
                return False
            paper.tags["is_blacklisted"] = "true"
            paper.tags["blacklist_reason"] = reason
            return True
        return False

    def remember(self, paper: PaperNode, status: str, reason: str = ""):
        key = paper_memory_key(paper)
        if not key:
            return
        self.records[key] = {
            "status": status,
            "reason": reason,
            "title": paper.title,
            "doi": normalize_doi_for_filter(paper.doi),
            "year": paper.year,
            "journal": paper.journal,
            "paper_domain": tag_get(paper.tags, "paper_domain", default=""),
            "enzyme_family": tag_get(paper.tags, "enzyme_family", default=""),
            "topic_domain": tag_get(paper.tags, "topic_domain", default=""),
            "review_status": paper.tags.get("review_status", ""),
            "updated_at": pd.Timestamp.now().isoformat(),
        }

    def migrate_legacy_other_blacklist(self) -> int:
        changed = 0
        for rec in self.records.values():
            if rec.get("status") == "blacklisted":
                reason = str(rec.get("reason", ""))
                if "Other_enzyme超过" in reason or "Other_enzyme overflow" in reason:
                    rec["status"] = "review_deferred"
                    rec["review_status"] = "deferred_legacy_other_overflow"
                    rec["updated_at"] = pd.Timestamp.now().isoformat()
                    changed += 1
        if changed:
            print(f"[Memory] Migrated {changed} legacy Other_enzyme quota blacklists to review_deferred")
        return changed


def load_text_cache(path: Path) -> Dict:
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def save_text_cache(path: Path, data: Dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def fetch_crossref_metadata(doi: str, timeout: float = 20.0) -> Dict[str, str]:
    doi = normalize_doi_for_filter(doi)
    if not doi:
        return {}
    try:
        resp = requests.get(f"https://api.crossref.org/works/{doi}", timeout=timeout)
        if not resp.ok:
            return {}
        msg = resp.json().get("message", {})
        subjects = msg.get("subject") or []
        containers = msg.get("container-title") or []
        journal = containers[0] if containers else ""
        abstract = re.sub(r"<[^>]+>", " ", msg.get("abstract") or "")
        return {
            "crossref_journal": journal,
            "crossref_subjects": " | ".join(subjects),
            "crossref_abstract": re.sub(r"\s+", " ", abstract).strip(),
        }
    except Exception:
        return {}


def fetch_openalex_metadata(doi: str, timeout: float = 20.0) -> Dict[str, str]:
    doi = normalize_doi_for_filter(doi)
    if not doi:
        return {}
    try:
        resp = requests.get(f"https://api.openalex.org/works/https://doi.org/{doi}", timeout=timeout)
        if not resp.ok:
            return {}
        data = resp.json()
        concepts = [c.get("display_name", "") for c in data.get("concepts", [])[:12] if c.get("display_name")]
        keywords = [k.get("display_name", "") for k in data.get("keywords", [])[:12] if k.get("display_name")]
        journal = ""
        source = ((data.get("primary_location") or {}).get("source") or {})
        if isinstance(source, dict):
            journal = source.get("display_name") or ""
        inverted = data.get("abstract_inverted_index") or {}
        abstract = ""
        if inverted:
            positions = []
            for word, idxs in inverted.items():
                for idx in idxs:
                    positions.append((idx, word))
            abstract = " ".join(word for _, word in sorted(positions))
        return {
            "openalex_journal": journal,
            "openalex_concepts": " | ".join(concepts),
            "openalex_keywords": " | ".join(keywords),
            "openalex_abstract": abstract,
        }
    except Exception:
        return {}


def extract_introduction_from_html(raw_html: str, max_chars: int = 7000) -> str:
    text = re.sub(r"(?is)<(script|style).*?>.*?</\1>", " ", raw_html or "")
    text = re.sub(r"(?i)</(p|div|section|h[1-4])>", "\n", text)
    text = re.sub(r"(?s)<[^>]+>", " ", text)
    text = html_lib.unescape(text)
    text = re.sub(r"[ \t\r\f\v]+", " ", text)
    text = re.sub(r"\n\s+", "\n", text)
    match = re.search(
        r"(?is)\b(1\.\s*)?(introduction|background)\b(.{500,}?)\n\s*(2\.|results|discussion|materials|methods|conclusion|references)\b",
        text,
    )
    if match:
        return re.sub(r"\s+", " ", match.group(0)).strip()[:max_chars]
    intro = re.search(r"(?is)\b(introduction|background)\b(.{500,})", text)
    if intro:
        return re.sub(r"\s+", " ", intro.group(0)).strip()[:max_chars]
    return ""


def fetch_introduction_snippet(paper: PaperNode, timeout: float = 25.0) -> str:
    """Best-effort HTML DOI landing-page introduction fetch; PDF parsing is a future layer."""
    candidates = []
    doi = normalize_doi_for_filter(paper.doi)
    if doi:
        candidates.append(f"https://doi.org/{doi}")
    if paper.url:
        candidates.append(paper.url)
    headers = {"User-Agent": "LiteratureAgent/1.0 (+metadata enrichment)"}
    for url in candidates:
        try:
            resp = requests.get(url, timeout=timeout, headers=headers, allow_redirects=True)
            ctype = resp.headers.get("content-type", "").lower()
            if resp.ok and "html" in ctype:
                intro = extract_introduction_from_html(resp.text)
                if intro:
                    return intro
        except Exception:
            continue
    return ""


def enrich_external_metadata(papers: List[PaperNode], cache_path: Path, enabled: bool = True):
    if not enabled:
        return
    cache = load_text_cache(cache_path)
    changed = False
    for i, paper in enumerate(papers, 1):
        doi = normalize_doi_for_filter(paper.doi)
        if not doi:
            continue
        key = doi.lower()
        meta = cache.get(key)
        if meta is None:
            meta = {}
        if not meta:
            print(f"  [{i}/{len(papers)}] 补充Crossref/OpenAlex元数据: {paper.get_display_title(50)}")
            meta = {}
            meta.update(fetch_crossref_metadata(doi))
            meta.update(fetch_openalex_metadata(doi))
            cache[key] = meta
            changed = True
            time.sleep(0.2)
        elif should_replace_journal(paper.journal) and not (meta.get("crossref_journal") or meta.get("openalex_journal")):
            refresh = {}
            refresh.update(fetch_crossref_metadata(doi))
            refresh.update(fetch_openalex_metadata(doi))
            if refresh:
                meta.update({k: v for k, v in refresh.items() if v})
                cache[key] = meta
                changed = True
                time.sleep(0.2)
        if meta:
            keywords = " | ".join(
                v for k, v in meta.items()
                if k.endswith("subjects") or k.endswith("concepts") or k.endswith("keywords")
                if v
            )
            if keywords:
                paper.tags["metadata_keywords"] = keywords
            for journal_key in ("crossref_journal", "openalex_journal"):
                journal = (meta.get(journal_key) or "").strip()
                if journal and should_replace_journal(paper.journal):
                    paper.journal = journal
                    paper.tags["journal_source"] = journal_key
                    break
            for abstract_key in ("crossref_abstract", "openalex_abstract"):
                if is_empty_tag(paper.abstract) and meta.get(abstract_key):
                    paper.abstract = meta[abstract_key]
    if changed:
        save_text_cache(cache_path, cache)


def apply_annotation_quality_filter(
    papers: List[PaperNode],
    target_non_other: int = 100,
    max_other: int = 20,
) -> Tuple[List[PaperNode], List[PaperNode]]:
    """Filter hard rejects while keeping Other_enzyme overflow as review-deferred, not blacklist."""
    kept_non_other: List[PaperNode] = []
    other_candidates: List[PaperNode] = []
    blacklisted: List[PaperNode] = []

    for paper in papers:
        family = normalize_enzyme_family(tag_get(paper.tags, "enzyme_family"), paper.title, paper.abstract)
        paper.tags["enzyme_family"] = family
        reject, reason = should_blacklist_paper(paper)
        if reject:
            paper.tags["is_blacklisted"] = "true"
            paper.tags["blacklist_reason"] = reason
            blacklisted.append(paper)
            continue
        if is_other_enzyme_family(family):
            other_candidates.append(paper)
        else:
            kept_non_other.append(paper)

    other_candidates.sort(key=rank_paper_for_quota, reverse=True)
    kept_other = other_candidates[:max_other]
    for paper in other_candidates[max_other:]:
        paper.tags["review_status"] = "deferred_other_overflow"
        paper.tags["review_reason"] = f"Other_enzyme超过优先展示上限 {max_other}，未拉黑；需要酶家族精查"
        paper.tags["manual_review_priority"] = "high"

    kept = kept_non_other + kept_other + other_candidates[max_other:]
    if len(kept_non_other) < target_non_other:
        print(
            f"[Quality] 非Other可注释文献 {len(kept_non_other)} 篇，低于目标 {target_non_other}。"
            "建议提高 --breadth-limit/--max-depth 或降低 --relevance-threshold 后重跑 discovery。"
        )
    print(
        f"[Quality] 保留 {len(kept)} 篇：非Other {len(kept_non_other)}，"
        f"Other优先 {len(kept_other)}/{len(other_candidates)}，"
        f"Other待精查 {max(0, len(other_candidates) - len(kept_other))}，黑名单 {len(blacklisted)}。"
    )
    return kept, blacklisted


def write_blacklist_report(blacklisted: List[PaperNode], output_dir: Path):
    path = output_dir / "blacklisted_papers.md"
    lines = [
        "# 黑名单文献",
        "",
        f"共 {len(blacklisted)} 篇被排除出最终注释/聚类。",
        "",
    ]
    for paper in sorted(blacklisted, key=rank_paper_for_quota, reverse=True):
        lines.extend([
            f"## {paper.title}",
            "",
            f"- **Reason**: {paper.tags.get('blacklist_reason', '未说明')}",
            f"- **DOI**: {paper.doi if paper.doi else 'missing'}",
            f"- **Year**: {paper.year}",
            f"- **Journal**: {paper.journal}",
            f"- **Enzyme family**: {tag_get(paper.tags, 'enzyme_family')}",
            f"- **Relevance**: {paper.tags.get('relevance_score', '')}",
            f"- **Citations**: {paper.tags.get('citation_count', '')}",
            "",
        ])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[Quality] 黑名单报告已保存: {path}")


def write_deferred_review_report(papers: List[PaperNode], output_dir: Path):
    deferred = [p for p in papers if p.tags.get("review_status", "") == "deferred_other_overflow"]
    if not deferred:
        return
    deferred.sort(key=rank_paper_for_quota, reverse=True)
    rows = []
    lines = [
        "# Deferred Enzyme-Family Review",
        "",
        "These papers were not blacklisted. They exceeded the current Other_enzyme display quota and need enzyme-family resolution.",
        "",
    ]
    for paper in deferred:
        row = {
            "title": paper.title,
            "doi": normalize_doi_for_filter(paper.doi),
            "year": paper.year,
            "journal": paper.journal,
            "review_reason": paper.tags.get("review_reason", ""),
            "enzyme_family": tag_get(paper.tags, "enzyme_family"),
            "paper_domain": tag_get(paper.tags, "paper_domain"),
            "flavin_cofactor": tag_get(paper.tags, "flavin_cofactor", default=""),
            "criteria_status": paper.tags.get("criteria_status", ""),
            "criteria_score": paper.tags.get("criteria_score", ""),
            "search_track": tag_get(paper.tags, "search_track", default=""),
            "candidate_for_photoenzyme_repurposing": tag_get(paper.tags, "candidate_for_photoenzyme_repurposing", default=""),
            "relevance_score": paper.tags.get("relevance_score", ""),
            "citation_count": paper.tags.get("citation_count", ""),
            "abstract": paper.abstract,
        }
        rows.append(row)
        doi = row["doi"] or "missing"
        lines.extend([
            f"## {paper.title}",
            "",
            f"- DOI: {doi}",
            f"- Year: {paper.year}",
            f"- Journal: {paper.journal}",
            f"- Current family: {row['enzyme_family']}",
            f"- Reason: {row['review_reason']}",
            f"- Flavin: {row['flavin_cofactor']}",
            f"- Criteria: {row['criteria_status']} {row['criteria_score']}",
            "",
        ])
    md_path = output_dir / "deferred_enzyme_review.md"
    csv_path = output_dir / "deferred_enzyme_review.csv"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    pd.DataFrame(rows).to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"[Quality] Deferred enzyme review saved: {md_path} / {csv_path}")


def write_zotero_overlap_report(
    papers: List[PaperNode],
    seed_index: Dict[str, Set[str]],
    library_index: Dict[str, Set[str]],
    output_dir: Path,
):
    path = output_dir / "zotero_overlap_report.md"
    status_counts = Counter(p.tags.get("zotero_status", "") or "unknown" for p in papers)
    lines = [
        "# Zotero Overlap Report",
        "",
        f"- Seed CSV index: DOI {len(seed_index.get('dois', set()))}, title {len(seed_index.get('titles', set()))}",
        f"- Full Zotero index: DOI {len(library_index.get('dois', set()))}, title {len(library_index.get('titles', set()))}, journal {len(library_index.get('journals', set()))}",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in status_counts.most_common():
        lines.append(f"- **{status}**: {count}")
    lines.extend(["", "## New Discoveries", ""])
    new_papers = [p for p in papers if p.tags.get("zotero_status") == "new_discovery"]
    for paper in sorted(new_papers, key=rank_paper_for_quota, reverse=True)[:200]:
        doi = normalize_doi_for_filter(paper.doi)
        doi_text = f"[{doi}](https://doi.org/{doi})" if doi else "missing"
        lines.extend([
            f"### {paper.title}",
            f"- DOI: {doi_text}",
            f"- Year: {paper.year}",
            f"- Journal: {paper.journal}",
            f"- Domain: {tag_get(paper.tags, 'paper_domain')}",
            f"- Field match: {paper.tags.get('field_match_status', '')}",
            f"- Criteria: {paper.tags.get('criteria_status', '')} {paper.tags.get('criteria_score', '')}; missing={paper.tags.get('criteria_missing', '')}",
            f"- Enzyme family: {tag_get(paper.tags, 'enzyme_family')}",
            "",
        ])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[Zotero] Overlap report saved: {path}")


def paper_to_dashboard_record(paper: PaperNode) -> Dict:
    return {
        "id": paper.id,
        "title": paper.title,
        "authors": paper.authors,
        "year": paper.year,
        "doi": paper.doi,
        "journal": paper.journal,
        "url": paper.url or (f"https://doi.org/{paper.doi}" if paper.doi else ""),
        "abstract": paper.abstract,
        "enzyme_family": tag_get(paper.tags, "enzyme_family"),
        "paper_domain": tag_get(paper.tags, "paper_domain"),
        "search_track": tag_get(paper.tags, "search_track"),
        "photoenzymatic_application_status": tag_get(paper.tags, "photoenzymatic_application_status"),
        "new_to_nature_reaction": tag_get(paper.tags, "new_to_nature_reaction"),
        "flavin_cofactor": tag_get(paper.tags, "flavin_cofactor"),
        "characterized_enzyme_evidence": tag_get(paper.tags, "characterized_enzyme_evidence"),
        "candidate_for_photoenzyme_repurposing": tag_get(paper.tags, "candidate_for_photoenzyme_repurposing"),
        "structure_similarity_hint": tag_get(paper.tags, "structure_similarity_hint"),
        "topic_domain": tag_get(paper.tags, "topic_domain"),
        "reaction_type": tag_get(paper.tags, "reaction_type"),
        "mechanism": tag_get(paper.tags, "radical_mechanism"),
        "bond": tag_get(paper.tags, "bond_formation"),
        "substrate": tag_get(paper.tags, "substrate_class"),
        "opportunity": tag_get(paper.tags, "enzyme_integration_opportunity"),
        "hydrolase_relevance": tag_get(paper.tags, "lipase_or_hydrolase_relevance"),
        "review_priority": tag_get(paper.tags, "manual_review_priority"),
        "zotero_status": paper.tags.get("zotero_status", ""),
        "metadata_keywords": tag_get(paper.tags, "metadata_keywords"),
        "relevance_score": paper.tags.get("relevance_score", ""),
        "citation_count": paper.tags.get("citation_count", ""),
        "discovery_depth": paper.tags.get("discovery_depth", ""),
        "criteria_status": paper.tags.get("criteria_status", ""),
        "criteria_score": paper.tags.get("criteria_score", ""),
        "criteria_missing": paper.tags.get("criteria_missing", ""),
        "field_match_status": paper.tags.get("field_match_status", ""),
        "criteria_field_hits": paper.tags.get("criteria_field_hits", ""),
        "criteria_enzyme_feature_hits": paper.tags.get("criteria_enzyme_feature_hits", ""),
        "criteria_reaction_hits": paper.tags.get("criteria_reaction_hits", ""),
        "key_innovation": paper.tags.get("key_innovation", ""),
        "limitation": paper.tags.get("limitation", ""),
        "future_direction": paper.tags.get("future_direction", ""),
    }


def citation_network_to_dashboard(analyzer) -> Dict:
    if not analyzer or not getattr(analyzer, "G", None):
        return {"nodes": [], "edges": []}
    nodes = []
    for node, attr in analyzer.G.nodes(data=True):
        nodes.append({
            "id": node,
            "label": attr.get("title", node)[:80],
            "title": attr.get("title", ""),
            "year": attr.get("year", 0),
            "citations": attr.get("citation_count", 0) or 0,
            "is_seed": bool(attr.get("is_seed")),
            "relevance": attr.get("relevance", 0) or 0,
            "depth": attr.get("depth", 0),
            "venue": attr.get("venue", ""),
        })
    edges = [{"source": u, "target": v} for u, v in analyzer.G.edges()]
    return {"nodes": nodes, "edges": edges}


def write_literature_dashboard(papers: List[PaperNode], output_dir: Path, citation_analyzer=None):
    records = [paper_to_dashboard_record(p) for p in papers]
    network = citation_network_to_dashboard(citation_analyzer)
    data_json = json.dumps({"papers": records, "network": network}, ensure_ascii=False)
    html = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Literature Dashboard</title>
  <style>
    :root {{
      --bg: #f7f7f4; --panel: #ffffff; --ink: #202326; --muted: #667078;
      --line: #d9ddd8; --blue: #356a88; --green: #48735c; --amber: #9a6b2f;
      --red: #9b4b45; --chip: #eef1ed;
    }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: Arial, "Microsoft YaHei", sans-serif; background: var(--bg); color: var(--ink); }}
    header {{ padding: 18px 24px 14px; border-bottom: 1px solid var(--line); background: #fff; position: sticky; top: 0; z-index: 5; }}
    h1 {{ margin: 0 0 10px; font-size: 24px; font-weight: 700; }}
    .controls {{ display: grid; grid-template-columns: 1.5fr repeat(4, minmax(150px, 1fr)); gap: 10px; }}
    input, select {{ width: 100%; border: 1px solid var(--line); background: #fff; padding: 9px 10px; border-radius: 6px; font-size: 14px; }}
    main {{ display: grid; grid-template-columns: minmax(380px, 0.95fr) minmax(440px, 1.05fr); gap: 16px; padding: 16px 24px 24px; }}
    .stats {{ display: flex; gap: 10px; flex-wrap: wrap; margin-top: 10px; color: var(--muted); font-size: 13px; }}
    .stat {{ background: var(--chip); padding: 5px 8px; border-radius: 6px; }}
    .list {{ display: flex; flex-direction: column; gap: 10px; max-height: calc(100vh - 150px); overflow: auto; padding-right: 4px; }}
    .paper {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 12px; cursor: pointer; }}
    .paper.active {{ border-color: var(--blue); box-shadow: 0 0 0 2px rgba(53,106,136,.12); }}
    .title {{ font-weight: 700; line-height: 1.35; }}
    .meta {{ color: var(--muted); font-size: 12px; margin-top: 5px; }}
    .chips {{ display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px; }}
    .chip {{ background: var(--chip); border-radius: 999px; padding: 3px 7px; font-size: 12px; }}
    .chip.high {{ color: #fff; background: var(--green); }}
    .chip.other {{ color: #fff; background: var(--amber); }}
    .chip.bad {{ color: #fff; background: var(--red); }}
    .detail {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; min-height: 320px; padding: 16px; }}
    .detail h2 {{ margin: 0 0 8px; font-size: 20px; }}
    .kv {{ display: grid; grid-template-columns: 150px 1fr; gap: 8px 12px; font-size: 14px; margin-top: 14px; }}
    .kv div:nth-child(odd) {{ color: var(--muted); }}
    .abstract {{ margin-top: 14px; line-height: 1.55; color: #303438; white-space: pre-wrap; }}
    .tabs {{ display: flex; gap: 8px; margin-bottom: 10px; }}
    button {{ border: 1px solid var(--line); background: #fff; border-radius: 6px; padding: 8px 10px; cursor: pointer; }}
    button.active {{ background: var(--blue); color: #fff; border-color: var(--blue); }}
    #network {{ width: 100%; height: 420px; border: 1px solid var(--line); border-radius: 8px; background: #fbfbf9; }}
    .hidden {{ display: none; }}
    @media (max-width: 980px) {{ main {{ grid-template-columns: 1fr; }} .controls {{ grid-template-columns: 1fr 1fr; }} .list {{ max-height: 520px; }} }}
  </style>
</head>
<body>
  <header>
    <h1>Literature Dashboard</h1>
    <div class="controls">
      <input id="q" placeholder="搜索标题、DOI、摘要、关键词">
      <select id="zotero"></select>
      <select id="domain"></select>
      <select id="family"></select>
      <select id="reaction"></select>
      <select id="priority"></select>
      <select id="fieldMatch"></select>
    </div>
    <div class="stats" id="stats"></div>
  </header>
  <main>
    <section>
      <div class="list" id="list"></div>
    </section>
    <section>
      <div class="tabs">
        <button id="tabDetail" class="active">文献详情</button>
        <button id="tabNetwork">引文网络</button>
      </div>
      <div id="detail" class="detail"></div>
      <div id="networkPanel" class="hidden">
        <canvas id="network"></canvas>
        <div class="meta">圆点大小按引用数；蓝色为种子，绿色为发现文献。点击节点会尝试匹配详情。</div>
      </div>
    </section>
  </main>
<script>
const DATA = {data_json};
let papers = DATA.papers || [];
let filtered = papers.slice();
let selected = null;
const $ = id => document.getElementById(id);
function uniq(field) {{
  return [...new Set(papers.map(p => p[field]).filter(v => v && v !== 'N/A'))].sort();
}}
function fillSelect(id, label, field) {{
  const el = $(id);
  el.innerHTML = `<option value="">${{label}}</option>` + uniq(field).map(v => `<option>${{escapeHtml(v)}}</option>`).join('');
}}
function escapeHtml(s) {{ return String(s ?? '').replace(/[&<>"']/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c])); }}
function applyFilters() {{
  const q = $('q').value.toLowerCase().trim();
  const zotero = $('zotero').value, domain = $('domain').value, fam = $('family').value, rxn = $('reaction').value, pri = $('priority').value, fieldMatch = $('fieldMatch').value;
  filtered = papers.filter(p => {{
    const blob = [p.title,p.doi,p.abstract,p.metadata_keywords,p.enzyme_family,p.reaction_type,p.mechanism,p.search_track,p.flavin_cofactor,p.photoenzymatic_application_status,p.characterized_enzyme_evidence,p.structure_similarity_hint,p.criteria_status,p.field_match_status,p.criteria_field_hits,p.criteria_enzyme_feature_hits,p.criteria_reaction_hits].join(' ').toLowerCase();
    return (!q || blob.includes(q)) && (!zotero || p.zotero_status === zotero) && (!domain || p.paper_domain === domain) && (!fam || p.enzyme_family === fam) && (!rxn || p.reaction_type === rxn) && (!pri || p.review_priority === pri) && (!fieldMatch || p.field_match_status === fieldMatch);
  }});
  renderList();
  renderStats();
}}
function renderStats() {{
  const count = filtered.length;
  const other = filtered.filter(p => p.enzyme_family === 'Other_enzyme').length;
  const noenz = filtered.filter(p => p.enzyme_family === 'no enzyme').length;
  const high = filtered.filter(p => p.review_priority === 'high').length;
  const fresh = filtered.filter(p => p.zotero_status === 'new_discovery').length;
  const cpass = filtered.filter(p => p.criteria_status === 'criteria_pass').length;
  const cborder = filtered.filter(p => p.criteria_status === 'criteria_borderline').length;
  $('stats').innerHTML = [`显示 ${{count}} / ${{papers.length}}`, `New ${{fresh}}`, `criteria pass ${{cpass}}`, `borderline ${{cborder}}`, `Other ${{other}}`, `no enzyme ${{noenz}}`, `high priority ${{high}}`].map(x => `<span class="stat">${{x}}</span>`).join('');
}}
function chipClass(p) {{
  if (p.review_priority === 'high') return 'high';
  if (p.enzyme_family === 'Other_enzyme') return 'other';
  if (p.paper_domain === 'unrelated/blacklist' || p.paper_domain === 'organic synthesis no-enzyme') return 'bad';
  return '';
}}
function renderList() {{
  $('list').innerHTML = filtered.map((p, i) => `
    <article class="paper ${{selected && selected.id === p.id ? 'active' : ''}}" onclick="selectPaper(${{i}})">
      <div class="title">${{escapeHtml(p.title)}}</div>
      <div class="meta">${{p.year || ''}} · ${{escapeHtml(p.journal || '')}} · ${{escapeHtml(p.doi || 'no DOI')}}</div>
      <div class="chips">
        <span class="chip ${{chipClass(p)}}">${{escapeHtml(p.enzyme_family)}}</span>
        <span class="chip">${{escapeHtml(p.zotero_status || '')}}</span>
        <span class="chip">${{escapeHtml(p.paper_domain)}}</span>
        <span class="chip">${{escapeHtml(p.criteria_status || '')}}</span>
        <span class="chip">${{escapeHtml(p.field_match_status || '')}}</span>
        <span class="chip">${{escapeHtml(p.search_track || p.flavin_cofactor || '')}}</span>
        <span class="chip">${{escapeHtml(p.reaction_type)}}</span>
        <span class="chip">${{escapeHtml(p.review_priority)}}</span>
      </div>
    </article>`).join('');
  if (!selected && filtered.length) selectPaper(0, false);
}}
function selectPaper(i, rerender=true) {{
  selected = filtered[i];
  if (!selected) return;
  $('detail').innerHTML = `
    <h2>${{escapeHtml(selected.title)}}</h2>
    <div class="meta">${{escapeHtml(selected.authors || '')}} · ${{selected.year || ''}} · ${{escapeHtml(selected.journal || '')}}</div>
    <div class="kv">
      <div>DOI</div><div>${{selected.doi ? `<a href="https://doi.org/${{escapeHtml(selected.doi)}}" target="_blank">${{escapeHtml(selected.doi)}}</a>` : 'missing'}}</div>
      <div>Zotero状态</div><div>${{escapeHtml(selected.zotero_status || '')}}</div>
      <div>文献域</div><div>${{escapeHtml(selected.paper_domain)}}</div>
      <div>Criteria</div><div>${{escapeHtml(selected.criteria_status || '')}} · ${{escapeHtml(selected.criteria_score || '')}} · field=${{escapeHtml(selected.field_match_status || '')}} · missing=${{escapeHtml(selected.criteria_missing || '')}}</div>
      <div>Criteria命中</div><div>${{escapeHtml([selected.criteria_field_hits, selected.criteria_enzyme_feature_hits, selected.criteria_reaction_hits].filter(Boolean).join(' / '))}}</div>
      <div>Flavin任务分轨</div><div>${{escapeHtml(selected.search_track || '')}}</div>
      <div>Flavin辅因子</div><div>${{escapeHtml(selected.flavin_cofactor || '')}}</div>
      <div>Photoenzyme状态</div><div>${{escapeHtml(selected.photoenzymatic_application_status || '')}} · new-to-nature=${{escapeHtml(selected.new_to_nature_reaction || '')}}</div>
      <div>酶家族</div><div>${{escapeHtml(selected.enzyme_family)}}</div>
      <div>反应类型</div><div>${{escapeHtml(selected.reaction_type)}}</div>
      <div>机制</div><div>${{escapeHtml(selected.mechanism)}}</div>
      <div>底物/成键</div><div>${{escapeHtml(selected.substrate)}} / ${{escapeHtml(selected.bond)}}</div>
      <div>表征/候选价值</div><div>${{escapeHtml(selected.characterized_enzyme_evidence || '')}} · ${{escapeHtml(selected.candidate_for_photoenzyme_repurposing || '')}}</div>
      <div>结构相似线索</div><div>${{escapeHtml(selected.structure_similarity_hint || '')}}</div>
      <div>酶口袋机会</div><div>${{escapeHtml(selected.opportunity)}} · ${{escapeHtml(selected.hydrolase_relevance)}}</div>
      <div>API关键词</div><div>${{escapeHtml(selected.metadata_keywords)}}</div>
    </div>
    <div class="abstract">${{escapeHtml(selected.abstract || '')}}</div>`;
  if (rerender) renderList();
}}
function setup() {{
  fillSelect('domain', '全部文献域', 'paper_domain');
  fillSelect('zotero', '全部Zotero状态', 'zotero_status');
  fillSelect('family', '全部酶家族', 'enzyme_family');
  fillSelect('reaction', '全部反应', 'reaction_type');
  fillSelect('priority', '全部优先级', 'review_priority');
  fillSelect('fieldMatch', '全部field命中', 'field_match_status');
  ['q','zotero','domain','family','reaction','priority','fieldMatch'].forEach(id => $(id).addEventListener('input', applyFilters));
  $('tabDetail').onclick = () => showTab('detail');
  $('tabNetwork').onclick = () => showTab('network');
  applyFilters();
}}
function showTab(tab) {{
  $('tabDetail').classList.toggle('active', tab === 'detail');
  $('tabNetwork').classList.toggle('active', tab === 'network');
  $('detail').classList.toggle('hidden', tab !== 'detail');
  $('networkPanel').classList.toggle('hidden', tab !== 'network');
  if (tab === 'network') drawNetwork();
}}
function drawNetwork() {{
  const c = $('network'), ctx = c.getContext('2d');
  const rect = c.getBoundingClientRect(); c.width = rect.width * devicePixelRatio; c.height = rect.height * devicePixelRatio; ctx.scale(devicePixelRatio, devicePixelRatio);
  ctx.clearRect(0,0,rect.width,rect.height);
  const nodes = DATA.network.nodes || [], edges = DATA.network.edges || [];
  if (!nodes.length) {{ ctx.fillText('没有 citation network 数据。运行 --discover 并保留 citation_analyzer。', 20, 30); return; }}
  const pos = new Map();
  nodes.forEach((n,i) => {{
    const angle = 2*Math.PI*i/nodes.length;
    const r = Math.min(rect.width, rect.height)*0.38*(0.65 + 0.35*((n.depth||0)+1)/3);
    pos.set(n.id, [rect.width/2 + Math.cos(angle)*r, rect.height/2 + Math.sin(angle)*r]);
  }});
  ctx.strokeStyle = 'rgba(70,80,80,.18)';
  edges.forEach(e => {{ const a=pos.get(e.source), b=pos.get(e.target); if(!a||!b)return; ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke(); }});
  nodes.forEach(n => {{ const p=pos.get(n.id); const size=Math.max(3, Math.min(14, Math.log((n.citations||0)+1)*2)); ctx.beginPath(); ctx.fillStyle=n.is_seed?'#356a88':'#48735c'; ctx.arc(p[0],p[1],size,0,Math.PI*2); ctx.fill(); }});
}}
setup();
</script>
</body>
</html>"""
    path = output_dir / "literature_dashboard.html"
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"交互式文献看板已保存: {path}")


def export_homolog_candidate_seeds(papers: List[PaperNode], output_dir: Path) -> Path:
    """Write evidence-backed enzyme candidates for downstream homolog/tree curation."""
    rows = []
    for paper in papers:
        track = tag_get(paper.tags, "search_track", default="")
        candidate_rank = tag_get(paper.tags, "candidate_for_photoenzyme_repurposing", default="")
        criteria_status = paper.tags.get("criteria_status", "")
        family = tag_get(paper.tags, "enzyme_family", default="")
        flavin = tag_get(paper.tags, "flavin_cofactor", default="")
        evidence = tag_get(paper.tags, "characterized_enzyme_evidence", default="")
        structure_hint = tag_get(paper.tags, "structure_similarity_hint", default="")
        if family in {"", "N/A", "no enzyme", "Other_enzyme", "unclear"}:
            continue
        keep = (
            criteria_status in {"criteria_pass", "criteria_borderline"}
            or track == "B_characterized_flavin_enzyme_not_photoenzymatic"
            or candidate_rank in {"high", "medium"}
        )
        has_flavin_or_evidence = (
            flavin in {"FAD", "FMN", "flavin unspecified"}
            or "flavin" in family.lower()
            or evidence in {"purified enzyme kinetics", "activity assay only", "structure only"}
            or structure_hint in {"PDB available", "same fold/cofactor pocket", "similar substrate pocket"}
        )
        if not keep or not has_flavin_or_evidence:
            continue
        rows.append({
            "title": paper.title,
            "doi": normalize_doi_for_filter(paper.doi),
            "year": paper.year,
            "journal": paper.journal,
            "zotero_status": paper.tags.get("zotero_status", ""),
            "criteria_status": criteria_status,
            "criteria_score": paper.tags.get("criteria_score", ""),
            "criteria_missing": paper.tags.get("criteria_missing", ""),
            "search_track": track,
            "enzyme_family": family,
            "enzyme_name_or_target": tag_get(paper.tags, "enzyme_name_or_target", default=""),
            "flavin_cofactor": flavin,
            "reaction_type": tag_get(paper.tags, "reaction_type", default=""),
            "new_to_nature_reaction": tag_get(paper.tags, "new_to_nature_reaction", default=""),
            "characterized_enzyme_evidence": evidence,
            "candidate_for_photoenzyme_repurposing": candidate_rank,
            "structure_similarity_hint": structure_hint,
            "mechanistic_evidence": tag_get(paper.tags, "mechanistic_evidence", default=""),
            "manual_review_priority": tag_get(paper.tags, "manual_review_priority", default=""),
            "metadata_keywords": tag_get(paper.tags, "metadata_keywords", default=""),
            "abstract": paper.abstract,
        })
    rows.sort(
        key=lambda row: (
            row["criteria_status"] == "criteria_pass",
            row["candidate_for_photoenzyme_repurposing"] == "high",
            row["manual_review_priority"] == "high",
            row["year"] or 0,
        ),
        reverse=True,
    )
    path = output_dir / "homolog_candidate_seeds.csv"
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")
    print(f"[Criteria] Homolog candidate seed table saved: {path} ({len(rows)} rows)")
    return path


def extract_pdb_ids(text: str) -> List[str]:
    ids = []
    for match in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", text or ""):
        if any(c.isalpha() for c in match):
            ids.append(match.upper())
    return dedupe_preserve_order(ids)


def extract_uniprot_accessions(text: str) -> List[str]:
    pattern = r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:-[0-9]+)?|[A-Z0-9]{10})\b"
    return dedupe_preserve_order(re.findall(pattern, text or ""))


def clean_enzyme_candidate_name(name: str) -> str:
    value = re.sub(r"\s+", " ", str(name or "")).strip(" .,:;()[]")
    value = re.sub(r"^(engineered|direct|visible-light|photoexcited|mutant|variant)\s+", "", value, flags=re.I)
    value = re.sub(r"\s+(variants?|mutants?)$", "", value, flags=re.I)
    false_names = {
        "engineered", "centered", "direct", "visible", "protein", "enzyme",
        "flavin", "cofactor", "radical",
    }
    if not value or value.lower() in false_names:
        return ""
    if len(value) < 4:
        return ""
    enzyme_signal = any(
        token in value.lower()
        for token in [
            "ase", "enzyme", "reductase", "oxidase", "halogenase", "monooxygenase",
            "dehydrogenase", "photodecarboxylase", "photolyase", "cryptochrome",
            "nitroreductase", "fap", "fdh", "ered", "oye", "fmo", "bbe",
        ]
    )
    if not enzyme_signal:
        return ""
    return value


def extract_enzyme_name_candidates(paper: PaperNode) -> List[str]:
    names = []
    explicit = tag_get(paper.tags, "enzyme_name_or_target", default="")
    if explicit and explicit not in {"unclear", "N/A"}:
        for part in re.split(r"[|;,/]+", explicit):
            part = clean_enzyme_candidate_name(part)
            if part and part.lower() not in {"unclear", "unknown", "n/a"}:
                names.append(part)
    text = f"{paper.title}. {paper.abstract}"
    phrase_patterns = [
        r"\b(?:engineered\s+)?(?:flavin-dependent|flavin dependent)\s+[A-Za-z -]{3,45}ases?\b",
        r"\b(?:fatty acid photodecarboxylase|photodecarboxylase|ene-reductase|ene reductase|old yellow enzyme|flavin-dependent halogenase|flavin dependent halogenase|photolyase|cryptochrome|nitroreductase|berberine bridge enzyme|flavin monooxygenase|flavin reductase|glucose oxidase|D-amino acid oxidase)\b",
        r"\b(?:[A-Z][A-Za-z0-9]{0,8})?(?:FAP|FDH|ERED|OYE)\b",
    ]
    for pattern in phrase_patterns:
        flags = 0 if "FAP|FDH|ERED|OYE" in pattern else re.I
        for match in re.findall(pattern, text, flags=flags):
            cleaned = clean_enzyme_candidate_name(match)
            if cleaned:
                names.append(cleaned)
    return dedupe_preserve_order(names)[:8]


def make_uniprot_query(name: str, family: str, paper: PaperNode) -> str:
    terms = []
    for value in [name, family]:
        value = str(value or "").strip()
        if value and value not in {"Other_enzyme", "unclear", "no enzyme", "N/A"}:
            terms.append(value)
    if paper.title:
        source_match = re.search(r"\b(?:from|of|in)\s+([A-Z][a-z]+(?:\s+[a-z]+){0,2})", paper.title)
        if source_match:
            terms.append(source_match.group(1))
    base = " ".join(dedupe_preserve_order(terms))
    if not base:
        base = paper.title[:80]
    return f'({base}) AND (flavin OR FAD OR FMN OR flavoenzyme)'


def fetch_uniprot_candidates(query: str, size: int = 3, timeout: float = 25.0) -> List[Dict[str, str]]:
    params = {
        "query": query,
        "format": "json",
        "size": size,
        "fields": "accession,id,protein_name,organism_name,sequence,length,xref_pdb,reviewed",
    }
    try:
        resp = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"[UniProt] Query failed for {query[:80]}: {e}")
        return []
    results = []
    for item in data.get("results", []):
        seq_obj = item.get("sequence") or {}
        protein = item.get("proteinDescription", {})
        recommended = protein.get("recommendedName", {}).get("fullName", {}).get("value", "")
        if not recommended:
            recommended = protein.get("submissionNames", [{}])[0].get("fullName", {}).get("value", "")
        pdb_ids = []
        for xref in item.get("uniProtKBCrossReferences", []) or []:
            if xref.get("database") == "PDB" and xref.get("id"):
                pdb_ids.append(xref["id"])
        results.append({
            "uniprot_id": item.get("primaryAccession", ""),
            "uniprot_entry": item.get("uniProtkbId", ""),
            "protein_name": recommended,
            "organism": (item.get("organism") or {}).get("scientificName", ""),
            "reviewed": str(item.get("entryType", "")),
            "sequence_length": seq_obj.get("length", ""),
            "sequence": seq_obj.get("value", ""),
            "pdb_ids": " | ".join(dedupe_preserve_order(pdb_ids)),
        })
    return results


def export_enzyme_seed_candidates(
    papers: List[PaperNode],
    output_dir: Path,
    fetch_uniprot: bool = False,
    max_seeds: int = 100,
    uniprot_per_seed: int = 2,
) -> Tuple[Path, Optional[Path]]:
    rows = []
    fasta_records = []
    for paper in sorted(papers, key=rank_paper_for_quota, reverse=True):
        family = tag_get(paper.tags, "enzyme_family", default="")
        if family in {"", "N/A", "no enzyme", "Other_enzyme", "unclear"}:
            continue
        flavin = tag_get(paper.tags, "flavin_cofactor", default="")
        track = tag_get(paper.tags, "search_track", default="")
        candidate_rank = tag_get(paper.tags, "candidate_for_photoenzyme_repurposing", default="")
        evidence = tag_get(paper.tags, "characterized_enzyme_evidence", default="")
        text = f"{paper.title} {paper.abstract} {paper.tags.get('metadata_keywords', '')}"
        if not (
            flavin in {"FAD", "FMN", "flavin unspecified"}
            or "flavin" in text.lower()
            or track in {"A_used_flavin_photoenzyme_new_to_nature", "B_characterized_flavin_enzyme_not_photoenzymatic"}
            or candidate_rank in {"high", "medium"}
        ):
            continue
        names = extract_enzyme_name_candidates(paper) or [family]
        pdb_ids = extract_pdb_ids(text)
        accessions = extract_uniprot_accessions(text)
        for name in names:
            base_row = {
                "seed_priority": candidate_rank or tag_get(paper.tags, "manual_review_priority", default=""),
                "enzyme_name": name,
                "enzyme_family": family,
                "uniprot_id": " | ".join(accessions),
                "pdb_ids": " | ".join(pdb_ids),
                "doi": normalize_doi_for_filter(paper.doi),
                "title": paper.title,
                "year": paper.year,
                "journal": paper.journal,
                "search_track": track,
                "flavin_cofactor": flavin,
                "reaction_type": tag_get(paper.tags, "reaction_type", default=""),
                "new_to_nature_reaction": tag_get(paper.tags, "new_to_nature_reaction", default=""),
                "characterized_enzyme_evidence": evidence,
                "structure_similarity_hint": tag_get(paper.tags, "structure_similarity_hint", default=""),
                "criteria_status": paper.tags.get("criteria_status", ""),
                "manual_review_priority": tag_get(paper.tags, "manual_review_priority", default=""),
                "evidence_summary": paper.abstract[:700],
                "uniprot_query": make_uniprot_query(name, family, paper),
                "sequence_source": "",
                "sequence": "",
            }
            rows.append(base_row)
            if len(rows) >= max_seeds and not fetch_uniprot:
                break
        if len(rows) >= max_seeds and not fetch_uniprot:
            break

    if fetch_uniprot:
        expanded_rows = []
        for row in rows[:max_seeds]:
            hits = fetch_uniprot_candidates(row["uniprot_query"], size=uniprot_per_seed)
            if not hits:
                expanded_rows.append(row)
                continue
            for hit in hits:
                merged = dict(row)
                merged.update({
                    "uniprot_id": hit["uniprot_id"],
                    "uniprot_entry": hit["uniprot_entry"],
                    "uniprot_protein_name": hit["protein_name"],
                    "organism": hit["organism"],
                    "reviewed": hit["reviewed"],
                    "sequence_length": hit["sequence_length"],
                    "pdb_ids": hit["pdb_ids"] or row["pdb_ids"],
                    "sequence_source": "UniProt REST search",
                    "sequence": hit["sequence"],
                })
                expanded_rows.append(merged)
                if hit["sequence"]:
                    header = "|".join(
                        [
                            hit["uniprot_id"],
                            re.sub(r"\s+", "_", hit["protein_name"] or row["enzyme_name"])[:60],
                            f"doi={row['doi'] or 'NA'}",
                            f"family={row['enzyme_family'].replace(' ', '_')}",
                        ]
                    )
                    fasta_records.append(f">{header}\n{hit['sequence']}")
        rows = expanded_rows

    rows = rows[: max_seeds if not fetch_uniprot else max(max_seeds, len(rows))]
    csv_path = output_dir / "enzyme_seed_candidates.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False, encoding="utf-8-sig")
    fasta_path = None
    if fasta_records:
        fasta_path = output_dir / "enzyme_seed_candidates.fasta"
        fasta_path.write_text("\n".join(fasta_records) + "\n", encoding="utf-8")
    print(f"[Seeds] Enzyme seed candidate table saved: {csv_path} ({len(rows)} rows)")
    if fasta_path:
        print(f"[Seeds] UniProt FASTA saved: {fasta_path} ({len(fasta_records)} sequences)")
    elif fetch_uniprot:
        print("[Seeds] UniProt FASTA not written because no sequences were resolved")
    return csv_path, fasta_path


def export_sequence_seed_review(output_dir: Path) -> Optional[Path]:
    """Run the offline sequence seed resolver when candidate CSVs are available."""
    script_path = Path(__file__).resolve().parent / "scripts" / "sequence_seed_resolver.py"
    enzyme_csv = output_dir / "enzyme_seed_candidates.csv"
    homolog_csv = output_dir / "homolog_candidate_seeds.csv"
    if not script_path.exists() or not enzyme_csv.exists():
        return None
    import subprocess
    cmd = [
        sys.executable,
        str(script_path),
        "--enzyme-csv",
        str(enzyme_csv),
        "--outdir",
        str(output_dir),
        "--prefix",
        "sequence_seed",
    ]
    if homolog_csv.exists():
        cmd.extend(["--homolog-csv", str(homolog_csv)])
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"[SequenceSeed] Review export failed: {e}")
        return None
    return output_dir / "sequence_seed_review.html"


def build_extraction_context(paper: PaperNode) -> str:
    abstract = (paper.abstract or "").strip()
    if len(abstract) > 1600:
        abstract = abstract[:1600].rsplit(" ", 1)[0] + " ..."
    parts = [abstract]
    keywords = paper.tags.get("metadata_keywords", "")
    if keywords:
        if len(keywords) > 500:
            keywords = keywords[:500]
        parts.append(f"API keywords/subjects/concepts: {keywords}")
    if paper.journal:
        parts.append(f"Venue/journal: {paper.journal}")
    return "\n".join(part for part in parts if part)


@dataclass
class PaperNode:
    """知识图谱中的文献节点"""
    id: str
    title: str
    authors: str = ""
    year: int = 0
    doi: str = ""
    journal: str = ""
    abstract: str = ""
    url: str = ""
    file_attachments: str = ""
    link_attachments: str = ""
    # LLM提取的语义标签
    tags: Dict[str, str] = field(default_factory=dict)
    # 衍生属性
    relevance_score: float = 0.0  # 基于引用数和期刊影响力的综合得分
    
    def get_display_title(self, max_len: int = 50) -> str:
        t = self.title.strip()
        return t if len(t) <= max_len else t[:max_len-3] + "..."


def load_schema(schema_profile: str = "photo-radical", schema_file: Optional[str] = None) -> Dict[str, str]:
    """加载内置schema或外部JSON schema文件。"""
    if schema_file:
        with open(schema_file, "r", encoding="utf-8") as f:
            custom_schema = json.load(f)
        if not isinstance(custom_schema, dict):
            raise ValueError("schema文件必须是JSON对象，例如 {\"field\": \"description\"}")
        return {str(k): str(v) for k, v in custom_schema.items()}
    if schema_profile not in SCHEMA_PROFILES:
        raise ValueError(f"未知schema profile: {schema_profile}")
    return SCHEMA_PROFILES[schema_profile]


def read_task_prompt(task_prompt: Optional[str] = None, task_file: Optional[str] = None) -> str:
    """Read an optional high-level research objective for LLM extraction."""
    parts = []
    if task_prompt:
        parts.append(str(task_prompt).strip())
    if task_file:
        path = Path(task_file)
        parts.append(path.read_text(encoding="utf-8").strip())
    return "\n\n".join(part for part in parts if part)


def write_task_spec(output_dir: Path, schema_profile: str, task_prompt: str, search_queries: Optional[List[str]] = None) -> None:
    if not task_prompt:
        return
    lines = [
        "# Literature Agent Task Spec",
        "",
        f"- Schema profile: `{schema_profile}`",
        "",
        "## User Research Objective",
        "",
        task_prompt,
        "",
        "## Operational Interpretation",
        "",
        "- First collect papers where flavin/FAD/FMN-dependent enzymes have already been used for photoenzymatic or light-driven new-to-nature reactions.",
        "- Then collect characterized flavin enzymes with enzyme assays, structures, kinetics, or substrate scope that have not yet been applied in recent photoenzymatic synthesis.",
        "- Deprioritize genome-only annotations and papers with no enzyme evidence.",
        "- Treat DOI/PMID/title metadata as discovery seeds; use citation/reference expansion and explicit keyword searches to find missed landmarks.",
    ]
    if search_queries:
        lines.extend(["", "## Search Queries", ""])
        lines.extend(f"- {q}" for q in search_queries)
    (output_dir / "task_spec.md").write_text("\n".join(lines), encoding="utf-8")


def write_criteria_spec(output_dir: Path, criteria: ResearchCriteria, search_queries: Optional[List[str]] = None) -> None:
    if not criteria.enabled():
        return
    lines = [
        "# Literature Agent Criteria Spec",
        "",
        "## Hard Screening Anchors",
        "",
        f"- Field/domain: {criteria.field or 'not set'}",
        f"- Enzyme feature/cofactor: {criteria.enzyme_feature or 'not set'}",
        f"- Reaction: {criteria.reaction or 'not set'}",
        "",
        "## Final Goal",
        "",
        criteria.goal or "not set",
        "",
        "## Expanded Matching Terms",
        "",
    ]
    for label, raw_value in criteria.hard_fields():
        if raw_value:
            terms = expand_criteria_terms(raw_value, label)
            lines.append(f"### {label}")
            lines.append("")
            lines.append(", ".join(terms[:80]) if terms else "No terms generated.")
            lines.append("")
    if search_queries:
        lines.extend(["## Search Queries", ""])
        lines.extend(f"- {q}" for q in search_queries)
    (output_dir / "criteria_spec.md").write_text("\n".join(lines), encoding="utf-8")


def default_search_queries_for_task(schema_profile: str, task_prompt: str) -> List[str]:
    text = f"{schema_profile} {task_prompt}".lower()
    if "flavin" not in text and "fad" not in text and "fmn" not in text and "yellow enzyme" not in text:
        return []
    if schema_profile == "flavin-photoenzyme" and not any(
        k in text for k in ["photoenzym", "photoenzyme", "photobiocatal", "new-to-nature", "new to nature"]
    ):
        return [
            "\"flavin-dependent\" enzyme \"substrate scope\"",
            "\"FAD\" enzyme substrate scope kinetics structure",
            "\"FMN\" enzyme substrate scope kinetics structure",
            "\"flavoenzyme\" kinetics structure substrate",
            "\"old yellow enzyme\" kinetics structure",
            "\"ene-reductase\" enzyme substrate scope",
        ]
    return [
        "\"flavin\" photoenzymatic \"new-to-nature\"",
        "\"FAD\" \"photoenzymatic\" radical enzyme",
        "\"FMN\" \"photoenzymatic\" biocatalysis",
        "\"flavoenzyme\" \"visible light\" biocatalysis",
        "\"old yellow enzyme\" photoenzymatic radical",
        "\"ene-reductase\" \"photoenzymatic\"",
        "\"flavin-dependent\" enzyme \"substrate scope\"",
        "\"flavoenzyme\" kinetics structure substrate",
    ]


def read_zotero_csv(csv_path: str, encoding: str = "auto") -> Tuple[pd.DataFrame, str]:
    """读取Zotero CSV，兼容UTF-8/UTF-16/GB18030等导出。"""
    encodings = [encoding] if encoding != "auto" else ["utf-8-sig", "utf-8", "utf-16", "utf-16-le", "gb18030"]
    errors = []
    for enc in encodings:
        try:
            df = pd.read_csv(csv_path, encoding=enc)
            return df, enc
        except UnicodeError as e:
            errors.append(f"{enc}: {e}")
        except pd.errors.ParserError as e:
            errors.append(f"{enc}: {e}")
    raise ValueError("无法读取CSV编码，请用 --csv-encoding 指定。尝试结果:\n" + "\n".join(errors))


def configure_chinese_font():
    """尽量为matplotlib配置可用中文字体，避免图表中文变方块。"""
    font_files = [
        r"C:\Windows\Fonts\msyh.ttc",
        r"C:\Windows\Fonts\simhei.ttf",
        r"C:\Windows\Fonts\simsun.ttc",
        "/mnt/c/Windows/Fonts/msyh.ttc",
        "/mnt/c/Windows/Fonts/simhei.ttf",
        "/mnt/c/Windows/Fonts/simsun.ttc",
    ]
    for font_file in font_files:
        if Path(font_file).exists():
            fm.fontManager.addfont(font_file)
            prop = fm.FontProperties(fname=font_file)
            name = prop.get_name()
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return name

    candidates = [
        "Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC",
        "Source Han Sans SC", "WenQuanYi Micro Hei", "Arial Unicode MS"
    ]
    installed = {f.name for f in fm.fontManager.ttflist}
    for name in candidates:
        if name in installed:
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return name
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False
    return None


def normalize_api_key(api_key: str) -> str:
    """清洗并校验API key，避免requests在HTTP header编码阶段报latin-1错误。"""
    key = (api_key or "").strip().strip("'\"")
    if not key:
        return ""
    try:
        key.encode("ascii")
    except UnicodeEncodeError:
        raise ValueError(
            "API key中包含非ASCII字符。请检查 export KIMI_API_KEY=... 是否仍是示例占位符，"
            "例如不要使用“你的key”，也不要把中文说明一起复制进去。"
        )
    if key in {"你的key", "YOUR_KIMI_API_KEY", "sk-你的key"}:
        raise ValueError("API key仍是占位符，请替换成Kimi控制台中的真实key。")
    if not key.startswith("sk-"):
        print("警告: API key不是常见的 sk- 开头，请确认是否复制完整。")
    return key


def mask_api_key(api_key: str) -> str:
    key = api_key or ""
    if len(key) <= 10:
        return "*" * len(key)
    return f"{key[:6]}...{key[-4:]}"


def get_api_key_source(args_api_key: Optional[str]) -> str:
    if args_api_key:
        return "--api-key"
    if os.getenv("KIMI_API_KEY") is not None:
        return "KIMI_API_KEY"
    if os.getenv("MOONSHOT_API_KEY") is not None:
        return "MOONSHOT_API_KEY"
    return "none"


def print_auth_diagnostic(raw_key: str, source: str):
    key = raw_key or ""
    print(
        "API Key预检: "
        f"source={source}, len={len(key)}, ascii={key.isascii()}, "
        f"masked={mask_api_key(key.strip().strip(chr(39) + chr(34)))}"
    )
    if key and not key.isascii():
        bad_chars = sorted({c for c in key if not c.isascii()})
        preview = " ".join(f"U+{ord(c):04X}" for c in bad_chars[:6])
        print(f"API Key预检: 检测到非ASCII字符 {preview}")


class LLMExtractor:
    """调用LLM API进行文献语义提取"""
    
    def __init__(self, config: Dict, schema: Dict[str, str], task_prompt: str = ""):
        self.config = config
        self.config["api_key"] = normalize_api_key(config.get("api_key", ""))
        self.schema = schema
        self.task_prompt = (task_prompt or "").strip()
        self.headers = {
            "Authorization": f"Bearer {self.config['api_key']}",
            "Content-Type": "application/json"
        }

    @staticmethod
    def _parse_json_content(content: str) -> Dict[str, str]:
        content = (content or "").strip()
        if content.startswith("```"):
            content = re.sub(r"^```(?:json)?\s*", "", content, flags=re.I)
            content = re.sub(r"\s*```$", "", content)
        try:
            return json.loads(content)
        except json.JSONDecodeError:
            match = re.search(r"\{.*\}", content, flags=re.S)
            if match:
                return json.loads(match.group(0))
            raise

    def _chat_json(self, prompt: str, max_tokens: Optional[int] = None) -> Dict[str, str]:
        base_payload = {
            "model": self.config['model'],
            "messages": [{"role": "user", "content": prompt}],
            "temperature": self.config['temperature'],
        }
        model_name = str(self.config.get("model", "")).lower()
        if model_name.startswith("kimi-k2") and not self.config.get("enable_thinking", False):
            base_payload["thinking"] = {"type": "disabled"}

        requested_tokens = int(max_tokens or self.config["max_tokens"])
        token_budgets = [requested_tokens]

        attempts = []
        json_modes = [self.config.get("json_mode", True)]
        if self.config.get("json_mode", True):
            json_modes.append(False)
        for token_budget in token_budgets:
            for token_field in ("max_tokens",):
                for use_json_mode in json_modes:
                    payload = dict(base_payload)
                    payload[token_field] = token_budget
                    if use_json_mode:
                        payload["response_format"] = {"type": "json_object"}
                    attempts.append((token_field, token_budget, use_json_mode, payload))

        last_error = None
        for token_field, token_budget, use_json_mode, payload in attempts:
            resp = requests.post(
                f"{self.config['base_url']}/chat/completions",
                headers=self.headers,
                json=payload,
                timeout=(
                    self.config.get("connect_timeout", 15),
                    self.config.get("read_timeout", 360),
                )
            )
            if not resp.ok:
                last_error = f"{resp.status_code} token_field={token_field} tokens={token_budget} json_mode={use_json_mode}: {resp.text[:300]}"
                if resp.status_code in (401, 403):
                    break
                continue

            data = resp.json()
            choice = (data.get("choices") or [{}])[0]
            message = choice.get("message") or {}
            content = message.get("content") or ""
            if not content.strip():
                reasoning_len = len(message.get("reasoning_content") or "")
                last_error = (
                    f"empty message.content token_field={token_field} tokens={token_budget} "
                    f"json_mode={use_json_mode} finish_reason={choice.get('finish_reason')} "
                    f"reasoning_content_chars={reasoning_len}"
                )
                continue
            try:
                return self._parse_json_content(content)
            except Exception as e:
                last_error = (
                    f"JSON parse failed token_field={token_field} tokens={token_budget} "
                    f"json_mode={use_json_mode}: {e}; content={content[:240]}"
                )
                continue

        raise ValueError(last_error or "no parseable JSON response")
        
    def extract(self, title: str, abstract: str, journal: str = "") -> Dict[str, str]:
        """对单篇文献进行结构化提取"""
        fallback_tags = fallback_extract(title, abstract, self.schema)
        if not abstract or len(abstract) < 50:
            # 摘要太短或缺失时，至少用标题/期刊关键词给出可聚类的粗标签。
            return fallback_tags
            
        schema_desc = "\n".join([f"- {k}: {v}" for k, v in self.schema.items()])
        task_block = ""
        if self.task_prompt:
            task_block = f"""

用户当前研究目标（最高优先级，请据此解释长难句任务并服务于检索/分类）:
{self.task_prompt}

请特别区分：
- A类：已经把flavin/FAD/FMN相关酶用于photoenzymatic或光驱动new-to-nature合成反应的报道。
- B类：已有酶学表征、结构、动力学或底物谱证据的flavin酶，但尚未被近年photoenzymatic合成领域应用的候选酶。
- C类：只是基因组注释、无实验酶学证据，除非非常关键，否则人工复核优先级应降低。
"""
        
        prompt = f"""你是一位photoenzymatic catalysis、radical chemistry、organic photochemistry和酶工程交叉领域的文献分析专家。请仔细阅读以下文献信息，并按schema提取结构化知识。

分析优先级：
1. 先判断文章属于photoenzymatic radical、photoredox radical reaction还是更广义的radical chemical reaction。
2. 如果尚未接入酶体系，不要强行归为酶催化；请标记为“仅光化学/自由基反应”或“可设想接入酶体系”。
3. lipase/alpha-beta hydrolase只作为潜在关联或人工筛选线索，不作为主分类前提。
4. enzyme_family必须尽量查清具体家族；不能把普通英文词或泛泛的Other_enzyme当作结果。若文献没有酶体系，请填"no enzyme"；若无法澄清且不适合注释，请is_blacklisted=true并说明原因。
{task_block}

文献标题: {title}
期刊: {journal}
摘要: {abstract}

提取Schema（请严格按JSON格式返回，只返回JSON，不要markdown代码块）:
{schema_desc}

注意：
1. 分类字段请优先根据标题和摘要做最合理的粗粒度判断；只有完全没有线索时才填"unclear"，尽量不要填"N/A"
2. key_innovation/limitation/future_direction必须用中文一句话概括
3. 多值字段用" | "分隔，如"C-C | C-N"
4. is_blacklisted必须返回true或false；blacklist_reason在false时可填空字符串
5. 必须返回合法JSON，不要任何额外解释"""

        try:
            result = self._chat_json(prompt)
            for k in self.schema.keys():
                if k not in result:
                    result[k] = "unclear" if k in CLASSIFICATION_FIELDS else "N/A"
            return merge_with_fallback(result, fallback_tags, self.schema)
            if False:
                attempts = []
            token_fields = ["max_completion_tokens", "max_tokens"]
            json_modes = [self.config.get("json_mode", True)]
            if self.config.get("json_mode", True):
                json_modes.append(False)
            for token_field in token_fields:
                for use_json_mode in json_modes:
                    payload = dict(base_payload)
                    payload[token_field] = self.config["max_tokens"]
                    if use_json_mode:
                        payload["response_format"] = {"type": "json_object"}
                    attempts.append((token_field, use_json_mode, payload))

            last_error = None
            for token_field, use_json_mode, payload in attempts:
                resp = requests.post(
                    f"{self.config['base_url']}/chat/completions",
                    headers=self.headers,
                    json=payload,
                    timeout=(
                        self.config.get("connect_timeout", 15),
                        self.config.get("read_timeout", 360),
                    )
                )
                if resp.ok:
                    break
                last_error = f"{resp.status_code} token_field={token_field} json_mode={use_json_mode}: {resp.text[:300]}"
                if resp.status_code in (401, 403):
                    break
            else:
                resp.raise_for_status()

            if not resp.ok:
                print(f"  LLM请求失败详情: {last_error}")
                resp.raise_for_status()

            resp.raise_for_status()
            data = resp.json()
            content = data['choices'][0]['message']['content']
            result = json.loads(content)
            # 确保所有key存在
            for k in self.schema.keys():
                if k not in result:
                    result[k] = "unclear" if k in CLASSIFICATION_FIELDS else "N/A"
            return merge_with_fallback(result, fallback_tags, self.schema)
            
        except requests.exceptions.ReadTimeout as e:
            print(
                f"  LLM提取超时: 等待模型输出超过 {self.config.get('read_timeout', 360)} 秒。"
                "可调高 --read-timeout，或降低 --max-tokens/缩短摘要输入。"
            )
            return fallback_tags
        except requests.exceptions.ConnectTimeout as e:
            print(
                f"  LLM连接超时: 连接API超过 {self.config.get('connect_timeout', 15)} 秒。"
                "请检查网络、base-url或代理。"
            )
            return fallback_tags
        except Exception as e:
            detail = ""
            if "resp" in locals():
                try:
                    detail = f" | 响应: {resp.text[:500]}"
                except Exception:
                    pass
            print(f"  LLM提取失败: {e}{detail}")
            return fallback_tags

    def resolve_other_enzyme(self, paper: "PaperNode") -> Dict[str, str]:
        """Second-pass LLM check for Other_enzyme/unclear enzyme family."""
        fallback_tags = fallback_extract(paper.title, paper.abstract, self.schema)
        introduction = paper.tags.get("introduction_snippet", "")
        prompt = f"""你是一位酶催化、光酶催化和有机光氧化还原文献筛选专家。请只解决 enzyme_family 归类问题，并判断这篇文献是否适合保留在可注释文献集中。

任务：
1. 彻底判断文献是否真的涉及酶/蛋白/生物催化体系。
2. 如果涉及，请给出具体 enzyme_family，例如 P450、ERED/OYE、FAP、flavoprotein、nonheme iron enzyme、lipase/esterase、alpha/beta-hydrolase、transaminase、dehydrogenase、halogenase、oxidase、reductase 等。
3. 如果不涉及酶，请填 "no enzyme"，不要填 Other_enzyme。
4. 如果只有泛泛的 enzyme、biocatalysis，但无法从标题/摘要判断具体家族，请填 "Other_enzyme"，并解释为什么仍值得保留或为什么应拉黑。
5. release/base/showcase/disease/phase 等普通英文词不是酶名，不能作为 enzyme_family。
6. 对无 DOI、非化学/非酶、无法注释、或 Other_enzyme 且无法澄清的文献，is_blacklisted=true。

标题: {paper.title}
期刊: {paper.journal}
年份: {paper.year}
DOI: {paper.doi}
API关键词/学科/概念: {paper.tags.get('metadata_keywords', '')}
摘要: {paper.abstract}
Introduction片段（若抓取成功）: {introduction}

请只返回合法JSON，字段如下：
{{
  "enzyme_family": "...",
  "paper_domain": "photoenzymatic/enzymatic|enzymatic non-photo|organic photoredox/radical no-enzyme|organic synthesis no-enzyme|unrelated/blacklist",
  "photoenzymatic_status": "...",
  "enzyme_integration_opportunity": "high|medium|low|unclear",
  "is_blacklisted": "true|false",
  "blacklist_reason": "..."
}}"""
        try:
            result = self._chat_json(prompt, max_tokens=min(1200, self.config["max_tokens"]))
            result["enzyme_family"] = normalize_enzyme_family(
                result.get("enzyme_family", ""),
                paper.title,
                paper.abstract,
            )
            return merge_with_fallback(result, fallback_tags, self.schema)
        except Exception as e:
            print(f"  Other_enzyme二次精查失败: {e}")
            return fallback_tags

    def resolve_flavin_enzyme_family(self, paper: "PaperNode") -> Dict[str, str]:
        """Focused low-token resolver for flavin enzyme family and seed usefulness."""
        fallback_tags = fallback_extract(paper.title, paper.abstract, self.schema)
        introduction = paper.tags.get("introduction_snippet", "")
        prompt = f"""You are resolving enzyme identity for a flavin/photoenzymatic literature triage table.

Return only JSON. Do not blacklist papers just because the enzyme family is unclear.

Rules:
1. If the paper mentions flavin, FAD, FMN, flavoenzyme, OYE, ERED, FAP, FDH, photolyase, cryptochrome, nitroreductase, oxidase, reductase, dehydrogenase, monooxygenase, halogenase, or berberine bridge enzyme, assume an enzyme/protein is plausible unless the abstract clearly says it is small-molecule flavin photocatalysis with no enzyme.
2. Give the narrowest enzyme_family supported by title/abstract: ERED/OYE, flavin-dependent halogenase, FAP, photolyase/cryptochrome, flavin monooxygenase, flavin reductase/dehydrogenase, nitroreductase, berberine bridge enzyme, flavoenzyme oxidase, dehydrogenase, oxidase, reductase, or other specific family.
3. Use "Other_enzyme" only when no narrower family is supportable. Use "no enzyme" only when clearly no protein/enzyme system is involved.
4. candidate_for_photoenzyme_repurposing should be high/medium for characterized flavin enzymes with structure, PDB, purified enzyme, kinetics, activity assay, substrate scope, or mutagenesis evidence, especially if not already a photoenzymatic new-to-nature application.
5. is_blacklisted must be false unless the paper is clearly unrelated or has no useful chemical/enzyme evidence.

Title: {paper.title}
Journal: {paper.journal}
Year: {paper.year}
DOI: {paper.doi}
API keywords/subjects/concepts: {paper.tags.get('metadata_keywords', '')}
Current tags: enzyme_family={tag_get(paper.tags, 'enzyme_family')}; flavin={tag_get(paper.tags, 'flavin_cofactor', default='')}; track={tag_get(paper.tags, 'search_track', default='')}
Abstract: {paper.abstract}
Introduction snippet: {introduction}

JSON schema:
{{
  "enzyme_family": "...",
  "enzyme_name_or_target": "specific enzyme/protein name, accession/PDB clue, or unclear",
  "flavin_cofactor": "FAD|FMN|flavin unspecified|no flavin|unclear",
  "paper_domain": "flavin photoenzymatic application|characterized flavin enzyme non-photo|photoenzymatic non-flavin|organic photoredox/radical no-enzyme|unrelated/blacklist",
  "photoenzymatic_application_status": "applied to new-to-nature reaction|native/known enzymology only|candidate not yet photoenzymatic|no enzyme|unclear",
  "new_to_nature_reaction": "true|false|unclear",
  "characterized_enzyme_evidence": "purified enzyme kinetics|activity assay only|structure only|genome annotation only|no enzyme evidence|unclear",
  "candidate_for_photoenzyme_repurposing": "high|medium|low|no",
  "structure_similarity_hint": "same fold/cofactor pocket|homolog of applied enzyme|similar substrate pocket|PDB available|only sequence annotation|unclear",
  "manual_review_priority": "high|medium|low",
  "is_blacklisted": "true|false",
  "blacklist_reason": ""
}}"""
        try:
            result = self._chat_json(prompt, max_tokens=min(1600, self.config["max_tokens"]))
            result["enzyme_family"] = normalize_enzyme_family(
                result.get("enzyme_family", ""),
                paper.title,
                paper.abstract,
            )
            if str(result.get("is_blacklisted", "")).lower() == "true" and is_other_enzyme_family(result.get("enzyme_family", "")):
                result["is_blacklisted"] = "false"
                result["blacklist_reason"] = ""
                result["review_status"] = "needs_enzyme_family_review"
            return merge_with_fallback(result, fallback_tags, self.schema)
        except Exception as e:
            print(f"  Flavin enzyme family resolver failed: {e}")
            return fallback_tags


class ZoteroImporter:
    """从Zotero CSV导入文献"""
    
    REQUIRED_COLS = ['Title', 'Author', 'Publication Year', 'DOI', 'Abstract Note', 'Url', 'Publication Title']
    
    def __init__(self, csv_path: str, encoding: str = "auto"):
        self.csv_path = Path(csv_path)
        self.df, self.encoding = read_zotero_csv(csv_path, encoding)
        print(f"CSV编码: {self.encoding}")
        # 标准化列名（Zotero CSV列名有空格）
        self.df.columns = [c.strip() for c in self.df.columns]
        
    def parse(self) -> List[PaperNode]:
        papers = []
        seen_dois = set()
        for _, row in self.df.iterrows():
            title = str(row.get('Title', '')).strip()
            if not title or title.lower() in ('nan', 'none'):
                continue
                
            doi = str(row.get('DOI', '')).strip()
            # 过滤无效 DOI
            if doi.lower() in ('nan', 'none', ''):
                doi = ""
            # 去重
            if doi and doi in seen_dois:
                continue
            if doi:
                seen_dois.add(doi)
            
            pid = hashlib.md5(doi.encode() if doi else title.encode()).hexdigest()[:10]
            
            # 处理作者（Zotero CSV中可能是"Last, First; Last2, First2"）
            authors_raw = str(row.get('Author', '')).strip()
            authors = authors_raw.split(';')[0].strip() if ';' in authors_raw else authors_raw
            
            year = 0
            try:
                year = int(float(row.get('Publication Year', 0)))
            except:
                pass
                
            node = PaperNode(
                id=pid,
                title=title,
                authors=authors,
                year=year,
                doi=doi,
                journal=str(row.get('Publication Title', '')).strip(),
                abstract=str(row.get('Abstract Note', '')).strip(),
                url=str(row.get('Url', '')).strip(),
                file_attachments=str(row.get('File Attachments', '')).strip(),
                link_attachments=str(row.get('Link Attachments', '')).strip()
            )
            papers.append(node)
            
        print(f"成功导入 {len(papers)} 篇文献")
        return papers


class KnowledgeGraphBuilder:
    """构建文献-概念知识网络"""
    
    def __init__(self, papers: List[PaperNode]):
        self.papers = papers
        self.G = nx.Graph()
        self.concept_nodes = set()
        
    def build(self) -> nx.Graph:
        """构建二分网络：文献节点 ↔ 概念节点"""
        for p in self.papers:
            # 文献节点
            self.G.add_node(
                p.id, 
                type='paper', 
                title=p.get_display_title(40),
                year=p.year,
                journal=p.journal,
                label=p.get_display_title(30)
            )
            
            # 从tags中提取多值概念，建立连接
            concepts = self._extract_concepts(p.tags)
            for concept_type, concept_vals in concepts.items():
                for val in concept_vals:
                    if is_non_concept_tag(val):
                        continue
                    concept_id = f"{concept_type}:{val}"
                    if concept_id not in self.concept_nodes:
                        self.G.add_node(
                            concept_id, 
                            type='concept', 
                            category=concept_type,
                            label=val,
                            color=self._get_concept_color(concept_type, val)
                        )
                        self.concept_nodes.add(concept_id)
                    
                    # 文献-概念边
                    self.G.add_edge(p.id, concept_id, relation=concept_type)
                    
        return self.G
    
    def _extract_concepts(self, tags: Dict) -> Dict[str, List[str]]:
        """将LLM提取的标签拆分为可网络化的概念"""
        concepts = defaultdict(list)
        
        # 单值字段直接映射
        single_fields = [
            'topic_domain', 'search_track', 'photoenzymatic_status',
            'photoenzymatic_application_status', 'new_to_nature_reaction',
            'enzyme_category',
            'enzyme_family', 'radical_or_excited_state_mechanism',
            'radical_mechanism', 'radical_source', 'reaction_type',
            'stereochemistry', 'selectivity', 'has_mechanistic_study',
            'light_source', 'cofactor', 'catalyst_or_photosensitizer',
            'flavin_cofactor', 'characterized_enzyme_evidence',
            'candidate_for_photoenzyme_repurposing', 'structure_similarity_hint',
            'enzyme_integration_opportunity', 'lipase_or_hydrolase_relevance',
            'manual_review_priority',
        ]
        for f in single_fields:
            if f in tags and not is_non_concept_tag(tags[f]):
                concepts[f].append(tags[f])
                
        # 多值字段拆分（如bond_formation可能含"C-C | C-N"）
        multi_fields = ['bond_formation', 'substrate_class', 'mechanistic_methods', 'mechanistic_evidence']
        for f in multi_fields:
            if f in tags and not is_non_concept_tag(tags[f]):
                vals = re.split(r'[|;/,、]', tags[f])
                concepts[f].extend([v.strip() for v in vals if v.strip()])
                
        return dict(concepts)
    
    def _get_concept_color(self, cat: str, val: str) -> str:
        """根据概念类别返回颜色"""
        val_norm = val.replace(' ', '_').replace('/', '_')
        if val_norm in MORANDI_COLORS:
            return MORANDI_COLORS[val_norm]
        # 按类别返回默认色
        cat_colors = {
            'enzyme_family': '#D4C4B0', 'reaction_type': '#8B9DAF',
            'radical_mechanism': '#F0E6D3', 'bond_formation': '#A5B5C5',
            'radical_or_excited_state_mechanism': '#F0E6D3',
            'substrate_class': '#C9B8A8', 'cofactor': '#E8D5C4',
            'catalyst_or_photosensitizer': '#E8D5C4',
            'photoenzymatic_status': '#A8C4B8', 'topic_domain': '#C4B8A8',
        }
        return cat_colors.get(cat, '#CCCCCC')
    
    def export_html(self, output_path: str):
        """导出交互式HTML网络图（PyVis）"""
        if Network is None:
            print("跳过交互式知识网络: 未安装pyvis。可运行 pip install pyvis 后重新生成。")
            return
        net = Network(height="900px", width="100%", bgcolor="#fafafa", font_color="#333")
        net.barnes_hut(gravity=-3000, central_gravity=0.3, spring_length=120)
        
        # 添加节点
        for node, attr in self.G.nodes(data=True):
            if attr.get('type') == 'paper':
                size = 20
                color = '#5B7C99'  # 文献节点：深蓝灰
                title = f"{attr['title']} ({attr.get('year', '?')})"
            else:
                size = 15
                color = attr.get('color', '#999')
                title = f"{attr['category']}: {attr['label']}"
                
            net.add_node(
                node, 
                label=attr.get('label', node), 
                title=title,
                color=color,
                size=size,
                shape='dot' if attr.get('type') == 'paper' else 'box'
            )
            
        # 添加边
        for u, v, attr in self.G.edges(data=True):
            net.add_edge(u, v, title=attr.get('relation', ''))
            
        # 物理模拟参数
        net.set_options("""
        var options = {
          "physics": {
            "forceAtlas2Based": {
              "gravitationalConstant": -100,
              "centralGravity": 0.005,
              "springLength": 200,
              "springConstant": 0.18
            },
            "maxVelocity": 146,
            "solver": "forceAtlas2Based",
            "timestep": 0.35,
            "stabilization": {"iterations": 150}
          }
        }
        """)
        
        net.save_graph(output_path)
        print(f"交互式知识网络已保存: {output_path}")
        
    def analyze_clusters(self) -> Dict:
        """分析文献聚类（基于共享概念）"""
        # 构建文献-文献相似度矩阵（基于共同概念数）
        paper_ids = [p.id for p in self.papers]
        n = len(paper_ids)
        sim_matrix = np.zeros((n, n))
        
        for i, p1 in enumerate(self.papers):
            concepts1 = set()
            for neighbor in self.G.neighbors(p1.id):
                if self.G.nodes[neighbor].get('type') == 'concept':
                    concepts1.add(neighbor)
            
            for j, p2 in enumerate(self.papers):
                if i == j:
                    sim_matrix[i][j] = 1.0
                    continue
                concepts2 = set()
                for neighbor in self.G.neighbors(p2.id):
                    if self.G.nodes[neighbor].get('type') == 'concept':
                        concepts2.add(neighbor)
                        
                intersection = len(concepts1 & concepts2)
                union = len(concepts1 | concepts2)
                sim_matrix[i][j] = intersection / union if union > 0 else 0
                
        # 简单聚类：基于反应类型和酶家族
        clusters = defaultdict(list)
        for p in self.papers:
            rt = tag_get(p.tags, 'reaction_type', default='Unknown').replace(' ', '_')
            ef = tag_get(p.tags, 'enzyme_family', default='Unknown').replace(' ', '_')
            cluster_key = f"{rt}__{ef}"
            clusters[cluster_key].append(p)
            
        return {
            'similarity_matrix': sim_matrix,
            'paper_ids': paper_ids,
            'clusters': dict(clusters)
        }


class Visualizer:
    """生成静态可视化图表"""
    
    def __init__(self, papers: List[PaperNode], output_dir: Path):
        self.papers = papers
        self.output_dir = output_dir
        plt.style.use('seaborn-v0_8-whitegrid')
        configure_chinese_font()
        plt.rcParams['font.size'] = 10
        
    def plot_timeline(self):
        """时间线：年份 vs 反应类型"""
        fig, ax = plt.subplots(figsize=(12, 6))
        
        data = []
        for p in self.papers:
            if p.year > 2000:
                rt = tag_get(p.tags, 'reaction_type', default='Unknown')
                data.append({'year': p.year, 'reaction': rt, 'title': p.get_display_title(30)})
                
        df = pd.DataFrame(data)
        if df.empty:
            return
            
        # 统计每年各反应类型数量
        pivot = df.groupby(['year', 'reaction']).size().unstack(fill_value=0)
        
        # 堆叠条形图
        cmap = plt.get_cmap('viridis')
        denom = max(1, len(pivot.columns) - 1)
        colors = [cmap(i / denom) for i, _ in enumerate(pivot.columns)]
        pivot.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.7)
        
        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Number of Papers', fontsize=12)
        ax.set_title('Photoenzymatic Radical Literature Timeline', fontsize=14, pad=20)
        ax.legend(title='Reaction Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(self.output_dir / 'timeline.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("时间线图已保存")
        
    def plot_heatmap(self):
        """反应类型 × 酶家族 矩阵热图"""
        records = []
        for p in self.papers:
            rt = tag_get(p.tags, 'reaction_type', default='Unknown')
            ef = tag_get(p.tags, 'enzyme_family', default='Unknown')
            if rt != 'Unknown' and ef != 'Unknown':
                records.append((rt, ef))
                
        if not records:
            return
            
        df = pd.DataFrame(records, columns=['reaction', 'enzyme'])
        pivot = pd.crosstab(df['reaction'], df['enzyme'])
        
        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns)*1.2), max(6, len(pivot)*0.8)))
        sns.heatmap(pivot, annot=True, fmt='d', cmap='YlOrBr', 
                   linewidths=0.5, ax=ax, cbar_kws={'label': 'Paper Count'})
        ax.set_title('Reaction Type × Enzyme Family Matrix', fontsize=14, pad=20)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'reaction_enzyme_matrix.png', dpi=300)
        plt.close()
        print("反应-酶矩阵热图已保存")
        
    def plot_mechanism_sankey(self):
        """机制流向图（简化版桑基/箭头图）"""
        # 统计各机制-反应组合
        flows = Counter()
        for p in self.papers:
            mech = tag_get(p.tags, 'radical_mechanism', default='Unknown')
            rt = tag_get(p.tags, 'reaction_type', default='Unknown')
            if mech != 'Unknown' and rt != 'Unknown':
                flows[(mech, rt)] += 1
                
        if not flows:
            return
            
        # 绘制双向条形图作为Sankey替代
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        mechs = Counter([k[0] for k in flows.keys()])
        rts = Counter([k[1] for k in flows.keys()])
        
        ax1.barh(list(mechs.keys()), list(mechs.values()), color='#B8C5D0')
        ax1.set_title('Radical Mechanism Distribution', fontsize=12)
        ax1.set_xlabel('Count')
        
        ax2.barh(list(rts.keys()), list(rts.values()), color='#8B9DAF')
        ax2.set_title('Reaction Type Distribution', fontsize=12)
        ax2.set_xlabel('Count')
        
        plt.suptitle('Mechanism → Reaction Flow Overview', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'mechanism_flow.png', dpi=300)
        plt.close()
        print("机制流向图已保存")


class ReportGenerator:
    """生成Markdown综述报告"""
    
    def __init__(self, papers: List[PaperNode], output_dir: Path):
        self.papers = papers
        self.output_dir = output_dir
        
    def generate(self):
        """生成完整报告"""
        # 1. 主综述报告
        self._write_main_report()
        # 2.5 Other/no-enzyme 专项清单
        self._write_ambiguous_enzyme_report()
        # 3. 研究空白分析
        self._write_gap_analysis()
        
    def _write_main_report(self):
        """生成主Markdown报告"""
        lines = [
            "# 光酶自由基催化知识脉络图谱",
            "",
            f"> 生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
            f"> 文献总量: {len(self.papers)} 篇",
            "",
            "## 一、领域概览",
            "",
            "### 1.1 酶家族覆盖度",
            ""
        ]
        
        # 酶家族统计
        enzyme_counts = Counter()
        for p in self.papers:
            ef = tag_get(p.tags, 'enzyme_family', default='Unknown')
            if ef != 'Unknown':
                enzyme_counts[ef] += 1
                
        for enz, count in enzyme_counts.most_common():
            lines.append(f"- **{enz}**: {count} 篇")
            
        lines.extend([
            "",
            "### 1.2 反应类型分布",
            ""
        ])
        
        rt_counts = Counter()
        for p in self.papers:
            rt = tag_get(p.tags, 'reaction_type', default='Unknown')
            if rt != 'Unknown':
                rt_counts[rt] += 1
                
        for rt, count in rt_counts.most_common():
            lines.append(f"- **{rt}**: {count} 篇")
            
        lines.extend([
            "",
            "### 1.3 自由基机制分类",
            ""
        ])
        
        mech_counts = Counter()
        for p in self.papers:
            m = tag_get(p.tags, 'radical_mechanism', default='Unknown')
            if m != 'Unknown':
                mech_counts[m] += 1
                
        for m, count in mech_counts.most_common():
            lines.append(f"- **{m}**: {count} 篇")
            
        lines.extend([
            "",
            "## 二、时间演化趋势",
            "",
            "![Timeline](timeline.png)",
            "",
            "## 三、反应-酶映射矩阵",
            "",
            "![Matrix](reaction_enzyme_matrix.png)",
            "",
            "## 四、机制流向概览",
            "",
            "![Mechanism Flow](mechanism_flow.png)",
            "",
            "## 五、核心创新点提炼",
            ""
        ])
        
        # 按年份排序，列出每篇的核心创新
        sorted_papers = sorted(self.papers, key=lambda x: x.year, reverse=True)
        for p in sorted_papers:
            innovation = p.tags.get('key_innovation', 'N/A')
            if innovation != 'N/A':
                lines.append(f"### {p.get_display_title(60)} ({p.year})")
                lines.append(f"- **DOI**: {p.doi if p.doi else 'N/A'}")
                lines.append(f"- **创新**: {innovation}")
                lines.append(f"- **局限**: {p.tags.get('limitation', 'N/A')}")
                lines.append("")
                
        lines.extend([
            "---",
            "*Generated by Zotero Photoenzyme Knowledge Agent*"
        ])
        
        path = self.output_dir / "knowledge_report.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"主报告已保存: {path}")
        
    def _write_clustered_cards(self):
        """按主题聚类输出文献卡片（Obsidian-friendly）"""
        # 基于反应类型聚类
        clusters = defaultdict(list)
        for p in self.papers:
            rt = tag_get(p.tags, 'reaction_type', default='Uncategorized')
            clusters[rt].append(p)
            
        for cluster_name, papers in clusters.items():
            safe_name = cluster_name.replace(' ', '_').replace('/', '_')
            lines = [
                f"# {cluster_name} 文献集群",
                "",
                f"共 {len(papers)} 篇文献",
                ""
            ]
            
            for p in papers:
                lines.extend([
                    f"## {p.get_display_title(80)}",
                    "",
                    f"- **Authors**: {p.authors}",
                    f"- **Year**: {p.year}",
                    f"- **Journal**: {p.journal}",
                    f"- **DOI**: [{p.doi}](https://doi.org/{p.doi})" if p.doi else "",
                    f"- **URL**: {p.url}" if p.url else "",
                    "",
                    f"### 知识标签",
                    f"- 主题领域: {tag_get(p.tags, 'topic_domain')}",
                    f"- 文献域: {tag_get(p.tags, 'paper_domain')}",
                    f"- Flavin任务分轨: {tag_get(p.tags, 'search_track')}",
                    f"- Photoenzyme应用状态: {tag_get(p.tags, 'photoenzymatic_application_status')}",
                    f"- New-to-nature: {tag_get(p.tags, 'new_to_nature_reaction')}",
                    f"- Flavin辅因子: {tag_get(p.tags, 'flavin_cofactor')}",
                    f"- 酶/光酶状态: {tag_get(p.tags, 'enzyme_category')}",
                    f"- 酶家族: {tag_get(p.tags, 'enzyme_family')}",
                    f"- 酶表征证据: {tag_get(p.tags, 'characterized_enzyme_evidence')}",
                    f"- Photoenzyme候选价值: {tag_get(p.tags, 'candidate_for_photoenzyme_repurposing')}",
                    f"- 结构相似线索: {tag_get(p.tags, 'structure_similarity_hint')}",
                    f"- 自由基/激发态机制: {tag_get(p.tags, 'radical_mechanism')}",
                    f"- 自由基来源: {tag_get(p.tags, 'radical_source')}",
                    f"- 成键类型: {tag_get(p.tags, 'bond_formation')}",
                    f"- 底物: {tag_get(p.tags, 'substrate_class')}",
                    f"- 选择性/立体化学: {tag_get(p.tags, 'stereochemistry')}",
                    f"- 光源: {tag_get(p.tags, 'light_source')}",
                    f"- 催化剂/辅因子: {tag_get(p.tags, 'cofactor')}",
                    f"- 接入酶口袋机会: {tag_get(p.tags, 'enzyme_integration_opportunity')}",
                    f"- lipase/hydrolase关联: {tag_get(p.tags, 'lipase_or_hydrolase_relevance')}",
                    f"- 人工复核优先级: {tag_get(p.tags, 'manual_review_priority')}",
                    f"- API关键词/概念: {tag_get(p.tags, 'metadata_keywords')}",
                    "",
                    f"### 核心发现",
                    f"{p.tags.get('key_innovation', 'N/A')}",
                    "",
                    f"### 摘要",
                    f"> {p.abstract[:500]}..." if len(p.abstract) > 500 else f"> {p.abstract}",
                    "",
                    "---",
                    ""
                ])
                
            path = self.output_dir / f"cluster_{safe_name}.md"
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\n'.join(filter(None, lines)))
                
        print(f"已生成 {len(clusters)} 个主题聚类文件")

    def _write_ambiguous_enzyme_report(self):
        """输出Other/no-enzyme完整文献名，方便人工复核。"""
        rows = []
        for p in self.papers:
            family = tag_get(p.tags, 'enzyme_family')
            if is_other_enzyme_family(family) or family == 'no enzyme':
                rows.append(p)
        rows.sort(key=rank_paper_for_quota, reverse=True)
        lines = [
            "# Other / No-enzyme 文献复核清单",
            "",
            f"共 {len(rows)} 篇保留在最终集合中。",
            "",
        ]
        for p in rows:
            lines.extend([
                f"## {p.title}",
                "",
                f"- **Enzyme family**: {tag_get(p.tags, 'enzyme_family')}",
                f"- **Paper domain**: {tag_get(p.tags, 'paper_domain')}",
                f"- **DOI**: {p.doi if p.doi else 'missing'}",
                f"- **Year**: {p.year}",
                f"- **Journal**: {p.journal}",
                f"- **Reaction**: {tag_get(p.tags, 'reaction_type')}",
                f"- **Mechanism**: {tag_get(p.tags, 'radical_mechanism')}",
                f"- **Relevance**: {p.tags.get('relevance_score', '')}",
                f"- **Citations**: {p.tags.get('citation_count', '')}",
                "",
            ])
        path = self.output_dir / "ambiguous_enzyme_review.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"Other/no-enzyme复核清单已保存: {path}")
        
    def _write_gap_analysis(self):
        """研究空白与机会分析"""
        # 统计存在的组合
        existing = set()
        for p in self.papers:
            rt = tag_get(p.tags, 'reaction_type', default='')
            ef = tag_get(p.tags, 'enzyme_family', default='')
            bond = tag_get(p.tags, 'bond_formation', default='')
            if rt and ef:
                existing.add((rt, ef, bond))
                
        # 定义潜在的高价值组合（基于领域常识）
        high_value_gaps = [
            ("C-H_arylation", "ERED", "C-C"),
            ("decarboxylative_coupling", "P450", "C-C"),
            ("halogenation", "FAP", "C-X"),
            ("cycloaddition", "Aldolase", "C-C"),
            ("Giese_addition", "Lipase", "C-C")
        ]
        
        lines = [
            "# 研究空白与机会分析",
            "",
            "## 当前覆盖的组合",
            ""
        ]
        
        for rt, ef, bond in sorted(existing):
            lines.append(f"- ✅ {rt} + {ef} → {bond}")
            
        lines.extend([
            "",
            "## 潜在高价值空白（文献中未覆盖或极少覆盖）",
            ""
        ])
        
        for combo in high_value_gaps:
            status = "✅ 已有覆盖" if combo in existing else "🔴 研究空白"
            lines.append(f"- {status}: {combo[0]} + {combo[1]} → {combo[2]}")
            
        lines.extend([
            "",
            "## 机制研究深度分析",
            ""
        ])
        
        mechanistic = [
            p for p in self.papers
            if not is_empty_tag(tag_get(p.tags, 'mechanistic_methods', default=''))
            or not is_empty_tag(tag_get(p.tags, 'mechanistic_evidence', default=''))
            or p.tags.get('has_mechanistic_study', '') == '是'
        ]
        lines.append(f"- 含机制研究的文献: {len(mechanistic)}/{len(self.papers)} ({len(mechanistic)/len(self.papers)*100:.1f}%)")
        
        methods = Counter()
        for p in mechanistic:
            m = tag_get(p.tags, 'mechanistic_methods', default='')
            for method in re.split(r'[|;/,、]', m):
                method = method.strip()
                if not is_empty_tag(method):
                    methods[method] += 1
                
        lines.append("")
        lines.append("| 方法 | 使用次数 |")
        lines.append("|------|---------|")
        for m, c in methods.most_common():
            lines.append(f"| {m} | {c} |")
            
        lines.extend([
            "",
            "## 未来方向建议",
            ""
        ])
        
        # 从文献中提取future_direction
        directions = [
            p.tags.get('future_direction', '') for p in self.papers
            if not is_empty_tag(p.tags.get('future_direction', ''))
        ]
        for i, d in enumerate(set(directions), 1):
            lines.append(f"{i}. {d}")
            
        path = self.output_dir / "gap_analysis.md"
        with open(path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        print(f"空白分析已保存: {path}")


def main():
    parser = argparse.ArgumentParser(description='Zotero Photoenzyme Knowledge Agent')
    parser.add_argument('--input', '-i', required=True, help='Zotero导出的CSV文件路径')
    parser.add_argument('--output', '-o', default='./knowledge_base', help='输出目录')
    parser.add_argument('--csv-encoding', default='auto', help='CSV编码，默认auto；可指定utf-8-sig/utf-16/gb18030等')
    parser.add_argument('--use-llm', action='store_true', help='调用LLM API进行语义提取（需要配置API Key）')
    parser.add_argument('--api-key', help='LLM API Key（优先级高于环境变量KIMI_API_KEY/MOONSHOT_API_KEY）')
    parser.add_argument('--base-url', default=LLM_CONFIG['base_url'], help='OpenAI兼容API base URL')
    parser.add_argument('--model', default=LLM_CONFIG['model'], help='LLM模型名')
    parser.add_argument('--temperature', type=float, default=LLM_CONFIG['temperature'], help='模型温度；kimi-k2.6通常要求为1')
    parser.add_argument('--max-tokens', type=int, default=LLM_CONFIG['max_tokens'], help='每篇文献JSON输出的最大token数')
    parser.add_argument('--connect-timeout', type=float, default=LLM_CONFIG['connect_timeout'], help='API连接超时秒数')
    parser.add_argument('--read-timeout', type=float, default=LLM_CONFIG['read_timeout'], help='API等待模型输出的读取超时秒数')
    parser.add_argument('--no-json-mode', action='store_true', help='禁用response_format json_object，部分模型/端点不支持时使用')
    parser.add_argument('--schema-profile', choices=sorted(SCHEMA_PROFILES.keys()), default='photo-radical', help='内置提取schema')
    parser.add_argument('--schema-file', help='外部JSON schema文件；提供后覆盖 --schema-profile')
    parser.add_argument('--task-prompt', help='高层研究目标/长难句任务说明；会注入LLM提取提示词并写入task_spec.md')
    parser.add_argument('--task-file', help='从UTF-8文本文件读取高层研究目标；可与--task-prompt叠加')
    parser.add_argument('--criteria-file', help='JSON或key:value文本，定义field/enzyme_feature/reaction/goal四个筛选条件')
    parser.add_argument('--field', action='append', help='硬指标1：领域/field，可写长短语；支持photoenzymatic/photoenzyme等模糊同义词')
    parser.add_argument('--enzyme-feature', action='append', help='硬指标2：酶特征/辅因子，如flavin cofactor、flavin-dependent、FAD、FMN')
    parser.add_argument('--reaction', action='append', help='硬指标3：反应条件，可写长难句，如new-to-nature且生物体内罕见、化工常见')
    parser.add_argument('--goal', action='append', help='最终目标；用于任务说明、排序和人工复核，不作为必须直接命中的短关键词')
    parser.add_argument('--criteria-action', choices=['annotate', 'filter', 'llm-gate'], default='llm-gate', help='criteria处理方式：只标注、过滤掉fail、或仅用于减少LLM调用')
    parser.add_argument('--no-export-homolog-seeds', action='store_true', help='不导出homolog_candidate_seeds.csv')
    parser.add_argument('--fetch-uniprot-seeds', action='store_true', help='为enzyme_seed_candidates.csv候选酶查询UniProt并导出FASTA；需要网络')
    parser.add_argument('--max-enzyme-seeds', type=int, default=100, help='最多导出多少个enzyme seed候选')
    parser.add_argument('--uniprot-per-seed', type=int, default=2, help='每个候选酶名最多保留多少个UniProt命中')

    # Discovery 扩展参数
    parser.add_argument('--discover', action='store_true', help='启用 Semantic Scholar 文献发现与引用网络扩展')
    parser.add_argument('--ss-api-key', default=os.getenv('SEMANTIC_SCHOLAR_API_KEY'), help='Semantic Scholar API Key（可选；免费版通常无需）')
    parser.add_argument('--expand-mode', choices=['citations', 'references', 'both', 'related'], default='both', help='扩展方向')
    parser.add_argument('--max-depth', type=int, default=1, help='扩展深度（1=直接邻居，2=邻居的邻居）')
    parser.add_argument('--breadth-limit', type=int, default=20, help='每篇种子文献最多扩展多少邻居')
    parser.add_argument('--relevance-threshold', type=float, default=0.25, help='相关性阈值（0-1），低于此值的发现文献会被过滤')
    parser.add_argument('--search-query', action='append', help='额外关键词搜索查询（可多次使用）')
    parser.add_argument('--ss-rate-limit', type=float, default=3.0, help='Semantic Scholar request interval in seconds')
    parser.add_argument('--ss-max-retries', type=int, default=0, help='Semantic Scholar 429 retry count; default 0 avoids long waits')
    parser.add_argument('--ss-search-limit', type=int, default=25, help='每条Semantic Scholar关键词搜索最多取多少篇；降低可减少429')
    parser.add_argument('--ss-max-search-queries', type=int, default=4, help='最多执行多少条Semantic Scholar关键词搜索；0=不执行补充关键词搜索')
    parser.add_argument('--ss-continue-after-429', action='store_true', help='429后仍继续后续Semantic Scholar请求；默认首次429即停止SS请求')
    parser.add_argument('--export-citation-network', action='store_true', help='导出引用网络为 GEXF（Gephi）和边列表 CSV')
    parser.add_argument('--export-cluster-md', action='store_true', help='额外导出所有cluster_*.md文献卡片；默认只用HTML看板浏览')
    parser.add_argument('--no-dashboard', action='store_true', help='不生成literature_dashboard.html')
    parser.add_argument('--fetch-introduction-for-other', action='store_true', default=True, help='Other/unclear二次精查前尝试从DOI网页抓取introduction片段')
    parser.add_argument('--no-fetch-introduction-for-other', dest='fetch_introduction_for_other', action='store_false', help='禁用Other introduction网页抓取')
    parser.add_argument('--target-non-other', type=int, default=100, help='最终报告中预期的非Other/no-enzyme可注释文献下限')
    parser.add_argument('--max-other', type=int, default=20, help='最终报告中Other_enzyme最多保留多少篇')
    parser.add_argument('--disable-quality-filter', action='store_true', help='禁用Other_enzyme/no-enzyme黑名单与配额过滤')
    parser.add_argument('--resolve-other-with-llm', action='store_true', default=True, help='对Other/unclear/noisy enzyme family进行LLM二次精查')
    parser.add_argument('--no-resolve-other-with-llm', dest='resolve_other_with_llm', action='store_false', help='禁用Other enzyme LLM二次精查')
    parser.add_argument('--blacklist-memory', help='永久黑名单/注释记忆JSON路径；默认写入输出目录/literature_blacklist_memory.json')
    parser.add_argument('--metadata-cache', help='Crossref/OpenAlex元数据缓存JSON路径；默认写入输出目录/literature_metadata_cache.json')
    parser.add_argument('--no-external-metadata', action='store_true', help='禁用Crossref/OpenAlex关键词/摘要补充')
    parser.add_argument('--zotero-sqlite', help='可选：Zotero本地zotero.sqlite路径，用于标记already_in_zotero')
    parser.add_argument('--llm-scope', choices=['all', 'ambiguous', 'other-only', 'enzyme-family', 'criteria', 'criteria-ambiguous'], default='all', help='LLM处理范围；enzyme-family只澄清flavin/Other/no-enzyme家族；criteria-ambiguous先按硬指标筛，再只精查模糊项')
    parser.add_argument('--enable-kimi-thinking', action='store_true', help='启用Kimi K2 thinking；默认关闭以避免抽取任务超时')
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not args.blacklist_memory:
        args.blacklist_memory = str(output_dir / "literature_blacklist_memory.json")
    if not args.metadata_cache:
        args.metadata_cache = str(output_dir / "literature_metadata_cache.json")
    schema = load_schema(args.schema_profile, args.schema_file)
    task_prompt = read_task_prompt(args.task_prompt, args.task_file)
    criteria = build_research_criteria(args)
    criteria_prompt = criteria.to_prompt()
    if criteria_prompt:
        task_prompt = "\n\n".join(part for part in [task_prompt, criteria_prompt] if part)
        if args.llm_scope == "all" and args.criteria_action == "llm-gate":
            args.llm_scope = "criteria-ambiguous"
    auto_queries = default_search_queries_for_task(args.schema_profile, task_prompt)
    criteria_queries = default_search_queries_for_criteria(criteria)
    if criteria_queries:
        auto_queries = dedupe_preserve_order(auto_queries + criteria_queries)
    if auto_queries:
        existing_queries = set(args.search_query or [])
        args.search_query = list(args.search_query or [])
        for q in auto_queries:
            if q not in existing_queries:
                args.search_query.append(q)
        print(f"[Task] 根据研究目标补充 {len(auto_queries)} 条flavin/photoenzyme搜索式")
    write_task_spec(output_dir, args.schema_profile, task_prompt, args.search_query)
    write_criteria_spec(output_dir, criteria, args.search_query)
    memory = LiteratureMemory(Path(args.blacklist_memory))
    memory_migrated = memory.migrate_legacy_other_blacklist()
    if memory_migrated:
        memory.save()
    
    # 1. 导入
    print("=" * 60)
    print("STEP 1: 从Zotero CSV导入文献")
    print("=" * 60)
    importer = ZoteroImporter(args.input, args.csv_encoding)
    papers = importer.parse()
    zotero_seed_index = build_zotero_index_from_papers(papers)
    zotero_library_index = merge_zotero_indices(
        zotero_seed_index,
        build_zotero_index_from_sqlite(args.zotero_sqlite) if args.zotero_sqlite else {"dois": set(), "titles": set(), "journals": set()},
    )
    annotate_zotero_status(papers, zotero_seed_index, zotero_library_index)
    
    if not papers:
        print("未找到有效文献，请检查CSV文件")
        return

    initial_count = len(papers)
    papers = [p for p in papers if not memory.apply(p)]
    if len(papers) != initial_count:
        print(f"[Memory] 命中永久黑名单，跳过 {initial_count - len(papers)} 篇种子文献")

    # === Discovery 扩展（可选）===
    discovery_results = None
    citation_analyzer = None
    if args.discover:
        if not _DISCOVERY_AVAILABLE:
            print(f"\n[Discovery] 不可用: {_DISCOVERY_IMPORT_ERROR}")
            print("请确保 discovery.py 在同目录，且已安装 requests/networkx/numpy。")
            print("跳过 discovery，继续原有流程...\n")
        else:
            print("\n" + "=" * 60)
            print("STEP 1.5: 文献发现与引用网络扩展 (Semantic Scholar)")
            print("=" * 60)

            client = SemanticScholarClient(
                api_key=args.ss_api_key,
                rate_limit_delay=args.ss_rate_limit,
                max_retries=args.ss_max_retries,
                stop_on_rate_limit=not args.ss_continue_after_429,
            )
            discovery = LiteratureDiscovery(
                client,
                relevance_threshold=args.relevance_threshold,
                max_papers=max(500, args.target_non_other + args.max_other + 100),
            )

            # 扩展
            discovery.expand(
                papers,
                mode=args.expand_mode,
                max_depth=args.max_depth,
                breadth_limit=args.breadth_limit,
            )

            # 额外关键词搜索
            if args.search_query and args.ss_max_search_queries != 0:
                search_queries = args.search_query
                if args.ss_max_search_queries > 0:
                    search_queries = args.search_query[:args.ss_max_search_queries]
                    skipped_queries = len(args.search_query) - len(search_queries)
                    if skipped_queries > 0:
                        print(f"[Discovery] Skipping {skipped_queries} extra keyword searches due to --ss-max-search-queries={args.ss_max_search_queries}")
                for q in search_queries:
                    if getattr(client, "rate_limit_stopped", False):
                        print("[Discovery] Semantic Scholar rate-limit circuit is open; skipping remaining keyword searches.")
                        break
                    print(f"[Discovery] 执行补充搜索: {q}")
                    discovery.search(q, limit=args.ss_search_limit)
            elif args.search_query and args.ss_max_search_queries == 0:
                print("[Discovery] Supplemental Semantic Scholar keyword searches disabled by --ss-max-search-queries=0")

            # 构建引用网络
            citation_analyzer = CitationNetworkAnalyzer()
            citation_analyzer.build(discovery)

            # 补充种子的真实引用数（否则种子全显示 0）
            if getattr(client, "rate_limit_stopped", False):
                print("[Network] Skipping seed citation enrichment because Semantic Scholar is rate-limited.")
            else:
                citation_analyzer.enrich_seed_citations(discovery, client)

            # 打印网络统计
            stats = citation_analyzer.network_stats()
            print(f"\n[Discovery] 网络统计: {stats['total_nodes']} 节点, {stats['total_edges']} 边")
            print(f"[Discovery]   种子: {stats['seed_nodes']} | 新发现: {stats['discovered_nodes']}")

            # 将发现的文献转换为 PaperNode 并入主列表
            seed_dois = {p.doi for p in papers if p.doi}
            new_nodes_count = 0
            discovered_candidates = sorted(
                discovery.discovered.values(),
                key=lambda dp: (dp.relevance_score, dp.citation_count, dp.year),
                reverse=True,
            )
            for dp in discovered_candidates:
                if not normalize_doi_for_filter(dp.doi):
                    continue
                if dp.doi in seed_dois:
                    continue
                # 避免重复 DOI
                if any(p.doi == dp.doi for p in papers):
                    continue

                node = PaperNode(
                    id=dp.doi.replace("/", "_").replace(".", "_"),
                    title=dp.title,
                    authors=", ".join(dp.authors[:3]) + (" et al." if len(dp.authors) > 3 else ""),
                    year=dp.year,
                    doi=dp.doi,
                    journal=dp.venue or "Unknown",
                    abstract=dp.abstract or "",
                    url=f"https://doi.org/{dp.doi}" if dp.doi else "",
                    file_attachments="",
                    link_attachments="",
                )
                # 预置 discovery 标签，后续 LLM 会覆盖/补充
                node.tags = {
                    "relevance_score": f"{dp.relevance_score:.2f}",
                    "discovery_depth": str(dp.discovery_depth),
                    "citation_count": str(dp.citation_count),
                    "is_discovered": "true",
                    "metadata_keywords": " | ".join(dp.fields_of_study or []),
                }
                if memory.apply(node):
                    continue
                papers.append(node)
                new_nodes_count += 1

            print(f"[Discovery] 合并后总文献: {len(papers)} 篇（新增 {new_nodes_count} 篇）\n")
            discovery_results = discovery
            annotate_zotero_status(papers, zotero_seed_index, zotero_library_index)

    if not args.no_external_metadata:
        print("\n" + "=" * 60)
        print("STEP 1.7: Crossref/OpenAlex关键词与摘要补充")
        print("=" * 60)
        enrich_external_metadata(papers, Path(args.metadata_cache), enabled=True)

    if criteria.enabled():
        print("\n" + "=" * 60)
        print("STEP 1.8: Criteria triage (field/enzyme/reaction)")
        print("=" * 60)
        for p in papers:
            discovery_tags = dict(p.tags)
            fallback_tags = fallback_extract(p.title, build_extraction_context(p), schema)
            p.tags = {**discovery_tags, **fallback_tags}
        apply_research_criteria(papers, criteria)
        if args.criteria_action == "filter":
            papers = filter_papers_by_criteria(papers, output_dir)
            annotate_zotero_status(papers, zotero_seed_index, zotero_library_index)

    # 2. LLM语义提取（可选但强烈推荐）
    if args.use_llm:
        print("\n" + "=" * 60)
        print("STEP 2: LLM语义提取（每篇约3-5秒）")
        print("=" * 60)
        
        if args.api_key:
            LLM_CONFIG['api_key'] = args.api_key
        print_auth_diagnostic(LLM_CONFIG.get('api_key', ''), get_api_key_source(args.api_key))
        try:
            LLM_CONFIG['api_key'] = normalize_api_key(LLM_CONFIG.get('api_key', ''))
        except ValueError as e:
            print(f"API Key配置错误: {e}")
            return
        LLM_CONFIG['base_url'] = args.base_url
        LLM_CONFIG['model'] = args.model
        LLM_CONFIG['temperature'] = args.temperature
        LLM_CONFIG['max_tokens'] = args.max_tokens
        LLM_CONFIG['connect_timeout'] = args.connect_timeout
        LLM_CONFIG['read_timeout'] = args.read_timeout
        LLM_CONFIG['json_mode'] = not args.no_json_mode
        LLM_CONFIG['enable_thinking'] = args.enable_kimi_thinking
        if str(LLM_CONFIG['model']).lower().startswith("kimi-k2") and LLM_CONFIG['temperature'] != 0.6:
            print(f"提示: {LLM_CONFIG['model']} 仅允许 temperature=0.6，已自动覆盖 {LLM_CONFIG['temperature']} -> 0.6")
            LLM_CONFIG['temperature'] = 0.6
            
        if not LLM_CONFIG['api_key']:
            print("警告: 未配置API Key，请设置环境变量 KIMI_API_KEY/MOONSHOT_API_KEY 或使用 --api-key 参数")
            print("跳过LLM提取，使用基于关键词的fallback...")
            # Fallback: 基于关键词的正则匹配
            for p in papers:
                discovery_tags = dict(p.tags)
                p.tags = {**discovery_tags, **fallback_extract(p.title, build_extraction_context(p), schema)}
                if criteria.enabled():
                    apply_research_criteria_to_paper(p, criteria)
        else:
            extractor = LLMExtractor(LLM_CONFIG, schema, task_prompt=task_prompt)
            llm_candidates = []
            for p in papers:
                discovery_tags = dict(p.tags)
                fallback_tags = fallback_extract(p.title, build_extraction_context(p), schema)
                p.tags = {**discovery_tags, **fallback_tags}
                if criteria.enabled():
                    apply_research_criteria_to_paper(p, criteria)
                if should_run_llm_for_paper(p, fallback_tags, args.llm_scope):
                    llm_candidates.append(p)
            print(f"\n  LLM候选文献: {len(llm_candidates)}/{len(papers)} (scope={args.llm_scope})")
            for i, p in enumerate(llm_candidates, 1):
                print(f"  [{i}/{len(llm_candidates)}] 提取: {p.get_display_title(50)}")
                discovery_tags = dict(p.tags)
                if args.llm_scope == "enzyme-family":
                    extracted = extractor.resolve_flavin_enzyme_family(p)
                else:
                    extracted = extractor.extract(p.title, build_extraction_context(p), p.journal)
                p.tags = {**discovery_tags, **extracted}
                if criteria.enabled():
                    apply_research_criteria_to_paper(p, criteria)
                time.sleep(0.5)  # 避免rate limit

            if args.resolve_other_with_llm:
                candidates = [
                    p for p in papers
                    if p.tags.get("criteria_status", "") != "criteria_fail" and (
                    is_other_enzyme_family(tag_get(p.tags, 'enzyme_family'))
                    or tag_get(p.tags, 'paper_domain', default='') in {'', 'unclear', 'unrelated/blacklist'}
                    )
                ]
                print(f"\n  Other/no-enzyme候选二次精查: {len(candidates)} 篇")
                for j, p in enumerate(candidates, 1):
                    print(f"    [{j}/{len(candidates)}] 精查酶家族: {p.get_display_title(60)}")
                    if args.fetch_introduction_for_other and not p.tags.get("introduction_snippet"):
                        intro = fetch_introduction_snippet(p)
                        if intro:
                            p.tags["introduction_snippet"] = intro
                    discovery_tags = dict(p.tags)
                    if args.schema_profile == "flavin-photoenzyme":
                        refined = extractor.resolve_flavin_enzyme_family(p)
                    else:
                        refined = extractor.resolve_other_enzyme(p)
                    p.tags = {**discovery_tags, **refined}
                    if criteria.enabled():
                        apply_research_criteria_to_paper(p, criteria)
                    time.sleep(0.5)
    else:
        print("\n" + "=" * 60)
        print("STEP 2: 使用本地正则提取（无需API）")
        print("=" * 60)
        for p in papers:
            discovery_tags = dict(p.tags)
            p.tags = {**discovery_tags, **fallback_extract(p.title, build_extraction_context(p), schema)}
            if criteria.enabled():
                apply_research_criteria_to_paper(p, criteria)

    blacklisted_papers: List[PaperNode] = []
    if not args.disable_quality_filter:
        print("\n" + "=" * 60)
        print("STEP 2.5: 注释质量过滤（Other/no enzyme/黑名单）")
        print("=" * 60)
        papers, blacklisted_papers = apply_annotation_quality_filter(
            papers,
            target_non_other=args.target_non_other,
            max_other=args.max_other,
        )
        for paper in papers:
            memory.remember(paper, "kept", "")
        for paper in blacklisted_papers:
            memory.remember(paper, "blacklisted", paper.tags.get("blacklist_reason", ""))
        memory.save()
        write_blacklist_report(blacklisted_papers, output_dir)
        write_deferred_review_report(papers, output_dir)
    annotate_zotero_status(papers, zotero_seed_index, zotero_library_index)
    write_zotero_overlap_report(papers, zotero_seed_index, zotero_library_index, output_dir)
    if not args.no_export_homolog_seeds:
        export_homolog_candidate_seeds(papers, output_dir)
        export_enzyme_seed_candidates(
            papers,
            output_dir,
            fetch_uniprot=args.fetch_uniprot_seeds,
            max_seeds=args.max_enzyme_seeds,
            uniprot_per_seed=args.uniprot_per_seed,
        )
        export_sequence_seed_review(output_dir)

    # 3. 构建知识网络
    print("\n" + "=" * 60)
    print("STEP 3: 构建知识图谱")
    print("=" * 60)
    kg_builder = KnowledgeGraphBuilder(papers)
    G = kg_builder.build()
    kg_builder.export_html(str(output_dir / "knowledge_graph.html"))
    
    # 4. 可视化
    print("\n" + "=" * 60)
    print("STEP 4: 生成可视化")
    print("=" * 60)
    viz = Visualizer(papers, output_dir)
    viz.plot_timeline()
    viz.plot_heatmap()
    viz.plot_mechanism_sankey()
    
    # 5. 生成报告
    print("\n" + "=" * 60)
    print("STEP 5: 生成Markdown报告")
    print("=" * 60)
    reporter = ReportGenerator(papers, output_dir)
    reporter.generate()
    if args.export_cluster_md:
        reporter._write_clustered_cards()
    if not args.no_dashboard:
        write_literature_dashboard(papers, output_dir, citation_analyzer)

    # === Discovery 扩展报告（可选）===
    if args.discover and _DISCOVERY_AVAILABLE and discovery_results and citation_analyzer:
        print("\n" + "=" * 60)
        print("STEP 6: 生成 Discovery 引用网络分析报告")
        print("=" * 60)
        disc_reporter = DiscoveryReportGenerator(discovery_results, citation_analyzer, output_dir)
        disc_reporter.generate()

        if args.export_citation_network:
            citation_analyzer.export_gexf(str(output_dir / "citation_network.gexf"))
            citation_analyzer.export_edgelist(str(output_dir / "citation_edges.csv"))
        citation_analyzer.export_html(str(output_dir / "citation_network.html"))
    
    print("\n" + "=" * 60)
    print("完成！输出文件:")
    print("=" * 60)
    for f in sorted(output_dir.iterdir()):
        print(f"  - {f.name}")


def _keyword_in_text(text: str, keyword: str) -> bool:
    key = str(keyword or "").lower()
    if not key:
        return False
    if len(key) <= 4 and re.fullmatch(r"[a-z0-9]+", key):
        return re.search(rf"(?<![a-z0-9]){re.escape(key)}(?![a-z0-9])", text) is not None
    return key in text


def _detect_enzyme_family(text: str) -> str:
    """智能酶家族检测：优先精确匹配，其次用 *ase 后缀推断"""
    enzyme_context = any(k in text for k in [
        'enzyme', 'enzymatic', 'biocatal', 'biocatalysis', 'protein',
        'active site', 'mutant', 'variant', 'catalyzes', 'catalyses',
    ])
    # 精确匹配（高置信度）
    precise = [
        (['p450', 'cytochrome p450'], 'P450'),
        (['ered', 'ene-reductase', 'old yellow enzyme', 'oye', 'ene reductase'], 'ERED/OYE'),
        (['fap', 'fatty acid photodecarboxylase', 'photodecarboxylase', 'cvfap'], 'FAP'),
        (['flavin-dependent halogenase', 'flavin dependent halogenase', 'fdh ', 'fdhs'], 'flavin-dependent halogenase'),
        (['photolyase', 'cryptochrome'], 'photolyase/cryptochrome'),
        (['berberine bridge enzyme', 'bbe'], 'berberine bridge enzyme'),
        (['flavin monooxygenase', 'flavin-containing monooxygenase', 'fmo'], 'monooxygenase'),
        (['flavin reductase'], 'flavin reductase'),
        (['acyl-coa dehydrogenase', 'acyl coa dehydrogenase'], 'acyl-CoA dehydrogenase'),
        (['d-amino acid oxidase', 'd amino acid oxidase', 'daao'], 'D-amino acid oxidase'),
        (['glucose oxidase', 'gox'], 'glucose oxidase'),
        (['amine oxidase'], 'amine oxidase'),
        (['sulfide:quinone oxidoreductase', 'sulfide quinone oxidoreductase', 'sqr'], 'sulfide:quinone oxidoreductase'),
        (['luciferase'], 'luciferase'),
        (['aldolase', 'class i aldolase', 'class ii aldolase'], 'Aldolase'),
        (['lipase', 'esterase', 'triacylglycerol lipase'], 'lipase/esterase'),
        (['transaminase', 'aminotransferase'], 'transaminase'),
        (['dehydrogenase'], 'dehydrogenase'),
        (['monooxygenase'], 'monooxygenase'),
        (['peroxidase', 'catalase'], 'peroxidase'),
        (['nitroreductase'], 'nitroreductase'),
        (['halogenase'], 'halogenase'),
        (['synthase', 'synthetase'], 'synthase'),
        (['kinase', 'phosphatase'], 'kinase/phosphatase'),
        (['protease', 'peptidase'], 'protease'),
        (['lyase'], 'lyase'),
        (['hydrolase'], 'hydrolase'),
        (['oxidase'], 'oxidase'),
        (['reductase'], 'reductase'),
    ]
    for keywords, family in precise:
        if any(_keyword_in_text(text, k) for k in keywords):
            return family
    # 兜底：用正则抓 *ase 词
    import re
    matches = re.findall(r'\b(\w+ase)\b', text)
    if matches:
        cleaned = [
            m for m in matches
            if m not in FALSE_ASE_WORDS and len(m) >= 7
        ]
        if not cleaned or not enzyme_context:
            return 'no enzyme'
        counter = Counter(cleaned)
        most_common = counter.most_common(1)[0][0]
        return most_common
    return 'no enzyme' if not enzyme_context else 'Other_enzyme'


def fallback_extract(title: str, abstract: str, schema: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """无需LLM时的增强正则fallback提取。大幅扩展关键词覆盖，减少N/A。"""
    schema = schema or EXTRACTION_SCHEMA
    text = (title + " " + abstract).lower()
    tags = {}
    
    # ========== 酶家族 ==========
    enzyme_family = _detect_enzyme_family(text)
    if 'enzyme_family' in schema:
        tags['enzyme_family'] = enzyme_family
    
    # ========== 酶类别 (schema v2) ==========
    if 'enzyme_category' in schema:
        if any(k in text for k in ['photoenzym', 'photobiocatal', 'light-driven enzym']):
            tags['enzyme_category'] = '天然辅因子依赖型'
        elif any(k in text for k in ['artificial photoenzyme', 'synthetic cofactor', 'engineered photoenzyme']):
            tags['enzyme_category'] = '人工光酶(合成辅因子)'
        elif any(k in text for k in ['dual catalysis', 'photocatalysis and enzym', 'chemoenzymatic']):
            tags['enzyme_category'] = '协同光催化/酶催化'
        elif 'enzyme' in text or 'enzymatic' in text or 'biocatal' in text or enzyme_family not in {'no enzyme', 'Other_enzyme'}:
            tags['enzyme_category'] = '酶催化'
        else:
            tags['enzyme_category'] = '无酶小分子光催化' if any(k in text for k in ['photocatalysis', 'photoredox']) else 'unclear'
    
    # ========== 反应类型 ==========
    if any(k in text for k in ['semipinacol', 'c-c bond formation', 'c-c bond construction', 'carbon-carbon bond']):
        rt = 'C-C bond formation'
    elif 'azidooxygenation' in text:
        rt = 'C-N bond formation'
    elif 'decarboxyl' in text or 'decarboxylation' in text:
        rt = 'decarboxylative coupling' if 'photoenzymatic_status' in schema else 'decarboxylative_coupling'
    elif 'arylat' in text or 'c-h arylation' in text or 'c(sp3)-h arylation' in text:
        rt = 'arylation' if 'photoenzymatic_status' in schema else 'C-H_arylation'
    elif 'alkylat' in text or 'c-h alkylation' in text or 'hydroalkylation' in text:
        rt = 'alkylation' if 'photoenzymatic_status' in schema else 'C-H_alkylation'
    elif any(k in text for k in ['halogenat', 'chlorinat', 'brominat', 'fluorinat', 'iodinat']):
        rt = 'halogenation'
    elif 'cycloaddition' in text or '[2+2]' in text or '[4+2]' in text:
        rt = 'cycloaddition'
    elif 'giese' in text or 'giese addition' in text:
        rt = 'Giese addition' if 'photoenzymatic_status' in schema else 'Giese_addition'
    elif 'minisci' in text:
        rt = 'Minisci reaction'
    elif 'cross-coupling' in text or 'cross coupling' in text:
        rt = 'cross-coupling'
    elif 'hydroamination' in text or 'carboamination' in text:
        rt = 'hydroamination/carboamination'
    elif 'cyclization' in text or 'carbocyclization' in text:
        rt = 'cyclization'
    elif 'rearrangement' in text or 'migration' in text or 'rearrang' in text:
        rt = 'rearrangement/migration'
    elif 'functionalization' in text or 'c-h functionalization' in text:
        rt = 'C-H functionalization'
    else:
        rt = 'other'
    tags['reaction_type'] = rt
    
    # ========== 自由基/激发态机制 ==========
    if 'photoredox' in text or 'photo-redox' in text or 'photoinduced electron transfer' in text:
        mech = '光氧化还原(PET)'
        mech_v2 = 'PET/photoredox'
    elif any(k in text for k in ['energy transfer', 'triplet sensitizer', 'triplet-triplet']):
        mech = '能量转移(EnT)'
        mech_v2 = 'EnT'
    elif 'hat' in text or 'hydrogen atom transfer' in text:
        mech = '氢原子转移(HAT)'
        mech_v2 = 'HAT'
    elif 'eda complex' in text or 'electron donor acceptor' in text:
        mech = 'EDA复合物'
        mech_v2 = 'EDA complex'
    elif 'single electron transfer' in text or ' set ' in text:
        mech = '单电子转移(SET)'
        mech_v2 = 'SET'
    elif 'lmct' in text:
        mech = '配体到金属电荷转移(LMCT)'
        mech_v2 = 'LMCT'
    elif 'direct excitation' in text or 'directly excited' in text:
        mech = '直接激发'
        mech_v2 = 'direct substrate excitation'
    elif 'cofactor excitation' in text or 'flavin excitation' in text or 'fad excitation' in text:
        mech = '辅因子激发'
        mech_v2 = 'cofactor excitation'
    else:
        mech = 'unclear'
        mech_v2 = 'unclear'
    tags['radical_mechanism'] = mech
    if 'radical_or_excited_state_mechanism' in schema:
        tags['radical_or_excited_state_mechanism'] = mech_v2
    
    # ========== 自由基来源 ==========
    if 'radical_source' in schema:
        if any(k in text for k in ['carboxylic acid', 'decarboxylat', 'carboxylate']):
            tags['radical_source'] = 'carboxylic acid/decarboxylation'
        elif any(k in text for k in ['alkyl halide', 'bromide', 'chloride', 'iodide']):
            tags['radical_source'] = 'alkyl halide'
        elif 'sulfonamide' in text:
            tags['radical_source'] = 'sulfonamide'
        elif 'pyridinium' in text or 'pyridin' in text:
            tags['radical_source'] = 'pyridinium salt'
        elif 'diaz' in text or 'azo' in text:
            tags['radical_source'] = 'diazo/azo precursor'
        elif 'peroxide' in text or 'perox' in text:
            tags['radical_source'] = 'peroxide'
        elif 'amine' in text or 'amino' in text:
            tags['radical_source'] = 'amine'
        elif 'alcohol' in text or 'aldehyde' in text:
            tags['radical_source'] = 'alcohol/aldehyde'
        elif 'alkene' in text or 'olefin' in text:
            tags['radical_source'] = 'alkene'
        elif 'enzyme cofactor' in text or 'flavin' in text or 'fad' in text:
            tags['radical_source'] = 'enzyme cofactor'
        else:
            tags['radical_source'] = 'unclear'
    
    # ========== 成键类型 ==========
    if 'bond_formation' in schema:
        bonds = []
        if 'c-c' in text or 'c(sp3)-c' in text or 'c(sp2)-c' in text or 'carbon-carbon' in text:
            bonds.append('C-C')
        if 'c-n' in text or 'c(sp3)-n' in text or 'amine' in text or 'amid' in text or 'carbon-nitrogen' in text:
            bonds.append('C-N')
        if 'c-o' in text or 'hydroxy' in text or 'alcohol' in text or 'ether' in text or 'carbon-oxygen' in text:
            bonds.append('C-O')
        if 'c-s' in text or 'thiol' in text or 'thio' in text or 'carbon-sulfur' in text:
            bonds.append('C-S')
        if any(k in text for k in ['c-f', 'c-cl', 'c-br', 'c-i', 'halogenat', 'fluorinat']):
            bonds.append('C-X')
        if 'c-h' in text or 'c(sp3)-h' in text or 'carbon-hydrogen' in text:
            bonds.append('C-H')
        if 'c-b' in text or 'boron' in text or 'boryl' in text:
            bonds.append('C-B')
        if 'c-si' in text or 'silicon' in text or 'silyl' in text:
            bonds.append('C-Si')
        if 'skeletal' in text or 'ring expansion' in text or 'ring contraction' in text:
            bonds.append('skeletal rearrangement')
        tags['bond_formation'] = ' | '.join(bonds) if bonds else 'unclear'
    
    # ========== 底物类别 ==========
    if 'substrate_class' in schema:
        subs = []
        if 'alkene' in text or 'olefin' in text or 'alkenyl' in text:
            subs.append('alkene')
        if 'pyridine' in text or 'pyridinium' in text or 'heteroarene' in text:
            subs.append('pyridine/heteroarene')
        if 'arene' in text or 'aromatic' in text or 'aryl' in text:
            subs.append('arene')
        if 'carboxylic acid' in text or 'carboxylate' in text:
            subs.append('carboxylic acid')
        if 'sulfonamide' in text:
            subs.append('sulfonamide')
        if 'amine' in text or 'amino' in text:
            subs.append('amine')
        if 'carbonyl' in text or 'aldehyde' in text or 'ketone' in text:
            subs.append('carbonyl compound')
        if 'imine' in text:
            subs.append('imine')
        if 'alkyl halide' in text or 'bromide' in text or 'chloride' in text:
            subs.append('alkyl halide')
        if 'drug' in text or 'pharmaceutical' in text:
            subs.append('drug-like molecule')
        if 'amino acid' in text or 'peptide' in text:
            subs.append('amino acid derivative')
        if 'lipid' in text or 'ester' in text or 'fatty acid' in text:
            subs.append('lipid/ester')
        tags['substrate_class'] = ' | '.join(subs) if subs else 'other'
    
    # ========== 选择性 ==========
    if any(k in schema for k in ['selectivity', 'stereochemistry']):
        sel = []
        if 'enantioselective' in text or 'enantiomeric' in text or 'ee ' in text:
            sel.append('enantioselective')
        if 'diastereoselective' in text or 'diastereomer' in text or 'dr ' in text:
            sel.append('diastereoselective')
        if 'regioselective' in text or 'regioisomer' in text:
            sel.append('regioselective')
        if 'chemoselective' in text:
            sel.append('chemoselective')
        if 'site-selective' in text or 'site selective' in text:
            sel.append('site-selective')
        if 'kinetic resolution' in text:
            sel.append('kinetic resolution')
        if 'stereocontrol' in text or 'stereoselective' in text:
            sel.append('stereoselective')
        val = ' | '.join(sel) if sel else ('stereocontrol unclear' if 'enzyme' in text or 'selectiv' in text else 'no selectivity focus')
        tags['selectivity'] = val
        tags['stereochemistry'] = val
    
    # ========== 机制研究证据 ==========
    if 'mechanistic_evidence' in schema or 'mechanistic_methods' in schema:
        mechs = []
        if any(k in text for k in ['epr', 'esr', 'electron paramagnetic', 'electron spin resonance']):
            mechs.append('EPR/ESR')
        if 'transient absorption' in text or 'ta spectroscopy' in text:
            mechs.append('transient absorption')
        if 'stern-volmer' in text or 'quenching' in text:
            mechs.append('Stern-Volmer')
        if 'quantum yield' in text or 'quantum efficiency' in text:
            mechs.append('quantum yield')
        if 'radical clock' in text:
            mechs.append('radical clock')
        if 'isotope label' in text or 'deuterium' in text or 'kinetic isotope' in text:
            mechs.append('isotope labeling')
        if 'dft' in text or 'density functional' in text:
            mechs.append('DFT')
        if 'crystal structure' in text or 'x-ray' in text or 'crystallography' in text:
            mechs.append('crystal structure')
        if 'md simulation' in text or 'molecular dynamics' in text:
            mechs.append('MD')
        if 'control experiment' in text:
            mechs.append('control experiments')
        if 'computational' in text or 'theoretical' in text:
            mechs.append('theoretical/computational')
        tags['mechanistic_evidence'] = ' | '.join(mechs) if mechs else 'none mentioned'
        tags['mechanistic_methods'] = tags['mechanistic_evidence']
    
    # ========== 光源 ==========
    if 'light_source' in schema:
        if 'blue led' in text or 'blue light' in text or '450 nm' in text:
            tags['light_source'] = 'blue LED'
        elif 'green led' in text or 'green light' in text or '520 nm' in text:
            tags['light_source'] = 'green LED'
        elif 'uv' in text or 'ultraviolet' in text:
            tags['light_source'] = 'UV'
        elif 'solar' in text or 'sunlight' in text:
            tags['light_source'] = 'solar'
        elif 'visible light' in text or 'white led' in text:
            tags['light_source'] = 'visible light'
        elif 'led' in text or 'light emitting diode' in text:
            tags['light_source'] = 'LED'
        else:
            tags['light_source'] = 'unclear'
    
    # ========== 催化剂/光敏剂/辅因子 ==========
    if 'catalyst_or_photosensitizer' in schema or 'cofactor' in schema or 'cofactor_or_photosensitizer' in schema:
        cats = []
        if any(k in text for k in ['flavin', 'fad', 'fmn']):
            cats.append('flavin/FAD/FMN')
        if 'nad' in text or 'nadh' in text or 'nadp' in text:
            cats.append('NAD(P)H')
        if 'ru(' in text or 'ruthenium' in text or 'bpy' in text:
            cats.append('Ru')
        if 'ir(' in text or 'iridium' in text or 'ppy' in text:
            cats.append('Ir')
        if 'organic photocatalyst' in text or 'eosin' in text or 'rose bengal' in text or 'acridinium' in text:
            cats.append('organic photocatalyst')
        if 'quinolinone' in text or 'quinolone' in text:
            cats.append('quinolinone')
        if 'palladium' in text or 'pd(' in text:
            cats.append('palladium')
        if 'enzyme active site' in text or 'active site' in text:
            cats.append('enzyme active site')
        if 'no catalyst' in text or 'catalyst-free' in text:
            cats.append('no catalyst')
        cat_str = ' | '.join(cats) if cats else 'unclear'
        tags['catalyst_or_photosensitizer'] = cat_str
        tags['cofactor'] = cat_str
        tags['cofactor_or_photosensitizer'] = cat_str
        if 'flavin_cofactor' in schema:
            if 'fad' in text:
                tags['flavin_cofactor'] = 'FAD'
            elif 'fmn' in text:
                tags['flavin_cofactor'] = 'FMN'
            elif 'flavin' in text or 'flavoenzyme' in text:
                tags['flavin_cofactor'] = 'flavin unspecified'
            else:
                tags['flavin_cofactor'] = 'no flavin'
    
    # ========== 核心创新点 / 局限 / 未来方向 ==========
    if 'key_innovation' in schema:
        # 尝试从标题中提取核心创新
        tags['key_innovation'] = f"基于标题推断: {title[:60]}..."
    if 'limitation' in schema:
        tags['limitation'] = '需LLM进一步提取'
    if 'future_direction' in schema:
        tags['future_direction'] = '需LLM进一步提取'
    
    # ========== topic_domain & photoenzymatic_status (schema v2) ==========
    if 'topic_domain' in schema:
        if any(k in text for k in ['photoenzym', 'photobiocatal', 'enzyme photocatalysis']):
            tags['topic_domain'] = 'photoenzymatic radical'
        elif any(k in text for k in ['photoredox', 'visible light', 'photocatalysis', 'energy transfer', 'photoinduced']):
            tags['topic_domain'] = 'photoredox radical reaction'
        elif any(k in text for k in ['radical', 'radicals', 'minisci', 'giese']):
            tags['topic_domain'] = 'radical chemical reaction'
        elif any(k in text for k in ['lipase', 'esterase', 'hydrolase', 'calb', 'candida antarctica']):
            tags['topic_domain'] = 'enzymatic catalysis'
        else:
            tags['topic_domain'] = 'other'

    flavin_enzyme_family_signal = enzyme_family not in {'no enzyme', 'Other_enzyme'} or any(k in text for k in [
        'ene-reductase', 'ene reductase', 'ered', 'old yellow enzyme', 'oye',
        'fatty acid photodecarboxylase', 'photodecarboxylase', 'fap',
        'flavin-dependent halogenase', 'flavin dependent halogenase', 'fdh',
        'flavin monooxygenase', 'nitroreductase', 'flavin reductase',
        'flavoenzyme', 'oxidase', 'reductase', 'dehydrogenase', 'halogenase',
    ])

    if 'paper_domain' in schema:
        has_enzyme_signal = any(k in text for k in [
            'photoenzym', 'photobiocatal', 'enzyme', 'enzymatic', 'biocatal',
            'active site', 'protein', 'flavoenzyme',
        ]) or flavin_enzyme_family_signal
        has_photo_radical_signal = any(k in text for k in [
            'photoredox', 'visible light', 'photoinduced', 'photocatalysis',
            'radical', 'minisci', 'giese', 'energy transfer', 'hat',
        ])
        has_organic_synthesis_signal = any(k in text for k in [
            'organic synthesis', 'synthetic method', 'total synthesis',
            'cross-coupling', 'functionalization', 'alkylation', 'arylation',
        ])
        if any(k in text for k in ['photoenzym', 'photobiocatal']):
            tags['paper_domain'] = 'photoenzymatic/enzymatic'
        elif has_enzyme_signal:
            tags['paper_domain'] = 'enzymatic non-photo' if not has_photo_radical_signal else 'photoenzymatic/enzymatic'
        elif has_photo_radical_signal and not flavin_enzyme_family_signal:
            tags['paper_domain'] = 'organic photoredox/radical no-enzyme'
        elif has_organic_synthesis_signal:
            tags['paper_domain'] = 'organic synthesis no-enzyme'
        else:
            tags['paper_domain'] = 'unrelated/blacklist'

    if 'photoenzymatic_application_status' in schema:
        has_flavin = tags.get('flavin_cofactor') in {'FAD', 'FMN', 'flavin unspecified'} or any(k in text for k in ['flavin', 'fad', 'fmn', 'flavoenzyme'])
        has_photoenzyme = any(k in text for k in ['photoenzym', 'photobiocatal', 'light-driven enzym', 'photoinduced enzym'])
        has_new_to_nature = any(k in text for k in [
            'new-to-nature', 'new to nature', 'new-to-nature biocatalysis',
            'non-natural', 'non natural', 'nonnatural', 'unnatural reaction',
            'abiological', 'non-native', 'non native', 'non-native transformation',
            'non-native reaction', 'noncanonical', 'non-canonical',
            'currently unknown in nature', 'unknown in nature', 'not known in nature',
            'not found in nature', 'promiscuous radical',
        ])
        has_characterized = any(k in text for k in [
            'purified', 'kinetic', 'kcat', 'km ', 'substrate scope', 'crystal structure',
            'pdb', 'activity assay', 'enzyme assay', 'mutagenesis', 'variant', 'active site',
            'mechanistic studies', 'catalyse', 'catalyze', 'catalysed', 'catalyzed',
        ])
        if has_photoenzyme and has_flavin:
            tags['photoenzymatic_application_status'] = 'applied to new-to-nature reaction' if has_new_to_nature else 'candidate not yet photoenzymatic'
            tags['new_to_nature_reaction'] = 'true' if has_new_to_nature else 'unclear'
            tags['search_track'] = 'A_used_flavin_photoenzyme_new_to_nature'
            tags['paper_domain'] = 'flavin photoenzymatic application'
        elif has_flavin and (has_characterized or flavin_enzyme_family_signal):
            tags['photoenzymatic_application_status'] = 'candidate not yet photoenzymatic'
            tags['new_to_nature_reaction'] = 'true' if has_new_to_nature else 'false'
            tags['search_track'] = 'B_characterized_flavin_enzyme_not_photoenzymatic'
            tags['paper_domain'] = 'characterized flavin enzyme non-photo'
        elif has_photoenzyme:
            tags['photoenzymatic_application_status'] = 'applied to new-to-nature reaction' if has_new_to_nature else 'unclear'
            tags['new_to_nature_reaction'] = 'true' if has_new_to_nature else 'unclear'
            tags['search_track'] = 'C_photoenzyme_non_flavin'
            tags['paper_domain'] = 'photoenzymatic non-flavin'
        elif has_photo_radical_signal and not (has_flavin and flavin_enzyme_family_signal):
            tags['photoenzymatic_application_status'] = 'no enzyme'
            tags['new_to_nature_reaction'] = 'false'
            tags['search_track'] = 'D_organic_photochemistry_no_enzyme'
            tags['paper_domain'] = 'organic photoredox/radical no-enzyme'
        else:
            tags['photoenzymatic_application_status'] = 'unclear'
            tags['new_to_nature_reaction'] = 'unclear'
            tags['search_track'] = 'E_background_or_unrelated'

        if 'characterized_enzyme_evidence' in schema:
            if any(k in text for k in ['kcat', 'km ', 'steady-state kinetic', 'kinetic parameter']):
                tags['characterized_enzyme_evidence'] = 'purified enzyme kinetics'
            elif any(k in text for k in ['activity assay', 'enzyme assay', 'substrate scope']):
                tags['characterized_enzyme_evidence'] = 'activity assay only'
            elif any(k in text for k in ['crystal structure', 'pdb', 'x-ray', 'crystallography']):
                tags['characterized_enzyme_evidence'] = 'structure only'
            elif any(k in text for k in ['genome', 'genomic', 'annotation']):
                tags['characterized_enzyme_evidence'] = 'genome annotation only'
            else:
                tags['characterized_enzyme_evidence'] = 'unclear'
        if 'candidate_for_photoenzyme_repurposing' in schema:
            if tags.get('search_track') == 'B_characterized_flavin_enzyme_not_photoenzymatic':
                tags['candidate_for_photoenzyme_repurposing'] = 'high' if any(k in text for k in ['crystal structure', 'pdb', 'substrate scope', 'kinetic']) else 'medium'
            elif tags.get('search_track') == 'A_used_flavin_photoenzyme_new_to_nature':
                tags['candidate_for_photoenzyme_repurposing'] = 'no'
            else:
                tags['candidate_for_photoenzyme_repurposing'] = 'low'
        if 'structure_similarity_hint' in schema:
            if any(k in text for k in ['crystal structure', 'pdb', 'x-ray', 'crystallography']):
                tags['structure_similarity_hint'] = 'PDB available'
            elif any(k in text for k in ['homolog', 'homologue', 'same fold', 'structural similarity']):
                tags['structure_similarity_hint'] = 'same fold/cofactor pocket'
            elif 'active site' in text or 'substrate pocket' in text:
                tags['structure_similarity_hint'] = 'similar substrate pocket'
            else:
                tags['structure_similarity_hint'] = 'unclear'

    if 'photoenzymatic_status' in schema:
        if any(k in text for k in ['photoenzym', 'photobiocatal']):
            tags['photoenzymatic_status'] = '已接入photoenzymatic'
        elif any(k in text for k in ['enzyme', 'enzymatic', 'biocatal', 'ase ']):
            tags['photoenzymatic_status'] = '仅酶催化背景'
        elif any(k in text for k in ['photoredox', 'visible light', 'photoinduced', 'radical']):
            tags['photoenzymatic_status'] = '仅光化学/自由基反应'
        else:
            tags['photoenzymatic_status'] = 'unclear'

    if 'enzyme_integration_opportunity' in schema:
        if any(k in text for k in ['enantioselective', 'stereoselective', 'site-selective', 'regioselective', 'drug', 'natural product']):
            tags['enzyme_integration_opportunity'] = 'medium'
        elif any(k in text for k in ['enzyme', 'active site', 'cofactor', 'flavin', 'fad', 'fmn', 'ase ']):
            tags['enzyme_integration_opportunity'] = 'high'
        else:
            tags['enzyme_integration_opportunity'] = 'unclear'

    if 'lipase_or_hydrolase_relevance' in schema:
        if any(k in text for k in ['lipase', 'esterase', 'hydrolase', 'calb', 'candida antarctica']):
            tags['lipase_or_hydrolase_relevance'] = 'direct'
        elif any(k in text for k in ['enzyme', 'active site', 'substrate recognition', 'stereoselective']):
            tags['lipase_or_hydrolase_relevance'] = 'stereocontrol analogy'
        elif any(k in text for k in ['hydrophobic', 'lipid', 'ester', 'interfacial']):
            tags['lipase_or_hydrolase_relevance'] = 'hydrophobic-substrate analogy'
        else:
            tags['lipase_or_hydrolase_relevance'] = 'no obvious relevance'

    if 'manual_review_priority' in schema:
        if tags.get('photoenzymatic_status') == '已接入photoenzymatic':
            tags['manual_review_priority'] = 'high'
        elif tags.get('enzyme_integration_opportunity') in ('high', 'medium'):
            tags['manual_review_priority'] = 'medium'
        else:
            tags['manual_review_priority'] = 'low'

    if 'confidence' in schema:
        tags['confidence'] = '0.55'
    
    if 'has_mechanistic_study' in schema:
        tags['has_mechanistic_study'] = '是' if len(mechs) > 0 else '否'
    
    # 其他未覆盖字段默认
    for k in schema.keys():
        if k not in tags:
            tags[k] = 'unclear' if k in CLASSIFICATION_FIELDS else 'N/A'
            
    return tags


if __name__ == "__main__":
    main()
