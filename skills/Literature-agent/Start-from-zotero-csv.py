#!/usr/bin/env python3
"""
Zotero Photoenzymatic Radical Knowledge Agent
==============================================
从Zotero导出的CSV出发，构建光酶/自由基文献的知识脉络。

Usage:
    1. Zotero中选中你的20-30篇文献 → 右键 "导出条目" → 格式选CSV → 保存为 zotero_library.csv
    2. python zotero_photoenzyme_agent.py --input zotero_library.csv --output ./knowledge_base/
"""

import os
import sys
import json
import argparse
import re
import time
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


# ==================== 配置区 ====================
LLM_CONFIG = {
    "api_key": os.getenv("KIMI_API_KEY") or os.getenv("MOONSHOT_API_KEY") or "",
    "base_url": os.getenv("MOONSHOT_BASE_URL") or "https://api.moonshot.cn/v1",
    "model": "kimi-k2.6",  
    "temperature": 1,
    "max_tokens": 3600,
    "timeout": 180,
    "json_mode": True
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
        "key_innovation": "核心创新点(中文一句话)",
        "limitation": "局限性/挑战(中文一句话)",
        "future_direction": "未来方向或可迁移启发(中文一句话)",
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


class LLMExtractor:
    """调用LLM API进行文献语义提取"""
    
    def __init__(self, config: Dict, schema: Dict[str, str]):
        self.config = config
        self.schema = schema
        self.headers = {
            "Authorization": f"Bearer {config['api_key']}",
            "Content-Type": "application/json"
        }
        
    def extract(self, title: str, abstract: str, journal: str = "") -> Dict[str, str]:
        """对单篇文献进行结构化提取"""
        if not abstract or len(abstract) < 50:
            # 摘要太短或缺失，返回空标签
            return {k: "N/A" for k in self.schema.keys()}
            
        schema_desc = "\n".join([f"- {k}: {v}" for k, v in self.schema.items()])
        
        prompt = f"""你是一位photoenzymatic catalysis、radical chemistry、organic photochemistry和酶工程交叉领域的文献分析专家。请仔细阅读以下文献信息，并按schema提取结构化知识。

分析优先级：
1. 先判断文章属于photoenzymatic radical、photoredox radical reaction还是更广义的radical chemical reaction。
2. 如果尚未接入酶体系，不要强行归为酶催化；请标记为“仅光化学/自由基反应”或“可设想接入酶体系”。
3. lipase/alpha-beta hydrolase只作为潜在关联或人工筛选线索，不作为主分类前提。

文献标题: {title}
期刊: {journal}
摘要: {abstract}

提取Schema（请严格按JSON格式返回，只返回JSON，不要markdown代码块）:
{schema_desc}

注意：
1. 如果信息不足，填"N/A"或"未明确"
2. key_innovation/limitation/future_direction必须用中文一句话概括
3. 多值字段用" | "分隔，如"C-C | C-N"
4. 必须返回合法JSON，不要任何额外解释"""

        try:
            base_payload = {
                "model": self.config['model'],
                "messages": [{"role": "user", "content": prompt}],
                "temperature": self.config['temperature'],
            }

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
                    timeout=self.config.get("timeout", 60)
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
                    result[k] = "N/A"
            return result
            
        except Exception as e:
            detail = ""
            if "resp" in locals():
                try:
                    detail = f" | 响应: {resp.text[:500]}"
                except Exception:
                    pass
            print(f"  LLM提取失败: {e}{detail}")
            return {k: "N/A" for k in self.schema.keys()}


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
        for _, row in self.df.iterrows():
            title = str(row.get('Title', '')).strip()
            if not title or title == 'nan':
                continue
                
            doi = str(row.get('DOI', '')).strip()
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
                    if val in ('N/A', '未明确', '无', ''):
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
        single_fields = ['enzyme_category', 'enzyme_family', 'radical_mechanism', 
                        'reaction_type', 'stereochemistry', 'has_mechanistic_study',
                        'light_source', 'cofactor']
        for f in single_fields:
            if f in tags and tags[f] not in ('N/A', '未明确', ''):
                concepts[f].append(tags[f])
                
        # 多值字段拆分（如bond_formation可能含"C-C | C-N"）
        multi_fields = ['bond_formation', 'substrate_class', 'mechanistic_methods']
        for f in multi_fields:
            if f in tags and tags[f] not in ('N/A', '未明确', '', '无'):
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
            'substrate_class': '#C9B8A8', 'cofactor': '#E8D5C4'
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
            rt = p.tags.get('reaction_type', 'Unknown').replace(' ', '_')
            ef = p.tags.get('enzyme_family', 'Unknown').replace(' ', '_')
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
                rt = p.tags.get('reaction_type', 'Unknown')
                data.append({'year': p.year, 'reaction': rt, 'title': p.get_display_title(30)})
                
        df = pd.DataFrame(data)
        if df.empty:
            return
            
        # 统计每年各反应类型数量
        pivot = df.groupby(['year', 'reaction']).size().unstack(fill_value=0)
        
        # 堆叠条形图
        colors = [MORANDI_COLORS.get(c.replace(' ', '_'), '#999') for c in pivot.columns]
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
            rt = p.tags.get('reaction_type', 'Unknown')
            ef = p.tags.get('enzyme_family', 'Unknown')
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
            mech = p.tags.get('radical_mechanism', 'Unknown')
            rt = p.tags.get('reaction_type', 'Unknown')
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
        # 2. 按聚类的文献卡片
        self._write_clustered_cards()
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
            ef = p.tags.get('enzyme_family', 'Unknown')
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
            rt = p.tags.get('reaction_type', 'Unknown')
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
            m = p.tags.get('radical_mechanism', 'Unknown')
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
            rt = p.tags.get('reaction_type', 'Uncategorized')
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
                    f"- 酶类别: {p.tags.get('enzyme_category', 'N/A')}",
                    f"- 酶家族: {p.tags.get('enzyme_family', 'N/A')}",
                    f"- 自由基机制: {p.tags.get('radical_mechanism', 'N/A')}",
                    f"- 成键类型: {p.tags.get('bond_formation', 'N/A')}",
                    f"- 底物: {p.tags.get('substrate_class', 'N/A')}",
                    f"- 立体化学: {p.tags.get('stereochemistry', 'N/A')}",
                    f"- 光源: {p.tags.get('light_source', 'N/A')}",
                    f"- 辅因子: {p.tags.get('cofactor', 'N/A')}",
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
        
    def _write_gap_analysis(self):
        """研究空白与机会分析"""
        # 统计存在的组合
        existing = set()
        for p in self.papers:
            rt = p.tags.get('reaction_type', '')
            ef = p.tags.get('enzyme_family', '')
            bond = p.tags.get('bond_formation', '')
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
        
        mechanistic = [p for p in self.papers 
                      if p.tags.get('has_mechanistic_study', '') == '是']
        lines.append(f"- 含机制研究的文献: {len(mechanistic)}/{len(self.papers)} ({len(mechanistic)/len(self.papers)*100:.1f}%)")
        
        methods = Counter()
        for p in mechanistic:
            m = p.tags.get('mechanistic_methods', '')
            for method in re.split(r'[|;/,、]', m):
                methods[method.strip()] += 1
                
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
        directions = [p.tags.get('future_direction', '') for p in self.papers 
                     if p.tags.get('future_direction', '') not in ('N/A', '未明确', '')]
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
    parser.add_argument('--temperature', type=float, default=LLM_CONFIG['temperature'], help='结构化提取建议0.0-0.3')
    parser.add_argument('--max-tokens', type=int, default=LLM_CONFIG['max_tokens'], help='每篇文献JSON输出的最大token数')
    parser.add_argument('--no-json-mode', action='store_true', help='禁用response_format json_object，部分模型/端点不支持时使用')
    parser.add_argument('--schema-profile', choices=sorted(SCHEMA_PROFILES.keys()), default='photo-radical', help='内置提取schema')
    parser.add_argument('--schema-file', help='外部JSON schema文件；提供后覆盖 --schema-profile')
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    schema = load_schema(args.schema_profile, args.schema_file)
    
    # 1. 导入
    print("=" * 60)
    print("STEP 1: 从Zotero CSV导入文献")
    print("=" * 60)
    importer = ZoteroImporter(args.input, args.csv_encoding)
    papers = importer.parse()
    
    if not papers:
        print("未找到有效文献，请检查CSV文件")
        return
        
    # 2. LLM语义提取（可选但强烈推荐）
    if args.use_llm:
        print("\n" + "=" * 60)
        print("STEP 2: LLM语义提取（每篇约3-5秒）")
        print("=" * 60)
        
        if args.api_key:
            LLM_CONFIG['api_key'] = args.api_key
        LLM_CONFIG['base_url'] = args.base_url
        LLM_CONFIG['model'] = args.model
        LLM_CONFIG['temperature'] = args.temperature
        LLM_CONFIG['max_tokens'] = args.max_tokens
        LLM_CONFIG['json_mode'] = not args.no_json_mode
            
        if not LLM_CONFIG['api_key']:
            print("警告: 未配置API Key，请设置环境变量 KIMI_API_KEY/MOONSHOT_API_KEY 或使用 --api-key 参数")
            print("跳过LLM提取，使用基于关键词的fallback...")
            # Fallback: 基于关键词的正则匹配
            for p in papers:
                p.tags = fallback_extract(p.title, p.abstract, schema)
        else:
            extractor = LLMExtractor(LLM_CONFIG, schema)
            for i, p in enumerate(papers, 1):
                print(f"  [{i}/{len(papers)}] 提取: {p.get_display_title(50)}")
                p.tags = extractor.extract(p.title, p.abstract, p.journal)
                time.sleep(0.5)  # 避免rate limit
    else:
        print("\n" + "=" * 60)
        print("STEP 2: 使用本地正则提取（无需API）")
        print("=" * 60)
        for p in papers:
            p.tags = fallback_extract(p.title, p.abstract, schema)
            
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
    
    print("\n" + "=" * 60)
    print("完成！输出文件:")
    print("=" * 60)
    for f in sorted(output_dir.iterdir()):
        print(f"  - {f.name}")


def fallback_extract(title: str, abstract: str, schema: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """无需LLM时的正则fallback提取"""
    schema = schema or EXTRACTION_SCHEMA
    text = (title + " " + abstract).lower()
    tags = {}
    
    # 酶家族
    if any(k in text for k in ['p450', 'cytochrome']):
        tags['enzyme_family'] = 'P450'
    elif any(k in text for k in ['ered', 'ene-reductase', 'old yellow enzyme', 'oye']):
        tags['enzyme_family'] = 'ERED'
    elif any(k in text for k in ['fap', 'fatty acid photodecarboxylase', 'photodecarboxylase']):
        tags['enzyme_family'] = 'FAP'
    elif any(k in text for k in ['aldolase', 'class i aldolase']):
        tags['enzyme_family'] = 'Aldolase'
    elif any(k in text for k in ['lipase', 'esterase']):
        tags['enzyme_family'] = 'Lipase'
    else:
        tags['enzyme_family'] = 'Other_enzyme' if 'enzyme_family' in schema else 'N/A'
        
    # 反应类型
    if 'decarboxyl' in text:
        tags['reaction_type'] = 'decarboxylative_coupling'
    elif 'arylat' in text or 'c-h arylation' in text:
        tags['reaction_type'] = 'C-H_arylation'
    elif 'alkylat' in text or 'c-h alkylation' in text:
        tags['reaction_type'] = 'C-H_alkylation'
    elif 'halogen' in text or 'chlorinat' in text or 'brominat' in text:
        tags['reaction_type'] = 'halogenation'
    elif 'cycloaddition' in text or 'cycloadd' in text:
        tags['reaction_type'] = 'cycloaddition'
    elif 'giese' in text:
        tags['reaction_type'] = 'Giese_addition'
    else:
        tags['reaction_type'] = 'other'
        
    # 机制
    if 'photoredox' in text or 'photo-redox' in text:
        tags['radical_mechanism'] = '光氧化还原(PET)'
        tags['radical_or_excited_state_mechanism'] = 'PET/photoredox'
    elif 'energy transfer' in text or 'triplet' in text:
        tags['radical_mechanism'] = '能量转移(EnT)'
        tags['radical_or_excited_state_mechanism'] = 'EnT'
    elif 'hat' in text or 'hydrogen atom transfer' in text:
        tags['radical_mechanism'] = '氢原子转移(HAT)'
        tags['radical_or_excited_state_mechanism'] = 'HAT'
    elif 'eda' in text:
        tags['radical_mechanism'] = 'EDA复合物'
        tags['radical_or_excited_state_mechanism'] = 'EDA complex'
    else:
        tags['radical_mechanism'] = '直接激发辅因子'
        tags['radical_or_excited_state_mechanism'] = 'unclear'

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

    if 'photoenzymatic_status' in schema:
        if any(k in text for k in ['photoenzym', 'photobiocatal']):
            tags['photoenzymatic_status'] = '已接入photoenzymatic'
        elif any(k in text for k in ['enzyme', 'enzymatic', 'biocatal']):
            tags['photoenzymatic_status'] = '仅酶催化背景'
        elif any(k in text for k in ['photoredox', 'visible light', 'photoinduced', 'radical']):
            tags['photoenzymatic_status'] = '仅光化学/自由基反应'
        else:
            tags['photoenzymatic_status'] = '未明确'

    if 'enzyme_integration_opportunity' in schema:
        if any(k in text for k in ['enantioselective', 'stereoselective', 'site-selective', 'regioselective', 'drug', 'natural product']):
            tags['enzyme_integration_opportunity'] = 'medium'
        elif any(k in text for k in ['enzyme', 'active site', 'cofactor', 'flavin', 'fad', 'fmn']):
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
        tags['confidence'] = '0.35'
        
    # 其他字段默认
    for k in schema.keys():
        if k not in tags:
            tags[k] = 'N/A'
            
    return tags


if __name__ == "__main__":
    main()
