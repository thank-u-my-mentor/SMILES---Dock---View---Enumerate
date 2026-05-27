#python tree_pipeline.py representatives_aln.fasta
import json
import csv
import os
import subprocess
import shutil
import re
import time
from pathlib import Path
from typing import List, Optional, Set, Dict
from collections import defaultdict
from dataclasses import dataclass
import argparse
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

IQTREE = shutil.which("iqtree") or shutil.which("iqtree2")

class Config:
    BASE_DIR = Path(".")
    OUTPUT_DIR = BASE_DIR / "tree_analysis"
    CLUSTALO = shutil.which("clustalo") or "/mnt/l/WSL/softwares/conda_envs/md/bin/clustalo"
    FASTTREE = shutil.which("FastTree") or shutil.which("fasttree")
    KIMI_API_KEY = os.getenv("KIMI_API_KEY")
    # 不再强制要求，因为 Tree pipeline 直接从 SSN CSV 读 taxonomy
    if not KIMI_API_KEY:
        print("  ⚠ 警告: KIMI_API_KEY 未设置，Unknown 序列将无法 Kimi 补全")
    KIMI_MODEL = "moonshot-v1-128k"
    MSA_THREADS = 16

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
    is_core: bool = False
    community: int = -1

    def __post_init__(self):
        self.length = len(self.sequence)


class PhylogeneticTreePipeline:
    def __init__(self, input_fasta: Path, output_dir: Optional[Path] = None,
                 core_fasta: Optional[Path] = None):
        self.input_fasta = Path(input_fasta)
        self.output_dir = output_dir or Config.OUTPUT_DIR
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.core_fasta = Path(core_fasta) if core_fasta else None

        # 全局缓存（Kimi 补全备用）
        self.global_cache_file = Config.BASE_DIR.parent / "kingdom_global_cache.json"
        self.kingdom_cache: Dict = {}
        self.load_global_kingdom_cache()
        self.sequences: List[SequenceRecord] = []

    def _load_core_ids(self) -> Set[str]:
        if not self.core_fasta or not self.core_fasta.exists():
            return set()
        core_ids = set()
        with open(self.core_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seq_id = line[1:].split()[0]
                    core_ids.add(seq_id)
        print(f"  ✓ 加载核心序列 {len(core_ids)} 条")
        return core_ids

    def load_global_kingdom_cache(self):
        if self.global_cache_file.exists():
            try:
                with open(self.global_cache_file, 'r', encoding='utf-8') as f:
                    raw = json.load(f)
                self.kingdom_cache = {str(k).strip(): v for k, v in raw.items()}
                print(f"  ✓ 已加载全局超级字典 ({len(self.kingdom_cache)} 条)")
            except Exception:
                self.kingdom_cache = {}

    def save_global_kingdom_cache(self):
        try:
            with open(self.global_cache_file, 'w', encoding='utf-8') as f:
                json.dump(self.kingdom_cache, f, indent=2, ensure_ascii=False)
        except Exception:
            pass

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

    def _get_local_dict(self, organism_name: str) -> dict:
        if not organism_name or organism_name == "Unknown":
            return {"kingdom": "Unknown", "phylum_class": "Unknown"}
        std_name = self._standardize_organism_name(organism_name)
        cache_key = f"org:{std_name}"
        if cache_key in self.kingdom_cache:
            return self.kingdom_cache[cache_key]
        return {"kingdom": "Unknown", "phylum_class": "Unknown"}

    def _get_kimi_kingdom(self, organism_name: str) -> dict:
        if not organism_name or organism_name == "Unknown":
            return {"kingdom": "Unknown", "phylum_class": "Unknown"}

        std_name = self._standardize_organism_name(organism_name)
        cache_key = f"org:{std_name}"

        lower_org = organism_name.lower()
        if "metagenome" in lower_org:
            result = {"kingdom": "Metagenome", "phylum_class": "Environmental"}
            self.kingdom_cache[cache_key] = result
            self.save_global_kingdom_cache()
            print(f"    🌊 宏基因组标记: {organism_name}")
            return result

        local = self._get_local_dict(organism_name)
        if local["kingdom"] != "Unknown":
            return local

        for attempt in range(4):
            try:
                prompt = f"""你是一个专业的生物分类与文献专家。
物种名称: {organism_name}
请严格返回以下JSON格式，不要任何额外文字：
{{
  "kingdom": "Plant" | "Animal" | "Fungi" | "Bacteria" | "Archaea" | "Protist" | "Unknown",
  "phylum_class": "具体门或纲名称，如果不确定填 Unknown",
  "literature": "DOI或PubMed链接，如果没有可填空字符串"
}}"""
                response = requests.post(
                    "https://api.moonshot.cn/v1/chat/completions",
                    headers={"Authorization": f"Bearer {Config.KIMI_API_KEY}", "Content-Type": "application/json"},
                    json={"model": Config.KIMI_MODEL, "messages": [{"role": "user", "content": prompt}],
                          "temperature": 0.1, "max_tokens": 400},
                    timeout=40
                )
                if response.status_code == 429:
                    time.sleep(2 ** attempt + 1)
                    continue
                if response.status_code != 200:
                    return {"kingdom": "Unknown", "phylum_class": "Unknown"}

                content = response.json()["choices"][0]["message"]["content"].strip()
                if content.startswith("```"):
                    content = content.split("```")[1].replace("json", "").strip()
                info = json.loads(content)

                result = {
                    "organism": organism_name,
                    "kingdom": info.get("kingdom", "Unknown"),
                    "phylum_class": info.get("phylum_class", "Unknown")
                }
                if result["kingdom"] != "Unknown":
                    self.kingdom_cache[cache_key] = result
                    self.save_global_kingdom_cache()
                return result
            except Exception:
                if attempt < 3:
                    time.sleep(1.5)
        return {"kingdom": "Unknown", "phylum_class": "Unknown"}

    def parse_fasta(self) -> List[SequenceRecord]:
        """从 SSN 输出直接读取 taxonomy 和 core 标记"""
        print(f"\n📥 解析同源 FASTA: {self.input_fasta}")
        
        # 加载 SSN taxonomy + core 标记
        ssn_nodes_file = Path("ssn_nodes.csv")
        ssn_data = {}
        if ssn_nodes_file.exists():
            with open(ssn_nodes_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    ssn_data[row['id']] = {
                        'kingdom': row.get('kingdom', 'Unknown'),
                        'organism': row.get('organism', 'Unknown'),
                        'phylum_class': row.get('phylum_class', 'Unknown'),
                        'is_core': row.get('is_core', 'FALSE').upper() == 'TRUE'
                    }
            print(f"  ✓ 加载 SSN 数据 ({len(ssn_data)} 条)")
        else:
            print(f"  ⚠ 未找到 {ssn_nodes_file}")

        records = []
        current_id = None
        current_seq = []

        with open(self.input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id and current_seq:
                        rec = self._create_record(current_id, ''.join(current_seq), ssn_data)
                        records.append(rec)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_id and current_seq:
            rec = self._create_record(current_id, ''.join(current_seq), ssn_data)
            records.append(rec)

        # 统计
        kingdom_stats = defaultdict(int)
        core_count = 0
        for r in records:
            kingdom_stats[r.kingdom] += 1
            if r.is_core:
                core_count += 1
        print(f"  📊 Kingdom 分布: {dict(kingdom_stats)}")
        print(f"  ⭐ Core 序列: {core_count} 条")
        
        unknown_count = sum(1 for r in records if r.kingdom == 'Unknown')
        print(f"  🔍 Unknown 数量: {unknown_count}")
        
        # Kimi 补全（仅对 Unknown）
        if unknown_count > 0:
            print(f"\n🤖 对 {unknown_count} 个 Unknown 进行 Kimi 补全...")
            organism_to_records = defaultdict(list)
            for rec in records:
                if rec.kingdom == "Unknown" and rec.organism != "Unknown":
                    organism_to_records[rec.organism].append(rec)
            
            for org, recs in organism_to_records.items():
                info = self._get_kimi_kingdom(org)
                for rec in recs:
                    rec.kingdom = info["kingdom"]
                    rec.phylum_class = info["phylum_class"]

        return records
        
    def _create_record(self, header: str, sequence: str, ssn_data: dict) -> SequenceRecord:
        raw_id = header[1:].split()[0] if header.startswith('>') else header.split()[0]
        uniprot_id = self._clean_uniprot_id(raw_id)
        
        # 从 SSN 数据读取
        data = ssn_data.get(raw_id, {})
        kingdom = data.get('kingdom', 'Unknown')
        organism = data.get('organism', 'Unknown')
        phylum_class = data.get('phylum_class', 'Unknown')
        is_core = data.get('is_core', False)
        
        if is_core:
            print(f"    ⭐ Core 序列: {raw_id}")
        
        return SequenceRecord(
            id=raw_id,
            uniprot_id=uniprot_id,
            species_code="UNKNOWN",
            sequence=sequence,
            kingdom=kingdom,
            organism=organism,
            phylum_class=phylum_class,
            is_core=is_core
        )

    def run_msa(self, fasta_file: Path, aln_file: Path) -> bool:
        if not Path(Config.CLUSTALO).exists():
            print("  ⚠ ClustalO 未找到")
            return False
        cmd = [Config.CLUSTALO, "-i", str(fasta_file), "-o", str(aln_file),
               "--outfmt", "fasta", "--force", "--threads", str(Config.MSA_THREADS)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ⚠ ClustalO 失败: {result.stderr[:500]}")
            return False
        return True

    def build_tree(self, aln_file: Path, tree_file: Path) -> bool:
        if not IQTREE:
            return False
        print("\n🌳 使用 IQ-TREE 构建最大似然树...")
        # -m LG+F+R4 是蛋白序列推荐模型（Le+G4+F+R）
        # -B 1000 是 ultrafast bootstrap
        # -T AUTO 自动检测线程
        cmd = [
            IQTREE, "-s", str(aln_file),
            "-m", "LG+F+R4",
            "-B", "1000",
            "-T", "16",
            "--prefix", str(self.output_dir / "iqtree_run")
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            # IQ-TREE 输出到 prefix.treefile
            src = self.output_dir / "iqtree_run.treefile"
            if src.exists():
                shutil.copy(str(src), str(tree_file))
                print(f"  ✅ 系统发育树生成完成: {tree_file}")
                return True
        print(f"  ⚠ IQ-TREE 失败: {result.stderr[:500]}")
        return False

    def export_itol_files(self):
        """生成 iTOL 交互式着色文件（柔和 pastel 版）"""
        # Pastel 柔和配色：饱和度低，白色背景下更高级
        colors = {
            'Plant': '#A8D5BA',      # 淡薄荷
            'Animal': '#8FB8E6',     # 淡天蓝
            'Fungi': '#F4C2A1',      # 淡杏
            'Bacteria': '#E8A0B0',   # 淡玫瑰（避免刺眼大红）
            'Archaea': '#C5A5CF',    # 淡紫
            'Protist': '#F9E79F',    # 淡黄
            'Metagenome': '#D5D5D5', # 浅灰
            'Unknown': '#AAAAAA',    # 中灰
        }

        color_list = "\t".join(colors.values())
        label_list = "\t".join(colors.keys())
        shape_list = "\t".join(["1"] * len(colors))

        # 1. Kingdom 颜色条 —— 去掉所有非法参数，关闭 strip 文字
        with open(self.output_dir / "itol_kingdom_color_strip.txt", 'w', encoding='utf-8') as f:
            f.write("DATASET_COLORSTRIP\n")
            f.write("SEPARATOR TAB\n")
            f.write("DATASET_LABEL\tKingdom\n")
            f.write("COLOR\t#000000\n")
            f.write("STRIP_WIDTH\t25\n")
            f.write("MARGIN\t5\n")
            f.write("BORDER_WIDTH\t1\n")
            f.write("BORDER_COLOR\t#000000\n")
            f.write("SHOW_STRIP_LABELS\t0\n")        # ← 关键：不显示条带文字，避免重复堆叠
            f.write("LEGEND_TITLE\tKingdom\n")
            f.write(f"LEGEND_SHAPES\t{shape_list}\n")
            f.write(f"LEGEND_COLORS\t{color_list}\n")
            f.write(f"LEGEND_LABELS\t{label_list}\n")
            f.write("DATA\n")
            for rec in self.sequences:
                color = colors.get(rec.kingdom, '#AAAAAA')
                # 第三列 label 仍要填，但 SHOW_STRIP_LABELS=0 时不会画在环上，只影响 legend
                f.write(f"{rec.id}\t{color}\t{rec.kingdom}\n")

        # 2. Phylum/Class 标签（保持不变）
        with open(self.output_dir / "itol_phylum_label.txt", 'w', encoding='utf-8') as f:
            f.write("LABELS\n")
            f.write("SEPARATOR TAB\n")
            f.write("DATA\n")
            for rec in self.sequences:
                f.write(f"{rec.id}\t{rec.phylum_class}\n")

        # 3. 核心序列高亮（BINARY 点状标记）
        core_records = [r for r in self.sequences if r.is_core]
        if core_records:
            with open(self.output_dir / "itol_core_highlight.txt", 'w', encoding='utf-8') as f:
                f.write("DATASET_BINARY\n")
                f.write("SEPARATOR TAB\n")
                f.write("DATASET_LABEL\tCore Sequences\n")
                f.write("COLOR\t#FF0000\n")
                f.write("FIELD_SHAPES\t2\n")        # 2 = 圆形
                f.write("FIELD_LABELS\tcore\n")
                f.write("DATA\n")
                for rec in core_records:
                    f.write(f"{rec.id}\t1\n")

        print("  ✓ iTOL 文件已生成（pastel 配色 + 无重复标签）")
            
    def run_pipeline(self):
        print(f"\n{'#'*80}\n# Phylogenetic Tree Pipeline - SSN 代表序列专用版\n{'#'*80}\n")
        self.sequences = self.parse_fasta()

        seq_fasta = self.output_dir / "input_sequences.fasta"
        with open(seq_fasta, 'w') as f:
            for rec in self.sequences:
                f.write(f">{rec.id}\n{rec.sequence}\n")

        aln_file = self.output_dir / "alignment.fasta"
        msa_ok = self.run_msa(seq_fasta, aln_file)
        
        # 关键修复：MSA 失败直接报错退出，不再回退到未对齐序列
        if not msa_ok:
            raise RuntimeError(
                "ClustalO MSA 失败，无法构建可靠系统发育树。\n"
                "建议：1) 检查 ClustalO 安装；2) 改用 MAFFT: mafft --auto input.fa > aln.fa"
            )

        tree_file = self.output_dir / "phylogenetic_tree.newick"
        self.build_tree(aln_file, tree_file)

        self.export_itol_files()

        print(f"\n{'='*80}\n✅ 交互式系统发育树生成完成！\n"
              f"输出目录: {self.output_dir}\n"
              f"文件列表：\n"
              f"  - alignment.fasta (MSA 结果)\n"
              f"  - phylogenetic_tree.newick (带 bootstrap 的树)\n"
              f"  - itol_kingdom_color_strip.txt (Kingdom 着色)\n"
              f"  - itol_phylum_label.txt (门/纲标签)\n"
              f"  - itol_core_highlight.txt (核心序列高亮)\n"
              f"\n使用 iTOL (https://itol.embl.de/) 上传 .newick + 上述 txt 文件\n"
              f"即可获得带物种分类分色和 bootstrap 支持度的交互式系统发育树\n{'='*80}")


def main():
    parser = argparse.ArgumentParser(description='CHR 系统发育树流水线')
    parser.add_argument('input', help='输入 FASTA 文件（SSN 的 representatives.fasta）')
    parser.add_argument('--output', default=None, help='输出目录')
    args = parser.parse_args()

    pipeline = PhylogeneticTreePipeline(
        input_fasta=Path(args.input),
        output_dir=Path(args.output) if args.output else None
    )
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()
