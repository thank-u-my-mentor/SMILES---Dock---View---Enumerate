from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


END_OF_CHAIN = 0xFFFFFFFE
FREE_SECTOR = 0xFFFFFFFF
FAT_SECTOR = 0xFFFFFFFD
DIF_SECTOR = 0xFFFFFFFC


def u16(data: bytes, offset: int) -> int:
    return struct.unpack_from("<H", data, offset)[0]


def u32(data: bytes, offset: int) -> int:
    return struct.unpack_from("<I", data, offset)[0]


def i32(data: bytes, offset: int) -> int:
    return struct.unpack_from("<i", data, offset)[0]


@dataclass
class OleDirEntry:
    name: str
    obj_type: int
    start_sector: int
    size: int


class OleFile:
    """Minimal Compound File Binary reader sufficient for BIFF8 .xls files."""

    def __init__(self, path: Path):
        self.path = path
        self.data = path.read_bytes()
        if self.data[:8] != bytes.fromhex("D0 CF 11 E0 A1 B1 1A E1"):
            raise ValueError(f"{path} is not an OLE compound file")
        self.sector_size = 1 << u16(self.data, 30)
        self.mini_sector_size = 1 << u16(self.data, 32)
        self.num_fat_sectors = u32(self.data, 44)
        self.first_dir_sector = u32(self.data, 48)
        self.mini_cutoff = u32(self.data, 56)
        self.first_minifat_sector = u32(self.data, 60)
        self.num_minifat_sectors = u32(self.data, 64)
        self.first_difat_sector = u32(self.data, 68)
        self.num_difat_sectors = u32(self.data, 72)
        self.difat = self._load_difat()
        self.fat = self._load_fat()
        self.directory = self._load_directory()
        self.root = self.directory[0] if self.directory else None
        self.minifat = self._load_minifat()
        self.mini_stream = b""
        if self.root and self.root.start_sector not in (FREE_SECTOR, END_OF_CHAIN):
            self.mini_stream = self._read_regular_stream(self.root.start_sector, self.root.size)

    def sector(self, sid: int) -> bytes:
        start = (sid + 1) * self.sector_size
        return self.data[start : start + self.sector_size]

    def _load_difat(self) -> List[int]:
        entries: List[int] = []
        for i in range(109):
            sid = u32(self.data, 76 + 4 * i)
            if sid not in (FREE_SECTOR, END_OF_CHAIN):
                entries.append(sid)
        sid = self.first_difat_sector
        for _ in range(self.num_difat_sectors):
            if sid in (FREE_SECTOR, END_OF_CHAIN):
                break
            sec = self.sector(sid)
            for i in range((self.sector_size // 4) - 1):
                val = u32(sec, 4 * i)
                if val not in (FREE_SECTOR, END_OF_CHAIN):
                    entries.append(val)
            sid = u32(sec, self.sector_size - 4)
        return entries[: self.num_fat_sectors]

    def _load_fat(self) -> List[int]:
        fat: List[int] = []
        for sid in self.difat:
            sec = self.sector(sid)
            fat.extend(struct.unpack("<" + "I" * (self.sector_size // 4), sec))
        return fat

    def _chain(self, start: int, fat: Optional[List[int]] = None) -> List[int]:
        if fat is None:
            fat = self.fat
        chain: List[int] = []
        sid = start
        seen = set()
        while sid not in (END_OF_CHAIN, FREE_SECTOR) and sid < len(fat):
            if sid in seen:
                raise ValueError("OLE sector chain loop detected")
            seen.add(sid)
            chain.append(sid)
            sid = fat[sid]
        return chain

    def _read_regular_stream(self, start: int, size: int) -> bytes:
        chunks = [self.sector(sid) for sid in self._chain(start)]
        return b"".join(chunks)[:size]

    def _load_directory(self) -> List[OleDirEntry]:
        raw = self._read_regular_stream(self.first_dir_sector, 1 << 30)
        entries: List[OleDirEntry] = []
        for offset in range(0, len(raw), 128):
            ent = raw[offset : offset + 128]
            if len(ent) < 128:
                break
            name_len = u16(ent, 64)
            if name_len >= 2:
                name = ent[: name_len - 2].decode("utf-16le", errors="ignore")
            else:
                name = ""
            obj_type = ent[66]
            start_sector = u32(ent, 116)
            size = struct.unpack_from("<Q", ent, 120)[0]
            entries.append(OleDirEntry(name, obj_type, start_sector, int(size)))
        return entries

    def _load_minifat(self) -> List[int]:
        if self.first_minifat_sector in (FREE_SECTOR, END_OF_CHAIN) or self.num_minifat_sectors == 0:
            return []
        raw = b"".join(self.sector(sid) for sid in self._chain(self.first_minifat_sector))
        entries = list(struct.unpack("<" + "I" * (len(raw) // 4), raw))
        return entries

    def streams(self) -> List[str]:
        return [e.name for e in self.directory if e.obj_type == 2]

    def open_stream(self, name: str) -> bytes:
        wanted = name.lower()
        for entry in self.directory:
            if entry.obj_type == 2 and entry.name.lower() == wanted:
                if entry.size < self.mini_cutoff and self.minifat:
                    return self._read_mini_stream(entry.start_sector, entry.size)
                return self._read_regular_stream(entry.start_sector, entry.size)
        raise KeyError(name)

    def _read_mini_stream(self, start: int, size: int) -> bytes:
        chunks = []
        for sid in self._chain(start, self.minifat):
            begin = sid * self.mini_sector_size
            chunks.append(self.mini_stream[begin : begin + self.mini_sector_size])
        return b"".join(chunks)[:size]


@dataclass
class BiffRecord:
    record_id: int
    payload: bytes
    offset: int


def iter_biff_records(data: bytes, start: int = 0) -> Iterable[BiffRecord]:
    pos = start
    limit = len(data)
    while pos + 4 <= limit:
        rid = u16(data, pos)
        length = u16(data, pos + 2)
        payload_start = pos + 4
        payload_end = payload_start + length
        if payload_end > limit:
            break
        yield BiffRecord(rid, data[payload_start:payload_end], pos)
        pos = payload_end


class SegmentedCursor:
    def __init__(self, segments: List[bytes]):
        self.segments = segments
        self.seg_idx = 0
        self.pos = 0

    def _advance_segment(self) -> None:
        self.seg_idx += 1
        self.pos = 0

    def read(self, n: int) -> bytes:
        out = bytearray()
        while n > 0 and self.seg_idx < len(self.segments):
            seg = self.segments[self.seg_idx]
            avail = len(seg) - self.pos
            if avail <= 0:
                self._advance_segment()
                continue
            take = min(n, avail)
            out.extend(seg[self.pos : self.pos + take])
            self.pos += take
            n -= take
        if n:
            raise EOFError("unexpected end of SST data")
        return bytes(out)

    def read_chars(self, char_count: int, is_16bit: bool) -> str:
        parts: List[str] = []
        remaining = char_count
        width = 2 if is_16bit else 1
        while remaining > 0:
            if self.seg_idx >= len(self.segments):
                raise EOFError("unexpected end of SST string")
            seg = self.segments[self.seg_idx]
            avail = len(seg) - self.pos
            if avail <= 0:
                self._advance_segment()
                if self.seg_idx >= len(self.segments):
                    raise EOFError("unexpected end of SST continuation")
                flags = self.read(1)[0]
                width = 2 if (flags & 0x01) else 1
                continue
            chars_here = min(remaining, avail // width)
            if chars_here == 0:
                # Very defensive handling for an odd byte at a BIFF boundary.
                raw = bytearray(self.read(avail))
                self._advance_segment()
                flags = self.read(1)[0]
                width = 2 if (flags & 0x01) else 1
                needed = width - len(raw)
                raw.extend(self.read(needed))
                parts.append(bytes(raw).decode("utf-16le" if width == 2 else "latin1", errors="replace"))
                remaining -= 1
                continue
            raw = self.read(chars_here * width)
            parts.append(raw.decode("utf-16le" if width == 2 else "latin1", errors="replace"))
            remaining -= chars_here
        return "".join(parts)


def parse_sst(records: List[BiffRecord]) -> List[str]:
    sst_index = next((i for i, rec in enumerate(records) if rec.record_id == 0x00FC), None)
    if sst_index is None:
        return []
    first = records[sst_index].payload
    if len(first) < 8:
        return []
    unique_count = u32(first, 4)
    segments = [first[8:]]
    i = sst_index + 1
    while i < len(records) and records[i].record_id == 0x003C:
        segments.append(records[i].payload)
        i += 1
    cur = SegmentedCursor(segments)
    strings: List[str] = []
    for _ in range(unique_count):
        try:
            cch = u16(cur.read(2), 0)
            flags = cur.read(1)[0]
            is_16bit = bool(flags & 0x01)
            has_ext = bool(flags & 0x04)
            has_rich = bool(flags & 0x08)
            rich_runs = u16(cur.read(2), 0) if has_rich else 0
            ext_size = u32(cur.read(4), 0) if has_ext else 0
            text = cur.read_chars(cch, is_16bit)
            if rich_runs:
                cur.read(4 * rich_runs)
            if ext_size:
                cur.read(ext_size)
            strings.append(text)
        except EOFError:
            break
    return strings


def decode_biff_string(raw: bytes, offset: int = 0) -> str:
    if offset + 3 > len(raw):
        return ""
    cch = u16(raw, offset)
    flags = raw[offset + 2]
    pos = offset + 3
    has_ext = bool(flags & 0x04)
    has_rich = bool(flags & 0x08)
    if has_rich:
        pos += 2
    if has_ext:
        pos += 4
    width = 2 if (flags & 0x01) else 1
    text_bytes = raw[pos : pos + cch * width]
    return text_bytes.decode("utf-16le" if width == 2 else "latin1", errors="replace")


def decode_rk(raw_rk: int) -> float:
    mult = raw_rk & 0x01
    is_int = raw_rk & 0x02
    if is_int:
        value = raw_rk >> 2
        if value & (1 << 29):
            value -= 1 << 30
        result = float(value)
    else:
        packed = struct.pack("<Q", (raw_rk & 0xFFFFFFFC) << 32)
        result = struct.unpack("<d", packed)[0]
    if mult:
        result /= 100.0
    return result


def clean_cell(value: Any) -> Any:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def read_xls(path: Path) -> Dict[str, List[List[Any]]]:
    ole = OleFile(path)
    stream_name = "Workbook" if "Workbook" in ole.streams() else "Book"
    workbook = ole.open_stream(stream_name)
    records = list(iter_biff_records(workbook))
    sst = parse_sst(records)

    sheets: List[Tuple[str, int]] = []
    for rec in records:
        if rec.record_id != 0x0085 or len(rec.payload) < 8:
            continue
        sheet_offset = u32(rec.payload, 0)
        name_len = rec.payload[6]
        name_flags = rec.payload[7]
        width = 2 if (name_flags & 0x01) else 1
        name = rec.payload[8 : 8 + name_len * width].decode(
            "utf-16le" if width == 2 else "latin1", errors="replace"
        )
        sheets.append((name, sheet_offset))

    result: Dict[str, List[List[Any]]] = {}
    for sheet_name, sheet_offset in sheets:
        cells: Dict[Tuple[int, int], Any] = {}
        pending_formula_cell: Optional[Tuple[int, int]] = None
        for rec in iter_biff_records(workbook, sheet_offset):
            rid = rec.record_id
            p = rec.payload
            if rid == 0x000A:
                break
            if rid == 0x00FD and len(p) >= 10:  # LABELSST
                row, col = u16(p, 0), u16(p, 2)
                idx = u32(p, 6)
                cells[(row, col)] = sst[idx] if idx < len(sst) else f"#SST[{idx}]"
            elif rid == 0x0204 and len(p) >= 8:  # LABEL
                row, col = u16(p, 0), u16(p, 2)
                cells[(row, col)] = decode_biff_string(p, 6)
            elif rid == 0x0203 and len(p) >= 14:  # NUMBER
                row, col = u16(p, 0), u16(p, 2)
                cells[(row, col)] = clean_cell(struct.unpack_from("<d", p, 6)[0])
            elif rid == 0x027E and len(p) >= 10:  # RK
                row, col = u16(p, 0), u16(p, 2)
                cells[(row, col)] = clean_cell(decode_rk(u32(p, 6)))
            elif rid == 0x00BD and len(p) >= 6:  # MULRK
                row, first_col = u16(p, 0), u16(p, 2)
                last_col = u16(p, len(p) - 2)
                pos = 4
                for col in range(first_col, last_col + 1):
                    if pos + 6 > len(p) - 2:
                        break
                    rk = u32(p, pos + 2)
                    cells[(row, col)] = clean_cell(decode_rk(rk))
                    pos += 6
            elif rid == 0x0201 and len(p) >= 8:  # BLANK
                row, col = u16(p, 0), u16(p, 2)
                cells.setdefault((row, col), "")
            elif rid == 0x00BE and len(p) >= 6:  # MULBLANK
                row, first_col = u16(p, 0), u16(p, 2)
                last_col = u16(p, len(p) - 2)
                for col in range(first_col, last_col + 1):
                    cells.setdefault((row, col), "")
            elif rid == 0x0006 and len(p) >= 14:  # FORMULA
                row, col = u16(p, 0), u16(p, 2)
                raw_result = p[6:14]
                pending_formula_cell = None
                if raw_result[:2] == b"\xff\xff":
                    pending_formula_cell = (row, col)
                else:
                    try:
                        cells[(row, col)] = clean_cell(struct.unpack("<d", raw_result)[0])
                    except struct.error:
                        pass
            elif rid == 0x0207 and pending_formula_cell is not None:  # STRING
                cells[pending_formula_cell] = decode_biff_string(p, 0)
                pending_formula_cell = None

        if cells:
            max_row = max(r for r, _ in cells)
            max_col = max(c for _, c in cells)
            table: List[List[Any]] = []
            for r in range(max_row + 1):
                row_values = [cells.get((r, c), "") for c in range(max_col + 1)]
                while row_values and row_values[-1] == "":
                    row_values.pop()
                table.append(row_values)
            result[sheet_name] = table
        else:
            result[sheet_name] = []
    return result


def dump_xls(path: Path, max_rows: int) -> None:
    sheets = read_xls(path)
    for name, rows in sheets.items():
        print(f"## sheet={name!r} rows={len(rows)}")
        for i, row in enumerate(rows[:max_rows], start=1):
            line = f"{i:03d}: " + "\t".join("" if v is None else str(v) for v in row)
            print(line.encode("utf-8", errors="replace").decode("utf-8"))


def normalize_text(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("\xa0", " ").strip()


def extract_order_mapping(path: Path) -> List[Dict[str, str]]:
    sheets = read_xls(path)
    rows = next(iter(sheets.values()))
    header_idx = None
    for i, row in enumerate(rows):
        norm = [normalize_text(v) for v in row]
        if "*基因名称" in norm and "uniprot ID" in norm:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find order-table header row with *基因名称 and uniprot ID")
    header = [normalize_text(v) for v in rows[header_idx]]
    gene_col = header.index("*基因名称")
    uniprot_col = header.index("uniprot ID")
    number_col = header.index("#") if "#" in header else 0

    accession_re = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[A-Z0-9]{7}|[A-Z0-9]{6,10})$")
    mapping: List[Dict[str, str]] = []
    seen = set()
    for idx, row in enumerate(rows[header_idx + 1 :], start=header_idx + 2):
        gene = normalize_text(row[gene_col] if gene_col < len(row) else "")
        uniprot = normalize_text(row[uniprot_col] if uniprot_col < len(row) else "")
        serial = normalize_text(row[number_col] if number_col < len(row) else "")
        if not gene or not uniprot:
            continue
        uniprot = re.split(r"[\s,;]+", uniprot)[0].strip()
        if not accession_re.match(uniprot):
            continue
        key = (uniprot, gene)
        if key in seen:
            continue
        seen.add(key)
        mapping.append({"excel_row": str(idx), "serial": serial, "gene_name": gene, "uniprot": uniprot})
    return mapping


def print_order_mapping(path: Path) -> None:
    mapping = extract_order_mapping(path)
    print("excel_row\tserial\tgene_name\tuniprot")
    for row in mapping:
        print("\t".join(row[k] for k in ["excel_row", "serial", "gene_name", "uniprot"]))


def extract_dashboard_data(text: str) -> Tuple[Dict[str, Any], Tuple[int, int]]:
    marker = "const DATA = "
    start = text.find(marker)
    if start < 0:
        raise ValueError("Could not find dashboard DATA block")
    json_start = start + len(marker)
    end = text.find(";\nconst candidates", json_start)
    if end < 0:
        end = text.find(";\r\nconst candidates", json_start)
    if end < 0:
        raise ValueError("Could not find end of dashboard DATA block")
    return json.loads(text[json_start:end]), (json_start, end)


def dashboard_css_patch() -> str:
    return r"""
    .pdb small { display:block; color:var(--muted); font-weight:400; font-size:12px; margin-top:2px; }
    .identity-summary { display:grid; grid-template-columns:repeat(5,minmax(92px,1fr)); gap:8px; margin:8px 0 10px; }
    .identity-summary .stat { text-align:center; }
    .identity-heatmap-wrap { overflow:auto; max-height:460px; border:1px solid var(--line); border-radius:8px; background:#fff; margin-top:8px; }
    table.identity-heatmap { border-collapse:separate; border-spacing:0; font-size:11px; min-width:max-content; }
    .identity-heatmap th,.identity-heatmap td { border-right:1px solid rgba(255,255,255,.42); border-bottom:1px solid rgba(255,255,255,.42); min-width:28px; height:26px; text-align:center; }
    .identity-corner { position:sticky; top:0; left:0; z-index:4; min-width:120px; background:#fff; color:var(--muted); border-right:1px solid var(--line); border-bottom:1px solid var(--line); }
    .identity-col-head { position:sticky; top:0; z-index:3; height:120px; min-width:28px; vertical-align:bottom; background:#fff; border-bottom:1px solid var(--line); }
    .identity-col-head > div { writing-mode:vertical-rl; transform:rotate(180deg); white-space:nowrap; max-height:112px; overflow:hidden; padding:4px 2px; color:var(--muted); font-weight:400; }
    .identity-row-head { position:sticky; left:0; z-index:2; min-width:120px; max-width:120px; background:#fff; text-align:right; padding:0 6px; border-right:1px solid var(--line); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; color:var(--muted); font-weight:400; }
    .identity-cell { cursor:pointer; color:transparent; }
    .identity-cell.show-number { color:#202326; font-size:10px; }
    .identity-cell.high-number { color:#fff; }
    .identity-cell:hover { outline:2px solid #202326; outline-offset:-2px; }
    .identity-cell.diagonal { box-shadow:inset 0 0 0 1px rgba(0,0,0,.22); }
    .identity-cell.focus { box-shadow:inset 0 0 0 2px #202326; }
    .identity-pair-detail { display:grid; grid-template-columns:128px 1fr; gap:6px 10px; margin:10px 0 12px; font-size:13px; line-height:1.35; border:1px solid var(--line); border-radius:8px; padding:10px; background:#fff; }
    .identity-pair-detail div:nth-child(odd) { color:var(--muted); }
    .identity-legend { display:grid; grid-template-columns:36px 1fr 42px; gap:8px; align-items:start; font-size:12px; color:var(--muted); margin-top:8px; }
    .identity-ramp-wrap { position:relative; padding-bottom:22px; }
    .identity-ramp { height:12px; border-radius:999px; border:1px solid var(--line); background:linear-gradient(90deg,#C94C4C 0%,#ECA76A 25%,#F2E7A6 50%,#9BCB9C 75%,#63B6C2 90%,#3F77B5 100%); }
    .identity-tick { position:absolute; top:15px; transform:translateX(-50%); white-space:nowrap; font-size:11px; }
    .identity-tick::before { content:""; display:block; width:1px; height:6px; margin:0 auto 2px; background:var(--line); }
"""


def dashboard_js_patch() -> str:
    return r"""
function candidateByPdb(pdb){
  return candidates.find(c=>c.pdb_id===pdb) || null;
}
function displayLabel(c){
  if(!c) return "";
  return c.gene_name || c.gene_label || c.pdb_id || "";
}
function labelWithPdb(c){
  const label=displayLabel(c);
  if(!c || !c.pdb_id || label===c.pdb_id) return label;
  return `${label} (${c.pdb_id})`;
}
function labelForPdb(pdb){
  return labelWithPdb(candidateByPdb(pdb)) || pdb || "";
}
const IDENTITY_STOPS=[
  [0,[201,76,76]],[25,[236,167,106]],[50,[242,231,166]],
  [75,[155,203,156]],[90,[99,182,194]],[100,[63,119,181]]
];
function mixColor(a,b,t){return Math.round(a+(b-a)*t);}
function identityColor(v){
  const x=Math.max(0,Math.min(100,Number(v)||0));
  for(let i=1;i<IDENTITY_STOPS.length;i++){
    const [rv,rc]=IDENTITY_STOPS[i], [lv,lc]=IDENTITY_STOPS[i-1];
    if(x<=rv){
      const t=(x-lv)/Math.max(1,rv-lv);
      return `rgb(${mixColor(lc[0],rc[0],t)}, ${mixColor(lc[1],rc[1],t)}, ${mixColor(lc[2],rc[2],t)})`;
    }
  }
  return "rgb(63, 119, 181)";
}
function identityPctValue(p){
  if(!p) return null;
  const direct=num(p.identity_percent);
  if(direct!==null) return direct;
  const frac=num(p.identity);
  return frac===null?null:frac*100;
}
function identityPairKey(a,b){
  return [a,b].sort().join("__");
}
const identityPairMap=new Map();
identityPairs.forEach(p=>{
  if(p.pdb_id_a && p.pdb_id_b) identityPairMap.set(identityPairKey(p.pdb_id_a,p.pdb_id_b),p);
});
function identityFor(a,b){
  if(a===b) return 100;
  return identityPctValue(identityPairMap.get(identityPairKey(a,b)));
}
function selectedIdentityIds(anchorId){
  const order=candidates.map(c=>c.pdb_id);
  if(anchorId && !order.includes(anchorId)) order.unshift(anchorId);
  return order;
}
function identityStats(ids){
  const vals=[];
  for(let i=0;i<ids.length;i++){
    for(let j=i+1;j<ids.length;j++){
      const v=identityFor(ids[i],ids[j]);
      if(v!==null) vals.push(v);
    }
  }
  const high90=vals.filter(v=>v>=90).length;
  const high70=vals.filter(v=>v>=70).length;
  const avg=vals.length?vals.reduce((a,b)=>a+b,0)/vals.length:0;
  return {pairs:vals.length, high90, high70, avg};
}
function identityPairDetailHtml(a,b){
  const ca=candidateByPdb(a), cb=candidateByPdb(b), v=identityFor(a,b);
  const value=v===null?"n/a":`${v.toFixed(2)}%`;
  return `<div>Pair</div><div><b>${esc(labelWithPdb(ca))} vs ${esc(labelWithPdb(cb))}</b></div>
    <div>Identity</div><div><span class="chip">${esc(value)}</span></div>
    <div>A UniProt</div><div>${esc(ca?.uniprot_ids_primary || ca?.uniprot_ids || "")}</div>
    <div>B UniProt</div><div>${esc(cb?.uniprot_ids_primary || cb?.uniprot_ids || "")}</div>
    <div>A organism</div><div>${esc(ca?.source_organism_corrected || "")}</div>
    <div>B organism</div><div>${esc(cb?.source_organism_corrected || "")}</div>`;
}
function renderIdentityPairDetail(a,b){
  const el=document.getElementById("identityPairDetail");
  if(el) el.innerHTML=identityPairDetailHtml(a,b);
}
function identityHeatmapHtml(anchorId){
  const ids=selectedIdentityIds(anchorId);
  const stats=identityStats(ids);
  const head=`<tr><th class="identity-corner">${ids.length} seq</th>`+ids.map(id=>{
    const c=candidateByPdb(id);
    return `<th class="identity-col-head" title="${esc(labelWithPdb(c))}"><div>${esc(displayLabel(c))}</div></th>`;
  }).join("")+`</tr>`;
  const rows=ids.map(a=>{
    const ca=candidateByPdb(a);
    const cells=ids.map(b=>{
      const v=identityFor(a,b);
      const text=v===null?"":(v>=70?v.toFixed(0):"");
      const bg=v===null?"#f1f1ec":identityColor(v);
      const cls=["identity-cell",a===b?"diagonal":"",a===anchorId||b===anchorId?"focus":"",text?"show-number":"",v>=94?"high-number":""].filter(Boolean).join(" ");
      return `<td class="${cls}" style="background:${bg}" data-a="${esc(a)}" data-b="${esc(b)}" title="${esc(labelForPdb(a))} vs ${esc(labelForPdb(b))}: ${v===null?"n/a":v.toFixed(2)+"%"}">${esc(text)}</td>`;
    }).join("");
    return `<tr><th class="identity-row-head" title="${esc(labelWithPdb(ca))}">${esc(displayLabel(ca))}</th>${cells}</tr>`;
  }).join("");
  setTimeout(()=>{
    document.querySelectorAll(".identity-cell").forEach(td=>td.addEventListener("click",()=>{
      renderIdentityPairDetail(td.dataset.a,td.dataset.b);
    }));
    renderIdentityPairDetail(anchorId,anchorId);
  },0);
  return `<div class="identity-summary">
      <span class="stat">Sequences ${ids.length}</span>
      <span class="stat">Pairs ${stats.pairs}</span>
      <span class="stat">>=90% ${stats.high90}</span>
      <span class="stat">>=70% ${stats.high70}</span>
      <span class="stat">mean ${stats.avg.toFixed(1)}%</span>
    </div>
    <div id="identityPairDetail" class="identity-pair-detail"></div>
    <div class="identity-heatmap-wrap"><table class="identity-heatmap">${head}${rows}</table></div>
    <div class="identity-legend"><span>low</span><div class="identity-ramp-wrap"><div class="identity-ramp"></div>
      <span class="identity-tick" style="left:0%">0%</span><span class="identity-tick" style="left:50%">50%</span>
      <span class="identity-tick" style="left:75%">75%</span><span class="identity-tick" style="left:90%">90%</span>
      <span class="identity-tick" style="left:100%">100%</span></div><span>high</span></div>`;
}
"""


def rewrite_dashboard_js(text: str) -> str:
    replacements = [
        (
            '    const blob=[c.pdb_id,c.uniprot_ids_primary,c.source_organism_corrected,c.expression_system_corrected,c.fluorine_comp_ids,c.metal_comp_ids,c.risk_notes,c.doi,c.nearest_F_atom,c.nearest_metal_atom].join(" ").toLowerCase();',
            '    const blob=[c.gene_name,c.gene_label,c.pdb_id,c.uniprot_ids_primary,c.source_organism_corrected,c.expression_system_corrected,c.fluorine_comp_ids,c.metal_comp_ids,c.risk_notes,c.doi,c.nearest_F_atom,c.nearest_metal_atom].join(" ").toLowerCase();',
        ),
        (
            '<input id="q" placeholder="Search PDB, UniProt, organism, ligand, metal, note">',
            '<input id="q" placeholder="Search gene name, PDB, UniProt, organism, ligand, metal, note">',
        ),
        (
            '<div class="row"><div class="pdb">${esc(c.pdb_id)}</div><span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></div>\n      <div class="meta">${esc(c.source_organism_corrected||"unknown")} · ${esc(c.uniprot_ids_primary||"no UniProt")} · ${dist(c.nearest_F_metal_A)}</div>',
            '<div class="row"><div class="pdb">${esc(displayLabel(c))}<small>${esc(c.pdb_id)}</small></div><span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></div>\n      <div class="meta">${esc(c.source_organism_corrected||"unknown")} · ${esc(c.uniprot_ids_primary||"no UniProt")} · ${dist(c.nearest_F_metal_A)}</div>',
        ),
        (
            '    <div>PDB</div><div><a href="${esc(c.pdb_url)}" target="_blank">${esc(c.pdb_id)}</a></div>',
            '    <div>Gene name</div><div><b>${esc(displayLabel(c))}</b></div>\n    <div>PDB</div><div><a href="${esc(c.pdb_url)}" target="_blank">${esc(c.pdb_id)}</a></div>',
        ),
        (
            '<div>UniProt</div><div>${esc(c.uniprot_ids_primary || c.uniprot_ids || "")}</div>',
            '<div>UniProt</div><div>${esc(c.uniprot_ids_primary || c.uniprot_ids || "")}</div>\n    <div>Order table row</div><div>${esc(c.order_excel_row || "")}</div>',
        ),
        (
            '  return `<h3>Browser note for ${esc(c.pdb_id)}</h3><textarea id="userNote" placeholder="Add your plasmid-ordering note here. Saved only in this browser.">${esc(value)}</textarea><div class="hint">This note is stored in browser localStorage. It does not modify the CSV/XLSX on disk.</div>`;',
            '  return `<h3>Browser note for ${esc(labelWithPdb(c))}</h3><textarea id="userNote" placeholder="Add your plasmid-ordering note here. Saved only in this browser.">${esc(value)}</textarea><div class="hint">This note is stored in browser localStorage. It does not modify the CSV/XLSX on disk.</div>`;',
        ),
        (
            '  return rows.map(r=>`<h3>${esc(r.scene||c.pdb_id)}</h3><a href="${esc(r.png_relpath)}" target="_blank"><img class="scene-img" src="${esc(r.png_relpath)}" loading="lazy" alt="${esc(c.pdb_id)} scene"></a><div class="meta">${esc(r.fluorine_comp_ids||"")} / ${esc(r.metal_comp_ids||"")} · ${esc(r.nearest_F_metal_A||"")} A</div>`).join("");',
            '  return rows.map(r=>`<h3>${esc(displayLabel(c))} scene</h3><a href="${esc(r.png_relpath)}" target="_blank"><img class="scene-img" src="${esc(r.png_relpath)}" loading="lazy" alt="${esc(labelWithPdb(c))} scene"></a><div class="meta">${esc(r.scene||c.pdb_id)} · ${esc(r.fluorine_comp_ids||"")} / ${esc(r.metal_comp_ids||"")} · ${esc(r.nearest_F_metal_A||"")} A</div>`).join("");',
        ),
        (
            '    return `<div class="pair ${high}"><b>${esc(c.pdb_id)} vs ${esc(other)}</b><br>${esc(p.identity_percent||pct(p.identity))}</div>`;',
            '    return `<div class="pair ${high}"><b>${esc(labelWithPdb(c))} vs ${esc(labelForPdb(other))}</b><br>${esc(p.identity_percent||pct(p.identity))}</div>`;',
        ),
        (
            '  $("detail").innerHTML=`<h2>${esc(c.pdb_id)} <span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></h2>\n    <div class="meta">${esc(c.source_organism_corrected||"unknown")} · ${esc(c.expression_system_corrected||"unknown")} · ${dist(c.nearest_F_metal_A)}</div>',
            '  $("detail").innerHTML=`<h2>${esc(displayLabel(c))} <span class="chip ${esc(c.final_decision)}">${esc(c.final_decision)}</span></h2>\n    <div class="meta">${esc(c.pdb_id)} · ${esc(c.source_organism_corrected||"unknown")} · ${esc(c.expression_system_corrected||"unknown")} · ${dist(c.nearest_F_metal_A)}</div>',
        ),
    ]
    for old, new in replacements:
        if old not in text:
            raise ValueError(f"Could not find expected JS/HTML snippet for replacement: {old[:80]}")
        text = text.replace(old, new)

    anchor = """function entityBucket(c){
  const n=parseInt(c.protein_entity_count_corrected||0,10);
  if(n>1) return "complex";
  if(n===1) return "single protein";
  return "unknown";
}
"""
    if dashboard_js_patch().strip() not in text:
        text = text.replace(anchor, anchor + dashboard_js_patch())

    old_identity = """function identityTab(c){
  const rows=pairsFor(c.pdb_id);
  if(!rows.length) return "<p class='meta'>No pairwise identity rows.</p>";
  return rows.map(p=>{
    const other=p.pdb_id_a===c.pdb_id?p.pdb_id_b:p.pdb_id_a;
    const high=(p.identity||0)>=0.9?"high":"";
    return `<div class="pair ${high}"><b>${esc(labelWithPdb(c))} vs ${esc(labelForPdb(other))}</b><br>${esc(p.identity_percent||pct(p.identity))}</div>`;
  }).join("");
}
"""
    new_identity = """function identityTab(c){
  const rows=pairsFor(c.pdb_id);
  const pairRows=rows.length?rows.map(p=>{
    const other=p.pdb_id_a===c.pdb_id?p.pdb_id_b:p.pdb_id_a;
    const high=(p.identity||0)>=0.9?"high":"";
    return `<div class="pair ${high}"><b>${esc(labelWithPdb(c))} vs ${esc(labelForPdb(other))}</b><br>${esc(p.identity_percent||pct(p.identity))}</div>`;
  }).join(""):"<p class='meta'>No pairwise identity rows.</p>";
  return `<h3>Identity matrix</h3>${identityHeatmapHtml(c.pdb_id)}<h3>Pairs for ${esc(displayLabel(c))}</h3>${pairRows}`;
}
"""
    if old_identity not in text:
        raise ValueError("Could not find identityTab function for replacement")
    text = text.replace(old_identity, new_identity)
    return text


def update_dashboard(input_html: Path, order_xls: Path, output_html: Path) -> Dict[str, Any]:
    text = input_html.read_text(encoding="utf-8")
    data, (json_start, json_end) = extract_dashboard_data(text)
    order_rows = extract_order_mapping(order_xls)
    by_uniprot = {row["uniprot"]: row for row in order_rows}
    matched: List[Dict[str, str]] = []
    unmatched_order = set(by_uniprot)
    for candidate in data.get("candidates", []):
        uniprot = normalize_text(candidate.get("uniprot_ids_primary") or candidate.get("uniprot_ids") or "")
        uniprot = re.split(r"[\s,;]+", uniprot)[0] if uniprot else ""
        order = by_uniprot.get(uniprot)
        if order:
            candidate["gene_name"] = order["gene_name"]
            candidate["gene_label"] = order["gene_name"]
            candidate["order_excel_row"] = order["excel_row"]
            candidate["order_serial"] = order["serial"]
            matched.append({"pdb_id": candidate.get("pdb_id", ""), "uniprot": uniprot, "gene_name": order["gene_name"]})
            unmatched_order.discard(uniprot)
        else:
            candidate["gene_name"] = ""
            candidate["gene_label"] = candidate.get("pdb_id", "")
            candidate.pop("order_excel_row", None)
            candidate.pop("order_serial", None)
    data["order_gene_map"] = by_uniprot

    new_json = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    text = text[:json_start] + new_json + text[json_end:]
    css = dashboard_css_patch()
    if ".identity-heatmap-wrap" not in text:
        text = text.replace("    a { color:var(--blue); }\n", "    a { color:var(--blue); }\n" + css)
    text = rewrite_dashboard_js(text)
    text = re.sub(r'<section class="list" id="list">.*?</section>\s*<section class="detail" id="detail">.*?</section>', '<section class="list" id="list"></section>\n    <section class="detail" id="detail"></section>', text, count=1, flags=re.S)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(text, encoding="utf-8")
    return {
        "matched_count": len(matched),
        "order_count": len(order_rows),
        "candidate_count": len(data.get("candidates", [])),
        "matched": matched,
        "unmatched_order_uniprot": sorted(unmatched_order),
        "output": str(output_html),
    }


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump-xls", type=Path)
    parser.add_argument("--extract-order-mapping", type=Path)
    parser.add_argument("--update-dashboard", type=Path, help="Input dashboard HTML")
    parser.add_argument("--order-xls", type=Path)
    parser.add_argument("--output-html", type=Path)
    parser.add_argument("--max-rows", type=int, default=30)
    args = parser.parse_args()
    if args.dump_xls:
        dump_xls(args.dump_xls, args.max_rows)
    if args.extract_order_mapping:
        print_order_mapping(args.extract_order_mapping)
    if args.update_dashboard:
        if not args.order_xls or not args.output_html:
            raise SystemExit("--update-dashboard requires --order-xls and --output-html")
        summary = update_dashboard(args.update_dashboard, args.order_xls, args.output_html)
        print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

