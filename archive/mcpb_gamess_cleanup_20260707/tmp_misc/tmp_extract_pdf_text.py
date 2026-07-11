#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import zlib
from pathlib import Path


def unescape_pdf_string(raw: bytes) -> str:
    out = bytearray()
    i = 0
    while i < len(raw):
        c = raw[i]
        if c != 0x5C:
            out.append(c)
            i += 1
            continue
        i += 1
        if i >= len(raw):
            break
        e = raw[i]
        i += 1
        if e in b"nrtbf":
            out.extend({ord("n"): b"\n", ord("r"): b"\r", ord("t"): b"\t", ord("b"): b"\b", ord("f"): b"\f"}[e])
        elif e in b"()\\":
            out.append(e)
        elif 48 <= e <= 55:
            octal = bytes([e])
            for _ in range(2):
                if i < len(raw) and 48 <= raw[i] <= 55:
                    octal += bytes([raw[i]])
                    i += 1
            out.append(int(octal, 8))
        elif e in b"\r\n":
            while i < len(raw) and raw[i] in b"\r\n":
                i += 1
        else:
            out.append(e)
    return out.decode("latin1", errors="ignore")


def extract_pdf_strings(stream: bytes) -> list[str]:
    strings: list[str] = []
    i = 0
    while i < len(stream):
        if stream[i] != ord("("):
            i += 1
            continue
        depth = 1
        j = i + 1
        escaped = False
        raw = bytearray()
        while j < len(stream) and depth:
            c = stream[j]
            if escaped:
                raw.append(ord("\\"))
                raw.append(c)
                escaped = False
            elif c == ord("\\"):
                escaped = True
            elif c == ord("("):
                depth += 1
                raw.append(c)
            elif c == ord(")"):
                depth -= 1
                if depth:
                    raw.append(c)
            else:
                raw.append(c)
            j += 1
        if raw:
            text = unescape_pdf_string(bytes(raw))
            if text and not text.isspace():
                strings.append(text)
        i = max(j, i + 1)
    return strings


def decompress_streams(data: bytes) -> list[bytes]:
    streams: list[bytes] = []
    for m in re.finditer(rb"stream\r?\n", data):
        start = m.end()
        end = data.find(b"endstream", start)
        if end < 0:
            continue
        raw = data[start:end].strip(b"\r\n")
        for wbits in (zlib.MAX_WBITS, -zlib.MAX_WBITS):
            try:
                streams.append(zlib.decompress(raw, wbits))
                break
            except Exception:
                pass
    return streams


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf")
    parser.add_argument("--pattern", default="")
    args = parser.parse_args()
    data = Path(args.pdf).read_bytes()
    chunks = []
    for idx, stream in enumerate(decompress_streams(data)):
        strings = extract_pdf_strings(stream)
        if not strings:
            continue
        text = " ".join(s.strip() for s in strings if s.strip())
        if args.pattern and args.pattern.lower() not in text.lower():
            continue
        chunks.append((idx, text))
    for idx, text in chunks:
        print(f"\n--- STREAM {idx} ---")
        print(text[:12000])
    print(f"\nMATCHED_STREAMS {len(chunks)}")


if __name__ == "__main__":
    main()
