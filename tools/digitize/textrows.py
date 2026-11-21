"""Print the text layer of PDF table pages as rows, with fractions merged and bold flagged.

Run: ``uv run --group digitize python -I -m tools.digitize.textrows PDF PAGES [--region ...]``
where PAGES is a comma-separated list of 1-based page numbers.

Kodak's developer publications set table times such as 6¾ in two ways. Some (J-109) stack the
fraction: the text layer holds "61⁄" on one baseline and the denominator "4" slightly lower. Others
(J-78, J-24, J-86) write "61⁄2" in one word or "6 ½" as two words. This tool joins all three forms
into "6¾" so that each table cell is one token, groups words into rows by baseline, and appends
"^B" to a word whose first glyph is set in a bold or semibold font (J-109 marks nominal speeds and
J-86 marks suggested starting points in bold). Each word is printed with its left edge in PDF points
(origin at the top of the page, as pdfplumber reports it).

The output is a reading aid for transcribing tables into ``staged/developers/*.yaml``; every cell
taken from it was also checked against the rendered page (``python -m tools.digitize render``).
"""

from __future__ import annotations

import argparse
import re

FRACTIONS = {"1/2": "½", "1/4": "¼", "3/4": "¾"}
SLASH = "⁄"  # FRACTION SLASH, as Kodak's PDFs encode it
VULGAR = "½¼¾"


def _is_bold(fontname: str) -> bool:
    return any(k in fontname for k in ("Bold", "Semibold"))


def _fraction(whole: str, num: str, den: str) -> str:
    return whole + FRACTIONS.get(f"{num}/{den}", f" {num}/{den}")


def words(page) -> list[dict]:
    """Words of a pdfplumber page with fractions merged and a ``bold`` flag."""
    raw = page.extract_words(x_tolerance=1.5, keep_blank_chars=False)
    chars = page.chars
    for w in raw:
        first = [
            c
            for c in chars
            if abs(c["x0"] - w["x0"]) < 0.05 and w["top"] - 0.5 <= c["top"] <= w["bottom"]
        ]
        w["fontname"] = first[0]["fontname"] if first else ""
    used: set[int] = set()
    out: list[dict] = []
    for i, w in enumerate(raw):
        text = w["text"]
        inline = re.fullmatch(r"(\(?)(\d*)(\d)" + SLASH + r"(\d)(.*)", text)
        if inline:
            g = inline.groups()
            out.append(dict(w, text=g[0] + _fraction(g[1], g[2], g[3]) + g[4]))
            continue
        stacked = re.fullmatch(r"(\(?)(\d*)(\d)" + SLASH + r"(\D*)", text)
        if stacked:
            best = None
            for j, d in enumerate(raw):
                if j in used or j == i or not re.fullmatch(r"\d", d["text"]):
                    continue
                dy = d["top"] - w["top"]
                dx = d["x0"] - w["x1"]
                if 0.5 < dy < 9 and -8 < dx < 6 and (best is None or abs(dx) + dy < best[0]):
                    best = (abs(dx) + dy, j)
            if best is not None:
                used.add(best[1])
                g = stacked.groups()
                out.append(dict(w, text=g[0] + _fraction(g[1], g[2], raw[best[1]]["text"]) + g[3]))
                continue
            out.append(dict(w, text=text + "?"))  # unmatched stacked fraction: flag it
            continue
        out.append(dict(w, _index=i))
    result = [w for w in out if w.get("_index") not in used]
    for w in result:
        w["bold"] = _is_bold(w["fontname"])
    return result


def rows(page, tolerance: float = 2.5) -> list[tuple[float, list[dict]]]:
    """Words grouped into rows by baseline; a lone ½ ¼ ¾ is joined to the number before it."""
    grouped: list[list] = []
    for w in sorted(words(page), key=lambda w: (w["top"], w["x0"])):
        if grouped and abs(grouped[-1][0] - w["top"]) <= tolerance:
            grouped[-1][1].append(w)
        else:
            grouped.append([w["top"], [w]])
    result = []
    for top, ws in grouped:
        ws.sort(key=lambda w: w["x0"])
        merged: list[dict] = []
        for w in ws:
            if (
                merged
                and w["text"][:1] in VULGAR
                and re.fullmatch(r"\d+", merged[-1]["text"])
                and w["x0"] - merged[-1]["x1"] < 6
            ):
                merged[-1] = dict(merged[-1], text=merged[-1]["text"] + w["text"], x1=w["x1"])
            else:
                merged.append(w)
        result.append((top, merged))
    return result


def main(argv: list[str] | None = None) -> int:
    import pdfplumber

    ap = argparse.ArgumentParser(prog="python -I -m tools.digitize.textrows")
    ap.add_argument("pdf")
    ap.add_argument("pages", help="comma-separated 1-based page numbers")
    ap.add_argument("--region", help="top,bottom: only rows whose top lies in this range")
    args = ap.parse_args(argv)
    lo, hi = (float(v) for v in args.region.split(",")) if args.region else (-1e9, 1e9)
    with pdfplumber.open(args.pdf) as pdf:
        for p in (int(v) for v in args.pages.split(",")):
            print(f"===== page {p} =====")
            for top, ws in rows(pdf.pages[p - 1]):
                if lo <= top <= hi:
                    cells = "  ".join(
                        f"{w['text']}{'^B' if w['bold'] else ''}@{w['x0']:.0f}" for w in ws
                    )
                    print(f"{top:6.1f} | {cells}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
