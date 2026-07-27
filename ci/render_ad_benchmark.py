#!/usr/bin/env python3
"""Render dependency-free overlaid AD versus finite-difference timing bars."""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("svg_path", type=Path)
    parser.add_argument("--title", default="QMBED native AD vs central finite differences")
    args = parser.parse_args()

    with args.csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    times = [
        float(row[key])
        for row in rows
        for key in ("analytic_seconds", "finite_difference_seconds")
    ]
    minimum = 10 ** math.floor(math.log10(min(times)))
    maximum = 10 ** math.ceil(math.log10(max(times)))
    width = 1600
    left = 415
    plot_width = 880
    row_height = 54
    top = 130
    height = top + row_height * len(rows) + 95

    def x(value: float) -> float:
        fraction = (math.log10(value) - math.log10(minimum)) / (
            math.log10(maximum) - math.log10(minimum)
        )
        return left + plot_width * fraction

    output = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Inter,Arial,sans-serif;fill:#172033}",
        ".title{font-size:30px;font-weight:700}.subtitle{font-size:17px;fill:#536179}",
        ".label{font-size:16px;font-weight:600}.value{font-size:14px}.tick{font-size:14px;fill:#69758a}",
        "</style>",
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="42" y="48" class="title">{html.escape(args.title)}</text>',
        '<rect x="42" y="76" width="30" height="18" rx="4" fill="#c6ccd8"/>',
        '<text x="82" y="91" class="subtitle">central finite difference</text>',
        '<rect x="310" y="79" width="30" height="12" rx="4" fill="#2768d7"/>',
        '<text x="350" y="91" class="subtitle">native AD</text>',
        '<text x="1450" y="91" class="subtitle">FD / AD</text>',
    ]
    decade = minimum
    while decade <= maximum * (1.0 + 1.0e-12):
        coordinate = x(decade)
        output.extend(
            [
                f'<line x1="{coordinate:.1f}" y1="{top - 18}" x2="{coordinate:.1f}" y2="{height - 65}" stroke="#e6e9ef" stroke-width="1"/>',
                f'<text x="{coordinate:.1f}" y="{height - 38}" text-anchor="middle" class="tick">{decade:g} s</text>',
            ]
        )
        decade *= 10

    for index, row in enumerate(rows):
        center = top + index * row_height
        analytic = float(row["analytic_seconds"])
        finite = float(row["finite_difference_seconds"])
        fd_end = x(finite)
        ad_end = x(analytic)
        duration_end = max(fd_end, ad_end)
        duration_x = duration_end - 10 if duration_end > 1120 else duration_end + 10
        duration_anchor = "end" if duration_end > 1120 else "start"
        output.extend(
            [
                f'<text x="42" y="{center + 6}" class="label">{html.escape(row["case_id"].removeprefix("ad_").replace("_ground_energy", "").replace("_", " "))}</text>',
                f'<rect x="{left}" y="{center - 15}" width="{max(fd_end - left, 2):.1f}" height="30" rx="6" fill="#c6ccd8"/>',
                f'<rect x="{left}" y="{center - 8}" width="{max(ad_end - left, 2):.1f}" height="16" rx="5" fill="#2768d7"/>',
                f'<text x="{duration_x:.1f}" y="{center + 5}" text-anchor="{duration_anchor}" class="value">{analytic:.4g} / {finite:.4g} s</text>',
                f'<text x="1450" y="{center + 6}" class="label">{float(row["speedup"]):.2f}×</text>',
            ]
        )
    output.extend(
        [
            f'<text x="42" y="{height - 17}" class="subtitle">Bars share a logarithmic time axis. Labels show native AD / finite-difference median wall time; lower is better.</text>',
            "</svg>",
        ]
    )
    args.svg_path.write_text("\n".join(output) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
