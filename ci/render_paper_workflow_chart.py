#!/usr/bin/env python3
"""Render the paired paper-workflow medians as a standalone SVG chart."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path


CASE_ORDER = (
    "paper_mbl_shift_invert_l14",
    "paper_xxz_lanczos_quench_l16",
    "paper_floquet_heating_l9",
    "paper_spinful_hubbard_l8",
    "paper_interacting_ssh_l16",
    "paper_translation_xxz_l18",
    "paper_tfim_fidelity_l16",
    "paper_pxp_revival_l24",
    "paper_bose_hubbard_quench_l11",
    "paper_hubbard_current_l10",
    "paper_conb_dsf_l16",
    "paper_particle_addition_6x3",
)


@dataclass(frozen=True)
class TimingPair:
    case_id: str
    label: str
    python_ms: float
    julia_ms: float


def _read_language(path: Path, expected_language: str) -> dict[str, dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        if row.get("suite") != "paper":
            continue
        language = row.get("language")
        if language != expected_language:
            raise ValueError(
                f"{path} contains language {language!r}; expected {expected_language!r}"
            )
        case_id = row.get("case_id", "")
        if not case_id:
            raise ValueError(f"{path} contains a paper row without case_id")
        if case_id in selected:
            raise ValueError(f"{path} contains duplicate case_id {case_id!r}")
        if row.get("supported", "true").lower() != "true":
            raise ValueError(f"{path} marks {case_id!r} as unsupported")
        selected[case_id] = row
    return selected


def load_pairs(python_csv: Path, julia_csv: Path) -> list[TimingPair]:
    python_rows = _read_language(python_csv, "python")
    julia_rows = _read_language(julia_csv, "julia")
    expected = set(CASE_ORDER)
    for language, rows in (("Python", python_rows), ("Julia", julia_rows)):
        missing = expected - rows.keys()
        extra = rows.keys() - expected
        if missing or extra:
            details = []
            if missing:
                details.append(f"missing {', '.join(sorted(missing))}")
            if extra:
                details.append(f"unexpected {', '.join(sorted(extra))}")
            raise ValueError(f"{language} paper CSV: {'; '.join(details)}")

    pairs = []
    for case_id in CASE_ORDER:
        python = python_rows[case_id]
        julia = julia_rows[case_id]
        python_label = python.get("benchmark", case_id)
        julia_label = julia.get("benchmark", case_id)
        if python_label != julia_label:
            raise ValueError(
                f"benchmark label mismatch for {case_id}: "
                f"{python_label!r} != {julia_label!r}"
            )
        python_ms = float(python["median_seconds"]) * 1_000
        julia_ms = float(julia["median_seconds"]) * 1_000
        if not (math.isfinite(python_ms) and python_ms > 0):
            raise ValueError(f"invalid Python median for {case_id}: {python_ms}")
        if not (math.isfinite(julia_ms) and julia_ms > 0):
            raise ValueError(f"invalid Julia median for {case_id}: {julia_ms}")
        pairs.append(TimingPair(case_id, python_label, python_ms, julia_ms))
    return pairs


def _nice_axis_max(maximum: float, target_ticks: int = 5) -> tuple[float, float]:
    rough_step = maximum / target_ticks
    magnitude = 10 ** math.floor(math.log10(rough_step))
    normalized = rough_step / magnitude
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if normalized <= candidate:
            step = candidate * magnitude
            break
    axis_max = math.ceil(maximum / step) * step
    return axis_max, step


def _text(value: object) -> str:
    return html.escape(str(value), quote=True)


def _speedup_annotation(python_ms: float, julia_ms: float) -> tuple[str, str]:
    ratio = python_ms / julia_ms
    if ratio > 1.01:
        return f"Julia {ratio:.2f}× faster", "faster"
    if ratio < 0.99:
        return f"Julia {1 / ratio:.2f}× slower", "slower"
    return f"Julia ≈ parity ({ratio:.2f}×)", "parity"


def render_svg(
    pairs: list[TimingPair],
    *,
    run_id: str,
    candidate_repository: str,
    candidate_ref: str,
) -> str:
    if len(pairs) != len(CASE_ORDER):
        raise ValueError(
            f"expected {len(CASE_ORDER)} timing pairs, found {len(pairs)}"
        )

    width = 1180
    left = 335
    right = 225
    top = 150
    bottom = 70
    row_height = 58
    python_bar_height = 22
    julia_bar_height = 10
    plot_width = width - left - right
    annotation_x = left + plot_width + 18
    height = top + row_height * len(pairs) + bottom
    axis_max, tick_step = _nice_axis_max(
        max(max(pair.python_ms, pair.julia_ms) for pair in pairs)
    )

    short_ref = candidate_ref if len(candidate_ref) <= 16 else candidate_ref[:12]
    subtitle = (
        f"GitHub Actions run {run_id} · {candidate_repository}@{short_ref} · "
        "median wall time, lower is better"
    )
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="chart-title chart-desc">'
        ),
        '<title id="chart-title">Twelve paper-workflow end-to-end timings</title>',
        (
            '<desc id="chart-desc">Overlaid horizontal bars compare Python QuSpin '
            'and Julia QuSpin median wall times for twelve paper-shaped workflows. '
            'Each row directly labels whether Julia is faster, slower, or at parity.</desc>'
        ),
        '<rect width="100%" height="100%" rx="12" fill="#ffffff"/>',
        (
            '<text x="28" y="42" font-family="system-ui, sans-serif" '
            'font-size="24" font-weight="700" fill="#172033">'
            'Twelve paper-workflow end-to-end timings</text>'
        ),
        (
            '<text x="28" y="70" font-family="system-ui, sans-serif" '
            f'font-size="13" fill="#596579">{_text(subtitle)}</text>'
        ),
        '<rect x="28" y="92" width="18" height="12" rx="2" fill="#377eb8"/>',
        (
            '<text x="54" y="103" font-family="system-ui, sans-serif" '
            'font-size="13" fill="#344054">Python QuSpin</text>'
        ),
        '<rect x="174" y="92" width="18" height="12" rx="2" fill="#e76f51"/>',
        (
            '<text x="200" y="103" font-family="system-ui, sans-serif" '
            'font-size="13" fill="#344054">Julia QuSpin</text>'
        ),
    ]

    tick = 0.0
    while tick <= axis_max + tick_step / 10:
        x = left + (tick / axis_max) * plot_width
        lines.extend(
            [
                (
                    f'<line x1="{x:.2f}" y1="{top - 12}" x2="{x:.2f}" '
                    f'y2="{height - bottom + 4}" stroke="#e6e9ef" stroke-width="1"/>'
                ),
                (
                    f'<text x="{x:.2f}" y="{height - bottom + 28}" '
                    'text-anchor="middle" font-family="system-ui, sans-serif" '
                    f'font-size="12" fill="#667085">{tick:g}</text>'
                ),
            ]
        )
        tick += tick_step
    lines.append(
        (
            f'<text x="{left + plot_width / 2:.2f}" y="{height - 16}" '
            'text-anchor="middle" font-family="system-ui, sans-serif" '
            'font-size="13" fill="#475467">median wall time (ms)</text>'
        )
    )

    for index, pair in enumerate(pairs):
        row_top = top + index * row_height
        center = row_top + row_height / 2
        python_y = center - python_bar_height / 2
        julia_y = center - julia_bar_height / 2
        lines.append(
            (
                f'<text x="{left - 16}" y="{center + 5:.2f}" text-anchor="end" '
                'font-family="system-ui, sans-serif" font-size="13" '
                f'font-weight="600" fill="#283548">{_text(pair.label)}</text>'
            )
        )
        for language, value, y, bar_height, color, css_class in (
            (
                "Python",
                pair.python_ms,
                python_y,
                python_bar_height,
                "#377eb8",
                "python-bar",
            ),
            (
                "Julia",
                pair.julia_ms,
                julia_y,
                julia_bar_height,
                "#e76f51",
                "julia-bar",
            ),
        ):
            bar_width = (value / axis_max) * plot_width
            lines.append(
                (
                    f'<rect class="{css_class}" data-case-id="{_text(pair.case_id)}" '
                    f'data-center-y="{center:.2f}" '
                    f'x="{left}" y="{y:.2f}" width="{bar_width:.2f}" '
                    f'height="{bar_height}" rx="3" fill="{color}">'
                    f'<title>{language}: {value:.3f} ms</title></rect>'
                )
            )
        annotation, speed_class = _speedup_annotation(
            pair.python_ms,
            pair.julia_ms,
        )
        lines.extend(
            [
                (
                    f'<text class="speedup-label {speed_class}" '
                    f'data-case-id="{_text(pair.case_id)}" '
                    f'x="{annotation_x:.2f}" y="{center - 2:.2f}" '
                    'font-family="system-ui, sans-serif" font-size="13" '
                    f'font-weight="600" fill="#344054">{_text(annotation)}</text>'
                ),
                (
                    f'<text class="timing-label" '
                    f'data-case-id="{_text(pair.case_id)}" '
                    f'x="{annotation_x:.2f}" y="{center + 14:.2f}" '
                    'font-family="ui-monospace, SFMono-Regular, monospace" '
                    'font-size="11" fill="#667085">'
                    f'P {pair.python_ms:.1f} · J {pair.julia_ms:.1f} ms</text>'
                ),
            ]
        )

    lines.append('</svg>')
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("python_csv", type=Path)
    parser.add_argument("julia_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", default="local")
    parser.add_argument(
        "--candidate-repository", default="matrixlab-research/QuSpin.jl"
    )
    parser.add_argument("--candidate-ref", default="unknown")
    args = parser.parse_args()
    chart = render_svg(
        load_pairs(args.python_csv, args.julia_csv),
        run_id=args.run_id,
        candidate_repository=args.candidate_repository,
        candidate_ref=args.candidate_ref,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(chart)


if __name__ == "__main__":
    main()
