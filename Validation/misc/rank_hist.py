#!/usr/bin/env python3
"""Create side-by-side Plotly box plots from a pipeline count table."""

from __future__ import annotations

import argparse
import csv
import math
import random
from pathlib import Path
from typing import Iterable

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError as error:
    raise SystemExit(
        "rank_hist.py requires Plotly. Install it with: python -m pip install plotly"
    ) from error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one box plot per numeric column in a tab-delimited count table. "
            "A header line, when present, must begin with '#'."
        )
    )
    parser.add_argument("input", type=Path, help="Input count table.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("rank_hist.html"),
        help="Output HTML file (default: rank_hist.html).",
    )
    parser.add_argument(
        "--mode",
        choices=("absolute", "relative"),
        default="absolute",
        help="Plot raw counts or normalized values (default: absolute).",
    )
    parser.add_argument(
        "--relative-to",
        type=float,
        metavar="NUMBER",
        help=(
            "Divide every numeric column by NUMBER in relative mode. Without this "
            "option, each column is divided row-wise by the preceding numeric column; "
            "the first numeric column is divided by 1."
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open the generated figure in a browser after writing it.",
    )
    args = parser.parse_args()
    if args.relative_to is not None and args.relative_to == 0:
        parser.error("--relative-to cannot be zero")
    return args


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a tab-delimited table whose optional header starts with '#'."""
    header: list[str] | None = None
    rows: list[list[str]] = []

    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = next(csv.reader([line[1:]], delimiter="\t"))
                continue
            rows.append(next(csv.reader([line], delimiter="\t")))

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    n_columns = max(len(row) for row in rows)
    rows = [row + [""] * (n_columns - len(row)) for row in rows]
    if header is None:
        header = [""] * n_columns
    else:
        header += [""] * (n_columns - len(header))
    return header, rows


def as_number(value: str) -> float | None:
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def numeric_columns(rows: list[list[str]]) -> list[int]:
    """Return columns containing at least one numeric value.

    This deliberately excludes identifier columns such as '#Sample'.
    """
    return [
        column
        for column in range(len(rows[0]))
        if any(as_number(row[column]) is not None for row in rows)
    ]


def relative_values(
    rows: list[list[str]], columns: list[int], relative_to: float | None
) -> tuple[list[list[float]], list[list[str]]]:
    """Normalize values using a global baseline or the prior numeric column."""
    values = [[] for _ in columns]
    point_labels = [[] for _ in columns]
    for row in rows:
        previous: float | None = None
        for index, column in enumerate(columns):
            value = as_number(row[column])
            if value is None:
                previous = None
                continue
            denominator = relative_to if relative_to is not None else (previous if index else 1.0)
            if denominator not in (None, 0):
                ratio = value / denominator
                if math.isfinite(ratio):
                    values[index].append(ratio)
                    point_labels[index].append(row[0])
            previous = value
    return values, point_labels


def absolute_values(
    rows: list[list[str]], columns: list[int]
) -> tuple[list[list[float]], list[list[str]]]:
    values = [[] for _ in columns]
    point_labels = [[] for _ in columns]
    for row in rows:
        for index, column in enumerate(columns):
            value = as_number(row[column])
            if value is not None:
                values[index].append(value)
                point_labels[index].append(row[0])
    return values, point_labels


def column_label(header: list[str], column: int) -> str:
    return header[column] or " "


def add_boxplots(
    figure: go.Figure,
    values: Iterable[list[float]],
    point_labels: Iterable[list[str]],
    labels: Iterable[str],
) -> None:
    for index, (column_values, column_point_labels, label) in enumerate(
        zip(values, point_labels, labels), start=1
    ):
        show_pointcloud = len(column_values) < 100
        figure.add_trace(
            go.Box(
                x=[0.0] * len(column_values),
                y=column_values,
                boxpoints=False if show_pointcloud else "outliers",
                marker_color="#197278",
                line_color="#0b4f52",
                customdata=column_point_labels,
                hovertemplate="sample: %{customdata}<br>value: %{y:.4g}<extra></extra>",
                showlegend=False,
            ),
            row=1,
            col=index,
        )
        if show_pointcloud:
            jitter = random.Random(index)
            figure.add_trace(
                go.Scatter(
                    x=[jitter.uniform(-0.24, 0.24) for _ in column_values],
                    y=column_values,
                    mode="markers",
                    marker={"color": "#197278", "opacity": 0.72, "size": 7},
                    customdata=column_point_labels,
                    hovertemplate="sample: %{customdata}<br>value: %{y:.4g}<extra></extra>",
                    showlegend=False,
                ),
                row=1,
                col=index,
            )
        figure.update_xaxes(
            range=(-0.5, 0.5), showgrid=False, showticklabels=False, row=1, col=index
        )


def main() -> None:
    args = parse_args()
    header, rows = read_table(args.input)
    columns = numeric_columns(rows)
    if not columns:
        raise ValueError("The input table contains no numeric columns")

    if args.mode == "relative":
        values, point_labels = relative_values(rows, columns, args.relative_to)
        x_title = "Relative proportion"
    else:
        values, point_labels = absolute_values(rows, columns)
        x_title = "Absolute count"

    labels = [column_label(header, column) for column in columns]
    figure = make_subplots(
        rows=1,
        cols=len(columns),
        subplot_titles=labels,
        shared_yaxes=args.mode == "relative",
        horizontal_spacing=min(0.03, 0.25 / len(columns)),
    )
    add_boxplots(figure, values, point_labels, labels)
    figure.update_layout(
        title=f"{args.input.name}: {args.mode} count distributions",
        template="plotly_white",
        height=520,
        width=max(1100, 190 * len(columns)),
        margin={"b": 90},
    )
    figure.add_annotation(
        text=x_title,
        xref="paper",
        yref="paper",
        x=0.5,
        y=-0.18,
        showarrow=False,
        font={"size": 14},
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.output, include_plotlyjs=True, auto_open=args.show)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
