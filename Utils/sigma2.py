import argparse
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd


def load_sigma_values(
    input_csv: Union[str, Path],
    *,
    column: str = "sigma0_field",
    to_mpa: bool = True,
) -> np.ndarray:
    input_csv = Path(input_csv)
    df = pd.read_csv(input_csv)
    if column in df.columns:
        values = df[column]
    else:
        if df.shape[1] < 2:
            raise ValueError(f"CSV列数不足，无法自动选择sigma列: {input_csv}")
        values = df.iloc[:, 1]

    values = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if values.size == 0:
        raise ValueError(f"未读取到有效sigma数据: {input_csv}")

    if to_mpa:
        values = values / 1e6

    return values


def export_sigma_distribution(
    input_csv: Union[str, Path],
    output_csv: Union[str, Path],
    *,
    column: str = "sigma0_field",
    to_mpa: bool = True,
    method: str = "kde",
    points: int = 400,
    bins: int = 200,
    x_min: Optional[float] = None,
    x_max: Optional[float] = None,
) -> pd.DataFrame:
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)

    values = load_sigma_values(input_csv, column=column, to_mpa=to_mpa)

    if x_min is None:
        x_min = float(np.min(values))
    if x_max is None:
        x_max = float(np.max(values))
    if not (x_max > x_min):
        raise ValueError("x_max 必须大于 x_min")

    method = method.lower().strip()
    if method == "kde":
        try:
            from scipy.stats import gaussian_kde
        except Exception as e:
            raise RuntimeError("method=kde 需要安装 scipy") from e

        x = np.linspace(x_min, x_max, int(points))
        kde = gaussian_kde(values)
        pdf = kde(x)
        dx = float(x[1] - x[0]) if x.size > 1 else 0.0
        if dx > 0:
            cdf = np.concatenate(([0.0], np.cumsum((pdf[:-1] + pdf[1:]) * 0.5 * dx)))
        else:
            cdf = np.zeros_like(x)
    elif method == "hist":
        pdf, edges = np.histogram(values, bins=int(bins), range=(x_min, x_max), density=True)
        x = (edges[:-1] + edges[1:]) * 0.5
        dx = float(edges[1] - edges[0]) if edges.size > 1 else 0.0
        cdf = np.cumsum(pdf * dx) if dx > 0 else np.zeros_like(x)
    else:
        raise ValueError("method 仅支持 kde 或 hist")

    out = pd.DataFrame({"x_mpa": x, "pdf": pdf, "cdf": cdf})
    out.to_csv(output_csv, index=False, encoding="utf-8-sig")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="sigma0_distribution_sigma0_dist_0000.csv",
        help="输入CSV文件路径",
    )
    parser.add_argument(
        "--output",
        default="sigma_distribution_export.csv",
        help="输出CSV文件路径",
    )
    parser.add_argument(
        "--column",
        default="sigma0_field",
        help="sigma列名（默认 sigma0_field）。若找不到将自动用第2列",
    )
    parser.add_argument("--method", default="kde", choices=["kde", "hist"])
    parser.add_argument("--points", type=int, default=400, help="kde输出点数（小区间数）")
    parser.add_argument("--bins", type=int, default=200, help="hist分箱数（小区间数）")
    parser.add_argument("--x-min", type=float, default=None)
    parser.add_argument("--x-max", type=float, default=None)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    input_path = Path(args.input)
    if not input_path.is_absolute():
        candidate = script_dir / input_path
        if candidate.exists():
            input_path = candidate

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = script_dir / output_path

    values = load_sigma_values(input_path, column=args.column, to_mpa=True)
    out = export_sigma_distribution(
        input_path,
        output_path,
        column=args.column,
        method=args.method,
        points=args.points,
        bins=args.bins,
        x_min=args.x_min,
        x_max=args.x_max,
    )
    mean_value = float(np.mean(values)) if values.size else float("nan")
    std_value = float(np.std(values)) if values.size else float("nan")
    print(f"Wrote: {Path(output_path).resolve()}")
    print(f"Mean(x_mpa): {mean_value:.6f}")
    print(f"Std(x_mpa): {std_value:.6f}")


if __name__ == "__main__":
    main()
