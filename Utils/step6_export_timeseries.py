import csv
import glob
import os
import re


def _parse_fields(value: str):
    result = []
    if not value:
        return result
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            name, _label = part.split(":", 1)
            result.append(name.strip())
        else:
            result.append(part)
    return result

def _find_first_input(case_dir):
    for name in sorted(os.listdir(case_dir)):
        if name.endswith(".i"):
            return os.path.join(case_dir, name)
    return None

def _parse_file_base(input_path):
    if not input_path or not os.path.isfile(input_path):
        return None
    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "#" in line:
                line = line.split("#", 1)[0]
            if "file_base" not in line:
                continue
            m = re.search(r"file_base\s*=\s*['\"]([^'\"]+)['\"]", line)
            if m:
                return m.group(1).strip()
    return None

def _strip_placeholders(value):
    return re.sub(r"\$\{[^}]+\}", "", value)


def _reduce_array(arr, reduction: str):
    reduction = (reduction or "mean").lower()
    n = arr.GetNumberOfTuples()
    if n <= 0:
        return None
    if arr.GetNumberOfComponents() <= 0:
        return None
    if reduction == "mean":
        s = 0.0
        for i in range(n):
            s += float(arr.GetComponent(i, 0))
        return s / float(n)
    if reduction == "min":
        v = float(arr.GetComponent(0, 0))
        for i in range(1, n):
            v = min(v, float(arr.GetComponent(i, 0)))
        return v
    if reduction == "max":
        v = float(arr.GetComponent(0, 0))
        for i in range(1, n):
            v = max(v, float(arr.GetComponent(i, 0)))
        return v
    raise ValueError(f"unknown reduction: {reduction}")


def _iter_leaf_datasets(obj):
    if obj is None:
        return

    if hasattr(obj, "GetNumberOfPartitionedDataSets") and hasattr(obj, "GetPartitionedDataSet"):
        for i in range(obj.GetNumberOfPartitionedDataSets()):
            yield from _iter_leaf_datasets(obj.GetPartitionedDataSet(i))
        return

    if hasattr(obj, "GetNumberOfPartitions") and hasattr(obj, "GetPartition"):
        for i in range(obj.GetNumberOfPartitions()):
            yield from _iter_leaf_datasets(obj.GetPartition(i))
        return

    if hasattr(obj, "GetNumberOfBlocks") and hasattr(obj, "GetBlock"):
        for i in range(obj.GetNumberOfBlocks()):
            yield from _iter_leaf_datasets(obj.GetBlock(i))
        return

    if hasattr(obj, "GetCellData") and hasattr(obj, "GetPointData") and hasattr(obj, "GetFieldData"):
        yield obj


def _fielddata_scalar(obj, name: str):
    if obj is None or not hasattr(obj, "GetFieldData"):
        return None
    fd = obj.GetFieldData()
    if fd is None:
        return None
    arr = fd.GetArray(name)
    if arr is None:
        return None
    if arr.GetNumberOfTuples() < 1 or arr.GetNumberOfComponents() < 1:
        return None
    return float(arr.GetComponent(0, 0))


def _reduce_over_leaves(data, name: str, reduction: str, where: str):
    reduction = (reduction or "mean").lower()
    if reduction == "mean":
        s = 0.0
        n_total = 0
        for ds in _iter_leaf_datasets(data):
            arr = None
            if where == "cell":
                cd = ds.GetCellData()
                arr = cd.GetArray(name) if cd is not None else None
            elif where == "point":
                pd = ds.GetPointData()
                arr = pd.GetArray(name) if pd is not None else None
            if arr is None:
                continue
            n = arr.GetNumberOfTuples()
            if n <= 0 or arr.GetNumberOfComponents() <= 0:
                continue
            for i in range(n):
                s += float(arr.GetComponent(i, 0))
            n_total += n
        if n_total == 0:
            return None
        return s / float(n_total)

    if reduction == "min":
        v = None
        for ds in _iter_leaf_datasets(data):
            arr = None
            if where == "cell":
                cd = ds.GetCellData()
                arr = cd.GetArray(name) if cd is not None else None
            elif where == "point":
                pd = ds.GetPointData()
                arr = pd.GetArray(name) if pd is not None else None
            if arr is None:
                continue
            r = _reduce_array(arr, "min")
            if r is None:
                continue
            v = r if v is None else min(v, r)
        return v

    if reduction == "max":
        v = None
        for ds in _iter_leaf_datasets(data):
            arr = None
            if where == "cell":
                cd = ds.GetCellData()
                arr = cd.GetArray(name) if cd is not None else None
            elif where == "point":
                pd = ds.GetPointData()
                arr = pd.GetArray(name) if pd is not None else None
            if arr is None:
                continue
            r = _reduce_array(arr, "max")
            if r is None:
                continue
            v = r if v is None else max(v, r)
        return v

    raise ValueError(f"unknown reduction: {reduction}")


def _get_scalar_at_t(data, name: str, reduction: str):
    if data is None:
        return None, None

    v = _fielddata_scalar(data, name)
    if v is not None:
        return v, "field"
    for ds in _iter_leaf_datasets(data):
        v = _fielddata_scalar(ds, name)
        if v is not None:
            return v, "field"

    v = _reduce_over_leaves(data, name, reduction, "cell")
    if v is not None:
        return v, "cell"

    v = _reduce_over_leaves(data, name, reduction, "point")
    if v is not None:
        return v, "point"

    return None, None


def export_exodus_timeseries(
    exodus_path: str,
    fields,
    output_csv_path: str,
    reduction: str = "mean",
):
    from paraview.simple import OpenDataFile, UpdatePipeline
    from paraview import servermanager

    reader = OpenDataFile(exodus_path)
    if reader is None:
        raise RuntimeError(f"cannot open: {exodus_path}")

    time_steps = []
    if hasattr(reader, "TimestepValues") and reader.TimestepValues:
        time_steps = [float(t) for t in reader.TimestepValues]
    elif hasattr(reader, "TimestepValue"):
        time_steps = [float(reader.TimestepValue)]

    if not time_steps:
        raise RuntimeError(f"no timesteps: {exodus_path}")

    os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)) or ".", exist_ok=True)

    series = {name: [] for name in fields}
    with open(output_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time", *fields])
        for t in time_steps:
            UpdatePipeline(time=t, proxy=reader)
            data = servermanager.Fetch(reader)
            row = [t]
            for name in fields:
                v, _src = _get_scalar_at_t(data, name, reduction)
                if v is None:
                    row.append("")
                    series[name].append(None)
                else:
                    row.append(v)
                    series[name].append(v)
            w.writerow(row)
    return time_steps, series


def _pick_primary_exodus(case_dir: str, output_dir_name: str):
    input_path = _find_first_input(case_dir)
    file_base = _parse_file_base(input_path)
    if file_base:
        base_dir = os.path.dirname(file_base)
        base_name = os.path.basename(file_base)
        prefix = _strip_placeholders(base_name)
        output_dir = os.path.join(case_dir, base_dir)
        if os.path.isdir(output_dir):
            matches = []
            for root, _, files in os.walk(output_dir):
                for file in files:
                    if not file.endswith(".e"):
                        continue
                    if prefix and not file.startswith(prefix):
                        continue
                    matches.append(os.path.join(root, file))
            if matches:
                return max(matches, key=os.path.getmtime)

    out_candidates = sorted(glob.glob(os.path.join(case_dir, "*_out.e")))
    if out_candidates:
        return out_candidates[0]
    candidates = []
    for root, _, files in os.walk(case_dir):
        if output_dir_name and os.path.sep + output_dir_name in root:
            continue
        for file in files:
            if file.endswith(".e"):
                candidates.append(os.path.join(root, file))
    if candidates:
        return max(candidates, key=os.path.getmtime)
    return None


def _merge_by_row_index(case_series):
    max_len = 0
    for _case_name, times, values in case_series:
        max_len = max(max_len, len(times), len(values))

    header = []
    for case_name, _times, _values in case_series:
        header.extend([f"{case_name}_time", f"{case_name}_pellet_total_strain_energy"])

    rows = []
    for i in range(max_len):
        row = []
        for _case_name, times, values in case_series:
            row.append(times[i] if i < len(times) else "")
            row.append(values[i] if i < len(values) else "")
        rows.append(row)
    return header, rows


def main():
    project_base_dir = os.environ.get("PROJECT_BASE_DIR", "")
    studies_subdir = os.environ.get("STUDIES_SUBDIR", "parameter_studies")
    fields = _parse_fields(os.environ.get("PV_FIELDS_SINGLE2", ""))
    if not fields:
        raise SystemExit("PV_FIELDS_SINGLE2 is empty")

    reduction = os.environ.get("PV_TIMESERIES_REDUCTION", "mean")
    output_dir_name = os.environ.get("PV_TIMESERIES_OUTPUT_DIR", "post_results")
    studies_dir = os.path.join(project_base_dir, studies_subdir)

    if "pellet_total_strain_energy" not in fields:
        export_fields = [*fields, "pellet_total_strain_energy"]
    else:
        export_fields = fields

    print("============================================================")
    print("Step 6 全时间域标量导出")
    print(f"PROJECT_BASE_DIR: {project_base_dir}")
    print(f"STUDIES_SUBDIR:   {studies_subdir}")
    print(f"输出目录名:       {output_dir_name}")
    print(f"字段列表:         {', '.join(export_fields)}")
    print(f"场量归约方式:     {reduction}")
    print("============================================================")

    case_dirs = sorted(glob.glob(os.path.join(studies_dir, "case_*")))
    if not case_dirs:
        raise SystemExit(f"未找到案例目录: {studies_dir}/case_*")

    merged_case_series = []

    for idx, case_dir in enumerate(case_dirs, start=1):
        case_name = os.path.basename(case_dir)
        exodus_path = _pick_primary_exodus(case_dir, output_dir_name)
        if not exodus_path:
            print(f"[{idx}/{len(case_dirs)}] {case_name}: 未找到 .e 文件，跳过")
            continue

        out_dir = os.path.join(case_dir, output_dir_name)
        os.makedirs(out_dir, exist_ok=True)
        base = os.path.splitext(os.path.basename(exodus_path))[0]
        out_csv = os.path.join(out_dir, f"{base}_step6_timeseries.csv")

        print(f"[{idx}/{len(case_dirs)}] {case_name}")
        print(f"  exodus: {os.path.basename(exodus_path)}")
        print(f"  csv:    {os.path.relpath(out_csv, project_base_dir) if project_base_dir else out_csv}")

        time_steps, series = export_exodus_timeseries(
            exodus_path,
            export_fields,
            out_csv,
            reduction=reduction,
        )

        energy = series.get("pellet_total_strain_energy", [])
        merged_case_series.append((case_name, time_steps, energy))
        print(f"  timesteps: {len(time_steps)}")

    if merged_case_series:
        merged_dir = os.path.join(studies_dir, output_dir_name)
        os.makedirs(merged_dir, exist_ok=True)
        merged_path = os.environ.get(
            "PV_TIMESERIES_MERGED_CSV",
            os.path.join(merged_dir, "step6_merged_pellet_total_strain_energy.csv"),
        )
        header, rows = _merge_by_row_index(merged_case_series)
        with open(merged_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        print("------------------------------------------------------------")
        print(f"汇总CSV: {merged_path}")
        # print("说明: 第1-2列为第1个案例(time, pellet_total_strain_energy)，第2个案例从第3-4列开始，以此类推")
        print("------------------------------------------------------------")


if __name__ == "__main__":
    main()
