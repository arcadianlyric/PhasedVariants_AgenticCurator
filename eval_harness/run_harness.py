#!/usr/bin/env python3
"""
Deterministic evaluation harness for PhasedVariants AgenticCurator outputs.

The runner grades saved artifacts only. It does not call LLM APIs, download data,
or mutate project inputs.
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TASK_DIR = PROJECT_ROOT / "eval_harness" / "tasks"
DEFAULT_OUTPUT = PROJECT_ROOT / "eval_harness" / "results" / "harness_summary.json"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def iter_task_files(task_path: Path) -> Iterable[Path]:
    if task_path.is_file():
        yield task_path
        return
    yield from sorted(task_path.glob("*.json"))


def get_path(data: Any, dotted_path: str) -> Any:
    current = data
    for part in dotted_path.split("."):
        if isinstance(current, list):
            current = current[int(part)]
        elif isinstance(current, dict):
            current = current[part]
        else:
            raise KeyError(f"Cannot resolve {part!r} in {dotted_path!r}")
    return current


def normalize_text(value: Any) -> str:
    return str(value).casefold()


def grade_check(artifact: Dict[str, Any], check: Dict[str, Any]) -> Dict[str, Any]:
    check_type = check["type"]
    name = check["name"]
    observed = get_path(artifact, check["path"])

    if check_type == "json_path_min":
        threshold = check["min"]
        passed = float(observed) >= float(threshold)
        detail = f"observed={observed} min={threshold}"
    elif check_type == "json_path_equals":
        expected = check["value"]
        passed = observed == expected
        detail = f"observed={observed!r} expected={expected!r}"
    elif check_type == "text_contains_all":
        text = normalize_text(observed)
        missing = [term for term in check["terms"] if term.casefold() not in text]
        passed = not missing
        detail = "all terms present" if passed else f"missing={missing}"
    elif check_type == "text_absent_all":
        text = normalize_text(observed)
        present = [term for term in check["terms"] if term.casefold() in text]
        passed = not present
        detail = "all terms absent" if passed else f"present={present}"
    else:
        raise ValueError(f"Unsupported check type: {check_type}")

    return {
        "name": name,
        "type": check_type,
        "passed": passed,
        "detail": detail,
    }


def grade_task(task_file: Path) -> Dict[str, Any]:
    task = load_json(task_file)
    prediction_path = PROJECT_ROOT / task["prediction_path"]
    checks: List[Dict[str, Any]] = []

    if not prediction_path.exists():
        return {
            "id": task["id"],
            "gene": task.get("gene"),
            "passed": False,
            "score": 0.0,
            "prediction_path": str(prediction_path.relative_to(PROJECT_ROOT)),
            "checks": [{
                "name": "prediction_artifact_exists",
                "type": "file_exists",
                "passed": False,
                "detail": "prediction artifact not found",
            }],
        }

    artifact = load_json(prediction_path)
    for check in task["checks"]:
        checks.append(grade_check(artifact, check))

    passed_count = sum(1 for check in checks if check["passed"])
    total = len(checks)
    return {
        "id": task["id"],
        "gene": task.get("gene"),
        "passed": passed_count == total,
        "score": round(passed_count / total, 3) if total else 0.0,
        "prediction_path": str(prediction_path.relative_to(PROJECT_ROOT)),
        "checks": checks,
    }


def run(task_path: Path, output_path: Path, write_output: bool) -> Dict[str, Any]:
    task_files = list(iter_task_files(task_path))
    results = [grade_task(task_file) for task_file in task_files]
    passed_count = sum(1 for result in results if result["passed"])
    summary = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "task_count": len(results),
        "passed_count": passed_count,
        "failed_count": len(results) - passed_count,
        "accuracy": round(passed_count / len(results), 3) if results else 0.0,
        "results": results,
    }

    if write_output:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, ensure_ascii=False)

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run deterministic evaluation harness.")
    parser.add_argument("--tasks", default=str(DEFAULT_TASK_DIR), help="Task JSON file or directory.")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT), help="Output summary JSON path.")
    parser.add_argument("--no-write", action="store_true", help="Print summary without writing output JSON.")
    args = parser.parse_args()

    summary = run(Path(args.tasks), Path(args.output), write_output=not args.no_write)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if summary["failed_count"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
