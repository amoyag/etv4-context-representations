#!/usr/bin/env python3
"""Check evidence receipts and hashes; never executes a scientific stage.

Mechanical integrity is necessary, not a certification of scientific merit.
Missing/pending stages deliberately fail. Review actors/evidence are auditable
claims, not proof of independence inferred from a boolean.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    return json.loads(path.read_text())

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def confined(root, relative):
    rel = Path(relative)
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    path = (root / rel).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes workspace: {relative}")
    return path

def check_stage(root, stage, policies, seen=None):
    seen = set() if seen is None else set(seen)
    if stage in seen:
        raise ValueError(f"Cyclic stage dependency: {stage}")
    seen.add(stage)
    policy = policies[stage]
    for dep in policy["requires"]:
        check_stage(root, dep, policies, seen)
    receipt_path = root / "handoff" / f"{stage}.json"
    if not receipt_path.is_file():
        raise ValueError(f"Stage {stage}: no receipt; not executed/verified")
    receipt = read(receipt_path)
    if receipt.get("stage") != stage or receipt.get("status") != "complete":
        raise ValueError(f"Stage {stage}: incomplete")
    if not receipt.get("actor") or not receipt.get("completed_at"):
        raise ValueError(f"Stage {stage}: missing actor/time")
    entries = receipt.get("outputs", [])
    output_map = {x["path"]: x["sha256"] for x in entries}
    if len(output_map) != len(entries):
        raise ValueError(f"Stage {stage}: duplicate output")
    for name in policy["outputs"]:
        if name not in output_map:
            raise ValueError(f"Stage {stage}: required output not registered: {name}")
    for name, expected in output_map.items():
        path = confined(root, name)
        if not path.is_file() or path.stat().st_size == 0 or digest(path) != expected:
            raise ValueError(f"Stage {stage}: missing/empty/changed output: {name}")
    checks = {x["id"]: x for x in receipt.get("checks", [])}
    for name in policy["checks"]:
        check = checks.get(name, {})
        if check.get("passed") is not True or not check.get("evidence"):
            raise ValueError(f"Stage {stage}: check not passed: {name}")
        if not all(p in output_map for p in check["evidence"]):
            raise ValueError(f"Stage {stage}: unhashed check evidence: {name}")
    if stage == "E1":
        qualification = read(root / "environment/qualification.json")
        config = read(root / "design/execution-config.json")
        if qualification.get("passed") is not True or qualification.get("biological_forwards") != 0:
            raise ValueError("E1: no successful artificial-only numerical qualification")
        if set(qualification.get("models", {})) != {"esm2_650m", "profam1"} or config.get("status") != "frozen":
            raise ValueError("E1: missing core model or unfrozen execution configuration")
        if not qualification.get("resources", {}).get("fits_before_freeze"):
            raise ValueError("E1: resource estimate did not fit the recorded freeze deadline")
        for model, data in qualification["models"].items():
            if not data.get("errors") or config["models"][model]["eta"] != data["eta"]:
                raise ValueError("E1: missing or inconsistent numerical evidence")
            for layer, error in data["errors"].items():
                expected = max(1e-7, 10 * max(error["repeat_max"], error["cpu_mps_max"]))
                if data["eta"][layer] != expected:
                    raise ValueError("E1: tolerance differs from the frozen numerical rule")
        for array in qualification.get("raw_arrays", []):
            if output_map.get(array["path"]) != array["sha256"]:
                raise ValueError("E1: raw numerical evidence is not in the receipt")
        if qualification.get("raw_array_count", 0) != sum(qualification.get("forwards", {}).values()) or not qualification.get("raw_array_count"):
            raise ValueError("E1: no actual model forwards or inconsistent execution count")
        for item in config.get("frozen_files", []):
            path = confined(root, item["path"])
            if digest(path) != item["sha256"]:
                raise ValueError("E1: frozen execution dependency changed: " + item["path"])
    if stage == "E2":
        coverage = read(root / "results/coverage.json")
        if coverage.get("core_complete") is not True or coverage.get("missing_core") != []:
            raise ValueError("E2: incomplete core run matrix")
        if not isinstance(coverage.get("completed_runs"), int) or coverage["completed_runs"] < 1:
            raise ValueError("E2: no actual completed runs")
    return receipt

def check_final(root, policies):
    check_stage(root, "R0", policies)
    decision = read(root / "handoff/review-ready.json")
    if decision.get("status") != "ready_for_author_review":
        raise ValueError("No review-ready decision")
    if decision.get("open_critical") != 0 or decision.get("open_major") != 0:
        raise ValueError("Unresolved central review findings")
    if decision.get("contribution_supported") is not True:
        raise ValueError("Scientific contribution not supported")
    disposition = read(root / "reviews/disposition.json")
    if not isinstance(disposition.get("findings"), list):
        raise ValueError("Missing findings register")
    for finding in disposition["findings"]:
        if finding.get("severity") in ("critical", "major"):
            if finding.get("status") != "closed_verified" or not finding.get("verification"):
                raise ValueError("Central finding lacks verified closure")
    if decision.get("posted") is not False:
        raise ValueError("This goal's package must remain unposted")
    reviewers = decision.get("reviewers", [])
    for role in ("scientific", "evidence"):
        if not any(r.get("role") == role and r.get("actor") and r.get("report") for r in reviewers):
            raise ValueError(f"Missing identified {role} reviewer")
    final_receipt = read(root / "handoff/R0.json")
    final_paths = {x["path"] for x in final_receipt["outputs"]}
    for reviewer in reviewers:
        if reviewer["report"] not in final_paths:
            raise ValueError("Reviewer report is not in final hashed outputs")
    for path in (root / "paper/sections").glob("*.tex"):
        if re.search(r"SCAFFOLD|TODO|TBD|PLACEHOLDER|NOT EXECUTED", path.read_text(), re.I):
            raise ValueError(f"Unresolved manuscript marker: {path.name}")
    main = (root / "paper/main.tex").read_text()
    if "SCAFFOLD" in main or "NOT EXECUTED" in main:
        raise ValueError("Manuscript is still a scaffold")

def setup(root):
    for name in ("AGENTS.md", "README.md", "PLAN.md", "design/protocol-v1.md",
                 "design/stages.json", "WORK_STATE.json", "sources/manifest.json",
                 "environment/requirements.lock.txt", "paper/main.tex"):
        path = root / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing setup file: {name}")
    state = read(root / "WORK_STATE.json")
    if state["execution_activated"] or any(s["status"] != "pending" for s in state["stages"].values()):
        raise ValueError("Setup-only check requires all execution stages pending")
    if len(list((root / "prompts").glob("*.xml"))) != 7:
        raise ValueError("Expected coordinator plus six stage prompts")
    for prompt in (root / "prompts").glob("*.xml"):
        if ET.parse(prompt).getroot().tag != "goal_spec":
            raise ValueError(f"Bad goal specification: {prompt.name}")
    for source in read(root / "sources/manifest.json")["sources"]:
        path = Path(source["path"])
        if not path.is_file() or digest(path) != source["sha256"]:
            raise ValueError(f"Changed source: {path}")
    if read(root / "environment/setup-report.json").get("passed") is not True:
        raise ValueError("Environment imports/dependency checks failed")
    if read(root / "environment/template-check.json").get("passed") is not True:
        raise ValueError("Template compilation not checked")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--setup", action="store_true")
    choice.add_argument("--stage", choices=("E0", "E1", "E2", "W0", "W1", "R0"))
    choice.add_argument("--final", action="store_true")
    args = parser.parse_args()
    try:
        policies = read(ROOT / "design/stages.json")["stages"]
        if args.setup:
            setup(ROOT)
        elif args.final:
            check_final(ROOT, policies)
        else:
            check_stage(ROOT, args.stage, policies)
    except (ValueError, KeyError, FileNotFoundError, ET.ParseError) as exc:
        print(f"NOT READY: {exc}", file=sys.stderr)
        return 2
    print("PASS: mechanical evidence checks only; see scientific review for interpretation")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
