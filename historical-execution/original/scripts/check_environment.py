#!/usr/bin/env python3
"""Import-only environment inspection. Never loads weights or runs inference."""
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
for key, value in {
    "HF_HOME": str(ROOT / ".cache/huggingface"),
    "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1",
    "HF_HUB_DISABLE_TELEMETRY": "1", "WANDB_MODE": "disabled",
    "WANDB_DISABLED": "true", "OTEL_SDK_DISABLED": "true",
    "MPLCONFIGDIR": str(ROOT / ".cache/matplotlib"),
    "XDG_CACHE_HOME": str(ROOT / ".cache"),
    "PIP_CACHE_DIR": str(ROOT / ".cache/pip"),
    "NUMBA_CACHE_DIR": str(ROOT / ".cache/numba"),
    "DO_NOT_TRACK": "1",
}.items():
    os.environ[key] = value

def main():
    modules = {"profam": "profam", "esm": "fair-esm", "torch": "torch",
               "transformers": "transformers", "numpy": "numpy",
               "pandas": "pandas", "scipy": "scipy", "matplotlib": "matplotlib"}
    checks = {}
    for module, dist in modules.items():
        try:
            importlib.import_module(module)
            checks[module] = {"imported": True, "version": metadata.version(dist)}
        except Exception as exc:
            checks[module] = {"imported": False, "error": f"{type(exc).__name__}: {exc}"}
    pip = subprocess.run([sys.executable, "-m", "pip", "check"], capture_output=True, text=True)
    tools = {name: shutil.which(name) for name in ("latexmk", "lualatex", "bibtex", "pdftotext", "pdftoppm")}
    report = {"scope": "setup_imports_only", "python": platform.python_version(),
              "host": platform.node(), "executable": sys.executable,
              "platform": platform.platform(), "modules": checks,
              "pip_check_exit": pip.returncode, "pip_check": pip.stdout.strip(),
              "tools": tools, "model_weights_loaded": False, "inference_executed": False}
    if checks["torch"]["imported"]:
        import torch
        report["backends"] = {"mps_built": torch.backends.mps.is_built(),
                              "mps_available": torch.backends.mps.is_available(),
                              "cuda_available": torch.cuda.is_available()}
    report["passed"] = all(x["imported"] for x in checks.values()) and pip.returncode == 0 and all(tools.values())
    (ROOT / "environment/setup-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
