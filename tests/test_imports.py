"""Tests that the package imports only the standard library, and never the network modules.

Two checks: the import statements in the source (static), and the modules that
are actually loaded when every module of the package is imported in a fresh
interpreter (some standard library modules load the network modules themselves:
xml.sax.saxutils loads urllib.request, http.client and socket).
Each check is also run on a copy of the package with an import added, and
must find it.
"""
import ast
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = "smiles2pptx"
NETWORK = {"urllib", "http", "socket", "ssl", "ftplib", "smtplib", "poplib", "imaplib",
           "xmlrpc", "socketserver", "webbrowser"}
# Loaded modules that are allowed although they are under a network package: urllib.parse
# only works on strings, and pathlib imports it in Python 3.10 and 3.11.
ALLOWED_AT_RUNTIME = {"urllib", "urllib.parse"}


def static_problems(root: Path) -> list:
    """Return the disallowed imports in the package source, as (file, module, reason)."""
    problems = []
    for path in sorted((root / PACKAGE).rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                names = [node.module]
            else:
                continue
            for name in names:
                top = name.split(".")[0]
                rel = path.relative_to(root).as_posix()
                if top == PACKAGE:
                    continue
                if top in NETWORK:
                    problems.append((rel, name, "network"))
                elif top not in sys.stdlib_module_names:
                    problems.append((rel, name, "not in the standard library"))
    return problems


_LOADER = r"""
import importlib, json, os, sys
root, package = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
before = set(sys.modules)
names = []
for folder, dirs, files in os.walk(os.path.join(root, package)):
    dirs[:] = [d for d in dirs if d != "__pycache__"]
    rel = os.path.relpath(folder, root).replace(os.sep, ".")
    for f in sorted(files):
        if f.endswith(".py"):
            names.append(rel if f == "__init__.py" else rel + "." + f[:-3])
for name in sorted(names):
    importlib.import_module(name)
print(json.dumps(sorted(set(sys.modules) - before)))
"""


def runtime_problems(root: Path) -> list:
    """Return the disallowed modules that importing every module of the package loads."""
    # -I: ignore the environment and the user site directory; -S: no site-packages,
    # so a module from outside the standard library cannot be found at all.
    r = subprocess.run([sys.executable, "-I", "-S", "-c", _LOADER, str(root), PACKAGE],
                       capture_output=True, text=True, encoding="utf-8")
    if r.returncode:
        last = (r.stderr.strip().splitlines() or ["(no message)"])[-1]
        return [("import failed", last)]
    problems = []
    for name in json.loads(r.stdout):
        top = name.split(".")[0]
        if top == PACKAGE or name in ALLOWED_AT_RUNTIME:
            continue
        if top in NETWORK:
            problems.append((name, "network"))
        elif top not in sys.stdlib_module_names:
            problems.append((name, "not in the standard library"))
    return problems


def broken_copy(tmp_path: Path, line: str) -> Path:
    """Copy the package into `tmp_path`, append `line` to the pptx writer and return `tmp_path`."""
    shutil.copytree(ROOT / PACKAGE, tmp_path / PACKAGE, ignore=shutil.ignore_patterns("__pycache__"))
    target = tmp_path / PACKAGE / "writers" / "pptx.py"
    target.write_text(target.read_text(encoding="utf-8") + "\n" + line + "\n", encoding="utf-8", newline="\n")
    return tmp_path


def test_source_imports_only_the_standard_library():
    assert static_problems(ROOT) == []


def test_loaded_modules_are_only_the_standard_library():
    assert runtime_problems(ROOT) == []


@pytest.mark.parametrize("line, static_finds, runtime_finds", [
    ("import urllib.request", "urllib.request", "urllib"),
    ("import fontTools", "fontTools", "import failed"),
    # Allowed by the static check, but it loads urllib and socket: only the runtime check sees it.
    ("from xml.sax.saxutils import escape", None, "urllib"),
])
def test_checks_find_an_added_import(tmp_path, line, static_finds, runtime_finds):
    root = broken_copy(tmp_path, line)
    found_static = [name for _, name, _ in static_problems(root)]
    found_runtime = [p[0] for p in runtime_problems(root)]
    if static_finds:
        assert found_static == [static_finds]
    else:
        assert found_static == []
    assert any(name == runtime_finds or name.startswith(runtime_finds + ".") for name in found_runtime), found_runtime
