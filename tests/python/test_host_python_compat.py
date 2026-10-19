"""The helper scripts that run on the user's HOST Python (not in a container) must stay valid on
Python 3.8: no match statements, no `X | Y` outside annotations, no `zip(strict=)`, no `datetime.UTC`.

Scripts that only run inside containers (Python 3.12) may use newer syntax.
"""

import ast
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
HOST_SCRIPTS = ["make_samplesheet.py", "check_samplesheet.py", "replay_run.py", "reanalyse_runs.py"]
NEWER_ATTRIBUTES = {"UTC", "removeprefix", "removesuffix", "is_relative_to"}


def annotation_nodes(tree):
    """Ids of every node that sits inside an annotation (lazy with `from __future__ import annotations`)."""
    inside = set()
    for node in ast.walk(tree):
        annotations = []
        if isinstance(node, ast.arg) and node.annotation is not None:
            annotations.append(node.annotation)
        if isinstance(node, ast.AnnAssign):
            annotations.append(node.annotation)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.returns is not None:
            annotations.append(node.returns)
        for annotation in annotations:
            inside.update(id(n) for n in ast.walk(annotation))
    return inside


def test_host_scripts_parse_as_python_3_8():
    for name in HOST_SCRIPTS:
        # SyntaxError for a match statement, parenthesised context managers, ...
        ast.parse((BIN / name).read_text(), filename=name, feature_version=(3, 8))


def test_host_scripts_use_no_post_3_8_runtime_features():
    for name in HOST_SCRIPTS:
        tree = ast.parse((BIN / name).read_text(), filename=name)
        assert any(
            isinstance(node, ast.ImportFrom) and node.module == "__future__"
            for node in ast.walk(tree)
        ), f"{name}: needs `from __future__ import annotations` for X | Y annotations"
        lazy = annotation_nodes(tree)
        for node in ast.walk(tree):
            where = f"{name}:{getattr(node, 'lineno', '?')}"
            if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "zip":
                assert not any(k.arg == "strict" for k in node.keywords), f"{where}: zip(strict=)"
            if isinstance(node, ast.Attribute):
                assert node.attr not in NEWER_ATTRIBUTES, (
                    f"{where}: .{node.attr} needs Python >= 3.9"
                )
            if (
                isinstance(node, ast.BinOp)
                and isinstance(node.op, ast.BitOr)
                and id(node) not in lazy
            ):
                operands = (node.left, node.right)
                is_type_union = any(
                    isinstance(o, (ast.Subscript, ast.Name, ast.Attribute))
                    or (isinstance(o, ast.Constant) and o.value is None)
                    for o in operands
                ) and not any(
                    isinstance(o, ast.Constant) and isinstance(o.value, int) for o in operands
                )
                names_only = all(isinstance(o, (ast.Subscript, ast.Constant)) for o in operands)
                assert not (is_type_union and names_only), f"{where}: `X | Y` outside an annotation"
