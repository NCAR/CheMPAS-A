#!/usr/bin/env python3
"""Generate isolated TMT bindings without modifying the dependency checkout."""

import re
import sys
from pathlib import Path


def generate(source: Path, output: Path) -> None:
    headers = [source / "include" / path for path in (
        "musica/musica.h", "musica/micm.h", "musica/tuvx.h",
        "tmt/prepare.h", "tmt/fable.h",
    )]
    declarations = "\n".join(path.read_text() for path in headers)
    declarations = re.sub(r"/\*.*?\*/", "", declarations, flags=re.S)
    symbols = set(re.findall(r"\b(\w+)\s*\([^;{}]*\)\s*;", declarations))
    if not {"DeleteString", "CreateMicm", "CreateTuvx"} <= symbols:
        raise ValueError("Unrecognized TMT public headers")

    output.mkdir(parents=True, exist_ok=True)
    (output / "namespace.h").write_text(
        "/* Generated: keep TMT allocation and ownership separate from MUSICA. */\n"
        + "\n".join(f"#define {name} ChempasTmt_{name}" for name in sorted(symbols))
        + "\n"
    )
    for path in sorted((source / "fortran").rglob("*.F90")):
        text = path.read_text()
        # Only module declarations/imports change; public type and method names,
        # include paths, constants, and the numerical implementation stay intact.
        text = re.sub(
            r"(?im)^(\s*(?:end\s+module|module|use)\s+)musica(?=\b|_)",
            r"\1tmt", text,
        )

        def binding(match: re.Match) -> str:
            name = match[2]
            if name not in symbols:
                raise ValueError(f"{path}: undeclared C binding {name}")
            return f'{match[1]}"ChempasTmt_{name}"'

        text = re.sub(r'(?i)(bind\s*\(\s*c\s*,\s*name\s*=\s*)[\"\'](\w+)[\"\']', binding, text)
        destination = output / path.relative_to(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text)


if __name__ == "__main__":
    generate(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
