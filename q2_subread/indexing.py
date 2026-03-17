# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
import tempfile
from pathlib import Path

from q2_types.feature_data import DNAFASTAFormat
from subprocess import run, CalledProcessError

from q2_subread._types import SubreadIndexDirFmt


def _run_subread_buildindex(
    reference_fp: Path,
    F: bool,
    B: bool,
    M: int,
    f: int,
    c: bool,
) -> SubreadIndexDirFmt:
    index = SubreadIndexDirFmt()

    with tempfile.TemporaryDirectory(prefix="q2-subread-index-") as tmpdir:
        cmd = [
            "subread-buildindex",
            "-o",
            str(Path(tmpdir) / "subread-index"),
        ]

        if F:
            cmd.append("-F")
        if B:
            cmd.append("-B")
        if M is not None:
            cmd.extend(["-M", str(M)])
        if f is not None:
            cmd.extend(["-f", str(f)])
        if c:
            cmd.append("-c")

        cmd.append(str(reference_fp))

        try:
            run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
        except CalledProcessError as exc:
            detail = exc.stderr.strip() or exc.stdout.strip()
            raise RuntimeError(
                f"subread-buildindex failed with exit code {exc.returncode}: {detail}"
            ) from exc

        for index_fp in Path(tmpdir).glob("subread-index*"):
            if index_fp.suffix in (".array", ".tab", ".lowinf", ".reads"):
                shutil.move(index_fp, index.path / index_fp.name)

    return index


def build_index(
    reference: DNAFASTAFormat,
    full: bool = False,
    block: bool = True,
    memory: int = 8000,
    threshold: int = 100,
    color: bool = False,
) -> SubreadIndexDirFmt:
    return _run_subread_buildindex(
        reference.path,
        F=full,
        B=block,
        M=memory,
        f=threshold,
        c=color,
    )
