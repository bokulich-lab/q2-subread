# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import sys
from pathlib import Path

import qiime2 as q2
from q2_types.feature_data import DNAFASTAFormat
from qiime2.plugins import subread
from rachis.plugin.testing import TestPluginBase

from q2_subread._types import SubreadIndexDirFmt
from q2_subread.indexing import _run_subread_buildindex


class TestIndexing(TestPluginBase):
    package = "q2_subread.tests"

    def setUp(self):
        super().setUp()
        self._path = os.environ.get("PATH", "")
        env_bin = os.fspath(Path(sys.executable).resolve().parent)
        if env_bin not in self._path.split(os.pathsep):
            os.environ["PATH"] = (
                os.pathsep.join([env_bin, self._path]) if self._path else env_bin
            )
            self.addCleanup(os.environ.__setitem__, "PATH", self._path)
        self.reference = q2.Artifact.import_data(
            "FeatureData[Sequence]", self.get_data_path("ref.fasta")
        )

    def _assert_index_files(self, index, suffix_prefix):
        observed_paths = sorted(index.path.iterdir(), key=lambda path: path.name)
        observed = [path.name for path in observed_paths]
        self.assertEqual(
            observed,
            sorted(
                [
                    f"subread-index.00.{suffix_prefix}.array",
                    f"subread-index.00.{suffix_prefix}.tab",
                    "subread-index.lowinf",
                    "subread-index.reads",
                ]
            ),
        )
        for path in observed_paths:
            if path.name.endswith(".lowinf"):
                self.assertGreaterEqual(path.stat().st_size, 0)
            else:
                self.assertGreater(path.stat().st_size, 0)

    def test_build_index_end_to_end(self):
        (index,) = subread.methods.build_index(reference=self.reference)

        self._assert_index_files(index.view(SubreadIndexDirFmt), "b")

    def test_run_subread_buildindex_passes_optional_flags(self):
        reference = self.reference.view(DNAFASTAFormat)
        observed = _run_subread_buildindex(
            reference.path,
            F=True,
            B=True,
            M=1000,
            f=50,
            c=True,
        )

        self._assert_index_files(observed, "c")
