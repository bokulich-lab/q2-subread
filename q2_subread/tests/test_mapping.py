# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from pathlib import Path
from unittest.mock import call, patch

import pandas as pd
import qiime2 as q2
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2.plugins import subread
from rachis.plugin.testing import TestPluginBase

from q2_subread._types import SubreadIndexDirFmt
from q2_subread.mapping import _run_subread_align


class TestMapping(TestPluginBase):
    package = "q2_subread.tests"

    def setUp(self):
        super().setUp()
        self.reads = q2.Artifact.import_data(
            "SampleData[SequencesWithQuality]",
            self.get_data_path("reads"),
            "CasavaOneEightSingleLanePerSampleDirFmt",
        )
        self.ref_index = q2.Artifact.import_data(
            "SubreadIndex", self.get_data_path("subread-index")
        )
        read_dir = self.reads.view(SingleLanePerSampleSingleEndFastqDirFmt)
        self.read_lookup = read_dir.manifest.view(pd.DataFrame).to_dict(orient="index")

    def _assert_bam_outputs(self, result):
        observed = sorted(path.name for path in result.path.iterdir())
        self.assertEqual(
            observed,
            [
                "sample1_alignment.bam",
                "sample2_alignment.bam",
            ],
        )
        for path in result.path.iterdir():
            self.assertGreater(path.stat().st_size, 0)

    def test_map_reads_end_to_end(self):
        def _mock_align(**kwargs):
            kwargs["output_path"].write_bytes(b"bam-data")

        with patch(
            "q2_subread.mapping._run_subread_align",
            side_effect=_mock_align,
        ) as mock_align, patch(
            "q2_types.per_sample_sequences._formats.BAMFormat._validate_",
        ):
            (observed,) = subread.methods.map_reads(
                reads=self.reads,
                reference_index=self.ref_index,
                threads=4,
                indels=8,
                multi_mapping=True,
                max_alignments=5,
                experiment_type="dna-seq",
            )
            observed_dir = observed.view(BAMDirFmt)
            self._assert_bam_outputs(observed_dir)
            expected_index = Path(self.ref_index.view(SubreadIndexDirFmt).path) / (
                "subread-index"
            )
            self.assertEqual(mock_align.call_count, 2)
            mock_align.assert_has_calls(
                [
                    call(
                        index_basename=expected_index,
                        forward_fp=self.read_lookup["sample1"]["forward"],
                        reverse_fp=None,
                        output_path=mock_align.call_args_list[0].kwargs["output_path"],
                        threads=4,
                        indels=8,
                        multi_mapping=True,
                        max_alignments=5,
                        min_frag_length=50,
                        max_frag_length=600,
                        experiment_type="dna-seq",
                    ),
                    call(
                        index_basename=expected_index,
                        forward_fp=self.read_lookup["sample2"]["forward"],
                        reverse_fp=None,
                        output_path=mock_align.call_args_list[1].kwargs["output_path"],
                        threads=4,
                        indels=8,
                        multi_mapping=True,
                        max_alignments=5,
                        min_frag_length=50,
                        max_frag_length=600,
                        experiment_type="dna-seq",
                    ),
                ]
            )
            self.assertEqual(
                [
                    call_.kwargs["output_path"].name
                    for call_ in mock_align.call_args_list
                ],
                [
                    "sample1_alignment.bam",
                    "sample2_alignment.bam",
                ],
            )

    def test_map_reads_end_to_end_invalid_fragment_lengths(self):
        with self.assertRaisesRegex(
            ValueError,
            "min_frag_length cannot be greater than max_frag_length",
        ):
            subread.methods.map_reads(
                reads=self.reads,
                reference_index=self.ref_index,
                min_frag_length=601,
                max_frag_length=600,
            )

    def test_run_subread_align_builds_paired_end_command(self):
        with patch("q2_subread.mapping.run") as mock_run:
            _run_subread_align(
                index_basename=Path("/tmp/subread-index"),
                forward_fp=Path("/tmp/sample_R1.fastq.gz"),
                reverse_fp=Path("/tmp/sample_R2.fastq.gz"),
                output_path=Path("/tmp/sample_alignment.bam"),
                threads=8,
                indels=12,
                multi_mapping=True,
                max_alignments=6,
                min_frag_length=75,
                max_frag_length=500,
                experiment_type="dna-seq",
            )

        mock_run.assert_called_once_with(
            [
                "subread-align",
                "-i",
                "/tmp/subread-index",
                "-r",
                "/tmp/sample_R1.fastq.gz",
                "-o",
                "/tmp/sample_alignment.bam",
                "-T",
                "8",
                "-I",
                "12",
                "-t",
                "1",
                "-R",
                "/tmp/sample_R2.fastq.gz",
                "-d",
                "75",
                "-D",
                "500",
                "--multiMapping",
                "-B",
                "6",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    def test_run_subread_align_subprocess_error(self):
        with patch(
            "q2_subread.mapping.run",
            side_effect=subprocess.CalledProcessError(
                1,
                "subread-align",
                stderr="alignment failed",
            ),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "subread-align failed with exit code 1: alignment failed",
            ):
                _run_subread_align(
                    index_basename=Path("/tmp/subread-index"),
                    forward_fp=Path("/tmp/sample_R1.fastq.gz"),
                    reverse_fp=None,
                    output_path=Path("/tmp/sample_alignment.bam"),
                    threads=1,
                    indels=16,
                    multi_mapping=False,
                    max_alignments=3,
                )
