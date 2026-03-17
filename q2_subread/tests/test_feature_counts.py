# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import patch

import pandas as pd
import qiime2 as q2
from q2_types.per_sample_sequences import BAMDirFmt
from qiime2.plugins import subread
from rachis.plugin.testing import TestPluginBase

from q2_subread import _methods


class TestFeatureCounts(TestPluginBase):
    package = "q2_subread.tests"

    def setUp(self):
        super().setUp()
        self.loci = q2.Artifact.import_data(
            "GenomeData[Loci]", self.get_data_path("loci")
        )
        self.reads = q2.Artifact.import_data(
            "SampleData[SequencesWithQuality]",
            self.get_data_path("reads"),
            "CasavaOneEightSingleLanePerSampleDirFmt",
        )
        self.ref_index = q2.Artifact.load(self.get_data_path("ref-index.qza"))
        self.maps = q2.Artifact.load(self.get_data_path("maps.qza"))

    def test_count_features_end_to_end(self):
        (obs_ft,) = subread.methods.count_features(
            alignment_maps=self.maps,
            loci=self.loci,
            paired_end=False,
            reference_id="ref",
        )
        exp_ft = pd.read_csv(self.get_data_path("counts.csv"), index_col=0)

        pd.testing.assert_frame_equal(obs_ft.view(pd.DataFrame), exp_ft)

    def test_count_features_end_to_end_subprocess_error(self):
        with patch(
            "q2_subread._methods.run",
            side_effect=subprocess.CalledProcessError(
                1, "featureCounts", stderr="annotation parsing failed"
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "annotation parsing failed"):
                subread.methods.count_features(
                    alignment_maps=self.maps,
                    loci=self.loci,
                    paired_end=False,
                    reference_id="ref",
                )

    def test_count_features_end_to_end_no_reference(self):
        loci = q2.Artifact.import_data(
            "GenomeData[Loci]", self.get_data_path("loci-wrong")
        )
        with self.assertRaisesRegex(ValueError, "annotation 'ref' was not found"):
            subread.methods.count_features(
                alignment_maps=self.maps,
                loci=loci,
                paired_end=False,
                reference_id="ref",
            )

    def test_count_features_end_to_end_unexpected_sample(self):
        def _mock_run_featurecounts(*, output_path, **kwargs):
            output_path.write_text(
                "\t".join(
                    [
                        "Geneid",
                        "Chr",
                        "Start",
                        "End",
                        "Strand",
                        "Length",
                        "unexpected_sample",
                    ]
                )
                + "\n"
                + "\t".join(["geneA", "chrToy", "1", "29", "+", "29", "3"])
                + "\n"
            )

        with patch(
            "q2_subread._methods._run_featurecounts",
            side_effect=_mock_run_featurecounts,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "featureCounts returned an unexpected sample column",
            ):
                subread.methods.count_features(
                    alignment_maps=self.maps,
                    loci=self.loci,
                    paired_end=False,
                    reference_id="ref",
                )

    def test_collect_bam_paths_no_bam_paths(self):
        with self.assertRaisesRegex(
            ValueError,
            "No BAM files were found in the alignment map artifact.",
        ):
            _methods._collect_bam_paths(BAMDirFmt())
