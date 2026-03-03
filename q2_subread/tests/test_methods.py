# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from q2_subread._methods import (
    align_reads,
    feature_counts,
)
from q2_subread.plugin_setup import plugin


class _DummyCollection:
    def __init__(self, relpaths):
        self._relpaths = [Path(path) for path in relpaths]

    def iter_views(self, view_type):
        for relpath in self._relpaths:
            yield relpath, None


class _DummyAlignmentMaps:
    def __init__(self, root, relpaths):
        self.path = Path(root)
        self.bams = _DummyCollection(relpaths)


class _DummyLoci:
    def __init__(self, root, relpaths):
        self.path = Path(root)
        self.loci = _DummyCollection(relpaths)


class _DummyReads:
    def __init__(self, root, relpaths):
        self.path = Path(root)
        self.sequences = _DummyCollection(relpaths)


class _DummyReferenceSequences:
    def __init__(self, root, relpaths):
        self.path = Path(root)
        self.genomes = _DummyCollection(relpaths)


class AlignmentMethodTests(unittest.TestCase):
    def test_align_reads_builds_index_and_aligns_each_single_end_sample(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads_dir = root / "reads"
            reference_dir = root / "reference"
            reads_dir.mkdir()
            reference_dir.mkdir()

            (reads_dir / "sample-1_S1_L001_R1_001.fastq.gz").write_bytes(b"FASTQ")
            (reads_dir / "sample-2_S1_L001_R1_001.fastq.gz").write_bytes(b"FASTQ")
            (reference_dir / "chr1.fasta").write_text(">chr1\nACGT\n")
            (reference_dir / "chr2.fasta").write_text(">chr2\nTGCA\n")

            reads = _DummyReads(
                reads_dir,
                [
                    "sample-1_S1_L001_R1_001.fastq.gz",
                    "sample-2_S1_L001_R1_001.fastq.gz",
                ],
            )
            reference_sequences = _DummyReferenceSequences(
                reference_dir,
                ["chr1.fasta", "chr2.fasta"],
            )

            captured = {"commands": []}

            def _run(cmd, check, capture_output, text):
                captured["commands"].append(cmd)
                if cmd[0] == "subread-buildindex":
                    captured["reference_paths"] = cmd[3:]
                elif cmd[0] == "subread-align":
                    output_path = Path(cmd[cmd.index("-o") + 1])
                    output_path.write_bytes(b"BAM")
                return subprocess.CompletedProcess(cmd, 0, "", "")

            with patch("q2_subread._methods.subprocess.run", side_effect=_run):
                observed = align_reads(
                    reads,
                    reference_sequences,
                    threads=6,
                    indels=12,
                    multi_mapping=True,
                    max_alignments=5,
                )

            observed_paths = sorted(p.name for p in Path(observed.path).glob("*.bam"))
            self.assertEqual(observed_paths, ["sample-1.bam", "sample-2.bam"])
            self.assertEqual(captured["commands"][0][0], "subread-buildindex")
            self.assertEqual(captured["commands"][1][0], "subread-align")
            self.assertEqual(captured["commands"][2][0], "subread-align")
            self.assertEqual(
                captured["reference_paths"],
                [
                    str(reference_dir / "chr1.fasta"),
                    str(reference_dir / "chr2.fasta"),
                ],
            )
            self.assertIn("--multiMapping", captured["commands"][1])
            self.assertEqual(
                captured["commands"][1][captured["commands"][1].index("-B") + 1],
                "5",
            )
            self.assertEqual(
                captured["commands"][1][captured["commands"][1].index("-T") + 1],
                "6",
            )
            self.assertEqual(
                captured["commands"][1][captured["commands"][1].index("-I") + 1],
                "12",
            )

    def test_align_reads_passes_paired_end_fragment_constraints(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads_dir = root / "reads"
            reference_dir = root / "reference"
            reads_dir.mkdir()
            reference_dir.mkdir()

            (reads_dir / "sample-1_S1_L001_R1_001.fastq.gz").write_bytes(b"FASTQ")
            (reads_dir / "sample-1_S1_L001_R2_001.fastq.gz").write_bytes(b"FASTQ")
            (reference_dir / "chr1.fasta").write_text(">chr1\nACGT\n")

            reads = _DummyReads(
                reads_dir,
                [
                    "sample-1_S1_L001_R1_001.fastq.gz",
                    "sample-1_S1_L001_R2_001.fastq.gz",
                ],
            )
            reference_sequences = _DummyReferenceSequences(
                reference_dir,
                ["chr1.fasta"],
            )

            captured = {"commands": []}

            def _run(cmd, check, capture_output, text):
                captured["commands"].append(cmd)
                if cmd[0] == "subread-align":
                    output_path = Path(cmd[cmd.index("-o") + 1])
                    output_path.write_bytes(b"BAM")
                return subprocess.CompletedProcess(cmd, 0, "", "")

            with patch("q2_subread._methods.subprocess.run", side_effect=_run):
                observed = align_reads(
                    reads,
                    reference_sequences,
                    minimum_fragment_length=75,
                    maximum_fragment_length=800,
                )

            observed_paths = sorted(p.name for p in Path(observed.path).glob("*.bam"))
            self.assertEqual(observed_paths, ["sample-1.bam"])
            align_cmd = captured["commands"][1]
            self.assertEqual(align_cmd[0], "subread-align")
            self.assertEqual(align_cmd[align_cmd.index("-d") + 1], "75")
            self.assertEqual(align_cmd[align_cmd.index("-D") + 1], "800")
            self.assertIn("-R", align_cmd)

    def test_align_reads_raises_useful_error_on_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads_dir = root / "reads"
            reference_dir = root / "reference"
            reads_dir.mkdir()
            reference_dir.mkdir()

            (reads_dir / "sample-1_S1_L001_R1_001.fastq.gz").write_bytes(b"FASTQ")
            (reference_dir / "chr1.fasta").write_text(">chr1\nACGT\n")

            reads = _DummyReads(reads_dir, ["sample-1_S1_L001_R1_001.fastq.gz"])
            reference_sequences = _DummyReferenceSequences(
                reference_dir,
                ["chr1.fasta"],
            )

            def _run(cmd, check, capture_output, text):
                if cmd[0] == "subread-align":
                    raise subprocess.CalledProcessError(
                        1,
                        cmd,
                        stderr="index not found",
                    )
                return subprocess.CompletedProcess(cmd, 0, "", "")

            with patch("q2_subread._methods.subprocess.run", side_effect=_run):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "subread-align failed with exit code 1: index not found",
                ):
                    align_reads(reads, reference_sequences)

    def test_align_reads_rejects_inverted_fragment_bounds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads_dir = root / "reads"
            reference_dir = root / "reference"
            reads_dir.mkdir()
            reference_dir.mkdir()

            (reads_dir / "sample-1_S1_L001_R1_001.fastq.gz").write_bytes(b"FASTQ")
            (reference_dir / "chr1.fasta").write_text(">chr1\nACGT\n")

            reads = _DummyReads(reads_dir, ["sample-1_S1_L001_R1_001.fastq.gz"])
            reference_sequences = _DummyReferenceSequences(
                reference_dir,
                ["chr1.fasta"],
            )

            with self.assertRaisesRegex(
                ValueError,
                "minimum_fragment_length cannot be greater than maximum_fragment_length",
            ):
                align_reads(
                    reads,
                    reference_sequences,
                    minimum_fragment_length=500,
                    maximum_fragment_length=100,
                )


class FeatureCountsMethodTests(unittest.TestCase):
    def test_feature_counts_builds_expected_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam_dir = root / "bams"
            loci_dir = root / "loci"
            bam_dir.mkdir()
            loci_dir.mkdir()

            for name in ("sample-1.bam", "sample-2.bam"):
                (bam_dir / name).write_bytes(b"BAM")

            (loci_dir / "chr1.gff").write_text(
                "##gff-version 3\n"
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-1\n"
            )
            (loci_dir / "chr2.gff").write_text(
                "##gff-version 3\n"
                "chr2\tsrc\tgene\t50\t200\t.\t-\t.\tID=gene-2\n"
            )

            alignment_maps = _DummyAlignmentMaps(
                bam_dir,
                ["sample-1.bam", "sample-2.bam"],
            )
            loci = _DummyLoci(loci_dir, ["chr1.gff", "chr2.gff"])

            captured = {}

            def _run(cmd, check, capture_output, text):
                captured["cmd"] = cmd
                annotation_path = Path(cmd[cmd.index("-a") + 1])
                output_path = Path(cmd[cmd.index("-o") + 1])
                captured["annotation"] = annotation_path.read_text()
                output_path.write_text(
                    "# Program:featureCounts v2.0.0\n"
                    "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
                    f"{bam_dir / 'sample-1.bam'}\t{bam_dir / 'sample-2.bam'}\n"
                    "gene-1\tchr1\t1\t100\t+\t100\t7\t9\n"
                    "gene-2\tchr2\t50\t200\t-\t151\t3\t4\n"
                )
                return subprocess.CompletedProcess(cmd, 0, "", "")

            with patch("q2_subread._methods.subprocess.run", side_effect=_run):
                observed = feature_counts(
                    alignment_maps,
                    loci,
                    feature_type="gene",
                    id_attribute="ID",
                    strand_mode="reverse",
                    threads=4,
                    paired_end=True,
                )

            table = pd.DataFrame(
                observed.matrix_data.toarray().T,
                index=observed.ids(axis="sample"),
                columns=observed.ids(axis="observation"),
            ).astype(int)

            expected = pd.DataFrame(
                [[7, 3], [9, 4]],
                index=["sample-1", "sample-2"],
                columns=["gene-1", "gene-2"],
            )

            pd.testing.assert_frame_equal(table, expected)
            self.assertEqual(captured["cmd"][0], "featureCounts")
            self.assertIn("-p", captured["cmd"])
            self.assertEqual(captured["cmd"][captured["cmd"].index("-s") + 1], "2")
            self.assertEqual(captured["cmd"][captured["cmd"].index("-T") + 1], "4")
            self.assertEqual(captured["annotation"].count("##gff-version 3"), 1)
            self.assertIn("ID=gene-1", captured["annotation"])
            self.assertIn("ID=gene-2", captured["annotation"])

    def test_feature_counts_raises_useful_error_on_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam_dir = root / "bams"
            loci_dir = root / "loci"
            bam_dir.mkdir()
            loci_dir.mkdir()

            (bam_dir / "sample-1.bam").write_bytes(b"BAM")
            (loci_dir / "chr1.gff").write_text(
                "##gff-version 3\n"
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-1\n"
            )

            alignment_maps = _DummyAlignmentMaps(bam_dir, ["sample-1.bam"])
            loci = _DummyLoci(loci_dir, ["chr1.gff"])

            error = subprocess.CalledProcessError(
                1,
                ["featureCounts"],
                stderr="annotation parsing failed",
            )

            with patch(
                "q2_subread._methods.subprocess.run",
                side_effect=error,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "featureCounts failed with exit code 1: annotation parsing failed",
                ):
                    feature_counts(alignment_maps, loci)


class PluginSetupTests(unittest.TestCase):
    def test_plugin_metadata(self):
        self.assertEqual(plugin.name, "subread")
        self.assertEqual(plugin.package, "q2_subread")
        self.assertEqual(
            plugin.description,
            "QIIME 2 plugin for wrapping the Subread software package.",
        )

    def test_feature_counts_registered(self):
        self.assertEqual(len(plugin.methods), 2)
