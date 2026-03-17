# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
from pathlib import Path
from typing import Optional

import biom
import numpy as np
import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.genome_data import LociDirectoryFormat
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from q2_types.per_sample_sequences import BAMDirFmt, BAMFormat, FastqGzFormat
from subprocess import run, CalledProcessError

_FEATURECOUNTS_STRAND_MODE = {
    "unstranded": "0",
    "forward": "1",
    "reverse": "2",
}

_FEATURECOUNTS_ANNOTATION_COLUMNS = {
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length",
}


# def _parse_casava_fastq_filename(path: Path) -> tuple[str, str]:
#     parts = path.name.rsplit("_", maxsplit=4)
#     if len(parts) != 5:
#         raise ValueError(
#             f"Unexpected FASTQ filename for CASAVA-formatted reads: {path.name!r}."
#         )
#
#     sample_id, _barcode_id, _lane, read_number, _suffix = parts
#     return sample_id, read_number
#
#
# def _collect_single_end_reads(reads) -> list[tuple[str, Path]]:
#     collected = []
#
#     for relpath, _ in reads.sequences.iter_views(FastqGzFormat):
#         path = Path(reads.path) / relpath
#         sample_id, read_number = _parse_casava_fastq_filename(path)
#         if read_number != "R1":
#             continue
#         collected.append((sample_id, path))
#
#     if not collected:
#         raise ValueError("No forward FASTQ files were found in the reads artifact.")
#
#     return sorted(collected, key=lambda item: item[0])
#
#
# def _collect_paired_end_reads(reads) -> list[tuple[str, Path, Path]]:
#     grouped = {}
#
#     for relpath, _ in reads.sequences.iter_views(FastqGzFormat):
#         path = Path(reads.path) / relpath
#         sample_id, read_number = _parse_casava_fastq_filename(path)
#         grouped.setdefault(sample_id, {})[read_number] = path
#
#     if not grouped:
#         raise ValueError("No FASTQ files were found in the reads artifact.")
#
#     collected = []
#     for sample_id in sorted(grouped):
#         pair = grouped[sample_id]
#         if "R1" not in pair or "R2" not in pair:
#             raise ValueError(
#                 f"Sample {sample_id!r} does not have both forward and reverse reads."
#             )
#         collected.append((sample_id, pair["R1"], pair["R2"]))
#
#     return collected


# def _collect_reference_paths(
#     reference_sequences: GenomeSequencesDirectoryFormat,
# ) -> list[Path]:
#     reference_paths = [
#         Path(reference_sequences.path) / relpath
#         for relpath, _ in reference_sequences.genomes.iter_views(DNAFASTAFormat)
#     ]
#
#     if not reference_paths:
#         raise ValueError(
#             "No reference FASTA files were found in the reference genome artifact."
#         )
#
#     return reference_paths
#
#
# def _run_subread_buildindex(
#     index_basename: Path,
#     reference_fasta_paths: list[Path],
# ) -> None:
#     cmd = [
#         "subread-buildindex",
#         "-o",
#         str(index_basename),
#         *(str(path) for path in reference_fasta_paths),
#     ]
#
#     try:
#         subprocess.run(
#             cmd,
#             check=True,
#             capture_output=True,
#             text=True,
#         )
#     except subprocess.CalledProcessError as exc:
#         detail = exc.stderr.strip() or exc.stdout.strip()
#         raise RuntimeError(
#             f"subread-buildindex failed with exit code {exc.returncode}: {detail}"
#         ) from exc


# def _run_subread_align(
#     *,
#     index_basename: Path,
#     read_1_path: Path,
#     output_path: Path,
#     read_2_path: Optional[Path],
#     threads: int,
#     indels: int,
#     multi_mapping: bool,
#     max_alignments: int,
#     minimum_fragment_length: Optional[int] = None,
#     maximum_fragment_length: Optional[int] = None,
# ) -> None:
#     cmd = [
#         "subread-align",
#         "-i",
#         str(index_basename),
#         "-r",
#         str(read_1_path),
#         "-o",
#         str(output_path),
#         "-T",
#         str(threads),
#         "-I",
#         str(indels),
#     ]
#
#     if read_2_path is not None:
#         cmd.extend(
#             [
#                 "-R",
#                 str(read_2_path),
#                 "-d",
#                 str(minimum_fragment_length),
#                 "-D",
#                 str(maximum_fragment_length),
#             ]
#         )
#
#     if multi_mapping:
#         cmd.extend(["--multiMapping", "-B", str(max_alignments)])
#
#     try:
#         subprocess.run(
#             cmd,
#             check=True,
#             capture_output=True,
#             text=True,
#         )
#     except subprocess.CalledProcessError as exc:
#         detail = exc.stderr.strip() or exc.stdout.strip()
#         raise RuntimeError(
#             f"subread-align failed with exit code {exc.returncode}: {detail}"
#         ) from exc


def _run_featurecounts(
    annotation_path: Path,
    output_path: Path,
    bam_paths: list[Path],
    feature_type: str,
    id_attribute: str,
    strand_mode: str,
    threads: int,
    paired_end: bool,
) -> None:
    cmd = [
        "featureCounts",
        "-F",
        # featureCounts uses "GTF" for both GTF and compatible GFF inputs.
        "GTF",
        "-a",
        str(annotation_path),
        "-o",
        str(output_path),
        "-t",
        feature_type,
        "-g",
        id_attribute,
        "-s",
        _FEATURECOUNTS_STRAND_MODE[strand_mode],
        "-T",
        str(threads),
    ]

    if paired_end:
        cmd.append("-p")

    cmd.extend(str(path) for path in bam_paths)

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
            f"featureCounts failed with exit code {exc.returncode}: {detail}"
        ) from exc


# def _build_alignment_result_dir(
#     sample_ids: list[str],
#     aligner,
# ) -> BAMDirFmt:
#     result = BAMDirFmt()
#
#     for sample_id in sample_ids:
#         output_path = result.bams.path_maker(sample_id=sample_id)
#         aligner(sample_id, output_path)
#
#     return result


def _collect_bam_paths(alignment_maps: BAMDirFmt) -> tuple[list[str], list[Path]]:
    sample_ids = []
    bam_paths = []

    for relpath, _ in alignment_maps.bams.iter_views(BAMFormat):
        path = Path(alignment_maps.path) / relpath
        sample_ids.append(path.stem.replace("_alignment", ""))
        bam_paths.append(path)

    if not bam_paths:
        raise ValueError("No BAM files were found in the alignment map artifact.")

    return sample_ids, bam_paths


def _stage_annotation(loci: LociDirectoryFormat, reference_id: str) -> Path:
    loci_path = Path(loci.path)
    fp_matches = sorted(
        path
        for path in loci_path.rglob("*.gff")
        if path.is_file() and path.stem == reference_id
    )
    if len(fp_matches) == 1:
        annotation_path = fp_matches[0]
    else:
        raise ValueError(
            f"Reference annotation {reference_id!r} was not found in the loci artifact."
        )

    return annotation_path


def _parse_feature_counts(
    output_path: Path,
    sample_ids: list[str],
    bam_paths: list[Path],
) -> biom.Table:
    counts = pd.read_csv(output_path, sep="\t", comment="#")
    feature_ids = counts["Geneid"].astype(str)

    bam_path_lookup = {}
    for sample_id, bam_path in zip(sample_ids, bam_paths):
        bam_path_lookup[str(bam_path)] = sample_id
        bam_path_lookup[str(bam_path.resolve())] = sample_id
        bam_path_lookup[bam_path.name] = sample_id

    observed_sample_columns = [
        column
        for column in counts.columns
        if column not in _FEATURECOUNTS_ANNOTATION_COLUMNS
    ]

    rename_map = {}
    for column in observed_sample_columns:
        sample_id = bam_path_lookup.get(column)
        if sample_id is None:
            sample_id = bam_path_lookup.get(Path(column).name)
        if sample_id is None:
            raise RuntimeError(
                "featureCounts returned an unexpected sample column: " f"{column!r}."
            )
        rename_map[column] = sample_id

    count_frame = counts.loc[:, observed_sample_columns].rename(columns=rename_map)
    count_frame.index = feature_ids
    count_frame = count_frame.groupby(level=0).sum()
    count_frame = count_frame.reindex(columns=sample_ids, fill_value=0)
    count_frame = count_frame.apply(pd.to_numeric, errors="raise").astype(int)

    return biom.Table(
        count_frame.to_numpy(dtype=np.int64),
        observation_ids=count_frame.index.tolist(),
        sample_ids=count_frame.columns.tolist(),
    )


def count_features(
    alignment_maps: BAMDirFmt,
    loci: LociDirectoryFormat,
    reference_id: str,
    feature_type: str = "gene",
    id_attribute: str = "ID",
    strand_mode: str = "unstranded",
    threads: int = 1,
    paired_end: bool = True,
) -> biom.Table:
    sample_ids, bam_paths = _collect_bam_paths(alignment_maps)

    with tempfile.TemporaryDirectory(prefix="q2-subread-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        output_path = tmpdir_path / "feature-counts.tsv"

        annotation_path = _stage_annotation(loci, reference_id)
        _run_featurecounts(
            annotation_path=annotation_path,
            output_path=output_path,
            bam_paths=bam_paths,
            feature_type=feature_type,
            id_attribute=id_attribute,
            strand_mode=strand_mode,
            threads=threads,
            paired_end=paired_end,
        )

        return _parse_feature_counts(output_path, sample_ids, bam_paths)


# def align_reads(
#     reads: object,
#     reference_sequences: GenomeSequencesDirectoryFormat,
#     threads: int = 1,
#     indels: int = 16,
#     multi_mapping: bool = False,
#     max_alignments: int = 3,
#     minimum_fragment_length: int = 50,
#     maximum_fragment_length: int = 600,
# ) -> BAMDirFmt:
#     if minimum_fragment_length > maximum_fragment_length:
#         raise ValueError(
#             "minimum_fragment_length cannot be greater than maximum_fragment_length."
#         )
#
#     sequence_paths = [
#         Path(reads.path) / relpath
#         for relpath, _ in reads.sequences.iter_views(FastqGzFormat)
#     ]
#     has_reverse_reads = any(
#         _parse_casava_fastq_filename(path)[1] == "R2" for path in sequence_paths
#     )
#
#     if has_reverse_reads:
#         collected_reads = _collect_paired_end_reads(reads)
#         read_lookup = {
#             sample_id: (forward_path, reverse_path)
#             for sample_id, forward_path, reverse_path in collected_reads
#         }
#         sample_ids = [sample_id for sample_id, _, _ in collected_reads]
#     else:
#         collected_reads = _collect_single_end_reads(reads)
#         read_lookup = {
#             sample_id: (forward_path, None)
#             for sample_id, forward_path in collected_reads
#         }
#         sample_ids = [sample_id for sample_id, _ in collected_reads]
#
#     reference_paths = _collect_reference_paths(reference_sequences)
#
#     with tempfile.TemporaryDirectory(prefix="q2-subread-") as tmpdir:
#         tmpdir_path = Path(tmpdir)
#         index_basename = tmpdir_path / "reference-index"
#
#         _run_subread_buildindex(index_basename, reference_paths)
#
#         def _align(sample_id: str, output_path: Path) -> None:
#             forward_path, reverse_path = read_lookup[sample_id]
#             _run_subread_align(
#                 index_basename=index_basename,
#                 read_1_path=forward_path,
#                 output_path=output_path,
#                 read_2_path=reverse_path,
#                 threads=threads,
#                 indels=indels,
#                 multi_mapping=multi_mapping,
#                 max_alignments=max_alignments,
#                 minimum_fragment_length=minimum_fragment_length,
#                 maximum_fragment_length=maximum_fragment_length,
#             )
#
#         return _build_alignment_result_dir(sample_ids, _align)
