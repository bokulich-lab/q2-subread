# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
from pathlib import Path
from subprocess import run, CalledProcessError

import biom
import numpy as np
import pandas as pd
from q2_types.genome_data import LociDirectoryFormat
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    BAMFormat,
)

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
