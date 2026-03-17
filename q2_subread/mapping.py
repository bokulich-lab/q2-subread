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
from subprocess import run, CalledProcessError
from typing import Optional, Union

import pandas as pd
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt, SequencesWithQuality, PairedEndSequencesWithQuality,
)
from q2_types.sample_data import SampleData

from q2_subread._types import SubreadIndexDirFmt


def _run_subread_align(
    index_basename: Path,
    forward_fp: Path,
    reverse_fp: Optional[Path],
    output_path: Path,
    threads: int,
    indels: int,
    multi_mapping: bool,
    max_alignments: int,
    min_fragment_length: Optional[int] = None,
    max_fragment_length: Optional[int] = None,
    experiment_type: str = "rna-seq",
) -> None:
    cmd = [
        "subread-align",
        "-i",
        str(index_basename),
        "-r",
        str(forward_fp),
        "-o",
        str(output_path),
        "-T",
        str(threads),
        "-I",
        str(indels),
    ]
    if experiment_type == "rna-seq":
        cmd.extend(["-t", "0"])
    elif experiment_type == "dna-seq":
        cmd.extend(["-t", "1"])

    if reverse_fp is not None:
        cmd.extend(
            [
                "-R",
                str(reverse_fp),
                "-d",
                str(min_fragment_length),
                "-D",
                str(max_fragment_length),
            ]
        )

    if multi_mapping:
        cmd.extend(["--multiMapping", "-B", str(max_alignments)])

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
            f"subread-align failed with exit code {exc.returncode}: {detail}"
        ) from exc


def _build_alignment_result_dir(
    sample_ids: list[str],
    aligner,
) -> BAMDirFmt:
    result = BAMDirFmt()
    with tempfile.TemporaryDirectory(prefix="q2-subread-align-") as tmpdir:
        for sample_id in sample_ids:
            output_path = Path(tmpdir) / f"{sample_id}_alignment.bam"
            aligner(sample_id, output_path)

            shutil.move(output_path, Path(result.path) / output_path.name)
    return result


def _map_reads(
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt,
        SingleLanePerSampleSingleEndFastqDirFmt,
    ],
    reference_index: SubreadIndexDirFmt,
    threads: int = 1,
    indels: int = 16,
    multi_mapping: bool = False,
    max_alignments: int = 3,
    min_fragment_length: int = 50,
    max_fragment_length: int = 600,
    experiment_type: str = "rna-seq",
) -> BAMDirFmt:
    if min_fragment_length > max_fragment_length:
        raise ValueError(
            "minimum_fragment_length cannot be greater than maximum_fragment_length."
        )

    read_lookup = reads.manifest.view(pd.DataFrame).to_dict(orient="index")
    sample_ids = list(read_lookup.keys())

    index_basename = Path(reference_index.path) / "subread-index"

    def _align(sample_id: str, output_path: Path) -> None:
        forward_path = read_lookup[sample_id].get("forward")
        reverse_path = read_lookup[sample_id].get("reverse")
        _run_subread_align(
            index_basename=index_basename,
            forward_fp=forward_path,
            reverse_fp=reverse_path,
            output_path=output_path,
            threads=threads,
            indels=indels,
            multi_mapping=multi_mapping,
            max_alignments=max_alignments,
            min_fragment_length=min_fragment_length,
            max_fragment_length=max_fragment_length,
            experiment_type=experiment_type,
        )

    return _build_alignment_result_dir(sample_ids, _align)


def map_reads(
    ctx,
    reads,
    reference_index,
    threads = 1,
    indels = 16,
    multi_mapping = False,
    max_alignments = 3,
    min_fragment_length = 50,
    max_fragment_length = 600,
    experiment_type = "rna-seq",
    num_partitions = 1
):
    if reads.type <= SampleData[SequencesWithQuality]:
        partition_action = ctx.get_action("demux", "partition_samples_single")
    elif reads.type <= SampleData[PairedEndSequencesWithQuality]:
        partition_action = ctx.get_action("demux", "partition_samples_paired")
    else:
        raise NotImplementedError(f"Unsupported reads type: {reads.type}")

    map_action = ctx.get_action("subread", "_map_reads")
    collate_action = ctx.get_action("assembly", "collate_alignments")

    (partitioned_seqs,) = partition_action(reads, num_partitions)

    all_maps = []
    for seq in partitioned_seqs.values():
        maps, = map_action(
            reads=seq, reference_index=reference_index, threads=threads, indels=indels,
            multi_mapping=multi_mapping, max_alignments=max_alignments,
            min_fragment_length=min_fragment_length, max_fragment_length=max_fragment_length,
            experiment_type=experiment_type
        )
        all_maps.append(maps)

    (collated_maps,) = collate_action(all_maps)

    return collated_maps
