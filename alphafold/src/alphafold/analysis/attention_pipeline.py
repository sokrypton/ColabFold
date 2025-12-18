import sys
import logging
import typing as T
import numpy as np

from pathlib import Path
from scipy.stats import zscore

from alphafold.analysis.utils import _map_indices_to_aligned
from alphafold.analysis import analyze_residue, plot_difference, process_attention


logger = logging.getLogger(__name__)


def run_pipeline(
    query_seq_path: str,
    query_attn_dir: str,
    query_name: str,
    target_seq_path: T.Optional[str] = None,
    target_attn_dir: T.Optional[str] = None,
    target_name: T.Optional[str] = None,
    alignment_path: T.Optional[str] = None,
    save_path: str = "attention_visualizations",
    query_highlight_indices: T.Optional[T.List[int]] = None,
    target_highlight_indices: T.Optional[T.List[int]] = None,
    query_highlight_color: str = "#AE0639",
    target_highlight_color: str = "#1f77b4",
) -> None:
    logger.info("Reading query sequence file: %s", query_seq_path)
    query_sequence = process_attention.read_sequence_file(sequence_file=query_seq_path)
    logger.info("Query sequence length: %d", len(query_sequence))

    query_pos_highlights = (
        list(query_highlight_indices) if query_highlight_indices else []
    )
    if query_pos_highlights:
        logger.info("Query highlight indices provided: %s", query_pos_highlights)

    logger.info("Processing attention data for query: %s", query_name)
    query_n = process_attention.get_n(folder_path=query_attn_dir)
    logger.info("Query attention n value: %d", query_n)

    query_attn_spectrum = process_attention.get_attention(
        folder_path=query_attn_dir, n=query_n
    )
    logger.info(
        "Retrieved query attention spectrum shape: %s", np.shape(query_attn_spectrum)
    )

    logger.info("Averaging and normalizing query attention")
    query_attn_avg = process_attention.average(attention_spectrum=query_attn_spectrum)
    logger.info("Query attention average: %s", query_attn_avg)

    query_attn_min_max = process_attention.min_max(data=query_attn_avg)
    logger.info(
        "Shape of query attention after min-max: %s", np.shape(query_attn_min_max)
    )
    logger.info(
        "Query min-max stats: min=%g, max=%g, mean=%g",
        float(np.min(query_attn_min_max)),
        float(np.max(query_attn_min_max)),
        float(np.mean(query_attn_min_max)),
    )

    query_zscores = zscore(query_attn_min_max)
    logger.info("Computed z-scores for query attention (len=%d)", query_zscores.size)
    logger.debug(
        "Query z-score sample: %s", np.array2string(query_zscores[:10], precision=3)
    )

    query_important_indices = analyze_residue.find_important(
        attention=query_attn_min_max, zscores=query_zscores
    )
    logger.info(
        "Query important indices detected (%d): %s",
        len(query_important_indices) if query_important_indices is not None else 0,
        query_important_indices,
    )

    if not target_name and not target_attn_dir:
        output_subdir = Path(save_path) / query_name
        output_subdir.mkdir(parents=True, exist_ok=True)
        logger.info("Output directory created: %s", output_subdir)

        logger.info("Generating plots and analysis for single protein: %s", query_name)
        analyze_residue.find_highest_attention(
            attention=query_attn_avg,
            sequence=query_sequence,
            output_dir=str(output_subdir),
            protein=query_name,
        )

        logger.info("Plotting attention scores for query: %s", query_name)
        plot_difference.plot_attention(
            attention_scores=query_attn_min_max,
            highlighted_scores=query_important_indices,
            protein_name=query_name,
            output_dir=str(output_subdir),
            sequence=query_sequence,
            query_highlight_positions=query_pos_highlights,
            query_highlight_color=query_highlight_color,
        )
        logger.info(
            "Saved attention plot: %s",
            output_subdir / f"{query_name}_average_attention.png",
        )

    if target_name and target_attn_dir:
        logger.info("Reading target sequence file: %s", target_seq_path)
        target_sequence = process_attention.read_sequence_file(
            sequence_file=target_seq_path
        )
        logger.info("Target sequence length: %d", len(target_sequence))

        target_pos_highlights = (
            list(target_highlight_indices) if target_highlight_indices else []
        )
        logger.info("Target highlight indices provided: %s", target_pos_highlights)

        if (len(query_sequence) != len(target_sequence)) and not alignment_path:
            logger.error(
                "Sequence length mismatch (Query: %d, Target: %d). Alignment path required.",
                len(query_sequence),
                len(target_sequence),
            )
            sys.exit(1)

        logger.info(
            "Processing attention data for target: %s, %s", target_name, target_attn_dir
        )
        target_n = process_attention.get_n(folder_path=target_attn_dir)
        logger.info("Target attention n value: %d", target_n)

        logger.info("Getting target attention spectrum")
        target_attn_spectrum = process_attention.get_attention(
            folder_path=target_attn_dir, n=target_n
        )
        logger.info(
            "Retrieved target attention spectrum shape: %s",
            np.shape(target_attn_spectrum),
        )

        logger.info("Averaging and normalizing target attention")
        target_attn_avg = process_attention.average(
            attention_spectrum=target_attn_spectrum
        )
        logging.info("Shape of target attention average: %s", np.shape(target_attn_avg))
        logger.info(
            "Target attention average stats: min=%g, max=%g, mean=%g",
            float(np.min(target_attn_avg)),
            float(np.max(target_attn_avg)),
            float(np.mean(target_attn_avg)),
        )

        output_subdir = Path(save_path) / f"{query_name}_{target_name}"
        output_subdir.mkdir(parents=True, exist_ok=True)
        logger.info("Output directory created: %s", output_subdir)

        if len(query_sequence) == len(target_sequence):
            target_attn_min_max = process_attention.min_max(data=target_attn_avg)
            logging.info(
                "Shape of target attention after min-max: %s",
                np.shape(target_attn_min_max),
            )
            target_zscores = zscore(target_attn_min_max)
            target_important_indices = analyze_residue.find_important(
                attention=target_attn_min_max, zscores=target_zscores
            )
            logger.info(
                "Target important indices detected (%d): %s",
                len(target_important_indices)
                if target_important_indices is not None
                else 0,
                target_important_indices,
            )

            query_blosum, target_blosum = analyze_residue.blosum_scores(
                sequence1=query_sequence, sequence2=target_sequence
            )
            logger.info(
                "Computed BLOSUM scores for query and target (arrays lengths: %d, %d)",
                len(query_blosum),
                len(target_blosum),
            )

            query_diff, target_diff = analyze_residue.calculate_differences(
                scores1=query_blosum,
                scores2=target_blosum,
                attention1=query_attn_min_max,
                attention2=target_attn_min_max,
                gaps1=[],
                gaps2=[],
                bool_alignment=False,
            )
            logger.info(
                "Calculated attention differences for unaligned sequences (len=%d)",
                len(query_diff),
            )

            for attn, seq, name in [
                (query_attn_avg, query_sequence, query_name),
                (target_attn_avg, target_sequence, target_name),
            ]:
                analyze_residue.find_highest_attention(
                    attention=attn,
                    sequence=seq,
                    output_dir=str(output_subdir),
                    protein=name,
                )
                logger.info(
                    "Ran find_highest_attention for %s (results in %s)",
                    name,
                    output_subdir,
                )

            plot_difference.plot_attention(
                attention_scores=query_attn_min_max,
                highlighted_scores=query_important_indices,
                protein_name=query_name,
                output_dir=str(output_subdir),
                sequence=query_sequence,
                query_highlight_positions=query_pos_highlights,
                query_highlight_color=query_highlight_color,
            )
            plot_difference.plot_attention(
                attention_scores=target_attn_min_max,
                highlighted_scores=target_important_indices,
                protein_name=target_name,
                output_dir=str(output_subdir),
                sequence=target_sequence,
                target_highlight_positions=target_pos_highlights,
                target_highlight_color=target_highlight_color,
            )
            plot_difference.plot_difference(
                attn_diff_scores=query_diff,
                protein_name=query_name,
                output_dir=str(output_subdir),
                sequence=query_sequence,
                query_highlight_positions=query_pos_highlights,
                query_highlight_color=query_highlight_color,
            )
            plot_difference.plot_difference(
                attn_diff_scores=target_diff,
                protein_name=target_name,
                output_dir=str(output_subdir),
                sequence=target_sequence,
                target_highlight_positions=target_pos_highlights,
                target_highlight_color=target_highlight_color,
            )

        if alignment_path and target_name:
            aligned_seq_query, aligned_seq_target = process_attention.read_alignment(
                protein1=query_name, protein2=target_name, alignment_path=alignment_path
            )

            (
                attn_query_aligned,
                attn_target_aligned,
                query_gaps,
                target_gaps,
            ) = process_attention.align_attention(
                attention1=query_attn_avg,
                attention2=target_attn_avg,
                sequence1=aligned_seq_query,
                sequence2=aligned_seq_target,
            )

            query_aligned_mm = process_attention.min_max(data=attn_query_aligned)
            target_aligned_mm = process_attention.min_max(data=attn_target_aligned)

            query_aligned_z = zscore(query_aligned_mm)
            target_aligned_z = zscore(target_aligned_mm)

            query_aligned_imp = analyze_residue.find_important(
                attention=query_aligned_mm, zscores=query_aligned_z
            )
            target_aligned_imp = analyze_residue.find_important(
                attention=target_aligned_mm, zscores=target_aligned_z
            )

            query_aligned_blo, target_aligned_blo = analyze_residue.blosum_scores(
                sequence1=aligned_seq_query, sequence2=aligned_seq_target
            )

            query_aligned_pos = _map_indices_to_aligned(
                aligned_seq=aligned_seq_query, indices=query_pos_highlights
            )
            target_aligned_pos = _map_indices_to_aligned(
                aligned_seq=aligned_seq_target, indices=target_pos_highlights
            )

            (
                diff_query_aligned,
                diff_target_aligned,
            ) = analyze_residue.calculate_differences(
                scores1=query_aligned_blo,
                scores2=target_aligned_blo,
                attention1=query_aligned_mm,
                attention2=target_aligned_mm,
                gaps1=query_gaps,
                gaps2=target_gaps,
                bool_alignment=True,
            )

            analyze_residue.find_highest_attention(
                attention=attn_query_aligned,
                sequence=aligned_seq_query,
                output_dir=str(output_subdir),
                protein=query_name,
            )
            analyze_residue.find_highest_attention(
                attention=attn_target_aligned,
                sequence=aligned_seq_target,
                output_dir=str(output_subdir),
                protein=target_name,
            )

            plot_difference.plot_attention(
                attention=query_aligned_mm,
                highlighted_scores=query_aligned_imp,
                protein_name=query_name,
                output_dir=str(output_subdir),
                sequence=aligned_seq_query,
                query_highlight_positions=query_aligned_pos,
                query_highlight_color=query_highlight_color,
            )
            plot_difference.plot_attention(
                attention=target_aligned_mm,
                highlighted_scores=target_aligned_imp,
                protein_name=target_name,
                output_dir=str(output_subdir),
                sequence=aligned_seq_target,
                target_highlight_positions=target_aligned_pos,
                target_highlight_color=target_highlight_color,
            )

            plot_difference.plot_difference(
                attn_diff_scores=diff_query_aligned,
                protein_name=query_name,
                output_dir=str(output_subdir),
                sequence=aligned_seq_query,
                query_highlight_positions=query_aligned_pos,
                query_highlight_color=query_highlight_color,
            )
            plot_difference.plot_difference(
                attn_diff_scores=diff_target_aligned,
                protein_name=target_name,
                output_dir=str(output_subdir),
                sequence=aligned_seq_target,
                target_highlight_positions=target_aligned_pos,
                target_highlight_color=target_highlight_color,
            )
