"""
Run an alphafold3-open model and write results in ColabFold's layout.
"""
import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

logger = logging.getLogger(__name__)


def per_residue_plddt(structure) -> List[float]:
    """Mean atom b-factor per residue, in chain order."""
    b = np.asarray(structure.atom_b_factor, dtype=float)
    chain = np.asarray(structure.chain_id)
    res = np.asarray(structure.res_id)
    if len(b) == 0:
        return []
    keys = np.asarray([f"{c}\t{r}" for c, r in zip(chain, res)])
    starts = np.flatnonzero(np.concatenate(([True], keys[1:] != keys[:-1])))
    return [float(v) for v in np.add.reduceat(b, starts) / np.diff(np.append(starts, len(b)))]


def scores_of(inference_result, ranking_score: float) -> Dict[str, Any]:
    from alphafold3.model import confidence_types

    full = confidence_types.StructureConfidenceFull.from_inference_result(inference_result)
    meta = inference_result.metadata
    pae = np.asarray(full.pae, dtype=float)
    scores = {
        "plddt": np.around(per_residue_plddt(inference_result.predicted_structure), 2).tolist(),
        "pae": np.around(pae, 2).tolist(),
        "max_pae": float(pae.max()),
        "ranking_score": float(ranking_score),
        "token_chain_ids": list(full.token_chain_ids),
    }
    for key, name in (("predicted_tm_score", "ptm"),
                      ("interface_predicted_tm_score", "iptm"),
                      ("has_clash", "has_clash"),
                      ("fraction_disordered", "fraction_disordered")):
        value = meta.get(key)
        if value is not None and np.isfinite(float(value)):
            scores[name] = float(value)
    return scores


def _print_line(tag: str, scores: Dict[str, Any], took: float) -> str:
    from colabfold.backend import metrics_line

    summary = {k: scores[k] for k in ("ptm", "iptm", "ranking_score") if k in scores}
    if scores["plddt"]:
        summary["mean_plddt"] = float(np.mean(scores["plddt"]))
    return f"{tag} took {took:.1f}s" + metrics_line(summary)


def rank_key(rank_by: str, scores: Dict[str, Any]) -> float:
    """The metric to rank on. ``auto`` is alphafold3's own ranking score."""
    if rank_by in ("auto", "ranking_score"):
        return float(scores["ranking_score"])
    if rank_by == "plddt":
        return float(np.mean(scores["plddt"])) if scores["plddt"] else float("nan")
    if rank_by == "multimer":
        rank_by = "iptm" if "iptm" in scores else "ptm"
    return float(scores.get(rank_by, scores["ranking_score"]))


def predict_structure(
    prefix: str,
    result_dir: Path,
    fold_input,
    model_runner,
    model_type: str,
    featurised_examples: List[Any],
    save_all: bool = False,
    rank_by: str = "auto",
    stop_at_score: float = 100.0,
    prediction_callback=None,
) -> Dict[str, Any]:
    import jax
    from alphafold3.model import post_processing

    ranked, best = [], -float("inf")
    for seed, example in zip(fold_input.rng_seeds, featurised_examples):
        if best >= stop_at_score:
            logger.info(f"stopping early, {rank_by} reached {best:.3f}")
            break
        start = time.time()
        result = model_runner.run_inference(example, jax.random.PRNGKey(seed))
        inference_results = model_runner.extract_inference_results(
            batch=example, result=result, target_name=fold_input.name
        )
        took = time.time() - start
        for sample, inference_result in enumerate(inference_results):
            processed = post_processing.post_process_inference_result(inference_result)
            tag = f"{model_type}_seed_{seed:03d}_sample_{sample}"
            scores = scores_of(inference_result, processed.ranking_score)
            ranked.append((rank_key(rank_by, scores), tag, processed, scores, took,
                           model_runner.extract_embeddings(result, len(scores["plddt"]))))
            best = max(best, ranked[-1][0])
            if prediction_callback is not None:
                prediction_callback(inference_result.predicted_structure, None,
                                    scores, example, (tag, False))

    logger.info(f"reranking models by '{rank_by}' metric")
    ranked.sort(key=lambda row: row[0], reverse=True)

    rank, metric, result_files = [], [], []
    for n, (_, tag, processed, scores, took, embeddings) in enumerate(ranked):
        new_tag = f"rank_{(n + 1):03d}_{tag}"
        rank.append(new_tag)
        metric.append(scores)
        logger.info(_print_line(new_tag, scores, took))
        cif = result_dir.joinpath(f"{prefix}_{new_tag}.cif")
        cif.write_bytes(processed.cif)
        result_files.append(cif)
        scores_file = result_dir.joinpath(f"{prefix}_scores_{new_tag}.json")
        scores_file.write_text(json.dumps(scores))
        result_files.append(scores_file)
        if save_all:
            full = result_dir.joinpath(f"{prefix}_confidences_{new_tag}.json")
            full.write_bytes(processed.structure_full_data_json)
            result_files.append(full)
        if embeddings:
            npz = result_dir.joinpath(f"{prefix}_representations_{new_tag}.npz")
            np.savez_compressed(npz, **embeddings)
            result_files.append(npz)
    return {"rank": rank, "metric": metric, "result_files": result_files}
