"""Shared ColabFold prediction options for the CAAT command-line wrappers."""

import argparse


def _model_order(value):
    try:
        order = [int(item) for item in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Expected comma-separated model numbers.") from exc
    if not order or any(item not in range(1, 6) for item in order):
        raise argparse.ArgumentTypeError("Model numbers must be between 1 and 5.")
    return order


def add_prediction_options(parser):
    """Register optional overrides; omitted values use batch.run's defaults."""
    group = parser.add_argument_group("additional ColabFold settings")
    options = [
        (("--num-recycles", "--num-recycle", "--num_recycles"), "num_recycles", int, "Number of prediction recycles (0 disables recycling)."),
        (("--recycle-early-stop-tolerance",), "recycle_early_stop_tolerance", float, "Recycling convergence tolerance."),
        (("--num-ensemble",), "num_ensemble", int, "Number of prediction ensembles."),
        (("--model-order",), "model_order", _model_order, "Comma-separated model order, e.g. 3,1,2."),
        (("--max-seq",), "max_seq", int, "Maximum number of MSA sequence clusters."),
        (("--max-extra-seq",), "max_extra_seq", int, "Maximum number of extra MSA sequences."),
        (("--custom-template-path",), "custom_template_path", str, "Custom template directory; requires --templates."),
        (("--stop-at-score",), "stop_at_score", float, "Stop prediction when this score is reached."),
    ]
    names = []
    for flags, dest, value_type, help_text in options:
        group.add_argument(*flags, dest=dest, type=value_type, default=argparse.SUPPRESS, help=help_text)
        names.append(dest)
    for flag, dest, choices in [
        ("--msa-mode", "msa_mode", ["mmseqs2_uniref_env", "mmseqs2_uniref_env_envpair", "mmseqs2_uniref", "single_sequence"]),
        ("--pair-mode", "pair_mode", ["unpaired", "paired", "unpaired_paired"]),
        ("--pair-strategy", "pairing_strategy", ["complete", "greedy"]),
        ("--rank", "rank_by", ["auto", "plddt", "ptm", "iptm", "multimer"]),
    ]:
        group.add_argument(flag, dest=dest, choices=choices, default=argparse.SUPPRESS)
        names.append(dest)
    for flag, dest, action in [
        ("--templates", "use_templates", "store_true"),
        ("--use-dropout", "use_dropout", "store_true"),
        ("--disable-cluster-profile", "use_cluster_profile", "store_false"),
        ("--save-recycles", "save_recycles", "store_true"),
        ("--save-all", "save_all", "store_true"),
    ]:
        group.add_argument(flag, dest=dest, action=action, default=argparse.SUPPRESS)
        names.append(dest)
    return tuple(names)


def prediction_kwargs(args, option_names):
    """Forward only explicitly supplied overrides to ColabFold."""
    return {name: getattr(args, name) for name in option_names if hasattr(args, name)}
