"""Optional RGI integration for the ColabFold diffusion notebooks."""

SUPPORTED_MODELS = frozenset(
    {
        "alphafold3",
        "openbind0",
        "openfold3",
        "boltz2",
        "protenix2",
        "rosettafold3",
        "chai1",
        "opendde",
        "esmfold2",
        "esmfold2_lm600m",
        "esmfold2_lm300m",
    }
)
ALPHAFOLD3_COLABFOLD_VERSION = "3.1.11"


def require_supported_model(model):
    if model not in SUPPORTED_MODELS:
        raise ValueError(
            f"RGI is not supported for {model!r}. Choose "
            f"{', '.join(sorted(SUPPORTED_MODELS))}, or turn use_rgi off."
        )
