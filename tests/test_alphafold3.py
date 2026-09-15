"""Backend routing and the parts of the AF3 path that do not need alphafold3."""
import json
from datetime import date
from pathlib import Path

import numpy as np
import pytest

from colabfold.backend import AF3_MODELS, get_backend, is_af3_model
from colabfold.batch import dropped_entities, warn_about_dropped_entities
from colabfold.input import classify_molecules


def test_af2_models_do_not_route_to_af3():
    pytest.importorskip("alphafold")
    for model_type in ["alphafold2", "alphafold2_ptm", "alphafold2_multimer_v3", "deepfold_v1"]:
        assert not is_af3_model(model_type)
        assert type(get_backend(model_type)).__name__ == "AF2Backend"


def test_af3_models_route_to_af3():
    for model_type in ["alphafold3", "protenix2", "boltz2", "chai1", "openfold3"]:
        assert is_af3_model(model_type)
    assert "protenix2" in AF3_MODELS


def test_af3_backend_says_how_to_install_itself():
    pytest.importorskip("colabfold.backend")
    try:
        import alphafold3  # noqa: F401
    except ModuleNotFoundError:
        with pytest.raises(RuntimeError, match="git\\+https://github.com/sokrypton/alphafold3"):
            get_backend("alphafold3")


def test_alphafold2_warns_about_entities_it_drops(caplog):
    protein, molecules = classify_molecules("MKVLLA:smiles|c1ccccc1|2:DNA|ACGT")
    queries = [("job", protein, None, molecules)]

    dropped = dropped_entities(queries, "alphafold2_ptm")
    assert dropped == {"job": {"SMILES": 2, "DNA": 1}}
    warn_about_dropped_entities(dropped, "alphafold2_ptm")
    assert "folds protein chains only" in caplog.text

    assert dropped_entities(queries, "protenix2") == {}


def test_smiles_keeps_its_case():
    protein, molecules = classify_molecules("MKVLLA:smiles|c1ccccc1:DNA|acgt")
    assert protein == ["MKVLLA"]
    assert [(m.name, seq) for m, seq, _ in molecules] == [("SMILES", "c1ccccc1"), ("DNA", "ACGT")]


def test_chain_ids_follow_the_alphabet():
    from colabfold.alphafold3.input import chain_id

    assert [chain_id(i) for i in [0, 1, 25, 26, 27]] == ["A", "B", "Z", "AA", "AB"]


def test_per_residue_plddt_averages_each_residue():
    from colabfold.alphafold3.predict import per_residue_plddt

    class Structure:
        atom_b_factor = np.array([90.0, 80.0, 70.0, 60.0, 50.0])
        chain_id = np.array(["A", "A", "A", "B", "B"])
        res_id = np.array([1, 1, 2, 1, 1])

    assert per_residue_plddt(Structure()) == [85.0, 70.0, 55.0]


def test_msa_coverage_encodes_a3m_insertions():
    from colabfold.alphafold3.plot import _encode, _rows

    rows = _rows(">q\nMKV\n>h\nMkKV\n")
    assert rows == ["MKV", "MkKV"]
    # the lower case insertion is dropped, so both rows are three columns long
    assert _encode(rows, 3).shape == (2, 3)


def test_paired_rows_get_a_species_id_alphafold3_will_pair_on():
    from colabfold.msa_pairing import check_pairing_preconditions, rewrite_paired_descriptions

    paired = [">101\nMKV\n>hit_a\nMKA\n", ">102\nGGG\n>hit_b\nGGA\n"]
    assert check_pairing_preconditions(paired, 2) == 2

    out = rewrite_paired_descriptions(paired, 2)
    first = [l for l in out[0].splitlines() if l.startswith(">")]
    second = [l for l in out[1].splitlines() if l.startswith(">")]
    assert first[1] == second[1], "row 1 of both chains must share a species id"
    assert first[0] != second[0] or first[0] == ">101", "the query must not be paired"

    import re
    regex = re.compile(r"(?:tr|sp)\|(?:[A-Z0-9]{6,10})(?:_\d+)?\|(?:[A-Z0-9]{1,10}_)(?P<SpeciesId>[A-Z0-9]{1,5})")
    assert regex.match(first[1][1:]), "alphafold3 must be able to parse the species id"


def test_pairing_refuses_when_a_chain_has_no_paired_block():
    from colabfold.msa_pairing import PairingError, check_pairing_preconditions

    with pytest.raises(PairingError, match="no paired block"):
        check_pairing_preconditions([">101\nMKV\n", ""], 2)
    with pytest.raises(PairingError, match="unequal depths"):
        check_pairing_preconditions([">101\nMKV\n>a\nMKA\n", ">102\nGGG\n"], 2)
    with pytest.raises(PairingError, match="one per chain"):
        check_pairing_preconditions([">101\nMKV\n"], 2)


def test_nothing_is_left_for_alphafold3_to_search():
    folding_input = pytest.importorskip("alphafold3.common.folding_input")
    from colabfold.alphafold3.input import no_af3_search, with_msas

    chains = [
        folding_input.ProteinChain(id="A", sequence="MKV", ptms=[], unpaired_msa=None,
                                   paired_msa=None, templates=None),
        folding_input.RnaChain(id="B", sequence="ACGU", modifications=[], unpaired_msa=None),
    ]
    fold_input = folding_input.Input(name="t", chains=chains, rng_seeds=[1])

    for out in (no_af3_search(fold_input), with_msas(fold_input, [">1\nMKV\n"], None)):
        for chain in out.chains:
            assert getattr(chain, "unpaired_msa", "") is not None
            assert getattr(chain, "paired_msa", "") is not None
            assert getattr(chain, "templates", ()) is not None


def test_build_fold_input_orders_chains_and_keeps_ligands():
    pytest.importorskip("alphafold3.common.folding_input")
    from colabfold.alphafold3.input import build_fold_input
    from colabfold.utils import MolType

    fold_input = build_fold_input(
        "job", ["MKV", "GGG"], [1, 2], None, None,
        molecules=[(MolType.SMILES, "c1ccccc1", 2), (MolType.DNA, "ACGT", 1)], seeds=[1, 2],
    )
    assert [(type(c).__name__, c.id) for c in fold_input.chains] == [
        ("ProteinChain", "A"), ("ProteinChain", "B"), ("ProteinChain", "C"),
        ("Ligand", "D"), ("Ligand", "E"), ("DnaChain", "F"),
    ]
    assert list(fold_input.rng_seeds) == [1, 2]
    assert fold_input.ligands[0].smiles == "c1ccccc1"


def _reference_kernel(q, k, v, mask_bias, nonbatched_bias, scale):
    """Stands in for tri_flash: same contract, plain jnp so it runs on CPU."""
    import jax
    import jax.numpy as jnp

    logits = jnp.einsum("bhqd,bhkd->bhqk", q, k) * scale
    if nonbatched_bias is not None:
        logits = logits + nonbatched_bias[None]
    logits = logits + mask_bias
    return jnp.einsum("bhqk,bhkd->bhqd", jax.nn.softmax(logits, axis=-1), v)


@pytest.mark.parametrize("lead", [(), (3,), (2, 3)])
def test_af3_attention_translation_matches_alphafold3(lead):
    jax = pytest.importorskip("jax")
    import jax.numpy as jnp

    from colabfold.alphafold3.attention import colabfold_attention

    rng = np.random.default_rng(0)
    seq_q, seq_k, heads, dim = 7, 13, 4, 16
    mk = lambda s: jnp.asarray(rng.normal(size=s), jnp.float32)
    q, k, v = mk(lead + (seq_q, heads, dim)), mk(lead + (seq_k, heads, dim)), mk(lead + (seq_k, heads, dim))
    bias = mk((heads, seq_q, seq_k))
    keep = rng.random((seq_k,)) > 0.3
    keep[0] = True
    mask = jnp.asarray(keep).reshape((1,) * (len(lead) + 2) + (seq_k,))
    scale = dim ** -0.5

    logits = jnp.einsum("...qhd,...khd->...hqk", q, k) * scale + bias
    want = jnp.einsum("...hqk,...khd->...qhd", jax.nn.softmax(jnp.where(mask, logits, -1e9), -1), v)

    got = colabfold_attention(q, k, v, mask=mask, bias=bias, scale=scale, kernel=_reference_kernel)
    assert got.shape == want.shape
    assert float(jnp.max(jnp.abs(got - want))) < 2e-5


def test_af3_attention_declines_a_per_row_bias():
    pytest.importorskip("jax")
    import jax.numpy as jnp

    from colabfold.alphafold3.attention import colabfold_attention

    q = jnp.zeros((2, 5, 4, 16))
    bias = jnp.zeros((2, 4, 5, 5))
    assert colabfold_attention(q, q, q, bias=bias, kernel=_reference_kernel) is None


@pytest.mark.parametrize("model_name", ["alphafold3", "protenix2", "openfold3", "boltz2"])
def test_featurises_against_the_real_package(model_name):
    pytest.importorskip("alphafold3.data.featurisation")
    pytest.importorskip("haiku")
    from colabfold.alphafold3.input import build_fold_input
    from colabfold.alphafold3.models import featurise

    fold_input = build_fold_input("t", ["MKVLLA"], [1], [">101\nMKVLLA\n>h\nMKVLLG\n"], None,
                                  seeds=[1])
    examples = featurise(fold_input, model_name, model_dir="/tmp/none")
    assert len(examples) == 1
    assert examples[0]["token_index"].shape[0] == 6


def test_ligands_become_tokens():
    pytest.importorskip("alphafold3.data.featurisation")
    from colabfold.alphafold3.input import build_fold_input
    from colabfold.alphafold3.models import featurise
    from colabfold.utils import MolType

    fold_input = build_fold_input("t", ["MKVLLA"], [1], [">101\nMKVLLA\n"], None,
                                  molecules=[(MolType.SMILES, "c1ccccc1", 1)], seeds=[1])
    examples = featurise(fold_input, "alphafold3", model_dir="/tmp/none")
    # six residues plus benzene's six carbons, which AlphaFold2 would have dropped
    assert examples[0]["token_index"].shape[0] == 12


def test_a_ccd_ligand_without_ideal_coordinates_featurises():
    pytest.importorskip("alphafold3.data.featurisation")
    from colabfold.alphafold3.input import build_fold_input
    from colabfold.alphafold3.models import featurise
    from colabfold.utils import MolType

    # TAC has no ideal coordinates, so alphafold3 dates its CCD entry to decide
    # whether to fall back to the model ones, and compares against this
    fold_input = build_fold_input("t", ["MKVLLA"], [1], [">101\nMKVLLA\n"], None,
                                  molecules=[(MolType.CCD, "TAC", 1)], seeds=[1])
    examples = featurise(fold_input, "alphafold3", model_dir="/tmp/none",
                         ref_max_modified_date=date.fromisoformat("2100-01-01"))
    assert examples[0]["token_index"].shape[0] > 6


def test_template_indices_compose_onto_the_seqres():
    """query -> hit -> SEQRES, checked by reading the residue letters back."""
    import dataclasses

    from colabfold.alphafold3.templates import map_hit_to_seqres, query_to_hit_map

    seqres = "MTTASPSQVRQNYHQDAEAAINRQINLELYASYVYLSMSYYFDRDDVALKNFAKYFLHQSHEE"
    template_part = seqres[20:30]
    query = "AAAAA" + template_part[:4] + "X" + template_part[4:] + "AAAAA"

    @dataclasses.dataclass
    class Hit:
        query: str
        hit_sequence: str
        indices_query: list
        indices_hit: list

    aligned_q = query[5:16]
    aligned_t = template_part[:4] + "-" + template_part[4:]
    indices_query, indices_hit, qi, hi = [], [], 5, 20
    for _, t in zip(aligned_q, aligned_t):
        indices_query.append(qi)
        qi += 1
        indices_hit.append(-1 if t == "-" else hi)
        hi += t != "-"
    hit = Hit(aligned_q, aligned_t, indices_query, indices_hit)

    hit_to_seqres = map_hit_to_seqres(hit.hit_sequence, seqres)
    mapping = {q: hit_to_seqres[h] for q, h in query_to_hit_map(hit, query).items()
               if h in hit_to_seqres}

    assert mapping == {5: 20, 6: 21, 7: 22, 8: 23, 10: 24, 11: 25, 12: 26, 13: 27, 14: 28, 15: 29}
    assert 9 not in mapping, "the query insertion has no template residue"
    assert all(query[q] == seqres[s] for q, s in mapping.items())


def test_attention_on_the_real_kernel_matches_fp32():
    jax = pytest.importorskip("jax")
    pytest.importorskip("alphafold.model.tri_flash")
    if jax.devices()[0].platform == "cpu":
        pytest.skip("needs a GPU")
    import jax.numpy as jnp

    from colabfold.alphafold3.attention import colabfold_attention

    rng = np.random.default_rng(0)
    seq, heads, dim = 64, 4, 32
    mk = lambda s: jnp.asarray(rng.normal(size=s) * 0.5, jnp.bfloat16)
    q, k, v = mk((1, seq, heads, dim)), mk((1, seq, heads, dim)), mk((1, seq, heads, dim))
    bias = jnp.asarray(rng.normal(size=(heads, seq, seq)) * 0.3, jnp.bfloat16)
    keep = rng.random((seq,)) > 0.25
    keep[0] = True
    mask = jnp.asarray(keep)[None, None, None, :]
    scale = dim ** -0.5

    logits = jnp.einsum("...qhd,...khd->...hqk", q.astype(jnp.float32), k.astype(jnp.float32))
    logits = jnp.where(mask, logits * scale + bias.astype(jnp.float32), -1e9)
    want = jnp.einsum("...hqk,...khd->...qhd", jax.nn.softmax(logits, -1), v.astype(jnp.float32))

    got = colabfold_attention(q, k, v, mask=mask, bias=bias, scale=scale)
    assert float(jnp.max(jnp.abs(got.astype(jnp.float32) - want))) < 5e-2


def test_af3_models_cite_their_weights():
    from colabfold.citations import af3_citations, citations

    assert af3_citations("alphafold2_ptm") == []
    assert af3_citations("openfold3") == ["Abramson2024", "OpenFold3", "OpenBind"]
    # the two Protenix releases have their own papers
    assert af3_citations("protenix1") == ["Abramson2024", "ProtenixV1"]
    assert af3_citations("protenix2") == ["Abramson2024", "ProtenixV2"]
    # every model must name the weights it runs on, and every key must resolve
    for model in ("openfold3", "openbind0", "protenix1", "protenix2", "boltz2", "chai1",
                  "intellifold2", "rosettafold3", "opendde", "esmfold2", "esmfold2_fast"):
        keys = af3_citations(model)
        assert len(keys) >= 2, f"{model} cites no weights"
        assert all(key in citations for key in keys)
    assert af3_citations("openbind0") == ["Abramson2024", "OpenFold3", "OpenBind"]


def test_a_json_input_keeps_the_msa_it_came_with():
    folding_input = pytest.importorskip("alphafold3.common.folding_input")
    from colabfold.alphafold3.input import msa_state, with_msas

    theirs = ">user\nMKV\n>hit\nMKA\n"
    chains = [
        folding_input.ProteinChain(id="A", sequence="MKV", ptms=[], unpaired_msa=theirs,
                                   paired_msa=""),
        folding_input.ProteinChain(id="B", sequence="MKV", ptms=[], unpaired_msa=None,
                                   paired_msa=None),
    ]
    fold_input = folding_input.Input(name="t", chains=chains, rng_seeds=[1])
    assert msa_state(fold_input) == "mixed"

    out = with_msas(fold_input, [">cf\nMKV\n"], None)
    assert out.protein_chains[0].unpaired_msa == theirs, "the file's own MSA must survive"
    assert out.protein_chains[1].unpaired_msa.startswith(">cf")


def test_copies_of_a_chain_share_one_msa():
    folding_input = pytest.importorskip("alphafold3.common.folding_input")
    from colabfold.alphafold3.input import with_msas

    chains = [folding_input.ProteinChain(id=i, sequence=s, ptms=[])
              for i, s in (("A", "MKV"), ("B", "MKV"), ("C", "AAW"))]
    fold_input = folding_input.Input(name="t", chains=chains, rng_seeds=[1])

    # the search returns one MSA per unique sequence, not one per chain
    out = with_msas(fold_input, [">cf\nMKV\n", ">cf\nAAW\n"], None)
    assert [c.unpaired_msa for c in out.protein_chains] == [
        ">cf\nMKV\n", ">cf\nMKV\n", ">cf\nAAW\n"]


def test_a_json_without_a_protein_still_makes_a_query(tmp_path):
    pytest.importorskip("alphafold3")
    from colabfold.alphafold3.input import polymer_lengths
    from colabfold.input import queries_from_af3_json

    path = tmp_path / "dna.json"
    path.write_text(json.dumps({
        "name": "dsDNA", "modelSeeds": [1], "dialect": "alphafold3", "version": 1,
        "sequences": [{"dna": {"id": "A", "sequence": "ACGTACGT"}},
                      {"dna": {"id": "B", "sequence": "ACGTACGT"}}],
    }))

    (name, sequence, _, extras), = queries_from_af3_json(path)
    assert name == "dsDNA"
    assert sequence == [], "nothing for the MSA search, but still a job"
    assert polymer_lengths(extras.fold_input) == [8, 8]


def test_af3_names_the_options_it_cannot_honour(caplog):
    pytest.importorskip("alphafold3")
    from colabfold.alphafold3.backend import AF3Backend
    from colabfold.backend import RunOptions

    backend = AF3Backend("openfold3")
    opts = RunOptions(model_type="openfold3", num_relax=1, rank_by="plddt")
    with caplog.at_level("WARNING"):
        backend.configure(opts, max_len=100, max_num=1, num_queries=1,
                          msa_mode="mmseqs2_uniref_env", is_complex=False, use_templates=False)
    assert "--amber" in caplog.text
    assert "--rank-by" not in caplog.text, "ranking is honoured now"

    caplog.clear()
    with caplog.at_level("WARNING"):
        backend.configure(RunOptions(model_type="openfold3"), max_len=100, max_num=1,
                          num_queries=1, msa_mode="mmseqs2_uniref_env", is_complex=False,
                          use_templates=False)
    assert "ignores" not in caplog.text, "defaults must not warn"


def test_multi_chain_coverage_plot_draws():
    pytest.importorskip("alphafold3.common.folding_input")
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from alphafold3.common import folding_input

    from colabfold.alphafold3.plot import plot_msa_coverage

    a3m = ">q\nMKVLLA\n>h1\nMKVLLG\n>h2\nMKVAAA\n"
    chains = [folding_input.ProteinChain(id=c, sequence="MKVLLA", ptms=[], unpaired_msa=a3m,
                                         paired_msa="") for c in ("A", "B")]
    fold_input = folding_input.Input(name="t", chains=chains, rng_seeds=[1])
    assert plot_msa_coverage(fold_input) is not None


def test_templates_reach_the_featurised_batch():
    pytest.importorskip("alphafold3.data.featurisation")
    pytest.importorskip("haiku")
    import numpy as np

    from alphafold3.common import folding_input

    from colabfold.alphafold3.input import build_fold_input
    from colabfold.alphafold3.models import featurise, mark_absl_flags_parsed

    mark_absl_flags_parsed()
    sequence = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    a3m = f">101\n{sequence}\n"

    def non_gap(templates):
        fold_input = build_fold_input("t", [sequence], [1], [a3m], None, seeds=[1],
                                      templates=templates)
        example = featurise(fold_input, "alphafold3", model_dir="/tmp/none")[0]
        return int(np.sum(example["template_aatype"] > 0))

    assert non_gap(None) == 0


def test_rank_key_picks_the_requested_metric():
    from colabfold.alphafold3.predict import rank_key

    scores = {"plddt": [80.0, 90.0], "ptm": 0.7, "iptm": 0.4, "ranking_score": 0.6}
    assert rank_key("auto", scores) == 0.6
    assert rank_key("ranking_score", scores) == 0.6
    assert rank_key("plddt", scores) == 85.0
    assert rank_key("ptm", scores) == 0.7
    assert rank_key("iptm", scores) == 0.4
    assert rank_key("multimer", scores) == 0.4  # iptm when there is one
    assert rank_key("multimer", {k: v for k, v in scores.items() if k != "iptm"}) == 0.7


def test_af3_only_warns_about_what_it_still_cannot_do(caplog):
    pytest.importorskip("alphafold3")
    from colabfold.alphafold3.backend import AF3Backend
    from colabfold.backend import RunOptions

    backend = AF3Backend("openfold3")
    # these are honoured now, so they must not be named
    opts = RunOptions(model_type="openfold3", rank_by="ptm", stop_at_score=90,
                      max_seq=256, save_pair_representations=True)
    with caplog.at_level("WARNING"):
        backend.configure(opts, max_len=100, max_num=1, num_queries=1,
                          msa_mode="mmseqs2_uniref_env", is_complex=False, use_templates=False)
    assert "ignores" not in caplog.text

    caplog.clear()
    with caplog.at_level("WARNING"):
        backend.configure(RunOptions(model_type="openfold3", initial_guess="x.pdb"),
                          max_len=100, max_num=1, num_queries=1,
                          msa_mode="mmseqs2_uniref_env", is_complex=False, use_templates=False)
    assert "--initial-guess" in caplog.text


def test_af3_ranking_is_not_forced_to_plddt():
    """The AlphaFold2 rule keys on the model name, which no alphafold3 model matches."""
    from colabfold.backend import is_af3_model

    for model_type in ("openfold3", "protenix2", "boltz2"):
        assert is_af3_model(model_type)
        assert "ptm" not in model_type and "multimer" not in model_type


def test_af3_weights_come_from_hugging_face(tmp_path):
    from colabfold.alphafold3.weights import model_dir_for, urls_for

    assert urls_for("protenix/protenix2.bin.zst", "sokrypton/af3-any-model") == [
        "https://huggingface.co/sokrypton/af3-any-model/resolve/main/"
        "protenix/protenix2.bin.zst"
    ]

    # beside AlphaFold2's params, and one directory per precision
    assert model_dir_for("protenix2", tmp_path) == tmp_path / "params/af3/protenix2"
    assert model_dir_for("protenix2", tmp_path, "int8") == tmp_path / "params/af3/protenix2-int8"


def test_a_dead_mirror_falls_back_to_the_next_url(tmp_path, monkeypatch):
    import colabfold.download as download

    tried = []

    class Response:
        headers = {"Content-Length": "5"}

        def raise_for_status(self):
            if "mirror" in tried[-1]:
                raise OSError("502 Bad Gateway")

        def iter_content(self, chunk_size):
            yield b"bytes"

    def get(url, **kwargs):
        tried.append(url)
        return Response()

    monkeypatch.setattr(download.requests, "get", get)
    dest = tmp_path / "weights.bin.zst"
    download.fetch_file(["https://mirror/x", "https://fallback/x"], dest, "test")

    assert tried == ["https://mirror/x", "https://fallback/x"]
    assert dest.read_bytes() == b"bytes"
    assert not list(tmp_path.glob("*.part")), "the partial file must not be left behind"
