import logging

from colabfold.batch import mk_hhsearch_db

HEADER = """data_{name}
#
loop_
_chem_comp.id
_chem_comp.type
LYS "L-peptide linking"
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 LYS
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . LYS A 1 1 ? 0.000 0.000 0.000 1.00 0.00 1 A 1
ATOM 2 C CA . LYS A 1 1 ? 1.458 0.000 0.000 1.00 0.00 1 A 1
ATOM 3 C C . LYS A 1 1 ? 2.009 1.420 0.000 1.00 0.00 1 A 1
ATOM 4 O O . LYS A 1 1 ? 1.251 2.390 0.000 1.00 0.00 1 A 1
"""

# like some PDB entries (e.g. 4udf, 1ejg), this one has an atom defined twice
DUPLICATE_ATOM = "ATOM 5 C C . LYS A 1 1 ? 2.100 1.500 0.000 1.00 0.00 1 A 1\n"


def test_mk_hhsearch_db_skips_unparsable_template(caplog, tmp_path):
    """A template that biopython cannot parse must not abort the run (#610)"""
    caplog.set_level(logging.WARNING)
    tmp_path.joinpath("good.cif").write_text(HEADER.format(name="good"))
    tmp_path.joinpath("broken.cif").write_text(
        HEADER.format(name="broken") + DUPLICATE_ATOM
    )

    mk_hhsearch_db(str(tmp_path))

    assert "broken.cif" in caplog.text
    a3m = tmp_path.joinpath("pdb70_a3m.ffdata").read_text()
    assert ">good_A" in a3m
    assert ">broken_A" not in a3m
