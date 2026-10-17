"""A branch that cannot be named fails the substituent it hangs from instead of dropping out of the name."""

import pytest

import openclatura
from openclatura import name_smiles, opsin_available, verify_with_opsin
from openclatura.graph_io import read_smiles
from openclatura.namer import name_subgraph

# Silicon is named, and a stannane or bismuthane parent is unaffected.
GUARDS = {
    "C[Si](C)(C)C=Cc1cccs1": "2-(2-(trimethylsilyl)ethenyl)thiophene",
    "C[Si](C)(C)C#Cc1ccccn1": "2-(2-(trimethylsilyl)ethynyl)pyridine",
    "C[Si](C)(C)CCOc1ccccn1": "2-((2-(trimethylsilyl)ethyl)oxy)pyridine",
    "C[Sn](C)(C)C#Cc1ccccc1": "trimethyl(2-phenylethynyl)stannane",
    "C[Bi](C)CCc1ccccc1": "dimethyl(2-phenylethyl)bismuthane",
    "CCOCCc1ccccn1": "2-(2-ethoxyethyl)pyridine",
}

# Bismuth and aluminium have no substituent prefix. These were named without
# the metal (2-ethylpyridine, 2-ethoxypyridine, ...), as the trimethylstannyl
# alkyne was named 2-ethynylpyridine (lamalab-org/openclatura#129).
UNNAMEABLE = [
    "C[Bi](C)CCc1ccccn1",
    "C[Bi](C)C#Cc1ccccn1",
    "C[Bi](C)CCOc1ccccn1",
    "C[Bi](C)CC[Si](C)(C)c1ccccn1",
    "C[Al](C)C=Cc1cccs1",
]


@pytest.mark.parametrize("smiles", sorted(GUARDS))
def test_branch_naming_guards(smiles):
    assert name_smiles(smiles) == GUARDS[smiles]


@pytest.mark.parametrize("smiles", UNNAMEABLE)
def test_unnamed_nested_branch_fails_the_name(smiles):
    result = openclatura.name(smiles)
    assert not result.ok
    assert "UnnamedAtomError" in result.error


@pytest.mark.parametrize("smiles", ["[2H]OCCc1ccccn1", "[2H]C([2H])([2H])Oc1ccccn1", "C[Si](C)(C)[2H]"])
def test_explicit_hydrogen_branch_does_not_fail_the_name(smiles):
    # An explicit hydrogen has no prefix; the component audit skips it as well.
    assert openclatura.name(smiles).ok


def test_chain_branch_with_unnamed_branch_has_no_name():
    # The ethyl branch on the pyridine carries a dimethylbismuthanyl group,
    # which has no prefix.
    mol = read_smiles("C[Bi](C)CCc1ccccn1")
    assert name_subgraph(mol, 4, set(range(5, 11)), upstream_atom=5) == ""


def test_heteroatom_branch_with_unnamed_ligand_has_no_name():
    # Not "hydroxy": the oxygen's ligand is there, it just has no name.
    mol = read_smiles("C[Bi](C)Oc1ccccn1")
    assert name_subgraph(mol, 3, set(range(4, 10)), upstream_atom=4) == ""


@pytest.mark.opsin
@pytest.mark.parametrize("smiles", sorted(GUARDS))
def test_branch_names_round_trip(smiles):
    if not opsin_available():
        pytest.skip("OPSIN and Java are required")
    check = verify_with_opsin(GUARDS[smiles], smiles)
    assert check.status == "matched", check.to_dict()
