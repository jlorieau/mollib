"""Pytests for the PDBReader."""

from math import sqrt

from mollib import Molecule
from mollib.core import MoleculeReader


def test_single_model():
    """Tests the loading of a single model into a molecule."""

    # Use the '2KXA' molecule as a test, which has 10 models. This test
    # is conducted by measuring the RMSD.
    model1 = Molecule('2kxa')  # by default, model 1 is loaded

    for model_id in range(1, 11):
        model = Molecule('2kxa', model_id=model_id)
        assert model.model_id == model_id

        # Calculate the RMSD
        count = 0
        msd = 0.0
        for atom1, atom2 in zip(model1.atoms, model.atoms):
            msd += sum(atom1.pos - atom2.pos) ** 2
            count += 1

        rmsd = sqrt(msd / (count - 1))

        # If it's the same model, they should have the same coordinates
        if model_id == 1:
            assert rmsd < 0.001
        # If it's a different model, the RMSD should be at least 0.2A
        else:
            assert rmsd > 0.2


def test_select_models():
    """Tests the loading of select models"""

    model_ids = [1, 5, 8]
    mr = MoleculeReader()
    models = mr.read(identifiers_or_files='2kxa', model_ids=model_ids)

    assert len(models) == 3

    for model_id, model in zip(model_ids, models):
        assert model.model_id == model_id
        assert len(list(model.atoms)) == len(list(models[0].atoms))


def test_all_models():
    """Tests the loading of all models."""

    # The '2kxa' structure should have 10 models
    mr = MoleculeReader()
    models = mr.read('2kxa')

    assert len(models) == 10
