"""Partial charge preservation is required before interpreting adsorption."""
import pytest
from cof_vision.cofdraw_science import import_cif, normalize_structure, structure_exports


def structure():
    return {'periodic': True, 'cell': [[10,0,0],[0,10,0],[0,0,10]], 'atoms': [
        {'element':'C', 'fractional':[.1,.2,.3], 'charge':.35},
        {'element':'N', 'fractional':[.24,.2,.3], 'charge':-.35}]}


def test_partial_charges_roundtrip_cif_and_repeated_view():
    result = structure_exports(structure(), repeats=[2,1,1])
    assert [a['charge'] for a in result['viewer']['atoms']] == [.35,-.35,.35,-.35]
    restored = import_cif(result['exports']['cif'])
    assert [a['charge'] for a in restored['atoms']] == [.35,-.35,.35,-.35]
    assert sum(a['charge'] for a in restored['atoms']) == pytest.approx(0)


def test_unknown_partial_charge_stays_unknown_not_zero():
    source = structure(); source['atoms'][1].pop('charge')
    restored = import_cif(structure_exports(source)['exports']['cif'])
    assert restored['atoms'][0]['charge'] == .35
    assert 'charge' not in restored['atoms'][1]
    source['atoms'][0]['charge'] = float('nan')
    with pytest.raises(ValueError, match='partial charge'):
        normalize_structure(source)
