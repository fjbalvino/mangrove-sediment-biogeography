"""Independently reconcile long-form evidence, matrices and summary counts."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def verify(directory):
    if not __debug__:
        raise RuntimeError('Verification requires Python without -O')
    directory = Path(directory)
    diag = json.loads((directory / 'diagnostico_y_validacion.json').read_text())
    long = pd.read_csv(directory / 'evidencia_por_regla_muestra.tsv.gz', sep='\t')
    summary = pd.read_csv(directory / 'resumen_reglas.tsv', sep='\t')
    scope = pd.read_csv(directory / 'alcance_65_procesos.tsv', sep='\t').fillna('')
    meta = pd.read_csv(directory / 'muestras_diseno_y_deteccion.tsv', sep='\t')
    ids = meta.run_accession.tolist()
    n = len(ids)
    assert n == diag['n_muestras'] and len(set(ids)) == n
    assert len(long) == 2 * n * (diag['n_reglas'] + diag['n_modulos_referencia'])
    assert not long.duplicated(['rule_id', 'run_accession', 'umbral_registros']).any()
    assert set(long.resultado_logico) <= {'TRUE', 'FALSE', 'UNKNOWN'}
    matrices = {}
    comparisons = 0
    for threshold in (1, 2):
        part = long[long.umbral_registros == threshold]
        wide = part.pivot(index='rule_id', columns='run_accession', values='resultado_logico')
        matrix = pd.read_csv(directory / f'matriz_97_reglas_ge{threshold}.tsv', sep='\t', index_col=0)
        matrices[threshold] = matrix
        assert list(matrix.columns) == ids
        assert len(matrix) == diag['n_reglas'] and matrix.index.is_unique
        expected = wide.loc[matrix.index, ids].apply(lambda s: s.map({'TRUE': 1., 'FALSE': 0., 'UNKNOWN': np.nan}))
        np.testing.assert_allclose(matrix.to_numpy(), expected.to_numpy(), equal_nan=True)
        counts = part.groupby(['rule_id', 'resultado_logico']).size().unstack(fill_value=0)
        for state, column in [('TRUE', 'n_cumple'), ('FALSE', 'n_no_cumple'), ('UNKNOWN', 'n_indeterminado')]:
            actual = summary[summary.umbral_registros == threshold].set_index('rule_id')[column]
            expected_counts = counts.reindex(columns=['TRUE', 'FALSE', 'UNKNOWN'], fill_value=0)[state]
            assert actual.index.is_unique
            pd.testing.assert_series_equal(actual.sort_index(), expected_counts.sort_index(), check_names=False)
        assert (counts.sum(axis=1) == n).all()
        representative = pd.read_csv(directory / f'matriz_65_reglas_representativas_ge{threshold}.tsv', sep='\t', index_col=0)
        assert list(representative.columns) == ids
        assert len(representative) == diag['n_procesos'] and representative.index.is_unique
        for row in scope.itertuples():
            if row.tipo == 'no_resuelto':
                assert (representative.loc[row.process_id] == 'requiere_otra_evidencia').all()
            else:
                observed = part[part.rule_id == row.regla_representativa].set_index('run_accession').loc[ids, 'estado']
                assert representative.loc[row.process_id].tolist() == observed.tolist()
            comparisons += 1
    assert matrices[1].index.equals(matrices[2].index)
    assert not ((matrices[2] == 1) & (matrices[1] != 1)).to_numpy().any()
    report = {'status': 'PASS', 'samples': n, 'rules': diag['n_reglas'],
              'reference_modules': diag['n_modulos_referencia'],
              'representative_row_comparisons': comparisons,
              'long_rows': len(long), 'scope': 'Internal output consistency; not biological validation'}
    (directory / 'verificacion_salidas.json').write_text(json.dumps(report, indent=2) + '\n')
    print('BIOGEOCHEMICAL_OUTPUTS_PASS ' + json.dumps(report))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True)
    try:
        verify(parser.parse_args().results)
    except Exception:
        print('BIOGEOCHEMICAL_OUTPUTS_FAIL')
        raise
