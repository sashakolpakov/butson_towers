"""Pure scalar restrictions through coprime period/cyclotomic field towers.

No right frame, Gram rescaling, or new source construction occurs here.
Every source entry and trace metric has a lazy Kronecker-factor definition.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import cache

import numpy as np
from sympy import Matrix

from experiments.scalar_lift.relative_scalar_restriction import edge, restriction_in_basis


@cache
def _local(lower_phase):
    if lower_phase == 1:
        return restriction_in_basis(12, 1, (0, 2, 1, 3))
    return edge(12, lower_phase)


def _json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass(frozen=True)
class PeriodRestrictionStage:
    field: object
    source: object
    lower_phase: int
    period_restricted: bool
    factor_order: str

    def __post_init__(self):
        if self.source.phase not in (6, 12) or self.lower_phase not in (1, 6, 12):
            raise ValueError('source and lower fields are outside the certified tower')
        if self.factor_order not in ('period-local', 'local-period'):
            raise ValueError('a literal factor order is required')

    @property
    def local(self):
        return _local(self.lower_phase)

    @property
    def width(self):
        return self.local.width*(self.field.degree if self.period_restricted else 1)

    @property
    def order(self):
        return self.source.order*self.width

    @property
    def coefficient_field_degree(self):
        return self.local.lower_degree*(1 if self.period_restricted else self.field.degree)

    def _coordinates(self, index):
        if not self.period_restricted:
            return 0, index
        if self.factor_order == 'period-local':
            return divmod(index, self.local.width)
        local, period = divmod(index, self.field.degree)
        return period, local

    def scalar_entry(self, source_exponent, row, column):
        """Coefficients in the displayed lower cyclotomic subfield.

        If L remains in the coefficient field, these coefficients are
        embedded as rational constants in L, where 1=-sum eta_j.
        """
        if not 0 <= row < self.width or not 0 <= column < self.width:
            raise ValueError('scalar indices outside the relative degree')
        pr, lr = self._coordinates(row)
        pc, lc = self._coordinates(column)
        action = self.local.actions.numerators[(12//self.source.phase)*int(source_exponent) % 12]
        numbers = action[lr, lc]*(pr == pc)
        return {'numerators': numbers.tolist(), 'denominator': self.local.actions.denominator,
                'coefficient_root_label': self.lower_phase,
                'period_coefficient': 'rational constant' if not self.period_restricted else None}

    def entry(self, row, column):
        if int(row) != row or int(column) != column or not 0 <= row < self.order or not 0 <= column < self.order:
            raise ValueError('matrix indices outside the stage order')
        i, a = divmod(int(row), self.width)
        j, b = divmod(int(column), self.width)
        return self.scalar_entry(self.source.entry(i, j), a, b)

    def scalar_table(self, source_exponent):
        local = self.local.action((12//self.source.phase)*int(source_exponent)).numerators
        if not self.period_restricted:
            return local
        identity = np.eye(self.field.degree, dtype=np.int64)
        return np.stack([np.kron(identity, local[:, :, k]) if self.factor_order == 'period-local'
                         else np.kron(local[:, :, k], identity)
                         for k in range(self.local.lower_degree)], axis=2)

    def metric_table(self):
        local = self.local.metric.numerators
        if not self.period_restricted:
            return local
        return np.stack([np.kron(self.field.trace_gram, local[:, :, k])
                         if self.factor_order == 'period-local'
                         else np.kron(local[:, :, k], self.field.trace_gram)
                         for k in range(self.local.lower_degree)], axis=2)

    def literal_hash(self):
        digest = hashlib.sha256()
        header = {'source_hash': self.source.source_hash, 'order': self.order,
                  'coefficient_field': [self.lower_phase,
                       None if self.period_restricted else self.field.conductor,
                       None if self.period_restricted else self.field.degree],
                  'coefficient_field_degree': self.coefficient_field_degree,
                  'source_phase': self.source.phase,
                  'action_denominator': self.local.actions.denominator,
                  'metric_denominator': self.local.metric.denominator}
        digest.update(json.dumps(header, sort_keys=True).encode())
        for exponent in self.source.used_alphabet:
            digest.update(np.asarray(self.scalar_table(exponent), dtype='<i8').tobytes())
        digest.update(np.asarray(self.metric_table(), dtype='<i8').tobytes())
        return digest.hexdigest()

    def metadata(self):
        field = self.field
        local = self.local
        return {'order': self.order, 'relative_width': self.width,
                'source_hash': self.source.source_hash, 'source_order': self.source.order,
                'coefficient_field': 'Q' if self.lower_phase == 1 and self.period_restricted else
                    (f'Q(zeta_{self.lower_phase})' if self.period_restricted else
                     ('L' if self.lower_phase == 1 else f'L(zeta_{self.lower_phase})')),
                'coefficient_field_degree': self.coefficient_field_degree,
                'period_restricted': self.period_restricted, 'literal_factor_order': self.factor_order,
                'literal_action_metric_sha256': self.literal_hash(),
                'local_action_table': local.actions.to_json(),
                'local_trace_metric': local.metric.to_json(),
                'local_inverse_trace_metric': local.inverse_metric.to_json(),
                'period_trace_metric_factor': ('ell I-h J' if self.period_restricted else 'coefficient field'),
                'period_inverse_trace_metric_factor': ('(I+h J)/ell' if self.period_restricted else None),
                'local_source_exponent_multiplier': 12//self.source.phase,
                'period_multiplication_sha256': field.proof()['multiplication_sha256'],
                'coefficient_embedding_when_L_remains': 'rational 1 = -sum eta_j',
                'weighted_Gram_identity': 'R (I_n tensor G^-1) R* = n (I_n tensor G^-1)',
                'Gram_scale': self.source.order, 'source_Gram_replayed': self.source.proof['full_Gram_verified'],
                'ordinary_Hadamard_Gram_or_flatness_claimed': False,
                'right_frame_inserted': False,
                'exact_local_restriction_definition_sha256': local.metadata()['exact_definition_sha256']}


def _verify_local_composition():
    """Independently compose the two small field edges, including traces."""
    first, second, combined = edge(12, 6), edge(6, 1), _local(1)
    for exponent in range(12):
        result = np.zeros((4, 4), dtype=np.int64)
        for i in range(2):
            for j in range(2):
                for b in range(2):
                    result[2*i:2*i+2, 2*j:2*j+2] += (
                        first.actions.numerators[exponent, i, j, b]
                        *second.actions.numerators[b, :, :, 0])
        if not np.array_equal(result, combined.actions.numerators[exponent, :, :, 0]):
            raise ValueError('the actual successive restriction actions differ')
    # Tr_12/Q = Tr_6/Q Tr_12/6 in the product power basis.
    if not np.array_equal(combined.metric.numerators[:, :, 0],
                          np.kron(2*np.eye(2, dtype=np.int64), second.metric.numerators[:, :, 0])):
        raise ValueError('the successive trace metric differs')
    return True


def period_towers(field, source):
    _verify_local_composition()
    first = [PeriodRestrictionStage(field, source, q, True, 'period-local') for q in (12, 6, 1)]
    second = [PeriodRestrictionStage(field, source, q, False, 'local-period') for q in (6, 1)]
    second.append(PeriodRestrictionStage(field, source, 1, True, 'local-period'))
    # A endpoint coordinates (period,local); B coordinates (local,period).
    r = field.degree
    permutation = np.array([j*4+a for a in range(4) for j in range(r)], dtype=np.int64)
    for exponent in source.used_alphabet:
        left = first[-1].scalar_table(exponent)[permutation][:, permutation]
        if not np.array_equal(left, second[-1].scalar_table(exponent)):
            raise ValueError('endpoint actions do not obey their exact perfect shuffle')
    if not np.array_equal(first[-1].metric_table()[permutation][:, permutation], second[-1].metric_table()):
        raise ValueError('endpoint trace metrics do not obey the exact perfect shuffle')
    # The auxiliary multiplication undergoes the very same permutation.
    if not np.array_equal(np.kron(field.multiplication, np.eye(4, dtype=np.int64))[permutation][:, permutation],
                          np.kron(np.eye(4, dtype=np.int64), field.multiplication)):
        raise ValueError('the auxiliary faithful action was dropped in the tower')
    return first, second


def tower_report(field, source):
    paths = period_towers(field, source)
    records = [[stage.metadata() for stage in path] for path in paths]
    initial = PeriodRestrictionStage(field, source, 12, False, 'local-period').metadata()
    hashes = [initial['literal_action_metric_sha256']] + [stage['literal_action_metric_sha256'] for path in records for stage in path]
    r = field.degree
    # A direct small trace-Gram proof fixes orientation of the inverse metric.
    G4 = Matrix(_local(1).metric.numerators[:, :, 0].tolist())
    actions = _local(1).actions.numerators[:, :, :, 0]
    small = Matrix.vstack(*(Matrix.hstack(*(Matrix(actions[(12//source.phase)*source.entry(i, j)].tolist())
                                            for j in range(source.order))) for i in range(source.order)))
    metric_inv = Matrix.diag(*([G4.inv()]*source.order))
    if small*metric_inv*small.T != source.order*metric_inv:
        raise ValueError('the complete source weighted rational Gram failed')
    for stage in (paths[0][-1], paths[1][-1]):
        off_period_column = stage.local.width if stage.factor_order == 'period-local' else 1
        if any(stage.scalar_entry(source.used_alphabet[0], 0, off_period_column)['numerators']):
            raise ValueError('the product-basis zero witness failed')
    return {'scope': 'Pure scalar restriction only; no right frame or Gram normalization',
            'source_id': source.id, 'source_hash': source.source_hash,
            'ambient_degree': 4*r, 'ambient_conductor_descriptor': [12, field.conductor],
            'initial_source': initial, 'paths': records, 'edge_relative_degrees': [[r, 2, 2], [2, 2, r]],
            'literal_state_count': len(hashes), 'literal_distinct_state_count': len(set(hashes)),
            'endpoint_coordinate_permutation': 'B(a,j) = A(j,a), index_A=4*j+a',
            'all_used_scalar_actions_and_metrics_shuffle_verified': True,
            'faithful_auxiliary_multiplication_shuffle_verified': True,
            'complete_source_weighted_Gram_verified': True,
            'source_Gram_scale_unchanged': source.order,
            'local_composition_all_twelve_roots_verified': True,
            'no_flat_endpoint_claimed': True,
            'product_basis_sign_endpoint_excluded': True,
            'zero_reason': 'Every source scalar is its local action tensor I_r, up to shuffle; r>1 leaves zero off-period blocks',
            'period_only_trace_normalization_cannot_change_source_action': True,
            'nonproduct_joint_bases_or_unequal_input_output_frames_excluded': False,
            'definition_sha256': _json_hash(records)}
