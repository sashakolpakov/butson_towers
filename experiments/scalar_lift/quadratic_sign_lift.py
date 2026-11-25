"""Exact fixed-coordinate sign lifting without a machine-word size limit.

The input specifies which of the four coordinates is occupied on every
multiplication orbit.  The remaining signs are binary variables x_i.
The imported exact mod-four elimination gives x=p+By over F_2.  We add
signs forced by extreme exact counts, derive the quadratic mod-eight
correlation equations on y, and solve them together with the complete
integer correlations and the two exact sign counts.

All positive outputs are independently expanded at every residue and
checked at every shift.  A solver timeout is recorded as unfinished.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import shutil
import subprocess
import tempfile
import json
import operator
from pathlib import Path
import time

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as frame


class InfeasibleAxes(ValueError):
    """The specified coordinate pattern fails a necessary exact condition."""


def _bits(value):
    while value:
        low = value & -value
        yield low.bit_length() - 1
        value ^= low


def _eliminate(rows, variables):
    """Solve packed binary rows, with the RHS in bit `variables`."""
    rows = list(rows)
    pivots = []
    for column in range(variables):
        found = next((i for i in range(len(pivots), len(rows))
                      if (rows[i] >> column) & 1), None)
        if found is None:
            continue
        at = len(pivots)
        rows[at], rows[found] = rows[found], rows[at]
        for i in range(len(rows)):
            if i != at and (rows[i] >> column) & 1:
                rows[i] ^= rows[at]
        pivots.append(column)
    mask = (1 << variables) - 1
    if any((row & mask) == 0 and (row >> variables) & 1 for row in rows):
        raise InfeasibleAxes("forced exact-count signs contradict the mod-four system")
    particular = sum(((rows[i] >> variables) & 1) << col
                     for i, col in enumerate(pivots))
    basis = []
    for free in range(variables):
        if free not in pivots:
            basis.append((1 << free) | sum(((rows[i] >> free) & 1) << col
                                           for i, col in enumerate(pivots)))
    return particular, basis


def sign_mask(system, assignment):
    """Recover the orbit-negative mask from a free-variable assignment."""
    if not 0 <= assignment < (1 << len(system['basis'])):
        raise ValueError("free-variable assignment out of range")
    result = system['particular']
    for i in _bits(assignment):
        result ^= system['basis'][i]
    return result


def prepare(axes, data, head_net=1, tail_net=1, first_negative=0,
            force_exact_count_extremes=True, fixed_signs=None):
    """Return a reduced affine sign system, exact counts and graph weights.

    Nonzero coordinate-zero orbits have total signed count `head_net`;
    all other orbits together have total signed count `tail_net`.
    Consequently the expanded coordinate sums are
    (1+m*head_net, (m//3)*tail_net, (m//3)*tail_net, (m//3)*tail_net).
    Orbit zero is +e0.  Orbit one must be coordinate zero, and its sign
    is fixed by `first_negative`; this is an explicitly stated gauge.
    """
    axes = tuple(axes)
    if any(isinstance(x, (bool, np.bool_)) for x in axes):
        raise ValueError("coordinate labels must be integers rather than Booleans")
    axes = tuple(operator.index(x) for x in axes)
    head_net, tail_net = operator.index(head_net), operator.index(tail_net)
    n = data['n']
    row_sums = [1+data['m']*head_net] + [data['m']//3*tail_net]*3
    if sum(value*value for value in row_sums) != data['p']:
        raise InfeasibleAxes("requested coordinate sums fail the exact norm condition")
    if len(axes) != n or any(not 0 <= x < 4 for x in axes) or axes[0] != 0:
        raise ValueError("one coordinate label in0..3 is required per orbit, with zero first")
    if axes[1] != 0:
        raise ValueError("the first nonzero orbit must be a head for this gauge")
    if first_negative not in (0, 1):
        raise ValueError("first_negative must be0 or1")
    heads = [i for i in range(1, n) if axes[i] == 0]
    tails = [i for i in range(1, n) if axes[i] != 0]
    wanted = []
    forced = {0: 0, 1: int(first_negative)}
    for index, bit in (fixed_signs or {}).items():
        index, bit = operator.index(index), operator.index(bit)
        if not 0 <= index < n or bit not in (0, 1):
            raise ValueError("fixed signs require in-range indices and binary values")
        if index in forced and forced[index] != bit:
            raise InfeasibleAxes("an additional sign contradicts the fixed gauge")
        forced[index] = bit
    for group, net in ((heads, head_net), (tails, tail_net)):
        if abs(net) > len(group) or (len(group) - net) % 2:
            raise InfeasibleAxes("coordinate count cannot have the requested signed sum")
        target = (len(group) - net) // 2
        wanted.append(target)
        known_negative = sum(forced[index] for index in group if index in forced)
        unknown = [index for index in group if index not in forced]
        needed = target-known_negative
        if not 0 <= needed <= len(unknown):
            raise InfeasibleAxes("fixed signs contradict an exact count")
        if force_exact_count_extremes and needed in (0, len(unknown)):
            bit = int(needed != 0)
            for index in unknown:
                forced[index] = bit
    try:
        system = frame.affine_sign_system(axes, data, head_net, tail_net, first_negative)
    except ValueError as exc:
        raise InfeasibleAxes(str(exc)) from exc
    old_particular, old_basis = system['particular'], system['basis']
    old_dimension = len(old_basis)
    extra = []
    for index, bit in forced.items():
        coefficients = sum(((vector >> index) & 1) << j
                           for j, vector in enumerate(old_basis))
        rhs = bit ^ ((old_particular >> index) & 1)
        extra.append(coefficients | (rhs << old_dimension))
    free_particular, free_basis = _eliminate(extra, old_dimension)
    particular = old_particular
    for j in _bits(free_particular):
        particular ^= old_basis[j]
    basis = []
    for vector in free_basis:
        image = 0
        for j in _bits(vector):
            image ^= old_basis[j]
        basis.append(image)
    system.update(particular=particular, basis=basis, axes=axes, data=data,
                  heads=heads, tails=tails, negative_counts=wanted,
                  head_net=int(head_net), tail_net=int(tail_net),
                  first_negative=int(first_negative), forced_signs=forced,
                  affine_dimension_before_extreme_counts=old_dimension)
    system['forms'] = [sum(((vector >> i) & 1) << j
                           for j, vector in enumerate(basis)) for i in range(n)]
    return system



def iter_head_sign_branches(axes, head_net=-5, first_negative=1):
    """Yield each complete head-sign assignment with the exact signed count.

    Orbit one is the fixed head.  For h=5,7,9 and head_net=-5,
    first_negative=1 this yields respectively1,6,28 branches.
    Pass each returned dictionary as `fixed_signs` to `prepare` or `solve_bits`.
    """
    heads = [i for i, axis in enumerate(axes) if i and axis == 0]
    if not heads or heads[0] != 1 or first_negative not in (0, 1):
        raise ValueError("head branching requires the stated first-head gauge")
    if abs(head_net) > len(heads) or (len(heads)-head_net) % 2:
        raise InfeasibleAxes("head count cannot realize the requested signed sum")
    target = (len(heads)-head_net)//2
    others = heads[1:]
    needed = target-first_negative
    if not 0 <= needed <= len(others):
        raise InfeasibleAxes("first-head gauge contradicts the exact head count")
    # Enumerating the smaller sign class avoids large intermediate tuples.
    selected_bit = int(needed <= len(others)-needed)
    selected_count = needed if selected_bit else len(others)-needed
    for selected in itertools.combinations(others, selected_count):
        chosen = set(selected)
        result = {1: int(first_negative)}
        result.update({i: selected_bit if i in chosen else 1-selected_bit for i in others})
        yield result


def correlations(system, negative_mask):
    """Evaluate the exact compressed integer correlations."""
    signs = np.array([1 - 2*((negative_mask >> i) & 1)
                      for i in range(len(system['axes']))], dtype=np.int64)
    return [int(diag + (signs @ weight @ signs)//2)
            for diag, weight in zip(system['diagonal'], system['weights'], strict=True)]


def quadratic_mod8(system):
    """Return all equations C_ell(p+By)/4=0 modulo two.

    Write s_i s_j = e_ij (1-2 L_ij(y)), where e_ij is +1 or-1 and
    L_ij is a binary linear form.  Modulo four,
    L_ij = sum a_u y_u - 2 sum_(u<v) a_u a_v y_u y_v.
    Each correlation is therefore quadratic modulo eight.  The affine
    mod-four equations make its constant and linear coefficients divisible
    by four, so the returned Boolean quadratic is exactly C_ell/4 mod2.

    `linear` is a packed bit mask; `quadratic` lists pairs (u,v) with u<v.
    """
    dimension = len(system['basis'])
    equations = []
    part, forms = system['particular'], system['forms']
    for diagonal, weight in zip(system['diagonal'], system['weights'], strict=True):
        constant = int(diagonal)
        linear_sums = [0]*dimension
        quadratic_rows = [0]*dimension
        for i in range(len(forms)):
            for j in range(i+1, len(forms)):
                w = int(weight[i, j])
                if not w:
                    continue
                signed_weight = w * (1-2*(((part >> i) ^ (part >> j)) & 1))
                constant += signed_weight
                form = forms[i] ^ forms[j]
                for u in _bits(form):
                    linear_sums[u] += signed_weight
                    if w & 1:
                        quadratic_rows[u] ^= form & ((1 << u)-1)
        if constant % 4 or any(value % 2 for value in linear_sums):
            raise ArithmeticError("affine mod-four reduction did not make division by four exact")
        equations.append({
            'constant': (constant//4) & 1,
            'linear': sum(((-value//2) & 1) << u for u, value in enumerate(linear_sums)),
            'quadratic': [(v, u) for u, row in enumerate(quadratic_rows) for v in _bits(row)],
        })
    return equations


def evaluate_quadratic(equation, assignment):
    value = equation['constant'] ^ ((equation['linear'] & assignment).bit_count() & 1)
    for u, v in equation['quadratic']:
        value ^= ((assignment >> u) & (assignment >> v)) & 1
    return value


def verify_states(states, data, head_net=None, tail_net=None):
    """Independently check every expanded periodic correlation and the sums."""
    states = list(map(int, states))
    b = frame.expand(states, data)
    values = [int(np.sum(b*np.roll(b, -shift, axis=1), dtype=np.int64))
              for shift in range(data['p'])]
    if values != [data['p']] + [0]*(data['p']-1):
        raise ValueError("candidate failed an expanded exact periodic correlation")
    sums = list(map(int, b.sum(axis=1)))
    if head_net is not None and tail_net is not None:
        target = [1+data['m']*head_net] + [data['m']//3*tail_net]*3
        if sums != target:
            raise ValueError("candidate failed the requested exact coordinate sums")
    return {'every_expanded_periodic_correlation_verified': True,
            'shift_count': data['p'], 'sequence_row_sums': sums,
            'sequence_support_sizes': list(map(int, np.count_nonzero(b, axis=1)))}


def _xor_equal(model, literals, rhs):
    if not literals:
        if rhs:
            model.add(False)
    elif rhs:
        model.add_bool_xor(literals)
    else:
        model.add_bool_xor([*literals, True])


def build_model(system, add_mod8=True, affine_links=True, exact_correlations=True):
    """Build early exact counts, MQs, and optional complete integer cuts."""
    from ortools.sat.python import cp_model
    model = cp_model.CpModel()
    n, dimension = len(system['axes']), len(system['basis'])
    signs = [model.new_bool_var(f'x{i}') for i in range(n)]
    free = [model.new_bool_var(f'y{j}') for j in range(dimension)] if affine_links else []
    if affine_links:
        for i, form in enumerate(system['forms']):
            _xor_equal(model, [signs[i]]+[free[u] for u in _bits(form)],
                       (system['particular'] >> i) & 1)
    else:
        for i, bit in system['forced_signs'].items():
            model.add(signs[i] == bit)
    for group, target in zip((system['heads'], system['tails']), system['negative_counts'], strict=True):
        model.add(sum(signs[i] for i in group) == target)
    pair_xors = {}
    if exact_correlations:
        for diagonal, weight in zip(system['diagonal'], system['weights'], strict=True):
            terms, absolute = [], int(diagonal)
            for i in range(n):
                for j in range(i+1, n):
                    w = int(weight[i, j])
                    if w:
                        absolute += w
                        if (i, j) not in pair_xors:
                            cut = model.new_bool_var(f'cut{i}_{j}')
                            _xor_equal(model, [signs[i], signs[j], cut], 0)
                            pair_xors[i, j] = cut
                        terms.append(w*pair_xors[i, j])
            if absolute % 2:
                raise InfeasibleAxes("unsigned correlation has odd parity")
            model.add(sum(terms) == absolute//2)
    equations = quadratic_mod8(system) if add_mod8 and affine_links else []
    products = {}
    for equation in equations:
        literals = [free[u] for u in _bits(equation['linear'])]
        for u, v in equation['quadratic']:
            if (u, v) not in products:
                z = model.new_bool_var(f'product{u}_{v}')
                model.add_multiplication_equality(z, [free[u], free[v]])
                products[u, v] = z
            literals.append(products[u, v])
        _xor_equal(model, literals, equation['constant'])
    error = model.validate()
    if error:
        raise ValueError(error)
    return model, signs, free, {
        'model_variables': len(model.proto.variables),
        'model_constraints': len(model.proto.constraints),
        'affine_dimension': dimension,
        'affine_dimension_before_extreme_counts': system['affine_dimension_before_extreme_counts'],
        'forced_sign_count': len(system['forced_signs']),
        'quadratic_mod8_equations': len(equations),
        'cut_variables': len(pair_xors), 'quadratic_product_variables': len(products),
        'exact_negative_counts': system['negative_counts'],
        'integer_correlation_constraints_in_model': bool(exact_correlations),
    }


def solve(axes, data, head_net=1, tail_net=1, first_negative=0, seconds=10,
          seed=0, all_solutions=False, max_solutions=1, add_mod8=True,
          affine_links=True, force_exact_count_extremes=True, exact_correlations=True,
          fixed_signs=None):
    """Run a bounded exact solver and verify every positive output independently."""
    from ortools.sat.python import cp_model
    if not exact_correlations and not (add_mod8 and affine_links):
        raise ValueError("MQ-only solving requires affine links and mod-eight equations")
    if seconds <= 0 or max_solutions < 0:
        raise ValueError("positive time limit and nonnegative solution limit are required")
    started = time.monotonic()
    try:
        system = prepare(axes, data, head_net, tail_net, first_negative,
                         force_exact_count_extremes, fixed_signs=fixed_signs)
    except InfeasibleAxes as exc:
        return {'solver_status': 'INFEASIBLE_ARITHMETIC', 'reason': str(exc),
                'complete_for_fixed_axes': True, 'solutions': [],
                'elapsed_seconds': time.monotonic()-started}
    model, signs, free, metadata = build_model(system, add_mod8, affine_links, exact_correlations)
    build_time = time.monotonic()-started
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = int(seed)
    solver.parameters.enumerate_all_solutions = bool(all_solutions or not exact_correlations)
    solutions = []
    limited = False
    relaxed_candidates = rejected_candidates = 0

    class Collector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            nonlocal limited, relaxed_candidates, rejected_candidates
            relaxed_candidates += 1
            negative_mask = sum(int(self.value(x)) << i for i, x in enumerate(signs))
            if any(correlations(system, negative_mask)):
                rejected_candidates += 1
                return
            states = [2*a+int(self.value(x)) for a, x in zip(system['axes'], signs, strict=True)]
            proof = verify_states(states, data, head_net, tail_net)
            solutions.append({'states': states, 'expanded_verification': proof})
            if not all_solutions or (max_solutions and len(solutions) >= max_solutions):
                limited = True
                self.stop_search()

    status = solver.solve(model, Collector())
    name = solver.status_name(status)
    return dict(metadata, solver_status=name,
                complete_for_fixed_axes=(not limited and status in (cp_model.OPTIMAL, cp_model.INFEASIBLE)),
                solutions=solutions, build_seconds=build_time,
                relaxed_candidates_checked=relaxed_candidates,
                nonzero_integer_correlations_rejected=rejected_candidates,
                solver_seconds=solver.wall_time, branches=solver.num_branches,
                conflicts=solver.num_conflicts, elapsed_seconds=time.monotonic()-started,
                scope='Fixed coordinate pattern and stated exact row-sum/gauge constraints only.')



def build_bits_binary(workdir):
    """Compile the small streaming Boolean-MQ kernel with its source hash."""
    source = Path(__file__).with_name('quadratic_sign_lift_bits.cpp')
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    identifier = hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    target = workdir/('quadratic-sign-bits-'+identifier)
    if not target.exists():
        compiler = shutil.which('c++')
        if compiler is None:
            raise RuntimeError('a C++17 compiler is required for the bit-sliced solver')
        subprocess.run([compiler, '-O3', '-std=c++17', str(source), '-o', str(target)],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return target


def solve_bits(axes, data, head_net=1, tail_net=1, first_negative=0,
               seconds=10, max_roots=100000, workdir=None, fixed_signs=None):
    """Enumerate MQ roots in native64-lane words, then test exact conditions.

    The forced extreme counts are included in the affine reduction. The
    remaining exact counts and all integer correlations are checked on
    each returned root, followed by an independent expanded verification.
    No positive result relies only on the modular equations.
    """
    started = time.monotonic()
    try:
        system = prepare(axes, data, head_net, tail_net, first_negative,
                         fixed_signs=fixed_signs)
    except InfeasibleAxes as exc:
        return {'solver_status': 'INFEASIBLE_ARITHMETIC', 'reason': str(exc),
                'complete_for_fixed_axes': True, 'solutions': [],
                'elapsed_seconds': time.monotonic()-started}
    dimension = len(system['basis'])
    if dimension > 60:
        raise ValueError('the streaming native kernel supports at most60 free variables')
    equations = quadratic_mod8(system)
    with tempfile.TemporaryDirectory(prefix='quadratic-sign-bits-') as temporary:
        directory = Path(temporary) if workdir is None else Path(workdir)
        directory.mkdir(parents=True, exist_ok=True)
        binary = build_bits_binary(directory)
        input_path, output_path = directory/'mq-input.txt', directory/'mq-result.json'
        lines = [f'{dimension} {len(equations)}']
        for equation in equations:
            rows = [0]*dimension
            for u, v in equation['quadratic']:
                rows[v] |= 1 << u
            lines += [f"{equation['constant']} {equation['linear']}",
                      ' '.join(map(str, rows))]
        input_path.write_text('\n'.join(lines)+'\n')
        subprocess.run([str(binary), str(input_path), str(output_path),
                        str(seconds), str(max_roots)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        native = json.loads(output_path.read_text())
    solutions = []
    passed_counts = 0
    for assignment in native['roots']:
        if any(evaluate_quadratic(e, assignment) for e in equations):
            raise ArithmeticError('native root failed independent quadratic evaluation')
        mask = sign_mask(system, assignment)
        if any(sum((mask >> i) & 1 for i in group) != target
               for group, target in zip((system['heads'], system['tails']),
                                        system['negative_counts'], strict=True)):
            continue
        passed_counts += 1
        if any(correlations(system, mask)):
            continue
        states = [2*a+((mask >> i) & 1) for i, a in enumerate(system['axes'])]
        proof = verify_states(states, data, head_net, tail_net)
        solutions.append({'states': states, 'expanded_verification': proof})
    return {
        'solver_status': 'COMPLETE' if native['complete'] else 'UNFINISHED',
        'complete_for_fixed_axes': native['complete'],
        'affine_dimension': dimension,
        'affine_dimension_before_extreme_counts': system['affine_dimension_before_extreme_counts'],
        'forced_sign_count': len(system['forced_signs']),
        'quadratic_mod8_equations': len(equations),
        'exact_negative_counts': system['negative_counts'],
        'mod8_roots': len(native['roots']), 'exact_count_survivors': passed_counts,
        'solutions': solutions, 'native': native,
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'Fixed coordinate pattern and stated exact row-sum/gauge constraints only.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axes', type=Path, required=True, help='JSON list or object with axes')
    parser.add_argument('--p', type=int, required=True)
    parser.add_argument('--m', type=int, required=True)
    parser.add_argument('--mu', type=int)
    parser.add_argument('--head-net', type=int, default=1)
    parser.add_argument('--tail-net', type=int, default=1)
    parser.add_argument('--first-negative', type=int, choices=(0, 1), default=0)
    parser.add_argument('--seconds', type=float, default=10)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--all-solutions', action='store_true')
    parser.add_argument('--max-solutions', type=int, default=1)
    parser.add_argument('--no-mod8', action='store_true')
    parser.add_argument('--no-affine-links', action='store_true')
    parser.add_argument('--bits', action='store_true', help='streaming bit-sliced MQ enumeration')
    parser.add_argument('--native-workdir', type=Path)
    parser.add_argument('--mq-only', action='store_true',
                        help='enumerate count-constrained MQ candidates and test exact correlations lazily')
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    record = json.loads(args.axes.read_text())
    axes = record['axes'] if isinstance(record, dict) else record
    data = frame.chart(args.p, args.m, args.mu)
    if args.bits:
        report = solve_bits(axes, data, args.head_net, args.tail_net, args.first_negative,
                            args.seconds, workdir=args.native_workdir)
    else:
        report = solve(axes, data, args.head_net, args.tail_net, args.first_negative,
                       args.seconds, args.seed, args.all_solutions, args.max_solutions,
                       not args.no_mod8, not args.no_affine_links,
                       exact_correlations=not args.mq_only)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('solver_status', 'elapsed_seconds', 'complete_for_fixed_axes')} |
                     {'solutions': len(report['solutions'])}))


if __name__ == '__main__':
    main()
