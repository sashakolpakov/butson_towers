"""Bounded exact SAT lifting of unrestricted periodic T-sequence signs.

Unsigned axes first pass the SL2 parity parametrization. The companion model
builder eliminates the modulo-four sign equations, without a machine-word
dimension limit. This backend encodes the resulting modulo-eight quadratics,
the four exact negative-sign counts, and all integer correlations using native
XOR constraints and Boolean binary adders. Every positive witness is expanded
independently at every shift. SAT/UNSAT applies only to the fixed axes and the
declared row sums and gauges; a bounded UNKNOWN excludes no construction.

The generic XOR and adder primitives are retained from the checked earlier
fixed-multiplier solver. This module has its own four-coordinate model and
does not enumerate an 80- or 86-dimensional affine space.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from math import gcd, isfinite
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


def _bits(value):
    while value:
        low = value & -value
        yield low.bit_length()-1
        value ^= low


class XorCNF:
    """Tseitin clauses, native XORs, and exact nonnegative weighted sums."""

    def __init__(self):
        self.variables = 0
        self.clauses = []
        self.xors = []
        self._and_cache = {}
        self._xor_cache = {}
        self._majority_cache = {}
        self.true = self.var()
        self.clause([self.true])

    def var(self):
        self.variables += 1
        return self.variables

    def clause(self, literals):
        literals = set(map(int, literals))
        if any(-v in literals for v in literals):
            return
        self.clauses.append(tuple(sorted(literals, key=lambda z: (abs(z), z))))

    def xor(self, literals, rhs):
        rhs = int(rhs)
        variables = set()
        for literal in literals:
            literal = int(literal)
            if literal < 0:
                rhs ^= 1
            v = abs(literal)
            if v in variables:
                variables.remove(v)
            else:
                variables.add(v)
        if not variables:
            if rhs:
                self.clause([])
            return
        if len(variables) == 1:
            v = next(iter(variables))
            self.clause([v if rhs else -v])
            return
        self.xors.append((tuple(sorted(variables)), rhs))

    def and_literal(self, a, b):
        if a == b:
            return a
        if a == -b or a == -self.true or b == -self.true:
            return -self.true
        if a == self.true:
            return b
        if b == self.true:
            return a
        key = tuple(sorted((a, b)))
        if key not in self._and_cache:
            c = self.var()
            self.clause([-c, a])
            self.clause([-c, b])
            self.clause([c, -a, -b])
            self._and_cache[key] = c
        return self._and_cache[key]

    def xor_literal(self, literals):
        variables, constant = set(), 0
        for literal in literals:
            if literal == self.true:
                constant ^= 1
                continue
            if literal == -self.true:
                continue
            if literal < 0:
                constant ^= 1
            v = abs(literal)
            if v in variables:
                variables.remove(v)
            else:
                variables.add(v)
        if not variables:
            return self.true if constant else -self.true
        key = tuple(sorted(variables))
        if len(key) == 1:
            value = key[0]
        else:
            if key not in self._xor_cache:
                value = self.var()
                self.xor([value, *key], 0)
                self._xor_cache[key] = value
            value = self._xor_cache[key]
        return -value if constant else value

    def majority_literal(self, a, b, c):
        if a == b or a == c:
            return a
        if b == c:
            return b
        if a == -b:
            return c
        if a == -c:
            return b
        if b == -c:
            return a
        values = (a, b, c)
        if self.true in values:
            x, y = [v for v in values if v != self.true]
            return -self.and_literal(-x, -y)
        if -self.true in values:
            x, y = [v for v in values if v != -self.true]
            return self.and_literal(x, y)
        key = tuple(sorted(values))
        if key not in self._majority_cache:
            out = self.var()
            for x, y in itertools.combinations(values, 2):
                self.clause([-x, -y, out])
                self.clause([x, y, -out])
            self._majority_cache[key] = out
        return self._majority_cache[key]

    def weighted_exact(self, terms, target):
        """Encode sum(weight * [literal]) == target by exact binary adders."""
        target = int(target)
        coefficients = Counter()
        for literal, weight in terms:
            literal, weight = int(literal), int(weight)
            if weight < 0:
                raise ValueError('negative weight')
            if literal == self.true:
                target -= weight
            elif literal != -self.true:
                coefficients[literal] += weight
        for variable in {abs(v) for v in coefficients}:
            constant = min(coefficients[variable], coefficients[-variable])
            target -= constant
            coefficients[variable] -= constant
            coefficients[-variable] -= constant
        terms = [(v, w) for v, w in coefficients.items() if w]
        total = sum(w for _, w in terms)
        if target < 0 or target > total:
            self.clause([])
            return
        if not terms:
            return
        divisor = 0
        for _, weight in terms:
            divisor = gcd(divisor, weight)
        if target % divisor:
            self.clause([])
            return
        terms = [(v, w//divisor) for v, w in terms]
        target //= divisor
        total //= divisor
        if target in (0, total):
            for literal, _ in terms:
                self.clause([literal if target else -literal])
            return
        columns = []
        for literal, weight in terms:
            for bit in range(weight.bit_length()):
                if (weight >> bit) & 1:
                    while len(columns) <= bit:
                        columns.append([])
                    columns[bit].append(literal)
        bit = 0
        while bit < len(columns) or (target >> bit):
            while len(columns) <= bit+1:
                columns.append([])
            column = columns[bit]
            while len(column) >= 3:
                a, b, c = column.pop(), column.pop(), column.pop()
                column.append(self.xor_literal([a, b, c]))
                columns[bit+1].append(self.majority_literal(a, b, c))
            if len(column) == 2:
                a, b = column.pop(), column.pop()
                column.append(self.xor_literal([a, b]))
                columns[bit+1].append(self.and_literal(a, b))
            expected = (target >> bit) & 1
            if column:
                self.clause([column[0] if expected else -column[0]])
            elif expected:
                self.clause([])
            bit += 1
            # Do not keep appending an empty carry column forever.
            if bit == len(columns)-1 and not columns[bit] and not (target >> bit):
                break

    def dimacs(self):
        lines = [f'p cnf {self.variables} {len(self.clauses)+len(self.xors)}']
        lines += [' '.join(map(str, row))+' 0' for row in self.clauses]
        for variables, rhs in self.xors:
            row = list(variables)
            if rhs == 0:
                row[0] = -row[0]
            lines.append('x'+' '.join(map(str, row))+' 0')
        return '\n'.join(lines)+'\n'

    def solve(self, seconds=5, seed=0, workdir=None):
        executable = shutil.which('cryptominisat5')
        if executable is None:
            raise RuntimeError('cryptominisat5 is required')
        if not isfinite(seconds) or seconds <= 0:
            raise ValueError('positive time limit required')
        if workdir is None:
            with tempfile.TemporaryDirectory(prefix='butson-xor-sat-') as directory:
                return self.solve(seconds, seed, directory)
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        path = workdir/'problem.cnf'
        path.write_text(self.dimacs())
        start = time.monotonic()
        command = [executable, '--verb', '0', '--threads', '1', '--random', str(int(seed) % 2147483647),
                   '--maxtime', str(seconds), str(path)]
        try:
            process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, timeout=seconds+1)
            output = process.stdout
        except subprocess.TimeoutExpired as exc:
            output = exc.stdout or ''
            if isinstance(output, bytes):
                output = output.decode(errors='replace')
            process = None
        (workdir/'solver.log').write_text(output)
        if re.search(r'^s SATISFIABLE\s*$', output, re.MULTILINE):
            status = 'SAT'
        elif re.search(r'^s UNSATISFIABLE\s*$', output, re.MULTILINE):
            status = 'UNSAT'
        elif process is None or 'INDETERMINATE' in output or 'UNKNOWN' in output:
            status = 'UNKNOWN'
        else:
            raise RuntimeError('unexpected solver output: '+output[:200])
        values = {}
        if status == 'SAT':
            for line in output.splitlines():
                if line.startswith('v '):
                    for word in line[2:].split():
                        literal = int(word)
                        if literal:
                            values[abs(literal)] = int(literal > 0)
        return dict(status=status, values=values, wall_seconds=time.monotonic()-start,
                    variables=self.variables, clauses=len(self.clauses), xors=len(self.xors),
                    process_wall_limit_reached=process is None,
                    problem_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def add_correlation(cnf, system, signs, index):
    """Encode C(s)=0 as the exact weighted cut sum (A-C(s))/2=A/2."""
    absolute = int(system['diagonal'][index])
    weights = system['weights'][index]
    terms = []
    for i in range(len(signs)):
        for j in range(i+1, len(signs)):
            weight = int(weights[i, j])
            if weight:
                if weight < 0:
                    raise ValueError('unsigned pair multiplicities must be nonnegative')
                absolute += weight
                # Use the already linked site signs, keeping every pair XOR
                # ternary instead of repeating long affine free-variable XORs.
                cut = cnf.xor_literal([signs[i], signs[j]])
                terms.append((cut, weight))
    if absolute % 2:
        raise ValueError('unsigned correlation parity failed')
    cnf.weighted_exact(terms, absolute//2)


def encode_sign_system(system, exact_correlations=True, mod8=True):
    """Return CNF plus free/site literals for a prepared four-count model.

    A site literal is its negative-sign bit x_i. The affine equations give
    x=p+By over F2. For a pair with unsigned multiplicity w, its signed
    contribution is w*(1-2*(x_i XOR x_j)). Consequently C(s)=0 is equivalent
    to an exact nonnegative weighted sum of pair XORs. The adder output bits
    are fixed to that integer target, with carries retained at every bit.
    """
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    cnf = XorCNF()
    free = [cnf.var() for _ in system['basis']]
    forms = system['forms']
    part = system['particular']

    def affine(form, constant=0):
        literal = cnf.xor_literal([free[i] for i in _bits(form)])
        return -literal if constant else literal

    signs = [affine(form, (part >> i) & 1) for i, form in enumerate(forms)]
    if len(system['groups']) != 4 or len(system['negative_counts']) != 4:
        raise ValueError('an unrestricted sign model requires four coordinate counts')
    for indices, count in zip(system['groups'], system['negative_counts'], strict=True):
        cnf.weighted_exact([(signs[i], 1) for i in indices], count)
    if mod8:
        for equation in lift.quadratic_mod8(system):
            literals = [free[u] for u in _bits(equation['linear'])]
            literals += [cnf.and_literal(free[u], free[v]) for u, v in equation['quadratic']]
            cnf.xor(literals, equation['constant'])
    if exact_correlations:
        for index in range(len(system['diagonal'])):
            add_correlation(cnf, system, signs, index)
    return cnf, free, signs


def solve(system, seconds=5, seed=0, workdir=None, exact_correlations=True,
          mod8=True, all_solutions=False, max_solutions=1):
    """Search a fixed unrestricted model within a bounded wall-time budget.

    By default all exact correlations are encoded immediately. The optional
    lazy mode starts from the count-constrained modular model and adds every
    violated exact correlation after each candidate. Each such cut is a
    necessary integer equation; it cannot discard a true solution. Returned
    SAT always means that a fully verified T-sequence was found, in both modes.
    All-solution enumeration blocks complete free-variable assignments, never
    auxiliary assignments. Only a final UNSAT establishes completeness.
    """
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    if not isfinite(seconds) or seconds <= 0 or max_solutions < 0:
        raise ValueError('positive finite seconds and a nonnegative solution limit are required')
    if not exact_correlations and not mod8:
        raise ValueError('lazy exact correlations require the modulo-eight model')
    if workdir is None:
        with tempfile.TemporaryDirectory(prefix='unrestricted-t-sign-sat-') as directory:
            return solve(system, seconds, seed, directory, exact_correlations,
                         mod8, all_solutions, max_solutions)
    started = time.monotonic()
    deadline = started+seconds
    cnf, free, signs = encode_sign_system(system, exact_correlations, mod8)
    build_seconds = time.monotonic()-started
    encoded = set(range(len(system['diagonal']))) if exact_correlations else set()
    solutions, attempts = [], []
    candidates = rejected = 0
    complete = False
    last_status = 'UNKNOWN'
    while time.monotonic() < deadline:
        remaining = deadline-time.monotonic()
        if remaining <= 0:
            break
        result = cnf.solve(remaining, int(seed)+len(attempts), workdir)
        model = result.pop('values')
        attempts.append(result)
        last_status = result['status']
        if last_status == 'UNSAT':
            complete = True
            break
        if last_status != 'SAT':
            break
        candidates += 1
        if not all(v in model for v in free):
            raise RuntimeError('solver omitted a free-variable value')
        assignment = sum(model[v] << i for i, v in enumerate(free))
        negative = lift.sign_mask(system, assignment)
        if any(sum((negative >> i) & 1 for i in group) != count
               for group, count in zip(system['groups'], system['negative_counts'], strict=True)):
            raise ArithmeticError('SAT encoding admitted an incorrect exact coordinate count')
        residuals = lift.correlations(system, negative)
        if any(residuals):
            if exact_correlations:
                raise ArithmeticError('SAT encoding admitted a nonzero exact correlation')
            violated = [i for i, value in enumerate(residuals) if value]
            if any(index in encoded for index in violated):
                raise ArithmeticError('SAT encoding violated a previously imposed exact correlation')
            for index in violated:
                add_correlation(cnf, system, signs, index)
                encoded.add(index)
            rejected += 1
            continue
        states = [2*a+((negative >> i) & 1) for i, a in enumerate(system['axes'])]
        proof = lift.verify_states(states, system['data'], system['target_row'])
        solutions.append(dict(states=states, assignment=assignment,
                              expanded_verification=proof))
        if not all_solutions or (max_solutions and len(solutions) >= max_solutions):
            break
        cnf.clause([-variable if model[variable] else variable for variable in free])
    return dict(
        status='SAT' if solutions else ('UNSAT' if complete else 'UNKNOWN'),
        solver_status=last_status, complete_for_fixed_axes=complete,
        affine_dimension=len(free), exact_negative_counts=list(system['negative_counts']),
        target_row=list(system['target_row']), solutions=solutions,
        exact_correlations_encoded_initially=bool(exact_correlations),
        exact_correlations_encoded_finally=len(encoded), modulo_eight_encoded=bool(mod8),
        candidates_checked=candidates, nonzero_correlations_rejected=rejected,
        build_seconds=build_seconds, elapsed_seconds=time.monotonic()-started,
        variables=cnf.variables, clauses=len(cnf.clauses), xors=len(cnf.xors),
        attempts=attempts,
        scope='Fixed coordinate pattern, target row sums and explicitly supplied sign gauges only.')


def solve_axes(axes, p_or_data, target_row, *, fixed_signs=None, **options):
    """Prepare and solve unrestricted axes; expose arithmetic infeasibility."""
    from experiments.scalar_lift import unrestricted_t_sign_model as lift
    started = time.monotonic()
    try:
        system = lift.prepare(axes, p_or_data, target_row, fixed_signs=fixed_signs)
    except lift.InfeasibleAxes as exc:
        return dict(status='UNSAT', solver_status='INFEASIBLE_ARITHMETIC',
                    reason=str(exc), complete_for_fixed_axes=True, solutions=[],
                    elapsed_seconds=time.monotonic()-started,
                    scope='Fixed coordinate pattern, target row sums and explicitly supplied sign gauges only.')
    return solve(system, **options)


def solve_rows(axes, data, target_row, seconds=120, seed=0, max_solutions=1,
               workdir=None, **options):
    """Campaign entrypoint with the common axes/data/row-sums interface."""
    return solve_axes(axes, data, target_row, seconds=seconds, seed=seed,
                      max_solutions=max_solutions, workdir=workdir, **options)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axes', type=Path, required=True)
    parser.add_argument('--p', type=int, required=True)
    parser.add_argument('--target-row', type=int, nargs=4, required=True)
    parser.add_argument('--seconds', type=float, default=60)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--workdir', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--lazy-correlations', action='store_true')
    parser.add_argument('--all-solutions', action='store_true')
    parser.add_argument('--max-solutions', type=int, default=1)
    args = parser.parse_args()
    record = json.loads(args.axes.read_text())
    axes = record['axes'] if isinstance(record, dict) else record
    result = solve_axes(axes, args.p, args.target_row, seconds=args.seconds,
                        seed=args.seed, workdir=args.workdir,
                        exact_correlations=not args.lazy_correlations,
                        all_solutions=args.all_solutions, max_solutions=args.max_solutions)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],
                          complete_for_fixed_axes=result['complete_for_fixed_axes'],
                          solutions=len(result['solutions']),
                          elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    main()
