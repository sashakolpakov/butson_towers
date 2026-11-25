"""CryptoMiniSat encoding of exact sign lifts, with native XOR constraints.

The finite Boolean encodings are checked separately on exhaustive controls.
A positive matrix witness is always expanded and checked over every shift.
A solver timeout says only that this particular search was unfinished.
"""
from __future__ import annotations

from collections import Counter
import itertools
import json
from math import gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


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
        if seconds <= 0:
            raise ValueError('positive time limit required')
        if workdir is None:
            with tempfile.TemporaryDirectory(prefix='butson-xor-sat-') as directory:
                return self.solve(seconds, seed, directory)
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        path = workdir/'problem.cnf'
        path.write_text(self.dimacs())
        start = time.monotonic()
        command = [executable, '--verb', '0', '--threads', '1', '--random', str(seed),
                   '--maxtime', str(seconds), str(path)]
        try:
            process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, timeout=seconds+10)
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
                    variables=self.variables, clauses=len(self.clauses), xors=len(self.xors))


def encode_sign_system(system, exact_correlations=True, mod8=True):
    from experiments.scalar_lift import quadratic_sign_lift as lift
    cnf = XorCNF()
    free = [cnf.var() for _ in system['basis']]
    forms = system['forms']
    part = system['particular']

    def affine(form, constant=0):
        literal = cnf.xor_literal([free[i] for i in lift._bits(form)])
        return -literal if constant else literal

    signs = [affine(form, (part >> i) & 1) for i, form in enumerate(forms)]
    for indices, count in zip((system['heads'], system['tails']), system['negative_counts'], strict=True):
        cnf.weighted_exact([(signs[i], 1) for i in indices], count)
    if mod8:
        for equation in lift.quadratic_mod8(system):
            literals = [free[u] for u in lift._bits(equation['linear'])]
            literals += [cnf.and_literal(free[u], free[v]) for u, v in equation['quadratic']]
            cnf.xor(literals, equation['constant'])
    if exact_correlations:
        for diagonal, weights in zip(system['diagonal'], system['weights'], strict=True):
            absolute = int(diagonal)
            terms = []
            for i in range(len(forms)):
                for j in range(i+1, len(forms)):
                    weight = int(weights[i, j])
                    if weight:
                        absolute += weight
                        cut = affine(forms[i] ^ forms[j], ((part >> i) ^ (part >> j)) & 1)
                        terms.append((cut, weight))
            if absolute % 2:
                raise ValueError('unsigned correlation parity failed')
            cnf.weighted_exact(terms, absolute//2)
    return cnf, free, signs


def solve(system, seconds=5, seed=0, workdir=None, exact_correlations=True, mod8=True):
    from experiments.scalar_lift import quadratic_sign_lift as lift
    cnf, free, _ = encode_sign_system(system, exact_correlations, mod8)
    result = cnf.solve(seconds, seed, workdir)
    model = result.pop('values')
    result.update(affine_dimension=len(free), exact_correlations_encoded=exact_correlations)
    if result['status'] == 'SAT':
        if not all(v in model for v in free):
            raise RuntimeError('solver omitted a free-variable value')
        assignment = sum(model[v] << i for i, v in enumerate(free))
        negative = lift.sign_mask(system, assignment)
        residuals = lift.correlations(system, negative)
        states = [2*a+((negative >> i) & 1) for i, a in enumerate(system['axes'])]
        result.update(assignment=assignment, states=states, correlations=residuals,
                      exact_zero=not any(residuals))
        if exact_correlations and any(residuals):
            raise ArithmeticError('SAT encoding admitted a nonzero exact correlation')
        if not any(residuals):
            result['expanded_verification'] = lift.verify_states(states, system['data'], system['head_net'], system['tail_net'])
    return result
