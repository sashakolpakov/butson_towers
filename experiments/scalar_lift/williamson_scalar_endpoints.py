"""Exact scalar restrictions of T-designs over rational cyclic group algebras.

This uses Q[C_w], a product of cyclotomic number fields, rather than calling
the weighted source a Butson matrix. The transported trace metric is the
identity in the inverse group basis; its actual used scalar images are
verified sign circulants. Large endpoints have literal entry/row recipes.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from functools import cache
import hashlib
import json
import math
import operator
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
def input_path(relative):
    path = ROOT / relative
    if path.exists():
        return path
    local = Path(__file__).resolve().with_name(Path(relative).name)
    nested = local.parent / "t127" / local.name
    return nested if not local.exists() and nested.exists() else local


INPUTS = input_path("experiments/scalar_lift/williamson_inputs.json")
STATUSES = input_path("experiments/scalar_lift/williamson_order_status.json")
T_INPUTS = {631: input_path("artifacts/hadamard-7572/t631.json"),
            127: input_path("artifacts/t-matrices-127/t127.json")}
# L(x)b: rows are the four coefficient sequences put into the GS array.
# Each entry is (variable index, coefficient sign).
L = (((0,1),(1,-1),(2,-1),(3,-1)),
     ((1,1),(0,1),(3,1),(2,-1)),
     ((2,1),(3,-1),(0,1),(1,1)),
     ((3,1),(2,1),(1,-1),(0,1)))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def array_hash(values):
    return digest(np.asarray(values, dtype=np.int8).tobytes())


def gs(rows):
    p = len(rows[0])
    index = (np.arange(p)[None,:]-np.arange(p)[:,None]) % p
    a,b,c,d = [np.asarray(row)[index] for row in rows]
    r = lambda matrix: matrix[:,::-1]
    return np.block([[a,r(b),r(c),r(d)],
                     [-r(b),a,-r(d.T),r(c.T)],
                     [-r(c),r(d.T),a,-r(b.T)],
                     [-r(d),-r(c.T),r(b.T),a]])


def expand_t(record):
    p,mu,m = record["p"],record["multiplier"],record["multiplier_order"]
    states = record["states"]
    require(p > 1 and m % 3 == 0 and (p-1) % m == 0
            and len(states) == 1+(p-1)//m
            and all(type(s) is int and 0 <= s < 8 for s in states), "invalid T orbit data")
    require(pow(mu,m,p) == 1 and all(pow(mu,d,p) != 1 for d in range(1,m) if m % d == 0),
            "incorrect T multiplier order")
    require(states[0] < 2, "origin channel is not fixed by the channel cycle")
    b = np.zeros((4,p),dtype=np.int64)
    b[0,0] = 1 if states[0] == 0 else -1
    seen,orbit_index = {0},1
    for representative in range(1,p):
        if representative in seen:
            continue
        axis,negative = divmod(states[orbit_index],2)
        point = representative
        for k in range(m):
            require(point not in seen, "overlapping T multiplier orbits")
            seen.add(point)
            channel = 0 if axis == 0 else 1+(axis-1+k)%3
            b[channel,point] = -1 if negative else 1
            point = point*mu % p
        require(point == representative, "T orbit did not close")
        orbit_index += 1
    require(orbit_index == len(states) and len(seen) == p, "incomplete T orbit cover")
    return b


def exact_index(value, bound):
    require(not isinstance(value, (bool,np.bool_)), "Boolean matrix index")
    value = operator.index(value)
    require(0 <= value < bound, "matrix index outside the endpoint")
    return value


def phi(n):
    result, left = n, n
    for p in range(2, math.isqrt(n)+1):
        if left % p == 0:
            result -= result // p
            while left % p == 0:
                left //= p
    if left > 1:
        result -= result // left
    return result


def algebra_description(w):
    divisors = [d for d in range(1,w+1) if w % d == 0]
    require(sum(phi(d) for d in divisors) == w, "cyclotomic degree sum failed")
    return {
        "algebra": f"Q[c]/(c^{w}-1)",
        "algebra_is_single_number_field": w == 1,
        "number_field_factors": [{"conductor":d,"field_degree":phi(d),
                                  "module_rank":1} for d in divisors],
        "rational_dimension": w,
        "coordinate_basis": "beta_j=c^(-j), 0<=j<w",
        "conjugation": "c -> c^(-1)",
        "trace_functional": f"tau(a)=[c^0]a=(1/{w})*sum_(d|{w}) Tr_Q(zeta_d)/Q(a_d)",
        "trace_metric": f"I_{w}",
        "metric_entry_formula": "G_ij=tau(conjugate(beta_i)*beta_j)=delta_ij",
        "multiplication_action": "rho(f)[r,s]=f_coefficients[(s-r) mod w]",
        "adjoint_identity": "rho(conjugate(f))=rho(f)^T",
        "endpoint_frame": f"I_{w}",
        "endpoint_frame_Gram_scale":1,
    }


def group_algebra_description(dimensions):
    dimensions = [operator.index(n) for n in dimensions]
    require(dimensions and all(n > 0 for n in dimensions), "invalid group dimensions")
    order,exponent = math.prod(dimensions),math.lcm(*dimensions)
    divisors = sorted({d for i in range(1,math.isqrt(exponent)+1) if exponent % i == 0
                       for d in (i,exponent//i)})
    counts,factors = {},[]
    for d in divisors:
        count = math.prod(math.gcd(d,n) for n in dimensions)
        count -= sum(counts[e] for e in counts if d % e == 0)
        counts[d] = count
        require(count % phi(d) == 0, "nonintegral cyclotomic field multiplicity")
        if count:
            factors.append({"conductor":d,"field_degree":phi(d),
                            "multiplicity":count//phi(d)})
    require(sum(x["field_degree"]*x["multiplicity"] for x in factors) == order,
            "group-algebra factor dimensions failed")
    return {"algebra":"Q[G]","group":" x ".join(f"C_{n}" for n in dimensions),
            "group_dimensions":dimensions,"rational_dimension":order,
            "group_exponent":exponent,"cyclotomic_field_factors":factors,
            "conjugation":"g -> -g in the additive group",
            "coordinate_basis":"beta_g=c^(-g), indexed by the recorded mixed-radix group coordinates",
            "trace_functional":"tau(a)=[identity]a=(1/|G|)*sum_all_field_factors Tr(a_factor)",
            "trace_metric":f"I_{order}",
            "multiplication_action":"rho(f)[r,s]=f_coefficients[s-r], with componentwise group subtraction",
            "adjoint_identity":"rho(conjugate(f))=rho(f)^T",
            "endpoint_frame_Gram_scale":1}


def verify_t_rows(b):
    b = np.asarray(b)
    require(b.ndim == 2 and b.shape[0] == 4
            and np.issubdtype(b.dtype,np.integer)
            and set(map(int,np.unique(b))) <= {-1,0,1}, "invalid T entry alphabet")
    p = b.shape[1]
    require(np.all(np.sum(abs(b),axis=0) == 1), "T supports do not partition positions")
    correlations = [int(np.sum(b*np.roll(b,-s,axis=1))) for s in range(p)]
    require(correlations == [p]+[0]*(p-1), "complete periodic T Gram failed")
    return {"order":p,"actual_alphabet":[-1,0,1],
            "supports_partition_every_position":True,
            "all_periodic_shifts_checked":p,
            "summed_circulant_Gram_scale":p,
            "generating_rows_int8_sha256":array_hash(b),
            "row_sums":b.sum(axis=1).tolist()}


def coefficient_sequences(b):
    b = np.asarray(b,dtype=np.int64)
    out = np.zeros((4,4,b.shape[1]),dtype=np.int64)
    for i in range(4):
        for j in range(4):
            variable, sign = L[i][j]
            out[variable,i] += sign*b[j]
    return out


def verify_orthogonal_design(b):
    t_proof = verify_t_rows(b)
    p = len(b[0])
    coefficients = coefficient_sequences(b)
    checks = 0
    for a in range(4):
        for z in range(a,4):
            for shift in range(p):
                value = int(np.sum(coefficients[a]*np.roll(coefficients[z],-shift,axis=1)))
                if a != z:
                    value += int(np.sum(coefficients[z]*np.roll(coefficients[a],-shift,axis=1)))
                require(value == (p if a == z and shift == 0 else 0),
                        "orthogonal-design polynomial Gram coefficient failed")
                checks += 1
    matrices = np.array([gs(rows).astype(np.int8) for rows in coefficients])
    require(np.all(np.sum(abs(matrices),axis=0) == 1),
            "orthogonal-design coefficient supports are not a partition")
    for matrix in matrices:
        matrix.setflags(write=False)
    matrices.setflags(write=False)
    proof = {
        "T_source":t_proof,"orthogonal_design_order":4*p,
        "variable_count":4,"variable_weights":[p]*4,
        "source_identity":"D(x)D(x)^T=p*(x0^2+x1^2+x2^2+x3^2)*I_(4p)",
        "construction":"D(x)=GS(L(x)b), with the literal L table retained below",
        "L_variable_and_sign":L,
        "GS_convention":"The standalone gs function in this recipe; reversal has J_uv=1 when u+v=p-1",
        "all_polynomial_periodic_Gram_coefficients_checked":checks,
        "full_polynomial_Gram_verified":True,
        "all_circulant_blocks_commute":True,
        "coefficient_matrices_entry_alphabet":[-1,0,1],
        "coefficient_matrices_disjoint_and_cover_every_entry":True,
        "coefficient_matrices_int8_sha256":[array_hash(x) for x in matrices],
        "full_dense_coefficient_matrix_products_performed":False,
    }
    return matrices, proof


@cache
def verified_t_design(p):
    path = T_INPUTS[p]
    record = json.loads(path.read_text())
    b = expand_t(record)
    if "b_sequences" in record:
        require(np.array_equal(b,np.asarray(record["b_sequences"])), "stored T rows disagree")
    matrices, proof = verify_orthogonal_design(b)
    proof.update({"literal_T_source":str(path.relative_to(ROOT)),
                  "literal_T_source_sha256":digest(path.read_bytes()),
                  "orbit_states":record["states"],
                  "multiplier":record["multiplier"],
                  "multiplier_order":record["multiplier_order"]})
    return matrices, proof


def verify_williamson(rows):
    rows = np.asarray(rows)
    require(rows.ndim == 2 and rows.shape[0] == 4
            and np.issubdtype(rows.dtype,np.integer)
            and np.all(abs(rows) == 1), "Williamson rows must have four sign words")
    rows = rows.astype(np.int64)
    w = rows.shape[1]
    indices = (np.arange(w)[None,:]-np.arange(w)[:,None]) % w
    matrices = rows[:,indices]
    require(np.array_equal(matrices,matrices.transpose(0,2,1)), "Williamson polynomial is not self-adjoint")
    for i in range(4):
        for j in range(i):
            require(np.array_equal(matrices[i]@matrices[j].T,
                                   matrices[j]@matrices[i].T), "Williamson amicability failed")
    gram = sum(a@a.T for a in matrices)
    require(np.array_equal(gram,4*w*np.eye(w,dtype=np.int64)), "Williamson summed Gram failed")
    # This is the exact source relation in Q[C_w], independent of matrix Gram code.
    convolution = np.zeros(w,dtype=np.int64)
    for row in rows:
        for i in range(w):
            convolution += int(row[i])*np.roll(row,i)
    require(convolution.tolist() == [4*w]+[0]*(w-1), "cyclic-algebra norm relation failed")
    matrices = matrices.astype(np.int8)
    matrices.setflags(write=False)
    return matrices, {
        "order":w,"first_rows_pm":["".join("+" if x == 1 else "-" for x in row) for row in rows],
        "actual_used_scalar_polynomials":"f_j(c)=sum_k first_rows[j][k]*c^k",
        "all_f_j_self_adjoint":True,"all_pairwise_amicability_identities_verified":True,
        "full_integer_Gram_verified":True,"cyclic_algebra_norm_relation_verified":True,
        "exact_norm_relation":f"sum_j f_j(c)^2={4*w} in Q[c]/(c^{w}-1)",
        "row_sums":rows.sum(axis=1).tolist(),
        "first_rows_int8_sha256":array_hash(rows),
        "entry_images_are_literal_sign_circulants":True,
    }


@cache
def verified_williamson(w):
    record = json.loads(INPUTS.read_text())
    entry = next((entry for entry in record["entries"] if entry["order"] == w),None)
    require(entry is not None, "no retained literal Williamson source at this order")
    require(len(entry["first_rows_pm"]) == 4
            and all(isinstance(word,str) and len(word) == w and set(word) <= {"+","-"}
                    for word in entry["first_rows_pm"]), "invalid literal Williamson sign word")
    rows = np.array([[1 if x == "+" else -1 for x in word]
                     for word in entry["first_rows_pm"]],dtype=np.int64)
    require(rows.shape == (4,w), "wrong length of literal Williamson source")
    matrices, proof = verify_williamson(rows)
    require(proof["first_rows_int8_sha256"] == entry["first_rows_int8_row_major_sha256"],
            "Williamson source hash mismatch")
    proof.update({"input_file":str(INPUTS.relative_to(ROOT)),
                  "input_file_sha256":digest(INPUTS.read_bytes()),
                  "primary_source":record["source_reference"]})
    return matrices, proof


@dataclass
class ScalarEndpoint:
    X: np.ndarray = field(repr=False)
    A: np.ndarray = field(repr=False)

    @property
    def order(self):
        return self.X.shape[1]*self.A.shape[1]

    def entry(self,row,column):
        row,column = exact_index(row,self.order),exact_index(column,self.order)
        i,a = divmod(row,self.A.shape[1]); j,b = divmod(column,self.A.shape[1])
        return sum(int(self.X[k,i,j])*int(self.A[k,a,b]) for k in range(4))

    def row(self,row):
        row = exact_index(row,self.order)
        i,a = divmod(row,self.A.shape[1])
        return sum(self.X[k,i,:,None]*self.A[k,a,None,:] for k in range(4)).reshape(-1)


@dataclass
class ImplicitScalarEndpoint:
    X: np.ndarray = field(repr=False)
    source: object = field(repr=False)

    @property
    def order(self):
        return self.X.shape[1]*self.source.order

    def entry(self,row,column):
        row,column = exact_index(row,self.order),exact_index(column,self.order)
        i,a = divmod(row,self.source.order); j,b = divmod(column,self.source.order)
        k = int(np.flatnonzero(self.X[:,i,j])[0])
        return int(self.X[k,i,j])*int(self.source.entry(k,a,b))

    def row(self,row):
        row = exact_index(row,self.order)
        i,a = divmod(row,self.source.order)
        signs = self.X[:,i].sum(axis=0,dtype=np.int8)
        variables = np.argmax(abs(self.X[:,i]),axis=0)
        source_rows = np.array([self.source.row(k,a) for k in range(4)],dtype=np.int8)
        return (signs[:,None]*source_rows[variables]).reshape(-1)


def endpoint(p,w):
    x,_ = verified_t_design(p)
    a,_ = verified_williamson(w)
    return ScalarEndpoint(x,a)


def row_controls(matrix):
    w = matrix.A.shape[1] if isinstance(matrix,ScalarEndpoint) else matrix.source.order
    selected = sorted({0,1,w-1,w,matrix.order//2,matrix.order-1})
    packed = []
    for i in selected:
        row = matrix.row(i)
        require(row.shape == (matrix.order,) and np.all(abs(row) == 1), "flat endpoint entry test failed")
        for j in selected:
            require(matrix.entry(i,j) == int(row[j]), "entry/row interfaces disagree")
        packed.append(int.from_bytes(np.packbits(row<0,bitorder="little").tobytes(),"little"))
    for i in range(len(selected)):
        for j in range(i):
            require((packed[i]^packed[j]).bit_count() == matrix.order//2,
                    "literal endpoint row inner product failed")
    return {"row_indices":selected,"all_selected_full_rows_are_signs":True,
            "all_selected_row_pairs_checked_by_exact_bit_counts":True,
            "complete_Gram_proof":"Full polynomial OD identity and complete exact Williamson norm relation; row checks are interface controls."}


def certificate(p=631,orders=None):
    record = json.loads(INPUTS.read_text())
    if orders is None:
        orders = [item["order"] for item in record["entries"]]
    status = json.loads(STATUSES.read_text()) if STATUSES.exists() else {}
    # Permit a descriptive envelope, with actual per-endpoint records under entries.
    status = status.get("entries",status)
    x,t_proof = verified_t_design(p)
    outputs = []
    for w in orders:
        a,source = verified_williamson(w)
        literal = ScalarEndpoint(x,a)
        known = status.get(f"{p}:{w}")
        used_polynomials = sorted({
            "".join("+" if sign*coefficient == 1 else "-" for coefficient in a[j,0])
            for j in range(4) for sign in set(map(int,np.unique(x[j]))) if sign
        })
        output = {"p":p,"w":w,"output_order":literal.order,
            "Williamson_source":source,"scalar_algebra":algebra_description(w),
            "weighted_source_order":4*p,"weighted_source_Gram_scale":4*p*w,
            "weighted_source_is_Butson":False,
            "weighted_source_used_alphabet":"{+f_j(c),-f_j(c):0<=j<4}, with literal polynomials recorded above",
            "actual_used_source_polynomials_pm":used_polynomials,
            "actual_used_source_polynomial_count":len(used_polynomials),
            "source_identity":f"D(f)D(f)*={p}*sum_j f_j(c)^2*I={4*p*w}I",
            "rational_dimension_after_restriction":4*p*w,
            "endpoint_Gram_scale":4*p*w,"endpoint_entry_alphabet":[-1,1],
            "endpoint_recipe":"H=sum_(j=0)^3 X_j tensor rho(f_j); X_j are the four disjoint GS orthogonal-design coefficients",
            "complete_algebraic_endpoint_Gram_verified":True,
            "all_used_scalar_images_verified_sign_circulants":True,
            "literal_entry_and_row_recipe_available":True,
            "dense_endpoint_allocated":False,
            "complete_dense_endpoint_Gram_performed":False,
            "literal_row_controls":row_controls(literal),
            "reported_order_status":known}
        outputs.append(output)
    return {"schema":"complete-cyclotomic-algebra-scalar-endpoints-v1",
        "T_design":t_proof,"outputs":outputs,
        "all_outputs_are_completed_literal_sign_matrix_recipes":True,
        "restriction_scope":"Restriction of scalars from the finite etale algebra Q[C_w], explicitly factored into cyclotomic number fields; weighted source entries are not called Butson roots.",
        "implementation_sha256":digest(Path(__file__).read_bytes()),
        "reconstruction_dependencies":["Python>=3.10","NumPy","neighboring literal T and Williamson input JSON files"],
        "Williamson_input_sha256":digest(INPUTS.read_bytes()),
        "status_record":str(STATUSES.relative_to(ROOT)) if STATUSES.exists() else None,
        "status_record_sha256":digest(STATUSES.read_bytes()) if STATUSES.exists() else None}


def family_certificate():
    try:
        from experiments.scalar_lift import williamson_spence_family as family
    except ModuleNotFoundError:
        import williamson_spence_family as family
    status_path = input_path("experiments/scalar_lift/williamson_family_order_status.json")
    statuses = json.loads(status_path.read_text()) if status_path.exists() else {}
    statuses = statuses.get("entries",statuses)
    x,t_proof = verified_t_design(631)
    source_proof = family.certificate()
    outputs = []
    for source_record in source_proof["first_expanded_members"]:
        q,r = source_record["q"],source_record["r"]
        source = family.SpenceSource(q,r)
        literal = ImplicitScalarEndpoint(x,source)
        rows = source.coefficient_rows().astype(np.int8)
        actual_polynomials = sorted({
            "".join("+" if sign*a == 1 else "-" for a in row)
            for row in rows for sign in (-1,1)})
        outputs.append({
            "q":q,"r":r,"source_order":source.order,
            "output_order":literal.order,
            "scalar_algebra":group_algebra_description(source.group_dimensions),
            "weighted_source_order":2524,
            "weighted_source_Gram_scale":literal.order,
            "weighted_source_is_Butson":False,
            "actual_used_source_polynomials_pm":actual_polynomials,
            "polynomial_coordinate_order":"mixed-radix coordinates on the recorded finite abelian group",
            "all_source_group_norm_coefficients_verified":True,
            "source_norm_certificate":source_record,
            "endpoint_Gram_scale":literal.order,
            "actual_endpoint_entry_alphabet":[-1,1],
            "complete_algebraic_endpoint_Gram_verified":True,
            "literal_endpoint_recipe":"sum_j X_j tensor A_j; source entries and rows are the exact SpenceSource recursion",
            "literal_entry_and_row_recipe_available":True,
            "literal_row_controls":row_controls(literal),
            "dense_endpoint_allocated":False,
            "complete_dense_endpoint_Gram_performed":False,
            "reported_order_status":statuses.get(f"family:{q}:{r}")})
    return {"schema":"complete-spence-cyclotomic-algebra-endpoints-v1",
            "T_design":t_proof,"family_source_proof":source_proof,
            "outputs":outputs,
            "all_outputs_are_completed_literal_sign_matrix_recipes":True,
            "implementation_sha256":digest(Path(__file__).read_bytes()),
            "source_recurrence_implementation_sha256":digest(Path(family.__file__).read_bytes()),
            "status_record_sha256":digest(status_path.read_bytes()) if status_path.exists() else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p",type=int,choices=sorted(T_INPUTS),default=631)
    parser.add_argument("--orders",type=int,nargs="+")
    parser.add_argument("--families",action="store_true",help="verify the five infinite Spence source families and their first T631 scalar endpoints")
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    require(not args.families or (args.p == 631 and args.orders is None),
            "--families uses the T631 design and its five recorded first expanded sources")
    report = family_certificate() if args.families else certificate(args.p,args.orders)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"orders":[x["output_order"] for x in report["outputs"]],
                      "complete_Gram_proofs":True,"report":str(args.output)}))


if __name__ == "__main__":
    main()
