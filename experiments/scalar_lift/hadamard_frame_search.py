"""Active exact frame search: SL2 parity axes, finite mod4/mod8 sign lifts.

See docs/HADAMARD_FRAME_SEARCH.md. The stochastic sampler currently targets
p=631; supplied maps can add other validated small-field charts. Finite axis-pool lifting supports prime twisted charts with at
most 63 orbit variables. No ignored research module or generated table is
required. All native candidates receive complete integer endpoint checks.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
NATIVE_SOURCE = Path(__file__).with_suffix(".cpp")
TABLE631_SHA256 = 'b6c7507dc7efb2d1cedb39f92828dced72fbd69dbb6183157c61c216a3227729'
D4 = np.array(((0,1,1,1),(-1,0,1,-1),(-1,-1,0,1),(-1,1,-1,0)), dtype=np.int64)
# Exact inverse Fourier/color basis. Its complete quadratic identity is checked
# by validate_parity_map(), rather than trusted from a numerical eigensystem.
PARITY_INVERSE_COLUMNS_631 = (
    425404997632973157, 932445608469002174, 137406980501082596,
    16354118061632560, 1067712475950855374, 374351069947066572,
    852149619653713806, 8784224337604941, 834659803112503595,
    889565314010800593, 334094183686342775, 499650198083033278,
    343854616418592859, 633950041397048960, 526687327096792332,
    754821233956843194, 310989893490661847, 179087894977585784,
    8776703039435792, 532140989096581191, 725743236514946116,
    964213289079989703, 10443268505080454, 1002763388014730813,
    336828038264500834, 644947619654339737, 825608000286928087,
    597543140425153709, 883856819647208384, 686022321374657028,
    997440284354041446, 621012106088539574, 1146553316779962455,
    231337793886600323, 40891743351930134, 44545606742272530,
    765548807163580212, 989943777445285922, 892657040088779597,
    879730336583006681, 282892711115313836, 122992105876703925,
    105231803364164325, 449357232292679477, 1052347590734575115,
    462155496618175419, 929862604395215579, 395947605640155305,
    79499731132077249, 63766604014704171, 65164032196656931,
    1145405152020841752, 492111612335104051, 338244334427540870,
    295750111248457326, 105287807578252276, 205011627573678042,
    187582493781427066, 798421410602501402, 571905232869284621,
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def binary_linear(columns, bits):
    value = 0
    while bits:
        bit = bits & -bits
        value ^= columns[bit.bit_length()-1]
        bits ^= bit
    return value


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            bit = value.bit_length()-1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def field_multiply(a, b, degree, modulus):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        if a & (1<<degree):
            a ^= modulus
        b >>= 1
    return value


def chart(p, m, mu=None):
    require(p > 2 and all(p % d for d in range(2, int(p**.5)+1)), "p must be prime")
    require(m > 0 and m % 3 == 0 and (p-1) % m == 0, "invalid multiplier order")
    if mu is None:
        factors = [d for d in range(2,p) if (p-1)%d == 0
                   and all(d%e for e in range(2,int(d**.5)+1))]
        primitive = next(g for g in range(2,p)
                         if all(pow(g,(p-1)//q,p) != 1 for q in factors))
        mu = pow(primitive,(p-1)//m,p)
    require(pow(mu,m,p) == 1 and all(pow(mu,d,p) != 1 for d in range(1,m) if m%d == 0),
            "mu does not have the specified order")
    index = np.zeros(p,dtype=np.int64)
    phase = np.zeros(p,dtype=np.int64)
    seen, orbits = {0}, [[(0,0)]]
    for x in range(1,p):
        if x in seen:
            continue
        orbit, y = [], x
        for k in range(m):
            require(y not in seen, "overlapping orbit")
            seen.add(y)
            index[y], phase[y] = len(orbits), k%3
            orbit.append((y,k%3))
            y = y*mu%p
        require(y == x, "orbit failed to close")
        orbits.append(orbit)
    shifts, covered = [], {0}
    for s in range(1,p):
        if s in covered:
            continue
        shifts.append(s)
        v = s
        for _ in range(m):
            covered.update((v,-v%p))
            v = v*mu%p
    rows = []
    for s in shifts:
        row = defaultdict(int)
        for x in range(p):
            y = (x+s)%p
            i,j,delta = int(index[x]),int(index[y]),int(phase[y]-phase[x])%3
            if i > j:
                i,j,delta = j,i,-delta%3
            if i == 0:
                delta = 0
            row[i,j,delta] += 1
        require(sum(row.values()) == p, "coefficient count failed")
        rows.append(dict(row))
    return dict(p=p,m=m,mu=mu,n=len(orbits),index=index,phase=phase,
                orbits=orbits,shifts=shifts,rows=rows)


def expand(states, data):
    require(len(states) == data['n'] and all(type(s) is int and 0<=s<8 for s in states)
            and states[0] == 0, "invalid exact integer states")
    states = np.asarray(states,dtype=np.int64)[data['index']]
    axis = states//2
    coordinate = np.where(axis == 0,0,1+(axis-1+data['phase'])%3)
    b = np.zeros((4,data['p']),dtype=np.int64)
    b[coordinate,np.arange(data['p'])] = 1-2*(states%2)
    return b


def write_tables(path, data, head_net=1, tail_net=1):
    n,L = data['n'],len(data['rows'])
    require(n <= 63, "native finite lift requires at most63 orbit variables")
    H = np.zeros((n,8,L),dtype=np.int16)
    K = np.zeros((n,n,8,8,L),dtype=np.int16)
    axis,sign = np.arange(8)//2,1-2*(np.arange(8)%2)
    for ell,row in enumerate(data['rows']):
        for (i,j,delta),count in row.items():
            rotated = np.where(axis == 0,0,1+(axis-1+delta)%3)
            if i == j:
                H[i,:,ell] += count*(axis == rotated)
            else:
                values = count*(axis[:,None] == rotated[None,:])*sign[:,None]*sign[None,:]
                K[i,j,:,:,ell] += values
                K[j,i,:,:,ell] += values.T
    with Path(path).open('wb') as output:
        np.array([data['p'],data['m'],data['mu'],n,L,1+data['m']*head_net,
                  data['m']//3*tail_net],dtype=np.int32).tofile(output)
        H.tofile(output)
        K.tofile(output)
    if (data['p'],data['m'],data['mu'],head_net,tail_net) == (631,21,458,1,1):
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == TABLE631_SHA256,
                "generated table differs from successful exact631 table")


def write_pool(path, data, patterns, head_net=1, tail_net=1, first_negative=0):
    require(all(len(a) == data['n'] and a[:2] == [0,0]
                and all(type(v) is int and 0<=v<4 for v in a) for a in patterns),
            "axis pools require exact integer colors0..3 and origin/first-head0")
    with Path(path).open('wb') as output:
        np.array([0x41584c34,data['p'],data['m'],data['n'],len(patterns),
                  head_net,tail_net,first_negative,24],dtype=np.int32).tofile(output)
        np.array(patterns,dtype=np.uint8).tofile(output)


def build_native(workdir):
    compiler = shutil.which('c++') or shutil.which('clang++') or shutil.which('g++')
    require(compiler is not None, "a C++17 compiler is required")
    binary = Path(workdir)/'hadamard_frame_search'
    subprocess.run([compiler,'-O3','-std=c++17',str(NATIVE_SOURCE),'-o',str(binary)],check=True)
    return binary


def builtin_parity_map():
    return dict(p=631,m=21,mu=458,field_degree=15,field_modulus=0xe7d5,
                inverse_columns=list(PARITY_INVERSE_COLUMNS_631),head_net=1,tail_net=1,
                fixed_first_head_negative=False)


def write_parity_map(path, config):
    degree = config['field_degree']
    Path(path).write_bytes(struct.pack('=7i',0x46503135,3,config['field_modulus'],degree,
                                     int(config['fixed_first_head_negative']),config['head_net'],config['tail_net'])
                          +b''.join(struct.pack('=QQ',column&((1<<64)-1),column>>64)
                                    for column in config['inverse_columns']))


def validate_parity_map(data=None, config=None):
    """Prove the full quadratic identity, not only random samples."""
    config = config or builtin_parity_map()
    data = data or chart(config['p'],config['m'],config['mu'])
    degree = config['field_degree']
    require(type(degree) is int and 1<=degree<=25,"native sampler requires field degree1..25")
    require(type(config['field_modulus']) is int
            and config['field_modulus'].bit_length() == degree+1
            and config['field_modulus']%2 == 1,"invalid binary field modulus")
    require((data['p'],data['m'],data['mu']) == (config['p'],config['m'],config['mu']),
            "parity map belongs to a different multiplier chart")
    variables = 4*degree
    columns = config['inverse_columns']
    require(len(columns) == variables and all(type(c) is int and 0<=c<(1<<variables) for c in columns)
            and data['n'] == 2*degree+1 and len(data['rows']) == degree,
            "map dimensions do not match a single small SL2 factor")
    require(binary_rank(columns) == variables, "singular color map")
    bh,bt = 1+data['m']*config['head_net'],data['m']//3*config['tail_net']
    require(bh*bh+3*bt*bt == data['p'], "DC norm condition failed")
    pairs = [(u,v) for u in range(variables) for v in range(u)]
    inputs = [0]+[1<<u for u in range(variables)]+[(1<<u)|(1<<v) for u,v in pairs]
    colors = [binary_linear(columns,value) for value in inputs]
    axes = np.array([[0]+[(value>>(2*j))&3 for j in range(2*degree)] for value in colors])
    values = np.zeros((len(inputs),degree),dtype=np.uint8)
    for ell,row in enumerate(data['rows']):
        for (i,j,delta),count in row.items():
            if count%2:
                rotated = np.where(axes[:,j] == 0,0,1+(axes[:,j]-1+delta)%3)
                values[:,ell] ^= (axes[:,i] == rotated)
    parity = [sum(int(bit)<<ell for ell,bit in enumerate(row)) for row in values]
    require(all(value == parity[0] for value in parity[1:1+variables]), "unexpected linear term")
    def polar(u,v):
        if u < v:
            u,v = v,u
        return parity[1+variables+pairs.index((u,v))]^parity[1+u]^parity[1+v]^parity[0]
    determinant_map = [polar(3*degree+j,0) for j in range(degree)]
    require(binary_rank(determinant_map) == degree and determinant_map[0] == parity[0],
            "determinant output map is singular or has wrong constant")
    for position,(u,v) in enumerate(pairs):
        term = 0
        if (3*degree<=u<4*degree and 0<=v<degree) or (2*degree<=u<3*degree and degree<=v<2*degree):
            term = field_multiply(1<<(u%degree),1<<(v%degree),degree,config['field_modulus'])
        observed = parity[1+variables+position]^parity[1+u]^parity[1+v]^parity[0]
        require(observed == binary_linear(determinant_map,term), "parity quadratic coefficient failed")
    return dict(color_map_rank=variables,determinant_map_rank=degree,
                full_quadratic_coefficients_checked=len(inputs),
                exact_identity='Cabs mod2=L(p*b1+q*b0+1), L invertible')


def gram(rows, order):
    for i,row in enumerate(rows):
        require(all((row^rows[j]).bit_count() == order//2 for j in range(i)),
                f"integer Gram failed at row{i}")


def packed_rows(matrix):
    require(np.all((matrix == 1)|(matrix == -1)), "non-sign endpoint entry")
    bits = np.packbits(matrix<0,axis=1,bitorder='little')
    return bits,[int.from_bytes(row.tobytes(),'little') for row in bits]


def gs(sequences):
    p = sequences.shape[1]
    indices = np.arange(p)
    A,B,C,D = [row[(indices[None,:]-indices[:,None])%p]
               for row in sequences.astype(np.int8)]
    rev = np.arange(p-1,-1,-1)
    return np.block([[A,B[:,rev],C[:,rev],D[:,rev]],
                     [-B[:,rev],A,-D.T[:,rev],C.T[:,rev]],
                     [-C[:,rev],D.T[:,rev],A,-B.T[:,rev]],
                     [-D[:,rev],-C.T[:,rev],B.T[:,rev],A]])


def verify_endpoint(states, data, output):
    """All-shift coefficients, complete P/Q/cross and complete endpoint Grams."""
    started = time.monotonic()
    output = Path(output)
    require(not output.exists(), "endpoint output already exists")
    # Source BH(3,6): -1 on diagonal, zeta elsewhere. Products are reduced
    # exactly in Z[zeta]/(zeta^2-zeta+1), with conjugate(zeta)=1-zeta.
    source = [[(-1,0) if i == j else (0,1) for j in range(3)] for i in range(3)]
    for i,j in itertools.product(range(3),repeat=2):
        total = [0,0]
        for k in range(3):
            a,b = source[i][k]
            c,d = source[j][k]
            c,d = c+d,-d
            total[0] += a*c-b*d
            total[1] += a*d+b*c+b*d
        require(total == ([3,0] if i == j else [0,0]),"exact source Gram failed")
    b = expand(states,data)
    a = D4@b
    p,d = data['p'],4*data['p']
    require(np.all(a*b == 0) and np.all(a*a+b*b == 1), "coefficient entry test failed")
    for shift in range(p):
        ar,br = np.roll(a,-shift,axis=1),np.roll(b,-shift,axis=1)
        observed = (int(np.sum(a*ar)),int(np.sum(b*br)),int(np.sum(a*br+b*ar)))
        require(observed == ((3*p,p,0) if shift == 0 else (0,0,0)), "coefficient Gram failed")
    P,Q = gs(a-b),gs(a+b)
    _,pr = packed_rows(P)
    _,qr = packed_rows(Q)
    gram(pr,d)
    gram(qr,d)
    r_numerators = {d//2}
    for i in range(d):
        require((pr[i]^qr[i]).bit_count() == d//4, "cross diagonal failed")
        for j in range(i):
            distance1,distance2 = (pr[i]^qr[j]).bit_count(),(qr[i]^pr[j]).bit_count()
            require(distance1+distance2 == d,"cross identity failed")
            r_numerators.update((d-2*distance1,d-2*distance2))
    H = np.block([[-P,Q,Q],[Q,-P,Q],[Q,Q,-P]])
    bits,rows = packed_rows(H)
    gram(rows,3*d)
    output.mkdir(parents=True)
    np.savez_compressed(output/f'h{3*d}.npz',negative_entry_bits=bits,
                        order=np.array(3*d),bitorder=np.array('little'))
    report = dict(p=p,m=data['m'],mu=data['mu'],states=states,endpoint_order=3*d,
                  source_exponents=[[3,1,1],[1,3,1],[1,1,3]],source_gram_scale=3,
                  actual_source_exponents=[1,3],actual_coefficient_alphabet=[-1,0,1],
                  source_field='Q(zeta_6)/Q',field_degree=2,module_rank=d//2,
                  rational_dimension=d,trace_metric=[[2,1],[1,2]],
                  incoming_form='h(u,v)=u^T*v/2-(u^T*(2R-I)*v)*sqrt(-3)/6; Tr h=u^T*v',
                  action='R=Q*P^T/d; R^T=I-R; R^2-R+I=0',frame='F=P; R*F=Q',
                  actual_R_entry_numerators=sorted(r_numerators),actual_R_entry_denominator=d,
                  actual_intermediate_alphabet='{0,-1} union actual_R_entry_numerators/d',
                  frame_gram_scale=d,endpoint_gram_scale=3*d,actual_alphabet=[-1,1],
                  all_coefficient_pair_cross_and_endpoint_integer_Grams_verified=True,
                  packed_sha256=hashlib.sha256(bits.tobytes()).hexdigest(),
                  new_order_claimed=False,new_equivalence_class_claimed=False,
                  elapsed_seconds=time.monotonic()-started)
    (output/'endpoint.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def affine_sign_system(axes, data, head_net=1, tail_net=1, first_negative=0):
    """Independent Python preparation for finite control comparisons."""
    n,L = data['n'],len(data['rows'])
    weights = np.zeros((L,n,n),dtype=np.int32)
    diagonal = np.zeros(L,dtype=np.int32)
    for ell,row in enumerate(data['rows']):
        for (i,j,delta),count in row.items():
            rotated = 0 if axes[j] == 0 else 1+(axes[j]-1+delta)%3
            if axes[i] == rotated:
                if i == j:
                    diagonal[ell] += count
                else:
                    weights[ell,i,j] += count
                    weights[ell,j,i] += count
    absolute = diagonal+weights.sum(axis=(1,2))//2
    require(np.all(absolute%2 == 0), "axis parity failed")
    head = sum(1<<i for i in range(1,n) if axes[i] == 0)
    tail = ((1<<n)-2)^head
    rows = []
    for ell in range(L):
        degree = weights[ell].sum(axis=1)
        rows.append(sum(int(v%2)<<i for i,v in enumerate(degree))
                    | (int(absolute[ell]//2%2)<<n))
    rows += [1,2|(int(first_negative)<<n),head|(((head.bit_count()-head_net)//2%2)<<n),
             tail|(((tail.bit_count()-tail_net)//2%2)<<n)]
    pivot = []
    for col in range(n):
        position = next((i for i in range(len(pivot),len(rows)) if rows[i]>>col&1),None)
        if position is None:
            continue
        at = len(pivot)
        rows[at],rows[position] = rows[position],rows[at]
        for i in range(len(rows)):
            if i != at and rows[i]>>col&1:
                rows[i] ^= rows[at]
        pivot.append(col)
    require(not any((row&((1<<n)-1)) == 0 and row>>n&1 for row in rows),
            "inconsistent sign control")
    particular = sum(((rows[i]>>n)&1)<<col for i,col in enumerate(pivot))
    basis = []
    for free in range(n):
        if free not in pivot:
            basis.append((1<<free)|sum(((rows[i]>>free)&1)<<col for i,col in enumerate(pivot)))
    return dict(particular=particular,basis=basis,weights=weights,diagonal=diagonal)


def affine_signs(axes, data, head_net=1, tail_net=1, first_negative=0):
    system = affine_sign_system(axes,data,head_net,tail_net,first_negative)
    particular,basis = system['particular'],system['basis']
    weights,diagonal = system['weights'],system['diagonal']
    n,L = data['n'],len(data['rows'])
    require(len(basis)<=20,"control enumeration is too large")
    masks = [particular]
    for vector in basis:
        masks += [value^vector for value in masks]
    signs = 1-2*((np.array(masks,dtype=np.uint64)[:,None]
                  >>np.arange(n,dtype=np.uint64))&1).astype(np.int32)
    correlations = np.array([diagonal[ell]+np.einsum('bi,ij,bj->b',signs,weights[ell],signs)//2
                             for ell in range(L)]).T
    return masks,correlations


def self_test(workdir, binary):
    """Complete fixed-k small chart and complete high-bit target filter control."""
    audit = validate_parity_map()
    with tempfile.TemporaryDirectory(prefix='frame-controls-',dir=workdir) as folder:
        folder = Path(folder)
        small = chart(19,3)
        patterns,expected = [],set()
        # First head fixed positive; exactly three of six nonorigin axes are head.
        for extra_heads in itertools.combinations(range(2,7),2):
            tails = [i for i in range(2,7) if i not in extra_heads]
            for colors in itertools.product((1,2,3),repeat=3):
                axes = [0]*7
                for i,color in zip(tails,colors):
                    axes[i] = color
                unsigned = expand([2*a for a in axes],small)
                if any(int(np.sum(unsigned*np.roll(unsigned,-s,axis=1)))%2
                       for s in range(1,19)):
                    continue
                patterns.append(axes)
                for negative in range(1<<5):
                    states = [2*a+((negative>>(i-2))&1 if i>=2 else 0)
                              for i,a in enumerate(axes)]
                    b = expand(states,small)
                    if b.sum(axis=1).tolist() != [4,1,1,1]:
                        continue
                    if all(int(np.sum(b*np.roll(b,-s,axis=1))) == 0 for s in range(1,19)):
                        expected.add(tuple(states))
        require(len(patterns) == 90 and len(expected) == 150,"independent T19 control changed")
        table,pool,result = folder/'table19.bin',folder/'axes19.bin',folder/'result19.json'
        write_tables(table,small)
        write_pool(pool,small,patterns)
        subprocess.run([str(binary),'lift',str(table),str(pool),str(result),'0','1','1','1'],
                       check=True,stdout=subprocess.DEVNULL)
        native = json.loads(result.read_text())
        got = {tuple(2*a+((mask>>i)&1) for i,a in enumerate(patterns[row['pattern']]))
               for row in native['patterns'] for mask in row['negative_masks']}
        require(got == expected,"native finite lift differs from complete T19 brute force")
        control_endpoint = verify_endpoint(list(sorted(expected)[0]),small,folder/'endpoint228')
        data = chart(631,21,458)
        planted = json.loads((ROOT/'artifacts/hadamard-7572/t631.json').read_text())['states']
        axes = [s//2 for s in planted]
        table,pool,result = folder/'table631.bin',folder/'axes631.bin',folder/'result631.json'
        write_tables(table,data)
        write_pool(pool,data,[axes])
        subprocess.run([str(binary),'lift',str(table),str(pool),str(result),'0','1','1','1'],
                       check=True,stdout=subprocess.DEVNULL)
        native = json.loads(result.read_text())['patterns'][0]
        masks,correlations = affine_signs(axes,data)
        expected_mod8 = {x for x,c in zip(masks,correlations) if np.all(c%8 == 0)}
        require(set(native['mod8_negative_masks']) == expected_mod8,
                "high-bit mod8 filter discarded or added a sign vector")
        planted_mask = sum((s%2)<<i for i,s in enumerate(planted))
        require(planted_mask in native['negative_masks'],"planted exact T631 solution was lost")
        port = json.loads((ROOT/'certificates/scalar-lift/hadamard-frame-p271.json').read_text())
        port_audit = validate_parity_map(config=port)
        port_data = chart(port['p'],port['m'],port['mu'])
        table,columns,result = folder/'table271.bin',folder/'columns271.bin',folder/'result271.json'
        write_tables(table,port_data,port['head_net'],port['tail_net'])
        write_parity_map(columns,port)
        subprocess.run([str(binary),'search',str(table),str(columns),str(result),'3','20261124271'],
                       check=True,stdout=subprocess.DEVNULL)
        repeat = json.loads(result.read_text())
        require('states' in repeat,"independent prime repeat did not recover T271")
        repeat_endpoint = verify_endpoint(repeat['states'],port_data,folder/'endpoint3252')
        require(repeat_endpoint['packed_sha256'] == port['known_order_positive_repeat']['endpoint']['packed_sha256'],
                "independent prime replay differs from its recorded exact matrix")
        report = dict(parity_map=audit,complete_T19_axis_patterns=90,
                      small_control_scope='p19,m3,k3,origin/first-head positive,DC(4,1,1,1)',
                      complete_T19_solution_set_size=150,small_endpoint_order=228,
                      small_endpoint_integer_Gram_verified=control_endpoint[
                          'all_coefficient_pair_cross_and_endpoint_integer_Grams_verified'],
                      target_affine_masks_checked=len(masks),
                      target_complete_mod8_survivor_set_size=len(expected_mod8),
                      planted_T631_retained=True,
                      independent_prime_repeat_parity_map=port_audit,
                      independent_prime_repeat_endpoint_order=3252,
                      independent_prime_repeat_full_integer_Gram_verified=True,
                      independent_prime_repeat_packed_sha256=repeat_endpoint['packed_sha256'])
    (Path(workdir)/'self-test.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--seed',type=int,default=1)
    parser.add_argument('--heads',help='nonorigin head counts, or all; default5,7,9 for631')
    parser.add_argument('--parity-map',type=Path,help='validated small-field determinant/color map JSON')
    parser.add_argument('--headless',action='store_true',help='conditionally enumerate zero-head axes, with first tail positive')
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--lift-pool',type=Path,help='JSON finite axis pool with p,m,patterns and DC metadata')
    args = parser.parse_args()
    require(np.isfinite(args.seconds) and args.seconds>0,"seconds must be positive")
    # Bind the actual source inputs before compilation and the long native run.
    source_hashes = dict(native_source_sha256=hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest(),
                         runner_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.workdir.mkdir(parents=True,exist_ok=True)
    binary = build_native(args.workdir)
    if args.self_test:
        print(json.dumps(self_test(args.workdir,binary),indent=2))
        return
    if args.lift_pool:
        pool = json.loads(args.lift_pool.read_text())
        data = chart(pool['p'],pool['m'],pool.get('mu'))
        patterns = [row['axes'] if isinstance(row,dict) else row for row in pool['patterns']]
        table,axes = args.workdir/'tables.bin',args.workdir/'axes.bin'
        write_tables(table,data,pool['head_net'],pool['tail_net'])
        write_pool(axes,data,patterns,pool['head_net'],pool['tail_net'],int(pool['fixed_first_head_negative']))
        subprocess.run([str(binary),'lift',str(table),str(axes),str(args.workdir/'lifts.json'),
                        '0','1','1','1'],check=True)
        found = json.loads((args.workdir/'lifts.json').read_text())
        if 'states' in found:
            print(json.dumps(verify_endpoint(found['states'],data,args.workdir/'verified'),indent=2))
        return
    config = json.loads(args.parity_map.read_text()) if args.parity_map else builtin_parity_map()
    data = chart(config['p'],config['m'],config['mu'])
    audit = validate_parity_map(data,config)
    table,columns = args.workdir/'tables.bin',args.workdir/'columns.bin'
    write_tables(table,data,config['head_net'],config['tail_net'])
    write_parity_map(columns,config)
    result = args.workdir/f'seed-{args.seed}.json'
    require(not result.exists(),"result already exists; choose a new seed or workdir")
    mode = 'headless' if args.headless else 'search'
    command = [str(binary),mode,str(table),str(columns),str(result),str(args.seconds),str(args.seed)]
    heads = args.heads or ('5,7,9' if data['p'] == 631 else 'all')
    if args.headless:
        require(config['head_net'] == 0 and config['tail_net']>0
                and not config['fixed_first_head_negative'],"headless mode requires zero head net and positive tails")
        require(args.heads in (None,'all','0'),"headless mode fixes the head count to zero")
        command.append('search')
    elif heads != 'all':
        command.extend(['0',heads])
    subprocess.run(command,check=True)
    found = json.loads(result.read_text())
    report = dict(parity_map_validation=audit,search_result=str(result),
                  source_capture_phase='before_native_build_and_search',**source_hashes)
    if args.headless:
        report['axis_generator'] = 'all admissible w for each uniformly sampled nonzero SL2 first row'
        report['axis_sampling_scope'] = 'uniform visit intensity; emitted axes are not independent draws'
        report['gauge'] = 'bprime(x)=C^k*b(r*x), positive b(r) chosen; origin and positive tail DC unchanged'
    if 'states' in found:
        report['endpoint'] = verify_endpoint(found['states'],data,args.workdir/f'seed-{args.seed}-verified')
    else:
        report['endpoint_found'] = False
    (args.workdir/f'seed-{args.seed}-run.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
