"""Native exact SL2 parity sampler with a serializable deterministic stream.

Only field arithmetic, binary coordinate conversion and head-count rejection
run natively. The established Python normalization transports orbit phases.
A parity word supplies coordinate labels, not a completed signed T-design.
"""
from __future__ import annotations

import ctypes
import hashlib
from pathlib import Path
import shutil
import subprocess

import numpy as np

from experiments.scalar_lift import hadamard_frame_search as base
from experiments.scalar_lift.large_parity_map import BinaryField, normalize_head

MASK64=(1<<64)-1


def build_native(directory):
    """Compile a source-addressed library; return its loaded interface."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).with_suffix('.cpp')
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    snapshot=directory/f'exact-parity-sampler-{digest[:16]}.cpp'
    library=snapshot.with_suffix('.dylib')
    if not snapshot.exists():snapshot.write_bytes(source.read_bytes())
    if not library.exists():
        compiler=shutil.which('c++') or shutil.which('clang++')
        if compiler is None:raise RuntimeError('native parity sampling requires C++17')
        temporary=library.with_suffix('.building.dylib')
        subprocess.run([compiler,'-O3','-std=c++17','-shared','-fPIC',str(snapshot),'-o',str(temporary)],check=True)
        temporary.replace(library)
    native=ctypes.CDLL(str(library.resolve()))
    u64=ctypes.c_uint64;uint=ctypes.c_uint
    words=np.ctypeslib.ndpointer(np.uint64,flags='C_CONTIGUOUS')
    native.exact_parity_field.argtypes=[uint,u64,u64,u64,ctypes.POINTER(u64),ctypes.POINTER(u64)]
    native.exact_parity_field.restype=ctypes.c_int
    native.exact_parity_colors.argtypes=[uint,u64,words,u64,u64,u64,words]
    native.exact_parity_colors.restype=ctypes.c_int
    native.exact_parity_sample.argtypes=[uint,u64,words,np.ctypeslib.ndpointer(np.uint8,flags='C_CONTIGUOUS'),ctypes.c_int,u64,ctypes.POINTER(u64),ctypes.POINTER(u64),words]
    native.exact_parity_sample.restype=ctypes.c_int
    return native


class NativeDeterminantSampler:
    """The DeterminantSampler interface with an independent xoshiro256** stream.

    Resume states are explicitly tagged; they cannot be mistaken for Python
    Random states. The wrapper's exact first-head normalization is unchanged.
    """
    def __init__(self,config,seed=None,*,workdir):
        self.config=config;self.degree=int(config['field_degree'])
        self.field=BinaryField(self.degree,int(config['field_modulus']))
        if not 1<=self.degree<=63 or not self.field.irreducible():
            raise ValueError('native sampler requires an irreducible binary field of degree1..63')
        self.columns=config['inverse_columns']
        if len(self.columns)!=4*self.degree or any(type(x) is not int or not 0<=x<1<<(4*self.degree) for x in self.columns):
            raise ValueError('native sampler requires a full bounded binary color map')
        self.packed_columns=np.array([[(value>>(64*k))&MASK64 for k in range(4)] for value in self.columns],dtype=np.uint64)
        self.native=build_native(workdir)
        seed_state=(0 if seed is None else int(seed))&MASK64
        seed_words=[]
        for _ in range(4):
            seed_state=(seed_state+0x9e3779b97f4a7c15)&MASK64
            z=seed_state;z=((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK64
            z=((z^(z>>27))*0x94d049bb133111eb)&MASK64
            seed_words.append(z^(z>>31))
        self.state=(ctypes.c_uint64*4)(*seed_words)
        self.draw_count=ctypes.c_uint64(0)
        self.rng=self;self._chart=None
        self.nonorigin_orbits=2*self.degree

    @property
    def draws(self):return self.draw_count.value

    @draws.setter
    def draws(self,value):
        if not 0<=value<=MASK64:raise ValueError('native draw count out of range')
        self.draw_count.value=value

    def getstate(self):return ('exact-parity-xoshiro256ss-v1',tuple(self.state))

    def setstate(self,value):
        if not isinstance(value,(list,tuple)) or len(value)!=2 or value[0]!='exact-parity-xoshiro256ss-v1' or not isinstance(value[1],(list,tuple)) or len(value[1])!=4 or any(type(v) is not int or not 0<=v<=MASK64 for v in value[1]) or not any(value[1]):
            raise ValueError('incorrect native sampler resume state')
        self.state[:]=value[1]

    def field_operation(self,a,b):
        product=ctypes.c_uint64();reciprocal=ctypes.c_uint64()
        if not 0<=a<1<<self.degree or not 0<=b<1<<self.degree:raise ValueError('field input out of range')
        error=self.native.exact_parity_field(self.degree,self.field.modulus,a,b,ctypes.byref(product),ctypes.byref(reciprocal))
        if error:raise ArithmeticError(f'native field operation failed: {error}')
        return product.value,reciprocal.value

    def color_parameters(self,p,q,w):
        if not all(type(v) is int and 0<=v<1<<self.degree for v in (p,q,w)) or not p|q:
            raise ValueError('bounded nonzero determinant first row required')
        output=np.zeros(4,dtype=np.uint64)
        error=self.native.exact_parity_colors(self.degree,self.field.modulus,self.packed_columns,p,q,w,output)
        if error:raise ArithmeticError(f'native coordinate conversion failed: {error}')
        return sum(int(word)<<(64*k) for k,word in enumerate(output))

    def head_count(self,color):
        mask=((1<<(2*self.nonorigin_orbits))-1)//3
        return self.nonorigin_orbits-((color|(color>>1))&mask).bit_count()

    def _sample(self,heads,first_head,max_draws):
        allowed=np.ones(self.nonorigin_orbits+1,dtype=np.uint8) if heads is None else np.zeros(self.nonorigin_orbits+1,dtype=np.uint8)
        if heads is not None:
            for h in heads:
                if type(h) is not int or not 0<=h<=self.nonorigin_orbits:raise ValueError('head count out of range')
                allowed[h]=1
        if not np.any(allowed):return None
        if max_draws is None:max_draws=MASK64
        if type(max_draws) is not int or not 0<=max_draws<=MASK64:raise ValueError('draw budget out of range')
        output=np.zeros(4,dtype=np.uint64)
        error=self.native.exact_parity_sample(self.degree,self.field.modulus,self.packed_columns,allowed,int(first_head),max_draws,self.state,ctypes.byref(self.draw_count),output)
        if error==1:return None
        if error:raise ArithmeticError(f'native sampler failed: {error}')
        return sum(int(word)<<(64*k) for k,word in enumerate(output))

    def draw_packed(self):return self._sample(None,False,1)

    def next_axes(self,heads=None,first_head=False,normalize_first_head=False,max_draws=None):
        color=self._sample(heads,first_head and not normalize_first_head,max_draws)
        if color is None:return None
        axes=[0]+[(color>>(2*j))&3 for j in range(self.nonorigin_orbits)]
        if normalize_first_head:
            if self._chart is None:
                c=self.config;self._chart=base.chart(c['p'],c['m'],c['mu'])
            axes,_scale=normalize_head(axes,self._chart)
        return axes
