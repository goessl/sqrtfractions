from math import sqrt, prod, sumprod
from collections import defaultdict
from operator import mul
from sympy.ntheory import factorint, primefactors
import sympy as sp
from itertools import repeat, chain
from functools import total_ordering, reduce, cached_property
from collections.abc import KeysView, ValuesView, ItemsView
from random import randint
from typing import Union, Any



@total_ordering
class SumOfSqrts:
    #creation
    def __init__(self, n:dict[int,int]|int=0):
        if isinstance(n, int):
            n = {1:n}
        elif isinstance(n, dict):
            if not all(isinstance(k, int) and isinstance(v, int) \
                    for k, v in n.items()):
                raise TypeError('Radicands and factors must be integers.')
            if not all(k>=0 for k in n.keys()):
                raise ValueError('Radicands must be equal to or greater than zero.')
        else:
            raise TypeError('Numerator must be an integer or a dictionary.')
        
        _n = defaultdict(int)
        for k, v in n.items():
            f = factorint(k)
            s = prod(p**(e//2) for p, e in f.items())
            k = prod(p**(e%2) for p, e in f.items())
            _n[k] += s * v
        _n.pop(0, None) #remove possible v*sqrt(0) item
        _n = {k:v for k, v in _n.items() if v}
        
        self.n = dict(sorted(_n.items()))
    
    @staticmethod
    def _create_directy(n:dict) -> 'SumOfSqrts':
        s = object.__new__(SumOfSqrts)
        s.n = n
        return s
    
    @staticmethod
    def random(N:int=10, precision:int=20) -> 'SumOfSqrts':
        n = {n:randint(-precision, +precision) for n in range(1, N+1)}
        return SumOfSqrts(n)
    
    
    
    def validate(self) -> bool:
        for k, v in self.items():
            if not isinstance(k, int) or not isinstance(v, int) \
                    or not v or not k>0 \
                    or any(v>=2 for k, v in factorint(k).items()):
                return False
        return True
    
    
    
    #conversion
    @cached_property
    def _float(self) -> float:
        return float(sumprod(self.values(), map(sqrt, self.keys())))
    
    def __float__(self) -> float:
        return self._float
    
    def is_integer(self) -> bool:
        return set(self.keys()) <= {1}
    
    def __int__(self) -> int:
        if self.is_integer():
            return self.n.get(1, 0)
        else:
            raise ValueError('doesn\'t represent an integer')
    
    def sympify(self) -> sp.Expr:
        return sp.sympify(sumprod(self.values(), map(sp.sqrt, self.keys())))
    
    
    
    #container
    def __len__(self) -> int:
        return len(self.n)
    
    def keys(self) -> KeysView[int]:
        return self.n.keys()
    
    def values(self) -> ValuesView[int]:
        return self.n.values()
    
    def items(self) -> ItemsView[int, int]:
        return self.n.items()
    
    
    
    #ordering
    def __eq__(self, other:Any) -> bool: #other:Union['SumOfSqrts', int]
        if isinstance(other, SumOfSqrts):
            return self.n == other.n
        elif isinstance(other, int):
            try:
                return int(self) == other
            except ValueError:
                return False
        return NotImplemented
    
    def __abs__(self) -> 'SumOfSqrts':
        return self if self>=0 else -self
    
    def __lt__(self, other:Any) -> bool: #other:Union['SumOfSqrts', int]
        #https://math.stackexchange.com/a/1076510
        if isinstance(other, (SumOfSqrts, int)):
            l = self - other
            
            while not l.is_integer():
                p = max(chain(*(primefactors(k) for k in l.keys())))
                
                r = SumOfSqrts({k//p:-v for k, v in l.items() if k%p==0})
                l = SumOfSqrts({k:v for k, v in l.items() if k%p!=0})
                #https://math.stackexchange.com/a/2347212
                l, r = l*abs(l), r*abs(r)*p
                
                l -= r
            
            return int(l) < 0
        return NotImplemented
    
    
    
    #arithmetic
    def __neg__(self) -> 'SumOfSqrts':
        return SumOfSqrts._create_directy({k:-v for k, v in self.items()})
    
    def __add__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        if isinstance(other, int):
            other = SumOfSqrts(other)
        if isinstance(other, SumOfSqrts):
            n = defaultdict(int, self.n)
            for k, v in other.items():
                n[k] += v
                if not n[k]:
                    del n[k]
            return SumOfSqrts._create_directy(n)
        return NotImplemented
    __radd__ = __add__
    
    def __iadd__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        if isinstance(other, int):
            other = SumOfSqrts(other)
        if isinstance(other, SumOfSqrts):
            for k, v in other.items():
                self.n[k] = self.n.get(k, 0) + v
                if not self.n[k]:
                    del self.n[k]
            return self
        return NotImplemented
    
    
    def __sub__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        return self + (-other)
    
    def __rsub__(self, other:Any) -> 'SumOfSqrts': #other:int
        return (-self) + other
    
    def __isub__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        self += -other
        return self
    
    
    def __mul__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        if isinstance(other, SumOfSqrts):
            n = defaultdict(int)
            for ki, vi in self.items():
                for kj, vj in other.items():
                    n[ki*kj] += vi * vj
            return SumOfSqrts(n)
        elif isinstance(other, int):
            if not other:
                return SumOfSqrts()
            return SumOfSqrts._create_directy({k:v*other for k, v in self.items()})
        return NotImplemented
    __rmul__ = __mul__
    
    def __imul__(self, other:Any) -> 'SumOfSqrts': #other:Union['SumOfSqrts', int]
        if isinstance(other, int):
            if not other:
                self.n.clear()
                return self
            for k in self.keys():
                self.n[k] *= other
            return self
        return NotImplemented
    
    
    def __floordiv__(self, other:Any) -> 'SumOfSqrts': #other:int
        if isintance(other, int):
            return SumOfSqrts._create_directy(
                    {k:v//other for k, v in self.items() if v//other})
        return NotImplemented
    
    def __ifloordiv__(self, other:Any) -> 'SumOfSqrts': #other:int
        if isintance(other, int):
            for k in tuple(self.keys()):
                self.n[k] //= other
                if not self.n[k]:
                    del self.n[k]
            return self
        return NotImplemented
    
    
    def __mod__(self, other:Any) -> 'SumOfSqrts': #other:int
        if isintance(other, int):
            return SumOfSqrts._create_directy(
                    {k:v%other for k, v in self.items() if v%other})
        return NotImplemented
    
    def __imod__(self, other:Any) -> 'SumOfSqrts': #other:int
        if isintance(other, int):
            for k in tuple(self.keys()):
                self.n[k] %= other
                if not self.n[k]:
                    del self.n[k]
            return self
        return NotImplemented
    
    
    def __pow__(self, other:Any) -> 'SumOfSqrts': #other:int
        if isinstance(other, int):
            if other >= 0:
                return reduce(mul, repeat(self, other), SumOfSqrts(1))
            else:
                raise ValueError('Exponent must be non-negative.')
        return NotImplemented
    
    
    
    #IO
    def __repr__(self) -> str:
        n = [f'{n:+}{chr(0x221A)}{r}' for r, n in self.items()]
        if len(n) <= 1: #no parentheses needed
            return n[0] if n else '0'
        return '(' + ''.join(n) + ')'
    
    def _repr_latex_(self) -> str:
        n = [f'{n:+d}\\sqrt{{{r}}}' for r, n in self.items()]
        return '$' + (''.join(n) if n else '0') + '$'
