"""Evaluate frozen KO rules at community level with explicit unknown values.

Import parse, evaluate, kos and self_tests from this module.
Execute this file to run the logic checks.
No network calls. No module direction or same-genome completeness is inferred.
"""
from pathlib import Path
import argparse, hashlib, json, re, sys, platform
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parent

def split_top(s,mode):
    depth=0; start=0; out=[]
    for i,c in enumerate(s):
        if c=='(': depth+=1
        elif c==')': depth-=1
        if depth<0: raise ValueError('Unbalanced parentheses')
        if depth==0 and ((mode=='space' and c.isspace()) or c==mode):
            if s[start:i].strip():out.append(s[start:i].strip())
            start=i+1
    if depth:raise ValueError('Unbalanced parentheses')
    if s[start:].strip():out.append(s[start:].strip())
    return out

def unwrap(s):
    while s.startswith('(') and s.endswith(')'):
        depth=0; closed=False
        for i,c in enumerate(s):
            depth+=(c=='(')-(c==')')
            if depth==0:
                closed=(i==len(s)-1);break
        if not closed:break
        s=s[1:-1].strip()
    return s

def parse(s,modules=None,stack=()):
    s=unwrap(s.strip())
    if not s:raise ValueError('Empty expression')
    # KEGG: top-level pathway blocks are each parenthesized, then ANDed.
    parts=split_top(s,'space')
    if len(parts)>1:return ('and',[parse(p,modules,stack) for p in parts])
    parts=split_top(s,',')
    if len(parts)>1:return ('or',[parse(p,modules,stack) for p in parts])
    depth=0; start=0; sign='+'; terms=[]
    for i,c in enumerate(s):
        depth+=(c=='(')-(c==')')
        if depth==0 and c in '+-':
            if i>start:terms.append((sign,s[start:i]))
            elif i!=0:raise ValueError('Consecutive operators')
            sign=c;start=i+1
    if start>0:
        terms.append((sign,s[start:]))
        children=[('optional',[parse(v,modules,stack)]) if op=='-' else parse(v,modules,stack) for op,v in terms]
        return ('and',children)
    if re.fullmatch(r'K\d{5}',s):return ('ko',s)
    if re.fullmatch(r'M\d{5}',s):
        if modules is None or s not in modules:raise ValueError('Unresolved module '+s)
        if s in stack:raise ValueError('Recursive module '+s)
        return parse(modules[s],modules,stack+(s,))
    raise ValueError('Unsupported token: '+s)

def kos(tree,include_optional=False):
    op,x=tree
    if op=='ko':return {x}
    if op=='optional' and not include_optional:return set()
    return set().union(*(kos(c,include_optional) for c in x))

def evaluate(tree,values,n):
    op,x=tree
    if op=='ko':return values.get(x,np.full(n,2,dtype=np.int8))
    if op=='optional':return np.ones(n,dtype=np.int8)
    v=np.stack([evaluate(c,values,n) for c in x])
    if op=='and':return np.where((v==0).any(axis=0),0,np.where((v==2).any(axis=0),2,1))
    if op=='or':return np.where((v==1).any(axis=0),1,np.where((v==2).any(axis=0),2,0))
    raise ValueError(op)

def self_tests():
    tests=[
        ('K00001+K00002-K00003',{'K00001':1,'K00002':1},1),
        ('K00001+K00002-K00003',{'K00001':1,'K00002':0},0),
        ('K00001,K00002',{'K00001':0},2),
        ('K00001,K00002',{'K00001':1},1),
        ('K00001 K00002',{'K00001':0},0),
        ('K00001 K00002',{'K00001':1},2),
        ('(K00001 K00002),K00003',{'K00001':0,'K00002':1,'K00003':1},1),
        ('K00001,K00002 K00003',{'K00001':1,'K00002':0,'K00003':0},0),
        ('(K00001+K00002-K00003)+(K00004,K00005)',{'K00001':1,'K00002':1,'K00005':1},1),
        ('K00001-(K00002+K00003)',{'K00001':1},1),
    ]
    for expr,d,expected in tests:
        actual=evaluate(parse(expr),{k:np.array([v]) for k,v in d.items()},1)[0]
        assert actual==expected,(expr,actual,expected)
    assert kos(parse('K00001+K00002-K00003'))=={'K00001','K00002'}
    assert kos(parse('K00001+K00002-K00003'),True)=={'K00001','K00002','K00003'}
    assert evaluate(parse('M00001',{'M00001':'K00001,K00002'}),{'K00002':np.array([1])},1)[0]==1
    for bad in ['K00001/K00002','K00001+','K00001++K00002','(K00001']:
        try:parse(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid expression accepted: '+bad)
    # Exhaustively compare all ternary inputs with independent scalar logic.
    import itertools
    for a,b,c in itertools.product((0,1,2),repeat=3):
        ab=0 if 0 in (a,b) else 2 if 2 in (a,b) else 1
        want=1 if 1 in (ab,c) else 2 if 2 in (ab,c) else 0
        actual=evaluate(parse('(K00001+K00002),K00003'),dict(zip(['K00001','K00002','K00003'],map(lambda x:np.array([x]),[a,b,c]))),1)[0]
        assert want==actual
    return len(tests)+27+7

if __name__ == "__main__":
    print(f"Rule-logic checks passed: {self_tests()}")
