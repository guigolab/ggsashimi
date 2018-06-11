#!/usr/bin/env python
import re
import importlib 
from collections import OrderedDict

sp = importlib.import_module('sashimi-plot')

def test_parse_coordinates():
    for s in ['chr1:1000-2000', 'chr1:1,000-2,000']:
        assert sp.parse_coordinates(s) == ("chr1", 999, 2000)

def test_count_operator():
    pos = 27037633
    _, start, end = sp.parse_coordinates('chr10:27035000-27050000')
    
    c = [0] * (end - start)
    j = OrderedDict()
    
    # soft clip
    new_pos = sp.count_operator('S', 1, pos, start, end, c, j)
    assert new_pos == pos
    assert all(v == 0 for v in c)
    assert len(j) == 0
    pos = new_pos

    # match
    new_pos = sp.count_operator('M', 42, pos, start, end, c, j)
    assert new_pos == pos+42
    p = [v == 1 for v in c[pos-start:pos-start+42]]
    assert len(p) == 42
    assert all(p)
    assert len(j) == 0
    pos = new_pos
    
    # skip
    new_pos = sp.count_operator('N', 2852, pos, start, end, c, j) 
    assert new_pos == pos+2852
    p = [v == 0 for v in c[pos-start:pos-start+2852]]
    assert len(p) == 2852
    assert all(p)
    assert len(j) == 1
    assert j[(pos,pos+2852)] == 1
    pos = new_pos

    # match
    new_pos = sp.count_operator('M', 58, pos, start, end, c, j) 
    assert new_pos == pos+58
    p = [v == 1 for v in c[pos-start:pos-start+58]]
    assert len(p) == 58
    assert all(p)
    assert j[(pos-2852,pos)] == 1

def test_flip_read():
    assert sp.flip_read('NONE', 4) == 0
    assert sp.flip_read('SENSE', 4) == 0
    assert sp.flip_read('MATE1_SENSE', 64) == 0
    assert sp.flip_read('MATE2_SENSE', 128) == 0
    assert sp.flip_read('ANTISENSE', 4) == 1
    assert sp.flip_read('MATE1_SENSE', 128) == 1
    assert sp.flip_read('MATE2_SENSE', 64) == 1

def test_intersect_introns():
    data = [
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27047991), 
        (27040713, 27044584), 
        (27044671, 27047991), 
        (27040713, 27047991), 
        (27040713, 27044584), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991), 
        (27040713, 27044584), 
        (27040713, 27047991), 
        (27044671, 27047991)
    ]
    
    i = list(sp.intersect_introns(data))
    assert len(i) == 2
    assert i == [(27040713, 27044584), (27044671, 27047991)]
