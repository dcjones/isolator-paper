#!/usr/bin/env python

'''
A simple GTF (Gene Transfer Format) parser.
'''

import re
from collections import namedtuple


attrib_pat = re.compile(
        r'\s*(\S+)\s+(\"([^\"]*)\"|([^\"\s]*))[\s*;]?')


row = namedtuple('row',
    '''
    seqname
    name
    feature
    start
    end
    score
    strand
    frame
    attributes
    ''')


def parse(f):
    '''
    Parse gtf rows from f, where f is something that generates strings (e.g.
    an open file).
    '''

    for line in f:
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == '#': continue

        fields = line.split('\t')
        if len(fields) != 9:
            raise Exception('Malformed GTF')

        attribs = attrib_pat.findall(fields[8])

        yield row(seqname = fields[0],
                  name    = fields[1],
                  feature = fields[2],
                  start   = int(fields[3]),
                  end     = int(fields[4]),
                  score   = fields[5],
                  strand  = fields[6],
                  frame   = int(fields[7]) if fields[7] != '.' else None,
                  attributes =
                    dict((mat[0], mat[2]) for mat in attribs))


