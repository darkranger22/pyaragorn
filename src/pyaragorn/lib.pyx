from libc.stdio cimport fopen, fdopen, fclose, stdout
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy

cimport libaragorn
from libaragorn cimport csw, data_set, gene

import os



cdef extern from * nogil:
    """
    void default_sw(csw* sw) {
        csw x = {
            {"tRNA", "tmRNA", "", "", "CDS", "overall"}, 
            NULL, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, STANDARD, 0, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
            0, METAZOAN_MT, 1, 0, 5, 5, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 
            3, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
            {0, 0, 0, 0, 0, 0}, 0, 0, 0, 0, NTAG, 10, 30, 
            {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
            0, 0, 0, 0, 0L, 100.0, 1.0, tRNAthresh, 4.0, 29.0, 26.0, 7.5, 8.0, 
            mtRNAtthresh, mtRNAdthresh, mtRNAdtthresh, -7.9, -6.0, tmRNAthresh, 
            14.0, 10.0, 25.0, 9.0, srpRNAthresh, CDSthresh, 
            {tRNAthresh, tmRNAthresh, srpRNAthresh, 0.0, CDSthresh}, 
            {
                45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 
                45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10, 65, 
                82, 65, 71, 79, 82, 78, 32, 118, 49, 46, 50, 46, 52, 49, 32, 
                32, 32, 68, 101, 97, 110, 32, 76, 97, 115, 108, 101, 116, 116, 
                10, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 
                45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10, 
                TERM
            }
        };
        memcpy(sw, &x, sizeof(csw));
    }
    """
    void default_sw(csw* sw)



cdef class Gene:
    cdef gene _gene

    @property
    def type(self):
        return ["tRNA", "tmRNA", "", "", "CDS"][<int> self._gene.genetype]

    @property
    def name(self):
        return self._gene.name

    @property
    def seq(self):
        return self._gene.seq

    @property
    def energy(self):
        return self._gene.energy



cdef class RNAFinder:
    cdef csw _sw

    def __init__(self):
        default_sw(&self._sw)
        self._sw.trna = True
        self._sw.tmrna = True
        self._sw.f = stdout
        self._sw.verbose = True
        
        self._sw.genespace = libaragorn.NT
        self._sw.genes = <gene*> calloc(self._sw.genespace, sizeof(gene))
        if self._sw.genes is NULL:
            raise MemoryError("failed to allocate data")

    def __dealloc__(self):
        free(self._sw.genes)

    def run(self, filename):

        cdef Gene     g
        cdef data_set ds
        cdef int      fd

        with open(filename) as f:
            fd = f.fileno()
            ds.bugmode = False
            ds.f = fdopen(fd, b"r")

            with nogil:
                libaragorn.bopt_fastafile(&ds, &self._sw)

        genes = []
        for i in range(89):
            g = Gene.__new__(Gene)
            memcpy(&g._gene, &self._sw.genes[i], sizeof(gene))
            genes.append(g)
        
        return genes