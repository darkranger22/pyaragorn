from cython.operator cimport postincrement, dereference
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.exc cimport PyErr_CheckSignals
from cpython.unicode cimport (
    PyUnicode_KIND,
    PyUnicode_DATA,
    PyUnicode_1BYTE_KIND,
    PyUnicode_GET_LENGTH,
    # PyUnicode_READ,
)

from libc.stdio cimport FILE, fopen, fdopen, fclose, fprintf, fputc, stdout, stderr
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy

cimport aragorn
from aragorn cimport csw, data_set, gene

import os

# --- Helpers ------------------------------------------------------------------

cdef extern from * nogil:
    Py_UCS4 PyUnicode_READ(int kind, const void* data, size_t pos)

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


cdef inline long int sq(data_set* d, long int pos) nogil:
    return (pos + d.psmax - 1) % d.psmax + 1

# --- Constants ----------------------------------------------------------------

cdef set _TRANSLATION_TABLES  = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 27)) | {29, 30} | {32, 33}


# --- Classes ------------------------------------------------------------------

cdef class Gene:

    cdef gene _gene
    cdef int  _genetic_code

    @staticmethod
    cdef Gene _new_gene(gene* _gene, int _genetic_code):
        cdef Gene obj

        if _gene.genetype == aragorn.tRNA:
            obj = TRNAGene.__new__(TRNAGene)
        elif _gene.genetype == aragorn.tmRNA:
            obj = TMRNAGene.__new__(TMRNAGene)
        else:
            raise NotImplementedError

        memcpy(&obj._gene, _gene, sizeof(gene))
        obj._genetic_code = _genetic_code
        return obj

    @property
    def type(self):
        return ["tRNA", "tmRNA", "", "", "CDS"][<int> self._gene.genetype]

    @property
    def start(self):
        return self._gene.start

    @property
    def stop(self):
        return self._gene.stop

    @property
    def length(self):
        return aragorn.seqlen(&self._gene)

    @property
    def sequence(self):
        cdef long  i
        cdef int   length = aragorn.seqlen(&self._gene)
        return ''.join([
            chr(aragorn.cpbase(self._gene.seq[i]))
            for i in range(length)
        ])

    @property
    def energy(self):
        cdef csw sw
        default_sw(&sw)
        return aragorn.nenergy(&self._gene, &sw)

    @property
    def strand(self):
        return -1 if self._gene.comp else +1




cdef class TRNAGene(Gene):

    @property
    def amino_acid(self):
        cdef csw sw
        cdef int* s = self._gene.seq + self._gene.anticodon
        (<int*> &sw.geneticcode)[0] = self._genetic_code
        if self._gene.cloop == 6:
            return (
                aragorn.aa(s - 1, &sw).decode('ascii'),
                aragorn.aa(s, &sw).decode('ascii'),
            )
        elif self._gene.cloop == 8:
            return (
                aragorn.aa(s, &sw).decode('ascii'),
                aragorn.aa(s + 1, &sw).decode('ascii')
            )
        else:
            return aragorn.aa(s, &sw).decode('ascii')

    @property
    def anticodon(self):
        cdef tuple c
        cdef int*  s = self._gene.seq + self._gene.anticodon
        if self._gene.cloop == 6:
            c = ( aragorn.cbase(s[0]), aragorn.cbase(s[1]) )
        elif self._gene.cloop == 8:
            c = ( aragorn.cbase(s[0]), aragorn.cbase(s[1]), aragorn.cbase(s[2]), aragorn.cbase(s[3]) )
        else:
            c = ( aragorn.cbase(s[0]), aragorn.cbase(s[1]), aragorn.cbase(s[2]) )
        return ''.join(map(chr, c))

    @property
    def anticodon_offset(self):
        cdef int x = 1 + self._gene.anticodon
        if self._gene.nintron > 0 and self._gene.intron <= self._gene.anticodon:
            x += self._gene.nintron
        return x


cdef class TMRNAGene(Gene):
    pass




cdef int[256] _map = [
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4, 5,-3,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-2,-4,-4, 0, 4, 1, 4,-4,-4, 2, 4, -4,-4, 4,-5, 4, 4,
    -4,-4,-4, 4, 4, 3, 3, 4, 4,-4, 4,-4,-4,-4,-4,-2,5,-4,0,4,1,4,
    -4,-4,2,4,-4,-4,4,-5,4,4,-4,-4,-4, 4,4,3,3,4,4,-4,
    4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4
]


cdef class Cursor:
    cdef object obj
    cdef const void*  data
    cdef int    kind
    cdef size_t length
    cdef size_t pos

    def __init__(self, obj):
        cdef const unsigned char[::1] view
        if isinstance(obj, str):
            self.kind = PyUnicode_KIND(obj)
            self.data = PyUnicode_DATA(obj)
            self.length = PyUnicode_GET_LENGTH(obj)
        else:
            view = obj
            self.kind = PyUnicode_1BYTE_KIND
            self.data = &view[0]
            self.length = view.shape[0]
        self.obj = obj
        self.pos = 0

    cdef int _forward(self, data_set* d) noexcept nogil:
        cdef Py_UCS4 x
        cdef int     base

        if self.pos >= self.length:
            return -1 #aragorn.TERM

        x = PyUnicode_READ(self.kind, self.data, self.pos)
        self.pos += 1

        if x >= 128:
            return 5 #aragorn.NOBASE

        base = _map[x]
        if base >= 0:
            d.ps += 1
            return base
        else:
            return 5 #aragorn.NOBASE

    cdef double _gc(self, data_set* d) noexcept nogil:
        cdef int  base
        cdef long i
        cdef long ngc  = 0

        for i in range(self.length):
            base = self._forward(d)
            if base == -1:
                break
            ngc += (base == <int> aragorn.base.Cytosine) or (base == <int> aragorn.base.Guanine)

        return <double> ngc / <double> d.ps


cdef class RNAFinder:
    cdef csw _sw

    def __init__(
        self,
        int translation_table = 1,
    ):
        default_sw(&self._sw)
        self._sw.trna = True
        self._sw.tmrna = True
        self._sw.f = stdout
        self._sw.verbose = False #True

        if translation_table not in _TRANSLATION_TABLES:
            raise ValueError(f"invalid translation table: {translation_table!r}")
        self._sw.geneticcode = translation_table

    def find_rna(self, object sequence):
        """Find RNA genes in the input DNA sequence.

        Arguments:
            sequence (`str` or buffer): The nucleotide sequence to process, 
                either as a string of nucleotides (upper- or lowercase), or
                as an object implementing the buffer protocol.

        Returns:
            `list` of `~pyaragorn.Gene`: A list of `Gene` corresponding to 
            RNA genes detected in the sequence according to the `RNAFinder` 
            parameters.

        """
        cdef Gene     g
        cdef data_set ds
        cdef int      n
        cdef int      nt
        cdef csw      sw
        cdef int*     vsort  = NULL
        cdef Cursor   cursor = Cursor(sequence)

        # copy parameters to ensure the `find_rna` method is re-entrant
        memcpy(&sw, &self._sw, sizeof(csw))

        try:
            with nogil:
                # allocate memory for the result genes
                sw.genespace = aragorn.NT
                sw.genes = <gene*> calloc(sw.genespace, sizeof(gene))
                if sw.genes is NULL:
                    raise MemoryError("failed to allocate memory")
                # detect RNA genes with the "batched" algorithm
                nt = self._bopt(cursor, &ds, &sw)
                # allocate array for sorting genes
                vsort = <int*> calloc(nt, sizeof(int))
                if vsort is NULL:
                    raise MemoryError("failed to allocate memory")
                # sort and threshold genes
                n = aragorn.gene_sort(&ds, nt, vsort, &sw)
            # recover genes
            genes = []
            for i in range(n):
                genes.append(Gene._new_gene(&sw.genes[vsort[i]], sw.geneticcode))
        finally:
            free(vsort)
            free(sw.genes)

        return genes

    cdef int _bopt(
        self,
        Cursor cursor,
        data_set* d,
        csw* sw
    ) except -1 nogil:
        # adapted from bopt_fastafile to use with our own
        cdef int nt
        cdef int seq[((2 * aragorn.LSEQ) + aragorn.WRAP) + 1]
        cdef int cseq[((2 * aragorn.LSEQ) + aragorn.WRAP) + 1]
        cdef int wseq[(2 * aragorn.WRAP) + 1]

        cdef long ngc
        cdef int base

        cdef long rewind
        cdef long drewind
        cdef long tmaxlen

        # compute width of sliding windows
        rewind = aragorn.MAXTAGDIST + 20
        if sw.trna or sw.mtrna:
            tmaxlen = aragorn.MAXTRNALEN + sw.maxintronlen
            if rewind < tmaxlen:
                rewind = tmaxlen
        if sw.tmrna:
            if rewind < aragorn.MAXTMRNALEN:
                rewind = aragorn.MAXTMRNALEN
        if sw.peptide:
            if sw.tagthresh >= 5 and rewind < aragorn.TSWEEP:
                rewind = aragorn.TSWEEP

        sw.loffset = rewind
        sw.roffset = rewind
        drewind = 2 * rewind

        # reinitialize dataset book-keeping 
        d.filepointer = 0
        d.ns = 0
        d.nf = 0
        d.nextseq = 0L
        d.nextseqoff = 0L
        d.seqstart = 0
        d.seqstartoff = 0
        d.ps = 0

        # count GC% and compute sequence length
        d.gc = cursor._gc(d)
        d.psmax = d.ps

        # reset dataset / cursor position
        d.ps = 0
        cursor.pos = 0

        # cleanly initialize gene array
        aragorn.init_gene(sw.genes, 0, aragorn.NT)

        cdef int i
        cdef bint flag
        cdef int length
        cdef int *s
        cdef int *sf
        cdef int *se
        cdef int *sc
        cdef int *swrap
        cdef long gap
        cdef long start
        cdef long vstart
        cdef long vstop
        cdef double sens
        cdef bint loop
        cdef bint NX
        cdef bint SH

        nt = 0
        flag = 0
        start = 1L

        loop = True
        NX = True
        SH = True

        se = seq
        if sw.linear:
            for i in range(rewind):
                postincrement(se)[0] = aragorn.NOBASE
            start -= rewind
        else:
            if d.psmax <= drewind:
                gap = drewind - d.psmax
                sc = se + gap
                while se < sc:
                    postincrement(se)[0] = aragorn.NOBASE

                swrap = wseq
                sc = se + d.psmax
                while se < sc:
                    se[0] = cursor._forward(d) #aragorn.move_forward(d)
                    postincrement(swrap)[0] = postincrement(se)[0]

                sc = swrap + gap
                while swrap < sc:
                    postincrement(swrap)[0] = aragorn.NOBASE

                swrap = wseq
                sc = swrap + d.psmax
                while swrap < sc:
                    postincrement(se)[0] = postincrement(swrap)[0]

                swrap = wseq
                sc = swrap + drewind
                while swrap < sc:
                    postincrement(se)[0] = postincrement(swrap)[0]

                sw.loffset = drewind
                sw.roffset = drewind
                start -= drewind
                flag = 1
                # goto SH
                loop = True
                SH = True
                NX = False

            else:
                swrap = wseq
                sc = seq + drewind
                while se < sc:
                    se[0] = cursor._forward(d) #aragorn.move_forward(d)
                    postincrement(swrap)[0] = postincrement(se)[0]

        # weird ass loop to emulate a GOTO
        while loop:

            # label NX: next
            sc = seq + aragorn.LSEQ
            if NX:
                while (se < sc):
                    postincrement(se)[0] = cursor._forward(d) #aragorn.move_forward(d)
                    if d.ps >= d.psmax:
                        if sw.linear:
                            for i in range(rewind):
                                postincrement(se)[0] = aragorn.NOBASE
                        else:
                            sc = wseq + drewind
                            swrap = wseq
                            while (swrap < sc):
                                postincrement(se)[0] = postincrement(swrap)[0]
                        flag = 1
                        SH = True
                        break

            # label SH: search
            if SH:
                length = <int> (se - seq)

                with gil:
                    PyErr_CheckSignals()

                # if (sw.verbose):
                #     vstart = sq(d, start + sw.loffset)
                #     vstop = sq(d, ((start + length) - sw.roffset) - 1)
                #     if (vstop < vstart):
                #         fprintf(stderr, "Searching from %ld to %ld\n", vstart, d.psmax)
                #         fprintf(stderr, "Searching from 1 to %ld\n", vstop)
                #     else:
                #         fprintf(stderr, "Searching from %ld to %ld\n", vstart, vstop)

                if (sw.both != 1):
                    sw.start = start
                    sw.comp = 0
                    nt = aragorn.tmioptimise(d, seq, length, nt, sw)

                if (sw.both > 0):
                    aragorn.sense_switch(seq, cseq, length)
                    sw.start = start + length
                    sw.comp = 1
                    nt = aragorn.tmioptimise(d, cseq, length, nt, sw)

                if not flag:
                    s = seq
                    sf = se - drewind
                    se = seq + drewind
                    while (s < se):
                        postincrement(s)[0] = postincrement(sf)[0]
                    start += length - drewind
                    # goto NX
                    NX = SH = loop = True
                    continue

                if nt < 1:
                    d.nf += 1
                if sw.maxintronlen > 0:
                    aragorn.remove_overlapping_trna(d, nt, sw)
                if sw.updatetmrnatags:
                    aragorn.update_tmrna_tag_database(sw.genes, nt, sw)

                # FIXME: here should sort genes and filter them with `gene_sort`
                # aragorn.batch_gene_set(d, nt, sw)

                # if sw.verbose:
                #     fprintf(stderr, "%s\nSearch Finished\n\n", d.seqname)

                d.ns += 1
                # exit loop
                loop = False

        return nt

    # if (d.ns > 1) and (sw.batch < 2):
    #     fprintf(f, ">end \t%d sequences", d.ns)
    #     if sw.trna or sw.mtrna:
    #         fprintf(f, " %d tRNA genes", sw.ngene[<int> aragorn.tRNA])
    #     if sw.tmrna:
    #         fprintf(f, " %d tmRNA genes", sw.ngene[<int> aragorn.tmRNA])
    #     if d.nf > 0:
    #         sens = (100.0 * (d.ns - d.nf)) / d.ns
    #         fprintf(f, ", nothing found in %d sequences, (%.2lf%% sensitivity)", d.nf, sens)
    #     fputc('\n', f)
    # if sw.updatetmrnatags:
    #     aragorn.report_new_tmrna_tags(sw)
