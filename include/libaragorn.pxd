from libc.stdio cimport FILE

cdef extern from "aragorn.h" nogil:

    const size_t STRLEN
    const size_t STRLENM1
    const size_t SHORTSTRLEN
    const size_t SHORTSTRLENM1

    const size_t NS

    const size_t MAXGCMOD
    const size_t MAMMAL_MT
    const size_t NGENECODE

    const size_t MININTRONLEN
    const size_t MAXINTRONLEN
    const size_t MINCTRNALEN
    const size_t MAXCTRNALEN
    const size_t MINTRNALEN
    const size_t MAXTRNALEN
    const size_t MAXETRNALEN

    const size_t NA
    const size_t ND
    const size_t NT
    const size_t NH
    const size_t NTH
    const size_t NC
    const size_t NGFT
    const size_t NTAG
    const size_t NTAGMAX

    cdef enum data_type:
        FASTA
        GENBANK

    cdef enum gencode:
        METAZOAN_MT
        STANDARD
        VERTEBRATE_MT

    cdef enum gene_type:
        noGENE
        tRNA
        tmRNA
        srpRNA
        rRNA
        CDS

    ctypedef struct annotated_gene:
        long int start
        long int stop
        int comp
        long int antistart
        long int antistop
        gene_type genetype
        bint pseudogene
        bint permuted
        bint detected
        char species[SHORTSTRLEN]

    ctypedef struct data_set:
        char filename[80]
        FILE *f
        char seqname[STRLEN]
        bint bugmode
        data_type datatype
        double gc
        long int filepointer
        long int ps
        long int psmax
        long int seqstart
        long int seqstartoff
        long int nextseq
        long int nextseqoff
        int ns
        int nf
        long int aseqlen
        int nagene[NS]
        annotated_gene gene[NGFT]

    ctypedef struct gene:
        char name[100]
        int seq[MAXTRNALEN + 1]
        int eseq[MAXETRNALEN + 1]
        int *ps
        int nbase
        int comp
        long int start
        long int stop
        int astem1
        int astem2
        int aatail
        int spacer1
        int spacer2
        int dstem
        int dloop
        int cstem
        int cloop
        int intron
        int nintron
        int anticodon
        int var
        int varbp
        int tstem
        int tloop
        int genetype
        double energy
        int asst
        int tps
        int tpe
        int annotation
        int annosc

    ctypedef struct csw:
        char genetypename[NS][10]
        FILE *f
        int batch
        int batchfullspecies
        int repeatsn
        int trna
        int tmrna
        int srprna
        int cds
        int mtrna
        int tvloop
        int cloop7
        int peptide
        gencode geneticcode
        int ngcmod
        int gcmod[MAXGCMOD]
        int gcfix
        int discrim
        int extastem
        int tarm
        int tagthresh
        int tarmlength
        int showconfig
        int libflag
        bint verbose
        int linear
        int both
        int reportpseudogenes
        int energydisp
        int secstructdisp
        int seqdisp
        int aataildisp
        int aataildiv
        int sp1max
        int sp2min
        int sp2max
        int mtxdetect
        int mtcdsscan
        int mtcompov
        int matchacceptor
        int maxintronlen
        int minintronlen
        int minintronlenreport
        int ioverlay
        int ifixedpos
        int ireportminintronlen
        int tmstrict
        int iamismatch
        int loffset
        int roffset
        long int start
        int comp
        gene* genes
        int genespace
        int srpspace
        int ngene[NS]
        int nps
        int annotated
        int dispmatch
        int updatetmrnatags
        int tagend
        int trnalenmisthresh
        int tmrnalenmisthresh
        int nagene[NS]
        int nafn[NS]
        int nafp[NS]
        int natfpd
        int natfptv
        int lacds
        int ldcds
        long int nabase
        double reportpsthresh
        double threshlevel
        double trnathresh
        double ttscanthresh
        double ttarmthresh
        double tdarmthresh
        double tastemthresh
        double tascanthresh
        double mttthresh
        double mtdthresh
        double mtdtthresh
        double mttarmthresh
        double mtdarmthresh
        double tmrnathresh
        double tmathresh
        double tmcthresh
        double tmcathresh
        double tmrthresh
        double srpthresh
        double cdsthresh
        double eref[NS]
        int tmrna_struct[200]


    void bopt_fastafile(data_set *d, csw *sw)