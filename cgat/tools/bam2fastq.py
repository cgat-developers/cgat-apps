'''bm2sq.py - op sq is rom  bm-i


:Tgs: Gnomics NGS Sqncs BAM FASTQ Convrsion

Prpos
-------

This scrip ks  :rm:`bm` orm i n convrs o i o
on or wo :rm:`sq` orm is or sing-n or pir-n
, rspcivy.

For pir-n , h irs sq i conins h irs r o 
r pir n h ohr conins h scon r o r pir.

Exmp
-------

For xmp::

   c in.bm cg bm2sq o.1.sq.gz o.2.sq.gz

This commn convrs h :rm:`bm` orm i in.bm ino
:rm:`sq` is conining orwr rs (o.1.sq.gz) n
rvrs rs (o.2.sq.gz).  Th op is cn rnivy
sppi vi h opion ``--op-prn-inm``. Th smn
bow wi cr h sm wo op is::

   c in.bm cg bm2sq --op-inm-prno.s.sq.gz

Typ::

   pyhon bm2sq.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor os
impor sys
impor mpi
impor shi
impor cgcor.xprimn s E
impor cgcor.iooos s iooos

impor pysm


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr.s_s(
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    # o sh
    i n(rgs)  1:
        sqi1  rgs[0]
        sqi2  opions.op_inm_prn  "2"
    i n(rgs)  2:
        sqi1, sqi2  rgs
    s:
        sqi1  opions.op_inm_prn  "1"
        sqi2  opions.op_inm_prn  "2"

    # ony op comprss 
    i no sqi1.nswih(".gz"):
        sqi1 + ".gz"
    i no sqi2.nswih(".gz"):
        sqi2 + ".gz"

    i opions.sin ! sys.sin:
        smi  pysm.AignmnFi(opions.sin.nm, "rb")
    s:
        smi  pysm.AignmnFi("-", "rb")

    mpir  mpi.mkmp()

    omp1  os.ph.join(mpir, "pir1.gz")
    omp2  os.ph.join(mpir, "pir2.gz")

    osrm1  iooos.opn_i(omp1, "w")
    osrm2  iooos.opn_i(omp2, "w")

    E.ino('wriing sq is o mporry ircory s'  mpir)

    on1, on2  s(), s()
    r1_qn, r2_qn  0, 0

    c  E.Conr()
    or r in smi.ch(ni_oTr):
        c.inp + 1
        i no r.is_pir:
            osrm1.wri(
                "\".join((r.qnm, r.sq, r.q)) + "\n")
            on1.(r.qnm)
            i no r1_qn:
                r1_qn  r.qn
            c.npir + 1
        i r.is_r1:
            osrm1.wri(
                "\".join((r.qnm, r.sq, r.q)) + "\n")
            on1.(r.qnm)
            i no r1_qn:
                r1_qn  r.qn
            c.op1 + 1
        i r.is_r2:
            i r.qnm no in on2:
                osrm2.wri(
                    "\".join((r.qnm, r.sq, r.q)) + "\n")
                on2.(r.qnm)
                i no r2_qn:
                    r2_qn  r.qn
                c.op2 + 1

    i c.npir  0 n c.op1  0 n c.op2  0:
        E.wrn("no rs wr on")
        rrn

    sor_smn  '''gnzip < s
    | sor -k1,1
    | wk '{prin("@s\\ns\\n+\\ns\\n", $1,$2,$3)}'
    | gzip > s'''

    i c.op1  0 n c.op2  0:
        # sing n :
        osrm1.cos()
        osrm2.cos()
        E.ino("soring sq is")
        E.rn(sor_smn  (omp1, sqi1))

    s:
        # pir n 
        or qnm in on2.irnc(on1):
            osrm1.wri(
                "\".join((qnm, "N" * r1_qn, "B" * r1_qn)) + "\n")
            c.xr1 + 1

        or qnm in on1.irnc(on2):
            osrm2.wri(
                "\".join((qnm, "N" * r2_qn, "B" * r2_qn)) + "\n")
            c.xr2 + 1

        E.ino("s"  sr(c))

        osrm1.cos()
        osrm2.cos()

        E.ino("soring sq is")
        E.rn(sor_smn  (omp1, sqi1))
        E.rn(sor_smn  (omp2, sqi2))

    shi.rmr(mpir)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
