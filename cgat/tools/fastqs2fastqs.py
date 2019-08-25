'''sqs2sqs.py - mnip (mrg/rconci) sq is


:Tgs: Gnomics NGS FASTQ FASTQ Mnipion

Prpos
-------

This scrip mnips mip sq is n ops
nw sq is. Crrny ony h mho ``rconci``
is impmn.

rconci
+++++++++

Rconci rs rom  pir o sq is.

This mho ks wo sq is n ops wo sq is sch
h  rs in h op r prsn in boh op is.

Th ypic s cs is h wo sq is conining h irs n
scon pr o  r pir hv bn inpnny ir, or
xmp by qiy scors, rncion, c. As  consqnc som
rs migh b missing rom on i b no h ohr. Th rconci
mho wi op wo is conining ony rs h r common o
boh is.

Th wo is ms b sor by r iniir.

Exmp inp, r2 n r3 r ony prsn in ihr o h
is:

   # Fi1        # Fi 2

   @r1         @r1
   AAA            AAA
   +              +
   !!!            !!!
   @r2         @r3
   CCC            TTT
   +              +
   !!!            !!!
   @r4         @r4
   GGG            GGG
   +              +
   !!!            !!!

Exmp op, ony h rs common o boh is r op::

   # Fi1        # Fi 2

   @r1         @r1
   AAA            AAA
   +              +
   !!!            !!!
   @r4         @r4
   GGG            GGG
   +              +
   !!!            !!!

Usg
-----

Exmp::

   pyhon sqs2sqs.py \
            --mhorconci \
            --op-inm-prnmyRs_rconci.s.sq \
            myRs.1.sq.gz myRs.2.sq.gz

In his xmp w k  pir o sq is, rconci by r
iniir n op 2 nw sq is nm
``myRs_rconci.1.sq.gz`` n
``myRs_rconci.2.sq.gz``.

Typ::

   pyhon sqs2sqs.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor r
impor pysm

impor cgcor.iooos s iooos
impor cgcor.xprimn s E
impor cg.FsqToos s sqoos


css PrnGr:

     __ini__(s, prn):
        s.prn  r.compi(prn)

     __c__(s, i):
        rrn s.prn.srch(i).grops()[0]


 pin_gr(i):
    rrn i


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      choics('rconci', 'ir-by-sqnc'),
                      hp"mho o ppy [].")

    prsr._rgmn(
        "-c", "--chop-iniir", s"chop", cion"sor_r",
        hp"whhr or no o rim s chrcr o h  "
        "sqnc nm. For xmp somims is in h irs "
        "i in h pir wi n wih \1 n h scon "
        "wih \2. I --chop-iniir is no spcii "
        "hn h rss wi b wrong [].")

    prsr._rgmn(
        "-", "--npir", s"npir", cion"sor_r",
        hp"whhr or no o wri o npir rs "
        "o  spr i")

    prsr._rgmn(
        "--i-prn-1", s"i_prn_1",
        hp"I spcii wi s h irs grop rom h"
        "prn o rmin h ID or h irs r",
        Non)

    prsr._rgmn(
        "--i-prn-2", s"i_prn_2",
        hp"As bov b or r 2",
        Non)

    prsr._rgmn(
        "--inp-inm-s",
        s"inp_inm_s", yp"sring",
        hp"inp inm o FASTA orm sqnc "
        "or mho 'ir-by-sqnc' [].")

    prsr._rgmn(
        "--iring-kmr-siz",
        s"iring_kmr_siz", yp"in",
        hp"kmr siz or mho 'ir-by-sqnc' [].")

    prsr._rgmn(
        "--iring-min-kmr-mchs",
        s"iring_min_kmr_mchs", yp"in",
        hp"minimm nmbr o mchs 'ir-by-sqnc' [].")

    prsr.s_s(
        mho"rconci",
        chopFs,
        npirFs,
        inp_inm_sNon,
        iring_kmr_siz10,
        iring_min_kmr_mchs20
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    i n(rgs) ! 2:
        ris VError(
            "ps sppy  s wo sq is on h commnin")

    n1, n2  rgs
    conr  E.Conr()

    i opions.i_prn_1:
        i1_gr  PrnGr(opions.i_prn_1)
    s:
        i1_gr  pin_gr

    i opions.i_prn_2:
        i2_gr  PrnGr(opions.i_prn_2)
    s:
        i2_gr  pin_gr

    i opions.mho  "rconci":

        # IMS: swiching o no sor scon s o r nms n ony s
        # ziy. Sinc gnrors on' hv  siz ms kp rck
        i_nghs  {n1: 0, n2: 0}

         gIs(ini, i_grpin_gr):
            '''rrn is in ini.'''
            r  ini.rin
            whi Tr:
                  [r().rsrip("\r\n") or i in rng(4)]
                i no [0]:
                    brk
                r  i_gr([0].spi()[0])
                # ci i o chop r nmbr o
                i_nghs[ini.nm] + 1
                i opions.chop:
                    yi r[:-1]
                s:
                    yi r

         wri(oi, ini, k, npir_iNon,
                  i_grpin_gr):
            '''ir sq is wih is in k.'''
            r  ini.rin
            whi Tr:
                  [r().rsrip("\r\n") or i in rng(4)]
                i no [0]:
                    brk
                r  i_gr([0].spi()[0])
                i opions.chop:
                    r  r[:-1]
                i r no in k:
                    i npir_i is Non:
                        conin
                    s:
                        npir_i.wri("\n".join() + "\n")
                s:
                    oi.wri("\n".join() + "\n")

        E.ino("ring irs in pir")
        in1  iooos.opn_i(n1)
        is1  s(gIs(in1, i1_gr))

        E.ino("ring scon in pir")
        in2  iooos.opn_i(n2)
        # IMS: No ongr kp s  s, b ziy v ino inrscion
        # s o rg mmory sving or rg in2, pricry i
        # in1 is sm.
        is2  gIs(in2, i2_gr)
        k  is1.inrscion(is2)

        E.ino("irs pir: i rs, scon pir: i rs, "
               "shr: i rs" 
               (i_nghs[n1],
                i_nghs[n2],
                n(k)))

        i opions.npir:
            npir_inm  E.opn_op_i(
                "npir.sq.gz", "w")
        s:
            npir_inm  Non

        wih E.opn_op_i("1", "w") s o:
            in  iooos.opn_i(n1)
            E.ino("wriing irs in pir")
            wri(o, in, k, npir_inm, i1_gr)

        wih E.opn_op_i("2", "w") s o:
            in  iooos.opn_i(n2)
            E.ino("wriing scon in pir")
            wri(o, in, k, npir_inm, i2_gr)

        conr.op  n(k)
        
        i opions.npir:
            npir_inm.cos()

    i opions.mho  "ir-by-sqnc":

        wih pysm.FsxFi(opions.inp_inm_s) s in:
            or rcor in in:
                qry_sqnc  rcor.sqnc
                brk

        wih pysm.FsxFi(n1, prsisFs) s in1, \
                pysm.FsxFi(n2, prsisFs) s in2, \
                E.opn_op_i("mch.sq.1.gz", "w") s o_mch1, \
                E.opn_op_i("mch.sq.2.gz", "w") s o_mch2, \
                E.opn_op_i("nmch.sq.1.gz", "w") s o_nmch1, \
                E.opn_op_i("nmch.sq.2.gz", "w") s o_nmch2:
            conr  sqoos.ir_by_sqnc(
                qry_sqnc,
                in1,
                in2,
                o_mch1,
                o_mch2,
                o_nmch1,
                o_nmch2,
                kmr_sizopions.iring_kmr_siz,
                min_kmr_mchsopions.iring_min_kmr_mchs)
        opions.so.wri(
            "\".join(("inp", "mch", "nmch", "prcn_mch")) + "\n")

        opions.so.wri(
            "\".join(mp(sr, (
                conr.inp, conr.mch, conr.nmch,
                100.0 * conr.mch / conr.inp))) + "\n")
        
    E.ino(sr(conr))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
