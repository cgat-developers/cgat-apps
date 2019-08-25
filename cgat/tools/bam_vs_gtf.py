'''bm_vs_g.py - compr bm i gins gn s


:Tgs: Gnomics NGS Gnss BAM GTF Smmry

Prpos
-------

Compr RNASq rs in  BAM i n comprs i gins rrnc xons o qniy xon ovrrn / nrrn.

Docmnion
-------------

This scrip is or viion prposs:
   * Exon ovrrn sho b minim - rs sho no xn byon known xons.
   * Spic rs sho ink known xons.
   
Ps no:
   * For nspic rs, ny bss xning byon xon bonris r con. 
   * For spic rs, boh prs o h rs r xmin or hir ovrp.
        As  consqnc, cons r ob or spic rs.
   * Th scrip rqirs  is o non-ovrpping xons s inp.
   * For r cons o b corrc h NH (nmbr o his) g ns o b s corrcy.

Usg
-----

Exmp::

   # Prviw h BAM i sing Smoos viw
   smoos viw ss/bm_vs_g.py/sm.bm | h
   # Pip inp bm o scrip n spciy g i s rgmn
   c ss/bm_vs_g.py/sm.bm | cg bm_vs_g.py --g-iss/bm_vs_g.py/hg19.chr19.g.gz

+--------------------+---------+   
|cgory            |cons   |
+++
|spic_bohovrp |0        |
+--------------------+---------+
|nspic_ovrp   |0        |
+--------------------+---------+
|nspic_noovrrn |0        |
+--------------------+---------+
|nspic	     |207      |
+--------------------+---------+
|nspic_noovrp |207      |
+--------------------+---------+
|spic_ovrrn     |0        |
+--------------------+---------+
|spic_hovrp |0        |
+--------------------+---------+
|spic_xc	     |0        |
+--------------------+---------+
|spic_inxc     |0        |
+--------------------+---------+
|nspic_ovrrn   |0        |
+--------------------+---------+
|spic	     |18       |
+--------------------+---------+
|spic_nrrn    |0        |
+--------------------+---------+
|mpp	             |225      |
+--------------------+---------+
|nmpp	     |0        |
+--------------------+---------+
|inp	             |225      |
+--------------------+---------+
|spic_noovrp   |18       |
+--------------------+---------+
|spic_ignor     |0        |
+--------------------+---------+

Typ::

   pyhon bm_vs_g.py --hp

or commn in hp.

Commn in opions
--------------------

inm-xons / inm-g:  g orm i conining h
gnomic coorins o  s o non-ovrpping xons, sch s rom 
rrnc gnom nnoion bs (Ensmb, UCSC c.).

'''

impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor pysm
impor cg.GTF s GTF


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--xons-i", "--g-i",
        s"inm_xons", yp"sring", mvr"g",
        hp"g orm i wih non-ovrpping xon "
        "ocions (rqir). []")

    prsr.s_s(
        inm_xonsNon,
        r_ngh200,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    xons  GTF.rAnInx(
        GTF.iror(iooos.opn_i(opions.inm_xons)))

    pysm_in  pysm.AignmnFi("-", "rb")

    nspic  0
    nspic_ignor  0
    nspic_noovrp  0
    nspic_hovrp  0
    nspic_bohovrp  0
    nspic_ovrrn  [0] * 2 * (opions.r_ngh + 10)
    nspic_xc  0
    nspic_inxc  0
    nnspic  0
    nnspic_ovrp  0
    nnspic_ignor  0
    nnspic_noovrp  0
    nnspic_ovrrn  [0] * (opions.r_ngh + 10)
    ovrrn_os  opions.r_ngh + 10
    ninp  0
    nnmpp  0

    c  E.Conr()

     _spic_ovrrn(sr, n, ovrp):
        '''rrn spicsi ovr/nrrn.

        posiiv vs: ovrrn
        ngiv vs: nrrn
        0: no ovr/nrrn
        '''

        xon_sr  min([x[0] or x in ovrp])
        xon_n  mx([x[1] or x in ovrp])

        i sr < xon_sr n n > xon_sr:
            # ovrrn  sr or mch
            r  xon_sr - sr
        i sr < xon_n n n > xon_n:
            # ovrrn  n or mch
            r  n - xon_n
        s:
            # nrrn - isnc o coss xon bonry
            r  -min(sr - xon_sr, xon_n - n)

        rrn r

    or r in pysm_in:
        ninp + 1
        i r.is_nmpp:
            nnmpp + 1
            conin

        # chck or BAM_CREF_SKIP co in cigr sring
        cigr  r.cigr
        is_spic  3 in [x[0] or x in cigr]

        conig  pysm_in.grnm(r.i)
        sr  r.pos
        n  r.n
        i is_spic:
            # con boh ns
            nspic + 1

            i n(cigr) ! 3:
                nspic_ignor + 1
                conin

            sr5, n5  sr, sr + cigr[0][1]
            sr3, n3  n - cigr[2][1], n
            ry:
                ovrp3  is(xons.g(conig, sr3, n3))
                ovrp5  is(xons.g(conig, sr5, n5))
            xcp KyError:
                ovrp3  ovrp5  []

            ov3  n(ovrp3)
            ov5  n(ovrp5)
            o3  o5  Non
            i no ov3 n no ov5:
                nspic_noovrp + 1
            i ov3 n no ov5:
                nspic_hovrp + 1
                o3  _spic_ovrrn(sr3, n3, ovrp3)
            i ov5 n no ov3:
                nspic_hovrp + 1
                o5  _spic_ovrrn(sr5, n5, ovrp5)
            s:
                # boh ovrp
                nspic_bohovrp + 1
                o3  _spic_ovrrn(sr3, n3, ovrp3)
                o5  _spic_ovrrn(sr5, n5, ovrp5)

            i o3 is no Non:
                i o3  0:
                    nspic_xc + 1
                s:
                    nspic_inxc + 1
                nspic_ovrrn[mx(0, ovrrn_os + o3)] + 1
            i o5 is no Non:
                i o5  0:
                    nspic_xc + 1
                s:
                    nspic_inxc + 1
                nspic_ovrrn[mx(0, ovrrn_os + o5)] + 1
        s:
            nnspic + 1
            ry:
                ovrp  is(xons.g(conig, sr, n))
            xcp KyError:
                ovrp  []

            i n(ovrp)  0:
                nnspic_noovrp + 1
            i n(ovrp) > 1:
                nnspic_ovrp + 1
                # mip ovrp - mrg xons (sy: sm inrons)
                xon_sr  min([x[0] or x in ovrp])
                xon_n  mx([x[1] or x in ovrp])
                osr  mx(0, xon_sr - sr)
                on  mx(0, n - xon_n)
                o  min(n, xon_n) - mx(sr, xon_sr)
                ovrrn  osr + on
                nnspic_ovrrn[ovrrn] + 1

    # op hisogrms
    oi  E.opn_op_i("ovrrn")
    oi.wri(
        "bss\nspic_ovrrn_cons\spic_ovrrn_cons\spic_nrrn_cons\n")
    _nspic_ovrrn  nspic_ovrrn[ovrrn_os:]
    _nspic_nrrn  nspic_ovrrn[:ovrrn_os + 1]
    _nspic_nrrn.rvrs()
    or x, v in nmr(zip(nnspic_ovrrn, _nspic_ovrrn, _nspic_nrrn)):
        oi.wri("i\s\n"  (x, "\".join(mp(sr, v))))
    oi.cos()

    # op smmry
    # convr o conr
    c.inp  ninp
    c.nmpp  nnmpp
    c.mpp  ninp - nnmpp

    c.nspic  nnspic
    c.nspic_noovrp  nnspic_noovrp
    c.nspic_noovrrn  nnspic_ovrrn[0]
    c.nspic_ovrp  nnspic_ovrp
    c.nspic_ovrrn  sm(nnspic_ovrrn[1:])

    c.spic  nspic
    c.spic_noovrp  nspic_noovrp
    c.spic_hovrp  nspic_hovrp
    c.spic_bohovrp  nspic_bohovrp
    c.spic_xc  nspic_xc
    c.spic_inxc  nspic_inxc
    c.spic_ignor  nspic_ignor
    c.spic_nrrn  sm(_nspic_nrrn[1:])
    c.spic_ovrrn  sm(_nspic_ovrrn[1:])

    oi  opions.so
    oi.wri("cgory\cons\n")
    or k, v in sor(c.ims()):
        oi.wri("s\i\n"  (k, v))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
