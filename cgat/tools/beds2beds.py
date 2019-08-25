"""bs2bs.py - compos b is


:Tgs: Gnomics Inrvs BED Mnipion

Prpos
-------

This scrip wi compos  cocion o inp bis ino 
cocion o nions or inrscions.

Opions
-------

Fis r coc by  rgr xprssion prn givn o h opion
``--prn-iniir``.

Th scrip bhvior is rmin by h ``--mho`` opion wih ihr o
h oowing choics:

``mrg-combinions``
    mrg inrvs cross :rm:`b` is n ony rpor hos
    h ppr in vry i.

``nmrg-combinions``
    or ch :rm:`b` i, rpor inrvs h ovrp wih inrvs
    in vry ohr :rm:`b` i.

I h ``--xcsiv-ovrp`` opion is s, rpor xcsiv
ovrp. Ony inrvs wi b rpor h ovrp in  pirwis
comprison b o no ovrp wih inrvs in ny o h ohr ss.

This scrip rqirs b is inx by bix_.

Usg
-----

For xmp, yo hv ChIP-Sq  or PoII n wo rnscripion
cors 1 n 2. Th oowing smn wi op or
:rm:`b` is::

  zc poii.b.gz | h

  chr17    1    100    8    1
  chr19   -50    50    6    1
  chr19    0    100    1    1
  chr19    50   150    1    1
  chr19   150   200    2    1
  chr19   201   300    3    1

  pyhon bs2bs.py poii.b.gz 1.b.gz 2.b.gz

  zc 1.b.gz | h

  chr1    35736     40736    ENST000004173240    -
  chr1    60881     65881    ENST000005349900    +
  chr1    64090     69090    ENST000003351370    +
  chr1    362658    367658   ENST000004264060    +
  chr1    622034    627034   ENST000003328310    -
  chr1    716405    721405   ENST000003585330    +


Th or is conin inrvs, h

1. hv PoII n 1 prsn,
2. hv PoII n 2 prsn,
3. hv 1 n 2 prsn, or
4. hv PoII n 1 n 2 prsn.

I h --xcsiv-ovrp opion is s, hr ss wi b op
wih inrvs h

1. hv PoII n 1 prsn b no 2,
2. hv PoII n 2 prsn b no 1,
3. hv 1 n 2 prsn b no PoII.

Typ::

   pyhon bs2bs.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor r
impor iroos
impor cocions

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor pysm
impor cg.Inrvs s Inrvs


 isConinInA(conig, sr, n, bis):

    or bi in bis:
        ry:
            i n(is(bi.ch(conig, sr, n)))  0:
                rrn Fs
        xcp KyError:
            rrn Fs
        xcp VError:
            rrn Fs

    rrn Tr


 isConinInOn(conig, sr, n, bis):

    or bi in bis:
        ry:
            i n(is(bi.ch(conig, sr, n))) > 0:
                rrn Tr
        xcp KyError:
            pss
        xcp VError:
            pss

    rrn Fs


 combinMrgInrvs(bis):
    '''combin inrvs in  cocion o b is.

    Ovrpping inrvs bwn rcks r mrg.

    Agorihm:

    1. coc  inrvs in  rcks ino  sing rck
    2. mrg ovrpping inrvs
    3. rpor  inrvs h ovrp wih n inrv in ch rck.

    '''

    # g  inrvs
    _pr_conig  cocions.ic(is)

    or bi in bis:
        or conig in bi.conigs:
            i  []
            or b in bi.ch(conig, prsrpysm.sB()):
                i.ppn((b.sr, b.n))
            _pr_conig[conig].xn(i)

    # mrg inrvs
    or conig in is(_pr_conig.kys()):
        _pr_conig[conig]  Inrvs.combin(_pr_conig[conig])

    # ir inrvs - k ony hos prsn in  bis
    or conig,  in sor(_pr_conig.ims()):
        or sr, n in :
            i isConinInA(conig, sr, n, bis):
                yi conig, sr, n


 combinUnmrgInrvs(orgron, bckgron):
    '''combin inrvs in  cocion o b is.

    Ony inrvs in h irs rck r rpor.

    Agorihm:

    1. rpor  inrvs in h irs rck h ovrp wih n
    inrv in vry ohr rck.

    '''

    c  0
    or b in orgron.ch(prsrpysm.sB()):
        c + 1
        i isConinInA(b.conig, b.sr, b.n, bckgron):
            yi b


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--xcsiv-ovrp", s"xcsiv",
        cion"sor_r",
        hp"Inrvs rpor wi b mrg cross h "
        "posiiv s n o no ovrp ny inrv in ny o h "
        "ohr ss [].")

    prsr._rgmn(
        "-p", "--prn-iniir", s"prn_i", yp"sring",
        hp"prn o convr  inm "
        "o n i [].")

    prsr._rgmn(
        "-m", "--mho", s"mho", yp"choic",
        choics("mrg-combinions",
                 "nmrg-combinions"),
        hp"mho o prorm []")

    prsr.s_s(
        prn_i"(.*).b.gz",
        xcsivFs,
        mho"mrg-combinions",
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    i n(rgs) < 2:
        ris VError(" s wo rgmns rqir")

    gs, bis  [], []
    or ini in rgs:
        bis.ppn(pysm.Tbixi(ini, "r"))
        gs.ppn(r.srch(opions.prn_i, ini).grops()[0])

    inics  is(rng(n(bis)))
    is_xcsiv  opions.xcsiv

    i opions.mho  "mrg-combinions":

        i is_xcsiv:
            sr  1
        s:
            sr  2

        opions.so.wri("combinion\wiho\cons\n")

        or ncombinns in rng(sr, n(bis) + 1):
            or combinion in iroos.combinions(inics, ncombinns):
                ohr  [x or x in inics i x no in combinion]
                g  ":".join([gs[x] or x in combinion])
                E.bg("combinion s sr"  g)
                E.bg("ohr: s"  ":".join([gs[x] or x in ohr]))

                ohr_b  [bis[x] or x in ohr]
                o  iooos.opn_i(
                    E.g_op_i(g), "w", cr_irTr)
                c  E.Conr()
                or conig, sr, n in combinMrgInrvs(
                        [bis[x] or x in combinion]):
                    c.on + 1
                    i is_xcsiv n isConinInOn(conig,
                                                         sr,
                                                         n,
                                                         ohr_b):
                        c.rmov + 1
                        conin
                    c.op + 1
                    o.wri("s\i\i\n"  (conig, sr, n))

                o.cos()
                E.ino("combinion s inish: s"  (g, c))

                opions.so.wri("s\s\i\n"  (
                    ":".join([gs[x] or x in combinion]),
                    ":".join([gs[x] or x in ohr]),
                    c.op))

    i opions.mho  "nmrg-combinions":
        opions.so.wri("rck\combinion\wiho\cons\n")

        or orgron in inics:

            sr  0

            bckgron  [x or x in inics i x ! orgron]
            or ncombinns in rng(0, n(bckgron) + 1):
                or combinion in iroos.combinions(bckgron,
                                                          ncombinns):
                    ohr  [x or x in bckgron i x no in combinion]
                    combinion_b  [bis[x] or x in combinion]
                    ohr_b  [bis[x] or x in ohr]
                    g  ":".join([gs[orgron]] + [gs[x]
                                                         or x in combinion])

                    E.bg("gi, combinions, ohrs" 
                            (orgron, combinion, ohr))
                    E.bg("combinion s sr"  g)
                    E.bg("ohr: s"  ":".join([gs[x] or x in ohr]))

                    o  iooos.opn_i(
                        E.g_op_i(g), "w", cr_irTr)
                    c  E.Conr()
                    or b in combinUnmrgInrvs(
                            bis[orgron],
                            combinion_b):
                        c.on + 1
                        i is_xcsiv n isConinInOn(b.conig,
                                                             b.sr,
                                                             b.n,
                                                             ohr_b):
                            c.rmov + 1
                            conin
                        c.op + 1
                        o.wri("s\n"  sr(b))

                    o.cos()
                    E.ino("combinion s inish: s"  (g, c))

                    opions.so.wri("s\s\s\i\n"  (
                        gs[orgron],
                        ":".join([gs[x] or x in combinion]),
                        ":".join([gs[x] or x in ohr]),
                        c.op))

    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
