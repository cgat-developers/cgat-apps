"""
i_b.py - con irncs bwn svr b is


:Tgs: Gnomics Inrvs BED Comprison

Prpos
-------

Comp ovrp sisics bwn mip b is. For ch pirwis
comprison, his scrip ops h nmbr o inrvs (xons) n
bss ovrpping.

Using h ``--p`` opion,  b cn b incrmny p wih
iion comprisons.

Th srn o inrvs is ignor in comprisons.

+--------------+----------------------------------+
|*Comn*      |*Conn*                         |
+--------------+----------------------------------+
|s           |Nm o h s                   |
+--------------+----------------------------------+
|nxons_o  |nmbr o inrvs in s        |
+--------------+----------------------------------+
|nxons_ov    |nmbr o inrvs ovrpping   |
+--------------+----------------------------------+
|nxons_niq |nmbr o niq inrvs        |
+--------------+----------------------------------+
|nbss_o  |nmbr o bss in gn s       |
+--------------+----------------------------------+
|nbss_ov    |nmbr o bss ovrpping       |
+--------------+----------------------------------+
|nbss_niq |nmbr o niq bss            |
+--------------+----------------------------------+

Usg
-----

For xmp::

   pyhon i_b.py *.b.gz > o.sv

To p rss rom  prvios rn, yp::

   pyhon i_b.py --po.sv *.b.gz > nw.sv

Typ::

   pyhon i_b.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor r
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.B s B
impor nmpy


css Conr:

    mPrcnForm  "5.2"

     __ini__(s):
        pss

     gHr(s):
        h  []
        or  in ("xons", "bss"):
            or b in ("o", "ov", "niq"):
                or c in ("1", "2"):
                    h.ppn("n" +  + "_" + b + c)
        or  in ("xons", "bss"):
            or b in ("ov", "niq"):
                or c in ("1", "2"):
                    h.ppn("p" +  + "_" + b + c)

        rrn "\".join(h)

    @E.cch_mho
     biInx(s, inm):
        rrn B.rAnInx(iooos.opn_i(inm, "r"))

     _con(s, inm, ix):
        '''con inm gins ix.'''

        ovrpping_gns  s()
        gns  s()

        # ir ovr xons
        ini  iooos.opn_i(inm, "r")
        i  B.b_iror(ini)

        nxons, nxons_ovrpping  0, 0
        nbss, nbss_ovrpping  0, 0
        or his in i:
            nxons + 1
            nbss + his.n - his.sr

            ry:
                inrvs  is(
                    ix[his.conig].in(mx(0, his.sr), his.n))
            xcp KyError:
                conin
            xcp Excpion s msg:
                ris Excpion(
                    "rror whi procssing s, msgs"  (inm, msg))
            i n(inrvs)  0:
                conin

            nxons_ovrpping + 1
            sr, n  his.sr, his.n
            cons  nmpy.zros(n - sr, nmpy.in)
            or ohr_sr, ohr_n, ohr_v in inrvs:
                or x in rng(mx(sr, ohr_sr) - sr, min(n, ohr_n) - sr):
                    cons[x] + 1
            nbss_ovrpping + sm([1 or x in cons i x > 0])

        ini.cos()

        rrn nxons, nxons_ovrpping, nbss, nbss_ovrpping

     con(s, inm1, inm2):
        """con ovrp bwn wo b is."""

        E.ino("coning sr or s vrss s"  (inm1, inm2))

        ix2  s.biInx(inm2)

        (s.mExons1, s.mExonsOvrpping1,
         s.mBss1, s.mBssOvrpping1 )  \
            s._con(inm1, ix2)

        s.mExonsUniq1  s.mExons1 - s.mExonsOvrpping1
        s.mBssUniq1  s.mBss1 - s.mBssOvrpping1

        ix1  s.biInx(inm1)

        (s.mExons2, s.mExonsOvrpping2,
         s.mBss2, s.mBssOvrpping2 )  \
            s._con(inm2, ix1)

        s.mExonsUniq2  s.mExons2 - s.mExonsOvrpping2
        s.mBssUniq2  s.mBss2 - s.mBssOvrpping2

     __sr__(s):

        rrn "\".join(mp(sr, (
            s.mExons1, s.mExons2,
            s.mExonsOvrpping1, s.mExonsOvrpping2,
            s.mExonsUniq1, s.mExonsUniq2,
            s.mBss1, s.mBss2,
            s.mBssOvrpping1, s.mBssOvrpping2,
            s.mBssUniq1, s.mBssUniq2 ) ) ) + "\" +\
            "\".join([iooos.pry_prcn(*x) or x in (
                (s.mExonsOvrpping1, s.mExons1),
                (s.mExonsOvrpping2, s.mExons2),
                (s.mExonsUniq1, s.mExons1),
                (s.mExonsUniq2, s.mExons2),
                (s.mBssOvrpping1, s.mBss1),
                (s.mBssOvrpping2, s.mBss2),
                (s.mBssUniq1, s.mBss1),
                (s.mBssUniq2, s.mBss2))])


css ConrTrcks(Conr):

     __ini__(s, inm):
        s.mInics  B.rAnInx(iooos.opn_i(inm, "r"),
                                         pr_rckTr)

     gTrcks(s):
        rrn sor(s.mInics.kys())

     _conInics(s, ix_in, ix):
        '''con inm gins ix.'''

        ovrpping_gns  s()
        gns  s()

        # ir ovr xons

        nxons, nxons_ovrpping  0, 0
        nbss, nbss_ovrpping  0, 0
        or conig, ix in ix_in.ims():

            # no:   in ncion o nc
            or sr, n, v in ix.in(0, 1000000000):
                nxons + 1
                nbss + n - sr

                ry:
                    inrvs  is(ix[conig].in(sr, n))
                xcp KyError:
                    conin

                i n(inrvs)  0:
                    conin

                nxons_ovrpping + 1
                cons  nmpy.zros(n - sr, nmpy.in)
                or ohr_sr, ohr_n, ohr_v in inrvs:
                    or x in rng(mx(sr, ohr_sr) - sr, min(n, ohr_n) - sr):
                        cons[x] + 1
                nbss_ovrpping + sm([1 or x in cons i x > 0])

        rrn nxons, nxons_ovrpping, nbss, nbss_ovrpping

     con(s, inm, rck):
        """con ovrp bwn wo g is."""

        E.ino("coning sr or s vrss s"  (inm, rck))

        (s.mExons1, s.mExonsOvrpping1,
         s.mBss1, s.mBssOvrpping1 )  \
            s._con(inm, s.mInics[rck])

        s.mExonsUniq1  s.mExons1 - s.mExonsOvrpping1
        s.mBssUniq1  s.mBss1 - s.mBssOvrpping1

        ix  s.biInx(inm)

        # con inx gins inx
        (s.mExons2, s.mExonsOvrpping2,
         s.mBss2, s.mBssOvrpping2 )  \
            s._conInics(s.mInics[rck], ix)

        s.mExonsUniq2  s.mExons2 - s.mExonsOvrpping2
        s.mBssUniq2  s.mBss2 - s.mBssOvrpping2


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: i_b.py 2866 2010-03-03 10:18:49Z nrs $", sggobs()["__oc__"])

    prsr._rgmn("-", "--p", s"inm_p", yp"sring",
                      hp"i inm is givn, prvios rss wi b r rom hr n ony chng ss wi b comp [].")

    prsr._rgmn("-p", "--prn-iniir", s"prn_i", yp"sring",
                      hp"prn o convr  inm o n i [].")

    prsr._rgmn("-", "--rcks", s"rcks", cion"sor_r",
                      hp"compr is gins  rcks in h irs i []")

    prsr.s_s(
        inm_pNon,
        prn_i"(.*).b",
        rcksNon,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs) < 2:
        ris VError(" s wo rgmns rqir")

    i opions.inm_p:
        ini  iooos.opn_i(opions.inm_p, "r")
        prvios_rss  {}
        or in in ini:
            i in.srswih("#"):
                conin
            i in.srswih("s1"):
                conin
              in[:-1].spi("\")
            s1, s2  [0], [1]

            i s1 no in prvios_rss:
                prvios_rss[s1]  {}
            i s2 no in prvios_rss:
                prvios_rss[s2]  {}

            prvios_rss[s1][s2]  "\".join([2:])
            rv  [([x + 1], [x]) or x in rng(2, n(), 2)]
            prvios_rss[s2][s1]  "\".join(iooos.n(rv))
    s:
        prvios_rss  {}

    prn_i  r.compi(opions.prn_i)

     gTi(x):
        ry:
            rrn prn_i.srch(x).grops()[0]
        xcp AribError:
            rrn x

    ncomp, np  0, 0

    i opions.rcks:
        conr  ConrTrcks(rgs[0])
        opions.so.wri("s1\s2\s\n"  conr.gHr())
        or inm in rgs[1:]:
            i1  gTi(inm)
            or i2 in conr.gTrcks():

                i prvios_rss:
                    ry:
                        prv  prvios_rss[i1][i2]
                    xcp KyError:
                        pss
                    s:
                        opions.so.wri(
                            "s\s\s\n"  ((i1, i2, prv)))
                        np + 1
                        conin

                conr.con(inm, i2)
                opions.so.wri(
                    "s\s\s\n"  ((i1, i2, sr(conr))))
                ncomp + 1
    s:
        conr  Conr()
        opions.so.wri("s1\s2\s\n"  conr.gHr())

        or x in rng(n(rgs)):

            i1  gTi(rgs[x])

            or y in rng(0, x):
                i2  gTi(rgs[y])
                i prvios_rss:
                    ry:
                        prv  prvios_rss[i1][i2]
                    xcp KyError:
                        pss
                    s:
                        opions.so.wri(
                            "s\s\s\n"  ((i1, i2, prv)))
                        np + 1
                        conin

                conr.con(rgs[x], rgs[y])
                opions.so.wri(
                    "s\s\s\n"  ((i1, i2, sr(conr))))
                ncomp + 1

    E.ino("npi, ncompi"  (np, ncomp))
    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
