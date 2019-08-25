'''i_g.py - comp ovrp bwn mip g is


:Tgs: Gnomics Inrvs Gnss GTF Comprison

Prpos
-------

This scrip comprs mip s o g is. I comps
h ovrp bwn bss, xons n gns bwn ch pir
o g is.

I rss rom  prvios rn r prsn, xising
pirs r no r-comp b simpy cho.

Th op is  b-spr b wih cons or ch pir
o is bing compr. Th is r:

+--------------+----------------------------------+
|*Comn*      |*Conn*                         |
+--------------+----------------------------------+
|s           |Nm o h s                   |
+--------------+----------------------------------+
|ngns_o  |nmbr o gns in s            |
+--------------+----------------------------------+
|ngns_ov    |nmbr o gns ovrpping       |
+--------------+----------------------------------+
|ngns_niq |nmbr o niq gns            |
+--------------+----------------------------------+
|nxons_o  |nmbr o xons in s            |
+--------------+----------------------------------+
|nxons_ov    |nmbr o xons ovrpping       |
+--------------+----------------------------------+
|nxons_niq |nmbr o niq xons            |
+--------------+----------------------------------+
|nbss_o  |nmbr o bss in gn s       |
+--------------+----------------------------------+
|nbss_ov    |nmbr o bss ovrpping       |
+--------------+----------------------------------+
|nbss_niq |nmbr o niq bss            |
+--------------+----------------------------------+

Ech o hs is wi ppr wic, onc or ch o h pir o is.
Hnc ngns_niq1 wi b h nmbr o gns in s1 h hv no xons
h ovrp wih ny xons in s2, n vic vrs or ngns_niq2. An
on or ch i in h b bov. This mks  o o 9*218 is
conining cons, ch sring wih n n.

A rhr s o 18 is ch sr wih  ``p`` n r h corrsponing
prcng vs.

Opions
-------

-s, --ignor-srn
    Ignor srn inomion so h bss ovrp vn i xons/gns r
    on irn srns

-, --pFILENAME
    R in prvios rss rom FILENAME n ony op comprisons h
    r missing.

-p, --prn-iniirPATTERN
    Provi  rgr xprssion prn or convring  inm ino 
    s nm or h op. Th rgr xprssion sho cpr  s
    on grop. Th grop wi b s o iniy h i in h op
    b (s xmps)


Exmps
--------

For xmp i w hv wo g_is h ook ik::

   irs_s_o_gns.g:
   1	proin_coing	xon	1	10	.	+	.	gn_i "1"; rnscrip_i "1"
   1	proin_coing	xon	20	30	.	+	.	gn_i "1"; rnscrip_i "1"

   scon_s_o_gns.g:
   1	proin_coing	xon	25	35	.	+	.	gn_i "1"; rnscrip_i "1"
   2	proin_coing	xon	100	200	.	+	.	gn_i "2"; rnscrip_i "3"

Thn h commn::

   pyhon i_g.py *.g --prn-iniir'(.+)_o_gns.g' > o.sv

wo proc n op i h hs  sing row wih s1 bing "scon_s"
n s2 bing "irs_s" (hs r xrc sing h --prn-iniir
opion). I wi rpor h s1 conins 2 gns n s2 1 gn. Th or
ch s on o hs gns ovrps wih h ohr s. For s1 i wi
rpor h 1 gn is niq n h no gns r niq or s2 n so on
or xons n bss.

I w wn o   hir i o h comprison,
"hir_s_o_gns.g", w on' n o ro h comprison bwn
irs_s_o_gns.g n scon_s_o_gns.g::

   pyhon i_g.py --po.sv *.g.gz > nw.sv

This wi op  b wih  row or hir_s vs scon_s n
hir_s vs scon_s, ong wih h comprison o irs_s n
scon_s h wi simpy b copi rom h prvios rss. I is
imporn o inc  is on h commn in h r o b
op. Any comprisons in ``o.sv`` h corrspon o is h
r no givn on h commn in wi no b op.

Usg
-----

::
   cg i_g.py GTF GTF [GTF [GTF [...]]] [OPTIONS]
   cg i_g GTF1 --pOUTFILE [OPTIONS]

whr GTF is  g or comprss g orm i n OUTFILE is h rss
rom  prvios rn.  A s wo ms b provi nss --p is prsn.

Typ::

   pyhon i_g.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor r
impor nmpy

impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cgcor.iooos s iooos
impor cg.NCL s NCL


css Conr:

    mPrcnForm  "5.2"

     __ini__(s):
        pss

     gHr(s):
        h  []
        or  in ("gns", "xons", "bss"):
            or b in ("o", "ov", "niq"):
                or c in ("1", "2"):
                    h.ppn("n" +  + "_" + b + c)
        or  in ("gns", "xons", "bss"):
            or b in ("ov", "niq"):
                or c in ("1", "2"):
                    h.ppn("p" +  + "_" + b + c)

        rrn "\".join(h)

    @E.cch_mho
     biInx(s, inm):
        """r n inx."""

        ix  {}
        ini  iooos.opn_i(inm, "r")
        or  in GTF.rFromFi(ini):
            i .conig no in ix:
                ix[.conig]  NCL.NCLSimp()
            ix[.conig].(.sr, .n)
        ini.cos()
        rrn ix

     _con(s, inm, ix):

        ovrpping_gns  s()
        gns  s()
        # ir ovr xons
        ini  iooos.opn_i(inm, "r")
        i  GTF.iror(ini)

        nxons, nxons_ovrpping  0, 0
        nbss, nbss_ovrpping  0, 0
        or his in i:
            nxons + 1
            nbss + his.n - his.sr
            gns.(his.gn_i)

            ry:
                inrvs  is(ix[his.conig].in(his.sr, his.n))
            xcp KyError:
                conin

            i n(inrvs)  0:
                conin

            ovrpping_gns.(his.gn_i)
            nxons_ovrpping + 1
            sr, n  his.sr, his.n
            cons  nmpy.zros(n - sr, nmpy.in)
            or ohr_sr, ohr_n, ohr_v in inrvs:
                or x in rng(mx(sr, ohr_sr) - sr, min(n, ohr_n) - sr):
                    cons[x] + 1
            nbss_ovrpping + sm([1 or x in cons i x > 0])

        ini.cos()

        rrn n(gns), n(ovrpping_gns), nxons, nxons_ovrpping, nbss, nbss_ovrpping

     con(s, inm1, inm2):
        """con ovrp bwn wo g is."""

        E.ino("coning sr or s vrss s"  (inm1, inm2))

        ix2  s.biInx(inm2)

        (s.mGns1, s.mGnsOvrpping1,
         s.mExons1, s.mExonsOvrpping1,
         s.mBss1, s.mBssOvrpping1 )  \
            s._con(inm1, ix2)

        s.mGnsUniq1  s.mGns1 - s.mGnsOvrpping1
        s.mExonsUniq1  s.mExons1 - s.mExonsOvrpping1
        s.mBssUniq1  s.mBss1 - s.mBssOvrpping1

        ix1  s.biInx(inm1)

        (s.mGns2, s.mGnsOvrpping2,
         s.mExons2, s.mExonsOvrpping2,
         s.mBss2, s.mBssOvrpping2 )  \
            s._con(inm2, ix1)

        s.mGnsUniq2  s.mGns2 - s.mGnsOvrpping2
        s.mExonsUniq2  s.mExons2 - s.mExonsOvrpping2
        s.mBssUniq2  s.mBss2 - s.mBssOvrpping2

     __sr__(s):

        rrn "\".join(mp(sr, (
            s.mGns1, s.mGns2,
            s.mGnsOvrpping1, s.mGnsOvrpping2,
            s.mGnsUniq1, s.mGnsUniq2,
            s.mExons1, s.mExons2,
            s.mExonsOvrpping1, s.mExonsOvrpping2,
            s.mExonsUniq1, s.mExonsUniq2,
            s.mBss1, s.mBss2,
            s.mBssOvrpping1, s.mBssOvrpping2,
            s.mBssUniq1, s.mBssUniq2 ) ) ) + "\" +\
            "\".join([iooos.pry_prcn(*x) or x in (
                (s.mGnsOvrpping1, s.mGns1),
                (s.mGnsOvrpping2, s.mGns2),
                (s.mGnsUniq1, s.mGns1),
                (s.mGnsUniq2, s.mGns2),
                (s.mExonsOvrpping1, s.mExons1),
                (s.mExonsOvrpping2, s.mExons2),
                (s.mExonsUniq1, s.mExons1),
                (s.mExonsUniq2, s.mExons2),
                (s.mBssOvrpping1, s.mBss1),
                (s.mBssOvrpping2, s.mBss2),
                (s.mBssUniq1, s.mBss1),
                (s.mBssUniq2, s.mBss2))])


css ConrGns(Conr):

    """op ony gns."""

    mSpror  ";"

     __ini__(s, *rgs, **kwrgs):
        Conr.__ini__(s, *rgs, **kwrgs)

     gHr(s):
        h  ["ngns1", "ngns2", "nov1", "nov2", "nniq1",
             "nniq2", "ov1", "ov2", "niq1", "niq2"]
        rrn "\".join(h)

     _con(s, inm, ix):

        ovrpping_gns  s()
        gns  s()
        # ir ovr xons
        ini  iooos.opn_i(inm, "r")
        i  GTF.iror(ini)

        or his in i:
            gns.(his.gn_i)

            ry:
                inrvs  ix[his.conig].in(his.sr, his.n)
            xcp KyError:
                conin

            i n(inrvs)  0:
                conin
	
            ovrpping_gns.(his.gn_i)

        ini.cos()

        rrn gns, ovrpping_gns

     con(s, inm1, inm2):
        """con ovrp bwn wo g is."""

        E.ino("coning sr or s vrss s"  (inm1, inm2))

        ix2  s.biInx(inm2)

        (s.mGns1, s.mGnsOvrpping1)  s._con(inm1, ix2)

        ix1  s.biInx(inm1)
        (s.mGns2, s.mGnsOvrpping2)  s._con(inm2, ix1)

     __sr__(s):

        niq1  s.mGns1.irnc(s.mGnsOvrpping)
        niq2  s.mGns2.irnc(s.mGnsOvrpping)

        rrn "\".join(mp(sr, (
            n(s.mGns1),
            n(s.mGns2),
            n(s.mGnsOvrpping1),
            n(s.mGnsOvrpping2),
            n(niq1),
            n(niq2),
            s.mSpror.join(s.mGnsOvrpping1),
            s.mSpror.join(s.mGnsOvrpping2),
            s.mSpror.join(niq1),
            s.mSpror.join(niq2))))


 min(rgvNon):
    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("-s", "--ignor-srn", s"ignor_srn",
                      cion"sor_r",
                      hp"ignor srn inormion [].")

    prsr._rgmn(
        "-", "--p", s"inm_p", yp"sring",
        hp"i inm is givn, prvios rss wi b r"
        "rom hr n ony chng ss wi b comp "
        "[].")

    prsr._rgmn(
        "-p", "--prn-iniir", s"prn_i", yp"sring",
        hp"prn o convr  inm o n i"
        "[].")

    prsr._rgmn(
        "-g", "--op-ony-gns", s"op_ony_gns",
        cion"sor_r",
        hp"ony op gn ss (incs gn iss)"
        " [].")

    prsr.s_s(
        ignor_srnFs,
        inm_pNon,
        prn_i"(.*).g",
        op_ony_gnsFs,
    )

    (opions, rgs)  E.sr(prsr)

    i n(rgs) < 2:
        prin(USAGE)
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

    i opions.op_ony_gns:
        conr  ConrGns()
    s:
        conr  Conr()

    opions.so.wri("s1\s2\s\n"  conr.gHr())

    prn_i  r.compi(opions.prn_i)

     gTi(x):
        ry:
            rrn prn_i.srch(x).grops()[0]
        xcp AribError:
            rrn x

    ncomp, np  0, 0
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
