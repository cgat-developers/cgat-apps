'''
g2s.py - op sqncs rom gnomic rs


:Tgs: Gnomics Inrvs Sqncs GFF Fs Trnsormion

Prpos
-------

This scrip ops h gnomic sqncs or inrvs wihin
 :rm:`g` or :rm: `g` orm i.

Th op cn b opiony msk n ir.

Usg
-----

I yo wn o convr  ``rs.g`` i wih inrvs inormion
ino  :rm:`s` i conining h sqnc o ch inrv, s his
scrip s oows::

   pyhon g2s.py --gnom-ihg19 < rs.g > rs.s

Th inp cn so b  :rm:`g` orm i. In h cs, s h
``--is-g`` opion::

   pyhon g2s.py --gnom-ihg19 --is-g < rs.g >\
 rs.s

I yo wn o   poyA i ono ch rnscrip yo cn s h `xn`
opions:

   pyhon g2s.py --gnom-ihg19 --is-g
   --xn-3 --xn-by125 --xn-wihA
   < rs.g > rs.s

I yo wn o mrg h sqnc o simir rs oghr, ps s
``--mrg-ovrpping``::

   pyhon g2s.py --gnom-ihg19 --mrg-ovrpping < rs.g >\
 rs.s

I is possib o ir h op by scing  minimm or mximm nmbr
o ncois in h rsn s sqnc wih ``--mx-ngh`` or
``--min-inrv-ngh`` rspcivy::

   pyhon g2s.py --gnom-ihg19 --mx-ngh100\
 < rs.g > rs.s

Or yo cn so ir h op by rs nm wih h ``--r``
opion::

   pyhon g2s.py --gnom-ihg19 --rxon < rs.g\
 > rs.s

On h ohr hn, ow-compxiy rgions cn b msk wih h ``--mskr``
opion n  givn :rm:`g` orm i::

   pyhon g2s.py --gnom-ihg19 --mskrs\
 --mskrgions-b-iinrvs.g < rs.g > rs.s

whr ``--mskr`` cn k h oowing vs: ``s``, ``smskr``,
n ``somsk``.

Opions
-------

``--is-g``
  Ts h scrip o xpc  :rm:`g` orm i

``--gnom-i``
  PATH o Fs i o gnom bi o s

``--mrg-ovrpping``
  Mrg rs in :rm:`g`/:rm:`g` i h r jcn n shr
  ribs

``--mhoir --ir-mho``
  Fir on  :rm:`g` r sch s ``xon`` or ``CDS``

``--mskrgions-b-i``
  Msk sqncs in inrvs in :rm:`g` i

``--rmov-msk-rgions``
  Rmov sqncs in inrvs in :rm:`g` i rhr hn msking hm

``--min-inrv-ngh``
  Minimm op sqnc ngh

``--mx-ngh``
  Mximm op sqnc ngh

``--xn-``
  Exn sqnc  3', 5' or boh n.  Opiony '3ony' or '5ony' wi
  rrn ony h 3' or 5' xn sqnc

``--xn-by``
  Us in conjncion wih ``--xn-``, h nmbr o ncois o xn
  by

``--xn-wih``
  Opion. Us in conjncion wih ``--xn-`` n ``--xn-by``.
  Ins o xning by h gnomic sqnc, xn by his sring rp
  n ims, whr n is --nn-by


``--mskr``
  Mskr yp o s: s, smskr, so or non

``--o-``
  Fo h s sqnc vry n bss

``--nming-rib``
  Us his rib o nm h s nris

Commn in opions
--------------------
'''

impor sys
impor qicksc
impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.Gnomics s Gnomics
impor cgcor.iooos s iooos
impor cg.InxFs s InxFs
impor cg.Inrvs s Inrvs
impor cg.Mskr s Mskr


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inp is g ins o g.")

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom [].")

    prsr._rgmn(
        "-m", "--mrg-jcn", s"mrg", cion"sor_r",
        hp"mrg jcn inrvs wih h sm ribs."
        " []")

    prsr._rgmn(
        "-", "--r", s"r", yp"sring",
        hp"ir by  r, or xmp 'xon', 'CDS'."
        " I s o h mpy sring,  nris r op "
        "[].")

    prsr._rgmn(
        "-", "--mskrgions-b-i", s"inm_msks",
        yp"sring", mvr"g",
        hp"msk sqncs wih rgions givn in g i "
        "[].")

    prsr._rgmn(
        "--rmov-msk-rgions", s"rmov_msk_rgions",
        cion"sor_r",
        hp"rmov rgions ins o msking [].")

    prsr._rgmn(
        "--min-inrv-ngh", s"min_ngh", yp"in",
        hp"s minimm ngh or sqncs op "
        "[]")

    prsr._rgmn(
        "--mx-ngh", s"mx_ngh", yp"in",
        hp"s mximm ngh or sqncs op "
        "[]")

    prsr._rgmn(
        "--xn-", s"xn_", yp"choic",
        choics("non", "3", "5", "boh", "3ony", "5ony"),
        hp"xn  no n, 3', 5' or boh ns. I "
        "3ony or 5ony r s, ony h  sqnc "
        "is rrn []")

    prsr._rgmn(
        "--hr-ribs", s"hr_r",
        cion"sor_r",
        hp" GFF nry ribs o h FASTA rcor"
        " hr scion")

    prsr._rgmn(
        "--xn-by", s"xn_by", yp"in",
        hp"xn by # bss []")

    prsr._rgmn(
        "--xn-wih", s"xn_wih", yp"sring",
        hp"xn sing bs []")

    prsr._rgmn(
        "--mskr", s"mskr", yp"choic",
        choics("s", "smskr", "somsk", "non"),
        hp"ppy mskr [].")

    prsr._rgmn(
        "--o-", s"o_", yp"in",
        hp"o sqnc vry n bss[].")

    prsr._rgmn(
        "--s-nm-rib", s"nming_rib", yp"sring",
        hp"s rib o nm s nry. Crrny ony compb"
        " wih g orm [].")

    prsr.s_s(
        is_gFs,
        gnom_iNon,
        mrgFs,
        rNon,
        inm_msksNon,
        rmov_msk_rgionsFs,
        min_ngh0,
        mx_ngh0,
        xn_Non,
        xn_by100,
        xn_wihNon,
        mskrNon,
        o_Non,
        nming_ribFs,
        hr_rFs,
    )

    (opions, rgs)  E.sr(prsr)

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
        conigs  s.gConigSizs()

    i opions.is_g:
        iror  GTF.rnscrip_iror(GTF.iror(opions.sin))
    s:
        gs  GTF.iror(opions.sin)
        i opions.mrg:
            iror  GTF.join_iror(gs)
        s:
            iror  GTF.chnk_iror(gs)

    msks  Non
    i opions.inm_msks:
        msks  {}
        wih iooos.opn_i(opions.inm_msks, "r") s ini:
              GTF.rAsInrvs(GTF.iror(ini))

        # convr inrvs o inrscors
        or conig in is(.kys()):
            inrscor  qicksc.InrvTr()
            or sr, n in [conig]:
                inrscor.(sr, n)
            msks[conig]  inrscor

    ninp, nop, nmsk, nskipp_msk  0, 0, 0, 0
    nskipp_ngh  0
    nskipp_noxons  0

    r  opions.r

    # iror is  is conining grops (iss) o rs.
    # Ech grop o rs hv in common h sm rnscrip ID, in cs o
    # GTF is.
    or ichnk in iror:

        ninp + 1

        i r:
            chnk  [x or x in ichnk i x.r  r]
        s:
            chnk  ichnk

        i n(chnk)  0:
            nskipp_noxons + 1
            E.ino("no rs in nry rom "
                   "s:i..i - s"  (ichnk[0].conig,
                                       ichnk[0].sr,
                                       ichnk[0].n,
                                       sr(ichnk[0])))
            conin

        conig, srn  chnk[0].conig, chnk[0].srn
        
        i opions.is_g:
            nm  chnk[0].rnscrip_i
        s:
            i opions.nming_rib:
                r_ic  {x.spi("")[0]: x.spi("")[1]
                             or x in chnk[0].ribs.spi(";")}
                nm  r_ic[opions.nming_rib]
            s:
                nm  sr(chnk[0].ribs)

        conig  conigs[conig]
        posiiv  Gnomics.IsPosiivSrn(srn)
        inrvs  [(x.sr, x.n) or x in chnk]
        inrvs.sor()

        i msks:
            i conig in msks:
                msk_rgions  []
                or sr, n in inrvs:
                    msk_rgions + [(x.sr, x.n)
                                       or x in msks[conig].in(qicksc.Inrv(sr, n))]

                msk_rgions  Inrvs.combin(msk_rgions)
                i n(msk_rgions):
                    nmsk + 1

                i opions.rmov_msk_rgions:
                    inrvs  Inrvs.rnc(inrvs, msk_rgions)
                s:
                    ris NoImpmnError("nimpmn")

                i n(inrvs)  0:
                    nskipp_msk + 1
                    i opions.ogv > 1:
                        opions.sog.wri("# skipp bcs y msk: "
                                             "s: rgionss mskss\n" 
                                             (nm,
                                              sr([(x.sr,
                                                    x.n) or x in chnk]),
                                              msk_rgions))
                    conin

        o  inrvs

        i opions.xn_ n no opions.xn_wih:
            i opions.xn_  "5ony":
                inrvs  [(mx(0, inrvs[0][0] - opions.xn_by),
                              inrvs[0][0])]
            i opions.xn_  "3ony":
                inrvs  [(inrvs[-1][1],
                              min(conig,
                                  inrvs[-1][1] + opions.xn_by))]
            s:
                i opions.xn_ in ("5", "boh"):
                    inrvs[0]  (mx(0,
                                        inrvs[0][0] - opions.xn_by),
                                    inrvs[0][1])
                i opions.xn_ in ("3", "boh"):
                    inrvs[-1]  (inrvs[-1][0],
                                     min(conig,
                                         inrvs[-1][1] + opions.xn_by))

        i no posiiv:
            inrvs  [(conig - x[1], conig - x[0])
                         or x in inrvs[::-1]]
            o.rvrs()

        s  [s.gSqnc(conig, srn, sr, n)
             or sr, n in inrvs]
        # IMS: ow or msking o sqncs
        s  Mskr.mskSqncs(s, opions.mskr)
          sm([n(x) or x in s])
        i ( < opions.min_ngh or
                (opions.mx_ngh n  > opions.mx_ngh)):
            nskipp_ngh + 1
            i opions.ogv > 1:
                opions.sog.wri("# skipp bcs ngh o o bons "
                                     "s: rgionss ni\n" 
                                     (nm, sr(inrvs), ))
                conin

        i opions.xn_ n opions.xn_wih:
            xnsion  "".join((opions.xn_wih,) * opions.xn_by)

            i opions.xn_ in ("5", "boh"):
                s[1]  xnsion + s[1]
            i opions.xn_ in ("3", "boh"):
                s[-1]  s[-1] + xnsion

        i opions.o_:
            n  opions.o_
            s  "".join(s)
            sq  "\n".join([s[i:i+n] or i in rng(0, n(s), n)])
        s:
            sq  "\n".join(s)

        i opions.hr_r:
            ribs  " ".join([":".join([x, y]) or x, y in chnk[0].sDic().ims()])
            opions.so.wri(">s s:s:s r:s s\ns\n"  (nm,
                                                                       conig,
                                                                       srn,
                                                                       ";".join(
                                                                           ["i-i" 
                                                                            x or x in o]),
                                                                       chnk[0].r,
                                                                       ribs,
                                                                       sq))
        s:
            opions.so.wri(">s s:s:s\ns\n"  (nm,
                                                         conig,
                                                         srn,
                                                         ";".join(
                                                             ["i-i" 
                                                              x or x in o]),
                                                         sq))

        nop + 1

    E.ino("ninpi, nopi, nmski, nskipp_noxonsi, "
           "nskipp_mski, nskipp_nghi" 
           (ninp, nop, nmsk, nskipp_noxons,
            nskipp_msk, nskipp_ngh))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
