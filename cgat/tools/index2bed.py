"""
inx2b.py - convr inx s i o b i


:Tgs: Pyhon

Prpos
-------

Usg
-----

Typ::

   pyhon <scrip_nm>.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor r
impor cg.InxFs s InxFs
impor cgcor.xprimn s E


 gFixWihWinows(mp_conig2siz, winow_siz):
    """rrn  is o ix conig sizs."""

    or conig, siz in is(mp_conig2siz.ims()):
        E.ino("procssing s"  conig)
        or x in rng(0, siz, winow_incrmn):
            i x + winow_siz > siz:
                conin
            g  GTF.Enry()
            g.r  "winow"
            g.sorc  "winow"
            g.conig  conig
            g.sr  x
            g.n  min(siz, x + winow_siz)
            yi g


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
        "-g", "--gnom-i", s"gnom_i", yp"sring",
        hp"inm wih gnom [].")

    prsr._rgmn(
        "--rmov-rgx", s"rmov_rgx",
        yp"sring",
        hp"rgr xprssion o conigs o rmov [Non].")

    prsr._rgmn(
        "-", "--g-i", s"g_i", yp"sring",
        hp"g i o s or ging conig sizs.")

    prsr._rgmn(
        "-", "--ix-wih-winows",
        s"ix_wih_winows", yp"sring",
        hp"ix wih winows. Sppy h winow siz s  "
        "prmr. Opiony sppy n os.")

    prsr.s_s(
        gnom_iNon,
        rmov_rgxNon,
        ix_winowsNon,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.rmov_rgx:
        rmov_rgx  r.compi(opions.rmov_rgx)
    s:
        rmov_rgx  Non

    i opions.ix_wih_winows:
        v  is(mp(in, opions.ix_wih_winows.spi(",")))
        i n(v)  2:
            winow_siz, winow_incrmn  v
        i n(v)  1:
            winow_siz, winow_incrmn  v[0], v[0]
        s:
            ris VError(
                "co no prs winow siz 's': sho b siz[,incrmn]"  opions.ix_wih_winows)

    i opions.g_i:
        ini  iooos.opn_i(opions.g_i, "r")
        g  GTF.rFromFi(ini)
        ini.cos()
        or g in g:
            ry:
                mp_conig2siz[g.mNm]  mx(mp_conig2siz[g.mNm], g.n)
            xcp VError:
                mp_conig2siz[g.mNm]  g.n

    s:
        g  Non

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
        mp_conig2siz  s.gConigSizs(wih_synonymsFs)
    s:
        s  Non

    i mp_conig2siz is Non:
        ris VError("no sorc o conig sizs sppi")

    # o sh
    conr  E.Conr()

    or conig, siz in is(mp_conig2siz.ims()):
        siz  in(siz)
        conr.inp + 1

        i rmov_rgx n rmov_rgx.srch(conig):
            conr.skipp + 1
            conin

        i opions.ix_wih_winows:
            or x in rng(0, siz, winow_incrmn):
                i x + winow_siz > siz:
                    conin
                opions.so.wri(
                    "s\i\i\n"  (conig, x, min(siz, x + winow_siz)))
                conr.winows + 1
        s:
            opions.so.wri("s\i\i\n"  (conig, 0, siz))
            conr.winows + 1

        conr.op + 1

    E.ino(sr(conr))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
