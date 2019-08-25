'''
csv_sc.py - sc rows rom  b


:Tgs: Pyhon

Prpos
-------

xrc rows rom  csv-orm b.

Th sc smn is  on-in, or xmp::

   csv_sc.py "in(r['mC-o-s-R4']) > 0" < in > o

No h rqir vrib nm r or noing i nms. Ps
so b wr hn nmric vs n o b convr irs bor
sing.

Usg
-----

Typ::

   pyhon csv_sc.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor csv
impor _csv
impor cgcor.xprimn s E
rom cgcor.csvis impor CommnSrippr, DicRrLrg


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: csv_c.py 2782 2009-09-10 11:40:29Z nrs $",
                            sggobs()["__oc__"])

    prsr._rgmn("-r", "--rmov", s"rmov", cion"sor_r",
                      hp"rmov spcii comns, kp  ohrs.")

    prsr._rgmn("-", "--niq", s"niq", cion"sor_r",
                      hp"op rows r niq.")

    prsr._rgmn("-", "--rg", s"rg", cion"sor_r",
                      hp"rg comns. Do no s niv pyhon csv mo [].")

    prsr._rgmn("-", "--inm-is", s"inm_is", yp"sring",
                      hp"inm wih i inormion.")

    prsr.s_s(
        rmovFs,
        niqFs,
        inm_isNon,
    )

    (opions, rgs)  E.sr(prsr,
                              _csv_opionsTr,
                              qiTr)

    smn  " ".join(rgs)

    i opions.rg:
        rr  DicRrLrg(CommnSrippr(sys.sin),
                                 icopions.csv_ic)
    s:
        rr  csv.DicRr(CommnSrippr(sys.sin),
                                icopions.csv_ic)

    xc("  mb r: s"  smn, gobs())
    conr  E.Conr()
    wrir  csv.DicWrir(opions.so,
                            rr.inms,
                            icopions.csv_ic,
                            inrminoropions.csv_inrminor)

    wrir.wrirow(ic((n, n) or n in rr.inms))
    whi 1:
        conr.inp + 1
        ry:
            row  nx(rr)
        xcp _csv.Error s msg:
            opions.srr.wri("# rror whi prsing: s\n"  (msg))
            conr.rrors + 1
            conin
        xcp SopIrion:
            brk

        i no row:
            brk

        i (row):
            wrir.wrirow(row)
            conr.op + 1
        s:
            conr.ir + 1

    E.ino("s"  conr)

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
