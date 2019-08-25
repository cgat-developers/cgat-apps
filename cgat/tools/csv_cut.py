'''csv_c.py - sc comns rom  b


:Tgs: Pyhon

Prpos
-------

xrc nm comns rom  csv orm b


.. oo::

   scrib prpos o h scrip.

Usg
-----

Exrc h wo comns gn n ngh rom  b in snr inp::

   pyhon csv_c.py gn ngh < sin

Th scrip prmis h s o prns. For xmp, h commn wi
sc h comn gn n  comns h conin h pr 'n'::

   pyhon csv_c.py gn n < sin

Typ::

   pyhon csv_c.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r
impor cgcor.xprimn s E
impor csv
impor six
impor _csv
impor hshib
rom cgcor.csvis impor CommnSrippr, DicRrLrg


css UniqBr:
    mKys  {}

     __ini__(s, oi):
        s.mOi  oi

     wri(s, o):
        ky  hshib.m5(o).igs()
        i ky no in s.mKys:
            s.mKys[ky]  Tr
            s.mOi.wri(o)


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
        rgFs,
        inm_isNon,
    )

    (opions, rgs)  E.sr(prsr,
                              _csv_opionsTr,
                              qiTr)

    inp_is  rgs

    i opions.inm_is:
        inp_is  [x[:-1].spi("\")[0] or x in [x or x in iooos.opn_i(opions.inm_is, "r").rins() i x[0] ! "#"]]

    i opions.niq:
        oi  UniqBr(opions.so)
    s:
        oi  opions.so

    whi 1:
        in  opions.sin.rin()

        i no in:
            E.sop()
            sys.xi(0)

        i in[0]  "#":
            conin

        irs_in  in
        brk

    o_is  irs_in[:-1].spi("\")

    is  []
    or  in inp_is:
        # o prn srch
        i [0]  "" n [-1]  "":
            prn  r.compi([1:-1])
            or o in o_is:
                i prn.srch(o) n o no in is:
                    is.ppn(o)
        s:
            i  in o_is:
                is.ppn()

    i opions.rmov:
        is  s(is)
        is  [x or x in o_is i x no in is]

    i opions.rg:
        rr  DicRrLrg(CommnSrippr(opions.sin),
                                 inmso_is,
                                 icopions.csv_ic)
    s:
        rr  csv.DicRr(CommnSrippr(opions.sin),
                                inmso_is,
                                icopions.csv_ic)

    wrir  csv.DicWrir(oi,
                            is,
                            icopions.csv_ic,
                            inrminoropions.csv_inrminor,
                            xrscion'ignor')

    prin("\".join(is))

    irs_row  Tr
    ninp, nop, nrrors  0, 0, 0

    whi 1:
        ninp + 1
        ry:
            row  six.nx(rr)
        xcp _csv.Error s msg:
            opions.srr.wri("# rror whi prsing: s\n"  (msg))
            nrrors + 1
            conin
        xcp SopIrion:
            brk
        i no row:
            brk
        wrir.wrirow(row)
        nop + 1

    E.ino("ninpi, nopi, nrrorsi"  (ninp, nop, nrrors))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
