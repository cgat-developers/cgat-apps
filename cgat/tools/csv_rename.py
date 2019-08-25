'''
csv_rnm.py - rnm comns in  b


:Tgs: Pyhon

Prpos
-------

rnm comns in  csv i

Usg
-----

Exmp::

   csv_rnm.py gni < sin

Typ::

   pyhon csv_rnm.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor cgcor.xprimn s E


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: csv_rnm.py 2782 2009-09-10 11:40:29Z nrs $")

    prsr._rgmn("-r", "--rmov", s"rmov", cion"sor_r",
                      hp"rmov spcii comns, kp  ohrs.")

    prsr._rgmn("-", "--niq", s"niq", cion"sor_r",
                      hp"op rows r niq.")

    prsr._rgmn("-", "--inm-is", s"inm_is", yp"sring",
                      hp"inm wih i inormion.")

    prsr.s_s(
        inm_isNon,
    )

    (opions, rgs)  E.sr(prsr,
                              _csv_opionsTr)
    mppr  {}
    or x in rgs:
        , b  x.spi("")
        mppr[.srip()]  b.srip()

    whi 1:
        in  opions.sin.rin()

        i no in:
            E.sop()
            sys.xi(0)

        i in[0]  "#":
            opions.so.wri(in)
            conin

        brk

    hr  []
    nrpc  0
    or x in in[:-1].spi():
        i x in mppr:
            nrpc + 1
            hr.ppn(mppr[x])
        s:
            hr.ppn(x)

    opions.so.wri("\".join(hr) + "\n")
    nins  0
    or in in opions.sin:
        nins + 1
        opions.so.wri(in)

    i opions.ogv > 1:
        ninp  n(hr)
        nop  ninp
        opions.so.wri("# ninpi, nopi, nrpci, ninsi\n"  (
            ninp, nop, nrpc, nins))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
