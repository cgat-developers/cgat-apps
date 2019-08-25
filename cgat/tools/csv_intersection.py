'''
csv_inrscion.py - inrsc wo bs


:Tgs: Pyhon

Prpos
-------

.. oo::
   
   scrib prpos o h scrip.

Usg
-----

Exmp::

   pyhon csv_inrscion.py --hp

Typ::

   pyhon csv_inrscion.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
rom cgcor.csvis impor rTb
impor csv
impor hshib


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

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: csv_inrscion.py 2782 2009-09-10 11:40:29Z nrs $")

    prsr._rgmn("-", "--niq", s"niq", cion"sor_r",
                      hp"op rows r niq.")

    prsr.s_s(
        rmovFs,
        niqFs,
    )

    (opions, rgs)  E.sr(prsr, _csv_opionsTr)

    i n(rgs) ! 2:
        ris VError("ps spciy wo is o join")

    opions.inm1, opions.inm2  rgs

    b1  rTb(iooos.opn_i(opions.inm1, "r"))
    b2  rTb(iooos.opn_i(opions.inm2, "r"))

    i opions.niq:
        oi  UniqBr(sys.so)
    s:
        oi  opions.so

    # bi nw i is
    nw_is  []

    or x in opions.join_is1:
        nw_is.ppn(x)

    or x in is1:
        i x no in opions.join_is1:
            nw_is.ppn(x)
        i x no in opions.join_is2:
            nw_is.ppn(x)

        wrir  csv.DicWrir(oi,
                                is,
                                icopions.csv_ic,
                                inrminoropions.csv_inrminor,
                                xrscion'ignor')

    i n(ins) > 0:

        o_is  ins[0][:-1].spi("\")

        i opions.rmov:
            is  []
            or x in o_is:
                i x no in inp_is:
                    is.ppn(x)
        s:
            is  inp_is

        rr  csv.DicRr(ins,
                                icopions.csv_ic)

        prin("\".join(is))

        irs_row  Tr
        or row in rr:
            row  iooos.convrDicionry(row)
            wrir.wrirow(row)

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
