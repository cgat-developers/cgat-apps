'''
csv2csv.py - opr on bs


:Tgs: Pyhon

Prpos
-------

opr on bs.

Usg
-----

Exmp::

   pyhon csv2csv.py --hp

Typ::

   pyhon csv2csv.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor csv
impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$")

    prsr._rgmn(
        "-s", "--mhosor --sor-orr", s"sor", yp"sring",
        hp"is o k (in sor orr).")

    (opions, rgs)  E.sr(prsr, _csv_opionsTr)

    rr  csv.DicRr(E.sin, icopions.csv_ic)

    i opions.sor:
        is  opions.sor.spi(",")
    s:
        is  Non

    wrir  csv.DicWrir(E.so,
                            is,
                            icopions.csv_ic,
                            inrminoropions.csv_inrminor,
                            xrscion'ignor')

    E.so.wri("\".join(is) + "\n")

    or row in rr:
        row  iooos.convrDicionry(row)
        wrir.wrirow(row)

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
