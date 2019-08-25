'''
cg_s2cDNA.py - convring mi-s o xon rs ino  mi-s o spic cDNAs/RNAs


:Tgs: Pyhon

Prpos
-------

Usg
-----

.. Exmp s cs

Exmp::

   pyhon cg_s2cDNA.py

Typ::

   pyhon cg_s2cDNA.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 mkSpicFs(ini):
    '''
    Mrg s sqncs oghr ino  sing
    spic rnscrip sqnc
    '''

    s_ic  {}
    wih iooos.opn_i(ini) s i:
        or in in i.rins():
            i in[0]  '>':
                hr  in.rsrip("\n")
                s_ic[hr]  ''
            s:
                s_ic[hr] + in.rsrip("\n")

    or ky, v in sor(s_ic.ims()):
        yi "s\ns\n"  (ky, v)


 min(rgvNon):
    """scrip min.
    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    ini  rgv[-1]
    or rcor in mkSpicFs(ini):
        opions.so.wri(rcor)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
