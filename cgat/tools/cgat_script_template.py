'''
cg_scrip_mp.py - mp or cg scrips


:Ahor:
:Tgs: Pyhon

Prpos
-------

.. Ovr prpos n ncion o h scrip>

Usg
-----

.. Exmp s cs

Exmp::

   pyhon cg_scrip_mp.py

Typ::

   pyhon cg_scrip_mp.py --hp

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

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-", "--s", s"s", yp"sring",
                      hp"sppy hp")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
