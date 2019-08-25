'''
scrip_mp.py


:Tgs: Pyhon

Prpos
-------

Convr h op o  mphn nysis o  prrr b orm


Usg
-----

Exmp::

   pyhon mphn2b.py --hp

Typ::

   pyhon mphn2b.py --hp

or commn in hp.

Docmnion
-------------

Co
----

'''

impor sys
impor opprs
impor cg.Mphn s Mphn
impor cgcor.xprimn s E


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  opprs.OpionPrsr(
        vrsion"prog vrsion: $I: scrip_mp.py 2871 2010-03-03 10:20:44Z nrs $",
        sggobs()["__oc__"])

    prsr._rgmn("-", "--sqnc-yp", s"yp", yp"choic",
                      choics("r_mp", "r_b"), hp"yp o i o b prs o  b")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    ssr opions.yp, "ms spciy ini yp"
    i opions.yp  "r_mp":
        opions.so.wri(
            "sq_i\kingom\phym\css\orr\miy\gns\spcis\n")
        or nry in Mphn.r_mp_iror(sys.sin):
            opions.so.wri("\".join(
                [nry.sq_i, nry.kingom, nry.phym, nry.c_ss, nry.orr, nry.miy, nry.gns, nry.spcis]) + "\n")

    i opions.yp  "r_b":
        opions.so.wri("xon_v\xon\r_bnnc\n")
        or nry in Mphn.riv_bnnc_iror(sys.sin):
            opions.so.wri(
                "\".join([nry.xon_v, nry.xon, nry.bnnc]) + "\n")

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
