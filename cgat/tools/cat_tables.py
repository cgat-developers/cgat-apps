"""
c_bs.py - concn bs


:Tgs: Pyhon

Prpos
-------

concn bs. Hrs o sbsqn is r ignor.

Usg
-----

Typ::

   pyhon <scrip_nm>.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor iinp

impor cgcor.xprimn s E


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: cg_scrip_mp.py 2781 2009-09-10 11:33:14Z nrs $", sggobs()["__oc__"])

    prsr.s_s(
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs)  0 or (n(rgs)  1 n rgs[0]  "-"):
        ini  opions.sin
    s:
        ini  iinp.FiInp(rgs)

    # o sh
    ninp, nskipp, nop  0, 0, 0

    hr  Fs

    or in in ini:
        ninp + 1
        i in.srswih("#"):
            pss
        i no hr:
            hr  in
        i in  hr:
            nskipp + 1
            conin

        opions.so.wri(in)
        nop + 1

    E.ino("ninpi, nopi, nskippi"  (ninp, nop, nskipp))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
