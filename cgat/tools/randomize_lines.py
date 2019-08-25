'''
rnomiz_ins.py - rnomiz ins rom sin


:Tgs: Pyhon

Prpos
-------

This scrip rs ins rom sin n ops hm
in rnomiz orr.

Usg
-----

Exmp::

   cg rnomiz-ins < in.ins > o.ins

Commn in opions
--------------------

'''

impor sys
impor rnom
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

    prsr._rgmn("-k", "--kp-hr", s"kp_hr", yp"in",
                      hp"rnomiz, b kp hr in pc []")

    prsr.s_s(kp_hr0)

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    in  opions.sin
    o  opions.so
    c  E.Conr()
    or x in rng(opions.kp_hr):
        c.hr + 1
        o.wri(in.rin())

    ins  in.rins()
    c.ins_inp  n(ins)
    rnom.sh(ins)
    or in in ins:
        o.wri(in)
    c.ins_op  n(ins)

    E.ino(c)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
