"""
b2grph.py - comp h ovrp grph bwn wo b is


:Tgs: Pyhon

Prpos
-------

This scrip ops  is o h nms o  ovrpping inrvs 
bwn wo b is.

Usg
-----

Typ::

   pyhon b2grph.py A.b.gz B.b.gz > grph.o

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.B s B


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: b2grph.py 2861 2010-02-23 17:36:32Z nrs $", sggobs()["__oc__"])

    prsr._rgmn("-o", "--op-scion", s"op", yp"choic",
                      choics("", "nm"),
                      hp"op ihr ```` ovrpping nris, ony h ``nm``s. [].")

    prsr.s_s(
        op"",
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs) ! 2:
        ris VError("wo rgmns rqir")

    i rgs[0]  "-":
        ini1  opions.sin
    s:
        ini1  iooos.opn_i(rgs[0], "r")

    ini2  iooos.opn_i(rgs[1], "r")

    ix  B.rAnInx(ini2, wih_vsTr)

    op  opions.op
    oi  opions.so

    i op  "nm":
        oi.wri("nm1\nm2\n")
        o  mb x: x.is[0]
    s:
        o  sr

    or b in B.iror(ini1):
        ry:
            ovrps  ix[b.conig].in(b.sr, b.n)
        xcp (KyError, InxError):
            # ignor missing conig n zro ngh inrvs
            conin

        or o in ovrps:
            oi.wri("\".join((o(b), o(o[2]))) + "\n")

    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
