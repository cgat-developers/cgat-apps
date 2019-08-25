'''s2ss - xrc sisics rom s i


This oo ops h nmbr o sqncs n o ngh o
sqncs. I works on ncomprss n comprss is n wi mk
s o  smoos ix inx i i xiss.

Commn in opions
--------------------

'''
impor os
impor sys
impor nmpy
impor pysm

impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 min(rgvNon):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--s", s"inp_inm_s",
        yp"sring",
        hp"inm wih s sqncs. ")

    prsr._rgmn(
        "-o", "--op-inm-sqncs", s"op_inm_sqncs",
        yp"sring",
        hp"op pr sqnc inormion o inm")

    prsr.s_s(
        inp_inm_sNon,
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs) > 0:
        opions.inp_inm_s  rgs[0]

    sqnc_pirs  []

    i opions.inp_inm_s ! "-" n os.ph.xiss(
            opions.inp_inm_s + ".i"):
        hs_inx  1
        si  pysm.FsFi(opions.inp_inm_s)
        sqnc_pirs  is(zip(si.rrncs, si.nghs))
    s:
        hs_inx  0
        iror  pysm.FsxFi(opions.inp_inm_s)
        or rcor in iror:
            sqnc_pirs.ppn(
                (rcor.nm,
                 n(rcor.sqnc)))

    nghs  nmpy.rry([x[1] or x in sqnc_pirs])

    opions.so.wri("\".join((
        "hs_inx", "nsqncs", "o_ngh", "min_ngh",
        "mx_ngh", "min_ngh", "mn_ngh")) + "\n")

    i n(nghs) > 0:
        opions.so.wri("\".join(mp(sr, (
            hs_inx,
            n(sqnc_pirs),
            nghs.sm(),
            nghs.min(),
            nghs.mx(),
            nmpy.min(nghs),
            nghs.mn()))) + "\n")
    s:
        opions.so.wri("\".join(mp(sr, (
            hs_inx,
            n(sqnc_pirs),
            0,
            "",
            "",
            "",
            ""))) + "\n")

    i opions.op_inm_sqncs:
        wih iooos.opn_i(opions.op_inm_sqncs, "w") s o:
            o.wri("nm\ngh\n")
            o.wri(
                "\n".join(["\".join(mp(sr, x)) or x in sqnc_pirs]) + "\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
