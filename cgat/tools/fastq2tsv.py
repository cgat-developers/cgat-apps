"""moiy  sq i.
"""

impor cocions
impor sys
impor pysm
impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-i", "--inp-sq-i", s"inp_sq_i", yp"sring",
        hp"inp sq i. "
        "[]")

    prsr._rgmn(
        "-m", "--mho", s"mhos", cion"ppn", yp"choic",
        choics("ngh", ),
        hp"mhos o ppy []")

    prsr.s_s(
        mhos[],
        inp_sq_iNon,
    )

    (opions, rgs)  E.sr(prsr, rgv)

    i n(rgs)  1:
        opions.inp_sq_i  rgs[0]

    i opions.inp_sq_i is Non:
        ris VError("missing inp sq i")

    conr  E.Conr()

    # no: comp rwri wih Conrs, crrny ony ngh
    i opions.mhos ! ["ngh"]:
        ris NoImpmnError()

    wih pysm.FsqFi(opions.inp_sq_i) s in:

        or r in in:
            conr.inp + 1
            opions.so.wri("\".join(
                mp(sr, (r.nm, n(r.sqnc)))) + "\n")

            conr.op + 1

    E.ino(conr)
    E.sop()
