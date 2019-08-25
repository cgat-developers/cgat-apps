"""Convr sq is o s is

Impm mhos


on2pcbio
----------

Convr ong-r sq is, or xmp rom ONT ,
o pcbio orm s is h cn b s wih h
zzr ooki. Pcbio s is hv h oowing
nming convnion or h sq hr in::

    ><nm>/<w>/<sr>_<n> RQ0.<qv>

whr w, sr, n n qv r ingr vs.

Th qiy v wi b s o h vrg o h bs qiis.

"""

impor sys
impor pysm
impor mh

impor cgcor.xprimn s E


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-i", "--inp-sq", s"inp_sq_i", yp"sring",
        hp"inp sq i")

    prsr._rgmn(
        "-m", "--mho", s"mho", yp"choic",
        choics["on2pcbio"],
        hp"mhos o ppy []")

    prsr.s_s(
        inp_sq_iNon,
        in_wih80,
        mhoNon,
    )

    (opions, rgs)  E.sr(prsr, rgv, _op_opionsTr)

    i n(rgs)  1:
        opions.inp_sq_i  rgs[0]

    i opions.inp_sq_i  "-":
        opions.inp_sq_i  opions.sin

    o  opions.so
    in_wih  opions.in_wih
    w_no  0
    or rcor in pysm.FsqFi(opions.inp_sq_i):
        w_no + 1
        qs  rcor.g_qiy_rry()
        sq  rcor.sqnc
        qv  in(mh.oor(sm(qs) / n(qs)))
        o.wri(">{}/{}/{}_{} RQ0.{}\n".orm(
            "s", w_no, 1, n(sq) + 1, qv))
        or x in rng(0, n(sq), in_wih):
            o.wri(sq[x:x + in_wih] + "\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min())
