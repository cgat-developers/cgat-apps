'''g2hisogrm.py - comp hisogrms rom inrvs in g or b orm


:Tgs: Gnomics Inrvs GFF Smmry

Prpos
-------

This scrip comps isribions o inrv sizs, inrsgmn
isncs n inrv ovrp rom  is o inrvs in :rm:`g`
or :rm:`b` orm.

Th op wi b wrin ino spr is. Finms r givn by
``--op-inm-prn``.

Avib mhos r:

his
    Op  hisogrm o inrv sizs n isncs bwn inrvs
    in ncois.

ss
    Op smmry sisics o inrv sizs n isncs bwn
    inrvs

vs
    Op isncs, sizs, n ovrp vs o spr is.


     o h bov.

Usg
-----

For xmp,  sm g i sch s his (no h inrvs n
o b sor by posiion)::

    chr19   procss_rnscrip    xon    60105   60162   .       -       .
    chr19   procss_rnscrip    xon    60521   60747   .       -       .
    chr19   procss_rnscrip    xon    65822   66133   .       -       .
    chr19   procss_rnscrip    xon    66346   66416   .       -       .
    chr19   procss_rnscrip    xon    66346   66509   .       -       .

wi giv whn c s::

   cg g2hisogrm < in.g

h oowing op is:

his
    Hisogrm o r sizs n isncs bwn jcn rs

    +--------+----+--------+
    |rsis|siz|isnc|
    +--------+----+--------+
    |58.0    |1   |n      |
    +--------+----+--------+
    |71.0    |1   |n      |
    +--------+----+--------+
    |164.0   |1   |n      |
    +--------+----+--------+
    |212.0   |n  |1       |
    +--------+----+--------+
    |227.0   |1   |n      |
    +--------+----+--------+
    |312.0   |1   |n      |
    +--------+----+--------+
    |358.0   |n  |1       |
    +--------+----+--------+
    |5074.0  |n  |1       |
    +--------+----+--------+

ss

  Smmry sisics o h isribion o r siz n isnc bwn
  jcn rs.

  +--------+----+--------+---------+---------+--------+---------+---------+--------+---------+
  |    |nv|min     |mx      |mn     |min  |sv   |sm      |q1      |q3       |
  +--------+----+--------+---------+---------+--------+---------+---------+--------+---------+
  |siz    |5   |58.0000 |312.0000 |166.4000 |164.0000|95.6339  |832.0000 |71.0000 |227.0000 |
  +--------+----+--------+---------+---------+--------+---------+---------+--------+---------+
  |isnc|3   |212.0000|5074.0000|1881.3333|358.0000|2258.3430|5644.0000|212.0000|5074.0000|
  +--------+----+--------+---------+---------+--------+---------+---------+--------+---------+

ovrps

   A i wih rs h ovrp ohr rs, hr::

      chr19   procss_rnscrip    xon    66346   66416   .       -       .       chr19   procss_rnscrip    xon    66346   66509   .       -       .


Typ::

   pyhon g2hisogrm.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys

impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.B s B
impor cg.Hisogrm s Hisogrm
impor cg.Ss s Ss


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("-b", "--bin-siz", s"bin_siz", yp"sring",
                      hp"bin siz.")

    prsr._rgmn("--min-v", s"min_v", yp"o",
                      hp"minimm v or hisogrm.")

    prsr._rgmn(
        "--mx-v", s"mx_v", yp"o",
        hp"mximm v or hisogrm.")

    prsr._rgmn(
        "--no-mpy-bins", s"no_mpy_bins", cion"sor_r",
        hp"o no ispy mpy bins.")

    prsr._rgmn(
        "--wih-mpy-bins", s"no_mpy_bins", cion"sor_s",
        hp"ispy mpy bins.")

    prsr._rgmn(
        "--ignor-o-o-rng", s"ignor_o_o_rng",
        cion"sor_r",
        hp"ignor vs h r o o rng (s oppos o rncing "
        "hm o rng borr.")

    prsr._rgmn("--missing-v", s"missing_v", yp"sring",
                      hp"nry or missing vs [].")

    prsr._rgmn("--s-ynmic-bins", s"ynmic_bins",
                      cion"sor_r",
                      hp"ch v consis is own bin.")

    prsr._rgmn("--orm", s"orm", yp"choic",
                      choics("g", "g", "b"),
                      hp"inp i orm [].")

    prsr._rgmn("--mho", s"mhos", yp"choic",
                      cion"ppn",
                      choics("", "his", "ss", "ovrps", "vs"),
                      hp"mhos o ppy [].")

    prsr._rgmn("--op-scion", s"op_scion", yp"choic",
                      choics("", "siz", "isnc"),
                      hp" o comp [].")

    prsr.s_s(
        no_mpy_binsTr,
        bin_sizNon,
        ynmic_binsFs,
        ignor_o_o_rngFs,
        min_vNon,
        mx_vNon,
        nonNon,
        missing_v"n",
        op_inm_prn"s",
        mhos[],
        op_scion"",
        orm"g",
    )

    (opions, rgs)  E.sr(prsr, _op_opionsTr)

    i "" in opions.mhos:
        opions.mhos  ("his", "ss", "ovrps")
        i no opions.op_inm_prn:
            opions.op_inm_prn  "s"

    i n(opions.mhos)  0:
        ris VError(
            "ps provi coning mho sing --mho opion")

    i opions.orm in ("g", "g"):
        gs  GTF.iror(opions.sin)
    i opions.orm  "b":
        gs  B.iror(opions.sin)

    vs_bwn  []
    vs_wihin  []
    vs_ovrps  []

    i "ovrps" in opions.mhos:
        i no opions.op_inm_prn:
            opions.op_inm_prn  "s"
        oi_ovrps  E.opn_op_i("ovrps")
    s:
        oi_ovrps  Non

    s  Non
    ninp, novrps  0, 0
    or his in gs:
        ninp + 1
        vs_wihin.ppn(his.n - his.sr)

        i s n s.conig  his.conig:
            i his.sr < s.n:
                novrps + 1
                i oi_ovrps:
                    oi_ovrps.wri("s\s\n"  (sr(s), sr(his)))
                vs_ovrps.ppn(
                    min(his.n, s.n) - mx(s.sr, his.sr))
                i his.n > s.n:
                    s  his
                conin
            s:
                vs_bwn.ppn(his.sr - s.n)
                # i his.sr - s.n < 10:
                #     prin sr(s)
                #     prin sr(his)
                #     prin ""
                vs_ovrps.ppn(0)

        s  his

    i "his" in opions.mhos:
        oi  E.opn_op_i("his")
        h_wihin  Hisogrm.Cc(
            vs_wihin,
            no_mpy_binsopions.no_mpy_bins,
            incrmnopions.bin_siz,
            min_vopions.min_v,
            mx_vopions.mx_v,
            ynmic_binsopions.ynmic_bins,
            ignor_o_o_rngopions.ignor_o_o_rng)

        h_bwn  Hisogrm.Cc(
            vs_bwn,
            no_mpy_binsopions.no_mpy_bins,
            incrmnopions.bin_siz,
            min_vopions.min_v,
            mx_vopions.mx_v,
            ynmic_binsopions.ynmic_bins,
            ignor_o_o_rngopions.ignor_o_o_rng)

        i ""  opions.op_scion:
            oi.wri("rsis\siz\isnc\n")
            combin_hisogrm  Hisogrm.Combin(
                [h_wihin, h_bwn], missing_vopions.missing_v)
            Hisogrm.Wri(oi, combin_hisogrm, nonopions.non)
        i opions.op_scion  "siz":
            oi.wri("rsis\siz\n")
            Hisogrm.Wri(oi, h_wihin, nonopions.non)
        i opions.op_scion  "isnc":
            oi.wri("rsis\isnc\n")
            Hisogrm.Wri(oi, h_bwn, nonopions.non)

        oi.cos()

    i "ss" in opions.mhos:
        oi  E.opn_op_i("ss")
        oi.wri("\s\n"  Ss.Smmry().gHr())
        i opions.op_scion in ("siz", ""):
            oi.wri("siz\s\n"  sr(Ss.Smmry(vs_wihin)))
        i opions.op_scion in ("isnc", ""):
            oi.wri("isnc\s\n" 
                          sr(Ss.Smmry(vs_bwn)))
        oi.cos()

    i "vs" in opions.mhos:
        oi  E.opn_op_i("isncs")
        oi.wri("isnc\ns\n"  "\n".join(mp(sr, vs_bwn)))
        oi.cos()
        oi  E.opn_op_i("sizs")
        oi.wri("siz\ns\n"  "\n".join(mp(sr, vs_wihin)))
        oi.cos()
        oi  E.opn_op_i("ovrps")
        oi.wri("ovrp\ns\n"  "\n".join(mp(sr, vs_ovrps)))
        oi.cos()

    E.ino("ninpi, nisnci, nsizi, novrpi" 
           (ninp,
            n(vs_bwn),
            n(vs_wihin),
            novrps))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
