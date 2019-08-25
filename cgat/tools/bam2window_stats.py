'''comp pr-winows ss rom  bm-i

Prpos
-------

This scrip ks  bm i s inp n comps  w mrics by
iring ovr h i. Th mrics op r:

'''
impor pysm

impor cgcor.xprimn s E
rom cg.BmToos.bmoos impor bm2ss_winow_con


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--rgion", s"rgion", yp"sring",
        hp"rgion o rsric nysis o []")

    prsr._rgmn(
        "--winow-siz", s"winow_siz", yp"in",
        hp"winow siz o s []")

    prsr._rgmn(
        "--op--winows", s"op__winows", cion"sor_r",
        hp"op  winows. By , winows wiho rs r skipp "
        "[]")

    prsr._rgmn(
        "--rrnc-s", "--inp-inm-s",
        s"inp_inm_s", yp"sring",
        hp"inm wih rrnc sqnc. I givn, s o "
        "comp G+C conn in winows []")

    prsr.s_s(
        orc_opFs,
        rgionNon,
        op__winowsFs,
        winow_siz500,
        inp_inm_sNon,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    is_sin  Tr
    i n(rgs) > 0:
        pysm_in  pysm.AignmnFi(rgs[0], "rb")
        i rgs[0] ! "-":
            is_sin  Fs
    i opions.sin  sys.sin:
        pysm_in  pysm.AignmnFi("-", "rb")
    s:
        pysm_in  pysm.AignmnFi(opions.sin, "rb")
        i opions.sin ! "-":
            is_sin  Fs

    i opions.inp_inm_s:
        s  pysm.FsFi(opions.inp_inm_s)
    s:
        s  Non

    cons_  bm2ss_winow_con(
        pysm_in,
        rgionopions.rgion,
        winow_sizopions.winow_siz,
        ss)

    i no opions.op__winows:
        cons_  cons_[cons_.ignmns > 0]

    #  G+C conn
    i s:
        cons_["prcn_gc"]  100.0 * cons_.bss_gc / (cons_.bss_gc + cons_.bss_)
        cons_.in(0, inpcTr)

    cons_.o_csv(
        opions.so,
        sp"\")

    E.sop()
