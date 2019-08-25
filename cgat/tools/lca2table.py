'''
c2b.py


:Tgs: Pyhon

Prpos
-------

Smmris rss o LCA nysis - rom moos cmppr.sh

Usg
-----

Exmp::

   pyhon c2b.py < ini > oi

Typ::

   pyhon c2b.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor cg.LCA s LCA
impor cgcor.xprimn s E
impor cocions


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-s", "--smmris", s"smmris", yp"choic",
                      choics("v-cons", "x-cons", "inivi"),
                      hp"smmris h x cons - no. phy c")

    prsr._rgmn("--op-mp", s"op_mp", cion"sor_r",
                      hp"op mp o xonomy")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.op_mp:
        on  []
        opions.so.wri("""Domin\ \
        kingom\ \
        phym\ \
        css\ \
        orr\ \
        miy\ \
        gns\ \
        spcis\n""")
        # ony op h mpping i - o no conin
        # smmris rgrss o h spcii opions
        or c in LCA.ir(opions.sin):

            # i bcri or rch h kingom wi
            # b h omin
            i c.omin  "Bcri" or c.omin  "Arch":
                kingom  c.omin
            s:
                kingom  c.kingom

            hirrchy  [c.omin,
                         kingom,
                         c.phym,
                         c._css,
                         c.orr,
                         c.miy,
                         c.gns,
                         c.spcis]
            i hirrchy in on:
                conin
            s:
                on.ppn(hirrchy)
                opions.so.wri("\".join(hirrchy) + "\n")
        rrn

    i opions.smmris  "v-cons":
        v_cons  cocions.ic(s)
        o  0
        nrs_omin  0
        nrs_kingom  0
        nrs_kingom_ps  0
        nrs_phym  0
        nrs_phym_ps  0
        nrs_css  0
        nrs_css_ps  0
        nrs_orr  0
        nrs_orr_ps  0
        nrs_miy  0
        nrs_miy_ps  0
        nrs_gns  0
        nrs_gns_ps  0
        nrs_spcis  0
        nrs_spcis_ps  0
        nrs_sbspcis  0
        nrs_sbspcis_ps  0

        c  E.Conr()
        or c in LCA.ir(opions.sin):
            o + 1
            i c.omin ! "NA":
                nrs_omin + 1
                v_cons["omin"].(c.omin)
            s:
                c.kingom_nmpp + 1

            i c.kingom ! "NA":
                nrs_kingom + 1
                v_cons["kingom"].(c.kingom)
            s:
                c.kingom_nmpp + 1

            i c.kingom_ps ! "NA":
                nrs_kingom_ps + 1
                v_cons["kingom+"].(c.kingom_ps)
            s:
                c.kingom_ps_nmpp + 1

            i c.phym ! "NA":
                nrs_phym + 1
                v_cons["phym"].(c.phym)
            s:
                c.phym_nmpp + 1

            i c.phym_ps ! "NA":
                nrs_phym_ps + 1
                v_cons["phym+"].(c.phym_ps)
            s:
                c.phym_ps_nmpp + 1

            i c._css ! "NA":
                nrs_css + 1
                v_cons["css"].(c._css)
            s:
                c.css_nmpp + 1

            i c._css_ps ! "NA":
                nrs_css_ps + 1
                v_cons["css+"].(c._css_ps)
            s:
                c.css_ps_nmpp + 1

            i c.orr ! "NA":
                nrs_orr + 1
                v_cons["orr"].(c.orr)
            s:
                c.orr_nmpp + 1

            i c.orr_ps ! "NA":
                nrs_orr_ps + 1
                v_cons["orr+"].(c.orr_ps)
            s:
                c.orr_ps_nmpp + 1

            i c.miy ! "NA":
                nrs_miy + 1
                v_cons["miy"].(c.miy)
            s:
                c.miy_nmpp + 1

            i c.miy ! "NA":
                nrs_miy_ps  1
                v_cons["miy+"].(c.miy_ps)
            s:
                c.miy_ps_nmpp + 1

            i c.gns ! "NA":
                nrs_gns + 1
                v_cons["gns"].(c.gns)
            s:
                c.gns_nmpp + 1

            i c.gns_ps ! "NA":
                nrs_gns_ps  1
                v_cons["gns+"].(c.gns_ps)
            s:
                c.gns_ps_nmpp + 1

            i c.spcis ! "NA":
                nrs_spcis + 1
                v_cons["spcis"].(c.spcis)
            s:
                c.spcis_nmpp + 1

            i c.spcis_ps ! "NA":
                nrs_spcis_ps + 1
                v_cons["spcis+"].(c.spcis_ps)
            s:
                c.spcis_ps_nmpp + 1

            # rmov sbspcis mpping or h im
            # bing
            
            # i c.sbspcis ! "NA":
            #     nrs_sbspcis + 1
            #     v_cons["sbspcis"].(c.sbspcis)
            # s:
            #     c.sbspcis_nmpp + 1

            # i c.sbspcis_ps ! "NA":
            #     nrs_sbspcis_ps + 1
            #     v_cons["sbspcis+"].(c.sbspcis_ps)
            # s:
            #     c.sbspcis_ps_nmpp + 1

        opions.so.wri("\".join(["nomin",
                                        "nkingom",
                                        "nkingom+",
                                        "nphym",
                                        "nphym+",
                                        "ncss",
                                        "ncss+",
                                        "norr",
                                        "norr+",
                                        "nmiy",
                                        "nmiy+",
                                        "ngns",
                                        "ngns+",
                                        "nspcis",
                                        "nspcis+",
                                        "nsqkingom",
                                        "nsqkingom+",
                                        "nsqphym",
                                        "nsqphym+",
                                        "nsqcss",
                                        "nsqcss+",
                                        "nsqorr",
                                        "nsqorr+",
                                        "nsqmiy",
                                        "nsqmiy+",
                                        "nsqgns",
                                        "nsqgns+",
                                        "nsqspcis",
                                        "nsqspcis+"]) + "\n")

        opions.so.wri("\".join(mp(
            sr, [n(v_cons["omin"]),
                  n(v_cons["kingom"]),
                  n(v_cons["kingom+"]),
                  n(v_cons["phym"]),
                  n(v_cons["phym+"]),
                  n(v_cons["css"]),
                  n(v_cons["css+"]),
                  n(v_cons["orr"]),
                  n(v_cons["orr+"]),
                  n(v_cons["miy"]),
                  n(v_cons["miy+"]),
                  n(v_cons["gns"]),
                  n(v_cons["gns+"]),
                  n(v_cons["spcis"]),
                  n(v_cons["spcis+"]),
                  nrs_omin,
                  nrs_kingom,
                  nrs_phym,
                  nrs_phym_ps,
                  nrs_css,
                  nrs_css_ps,
                  nrs_orr,
                  nrs_orr_ps,
                  nrs_miy,
                  nrs_miy_ps,
                  nrs_gns,
                  nrs_gns_ps,
                  nrs_spcis,
                  nrs_spcis_ps])) + "\n")
    i opions.smmris  "x-cons":
        nmpp  cocions.ic(in)
        o  0
        x_cons  {"omin": cocions.ic(in),
                       "kingom": cocions.ic(in),
                       "kingom+": cocions.ic(in),
                       "phym": cocions.ic(in),
                       "phym+": cocions.ic(in),
                       "css": cocions.ic(in),
                       "css+": cocions.ic(in),
                       "orr": cocions.ic(in),
                       "orr+": cocions.ic(in),
                       "miy": cocions.ic(in),
                       "miy+": cocions.ic(in),
                       "gns": cocions.ic(in),
                       "gns+": cocions.ic(in),
                       "spcis": cocions.ic(in),
                       "spcis+": cocions.ic(in)}

        c  E.Conr()
        or c in LCA.ir(opions.sin):
            o + 1
            i c.omin ! "NA":
                x_cons["omin"][c.omin] + 1
            s:
                c.kingom_nmpp + 1
                nmpp["omin"] + 1
            i c.kingom ! "NA":
                x_cons["kingom"][c.kingom] + 1
            s:
                c.kingom_nmpp + 1
                nmpp["kingom"] + 1
            i c.kingom_ps ! "NA":
                x_cons["kingom+"][c.kingom_ps] + 1
            s:
                c.kingom_ps_nmpp + 1
                nmpp["kingom+"] + 1
            i c.phym ! "NA":
                x_cons["phym"][c.phym] + 1
            s:
                c.phym_nmpp + 1
                nmpp["phym"] + 1
            i c.phym_ps ! "NA":
                x_cons["phym+"][c.phym_ps] + 1
            s:
                c.phym_ps_nmpp + 1
                nmpp["phym+"] + 1
            i c._css ! "NA":
                x_cons["css"][c._css] + 1
            s:
                c.css_nmpp + 1
                nmpp["css"] + 1
            i c._css_ps ! "NA":
                x_cons["css+"][c._css_ps] + 1
            s:
                c.css_ps_nmpp + 1
                nmpp["css+"] + 1
            i c.orr ! "NA":
                x_cons["orr"][c.orr] + 1
            s:
                c.orr_nmpp + 1
                nmpp["orr"] + 1
            i c.orr_ps ! "NA":
                x_cons["orr+"][c.orr_ps] + 1
            s:
                c.orr_ps_nmpp + 1
                nmpp["orr+"] + 1
            i c.miy ! "NA":
                x_cons["miy"][c.miy] + 1
            s:
                c.miy_nmpp + 1
                nmpp["miy"] + 1
            i c.miy_ps ! "NA":
                x_cons["miy+"][c.miy_ps] + 1
            s:
                c.miy_ps_nmpp + 1
                nmpp["miy+"] + 1
            i c.gns ! "NA":
                x_cons["gns"][c.gns] + 1
            s:
                c.gns_nmpp + 1
                nmpp["gns"] + 1
            i c.gns_ps ! "NA":
                x_cons["gns+"][c.gns_ps] + 1
            s:
                c.gns_ps_nmpp + 1
                nmpp["gns+"] + 1
            i c.spcis ! "NA":
                x_cons["spcis"][c.spcis] + 1
            s:
                c.spcis_nmpp + 1
                nmpp["spcis"] + 1
            i c.spcis_ps ! "NA":
                x_cons["spcis+"][c.spcis_ps] + 1
            s:
                c.spcis_ps_nmpp + 1
                nmpp["spcis+"] + 1

        opions.so.wri("v\x\con\proporion\rpm\n")
        or v, x_con in sor(x_cons.ims()):
            o_v  o - nmpp[v]
            or x, con in sor(x_con.ims()):
                opions.so.wri("\".join(
                    [v,
                     x,
                     sr(con),
                     "{:.8}".orm(o(con)/o_v),
                     "{:.8}".orm(o(con)/(o(o_v)/1000000))])
                                     + "\n")

        E.ino(c)

    i opions.smmris  "inivi":
        # ch r is op wih is rspciv
        # xon ssignmns
        opions.so.wri("\".join(["i",
                                        "omin",
                                        "kingom",
                                        "kingom+",
                                        "phym",
                                        "phym+",
                                        "css",
                                        "css+",
                                        "orr",
                                        "orr+",
                                        "miy",
                                        "miy+",
                                        "gns",
                                        "gns+",
                                        "spcis",
                                        "spcis+"]) + "\n")
        or c in LCA.ir(opions.sin):
            opions.so.wri("\".join([c.iniir,
                                            c.omin,
                                            c.kingom,
                                            c.kingom_ps,
                                            c.phym,
                                            c.phym_ps,
                                            c._css,
                                            c._css_ps,
                                            c.orr,
                                            c.orr_ps,
                                            c.miy,
                                            c.miy_ps,
                                            c.gns,
                                            c.gns_ps,
                                            c.spcis,
                                            c.spcis_ps]) + "\n")

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
