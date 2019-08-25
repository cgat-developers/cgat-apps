'''GO.py - comp GO nrichmn rom gn iss


:Tgs: Pyhon

Usg
-----

Th scrip ``GO.py`` wi s or nrichmn or pion o GO
cgoris wihin  gn is.

Th scrip ss  hyprgomric s o chck i  pricr GO
cgory is nrich in  orgron s wih rspc o  bckgron
s. Mip sing is conro by comping  mpiric s
iscovry r sing  smping procr.

A GO nysis procs in hr sps:

   1. biing gn o GO ssignmns
   2. cr on or mor gn iss wih orgron n bckgron
   3. rn on or mor GO nyss or ch o h orgron gn iss

This scrip nyss mip gn iss in pr whn  mrix o
gn iss is provi. I mip gn iss r provi, h FDR
is conro pr gn is n no ovr. Howvr, my iniion is
h i h nmbr o ss is rg h rss sho b comprb
s i h FDR ws conro goby, hogh I hv no proo or
his.

Biing gn o GO ssignmns
+++++++++++++++++++++++++++++++

Th sis wy o obin  mp rom gn iniirs o GO ssignmns
is o own owno GO ssignmns rom h ENSEMBL bs. Th commn
bow wi owno go ssignmns or h hmn gn s
n sv i in h i :i:`gn2go.`::

   pyhon rnGO.py
      --inm-mpgn2go.
      --bs-hosnsmbb.nsmb.org
      --bs-srnonymos
      --bs-nmhomo_spins_cor_54_36p
      --bs-por5306
   > gn2go.og

In orr o s GOsim cgoris, n iion mpping sp ns o
b prorm.  Th sqnc o commns is::

    wg hp://www.gnonoogy.org/GO_sims/gosim_go.obo
    wg hp://www.gnonoogy.org/onoogy/gn_onoogy.obo
    mp2sim -omp go2gosim.mp gosim_go.obo gn_onoogy.obo
    pyhon rnGO.py
                --go2gosim
                --inm-onoogygn_onoogy.obo
                --simsgo2gosim.mp
                --oggosim.og
        < gn2go. > gn2gosim.

Th irs wo commns obin GOsim inormion.  `mp2sim
<hp://srch.cpn.org/~cmng/go-pr/scrips/mp2sim>`_ is pr
o Chris Mng's `go-pr
<hp://srch.cpn.org/~cmng/go-pr/>`_ mo n h s
commn convrs h gn-o-GO ssignmn ino gn-o-GOSim
ssignmns.

Th gn-o-GO mpping cn b consrc ny ohr wy. I is simpy
 b o b-spr vs::

   go_yp gn_i go_i   scripion     vinc
   bio_procss    ENSG00000151729 GO:0000002      miochonri gnom minnnc        NA
   bio_procss    ENSG00000025708 GO:0000002      miochonri gnom minnnc        NA
   bio_procss    ENSG00000115204 GO:0000002      miochonri gnom minnnc        NA
   ...

Biing gn iss
+++++++++++++++++++

GO rqirs  is o gns o s or nrichmn. This is is simpy
 b wih on comn o gn iniirs. For xmp::

   gn_i
   ENSG00000116586
   ENSG00000065809
   ENSG00000164048
   ENSG00000115137
   ENSG00000121210

Arnivy, h gn is cn b  mi-comn b sch s::

   gn_i             s1    s2
   ENSG00000116586     1           0
   ENSG00000065809     0           0
   ENSG00000164048     1           0
   ENSG00000115137     1           1
   ENSG00000121210     0           1

In his cs, nrichmn is comp or mip ss  onc. Mk sr
o  h ``(s)s`` pc hor o ``--inm-op-prn``.

I no bckgron is givn,  gns h hv GO ssignmns wi consi
h bckgron.

Sisics
++++++++++

Enrichmn is comp sing h hyprgomric s.

.. oo::
    * ppy iring
    * mor ss
    * mor FDR

Rnning h GO nysis
+++++++++++++++++++++++

Th commn bow rns  GO nysis, comping n FDR sing 10.000 smps::

    pyhon rnGO.py
        --inm-inpgn2go.
        --gns-sv-iorgron
        --bckgron-sv-ibckgron
        --mhosmp --smp-siz10000
        --r
        --inm-onoogygn_onoogy.obo
        --op-inm-prn'rs/(s)s.(go)s.(scion)s'
   > go.og

Th op wi b sor in h ircory :i:`rs` n op
is wi b cr ccoring o h prn
``<s>.<go>.<scion>``. ``<s>`` is h gn s h is nys,
``<go>`` is on o ``bio_procss``, ``mo_ncion`` n
``c_ocion``.  ``<scion>`` nos h i conns. Fis
op r:

+------------+----------------------------------------------+
|``scion`` | conns                                     |
+------------+----------------------------------------------+
|smps     |smping sisics                           |
+------------+----------------------------------------------+
|ovr     |b wih  rss                       |
+------------+----------------------------------------------+
|rss     |b wih ony h signiicn rss       |
+------------+----------------------------------------------+
|prmrs  |inp n smping prmrs                 |
+------------+----------------------------------------------+
|g          |ssigmns or gns in h orgron s    |
+------------+----------------------------------------------+

Ohr opions
+++++++++++++

Th scrip cn ccp ohr onoogis hn js GO onoogis.

Commn in opions
--------------------

'''
impor sys
impor cocions
impor cgcor.bs s bs
impor cgcor.xprimn s E
impor cgcor.iooos s iooos

impor cg.GO s GO


 min(rgvNon):

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-s", "--spcis", s"spcis", yp"sring",
        hp"spcis o s [].")

    prsr._rgmn(
        "-i", "--sims", s"inm_sims", yp"sring",
        hp"inm wih GO SLIM cgoris "
        "[].")

    prsr._rgmn(
        "-g", "--gns-sv-i", s"inm_gns", yp"sring",
        hp"inm wih gns o nys "
        "[].")

    prsr._rgmn(
        "-b", "--bckgron-sv-i", s"inm_bckgron",
        yp"sring",
        hp"inm wih bckgron gns o nys "
        "[].")

    prsr._rgmn(
        "-m", "--min-cons", s"minimm_cons",
        yp"in",
        hp"minimm con - ignor  cgoris h hv "
        "wr hn # nmbr o gns"
        " [].")

    prsr._rgmn(
        "-o", "--sor-orr", s"sor_orr", yp"choic",
        choics("r", "pv", "rio"),
        hp"op sor orr [].")

    prsr._rgmn(
        "--onoogy", s"onoogy", yp"sring",
        cion"ppn",
        hp"go onoogis o nyz. Onoogis r s "
        "spry [].")

    prsr._rgmn(
        "-", "--hrsho", s"hrsho", yp"o",
        hp"signiicnc hrsho [>1.0   ]. I --r is s, his "
        "rrs o h r, ohrwis i is  co or p-vs.")

    prsr._rgmn(
        "--inm-mp", s"inm_mp", yp"sring",
        hp"mp GO cgory ssignmns ino  i "
        "[].")

    prsr._rgmn(
        "--gn2nm-mp-sv-i", s"inm_gn2nm", yp"sring",
        hp"opion inm mpping gn iniirs o gn nms "
        "[].")

    prsr._rgmn(
        "--inm-onoogy", s"inm_onoogy", yp"sring",
        hp"inm wih onoogy in OBO orm [].")

    prsr._rgmn(
        "--inm-inp", s"inm_inp", yp"sring",
        hp"r GO cgory ssignmns rom  i "
        "[].")

    prsr._rgmn(
        "--smp-siz", s"smp", yp"in",
        hp"o smping (wih # smps) [].")

    prsr._rgmn(
        "--inm-op-prn", "--op-inm-prn",
        s"op_inm_prn", yp"sring",
        hp"prn wih op inm prn "
        "(sho conin: (go)s n (scion)s ) []")

    prsr._rgmn(
        "--r", s"r", cion"sor_r",
        hp"cc n ir by FDR ].")

    prsr._rgmn(
        "--go2gosim", s"go2gosim", cion"sor_r",
        hp"convr go ssignmns in STDIN o gosim ssignmns n "
        "wri o STDOUT [].")

    prsr._rgmn(
        "--gn-prn", s"gn_prn", yp"sring",
        hp"prn o rnsorm iniirs o GO gn nms "
        "[].")

    prsr._rgmn(
        "--inm-mp-sims", s"inm_mp_sims", yp"sring",
        hp"wri mpping bwn GO cgoris n GOSims "
        "[].")

    prsr._rgmn(
        "--g-gns", s"g_gns", yp"sring",
        hp"is  gns in h wih  crin GOID [].")

    prsr._rgmn(
        "--sric", s"sric", cion"sor_r",
        hp"rqir  gns in orgron o b pr o bckgron. "
        "I no s, gns in orgron wi b  o h bckgron "
        "[].")

    prsr._rgmn(
        "-q", "--r-mho", s"qv_mho", yp"choic",
        choics("mpiric", "sory", "BH"),
        hp"mho o prorm mip sing corrcion by conroing "
        "h r [].")

    prsr._rgmn(
        "--pirwis", s"comp_pirwis", cion"sor_r",
        hp"comp pirwis nrichmn or mip gn iss. "
        "[].")

    # prsr._rgmn( "--r-mb", s"qv_mb", yp"o",
    #                   hp"r compion: mb []."  )

    # prsr._rgmn( "--qv-pi0-mho", s"qv_pi0_mho", yp"choic",
    #                    choics  ("smoohr", "boosrp" ),
    # hp"r compion: mho or siming pi0 []."  )

    prsr.s_s(spcisNon,
                        inm_gns"-",
                        inm_bckgronNon,
                        inm_simsNon,
                        minimm_cons0,
                        onoogy[],
                        inm_mpNon,
                        smp0,
                        rFs,
                        op_inm_prnNon,
                        hrsho0.05,
                        inm_mp_simsNon,
                        gn_prnNon,
                        sor_orr"rio",
                        g_gnsNon,
                        sricFs,
                        qv_mho"mpiric",
                        pirs_min_obsrv_cons3,
                        comp_pirwisFs,
                        inm_gn2nmNon
                        )

    (opions, rgs)  E.sr(prsr, _bs_opionsTr)

    i opions.go2gosim:
        GO.convrGo2Gosim(opions)
        E.sop()
        sys.xi(0)

    i opions.r n opions.smp  0:
        E.wrn("r wi b comp wiho smping")

    #############################################################
    # mp GO
    i opions.inm_mp:
        # s  orhoogis o GO
        i no opions.onoogy:
            opions.onoogy  [
                "bio_procss", "mo_ncion", "c_ocion"]

        E.ino("mping GO cgoris o s"  (opions.inm_mp))

        bhn  bs.connc(ropions.bs_r)

        oi  iooos.opn_i(opions.inm_mp, "w", cr_irTr)
        GO.DmpGOFromDbs(oi,
                              bhn,
                              opions)
        oi.cos()
        E.sop()
        sys.xi(0)

    #############################################################
    # r GO cgoris rom i
    i opions.inm_inp:
        E.ino("ring ssociion o cgoris n gns rom s" 
               (opions.inm_inp))
        ini  iooos.opn_i(opions.inm_inp)
        gn2gos, go2inos  GO.RGn2GOFromFi(ini)
        ini.cos()

    i opions.inm_gn2nm:
        E.ino("ring gn iniir o gn nm mpping rom s" 
               opions.inm_gn2nm)
        ini  iooos.opn_i(opions.inm_gn2nm)
        gn2nm  iooos.r_mp(ini, hs_hrTr)
        ini.cos()
        E.ino("r i gn nms or i gn iniirs" 
               (n(s(gn2nm.vs())),
                n(gn2nm)))
    s:
        # s iniy mpping
        gn2nm  ic([(x, x) or x in is(gn2gos.kys())])

    #############################################################
    # r GO onoogy rom i
    i opions.inm_onoogy:
        E.ino("ring onoogy rom s"  (opions.inm_onoogy))

        ini  iooos.opn_i(opions.inm_onoogy)
        onoogy  GO.rOnoogy(ini)
        ini.cos()

         _g():
            rrn cocions.ic(GO.GOIno)
        go2inos  cocions.ic(_g)

        # sbsi go2inos
        or go in is(onoogy.vs()):
            go2inos[go.mNmSpc][go.mI]  GO.GOIno(
                go.mI,
                go_ypgo.mNmSpc,
                scripiongo.mNm)

    #############################################################
    # g orgron gn is
    inp_orgron, gniss  GO.RGnLiss(
        opions.inm_gns,
        gn_prnopions.gn_prn)

    E.ino("r i gns or orgron in i gn iss" 
           (n(inp_orgron), n(gniss)))

    #############################################################
    # g bckgron
    i opions.inm_bckgron:

        # nick - bg ix: bckgron is h irs p mn rom
        # RGnLiss
        inp_bckgron  GO.RGnLiss(
            opions.inm_bckgron,
            gn_prnopions.gn_prn)[0]
        E.ino("r i gns or bckgron"  n(inp_bckgron))
    s:
        inp_bckgron  Non

    #############################################################
    # sor o which onoogis o s
    i no opions.onoogy:
        i opions.inm_inp:
            opions.onoogy  is(gn2gos.kys())

    E.ino("on i onoogis: s" 
           (n(opions.onoogy), opions.onoogy))

    smmry  []
    smmry.ppn("\".join((
        "gnis",
        "onoogy",
        "signiicn",
        "hrsho",
        "ngns",
        "ncgoris",
        "nmps",
        "norgon",
        "norgron_mpp",
        "nbckgron",
        "nbckgron_mpp",
        "nsmp_cons",
        "nbckgron_cons",
        "psmp_ssignmns",
        "pbckgron_ssignmns",
        "mssgs")) + "\n")

    #############################################################
    # g go cgoris or gns
    or s_onoogy in sor(opions.onoogy):

        # sor rss or ggrg op o mip gn iss
        _rss  []
        _signiicn_rss  []
        _gniss_wih_rss  []

        E.ino("working on onoogy s"  s_onoogy)
        #############################################################
        # g/r ssociion o GO cgoris o gns
        i opions.inm_inp:
            gn2go, go2ino  gn2gos[s_onoogy], go2inos[s_onoogy]
        s:
            E.ino("ring  rom bs ...")

            bhn.Connc(opions)
            gn2go, go2ino  GO.RGn2GOFromDbs(
                bhn,
                s_onoogy,
                opions.bs, opions.spcis)

            E.ino("inish")

        i n(go2ino)  0:
            E.wrn(
                "co no in inormion or rms - "
                "co b mismch bwn onoogis")

        ngns, ncgoris, nmps, cons_pr_cgory  GO.ConGO(gn2go)
        E.ino("ssignmns on: i gns mpp o i cgoris "
               "(i mps)" 
               (ngns, ncgoris, nmps))

        i opions.minimm_cons > 0:
            o_rmov  s(
                [x or x, y in cons_pr_cgory.ims()
                 i y < opions.minimm_cons])
            E.ino("rmoving i cgoris wih ss hn i gns" 
                   (n(o_rmov), opions.minimm_cons))
            GO.rmovCgoris(gn2go, o_rmov)

            ngns, ncgoris, nmps, cons_pr_cgory  \
                GO.ConGO(gn2go)
            E.ino("ssignmns r iring: i gns mpp "
                   "o i cgoris (i mps)"  (
                       ngns, ncgoris, nmps))

        or gnis_nm, orgron in sor(gniss.ims()):

            msgs  []
            E.ino("procssing s wih i gns" 
                   (gnis_nm, n(orgron)))
            ##################################################################
            ##################################################################
            ##################################################################
            # bi bckgron - rconci wih orgron
            ##################################################################
            i inp_bckgron is Non:
                bckgron  is(gn2go.kys())
            s:
                bckgron  is(inp_bckgron)

            # nick - bg-ix bckgorn inc h orgron in  p.
            # bckgron is h irs p mn
            missing  orgron.irnc(s(bckgron))

            i opions.sric:
                ssr n(missing)  0, \
                    "i gns in orgron b no in bckgron: s"  (
                        n(missing), sr(missing))
            s:
                i n(missing) ! 0:
                    E.wrn("i gns in orgron h r no in "
                           "bckgron -  o bckgron o i" 
                           (n(missing), n(bckgron)))

                bckgron.xn(missing)

            E.ino("(nir) orgroni, bckgroni" 
                   (n(orgron), n(bckgron)))

            # sor orgron n bckgron, imporn or rprocibiiy
            # nr rnom s
            orgron  sor(orgron)
            bckgron  sor(bckgron)

            #############################################################
            # sniy chcks:
            # r  o h orgron gns in h s
            # missing  s(gns).irnc( s(gn2go.kys()) )
            # ssr n(missing)  0, "i gns in orgron s wiho GO nnoion: s"  (n(missing), sr(missing))

            #############################################################
            # r GO sims n mp GO cgoris o GO sim cgoris
            i opions.inm_sims:
                go_sims  GO.GGOSims(
                    iooos.opn_i(opions.inm_sims, "r"))

                i opions.ogv > 1:
                    v  s()
                    or x in is(go_sims.vs()):
                        or xx in x:
                            v.(xx)
                    opions.sog.wri(
                        "# r go sims rom s: goi, simi\n" 
                        (opions.inm_sims,
                         n(go_sims),
                         n(v)))

                i opions.inm_mp_sims:
                    i opions.inm_mp_sims  "-":
                        oi  opions.so
                    s:
                        oi  iooos.opn_i(
                            opions.inm_mp_sims, "w")

                    oi.wri("GO\GOSim\n")
                    or go, go_sim in sor(is(go_sims.ims())):
                        oi.wri("s\s\n"  (go, go_sim))

                    i oi ! opions.so:
                        oi.cos()

                gn2go  GO.MpGO2Sims(gn2go, go_sims, onoogyonoogy)

                i opions.ogv > 1:
                    ngns, ncgoris, nmps, cons_pr_cgory  \
                        GO.ConGO(gn2go)
                    opions.sog.wri(
                        "# r go sim iring: i gns mpp o "
                        "i cgoris (i mps)\n"  (
                            ngns, ncgoris, nmps))

            #############################################################
            # Js mp o h gn is
            i opions.g_gns:
                g, bg, ng  [], [], []

                or gn, vv in is(gn2go.ims()):
                    or v in vv:
                        i v.mGOI  opions.g_gns:
                            i gn in gns:
                                g.ppn(gn)
                            i gn in bckgron:
                                bg.ppn(gn)
                            s:
                                ng.ppn(gn)

                # skip o nx GO css
                i no (bg or ng):
                    conin

                opions.so.wri(
                    "# gns in GO cgory s\n"  opions.g_gns)
                opions.so.wri("gn\s\n")
                or x in sor(g):
                    opions.so.wri("s\s\n"  ("g", x))
                or x in sor(bg):
                    opions.so.wri("s\s\n"  ("bg", x))
                or x in sor(ng):
                    opions.so.wri("s\s\n"  ("ng", x))

                E.ino("ngi, nbgi, nngi"  (n(g), n(bg), n(ng)))

                E.sop()
                sys.xi(0)

            #############################################################
            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'orgron',
                                     sgnis_nm)

            oi.wri("gn_i\ns\n"  ("\n".join(sor(orgron))))
            i opions.op_inm_prn:
                oi.cos()

            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'bckgron',
                                     sgnis_nm)

            # Jhro bg ix - s scion 'bi bckgron' or ssignmn
            oi.wri("gn_i\ns\n"  ("\n".join(sor(bckgron))))
            i opions.op_inm_prn:
                oi.cos()

            #############################################################
            # o h nysis
            go_rss  GO.AnysGO(gn2go, orgron, bckgron)

            i n(go_rss.mSmpGns)  0:
                E.wrn("s: no gns wih GO cgoris - nysis bor" 
                       gnis_nm)
                conin

            pirs  is(go_rss.mRss.ims())

            #############################################################
            # cc r or ch hypohsis
            i opions.r:
                rs, smps, mho  GO.compFDRs(go_rss,
                                                       orgron,
                                                       bckgron,
                                                       opions,
                                                       s_onoogy,
                                                       gn2go,
                                                       go2ino)
                or x, v in nmr(pirs):
                    v[1].mQV  rs[v[0]][0]
            s:
                rs, smps, mho  {}, {}, Non

            msgs.ppn("rs"  mho)

            i opions.sor_orr  "r":
                pirs.sor(kymb x: x[1].mQV)
            i opions.sor_orr  "rio":
                pirs.sor(kymb x: x[1].mRio)
            i opions.sor_orr  "pv":
                pirs.sor(kymb x: x[1].mPV)

            #############################################################
            #############################################################
            #############################################################
            # op h  rs
            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'ovr',
                                     sgnis_nm)

            GO.opRss(
                oi, pirs, go2ino, opions, rsrs, smpssmps)

            i opions.op_inm_prn:
                oi.cos()

            #############################################################
            #############################################################
            #############################################################
            # ir signiicn rss n op
            ir_pirs  GO.scSigniicnRss(pirs, rs, opions)

            nsc  n(ir_pirs)
            nsc_p  n([x or x in ir_pirs i x[1].mRio > 1])
            nsc_own  n(
                [x or x in ir_pirs i x[1].mRio < 1])

            ssr nsc_p + nsc_own  nsc

            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'rss',
                                     sgnis_nm)

            GO.opRss(oi,
                             ir_pirs,
                             go2ino,
                             opions,
                             rsrs,
                             smpssmps)

            i opions.op_inm_prn:
                oi.cos()

            #############################################################
            #############################################################
            #############################################################
            # sv rss or mi-gn-is nysis
            _rss.ppn(pirs)
            _signiicn_rss.ppn(ir_pirs)
            _gniss_wih_rss.ppn(gnis_nm)

            #############################################################
            #############################################################
            #############################################################
            # op prmrs
            ngns, ncgoris, nmps, cons_pr_cgory  \
                GO.ConGO(gn2go)

            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'prmrs',
                                     sgnis_nm)

            nbckgron  n(bckgron)
            i nbckgron  0:
                nbckgron  n(go_rss.mBckgronGns)

            oi.wri(
                "# inp go mppings or gn is 's' n cgory 's'\n" 
                (gnis_nm, s_onoogy))
            oi.wri("prmr\v\scripion\n")
            oi.wri("mpp_gns\i\mpp gns\n"  ngns)
            oi.wri(
                "mpp_cgoris\i\mpp cgoris\n"  ncgoris)
            oi.wri("mppings\i\mppings\n"  nmps)
            oi.wri("gns_in_g\i\gns in orgron\n" 
                          n(orgron))
            oi.wri(
                "gns_in_g_wih_ssignmn\i\gns in orgron wih GO ssignmns\n" 
                (n(go_rss.mSmpGns)))
            oi.wri(
                "gns_in_bg\i\inp bckgron\n"  nbckgron)
            oi.wri(
                "gns_in_bg_wih_ssignmn\i\gns in bckgron wih GO ssignmns\n"  (
                    n(go_rss.mBckgronGns)))
            oi.wri(
                "ssociions_in_g\i\ssociions in smp\n" 
                go_rss.mSmpConsTo)
            oi.wri(
                "ssociions_in_bg\i\ssociions in bckgron\n" 
                go_rss.mBckgronConsTo)
            oi.wri(
                "prcn_gns_in_g_wih_ssociion\s\prcn gns in smp wih GO ssignmns\n"  (
                    iooos.pry_prcn(n(go_rss.mSmpGns),
                                           n(orgron), "5.2")))
            oi.wri(
                "prcn_gns_in_bg_wih_ssociions\s\prcn gns bckgron wih GO ssignmns\n"  (
                    iooos.pry_prcn(n(go_rss.mBckgronGns),
                                           nbckgron, "5.2")))
            oi.wri(
                "signiicn\i\signiicn rss rpor\n"  nsc)
            oi.wri(
                "signiicn_p\i\signiicn p-rg rss rpor\n"  nsc_p)
            oi.wri(
                "signiicn_own\i\signiicn p-rg rss rpor\n"  nsc_own)
            oi.wri(
                "hrsho\6.4\signiicnc hrsho\n"  opions.hrsho)

            i opions.op_inm_prn:
                oi.cos()

            smmry.ppn("\".join(mp(sr, (
                gnis_nm,
                s_onoogy,
                nsc,
                opions.hrsho,
                ngns,
                ncgoris,
                nmps,
                n(orgron),
                n(go_rss.mSmpGns),
                nbckgron,
                n(go_rss.mBckgronGns),
                go_rss.mSmpConsTo,
                go_rss.mBckgronConsTo,
                iooos.pry_prcn(
                    n(go_rss.mSmpGns), n(orgron), "5.2"),
                iooos.pry_prcn(
                    n(go_rss.mBckgronGns), nbckgron, "5.2"),
                ",".join(msgs)))) + "\n")

            #############################################################
            #############################################################
            #############################################################
            # op h g prns
            oi  GO.gFiNm(opions,
                                     gos_onoogy,
                                     scion'wihgns',
                                     sgnis_nm)

            GO.opRss(oi, pirs, go2ino, opions,
                             rsrs,
                             smpssmps,
                             gn2gogn2go,
                             orgronorgron,
                             gn2nmgn2nm)

            i opions.op_inm_prn:
                oi.cos()

        i n(gniss) > 1:

            ###################################################################
            # op vrios smmry is
            # signiicn rss
            GO.opMipGnLisRss(_signiicn_rss,
                                             _gniss_wih_rss,
                                             s_onoogy,
                                             go2ino,
                                             opions,
                                             scion'signiicn')

            #  rss
            GO.opMipGnLisRss(_rss,
                                             _gniss_wih_rss,
                                             s_onoogy,
                                             go2ino,
                                             opions,
                                             scion'')

            i opions.comp_pirwis:
                GO.pirwisGOEnrichmn(_rss,
                                        _gniss_wih_rss,
                                        s_onoogy,
                                        go2ino,
                                        opions)

    oi_smmry  opions.so
    oi_smmry.wri("".join(smmry))

    E.sop()


i __nm__  "__min__":
    sys.xi(min())
