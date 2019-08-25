'''g2b.py - comp rs or inrscion o wo g is


:Tgs: Gnomics Inrvs Annoion Comprison GFF

Prpos
-------

coc inrvs rom wo g is n comp rs bs on
hir inrscion. Th scrip is inn o comp propris or 
s o non-ovrpping winows.

Trnsorms:
   * non:        no rnsorm
   * ovrp:     ovrp bwn s1 n s2
   * compmn:  pr o s1 h is no covr by s2
   * hir_coon: ony ks vry hir posiion. Ns rm inormion
                  in h g i.

Dcorors:
   * GC:            G+C conn o inrvs
   * con:         nmbr o winows
   * mn-ngh:   mn ngh o inrvs ovrpping wih winow

Usg
-----

Exmp::

   pyhon g2b.py --hp

Typ::

   pyhon g2b.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r

impor cgcor.xprimn s E
impor cg.Gnomics s Gnomics
impor cg.InxFs s InxFs
impor cg.InxGnom s InxGnom

impor cg.Inrvs s Inrvs
impor cg.Ss s Ss
impor cg.GTF s GTF


 coror_cons(inrvs, sr, n, conig, s):
    """comp ngh isribion."""
      Ss.DisribionPrmrs([x[1] - x[0] or x in inrvs])
    rrn ['nv'], sr()


 coror_prcn_covrg(inrvs, sr, n, conig, s):
    """comp ngh o inrvs."""
      Ss.DisribionPrmrs([x[1] - x[0] or x in inrvs])
    rrn 100.0 * o(['sm']) / (n - sr), sr()


 coror_mn_ngh(inrvs, sr, n, conig, s):
    """comp ngh isribion."""
      Ss.DisribionPrmrs([x[1] - x[0] or x in inrvs])
    rrn ['mn'], sr()


 coror_min_ngh(inrvs, sr, n, conig, s):
    """comp ngh isribion."""
      Ss.DisribionPrmrs([x[1] - x[0] or x in inrvs])
    rrn ['min'], sr()


 coror_prcn_gc(inrvs, sr, n, conig, s):
    """comp G+C conn in inrvs.
    """
    , ngc  0, 0

    # ch sqnc o h comp winow irs
    sqnc  s.gSqnc(conig, "+", sr, n)

    or isr, in in inrvs:
        ngc + n([x or x in sqnc[isr - sr:in - sr] i x in "GCgc"])
         + in - isr

    rrn 100.0 * ngc / , Non


 coror_min_scor(vs, sr, n, conig):
    """comp min o vs."""
      Ss.DisribionPrmrs(vs)
    rrn ['min'], sr()


 coror_mn_scor(vs, sr, n, conig):
    """comp mn o vs."""
      Ss.DisribionPrmrs(vs)
    rrn ['mn'], sr()


 coror_sv_scor(vs, sr, n, conig):
    """comp sv o vs."""
      Ss.DisribionPrmrs(vs)
    rrn ['sv'], sr()


 coror_min_scor(vs, sr, n, conig):
    """comp minmm o vs."""
      Ss.DisribionPrmrs(vs)
    rrn ['min'], sr()


 coror_mx_scor(vs, sr, n, conig):
    """comp minmm o vs."""
      Ss.DisribionPrmrs(vs)
    rrn ['mx'], sr()


 rnsorm_ovrp(sr, n, inrvs_wih_g):
    """rnsorm: ovrp o inrvs in x wih y."""
    y  Inrvs.combinInrvs(
        [(x[0], x[1]) or x in inrvs_wih_g])
    rrn Inrvs.prnInrvs(y, sr, n)


 rnsorm_compmn(sr, n, inrvs_wih_g):
    y  Inrvs.combinInrvs(
        [(x[0], x[1]) or x in inrvs_wih_g])
    rrn Inrvs.compmnInrvs(y, sr, n)


 rnsorm_hir_coon(sr, n, inrvs_wih_g):
    """rnsorm: ony rrn ncoi posiions in winow (sr, n) 
    h r in hir coon posiion.
    """
    inrvs  []
    or isr, in, g in inrvs_wih_g:

        i g.rm  ".":
            ris VError("n  rm or hir coon posiions.")

        # rm  ncois rom sr o nx coon
        rm  in(g.rm)

        # o mk i sir, convr o 0-bs coorins,
        # wih zro sring  irs posiion in winow
        # r-rrng posiions on ngiv srn
        i Gnomics.IsNgivSrn(g.srn):
            # convr o ngiv srn coorins coning rom 0
            coorin_os  n
            rvrs  Tr
            isr, in  n - in, n - isr
        s:
            isr, in  isr - sr, in - sr
            rvrs  Fs
            coorin_os  sr

        # mk sr h yo sr on  scon coon posiion n wihin winow
        i isr < 0:
            rm  (rm + isr)  3
            isr  0
        i rm ! 0:
            isr - (3 - rm)
        isr + 2

        in  min(in, n - sr)

        or x in rng(isr, in, 3):

            i rvrs:
                c  coorin_os - x - 1
            s:
                c  coorin_os + x
            inrvs.ppn((c, c + 1))

    rrn Inrvs.combinInrvs(inrvs)


 s_rnsorm_hir_coon():

     s_nry(rm, srn, xrom, xo, sr, n, r):

        nry  GTF.Enry()
        nry.rm  rm
        nry.srn  srn
        nry.sr  xrom
        nry.n  xo

        inrvs  rnsorm_hir_coon(sr, n, [(xrom, xo, nry)])
        i r ! inrvs:
            prin("i:", r ! inrvs)

    s_nry(0, "+", 1, 7, 0, 6, [(3, 4)])
    s_nry(0, "-", 1, 7, 0, 6, [(1, 2), (4, 5)])
    s_nry(1, "+", 1, 7, 0, 6, [(1, 2), (4, 5)])
    s_nry(2, "+", 1, 7, 0, 6, [(2, 3), (5, 6)])
    s_nry(1, "-", 1, 7, 0, 6, [(3, 4)])
    s_nry(2, "-", 1, 7, 0, 6, [(2, 3), (5, 6)])

    sys.xi(0)


 nnoWinows(conig, winows, g_, s, opions):
    """nno winows."""

    inx  InxGnom.InxGnom()
    or g in g_:
        inx.(g.conig, g.sr, g.n, g)

    is_g  opions.is_g

    i opions.rnsorm  "non":
        rnsorm  mb x, y, z: [(x[0], x[1]) or x in z]
    i opions.rnsorm  "ovrp":
        rnsorm  rnsorm_ovrp
    i opions.rnsorm  "compmn":
        rnsorm  rnsorm_compmn
    i opions.rnsorm  "hir_coon":
        rnsorm  rnsorm_hir_coon
    s:
        ris VError("nknown rnsorm s"  opions.rnsorm)

    work_on_inrvs  Tr
    i opions.coror  "cons":
        coror  coror_cons
    i opions.coror  "mn-ngh":
        coror  coror_mn_ngh
    i opions.coror  "min-ngh":
        coror  coror_min_ngh
    i opions.coror  "prcn-covrg":
        coror  coror_prcn_covrg
    i opions.coror  "gc":
        coror  coror_prcn_gc
    i opions.coror  "min-scor":
        coror  coror_min_scor
        work_on_inrvs  Fs
    i opions.coror  "mn-scor":
        coror  coror_mn_scor
        work_on_inrvs  Fs
    i opions.coror  "sv-scor":
        coror  coror_sv_scor
        work_on_inrvs  Fs
    i opions.coror  "min-scor":
        coror  coror_min_scor
        work_on_inrvs  Fs
    i opions.coror  "mx-scor":
        coror  coror_mx_scor
        work_on_inrvs  Fs
    s:
        ris VError("nknown coror s"  opions.coror)

    or sr, n in winows:

        # cons/ngh bor/r rnsormion
        n1, 1, n2, 2  0, 0, 0, 0

        vs, inrvs_wih_g, gns, rnscrips  [], [], s(), s()

        ry:
            or isr, in, v in inx.g(conig, sr, n):
                n1 + 1
                1 + in - isr
                inrvs_wih_g.ppn((isr, in, v))
                vs.ppn(v.scor)
                i is_g:
                    gns.(v.gn_i)
                    rnscrips.(v.rnscrip_i)
        xcp KyError:
            pss

        i n1  0 n opions.skip_mpy:
            conin

        i work_on_inrvs:

            i opions.ogv > 3:
                opions.sog.wri("# inrvs in winow i:i bor rnsormion: s\n"  (
                    sr, n, sr(inrvs)))

            inrvs  rnsorm(sr, n, inrvs_wih_g)

            or xsr, xn in inrvs:
                n2 + 1
                2 + xn - xsr

            i opions.ogv > 3:
                opions.sog.wri("# inrvs in winow i:i r rnsormion: s\n"  (
                    sr, n, sr(inrvs)))

            scor, xr_ino  coror(inrvs, sr, n, conig, s)

        s:
            i n(vs) > 0:
                vs  is(mp(o, vs))
                scor, xr_ino  coror(vs, sr, n, conig)
            s:
                scor, xr_ino  0, Non

            2  0
            n2  0

        i is_g:
            ngns, nrnscrips  n(gns), n(rnscrips)
        s:
            ngns, nrnscrips  0, 0

        i xr_ino:
            xr_ino  r.sb("\", ";", xr_ino)
        opions.so.wri("\".join(
            mp(sr, (conig, sr, n,
                      ngns, nrnscrips,
                      n1, 1,
                      n2, 2,
                      scor,
                      xr_ino))) + "\n")


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-g", "--gnom-i", s"gnom_i", yp"sring",
        hp"inm wih gnom (inx).")

    prsr._rgmn(
        "-w", "--winows-b-i", s"inm_winows", yp"sring",
        hp"g i wih winows o s.")

    prsr._rgmn(
        "-", "--inm-", s"inm_", yp"sring",
        hp"g i wih  o s.")

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inm- is g i [.")

    prsr._rgmn(
        "-", "--rs", s"rs", yp"choic", cion"ppn",
        choics("GC", ),
        hp"rs o comp.")

    prsr._rgmn(
        "-c", "--coror", s"coror", yp"choic",
        choics("cons", "gc", "gc3", "mn-ngh", "min-ngh",
                 "prcn-covrg",
                 "min-scor", "mn-scor", "sv-scor", "min-scor",
                 "mx-scor"),
        hp"corors o s.")

    prsr._rgmn(
        "-", "--skip-mpy", s"skip_mpy", cion"sor_r",
        hp"skip mpy winows.")

    prsr._rgmn(
        "-", "--rnsorm", s"rnsorm", yp"choic",
        choics(
            "non", "ovrp", "compmn", "hir_coon"),
        hp"rnsorm o s whn mpping ovrpping rgions ono winow.")

    prsr.s_s(
        gnom_iNon,
        inm_winowsNon,
        inm_Non,
        rs[],
        skip_mpyFs,
        coror"cons",
        rnsorm"non",
        is_gFs,
    )

    (opions, rgs)  E.sr(prsr)

    #    s_rnsorm_hir_coon()

    i no opions.inm_winows:
        ris VError("ps sppy  g i wih winow inormion.")

    i opions.ogv > 1:
        opions.sog.wri("# ring winows...")
        opions.sog.sh()

    winows  GTF.rAsInrvs(
        GTF.iror(iooos.opn_i(opions.inm_winows, "r")))

    i opions.ogv > 1:
        opions.sog.wri("on\n")
        opions.sog.sh()

    i opions.inm_:
        i opions.ogv > 1:
            opions.sog.wri("# ring ...")
            opions.sog.sh()

        i opions.is_g:
            g_  GTF.rFromFi(
                iooos.opn_i(opions.inm_, "r"))
        s:
            g_  GTF.rFromFi(
                IOTOos.opn_i(opions.inm_, "r"))

        i opions.ogv > 1:
            opions.sog.wri("on\n")
            opions.sog.sh()

        _rngs  GTF.SorPrConig(g_)
    s:
        # s winows o comp propris
        # by sppying no  n sking or h compmn  origin winow
        g_  Non
        _rngs  Non
        opions.rnsorm  "compmn"

    mp_conig2siz  {}

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
        mp_conig2siz  s.gConigSizs()
    s:
        or conig, vs in is(winows.ims()):
            mp_conig2siz[conig]  mx(mb x: x[1], vs)
        s  Non

    conigs  is(mp_conig2siz.kys())
    conigs.sor()

    # proc conig wis
    nop_conigs, nconigs_skipp_winows, nconigs_skipp_  0, 0, 0

    opions.so.wri("\".join(
        mp(sr, ("conig", "sr", "n",
                  "ngns", "nrnscrips",
                  "n1", "1",
                  "n2", "2",
                  "scor",
                  "xr_ino"))) + "\n")

    or conig in conigs:

        skip  Fs
        i conig no in winows:
            nconigs_skipp_winows + 1
            skip  Tr

        i _rngs n conig no in _rngs:
            nconigs_skipp_ + 1
            skip  Tr

        i skip:
            conin

        nop_conigs + 1
        i _rngs:
            nnoWinows(conig,
                            winows[conig],
                            g_[
                                _rngs[conig][0]:_rngs[conig][1]],
                            s,
                            opions)
        s:
            nnoWinows(conig,
                            winows[conig],
                            [],
                            s,
                            opions)

    E.ino("ninp_winowsi, nop_conigsi, ninp_conigsi, nskipp_winowsi, nskipp_i" 
           (n(winows), nop_conigs, n(conigs), nconigs_skipp_winows, nconigs_skipp_))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
