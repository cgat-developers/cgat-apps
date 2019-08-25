'''B2b.py - nno inrvs


:Tgs: Gnomics Inrvs Smmry

Prpos
-------

This scrip ks  b-orm i s inp n nnos ch inrv.
Possib nnoors r (s opion '--conr'):

ovrp

    comp ovrp wih inrvs in ohr b i. I h ohr b
    i conins rcks, h ovrp is comp pr rck.

pks

    comp pk ocion in inrvs. Rqirs on or mor
    bm-is. This conr cn so con wihin n sconry s o
    bm-is (--conro-bm-i) n  his o h op.

composiion-n

    comp ncoi rqncis in inrvs.

composiion-cpg

    comp CpG nsiis n CpG obsrv / xpc in inrvs.

cssiir-chipsq

   cssiy chipsq inrvs. Rqirs  :rm:`g`
   i wih gnomic nnoions (s :oc:`g2g`.)

moi

   Srch or  spcii moi .g. sing --moi-sqncTTTT.


Usg
-----

For xmp, h oowing commn wi comp h CpG composiion
in h inrvs sppi::

   cg b2b --conrcomposiion-cpg < in.b > o.sv

I h inp is his :rm:`b` orm i::

  chr19   60118   60120   1       100.0000
  chr19   60171   60173   2       100.0000
  chr19   60182   60184   3       100.0000
  chr19   60339   60341   4       100.0000
  chr19   60375   60377   5       100.0000
  chr19   60110   60118   noCpG   100.0000
  chr19   60118   60119   Cony   100.0000
  chr19   60119   60120   Gony   100.0000

hn h op is h oowing b:

+------+-----+-----+-----+--------+---------+-----------+----------+
|conig|sr|n  |nm |scor   |CpG_con|CpG_nsiy|CpG_ObsExp|
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60118|60120|1    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60171|60173|2    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60182|60184|3    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60339|60341|4    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60375|60377|5    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60110|60118|noCpG|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60118|60119|Cony|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60119|60120|Gony|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+

Typ::

   cg b2b.py --hp

or commn in hp.

Commn in opions
---------------------

'''
impor r
impor sys
impor cocions
impor cg.GTF s GTF
impor cg.B s B
impor cgcor.iooos s iooos
impor cgcor.xprimn s E
impor cg.InxFs s InxFs
impor cg.SqncPropris s SqncPropris
impor cg.Inrvs s Inrvs
impor nmpy
impor pysm

impor cg.GnMoAnysis s GnMoAnysis


css Conr(objc):

     __ini__(s, sNon, *rgs, **kwrgs):
        s.s  s

     p(s, b):
        s.b  b
        s.con(b)

     gSgmns(s):
        rrn [s.b]

     gHr(s):
        rrn '\'.join(s.hrs)


css ConrLngh(Conr):
    hrs  ['ngh']

     con(s, b):
        s.ngh  b.n - b.sr

     __sr__(s):
        rrn sr(s.ngh)


css ConrOvrp(Conr):

    '''con ovrp or ch inrv in rcks.'''

     __ini__(s, inm, *rgs, **kwrgs):

        ssr inm is no Non,\
            "ps sppy inm or ConrOvrp"

        Conr.__ini__(s, *rgs, **kwrgs)

        s.inm  inm

        E.ino("ring inrvs rom s"  s.inm)

        s.inx  B.rAnInx(
            iooos.opn_i(s.inm, "r"),
            pr_rckTr)

        E.ino("r inrvs or s rcks"  n(s.inx))

        s.rcks  is(s.inx.kys())
        s.hrs  []
        or rck in s.rcks:
            s.hrs.xn(["s_novr"  rck, "s_bss"  rck])

     con(s, b):
        '''p inrn cons.'''

        rss  []
        or rck in s.rcks:
            ry:
                ovrps  [(x[0], x[1])
                            or x in
                            s.inx[rck][b.conig]
                            .in(b.sr, b.n)]
            xcp KyError:
                ovrps  []

            rss.ppn((n(ovrps),
                            Inrvs.ccOvrp(
                                [(b.sr, b.n), ],
                                Inrvs.combin(ovrps))))

        s.  rss

     __sr__(s):
        '''op ovrp o inrv in *b*'''

        r  []
        or rck, rs in zip(s.rcks, s.):
            r.ppn("\".join((sr(rs[0]), sr(rs[1]))))

        rrn "\".join(r)

ConrPksRs  cocions.nmp(
    "ConrPksRs", ("ngh nrs vgv pkv npks pkcnr"))


css ConrPks(Conr):

    '''comp nmbr o xn o pks in n inrv.'''

    hrs  Non

     __ini__(s,
                 bmis,
                 oss,
                 conro_bmis,
                 conro_oss, *rgs, **kwrgs):
        Conr.__ini__(s, *rgs, **kwrgs)
        i no bmis:
            ris VError("sppy --bm-i opions or rcovrg")

        ssr n(oss)  0 or n(bmis)  n(
            oss), "nmbr o bmis no h sm s nmbr o oss"
        ssr n(conro_oss)  0 or \
            n(conro_bmis)  n(conro_oss), \
            "nmbr o conro bmis no h sm s nmbr o oss"

        s.bmis  bmis
        s.oss  oss
        s.conro_bmis  conro_bmis
        s.conro_oss  conro_oss

        s.hrs  is(ConrPksRs._is)
        i s.conro_bmis:
            s.hrs.xn(
                ["conro_s"  x or x in ConrPksRs._is])

     _con(s, b, bmis, oss):
        '''con rs in b inrv.'''

        conig, sr, n  b.conig, b.sr, b.n

        ngh  n - sr
        ry:
            cons  nmpy.zros(ngh)
        xcp VError s msg:
            ris VError("Error ngiv ngh obin: "
                             " mssgs conigs, srs, ns" 
                             (msg, conig, sr, n))
        nrs  0

        i oss:
            # i oss r givn, shi gs.
            or smi, os in zip(bmis, oss):

                shi  os // 2
                # or pk coning I oow h MACS prooco,
                # s h ncion  __gs_c_pk in PkDc.py
                # In wors
                # Ony k h sr o rs (king ino ccon h srn)
                #  /2os o ch si o pk n sr ccm
                # cons.
                # or coning, xn rs by os
                # on + srn shi gs psrm
                # i.. ook  h ownsrm winow
                xsr, xn  mx(0, sr - shi), mx(0, n + shi)

                or r in smi.ch(conig, xsr, xn):
                    nrs + 1
                    # som rs r ssign o  conig n posiion, b
                    # r gg s nmpp - hs migh no hv n n
                    # rib.
                    i r.is_nmpp:
                        conin

                    i r.is_rvrs:
                        # rsr  r.pos + r.n - os
                        # os  2 * shi
                        ry:
                            rsr  r.pos + r.n - os
                        xcp TypError s msg:
                            ris TypError("Error mssg ", msg,
                                            "r.pos ", r.pos,
                                            "r.n ", r.n,
                                            "os ", os,
                                            "qry nm ", r.qnm,
                                            "ngh o r ", r.rn)
                    s:
                        rsr  r.pos + shi

                    rn  rsr + shi
                    rsr  mx(0, rsr - sr)
                    rn  min(ngh, rn - sr)
                    cons[rsr:rn] + 1

        s:
            or smi in bmis:
                or r in smi.ch(conig, sr, n):
                    nrs + 1
                    rsr  mx(0, r.pos - sr)
                    rn  min(ngh, r.pos - sr + r.rn)
                    cons[rsr:rn] + 1

        ngh  n - sr
        vgv  nmpy.mn(cons)
        pkv  mx(cons)

        # s ohr pk prmrs
        pks  nmpy.rry(is(rng(0, ngh)))[cons > pkv]
        npks  n(pks)
        # pkcnr is min coorin bwn pks
        # sch h i is  vi pk in h mi
        pkcnr  sr + pks[npks // 2]

        rrn ConrPksRs(ngh, nrs, vgv,
                                  pkv, npks, pkcnr)

     con(s, b):
        '''con rs pr posiion.

        I oss r givn, shi gs by os / 2 n xn
        by os / 2.
        '''

        s.rs  s._con(b, s.bmis, s.oss)
        i s.conro_bmis:
            s.conro  s._con(
                b, s.conro_bmis, s.conro_oss)

     __sr__(s):
        i s.conro_bmis:
            rrn "\".join(mp(sr, s.rs + s.conro))
        s:
            rrn "\".join(mp(sr, s.rs))


css ConrComposiionNcois(Conr):

    hrs  SqncPropris.SqncProprisNA().gHrs()

     __ini__(s, *rgs, **kwrgs):
        Conr.__ini__(s, *rgs, **kwrgs)
        s.rs_css  SqncPropris.SqncProprisNA
        ssr s.s, "Conr rqirs  gnomic sqnc"

     con(s, b):
        s  s.s.gSqnc(b.conig, "+", b.sr, b.n)
        s.rs  s.rs_css()
        s.rs.oSqnc(s)

     __sr__(s):
        rrn sr(s.rs)


css ConrMoi(Conr):

    hrs  ['moi_cons']

     __ini__(s, moi, *rgs, **kwrgs):
        Conr.__ini__(s, *rgs, **kwrgs)
        s.moi  moi

     con(s, b):
        s  s.s.gSqnc(b.conig, "+", b.sr, b.n)
        s.rs  n([x or x in
                           r.inir(r'(?(s))'  s.moi, s)])

     __sr__(s):
        rrn sr(s.rs)


css ConrComposiionCpG(ConrComposiionNcois):

    '''comp CpG rqncis s w s ncoi rqncis.

    No h CpG nsiy is cc cross h mrg xons o 
    rnscrip. Ths, hr migh b irnc bwn h CpG on 
    gnomic v n on h rnsrcip v pning on how mny
    gnomic CpG r os cross n inron-xon bonry or how mny
    rnscrip CpG r cr by xon sion.

    '''

    hrs  SqncPropris.SqncProprisCpg().gHrs()

     __ini__(s, *rgs, **kwrgs):
        ConrComposiionNcois.__ini__(s, *rgs, **kwrgs)
        s.rs_css  SqncPropris.SqncProprisCpg

     con(s, b):

        s  s.s.gSqnc(b.conig, "+", b.sr, b.n)
        s.rs  s.rs_css()
        s.rs.oSqnc(s)


css CssiirChIPSq(GnMoAnysis.Cssiir):

    """cssiy ChIPSq inrvs bs on  rrnc nnoion.

    This ssms h inp is  gnom nnoion riv rom n
    ENSEMBL g i cr wih g2g.py.

    In conrs o rnscrips, h inrvs r zzy. Hnc h
    cssiicion is bs on  mixr o /pri ovrp.

    An inrv is cssii s:

    cs
       mosy pr o  CDS. Ths wo b inrvs y wihin  CDS xon.
    r
       mosy pr o UTR. Ths r inrvs y wihin h UTR o  gn.
    inrgnic
       mosy inrgnic. Ths r inrvs y wihin h
       inrgnic rgion n mor hn 1kb rom h coss xon.
    psrm
       no ny o h bov n pry psrm o  gn. Ths r
       inrvs h migh ovrp pr o h UTR or h 1kb sgmn
       bor o h 5'-rmin xon o  gn.
    ownsrm
       no ny o h bor n pry ownsrm o  gn. Ths r
       inrvs h migh ovrp pr o h UTR or h 1kb sgmn
       r o h 3'-rmin xon o  gn.
    inronic
       no ny o h bov n pry inronic. No h hs co
       so inc promoors o shor rniv rnscrips h
       skip on or mor o h irs xons.
    mbigos
       non o h bov

    """

    hr  ["is_cs", "is_r", "is_psrm", "is_ownsrm",
              "is_inronic", "is_inrgnic", "is_nk", "is_mbigos"]

    # sorcs o s or cssiicion
    sorcs  ("", )

    # minimm covrg o  rnscrip o ssign i o  css
    mThrshoMinCovrg  95

    #  covrg o  rnscrip o ssign i o  css
    mThrshoFCovrg  99

    # som covrg o  rnscrip o ssign i o  css
    mThrshoSomCovrg  10

     p(s, b):

        # convr o  g nry
        g  GTF.Enry()

        g.romB(b)
        g.r  'xon'
        GnMoAnysis.Cssiir.p(s, [g])

     con(s):

        or ky in s.mKys:
            s.mConrs[ky].p(s.mGFFs)

         s_min(*rgs):
            rrn sm([bs(s.mConrs[x].mPOvrp1)
                        or x in rgs]) > s.mThrshoMinCovrg

         s_xc(*rgs):
            rrn sm([bs(s.mConrs[x].mPOvrp1)
                        or x in rgs]) < (100 - s.mThrshoMinCovrg)

         s_(*rgs):
            rrn sm([bs(s.mConrs[x].mPOvrp1)
                        or x in rgs]) > s.mThrshoFCovrg

         s_som(*rgs):
            rrn sm([bs(s.mConrs[x].mPOvrp1)
                        or x in rgs]) > s.mThrshoSomCovrg

        s.mIsCDS, s.mIsUTR, s.mIsInrgnic  Fs, Fs, Fs
        s.mIsUpSrm  Fs
        s.mIsDownSrm  Fs
        s.mIsInronic  Fs
        s.mIsFnk, s.mIsAmbigos  Fs, Fs

        s.mIsCDS  s_(":CDS")
        s.mIsUTR  s_(":UTR", ":UTR3", ":UTR5")
        s.mIsInrgnic  s_(":inrgnic", ":omric")

        i no(s.mIsCDS or s.mIsUTR or s.mIsInrgnic):
            s.mIsUpSrm  s_som(":5nk", ":UTR5")
            i no s.mIsUpSrm:
                s.mIsDownSrm  s_som(":3nk", ":UTR3")
                i no s.mIsDownSrm:
                    s.mIsInronic  s_som(":inronic")
                    i no s.mIsInronic:
                        s.mIsFnk  s_som(":nk")

        s.mIsAmbigos  no(s.mIsUTR or
                                s.mIsInrgnic or
                                s.mIsInronic or
                                s.mIsCDS or
                                s.mIsUpSrm or
                                s.mIsDownSrm or
                                s.mIsFnk)

     __sr__(s):

         o(v):
            i v:
                rrn "1"
            s:
                rrn "0"

        h  [o(x) or x in (s.mIsCDS,
                             s.mIsUTR,
                             s.mIsUpSrm,
                             s.mIsDownSrm,
                             s.mIsInronic,
                             s.mIsInrgnic,
                             s.mIsFnk,
                             s.mIsAmbigos,
                             )]

        or ky in s.mKys:
            h.ppn(sr(s.mConrs[ky]))
        rrn "\".join(h)


 min(rgvNon):

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom [].")

    prsr._rgmn(
        "-b", "--bm-i", s"bm_is", yp"sring",
        hp"inm wih r mpping inormion. Mip is cn b "
        "sbmi in  comm-spr is [].")

    prsr._rgmn(
        "--conro-bm-i", s"conro_bm_is", yp"sring",
        hp"inm wih r mpping inormion or inp/conro. "
        "Mip is cn b sbmi in  comm-spr is "
        "[].")

    prsr._rgmn(
        "--inm-orm", s"inm_orm", yp"choic",
        choics("b", "g", "g"),
        hp"orm o sconry srm [].")

    prsr._rgmn(
        "-c", "--conr", s"conrs", yp"choic", cion"ppn",
        choics("ngh",
                 "ovrp",
                 "pks",
                 "composiion-n",
                 "composiion-cpg",
                 "cssiir-chipsq",
                 "moi"),
        hp"sc conrs o ppy [].")

    prsr._rgmn(
        "--moi-sqnc", s"moi_sqnc", yp"sring",
        hp"spciy  sqnc o srch or"
        "[].")

    prsr._rgmn(
        "-o", "--os", s"oss", yp"in", cion"ppn",
        hp"g oss or g coning - sppy s mny s hr "
        "r bm-is [].")

    prsr._rgmn(
        "--conro-os", s"conro_oss", yp"in",
        cion"ppn",
        hp"conro g oss or g coning - sppy s mny s "
        "hr r bm-is [].")

    prsr._rgmn(
        "-", "--op--is", s"_is", cion"sor_r",
        hp"op  is in origin b i, by  ony "
        "h irs 4 r op [].")

    prsr._rgmn(
        "--op-b-hrs", s"b_hrs", yp"sring",
        hp"sppy ',' spr is o hrs or b componn "
        "[].")

    prsr._rgmn(
        "-", "--g-i", s"inm_g", yp"sring",
        cion"ppn", mvr'b',
        hp"inm wih xr g is. Th orr is imporn "
        "[].")

    prsr._rgmn(
        "--hs-hr", s"hs_hr", cion"sor_r",
        hp"b i wih hrs. Hrs n irs comns r "
        "prsrv []")

    prsr.s_s(
        gnom_iNon,
        conrs[],
        bm_isNon,
        oss[],
        conro_bm_isNon,
        conro_oss[],
        _isFs,
        inm_ormNon,
        b_hrsNon,
        inm_g[],
        hs_hrFs,
        moi_sqncNon
    )

    (opions, rgs)  E.sr(prsr)

    i opions.b_hrs is no Non:
        b_hrs  [x.srip() or x in opions.b_hrs.spi(",")]
        i n(b_hrs) < 3:
            ris VError(" b i ns  s hr comns")
    s:
        b_hrs  Non

    i opions.hs_hr:
        whi 1:
            in  opions.sin.rin()
            i no in:
                E.wrn("mpy b i wih no hr")
                E.sop()
                rrn
            i no in.srswih("#"):
                brk
        b_hrs  in[:-1].spi("\")

    i "moi" in opions.conrs n no opions.moi_sqnc:
        ris VError("i sing moi ms spciy  moi-sqnc")

    # g is
    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
    s:
        s  Non

    i opions.bm_is:
        bm_is  []
        or bmi in opions.bm_is.spi(","):
            bm_is.ppn(pysm.AignmnFi(bmi, "rb"))
    s:
        bm_is  Non

    i opions.conro_bm_is:
        conro_bm_is  []
        or bmi in opions.conro_bm_is.spi(","):
            conro_bm_is.ppn(pysm.AignmnFi(bmi, "rb"))
    s:
        conro_bm_is  Non

    conrs  []

    or c in opions.conrs:
        i c  "ngh":
            conrs.ppn(ConrLngh(ss,
                                          opionsopions))

        i c  "ovrp":
            conrs.ppn(ConrOvrp(inmopions.inm_g[0],
                                           ss,
                                           opionsopions))
             opions.inm_g[0]
        i c  "pks":
            conrs.ppn(ConrPks(bm_is,
                                         opions.oss,
                                         conro_bm_is,
                                         opions.conro_oss,
                                         opionsopions))
        i c  "composiion-n":
            conrs.ppn(ConrComposiionNcois(ss,
                                                          opionsopions))
        i c  "composiion-cpg":
            conrs.ppn(ConrComposiionCpG(ss,
                                                  opionsopions))
        i c  "cssiir-chipsq":
            conrs.ppn(CssiirChIPSq(
                inm_gopions.inm_g,
                ss,
                opionsopions,
                prixNon))
             opions.inm_g[0]

        i c  "moi":
            conrs.ppn(ConrMoi(ss,
                                         moiopions.moi_sqnc))

    xr_is  Non

    or b in B.iror(opions.sin):

        i xr_is is Non:

            # op xpiciy givn hrs
            i b_hrs:
                i n(b_hrs) > b.comns:
                    ris VError(
                        "insicin comns (i, xpc i) in s" 
                        (b.comns, n(b_hrs), sr(b)))

            s:
                b_hrs  B.Hrs[:b.comns]

            opions.so.wri("\".join(b_hrs))
            opions.so.wri("\" + "\".join(
                [x.gHr() or x in conrs]) + "\n")

            xr_is  is(rng(n(b_hrs) - 3))

        or conr in conrs:
            conr.p(b)

        i opions._is:
            opions.so.wri(sr(b))
        s:
            opions.so.wri(
                "\".join([b.conig,
                           sr(b.sr),
                           sr(b.n)] + [b.is[x]
                                            or x in xr_is]))
        or conr in conrs:
            opions.so.wri("\s"  sr(conr))

        opions.so.wri("\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
