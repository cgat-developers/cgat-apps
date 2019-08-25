'''bm2bm.py - moiy bm is


Prpos
-------

This scrip rs  :rm:`bm` orm i rom sin, prorms n
cion (s mhos bow) hn ops  moii :rm:`bm`
orm i on so.

.. no::
   Yo n o rirc ogging inormion o  i (vi -L) or rn i o
   vi -v 0 in orr o g  vi sm/bm i.

Docmnion
-------------

Th scrip impmns h oowing mhos:

``s-nh``

   s h NH g. Som oos (bowi, bw) o no s h NH g.
   I s, his opion wi s h NH g (or mpp rs).
   This opion rqirs h bm/sm i o b sor by r nm.

``ns-nmpp_mpq``

   som oos s h mpping qiy o nmpp rs. This
   css  vioion in h Picr oos.

``ir``

   rmov ignmns bs on  vriy o gs. Th iring mho
   is rmin by h ``--ir-mho`` opion. Ths my b
   ``niq``, ``non-niq``, ``mpp``, ``NM`` or ``CM``.  I
   ``niq`` is s, ony niqy mpping rs wi b op. I
   ``non-niq`` is s hn ony mi-mpping rs wi b
   op. This mho irs chcks or h NH g - i s,  niq
   mch sho hv  mos NH1 his.  I no s, h mho chcks
   or BWA gs. Crrny i chcks i X0 is s (X0Nmbr o bs
   his on by BWA).  I ``mpp`` is givn, nmpp rs wi b
   rmov. I ``NM`` or ``CM`` is s, h ignmn o rs in wo
   sm is (inp n rrnc) is compr n ony rs wih 
   owr nmbr o mismchs in h inp compr o h rrnc
   sm i wi b kp. I ``CM`` is s, h coorspc mismch
   g (or ABI Soi rs) wi b s o con irncs o h
   rrnc sm i. By , h ``NM`` (nmbr o mismchs)
   g is s. Th g h is s ns o prsn in boh inp
   sm i n h rrnc sm i. I ``niq`` is givn his
   wi NOT rmov ny nmpp rs.  This cn b chiv by
   proviing h ``ir`` opion wic, onc ch wih ``mpp``
   n ``niq``.

   .. no::

      Th ir mhos cn' crrny combin wih ny o
      h ohr mhos - his is work in progrss.

``srip-sqnc``

   rmov h sqnc rom  rs in  bm-i. No h
   sripping h sqnc wi so rmov h qiy scors.
   Sripping is no rvrsib i h r nms r no niq.

``srip-qiy``

   rmov h qiy scors rom  rs in  bm-i.
   Sripping is no rvrsib i h r nms r no niq.

``s-sqnc``

   s h sqnc n qiy scors in h bm i o som mmy
   vs ('A' or sqnc, 'F' or qiy which is  vi scor in
   mos sq ncoings. Ncssry or som oos h cn no work
   wih bm-is wiho sqnc.

``nsrip``

    sqnc n qiy scors bck o  bm i. Rqirs 
   :rm:`sq` orm i wih h sqncs n qiy scors
   o insr.

``ns-nmpp-mpq``

   ss h mpping qiy o nmpp rs o 0.

``kp-irs-bs``

   kp ony h irs bs o rs so h r coning oos wi
   ony consir h irs bs in h cons

``ownsmp-sing``

   gnrs  ownsmp :rm:`bm` i by rnomy sbsmping
   rs rom  sing n :rm:`bm` i. Th ownsmping
   rins mimpping rs. Th s o his rqirs ownsmping
   prmr o b s n opiony rnoms.

``ownsmp-pir``

   gnrs  ownsmp :rm:`bm` i by rnomy sbsmping
   rs rom  pir n :rm:`bm` i. Th ownsmping
   rins mimpping rs. Th s o his rqirs ownsmping
   prmr o b s n opiony rnoms.

``-sqnc-rror``

     crin mon o rnom rror o r sqncs. This mho
   picks  crin proporion o posiions wihin  r's sqnc
   n rs h ncoi o  rnomy chosn rniv. Th
   mo is niv n ppis niorm probbiiis or posiions n
   ncois. Th mho os no p bs qiis, h
   ignmn n h NM g. As  rs, rror rs h r
   comp vi h NM g wi b nc. Th rror r is s
   by --rror-r.

By , h scrip works rom sin n ops o so.

Usg
-----

For xmp::

   cg bm2bm --mhoir --ir-mhompp < in.bm > o.bm

wi rmov  nmpp rs rom h bm-i.

Exmp or rnning ownsmp::

cg bm2bm --mhoownsmp-pir --ownsmp30000
--rnoms1 -L o.og < Pir.bm > o.bm

Typ::

   cg bm2bm --hp

or commn in hp.

Commn in opions
--------------------

'''

impor os
impor sys
impor mpi
impor shi
impor rnom
impor pysm
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor iroos
impor mh

rom cg.BmToos.bmoos impor bm2bm_ir_bm, SNH


css SbsBm(objc):

    ''' bs css or prorming ownsmping on sing n
    pir bm i

    A icionry o h r nms is m n hn n
    rry mching h nmbrs o rs in h icionry wi b
    cr. This rry conins 1s h mch h nmbr o rs
    o b ownsmp o n 0s o i o h rs o h rry.

    Th r is kp i h v in h icionry mchs 1. Th
    r is yi n ir ovr o proc h op bm.

    Th scrip wi hn mimpping rs.
    '''

     __ini__(s, ini, ownsmp, pir_nNon,
                 sing_nNon, rnom_sNon):

        s.pysm_in1, s.pysm_in2  iroos.(ini)
        s.ownsmp  ownsmp
        s.pir_n  pir_n
        s.sing_n  sing_n
        s.rnom_s  rnom_s

     is_o_rs(s, pirNon):

        '''
        This wi cr  icionry o niq rs in h bm
        '''
        r_is  []

        or r in s.pysm_in1:
            i pir is Tr:
                i r.is_propr_pir:
                    r_is.ppn(r.qnm)
            s:
                r_is.ppn(r.qnm)

        rrn sor(s(r_is))

     ownsmp_pir(s):

        '''
        This ncion wi ownsmp  pir bm i.
        I wi rin mimpping rs i hy hv no bn
        pr-ir
        '''
        i s.rnom_s is no Non:
            rnom.s(s.rnom_s)

        coc_is  s.is_o_rs(pirTr)
        r_is  rnom.smp(coc_is, s.ownsmp)

        i s.ownsmp  n(coc_is):
            E.wrn('''Th ownsmp rs is q o h
            nmbr o niq rs''')
            or r in s.pysm_in2:
                yi r

        # yi r i i is in r_is, i h r is mimpp
        # hn  mimpping rs wi b yi
        s:
            E.wrn('''Mimping rs hv bn c n hs wi
            b op o h in bm i''')
            or r in s.pysm_in2:
                i r.qnm in r_is:
                    yi r

     ownsmp_sing(s):

        '''
        This ncion wi ownsmp  sing bm i.
        I wi rin mimpping rs i no pr-ir
        '''
        i s.rnom_s is no Non:
            rnom.s(s.rnom_s)

        coc_is  s.is_o_rs(pirFs)
        r_is  rnom.smp(coc_is, s.ownsmp)

        i s.ownsmp  n(coc_is):
            E.wrn('''Th ownsmp rs is q o h
            nmbr o niq rs''')
            or r in s.pysm_in2:
                yi r

        # yi r i i is in r_is, i h rs is mimpp
        # hn  mimpping rs wi b yi
        s:
            E.wrn('''Mimping rs hv bn c n hs wi
            b op o h in bm i''')
            or r in s.pysm_in2:
                i r.qnm in r_is:
                    yi r


 procss_bm(ini, oi, opions):
    
    i "ir" in opions.mhos:
        i "rmov-is" in opions.ir_mhos or "kp-is" in opions.ir_mhos:

            i  ini.ch(ni_oTr)
            c  E.Conr()
            i "rmov-is" in opions.ir_mhos:
                or r in i:
                    c.inp + 1
                    i r.qry_nm in ir_qry_nms:
                        c.skipp + 1
                        conin
                    oi.wri(r)
                    c.op + 1
            i "kp-is" in opions.ir_mhos:
                or r in i:
                    c.inp + 1
                    i r.qry_nm no in ir_qry_nms:
                        c.skipp + 1
                        conin
                    oi.wri(r)
                    c.op + 1

            E.ino("cgory\cons\ns\n"  c.sTb())
        s:
            rmov_mismchs, coor_mismchs  Fs, Fs

            i "NM" in opions.ir_mhos:
                rmov_mismchs  Tr

            i "CM" in opions.ir_mhos:
                rmov_mismchs  Tr
                coor_mismchs  Tr

            i "min-ngh" in opions.ir_mhos n opions.minimm_r_ngh  0:
                ris VError("ps spciy --minimm-r-ngh whn sing "
                                 "--ir-mhomin-r-ngh")

            i "min-vrg-bs-qiy" in opions.ir_mhos n opions.minimm_vrg_bs_qiy  0:
                ris VError("ps spciy --min-vrg-bs-qiy whn "
                                 "sing --ir-mhomin-vrg-bs-qiy")

            i rmov_mismchs:
                i no opions.rrnc_bm:
                    ris VError(
                        "rqiring rrnc bm i or rmoving by "
                        "mismchs")

                pysm_r  pysm.AignmnFi(opions.rrnc_bm, "rb")
            s:
                pysm_r  Non

            # ir n gs r h opposi wy ron
            c  bm2bm_ir_bm(
                ini, oi, pysm_r,
                rmov_nonniq"niq" in opions.ir_mhos,
                rmov_niq"non-niq" in opions.ir_mhos,
                rmov_conigsNon,
                rmov_nmpp"mpp" in opions.ir_mhos,
                rmov_mismchsrmov_mismchs,
                ir_rror_ropions.rror_r,
                coor_mismchscoor_mismchs,
                minimm_r_nghopions.minimm_r_ngh,
                minimm_vrg_bs_qiyopions.minimm_vrg_bs_qiy)

            opions.sog.wri("cgory\cons\ns\n"  c.sTb())
    s:

        # s p h moiying irors
        i  ini.ch(ni_oTr)

         nop(x):
            rrn Non
        # ncion o chck i procssing sho sr
        pr_chck_  nop

        i "ns-nmpp-mpq" in opions.mhos:
             ns_nmpp_mpq(i):
                or r in i:
                    i r.is_nmpp:
                        r.mpq  0
                    yi r
            i  ns_nmpp_mpq(i)

        i "s-sqnc" in opions.mhos:
             s_sqnc(i):
                or r in i:
                    # cn' g  ngh o nmpp rs
                    i r.is_nmpp:
                        r.sq  "A"
                        r.q  "F"
                    s:
                        r.sq  "A" * r.inrr_ngh
                        r.q  "F" * r.inrr_ngh

                    yi r
            i  s_sqnc(i)

        i "srip-sqnc" in opions.mhos or "srip-qiy" in \
           opions.mhos:
             srip_sqnc(i):
                or r in i:
                    r.sq  Non
                    yi r

             chck_sqnc(rs):
                i rs[0].sq is Non:
                    rrn 'no sqnc prsn'
                rrn Non

             srip_qiy(i):
                or r in i:
                    r.q  Non
                    yi r

             chck_qiy(rs):
                i rs[0].q is Non:
                    rrn 'no qiy inormion prsn'
                rrn Non

             srip_mch(i):
                or r in i:
                    ry:
                        nm  r.op('NM')
                    xcp KyError:
                        nm  1
                    i nm  0:
                        r.sq  Non
                    yi r

            i opions.srip_mho  "":
                i "srip-sqnc" in opions.mhos:
                    i  srip_sqnc(i)
                    pr_chck_  chck_sqnc
                i "srip-qiy" in opions.mhos:
                    i  srip_qiy(i)
                    pr_chck_  chck_qiy
            i opions.srip_mho  "mch":
                i  srip_mch(i)

        i "nsrip" in opions.mhos:
             biRDicionry(inm):
                i no os.ph.xiss(inm):
                    ris OSError("i no on: s"  inm)
                sqi  pysm.FsxFi(inm)
                sq2sqnc  {}
                or x in sqi:
                    i x.nm in sq2sqnc:
                        ris VError(
                            "r s pic - cn no nsrip"  x.nm)

                    sq2sqnc[x.nm]  (x.sqnc, x.qiy)
                rrn sq2sqnc

            i no opions.sq_pir1:
                ris VError(
                    "ps sppy sq i(s) or nsripping")
            sq2sqnc1  biRDicionry(opions.sq_pir1)
            i opions.sq_pir2:
                sq2sqnc2  biRDicionry(opions.sq_pir2)

             nsrip_npir(i):
                or r in i:
                    r.sq, r.q  sq2sqnc1[r.qnm]
                    yi r

             nsrip_pir(i):
                or r in i:
                    i r.is_r1:
                        r.sq, r.q  sq2sqnc1[r.qnm]
                    s:
                        r.sq, r.q  sq2sqnc2[r.qnm]
                    yi r

            i opions.sq_pir2:
                i  nsrip_pir(i)
            s:
                i  nsrip_npir(i)

        i "s-nh" in opions.mhos:
            i  SNH(i)

        # kp irs bs o rs by chnging h cigrsring o
        # '1M' n, in rs mpping o h rvrs srn,
        # chngs h pos o n - 1
        # Ns o b rcor o mk i mor gnr
        # (s bs, mipoin, ..)
        i "kp_irs_bs" in opions.mhos:
             kp_irs_bs(i):
                or r in i:
                    i r.is_rvrs:
                        r.pos  r.n - 1
                        r.cigrsring  '1M'
                    i no r.is_nmpp:
                        r.cigrsring  '1M'
                    yi r
            i  kp_irs_bs(i)

        # r irs r n chck i procssing sho conin
        # ony possib whn no working rom sin
        # Rcoring: s cch o so o  pr-chck or
        # sin inp.
        i no ini.is_srm:
            # g irs r or chcking pr-coniions
            irs_rs  is(ini.h(1))

            msg  pr_chck_(irs_rs)
            i msg is no Non:
                i opions.orc:
                    E.wrn('proccssing conins, hogh: s'  msg)
                s:
                    E.wrn('procssing no sr: s'  msg)
                    rrn

        i "ownsmp-sing" in opions.mhos:

            i no opions.ownsmp:
                ris VError("Ps provi ownsmp siz")

            s:
                own  SbsBm(inii,
                                 ownsmpopions.ownsmp,
                                 pir_nNon,
                                 sing_nTr,
                                 rnom_sopions.rnom_s)
                i  own.ownsmp_sing()

        i "ownsmp-pir" in opions.mhos:

            i no opions.ownsmp:
                ris VError("Ps provi ownsmp siz")

            s:
                own  SbsBm(inii,
                                 ownsmpopions.ownsmp,
                                 pir_nTr,
                                 sing_nNon,
                                 rnom_sopions.rnom_s)
                i  own.ownsmp_pir()

        i "-sqnc-rror" in opions.mhos:
             _sqnc_rror(i):
                rror_r  opions.rror_r
                mp_nc2vr  {"A": "CGT",
                               "C": "AGT",
                               "G": "ACT",
                               "T": "ACG"}
                or r in i:
                    sqnc  is(r.qry_sqnc)
                    qs  r.qry_qiis
                    npos  in(mh.oor(n(sqnc) * rror_r))
                    posiions  rnom.smp(rng(n(sqnc)), npos)
                    or pos in posiions:
                        ry:
                              mp_nc2vr[sqnc[pos]]
                        xcp KyError:
                            conin
                        sqnc[pos]  [rnom.rnin(0, n() - 1)]

                    r.qry_sqnc  "".join(sqnc)
                    r.qry_qiis  qs
                    yi r

            i  _sqnc_rror(i)

        # conin procssing i n
        or r in i:
            oi.wri(r)

                    
 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mhos", s"mhos", yp"choic",
                      cion"ppn",
                      choics("ir",
                               "kp-irs-bs",
                               "s-nh",
                               "s-sqnc",
                               "srip-sqnc",
                               "srip-qiy",
                               "nsrip",
                               "ns-nmpp-mpq",
                               "ownsmp-sing",
                               "ownsmp-pir",
                               "-sqnc-rror"),
                      hp"mhos o ppy []")

    prsr._rgmn("--srip-mho", s"srip_mho", yp"choic",
                      choics("", "mch"),
                      hp"in which sqncs/qiis o srip. "
                      "mch mns h sripping ony ppis o nris "
                      "wiho mismchs (rqirs NM g o b prsn). "
                      "[]")

    prsr._rgmn("--ir-mho", s"ir_mhos",
                      cion"ppn", yp"choic",
                      choics('NM', 'CM',
                               "mpp", "niq", "non-niq",
                               "rmov-is",
                               "kp-is",
                               "rror-r",
                               "min-r-ngh",
                               "min-vrg-bs-qiy"),
                      hp"ir mho o ppy o rmov ignmns "
                      "rom  bm i. Mip mhos cn b sppi "
                      "[]")

    prsr._rgmn("--rrnc-bm-i", s"rrnc_bm",
                      yp"sring",
                      hp"bm-i o ir wih []")

    prsr._rgmn("--orc-op", s"orc", cion"sor_r",
                      hp"orc procssing. Som mhos sch "
                      "s srip/nsrip wi sop procssing i "
                      "hy hink i no ncssry "
                      "[]")

    prsr._rgmn("--op-sm", s"op_sm", cion"sor_r",
                      hp"op in sm orm []")

    prsr._rgmn(
        "--irs-sq-i", "-1", s"sq_pir1", yp"sring",
        hp"sq i wih r inormion or irs "
        "in pir or npir. Us or nsripping sqnc "
        "n qiy scors []")

    prsr._rgmn(
        "--scon-sq-i", "-2", s"sq_pir2", yp"sring",
        hp"sq i wih r inormion or scon "
        "in pir. Us or nsripping sqnc "
        "n qiy scors  []")

    prsr._rgmn(
        "--ownsmp", s"ownsmp",
        yp"in",
        hp"Nmbr o rs o ownsmp o")

    prsr._rgmn(
        "--inm-r-is", s"inm_r_is",
        yp"sring",
        hp"Finm wih is o rs o ir i 'kp-is' or 'rmov-is' "
        "ir mho is chosn []")

    prsr._rgmn(
        "--rror-r", s"rror_r",
        yp"o",
        hp"rror r o s s ir. Rs wih n rror r "
        "highr hn h hrsho wi b rmov []")

    prsr._rgmn(
        "--minimm-r-ngh", s"minimm_r_ngh",
        yp"in",
        hp"minimm r ngh whn iring []")

    prsr._rgmn(
        "--minimm-vrg-bs-qiy", s"minimm_vrg_bs_qiy",
        yp"o",
        hp"minimm vrg bs qiy whn iring []")

    prsr.s_s(
        mhos[],
        op_smFs,
        rrnc_bmNon,
        ir_mhos[],
        srip_mho"",
        orcFs,
        sq_pir1Non,
        sq_pir2Non,
        ownsmpNon,
        rnom_sNon,
        inm_r_isNon,
        rror_rNon,
        minimm_r_ngh0,
        minimm_vrg_bs_qiy0,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.sin ! sys.sin:
        bmi  opions.sin.nm
    i rgs:
        bmi  rgs[0]
        i n(rgs) > 1:
            ris VError("mip bm is provi in rgmns")
    s:
        bmi  "-"
        
    i "rmov-is" in opions.ir_mhos or "kp-is" in opions.ir_mhos:
        i "rmov-is" in opions.ir_mhos n "kp-is" in opions.ir_mhos:
            ris VError("i is no possib o spciy rmov-is n kp-is")

        wih iooos.opn_i(opions.inm_r_is) s in:
            ir_qry_nms  s([x.srip() or x in in.rins() i no x.srswih("#")])
        E.ino("r qry_sqnc ir is wih {} r nms".orm(n(ir_qry_nms)))

    i "rror-r" in opions.ir_mhos n no opions.rror_r:
        ris VError("iring by rror-r rqirs --rror-r o b s")

    i "-sqnc-rror" in opions.mhos n no opions.rror_r:
        ris VError("---rror-r rqirs --rror-r o b s")

    E.ino('procssing s'  bmi)
    i bmi ! "-" n iooos.is_mpy(bmi):
        E.wrn('ignoring mpy i s'  bmi)
        E.sop()
        rrn

    i opions.so ! sys.so:
        op_bmi  opions.so.nm
    s:
        op_bmi  "-"
        i opions.sog  sys.so:
            ris VError("rirc og-srm o i (--og) i oping o so")
        
    i opions.op_sm:
        op_mo  "wh"
    s:
        op_mo  "wb"

    # ring bm rom sin os no work wih ony h "r" g
    wih pysm.AignmnFi(bmi, "rb") s pysm_in:
        wih pysm.AignmnFi(op_bmi, op_mo,
                                 mppysm_in) s pysm_o:
            procss_bm(pysm_in, pysm_o, opions)

    # wri oor n op bnchmrk inormion.
    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
