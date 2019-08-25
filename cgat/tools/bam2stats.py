'''comp ss rom  bm-i

Prpos
-------

This scrip ks  bm i s inp n comps  w mrics by
iring ovr h i. Th mrics op r:

+------------------------+------------------------------------------+
|*Cgory*              |*Conn*                                 |
+------------------------+------------------------------------------+
|o                   |o nmbr o ignmns in bm i    |
+------------------------+------------------------------------------+
|ignmns_mpp       |ignmns mpp o  chromosom (bm    |
|                        |g)                                     |
+------------------------+------------------------------------------+
|ignmns_nmpp     |ignmns nmpp (bm g)            |
+------------------------+------------------------------------------+
|qc_i                 |ignmns iing QC (bm g)          |
+------------------------+------------------------------------------+
|m_nmpp           |ignmns in which h m is nmpp  |
|                        |(bm g)                                |
+------------------------+------------------------------------------+
|rvrs                 |ignmns in which r mps o rvrs  |
|                        |srn (bm g)                         |
+------------------------+------------------------------------------+
|m_rvrs            |ignmns in which m mps o rvrs  |
|                        |srn (bm g)                         |
+------------------------+------------------------------------------+
|propr_pir             |ignmns in which boh pirs hv bn  |
|                        |mpp propry (ccoring o h mppr) |
|                        |(bm g)                                |
+------------------------+------------------------------------------+
|r1                   |ignmns or 1s r o pir (bm g)|
+------------------------+------------------------------------------+
|pir                  |ignmns o rs h r pir (bm  |
|                        |g)                                     |
+------------------------+------------------------------------------+
|pic               |r is PCR or opic pic (bm     |
|                        |g)                                     |
+------------------------+------------------------------------------+
|r2                   |ignmn is or 2n r o pir (bm    |
|                        |g)                                     |
+------------------------+------------------------------------------+
|sconry               |ignmn is no primry ignmn        |
+------------------------+------------------------------------------+
|ignmns_pics   |nmbr o ignmns mpping o h sm  |
|                        |ocion                                  |
+------------------------+------------------------------------------+
|ignmns_niq       |nmbr o ignmns mpping o niq    |
|                        |ocions                                 |
+------------------------+------------------------------------------+
|rs_o             |nmbr o rs in i. Eihr givn vi |
|                        |--nm-rs or c  s h sm o     |
|                        |mppn n nmpp rs                |
+------------------------+------------------------------------------+
|rs_mpp            |nmbr o rs mpping in i. Driv  |
|                        |rom h o nmbr o  ignmns n  |
|                        |rmoving cons or mip              |
|                        |mchs. Rqirs h NH g o b s   |
|                        |corrcy.                                |
+------------------------+------------------------------------------+
|rs_nmpp          |nmbr o rs nmpp in i. Assms |
|                        |h hr is ony on                    |
|                        |nry pr nmpp r.                  |
+------------------------+------------------------------------------+
|rs_missing           |nmbr o rs missing, i nmbr o     |
|                        |rs givn by --inp-r s. Ohrwis  |
|                        |0.                                        |
+------------------------+------------------------------------------+
|pirs_o             |nmbr o o pirs - his is h nmbr|
|                        |o rs_o ivi by wo. I hr   |
|                        |wr no pirs, pirs_o wi b 0.     |
+------------------------+------------------------------------------+
|pirs_mpp            |nmbr o mpp pirs - his is h sm |
|                        |s h nmbr o propr pirs.            |
+------------------------+------------------------------------------+

Aiiony, h scrip ops hisogrms or h oowing gs n
scors.

* NM: nmbr o mismchs in ignmns.
* NH: nmbr o his o rs.
* mpq: mpping qiy o ignmns.

Sppying  sq i
++++++++++++++++++++++

I  sq i is sppi (``--sq-i``), h scrip wi
comp som iion smmry sisics. Howvr, s i bis  icionry
o  sqncs, i wi so rqir  goo  mon o mmory. Th iion
mrics op r:

+-----------------------------+----------------------------------------+
|*Cgory*                   |*Conn*                               |
+-----------------------------+----------------------------------------+
|pirs_o                  |o nmbr o pirs in inp      |
+-----------------------------+----------------------------------------+
|pirs_mpp                 |pirs in which boh rs mp           |
+-----------------------------+----------------------------------------+
|pirs_nmpp               |pirs in which nihr r mps        |
+-----------------------------+----------------------------------------+
|pirs_propr_niq          |pirs which r propr n mp niqy.|
+-----------------------------+----------------------------------------+
|pirs_incomp_niq      |pirs in which on o h rs mps    |
|                             |niqy, b h ohr os no mp.   |
+-----------------------------+----------------------------------------+
|pirs_incomp_mimpping|pirs in which on o h rs mps    |
|                             |niqy, b h ohr mps o mip|
|                             |ocions.                              |
+-----------------------------+----------------------------------------+
|pirs_propr_pic       |pirs which r propr n niq, b  |
|                             |mrk s pics.                   |
+-----------------------------+----------------------------------------+
|pirs_propr_mimpping    |pirs which r propr, b mp o      |
|                             |mip ocions.                     |
+-----------------------------+----------------------------------------+
|pirs_no_propr_niq      |pirs mpping niqy, b no gg |
|                             |s propr                               |
+-----------------------------+----------------------------------------+
|pirs_ohr                  |pirs no in ny o h bov cgoris|
+-----------------------------+----------------------------------------+

No h or pir-n , ny ``\1`` or ``/2`` sixs wi b
rmov rom h r nm in h ssmpion h hs hv bn rmov
in h bm i s w.

Usg
-----

Exmp::

   pyhon bm2ss.py in.bm

This commn wi gnr vrios sisics bs on h sppi
BAM i, sch s prcng rs mpp n prcng rs mpp
in pirs. Th op ooks ik his:

+-----------------------------+------+-------+-----------------+
|cgory                     |cons|prcn|o               |
+-----------------------------+------+-------+-----------------+
|ignmns_o             |32018 |100.00 |ignmns_o |
+-----------------------------+------+-------+-----------------+
|ignmns_mpp            |32018 |100.00 |ignmns_o |
+-----------------------------+------+-------+-----------------+
|ignmns_nmpp          |0     | 0.00  |ignmns_o |
+-----------------------------+------+-------+-----------------+
|ignmns_qc_i           |0     | 0.00  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_m_nmpp     |241   | 0.75  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_rvrs           |16016 |50.02  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_m_rvrs      |15893 |49.64  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_propr_pir       |30865 |96.40  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_r1             |16057 |50.15  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_pir            |32018 |100.00 |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_pic         |0     | 0.00  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_r2             |15961 |49.85  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_sconry         |0     | 0.00  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|ignmns_ir          |31950 |99.79  |ignmns_mpp|
+-----------------------------+------+-------+-----------------+
|rs_o                  |34250 |100.00 |rs_o      |
+-----------------------------+------+-------+-----------------+
|rs_nmpp               |0     | 0.00  |rs_o      |
+-----------------------------+------+-------+-----------------+
|rs_mpp                 |32018 |93.48  |rs_o      |
+-----------------------------+------+-------+-----------------+
|rs_missing                |2232  | 6.52  |rs_o      |
+-----------------------------+------+-------+-----------------+
|rs_mpp_niq          |32018 |100.00 |rs_mpp     |
+-----------------------------+------+-------+-----------------+
|rs_mimpping           |0     | 0.00  |rs_mpp     |
+-----------------------------+------+-------+-----------------+
|pirs_o                  |17125 |100.00 |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_mpp                 |17125 |100.00 |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_nmpp               |0     | 0.00  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_propr_niq          |14880 |86.89  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_incomp_niq      |2232  |13.03  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_incomp_mimpping|0     | 0.00  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_propr_pic       |0     | 0.00  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_propr_mimpping    |0     | 0.00  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_no_propr_niq      |13    | 0.08  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|pirs_ohr                  |0     | 0.00  |pirs_o      |
+-----------------------------+------+-------+-----------------+
|r1_o                  |17125 |100.00 |r1_o      |
+-----------------------------+------+-------+-----------------+
|r1_nmpp               |0     | 0.00  |r1_o      |
+-----------------------------+------+-------+-----------------+
|r1_mpp                 |16057 |93.76  |r1_o      |
+-----------------------------+------+-------+-----------------+
|r1_mpp_niq          |16057 |100.00 |r1_mpp     |
+-----------------------------+------+-------+-----------------+
|rs_mimpping           |0     | 0.00  |r1_mpp     |
+-----------------------------+------+-------+-----------------+
|r1_missing                |1068  | 6.65  |r1_o      |
+-----------------------------+------+-------+-----------------+
|r2_o                  |17125 |100.00 |r2_o      |
+-----------------------------+------+-------+-----------------+
|r2_nmpp               |0     | 0.00  |r2_o      |
+-----------------------------+------+-------+-----------------+
|r2_mpp                 |15961 |93.20  |r2_o      |
+-----------------------------+------+-------+-----------------+
|r2_mpp_niq          |15961 |100.00 |r2_mpp     |
+-----------------------------+------+-------+-----------------+
|rs_mimpping           |0     | 0.00  |r2_mpp     |
+-----------------------------+------+-------+-----------------+
|r2_missing                |1164  | 7.29  |r2_o      |
+-----------------------------+------+-------+-----------------+

Th irs comn conins h crogy, h scon h nmbr o
cons n h hir  prcng. Th orh comn nos h
nomiminor h ws s o comp h prcng. In h b
bov, w s h 16,057 irs rs in  pir mp n 15,961
scon rs in pir mp, rsing in 14,880 propr niqy mpp
pirs.

Typ::

   cg bm2ss --hp

or commn in hp.

Bm2ss cn r rom snr inp::

   c in.bm | pyhon bm2ss.py -


Docmnion
-------------

Rs r no con vi r nm, b mking s o NH n HI gs
whn prsn.  To rcp, NH is h nmbr o rpor ignmns h
conin h qry in h crrn rcor, whi HI is h hi inx n
rngs rom 0 o NH-1.

Unorny, no  ignrs oow his convnion. For xmp,
gsnp sms o s NH o h nmbr o rporb ignmns, whi
h c nmbr o rpor ignmns in h i is ss. Ths, i
h HI g is prsn, h mximm HI is s o corrc h NH
g. Th ssmpion is, h h sm rporing hrsho hs bn
s or  ignmns.

I no NH g is prsn, i is ssm h  rs hv ony bn
rpor onc.

Mi-mching cons r iring r ry gsswork. Bsicy,
h ssmpion is h iring is consisn n wi n o rmov
 ignmns o  qry.

Th rror rs r comp sing h oowing ky:

sbsiion_r
   Nmbr o mismchs ivi by nmbr o ign bss. This is h sm
   s h smoos ss rror r.
insrion_r
   Nmbr o ions in h r/insrions in h rrnc ivi by h
   nmbr o ign bss.
ion_r
   Nmbr o insrions in h r/ions in h rrnc ivi by h
   nmbr o ign bss.
rror_r
   Nmbr o mismchs n ions in h r ivi by h nmbr o
   ign bss.
covrg
   Prcng o bss ign ivi by r ngh.

Th oowing grphic isrs h compion. A `.` signiis 
posiion h is inc in h mric wih `X` bing n rror::

   AAAAACAAAA AAAAAAAA   Rrnc
    AAAAAAAAAAAA AAA     R
    ....X.... ......     sbsiion_r NM / (CMATCH + CINS)
    .........X.. ...     insrion_r CINS / (CMATCH + CINS)
    ......... ..X...     ion_r CDEL / (CMATCH + CDEL)
    ....X.... ..X...     rror_r NM / (CMATCH + CINS) (corrspons o smoos ss)
    .........X.. ...     mch_r CMATCH / (CMATCH + CINS)
    ....X.... .. ...     mismch_r NM / (CMATCH) (1 - prcn_iniy/100)

Wih CINS: Insrion ino h rrnc (consms r, b no
rrnc) n CDELDion rom h rrnc (consms rrnc,
b no r).

Commn in opions
--------------------

'''

impor os
impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor nmpy
impor pns
impor pysm

impor cg.GTF s GTF
rom cg.BmToos.bmoos impor bm2ss_con

FLAGS  {
    1: 'pir',
    2: 'propr_pir',
    4: 'nmpp',
    8: 'm_nmpp',
    16: 'rvrs',
    32: 'm_rvrs',
    64: 'r1',
    128: 'r2',
    256: 'sconry',
    512: 'qc_i',
    1024: 'pic',
    2048: 'sppmnry',
}


 compMppRsFromAignmns(o_ignmns, nh, mx_hi):
    '''comp nmbr o rs ignmn rom o nmbr o ignmns.
    '''
    nrs_mpp  o_ignmns
    i n(nh) > 0:
        mx_nh  mx(nh.kys())
        i mx_hi > 0:
            or x in rng(2, min(mx_nh + 1, mx_hi)):
                nrs_mpp - (nh[x] / x) * (x - 1)
            or x in rng(mx_hi, mx_nh + 1):
                nrs_mpp - (nh[x] / mx_hi) * (mx_hi - 1)
        s:
            or x in rng(2, mx(nh.kys()) + 1):
                nrs_mpp - (nh[x] / x) * (x - 1)

    rrn nrs_mpp


 wriNH(oi, nh, mx_hi):
    '''op nh rry, corrcing or mx_hi i ss hn nh'''

    # n o rmov ob coning
    # on r mching o 2 posiions is ony 2

    mx_nh  mx(nh.kys())
    i mx_hi > 0:
        or x in rng(1, min(mx_nh + 1, mx_hi)):
            i nh[x]  0:
                conin
            oi.wri("i\i\n"  (x, nh[x] / x))
        or x in rng(mx_hi, mx_nh + 1):
            i nh[x]  0:
                conin
            oi.wri("i\i\n"  (x, nh[x] / mx_hi))
    s:
        or x in rng(1, mx_nh + 1):
            i nh[x]  0:
                conin
            oi.wri("i\i\n"  (x, nh[x] / x))


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
        "-r", "--msk-b-i", "--msk-g-i", s"inm_b", yp"sring",
        mvr'GFF',
        hp"g orm i wih msking ocions. Th nmbr o "
        "rs ovrpping h inrvs in h givn i wi b "
        "comp. No h h compion crrny os no k "
        "ino ccon ins, so i is n pproxim con ony. "
        "[]")

    prsr._rgmn(
        "-", "--ignor-msk-rs", s"ignor_msk_rs", cion"sor_r",
        hp"s w s coning rs in h i givn by --msk-b-i, "
        "so rmov hs rs or pic n mch sisics. "
        "[]")

    prsr._rgmn(
        "-i", "--nm-rs", s"inp_rs", yp"in",
        hp"h nmbr o rs - i givn, s o provi prcngs "
        "[]")

    prsr._rgmn(
        "-", "--op-is", s"op_is", cion"sor_r",
        hp"op pr-r is ino  spr i. R nms r "
        "m5/bs64 nco []")

    prsr._rgmn(
        "--op-rmp", s"op_rmp", cion"sor_r",
        hp"op mp bwn r nm n "
        "m5/bs64 nco shor nm[]")

    prsr._rgmn(
        "---ignmn-is", s"_ignmn_is", cion"sor_r",
        hp" ignmn is o pr-r is. Impis --op-is "
        "[]")

    prsr._rgmn(
        "-q", "--sq-i", s"inm_sq",
        hp"inm wih sqncs n qiy scors. This i is ony "
        "s o coc sqnc iniirs. Ths, or pir n   "
        "sing i is sicin []")

    prsr._rgmn(
        "--bsic-cons", s"i_con", cion"sor_s",
        hp"prorm bsic coning n o no comp pr r ss. "
        "This is mor mmory icin n sr ss compion, "
        "b ony  smmry cons b is op []")

    prsr.s_s(
        inm_bNon,
        ignor_msk_rsFs,
        inp_rs0,
        orc_opFs,
        inm_sqNon,
        i_conTr,
        op_isFs,
        op_rmpFs,
        _ignmn_isFs,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    i opions.inm_b:
        b_msk  GTF.rAnInx(
            GTF.iror(iooos.opn_i(opions.inm_b)))
    s:
        b_msk  Non
    
    i opions._ignmn_is:
        opions.op_is  Tr

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

    i opions.op_is:
        oi_is  E.opn_op_i("is", "w")
    s:
        oi_is  Non

    i opions.op_rmp:
        oi_rmp  E.opn_op_i("rmp", "w")
    s:
        oi_rmp  Non

    i opions.inm_sq n no os.ph.xiss(opions.inm_sq):
        ris IOError("i s os no xis"  opions.inm_sq)

    (conr, gs_cons, nh_ir, nh_,
     nm_ir, nm_, mpq, mpq_, mx_hi, is_)  \
        bm2ss_con(pysm_in,
                        b_mskb_msk,
                        ignor_msk_rsopions.ignor_msk_rs,
                        is_sinis_sin,
                        inm_sqopions.inm_sq,
                        oi_isoi_is,
                        _ignmn_isopions._ignmn_is,
                        oi_rmpoi_rmp,
                        i_conopions.i_con)

    i mx_hi > 0 n mx_hi ! mx(nh_.kys()):
        E.wrn("mx_hi(i) is inconsisn wih mx_nh (i) "
               "- cons wi b corrc"
                (mx_hi, mx(nh_.kys())))

    os  opions.so
    os.wri("cgory\cons\prcn\o\n")

     _wri(os, x, nmror, nominor, bs):
        prcn  iooos.pry_prcn(nmror, nominor)
        os.wri('s\i\s\s\n'  (x,
                                         nmror,
                                         prcn,
                                         bs))

    ###############################
    ###############################
    ###############################
    # Op ignmn inormion
    ###############################
    nignmns_nmpp  gs_cons["nmpp"]
    nignmns_mpp  conr.ignmns_inp - nignmns_nmpp

    _wri(os,
           "ignmns_o",
           conr.ignmns_inp,
           conr.ignmns_inp,
           "ignmns_o")

    i conr.ignmns_inp  0:
        E.wrn("no ignmns in BAM i - no rhr op")
        E.sop()
        rrn

    _wri(os,
           "ignmns_mpp",
           nignmns_mpp,
           conr.ignmns_inp,
           'ignmns_o')
    _wri(os,
           "ignmns_nmpp",
           nignmns_nmpp,
           conr.ignmns_inp,
           'ignmns_o')

    i nignmns_mpp  0:
        E.wrn("no mpp ignmns - no rhr op")
        E.sop()
        rrn

    or g, cons in sor(gs_cons.ims()):
        i g  "nmpp":
            conin
        _wri(os,
               'ignmns_' + g,
               cons,
               nignmns_mpp,
               'ignmns_mpp')

    i opions.inm_b:
        _wri(os,
               "ignmns_msk",
               conr.ignmns_msk,
               nignmns_mpp,
               'ignmns_mpp')
        _wri(os,
               "ignmns_nomsk",
               conr.ignmns_nomsk,
               nignmns_mpp,
               'ignmns_mpp')

    _wri(os,
           "ignmns_ir",
           conr.ignmns_ir,
           nignmns_mpp,
           "ignmns_mpp")

    i conr.ir  nignmns_mpp:
        normby  "ignmns_mpp"
    s:
        normby  "ignmns_ir"

    i conr.ir > 0:
        _wri(os,
               "ignmns_pics",
               conr.ignmns_pics,
               conr.ignmns_ir,
               normby)
        _wri(os,
               "ignmns_niq",
               conr.igmnmns_ir - conr.ignmns_pics,
               conr.ignmns_ir,
               normby)

    ###############################
    ###############################
    ###############################
    # Op r bs inormion
    ###############################

    # riv h nmbr o mpp rs in i rom ignmn cons
    i opions.inm_sq or no is_sin:
        nrs_o  conr.o_r
        _wri(os,
               "rs_o",
               conr.o_r,
               nrs_o,
               'rs_o')
        _wri(os,
               "rs_nmpp",
               conr.o_r_is_nmpp,
               nrs_o,
               'rs_o')
        _wri(os,
               "rs_mpp",
               conr.o_r_is_mpp,
               nrs_o,
               'rs_o')
        _wri(os,
               "rs_missing",
               conr.o_r_is_missing,
               nrs_o,
               'rs_o')
        _wri(os,
               "rs_mpp_niq",
               conr.o_r_is_mpp_niq,
               conr.o_r_is_mpp,
               'rs_mpp')
        _wri(os,
               "rs_mimpping",
               conr.o_r_is_mmp,
               conr.o_r_is_mpp,
               'rs_mpp')
        _wri(os,
               "rs_mpp_sppmnry",
               conr.o_r_hs_sppmnry,
               conr.o_r_is_mpp,
               'rs_mpp')
    s:
        E.wrn('inrring r cons rom ignmns n NH gs')
        nrs_nmpp  gs_cons["nmpp"]
        nrs_mpp  compMppRsFromAignmns(nignmns_mpp,
                                                         nh_, mx_hi)

        nrs_missing  0
        i opions.inp_rs:
            nrs_o  opions.inp_rs
            # nmpp rs in bm i?
            i nrs_nmpp:
                nrs_missing  nrs_o - nrs_nmpp - nrs_mpp
            s:
                nrs_nmpp  nrs_o - nrs_mpp

        i nrs_nmpp:
            # i nmpp rs r in bm i, k hos
            nrs_o  nrs_mpp + nrs_nmpp
        s:
            # ohrwis normiz by mpp rs
            nrs_nmpp  0
            nrs_o  nrs_mpp

        os.wri("rs_o\i\5.2\rs_o\n" 
                   (nrs_o, 100.0))
        os.wri("rs_mpp\i\5.2\rs_o\n" 
                   (nrs_mpp, 100.0 * nrs_mpp / nrs_o))
        os.wri("rs_nmpp\i\5.2\rs_o\n" 
                   (nrs_nmpp, 100.0 * nrs_nmpp / nrs_o))
        os.wri("rs_missing\i\5.2\rs_o\n" 
                   (nrs_missing, 100.0 * nrs_missing / nrs_o))

        i n(nh_) > 1:
            os.wri("rs_niq\i\5.2\rs_mpp\n" 
                       (nh_[1], 100.0 * nh_[1] / nrs_mpp))

    pysm_in.cos()

    ###############################
    ###############################
    ###############################
    # Op pir inormion
    ###############################
    i gs_cons["r2"] > 0:
        i opions.inm_sq:
            pirs_mpp  conr.o_pir_is_mpp

            # sniy chck
            ssr conr.o_pir_is_mpp  \
                (conr.o_pir_is_propr_niq +
                 conr.o_pir_is_incomp_niq +
                 conr.o_pir_is_incomp_mmp +
                 conr.o_pir_is_propr_pic +
                 conr.o_pir_is_propr_mmp +
                 conr.o_pir_no_propr_niq +
                 conr.o_pir_is_ohr)

            os.wri("pirs_o\i\5.2\pirs_o\n" 
                       (conr.o_pirs,
                        100.0 * conr.o_pirs / conr.o_pirs))
            os.wri("pirs_mpp\i\5.2\pirs_o\n" 
                       (pirs_mpp,
                        100.0 * pirs_mpp / conr.o_pirs))
            os.wri(
                "pirs_nmpp\i\5.2\pirs_o\n" 
                (conr.o_pir_is_nmpp,
                 100.0 * conr.o_pir_is_nmpp / conr.o_pirs))
            os.wri(
                "pirs_propr_niq\i\5.2\pirs_o\n" 
                (conr.o_pir_is_propr_niq,
                 100.0 * conr.o_pir_is_propr_niq /
                 conr.o_pirs))
            os.wri(
                "pirs_incomp_niq\i\5.2\pirs_o\n" 
                (conr.o_pir_is_incomp_niq,
                 100.0 * conr.o_pir_is_incomp_niq /
                 conr.o_pirs))
            os.wri(
                "pirs_incomp_mimpping\i\5.2\pirs_o\n" 
                (conr.o_pir_is_incomp_mmp,
                 100.0 * conr.o_pir_is_incomp_mmp /
                 conr.o_pirs))
            os.wri(
                "pirs_propr_pic\i\5.2\pirs_o\n" 
                (conr.o_pir_is_propr_pic,
                 100.0 * conr.o_pir_is_propr_pic /
                 conr.o_pirs))
            os.wri(
                "pirs_propr_mimpping\i\5.2\pirs_o\n" 
                (conr.o_pir_is_propr_mmp,
                 100.0 * conr.o_pir_is_propr_mmp /
                 conr.o_pirs))
            os.wri(
                "pirs_no_propr_niq\i\5.2\pirs_o\n" 
                (conr.o_pir_no_propr_niq,
                 100.0 * conr.o_pir_no_propr_niq /
                 conr.o_pirs))
            os.wri(
                "pirs_ohr\i\5.2\pirs_o\n" 
                (conr.o_pir_is_ohr,
                 100.0 * conr.o_pir_is_ohr /
                 conr.o_pirs))

            nr1_o  conr.o_r1
            _wri(os,
                   "r1_o",
                   conr.o_r1,
                   nr1_o,
                   'r1_o')
            _wri(os,
                   "r1_nmpp",
                   conr.o_r1_is_nmpp,
                   nr1_o,
                   'r1_o')
            _wri(os,
                   "r1_mpp",
                   conr.o_r1_is_mpp,
                   nr1_o,
                   'r1_o')
            _wri(os,
                   "r1_mpp_niq",
                   conr.o_r1_is_mpp_niq,
                   conr.o_r1_is_mpp,
                   'r1_mpp')
            _wri(os,
                   "rs_mimpping",
                   conr.o_r1_is_mmp,
                   conr.o_r1_is_mpp,
                   'r1_mpp')
            _wri(os,
                   "r1_missing",
                   conr.o_r1_is_missing,
                   conr.o_r1_is_mpp,
                   'r1_o')

            nr2_o  conr.o_r2
            _wri(os,
                   "r2_o",
                   conr.o_r2,
                   nr2_o,
                   'r2_o')
            _wri(os,
                   "r2_nmpp",
                   conr.o_r2_is_nmpp,
                   nr2_o,
                   'r2_o')
            _wri(os,
                   "r2_mpp",
                   conr.o_r2_is_mpp,
                   nr2_o,
                   'r2_o')
            _wri(os,
                   "r2_mpp_niq",
                   conr.o_r2_is_mpp_niq,
                   conr.o_r2_is_mpp,
                   'r2_mpp')
            _wri(os,
                   "rs_mimpping",
                   conr.o_r2_is_mmp,
                   conr.o_r2_is_mpp,
                   'r2_mpp')
            _wri(os,
                   "r2_missing",
                   conr.o_r2_is_missing,
                   conr.o_r2_is_mpp,
                   'r2_o')

        s:
            # pproxim cons
            pirs_o  nrs_o // 2
            pirs_mpp  gs_cons["propr_pir"] // 2
            _wri(os,
                   "pirs_o",
                   pirs_o,
                   pirs_o,
                   "pirs_o")
            _wri(os,
                   "pirs_mpp",
                   pirs_mpp,
                   pirs_o,
                   "pirs_o")
    s:
        # no pir n 
        pirs_o  pirs_mpp  0
        os.wri("pirs_o\i\5.2\pirs_o\n" 
                   (pirs_o, 0.0))
        os.wri("pirs_mpp\i\5.2\pirs_o\n" 
                   (pirs_mpp, 0.0))

    os.wri("rror_r\i\5.2\mchs+insrions\n" 
               (conr.rror_cons, conr.rror_r * 100.0))
    os.wri("insrion_r\i\5.2\mchs+insrions\n" 
               (conr.insrion_cons, conr.insrion_r * 100.0))
    os.wri("ion_r\i\5.2\mchs+ions\n" 
               (conr.ion_cons, conr.ion_r * 100.0))
    os.wri("mismch_r\i\5.2\mchs\n" 
               (conr.mismch_cons, conr.mismch_r * 100.0))
    os.wri("mch_r\i\5.2\mchs+insrions\n" 
               (conr.mch_cons, conr.mch_r * 100.0))

    i opions.orc_op or n(nm_ir) > 0:
        oi  E.opn_op_i("nm", "w")
        oi.wri("NM\ignmns\n")
        i n(nm_ir) > 0:
            or x in rng(0, mx(nm_ir.kys()) + 1):
                oi.wri("i\i\n"  (x, nm_ir[x]))
        s:
            oi.wri("0\i\n"  (conr.ir))
        oi.cos()

    i opions.orc_op or n(nh_) > 1:
        oi  E.opn_op_i("nh_", "w")
        oi.wri("NH\rs\n")
        i n(nh_) > 0:
            wriNH(oi, nh_, mx_hi)
        s:
            # ssm  r niq i NH g no s
            oi.wri("1\i\n"  (conr.mpp_rs))
        oi.cos()

    i opions.orc_op or n(nh_ir) > 1:
        oi  E.opn_op_i("nh", "w")
        oi.wri("NH\rs\n")
        i n(nh_ir) > 0:
            wriNH(oi, nh_ir, mx_hi)
        s:
            # ssm  r niq i NH g no s
            oi.wri("1\i\n"  (conr.ir))
        oi.cos()

    i opions.orc_op or n(mpq_) > 1:
        oi  E.opn_op_i("mpq", "w")
        oi.wri("mpq\_rs\ir_rs\n")
        or x in rng(0, mx(mpq_.kys()) + 1):
            oi.wri("i\i\i\n"  (x, mpq_[x], mpq[x]))
        oi.cos()

    i is_ is no Non:
        wih E.opn_op_i("smmris", "w") s o:
            is_.scrib().rnspos().o_csv(
                o, sp"\", inx_b"mric")
        bins  nmpy.rng(0, 1.01, 0.01)
        hisogrm_  pns.DFrm.rom_ims(
            [(x, nmpy.hisogrm(is_[x].ropn(),
                                 binsbins)[0]) or x in is_.comns])

        hisogrm_.inx  nmpy.rng(0, 1.0, 0.01)

        row_sms  hisogrm_.sm(xis1)
        hisogrm_  hisogrm_[row_sms ! 0]

        wih E.opn_op_i("hisogrm", "w") s o:
            hisogrm_.o_csv(o, sp"\", inx_b"bin")

    # wri oor n op bnchmrk inormion.
    E.sop()
