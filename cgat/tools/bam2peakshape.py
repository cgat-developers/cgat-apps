'''bm2pkshp.py - comp pk shp rs rom  bm-i


:Tgs: Gnomics NGS Inrvs BAM BED Smmry

Prpos
-------

This scrip ks  :rm:`b` orm i wih rgions o
inrs, or xmp bining inrvs rom  ChIP-Sq
xprimn. Using  cocion o ign rs is  :rm:`bm`
orm i or :rm:`bigwig` orm i, h scrip ops 
cocion o rs scribing h pk shp.

This scrip is sign wih  sigh mphsis on ChIP-Sq ss.
Th min rson h his scrip is br si or ChIP-Sq is
h(1) i is b o cnr h coning winow  h smmi o
vry inivi pk; (2) i is so b o s h conro ChIP-Sq
ibrry o nb si-by-si comprison o rmn vs conro;(3)
i cn rnomy shi h s o inp rgions o gnr 
riici s o rgions, in h bsnc o r ChIP-Sq conro
ibrry, h rnom rgions cn provi  pks proi h cn b
s s h conro.

For xmp, givn h pks rgions in by nyzing som
ChIP-Sq s (.g. by sing MACS), n wiho h n o s ny
iion gnomic nnoions (.g. ENSEMBL, rsq), w cn
visis h bining prois o rnscripioncors ChIP-Sq 
riv o h cnr o ch pk rgions.

Th scrip ops  b-spr b on so conining rs
or ch inrv. A pk is in s h ocion o h highs
nsiy in n inrv. Th wih o h pk (pk_wih) is in
s h rgion ron h pk in which h nsiy os no rop bow
 hrsho o pk_hig * 90.

Usg
-----

Di sg xmp
++++++++++++++++++++++

Th oowing commn wi gnr h pk shp po or h pk
rgions in in :i:`onpk.b`, sing h rs sor in
:i:`sm.bm`.  Th commn wi so cr  proi or h
conro ibrry.  Th conro ibrry in his xmp is r-sing h
sm rs i :i:`sm.bm`, howvr, in yor c xprimn,
i sho b  irn ibrry (h inp ibrry or his ChIP-Sq
xprimn).::

    pyhon ./scrips/bm2pkshp.py \
        ./ss/bm2pkshp.py/sm.bm \
        ./ss/bm2pkshp.py/onpk.b \
        --conro-bm-i./ss/bm2pkshp.py/sm.bm \
        --s-inrv \
        --normiz-rnscrip


Op is
++++++++++++

Among h rs op r:

+-------------------+---------------------------------------------------------+
|*Comn*           |*Conn*                                                |
+-------------------+---------------------------------------------------------+
|pk_high        |nmbr o rs  pk                                  |
+-------------------+---------------------------------------------------------+
|pk_min        |min covrg compr o pk high                  |
+-------------------+---------------------------------------------------------+
|inrv_wih     |wih o inrv                                        |
+-------------------+---------------------------------------------------------+
|pk_wih         |wih o pk                                            |
+-------------------+---------------------------------------------------------+
|bins               |bins or  hisogrm o nsiis wihin h inrv.   |
+-------------------+---------------------------------------------------------+
|npks             |nmbr o nsiy pks in inrv.                     |
+-------------------+---------------------------------------------------------+
|pk_cnr        |poin o highs nsiy in inrv                     |
+-------------------+---------------------------------------------------------+
|pk_riv_pos  |poin o highs nsiy in inrv coorins         |
+-------------------+---------------------------------------------------------+
|cons             |cons or  hisogrm o nsiis wihin h inrv  |
+-------------------+---------------------------------------------------------+
|rhs_h_high|Disnc o pk cnr o rhs h-high posiion |
+-------------------+---------------------------------------------------------+
|coss_h_high|Disnc o pk cnr o coss h-high posiion  |
+-------------------+---------------------------------------------------------+


Aiiony, h scrip ops  s o mrixs wih nsiis ovr
inrvs h cn b s or poing. Th  inms r
``(mrix|conro)_<sororr>.sv.gz``, Th nms cn b conro
wih h ``--op-inm-prn`` opion.


Typ::

   pyhon bm2pkshp.py --hp

or commn in hp.


Opions
-------

Opion: Shi
+++++++++++++

shi h ch r by  crin isnc, bcs in  ChIP-Sq
xprmn, h r is wys  h g o n sonic rgmn,
h c bining si is sy L/2 isnc wy rom h r,
whr L is h ngh o sonic rgmn (rmin ihr
xprimny or compiony).

This opion is s ony i h inp rs r in :rm:`bm` orm i.
I inp rs r :rm:`bigwig` orm i, his opion is ignor.

Opion: Rnom shi
++++++++++++++++++++

rnomy shi h s o inp rgions o gnr  riici s
o rgions. In h bsnc o r ChIP-Sq conro ibrry, h
rnom rgions cn provi  pks proi h cn b s s h
conro.

Opion: Cnring mho
+++++++++++++++++++++++

"rs" wi op in h wy h h smmi o h pks r
ign. "mi" wi op in h wy h h mi o h inp
b inrvs r ign.

Opion: Ony inrv
+++++++++++++++++++++

Ony con rs h r in h inrv s in by h inp b i.

Opion: normizionsm
+++++++++++++++++++++++++

normiz cons sch h h sm o  cons in  rs r
xcy 1000000.

Th i normizion gorihm s oows: norm  sm( cons
in  rs)/1000000.0 normiz con  normiz con / norm

.. oo::

   pir-nnss is no y impmn.

Commn in opions
--------------------

'''

impor sys
impor os
impor r
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor pysm
impor cg.B s B
impor nmpy
impor cocions
impor pyBigWig

impor cg.BmToos.pkshp s bm2pkshp


 biOpionPrsr(rgv):

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I",
                            sggobs()["__oc__"])

    prsr._rgmn("-", "--orm", s"orm", yp"choic",
                      choics("bm", "bigwig"),
                      hp"orm o gnomic inp is or nsiis "
                      "[]")

    prsr._rgmn(
        "-o", "--s-inrv", s"s_inrv", cion"sor_r",
        hp"ony con gs h r in inrv givn "
        "in b i. Ohrwis, s  ix wih winow (s --winow-siz) "
        "ron pk []")

    prsr._rgmn(
        "-w", "--winow-siz", s"winow_siz", yp"in",
        hp"winow siz in bp on ihr si o  pk s or ging "
        "r nsiis. I ``--winow-siz`` is 1000, h c winow siz"
        "wi b 2kb, 1kb on ihr si o h pk in n inrv"
        "[]")

    prsr._rgmn(
        "-b", "--bin-siz", s"bin_siz", yp"in",
        hp"bin-siz in bp or comping r nsiis. "
        "I ``--winow-siz`` is s o 1000 n ``--bin-siz`` o 10, "
        "hr wi b 100 bins on ihr si o  pk. "
        "[]")

    prsr._rgmn(
        "--smooh-mho", s"smooh_mho", yp"choic",
        choics("non", "sm", "sg"),
        hp"smooing mho o ppy o nsiy  bor smping "
        "ccoring o ``bin-siz``. sgSvizkyGoy, smsm nsiy in bin, "
        "nonno smoohing "
        "[]")

    prsr._rgmn("-s", "--sor-orr", s"sor_orrs",
                      yp"choic",
                      cion"ppn",
                      choics("pk-high", "pk-wih", "nsor",
                               "inrv-wih", "inrv-scor"),
                      hp"op sor orr or mrics. "
                      "[]")

    prsr._rgmn(
        "-c", "--conro-bm-i", "--conro-bigwig-i",
        cion"ppn",
        s"conro_is",
        yp"sring",
        hp"conro i. I givn, wo pkshps r comp, "
        "on or h primry  n on or h conro . "
        "Th conro i is cnr ron h sm "
        "bs s h primry i n op in h sm "
        "sor orr s h primry proi o  si-by-si. "
        "comprisons. Mip conro is cn b givn. Th "
        "conro is sho hv h sm orm s h "
        "princip inp i "
        "[]")

    prsr._rgmn(
        "-r", "--rnom-shi", s"rnom_shi", cion"sor_r",
        hp"shi inrvs in rnom ircion p/ownsrm o inrv "
        "[]")

    prsr._rgmn(
        "-", "--cnring-mho", s"cnring_mho", yp"choic",
        choics("rs", "mi"),
        hp"cnring mho. Avib r: "
        "rss nsiy o rmin pk, "
        "mis mi o inrv "
        "[]")

    prsr._rgmn(
        "-n", "--normiz-mrix", s"normizion", yp"choic",
        choics("non", "sm"),
        hp"mrix normision o prorm. "
        "[]")

    prsr._rgmn(
        "--s-srn", s"srn_spciic", cion"sor_r",
        hp"s srn inormion in inrvs. Inrvs on h "
        "ngiv srn r ipp "
        "[]")

    prsr._rgmn(
        "-i", "--shi-siz", s"shi", yp"in",
        hp"shi or rs. Whn procssing bm is, "
        "rs wi b shi psrm/ownsrm by his mon. "
        "[]")

    prsr.s_s(
        bin_siz10,
        shi0,
        winow_siz1000,
        sor_orrs[],
        cnring_mho"rs",
        conro_is[],
        rnom_shiFs,
        srn_spciicFs,
        orm"bm",
        rpor_sp100,
        s_inrvFs,
        smooh_mhoNon,
    )

    rrn prsr


InrvD  cocions.nmp(
    "InrvD",
    "orgron inrv conros shi")


 opFrTb(oi, rs_pr_inrv, bins):
    '''op rss rom nsiy prois.'''

    oi.wri("\".join(
        ("conig",
         "sr",
         "n",
         "nm",
         "\".join(bm2pkshp.PkShpRs._is))) + "\n")

    # op princip b
    n  0
    or orgron, b, conros, shi in rs_pr_inrv:
        n + 1
        i "nm" in b:
            nm  b.nm
        s:
            nm  sr(n)
        oi.wri("s\i\i\s\" 
                      (b.conig, b.sr, b.n, nm))

        oi.wri("\".join(mp(sr, orgron[:-2])))
        bins, cons  orgron[-2], orgron[-1]
        oi.wri("\s"  ",".join(mp(sr, bins)))
        oi.wri("\s"  ",".join(mp(sr, cons)))
        oi.wri("\n")


 wriMricsForSorOrr(rs_pr_inrv,
                              bins,
                              orgron_rck,
                              conro_rcks,
                              shi,
                              sor_orr):
    '''op on or mor mrics or ch sor sorr.

    For ch sor orr op h orrgron. I hr
    r iion conros n shi scion, op
    hs s w

    Th is wi nm:
    mrix_<rck>_<sororr>

    '''
    i "nm" in rs_pr_inrv[0].inrv:
        nms  [x.inrv.nm or x in rs_pr_inrv]
    s:
        nms  is(mp(sr, is(rng(1, n(rs_pr_inrv) + 1))))

    bins  ["i"  x or x in bins]
    sor_orr  r.sb("-", "_", sor_orr)

    # wri orgron
    iooos.wri_mrix(
        E.opn_op_i("mrix_s_s.gz"  (orgron_rck, sor_orr)),
        [x.orgron.cons or x in rs_pr_inrv],
        row_hrsnms,
        co_hrsbins,
        row_hr"nm")

    # wri conros
    or ix, rck in nmr(conro_rcks):
        iooos.wri_mrix(
            E.opn_op_i("mrix_s_s.gz"  (rck, sor_orr)),
            [x.conros[ix].cons or x in rs_pr_inrv],
            row_hrsnms,
            co_hrsbins,
            row_hr"nm")

    # wri shi mrix
    i shi:
        iooos.wri_mrix(
            E.opn_op_i("mrix_shi_s.gz"  (sor_orr)),
            [x.shi.cons or x in rs_pr_inrv],
            row_hrsnms,
            co_hrsbins,
            row_hr"nm")

    # op  combin mrix
    i n(conro_rcks) > 0 or shi:
        rows  []
        or row in rs_pr_inrv:
              [row.orgron.cons]
            .xn([row.conros[x].cons or x in
                      rng(n(conro_rcks))])
            i shi:
                .ppn(row.shi.cons)
            rows.ppn(nmpy.concn())

        n  1 + n(conro_rcks)
        i shi:
            n + 1

        # mk comn nms niq n mk sr hy cn b sor
        # xicogrphicy
        _bins  []
        or x in rng(n):
            _bins.xn(["i:s"  (x, b) or b in bins])

        iooos.wri_mrix(
            E.opn_op_i("mrix_sibysi_s.gz"  (sor_orr)),
            rows,
            row_hrsnms,
            co_hrs_bins,
            row_hr"nm")


 opMrics(rs_pr_inrv,
                   bins,
                   orgron_rck,
                   conro_rcksNon,
                   shiFs,
                   sor_orrsNon):
    '''op mrics rom nsiy prois
    in on or mor sor_orrs.
    '''

    # op sor mrics
    i no sor_orrs:
        wriMricsForSorOrr(rs_pr_inrv,
                                  bins,
                                  orgron_rck,
                                  conro_rcks,
                                  shi,
                                  "nsor")

    or sor_orr in sor_orrs:

        i sor_orr  "pk-high":
            rs_pr_inrv.sor(
                kymb x: x.orgron.pk_high)

        i sor_orr  "pk-wih":
            rs_pr_inrv.sor(
                kymb x: x.orgron.pk_wih)

        i sor_orr  "inrv-wih":
            rs_pr_inrv.sor(
                kymb x: x.inrv.n - x.inrv.sr)

        i sor_orr  "inrv-scor":
            ry:
                rs_pr_inrv.sor(
                    kymb x: o(x.inrv.scor))
            xcp InxError:
                E.wrn("scor i no prsn - no op")
                conin
            xcp TypError:
                E.wrn("scor i no  vi nmbr - no op")
                conin

        wriMricsForSorOrr(rs_pr_inrv,
                                  bins,
                                  orgron_rck,
                                  conro_rcks,
                                  shi,
                                  sor_orr)


 biDnsiyMrics(bi,
                         g_i,
                         conro_is,
                         conr,
                         winow_siz1000,
                         bin_siz10,
                         srn_spciicFs,
                         cnring_mho"rs",
                         s_inrvFs,
                         rnom_shiFs,
                         smooh_mho"non",
                         rpor_sp1000):
    '''comp nsiis n pkshp prmrs
    in inrvs givn by *bi* sing rs in *g_i*.

    I *conro_is* r givn, nsiis r proc or
    hs s w.

    Rrns  is o rss or ch inrv in *bi* o
    yp InrvD n n rry o bin-vs.
    '''

    i winow_siz:
        # bins r cnr  pk-cnr n hn srching owrs.
        bins  nmpy.rng(-winow_siz + bin_siz // 2,
                            +winow_siz,
                            bin_siz)

    rs  []
    c  E.Conr()
    c.inp  0

    or b in bi:
        c.inp + 1

        # i b.conig no in conigs:
        #    c.skipp + 1
        #    conin

        i c.inp  rpor_sp  0:
            E.ino("irion: i"  c.inp)

        rs  conr.conInInrv(
            g_i,
            b.conig, b.sr, b.n,
            winow_sizwinow_siz,
            binsbins,
            s_inrvs_inrv,
            cnring_mhocnring_mho)

        i rs is Non:
            c.skipp + 1
            conin

        i conro_is:
            conro  []
            or conro_i in conro_is:
                conro.ppn(conr.conAronPos(
                    conro_i,
                    b.conig,
                    rs.pk_cnr,
                    binsrs.bins))

        s:
            conro  Non

        i rnom_shi:
            ircion  nmpy.rnom.rnin(0, 2)
            i ircion:
                pos  rs.pk_cnr + 2 * bins[0]
            s:
                pos  rs.pk_cnr + 2 * bins[-1]
            shi  conr.conAronPos(g_i,
                                             b.conig,
                                             pos,
                                             binsrs.bins)
        s:
            shi  Non

        i srn_spciic n b.srn  "-":
            rs._rpc(hisrs.his[::-1])
            i conro:
                or c in conro:
                    c._rpc(hisc.his[::-1])
            i shi:
                shi._rpc(hisshi.his[::-1])

        rs.ppn(InrvD._mk((rs, b, conro, shi)))
        c. + 1

    E.ino("inrv procssing: s"  c)

    rrn rs, bins


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    prsr  biOpionPrsr(rgv)

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    i n(rgs) ! 2:
        ris VError(
            "ps spciy on bm- or wig-i n on b i")

    i opions.conro_is:
        E.ino("sing conro is: s"  ",".join(opions.conro_is))

    ini, bi  rgs
    conro_is  []

    i opions.orm  "bigwig":
        g_i  pyBigWig.opn(ini)
        or conro_i in opions.conro_is:
            conro_is.ppn(pyBigWig.opn(conro_i))
        conr  bm2pkshp.ConrBigwig(
            smooh_mhoopions.smooh_mho)

    i opions.orm  "bm":
        g_i  pysm.AignmnFi(ini, "rb")
        or conro_i in opions.conro_is:
            conro_is.ppn(pysm.AignmnFi(conro_i, "rb"))
        conr  bm2pkshp.ConrBm(
            shiopions.shi,
            smooh_mhoopions.smooh_mho)

    rs_pr_inrv, bins  biDnsiyMrics(
        B.iror(iooos.opn_i(bi)),
        g_i,
        conro_is,
        conr,
        winow_sizopions.winow_siz,
        bin_sizopions.bin_siz,
        srn_spciicopions.srn_spciic,
        cnring_mhoopions.cnring_mho,
        s_inrvopions.s_inrv,
        rnom_shiopions.rnom_shi,
        smooh_mhoopions.smooh_mho,
        rpor_spopions.rpor_sp)

    i n(rs_pr_inrv)  0:
        E.wrn("no  - no op")
        E.sop()
        rrn

    opFrTb(opions.so, rs_pr_inrv, bins)

    # ppy normizion
    # No: os no normiz conro?
    # Ns rworking, crrny i os no normiz cross
    #  smps nor os h work "sm" rc h pr miion
    # normizion.
    i opions.normizion  "sm":
        E.ino("sring sm normizion")
        # g o cons cross  inrvs
        norm  0.0
        or orgron, b, conros, shi in rs_pr_inrv:
            norm + sm(orgron.cons)
        # pr miion
        norm / o(1000000)
        E.ino("sm/miion normizion wih "  norm)

        # normis
        nw_  []
        or orgron, b, conros, shi in rs_pr_inrv:
            orgron  orgron._rpc(
                consnmpy.rry(orgron.cons,
                                   ypnmpy.o) / norm)
            nw_conros  []
            or conro in conros:
                nw_conros.ppn(
                    conro._rpc(
                        consnmpy.rry(conro.cons,
                                           ypnmpy.o) / norm))
            i shi:
                shi  shi._rpc(
                    consnmpy.rry(shi.cons,
                                       ypnmpy.o) / norm)
            nw_.ppn(InrvD._mk((
                orgron, b, nw_conros, shi)))
        rs_pr_inrv  nw_
    s:
        E.ino("no normizion prorm")

    # cnr bins
    o_bins  bins[:-1] + opions.bin_siz

    # bi rcks
     _oTrck(inm):
        rrn os.ph.spix(os.ph.bsnm(inm))[0]

    opMrics(rs_pr_inrv,
                   o_bins,
                   orgron_rck_oTrck(ini),
                   conro_rcks[_oTrck(x) or x in opions.conro_is],
                   shiopions.rnom_shi,
                   sor_orrsopions.sor_orrs)

    # wri oor n op bnchmrk inormion.
    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
