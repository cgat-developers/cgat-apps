'''bm2gnproi.py - bi m-gn proi or  s o rnscrips/gns


:Tgs: Gnomics NGS Gnss Inrvs GTF BAM Smmry

Prpos
-------

This scrip ks  :rm:`g` orm i,  shor rs
:rm:`bm` orm i n comps m-gn prois ovr
vrios nnoions riv rom h :rm:`g` i.

A m-gn proi is n bsrc gnomic niy ovr which rs
sor in  :rm:`bm` orm i hv bn con. A m-gn
migh b n iiz kryoic gn (psrm, xonic sqnc,
ownsrm) or ny ohr gnomic nmrk o inrs sch s
rnscripion sr sis.

Th scrip cn b s o visiz bining prois o  chromin
mrk in gn bois, bining o rnscripion cors in promoors or
sqncing bis (.g. 3' bis) in RNA-Sq .

This scrip is sign wih  sigh mphsis on RNA-Sq
ss. For xmp, i ks cr o spic rs, by sing h
CIGAR sring in h BAM i o ccry in ign bss (i
h --bs-ccr is spcii, crrny his r is rn o
by  or sp rsons).

Arnivy, or h prpos o visizing bining prois o
rnscripion cor ChIP-Sq wiho h n o s ny gnomic
nnoions (ENSEMBL, or rsq), yo my so consir sing
:oc:`bm2pkshp`, which is sign wih  sigh mphsis on
Chip-Sq ss. For xmp, :oc:`bm2pkshp` is b o cnr
h coning winow o h smmi o vry inivi pk.
:oc:`bm2pkshp` is so b o: (1) po h conro ChIP-Sq
ibrry o nb si-by-si comprison; (2) rnomiz h givn
rgions o provi  smi-conro.

Usg
-----

Qick sr xmps
++++++++++++++++++++

Th oowing commn wi gnr h gn proi po simir o
Fig 1() in h pbish cg ppr, b sing  s s h is
mch smr n simpr hn h s s or pbishing h cg
ppr. ::

    pyhon ./scrips/bm2gnproi.py
        --bm-i./ss/bm2gnproi.py/mipRsSpicOAInronsAnSconExon.bm
        --g-i./ss/bm2gnproi.py/ongnWihoAnyCDS.g.gz
        --mhognproi
        --rporrgn

In h oowing,  sighy mor invov xmp wi s mor
rs o his scrip. Th oowing commn gnr h gn
proi showing bs ccrcy o psrm (500bp), xons, inrons n
ownsrm(500bp) o  gn mo rom som sr sppi RNA-Sq 
n gns. ::

    pyhon ./scrips/bm2gnproi.py
        --bm-i./rnsq.bm
        --g-i./gns.g.gz
        --mhognproiwihinrons
        --rporrgn
        --xnsion-psrm500
        --rsoion-psrm500
        --xnsion-ownsrm500
        --rsoion-ownsrm500

Th op wi conin r covrg ovr gns. Th proi wi
conin or spr sgmns:

1. h psrm rgion o  gn ( s o b 500bp ),
   (``--xnsion-psrm500``).

2. h rnscrib rgion o  gn. Th rnscrib rgion o vry gn wi
   b sc o 1000 bp (), shrinking ongr rnscrips n
   xpning shorr rnscrips.

3. h inronic rgions o  gn. Ths wi b sc o 1000b ().

4. h ownsrm rgion o  gn (s o b 500bp),
   (``--xnsion-ownsrm500``).



Di xpinion
+++++++++++++++++++++

Th :i:`bm2gnproi.py` scrip rs in  s o rnscrips
rom  :rm:`g` orm i. For ch rnscrip, ovrpping
rs rom h provi :rm:`bm` i r coc. Th cons
wihin h rnscrip r hn mpp ono h m-gn srcr n
cons r ggrg ovr  rnscrips in h :rm:`g` i.

:rm:`Bm` is n o b sor by coorin n inx.

A m-gn srcr hs wo componns - rgions o vrib siz,
sch s xons, inrons, c, which nvrhss hv  ix sr n
n coorin in  rnscrip. Th ohr componn r rgions o
ix wih, sch  rgions o  crin siz psrm or ownsrm
o  nmrk sch s  rnscripion sr si.

Th siz o h ormr css, rgions o vrib siz, cn b vri
wih ``--rsoion`` opions. For xmp, h opion
``--rsoion-psrm-r1000`` wi cr  m-gn wih 
1000bp psrm UTR rgion. UTRs h r rgr wi b comprss,
n UTRs h r smr, wi b srch o i h 1000bp
m-gn UTR rgion.

Th siz o ix-wih rgions cn b s wih ``--xnsion``
opions. For xmp, h opions ``--xnsion-psrm`` wi s
h siz o h psrm xnsion rgion o 1000bp. No h no
scing is rqir whn coning rs owrs h ix-wih
m-gn proi.

Typ::

   pyhon bm2gnproi.py --hp

or commn in hp.

Opions
-------

Th scrip provis  vriy o irn m-gn srcrs i..
gnprois, scb vi sing h opion: (``--mho``).

Prois
++++++++

Dirn prois r ccssib hrogh h ``--mho`` opion. Mip
mhos cn b ppi  h sm im. Whi ``psrm`` n ``ownsrm``
ypicy hv  ix siz, h ohr rgions sch s ``CDS``, ``UTR`` wi b
sc o  common siz.

rproi
    UPSTREAM - UTR5 - CDS - UTR3 - DOWNSTREAM
    gn mos wih UTR. Spr h coing scion rom h non-coing pr.

gnproi
    UPSTREAM - EXON - DOWNSTREAM
    simp xonic gn mos

gnproiwihinrons
    UPSTREAM - EXON - INTRON - DOWNSTREAM

    gn mos conining so inronic sqnc, ony corrc i
    s wih ``--s-bs-ccrcy`` opion.

sprxonproi
    UPSTREAM - FIRST EXON - EXON - LAST EXON - DOWNSTREAM

    gn mos wih h irs n s xons spr o rom 
    ohr xons.  Ony ppicb o gn mos wih > 1 xons.

sprxonproiwihinrons
    UPSTREAM - FIRST EXON - EXON - INTRON - LAST EXON - DOWNSTREAM

    gn mos wih irs n s xons spr o, n incs
     inrons oghr.  Excs gns wih < 2 xons n no inrons.

gnproibsoisncromhrprimn

    UPSTREAM - EXON (bso isnc, s bow) - INTRON (bso
    isnc, s bow) - DOWNSTREAM (h ownsrm o h xons)
    rgion, h scrip cons ovr h mRNA rnscrip ony, skipping
    inrons. Dsign o visiz h 3 prim bis in RNASq ,
    ony corrc i s oghr wih ``--s-bs-ccrcy`` opion.

    bso isnc: In orr o o visiz h 3 prim bis,
    gns r no sppos o b srch o q ngh s i i in
     ohr coning mhos. In his coning mho, w irs s
     ix ngh sing
    ``--xnsion-xons-bso-isnc-opoy``, h scrip wi
    iscr gns shorr hn his ix ngh. For gns (whn 
    h xons sich oghr) ongr hn his ix ngh, h
    scrip wi ony con ovr his ix ngh (  bso
    isnc ) rom hr prim n, ins o comprss h ongr
    gns. Sm gos or bso isnc inron coning.

ssproi
    UPSTREAM - DOWNSTREAM
    rnscripion sr/sop sis

inrvproi
    UPSTREAM - INTERVAL - DOWNSTREAM
    Simir o gnproi, b con ovr h comp spn o h gn
    (incing inrons).

mipoinproi
    UPSTREAM  - DOWNSTREAM
    ggrg ovr mipoin o gn mo


Normizion
+++++++++++++

Normizion cn b ppi in wo sgs o h compion.

Con vcor normizion
^^^^^^^^^^^^^^^^^^^^^^^^^^

Bor ing cons o h m-gn proi, h proi or h
inivi rnscrip cn b normiz. Wiho normizion, highy
xprss gns wi conrib mor o h m-gn proi hn
owy xprss gns.  Normizion cn ssr h ch gn
conribs n q mon.

Normizion is ppi o h vcor o r cons h is comp
or ch rnscrip. Normizion cn b ppi or h who
rnscrip (``o``) or on  pr sgmn bsis pning on h
conr. For xmp, in h gn conr, xons, psrm n
ownsrm sgmns cn b normiz inpnny.

Cons cn b normiz ihr by h mximm or h sm o 
cons in  sgmn or cross h who rnscrip. Normizion is
conro wih h commn in opion ``--normiz-rncrip``. Is
rgmns r:

* ``non``: no normizion
* ``sm``: sm o cons wihin  rgion (xons, psrm, ...).
  Th r nr h crv wi sm o 1 or ch rgion.
* ``mx``: mximm con wihin  rgion (xons,psrm, ...).
* ``o-sm``: sm o cons cross  rgions. Th r
  nr h crv wi sm o 1 or
  h comp rnscrip.
* ``o-mx``: mximm con cross  rgions.

Th opions bov conro h conribion o inivi rnscrips
o  m-gn proi n r hs si or xmp or RNA-Sq .

Th opions bov o no conro or irn r-phs or ny
oc biss. To compr m-gn prois bwn smps,
iion normizion is rqir.

M-gn proi normizion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To nb comprison bwn xprimns, h m-gn proi is
cn b normiz.  Normizion  proi cn hp compring h
shps o prois bwn irn xprimns inpnn o h
nmbr o rs or rnscrips s in h consrcion o h
m-gn proi.

M-gn proi normizion is conro vi h
``--normiz-proi`` opion. Possib normizion r:

* non: no normizion
* r: normiz sch h h r nr h m-gn proi is 1.
* cons: normiz by nmbr o rs (gns,ss) h hv bn con.
* bckgron: normiz wih bckgron (s bow).

A spci normizion is civ wih h ``bckgron`` opion.
Hr, h cons  h  n righ mos rgions r s o
sim  bckgron v or ch rnscrip. Th cons r hn
ivi by his bckgron-v. Th ssmpion is h h m-gn
mo is comp ovr  rg nogh r o inc gnomic
bckgron.

Gns vrss rnscrips
++++++++++++++++++++++++

Th  is o coc rs on  pr-rnscrip
v. Arnivy, h scrip cn mrg  rnscrips o  gn
ino  sing vir rnscrip. No h his vir rnscrip
migh no b  bioogicy psib rnscrip. I is sy br
o provi :i:`bm2gnproi.py` wih  s o rprsniv
rnscrips pr gn in orr o voi p-wighing gns wih
mip rnscrips.

Conro
+++++++

I conro is (chip-sq inp rcks) r sppi, cons in h
conro i cn b s o comp  o-chng.

B n wigg is
++++++++++++++++++++

Th nsiis cn b comp rom :rm:`b` or :rm:`wigg`
orm is. I  :rm:`b` orm i is sppi, i ms
b comprss wih n inx wih :i:`bix`.

.. no::

   Pir-nnss is ignor. Boh ns o  pir-n r r
   r iniviy.


Commn in opions
--------------------

'''

impor os
impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor pysm
impor cg.GTF s GTF
impor nmpy
impor pns
impor pyBigWig

rom cg.BmToos impor gnproi


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
        """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mhos", yp"choic",
                      cion"ppn",
                      choics("gnproi", "ssproi", "rproi",
                               "inrvproi", "mipoinproi",
                               "gnproiwihinrons",
                               "gnproibsoisncromhrprimn",
                               "sprxonproi",
                               "sprxonproiwihinrons",
                               ),
                      hp'conrs o s. Conrs scrib h '
                      'm-gn srcr o s. '
                      'No sing gnproiwihinrons, or '
                      'gnproibsoisncromhrprimn wi '
                      'omicy rn on h --s-bs-ccrcy opion'
                      '[].')

    prsr._rgmn("-b", "--bm-i", "--bi", "--bigwigi",
                      s"inis",
                      mvr"BAM",
                      yp"sring", cion"ppn",
                      hp"BAM/b/bigwig is o s. Do no mix "
                      "irn yps []")

    prsr._rgmn("-c", "--conro-bm-i", s"conrois",
                      mvr"BAM",
                      yp"sring", cion"ppn",
                      hp"conro/inp o s. Sho b o h sm "
                      "yp s h bm/b/bigwig i"
                      " []")

    prsr._rgmn("-g", "--g-i", s"gi", yp"sring",
                      mvr"GTF",
                      hp"GTF i o s. "
                      "[]")

    prsr._rgmn(
        "--normiz-rnscrip",
        s"rnscrip_normizion",
        yp"choic",
        choics("non", "mx", "sm", "o-mx", "o-sm"),
        hp"normizion o ppy on ch rnscrip "
        "proi bor ing o m-gn proi. "
        "[]")

    prsr._rgmn(
        "--normiz-proi",
        s"proi_normizions",
        yp"choic", cion"ppn",
        choics("", "non", "r", "cons", "bckgron"),
        hp"normizion o ppy on m-gn "
        "proi normizion. "
        "[]")

    prsr._rgmn(
        "-r", "--rporr", s"rporr", yp"choic",
        choics("gn", "rnscrip"),
        hp"rpor rss or gns or rnscrips."
        " Whn 'gns` is chosn, xons cross  rnscrips or"
        "  gn r mrg. Whn 'rnscrip' is chosn, cons r"
        " comp or ch rnscrip spry wih ch rnscrip"
        " conribing qy o h m-gn proi."
        " []")

    prsr._rgmn("-i", "--shi-siz", s"shis", yp"in",
                      cion"ppn",
                      hp"shi rs in :rm:`bm` orm i "
                      "bor comping nsiis (ChIP-Sq). "
                      "[]")

    prsr._rgmn("-", "--mrg-pirs", s"mrg_pirs",
                      cion"sor_r",
                      hp"mrg pirs in :rm:`bm` orm "
                      "i bor comping "
                      "nsiis (ChIP-Sq). "
                      "[]")

    prsr._rgmn("-", "--s-bs-ccrcy", s"bs_ccrcy",
                      cion"sor_r",
                      hp"comp nsiis wih bs ccrcy. Th  "
                      "is o ony s h sr n n o h ign rgion "
                      "(RNA-Sq) "
                      "[]")

    prsr._rgmn("-", "--xn", s"xns", yp"in",
                      cion"ppn",
                      hp"xn rs in :rm:`bm` orm i "
                      "(ChIP-Sq). "
                      "[]")

    prsr._rgmn("--rsoion-psrm", s"rsoion_psrm",
                      yp"in",
                      hp"rsoion o psrm rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-ownsrm", s"rsoion_ownsrm",
                      yp"in",
                      hp"rsoion o ownsrm rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-psrm-r",
                      s"rsoion_psrm_r",
                      yp"in",
                      hp"rsoion o psrm UTR rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-ownsrm-r",
                      s"rsoion_ownsrm_r",
                      yp"in",
                      hp"rsoion o ownsrm UTR rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-cs", s"rsoion_cs", yp"in",
                      hp"rsoion o cs rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-irs-xon", s"rsoion_irs",
                      yp"in",
                      hp"rsoion o irs xon in gn, in bp"
                      "[]")

    prsr._rgmn("--rsoion-s-xon", s"rsoion_s",
                      yp"in",
                      hp"rsoion o s xon in gn, in bp"
                      "[]")

    prsr._rgmn("--rsoion-inrons",
                      s"rsoion_inrons", yp"in",
                      hp"rsoion o inrons rgion in bp "
                      "[]")

    prsr._rgmn("--rsoion-xons-bso-isnc-opoy",
                      s"rsoion_xons_bso_isnc_opoy",
                      yp"in",
                      hp"rsoion o xons bso isnc "
                      "opoy in bp "
                      "[]")

    prsr._rgmn("--rsoion-inrons-bso-isnc-opoy",
                      s"rsoion_inrons_bso_isnc_opoy",
                      yp"in",
                      hp"rsoion o inrons bso isnc "
                      "opoy in bp "
                      "[]")

    prsr._rgmn("--xnsion-xons-bso-isnc-opoy",
                      s"xnsion_xons_bso_isnc_opoy",
                      yp"in",
                      hp"xnsion or xons rom h bso "
                      "isnc rom h opoy in bp "
                      "[]")

    prsr._rgmn(
        "--xnsion-inrons-bso-isnc-opoy",
        s"xnsion_inrons_bso_isnc_opoy", yp"in",
        hp"xnsion or inrons rom h bso isnc rom "
        "h opoy in bp []")

    prsr._rgmn(
        "--xnsion-psrm", s"xnsion_psrm", yp"in",
        hp"xnsion psrm rom h irs xon in bp"
        "[]")

    prsr._rgmn(
        "--xnsion-ownsrm", s"xnsion_ownsrm", yp"in",
        hp"xnsion ownsrm rom h s xon in bp"
        "[]")

    prsr._rgmn(
        "--xnsion-inwr", s"xnsion_inwr", yp"in",
        hp"xnsion inwr rom  TSS sr si in bp"
        "[]")

    prsr._rgmn(
        "--xnsion-owr", s"xnsion_owr", yp"in",
        hp"xnsion owr rom  TSS sr si in bp"
        "[]")

    prsr._rgmn("--sc-nk-ngh", s"sc_nks", yp"in",
                      hp"sc nks o (ingr mips o) gn ngh"
                      "[]")

    prsr._rgmn(
        "--conro-cor", s"conro_cor", yp"o",
        hp"cor or normizing conro n orgron . "
        "Comp rom  i no s. "
        "[]")

    prsr._rgmn("--op--prois", s"op__prois",
                      cion"sor_r",
                      hp"kp inivi prois or ch "
                      "rnscrip n op. "
                      "[]")

    prsr._rgmn("--cons-sv-i", s"inp_inm_cons",
                      yp"sring",
                      hp"inm wih con  or ch rnscrip. "
                      "Us his ins "
                      "o rcomping h proi. Us or poing h "
                      "m-gn proi "
                      "rom prviosy comp cons "
                      "[]")

    prsr._rgmn(
        "--bckgron-rgion-bins",
        s"bckgron_rgion_bins",
        yp"in",
        hp"nmbr o bins on ihr n o h proi "
        "o b consir or bckgron m-gn normizion "
        "[]")

    prsr._rgmn("--op-rs",
                      s"rsoion_imgs", yp"in",
                      hp"h op pi or h igr po - wi  o "
                      "[]")

    prsr._rgmn("--img-orm", s"img_orm", yp"sring",
                      hp"Th op orm or h igr po - s o "
                      "[]")                      

    prsr.s_s(
        rmov_rnFs,
        ignor_pirsFs,
        orc_opFs,
        bin_siz10,
        xns[],
        shis[],
        sor[],
        rporr"rnscrip",
        rsoion_cs1000,
        rsoion_inrons1000,
        # 3kb is  goo bnc o sing ong nogh 3 prim bis n no omi
        # oo mny gns. Tim 31h Ag 2013
        rsoion_xons_bso_isnc_opoy3000,
        # inrons is ony or ssss h nois v, hs o on n  ong
        # rgion,  ong rgion hs h si c o omi mor gns. Tim
        # 31h Ag 2013
        rsoion_inrons_bso_isnc_opoy500,
        # xnsion cn simpy js b h sm s rsoion
        xnsion_xons_bso_isnc_opoy3000,
        xnsion_inrons_bso_isnc_opoy500,
        rsoion_psrm_r1000,
        rsoion_ownsrm_r1000,
        rsoion_psrm1000,
        rsoion_ownsrm1000,
        rsoion_irs1000,
        rsoion_s1000,
        # mn ngh o rnscrips: bo 2.5 kb
        xnsion_psrm2500,
        xnsion_ownsrm2500,
        xnsion_inwr3000,
        xnsion_owr3000,
        poTr,
        mhos[],
        inis[],
        conrois[],
        giNon,
        proi_normizions[],
        rnscrip_normizionNon,
        sc_nks0,
        mrg_pirsFs,
        min_insr_siz0,
        mx_insr_siz1000,
        bs_ccrcyFs,
        mrix_orm"sing",
        conro_corNon,
        op__proisFs,
        bckgron_rgion_bins10,
        inp_inm_consNon,
        rsoion_imgsNon,
        img_orm"png",
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    # Kp or bckwrs compbiiy
    i n(rgs)  2:
        ini, g  rgs
        opions.inis.ppn(ini)
        opions.gi  g

    i no opions.gi:
        ris VError("no GTF i spcii")

    i opions.gi  "-":
        opions.gi  opions.sin
    s:
        opions.gi  iooos.opn_i(opions.gi)

    i n(opions.inis)  0:
        ris VError("no bm/wig/b is spcii")

    or mhosRqirsBsAccrcy in [
            "gnproiwihinrons",
            "gnproibsoisncromhrprimn",
    ]:
        # I yo impmn ny mhos h yo o no wn h
        # spic o inrons or xons ppr o b covr by
        # non-xisn rs, i is br yo  hos mhos impy
        # --bs-ccrrcy by  hm hr.
        i mhosRqirsBsAccrcy in opions.mhos:
            opions.bs_ccrcy  Tr

    i opions.rporr  "gn":
        g_iror  GTF._gn_iror(GTF.iror(opions.gi))
    i opions.rporr  "rnscrip":
        g_iror  GTF.rnscrip_iror(GTF.iror(opions.gi))

    # Sc rngconr bs on i yp
    i n(opions.inis) > 0:
        i opions.inis[0].nswih(".bm"):
            bmis  [pysm.AignmnFi(x, "rb") or x in opions.inis]

            i opions.conrois:
                conrois  [pysm.AignmnFi(x, "rb")
                                or x in opions.conrois]
            s:
                conrois  Non

            orm  "bm"
            i opions.mrg_pirs:
                rng_conr  gnproi.RngConrBAM(
                    bmis,
                    shisopions.shis,
                    xnsopions.xns,
                    mrg_pirsopions.mrg_pirs,
                    min_insr_sizopions.min_insr_siz,
                    mx_insr_sizopions.mx_insr_siz,
                    conroisconrois,
                    conro_coropions.conro_cor)

            i opions.shis or opions.xns:
                rng_conr  gnproi.RngConrBAM(
                    bmis,
                    shisopions.shis,
                    xnsopions.xns,
                    conroisconrois,
                    conro_coropions.conro_cor)

            i opions.bs_ccrcy:
                rng_conr  gnproi.RngConrBAMBsAccrcy(
                    bmis,
                    conroisconrois,
                    conro_coropions.conro_cor)
            s:
                rng_conr  gnproi.RngConrBAM(
                    bmis,
                    conroisconrois,
                    conro_coropions.conro_cor)

        i opions.inis[0].nswih(".b.gz"):
            bis  [pysm.Tbixi(x) or x in opions.inis]

            i opions.conrois:
                conrois  [pysm.Tbixi(x)
                                or x in opions.conrois]
            s:
                conrois  Non

            rng_conr  gnproi.RngConrB(
                bis,
                conroisconrois,
                conro_coropions.conro_cor)

        i opions.inis[0].nswih(".bw"):
            wigis  [pyBigWig.opn(x) or x in opions.inis]
            rng_conr  gnproi.RngConrBigWig(wigis)

        s:
            ris NoImpmnError(
                "cn' rmin i yp or s"  sr(opions.inis))

    conrs  []
    or mho in opions.mhos:
        i mho  "rproi":
            conrs.ppn(gnproi.UTRConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_psrm_r,
                opions.rsoion_cs,
                opions.rsoion_ownsrm_r,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm,
            ))

        i mho  "gnproi":
            conrs.ppn(gnproi.GnConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_cs,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm,
                opions.sc_nks))

        i mho  "gnproiwihinrons":
            conrs.ppn(gnproi.GnConrWihInrons(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_cs,
                opions.rsoion_inrons,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm,
                opions.sc_nks))

        i mho  "gnproibsoisncromhrprimn":
            # opions.xnsion_xons_bso_isnc_osrsi,
            # opions.xnsion_inrons_bso_isnc_osrsi,
            # Tim 31h Ag 2013:  possib r or r,  i iv prim
            # bis is o yor inrs.
            # (yo n o cr nohr css). I is no vry iic o
            # riv rom his css, b is no impmn y
            # This r r is sighy irn h TSS proi
            # ry impmn, bcs in his r r inrons r
            # skipp,
            conrs.ppn(
                gnproi.GnConrAbsoDisncFromThrPrimEn(
                    rng_conr, opions.rsoion_psrm,
                    opions.rsoion_ownsrm,
                    opions.rsoion_xons_bso_isnc_opoy,
                    opions.rsoion_inrons_bso_isnc_opoy,
                    opions.xnsion_psrm,
                    opions.xnsion_ownsrm,
                    opions.xnsion_xons_bso_isnc_opoy,
                    opions.xnsion_inrons_bso_isnc_opoy,
                    opions.sc_nks))

        i mho  "ssproi":
            conrs.ppn(gnproi.TSSConr(
                rng_conr,
                opions.xnsion_owr,
                opions.xnsion_inwr))

        i mho  "inrvproi":
            conrs.ppn(gnproi.RgionConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_cs,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm))

        i mho  "mipoinproi":
            conrs.ppn(gnproi.MipoinConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm))

        #  nw mho o spi 1s n s xons o
        # rqirs  rprsniv rnscrip or rch gn
        # g sho b sor gn-posiion
        i mho  "sprxonproi":
            conrs.ppn(gnproi.SprExonConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_irs,
                opions.rsoion_s,
                opions.rsoion_cs,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm))

        i mho  "sprxonproiwihinrons":
            conrs.ppn(gnproi.SprExonWihInronConr(
                rng_conr,
                opions.rsoion_psrm,
                opions.rsoion_irs,
                opions.rsoion_s,
                opions.rsoion_cs,
                opions.rsoion_inrons,
                opions.rsoion_ownsrm,
                opions.xnsion_psrm,
                opions.xnsion_ownsrm))

    # s normizion
    or c in conrs:
        c.sNormizion(opions.rnscrip_normizion)
        i opions.op__prois:
            c.sOpProis(iooos.opn_i(E.g_op_i(c.nm) +
                                                  ".prois.sv.gz", "w"))

    i opions.inp_inm_cons:
        # r cons rom i
        E.ino("ring cons rom s"  opions.inp_inm_cons)
        _cons  pns.r_csv(
            iooos.opn_i(opions.inp_inm_cons),
            sp'\', hr0, inx_co0)

        i n(conrs) ! 1:
            ris NoImpmnError(
                'coning rom mrix ony impmn or 1 conr.')
        # bi conr bs on rrnc conr
        conr  gnproi.UnsgmnConr(conrs[0])
        conrs  [conr]
        gnproi.conFromCons(conrs, _cons)

    s:
        E.ino("sring coning wih i conrs"  n(conrs))
        r_nms  gnproi.conFromGTF(conrs,
                                                 g_iror)

    # op mrics
    i no opions.proi_normizions:
        opions.proi_normizions.ppn("non")
    i "" in opions.proi_normizions:
        opions.proi_normizions  ["non",
                                          "r",
                                          "cons",
                                          "bckgron"]

    or mho, conr in zip(opions.mhos, conrs):
        prois  []
        or norm in opions.proi_normizions:
            # bi mrix, ppy normizion
            proi  conr.gProi(
                normiznorm,
                bckgron_rgion_binsopions.bckgron_rgion_bins)
            prois.ppn(proi)

        or x in rng(1, n(prois)):
            ssr prois[0].shp  prois[x].shp

        # bi  sing mrix o  prois or op
        mrix  nmpy.concn(prois)
        mrix.shp  n(prois), n(prois[0])
        mrix  mrix.rnspos()

        wih iooos.opn_i(E.g_op_i(conr.nm) +
                               ".mrix.sv.gz", "w") s oi:
            oi.wri("bin\rgion\rgion_bin\s\n"  "\".join(
                opions.proi_normizions))
            is  []
            bins  []
            or i, nbins in zip(conr.is, conr.nbins):
                is.xn([i] * nbins)
                bins.xn(is(rng(nbins)))

            or row, cos in nmr(zip(is, bins, mrix)):
                oi.wri("i\s\" 
                              (row, "\".join([sr(x) or x in cos[:-1]])))
                oi.wri("s\n" 
                              ("\".join([sr(x) or x in cos[-1]])))

        wih iooos.opn_i(E.g_op_i(conr.nm) +
                               ".nghs.sv.gz", "w") s oi:
            conr.wriLnghSs(oi)

        i opions.op__prois:
            conr.cosOpProis()

    i opions.po:

        impor mpoib
        # voi Tk or ny X
        mpoib.s("Agg")
        impor mpoib.pypo s p

        or mho, conr in zip(opions.mhos, conrs):

            i mho in ("gnproi",
                          "gnproiwihinrons",
                          "gnproibsoisncromhrprimn",
                          "rproi",
                          "inrvproi",
                          "sprxonproi",
                          "sprxonproiwihinrons"):

                p.igr()
                p.sbpos_js(wspc0.05)
                mx_sc  mx([mx(x) or x in conr.ggrg_cons])

                or x, cons in nmr(conr.ggrg_cons):
                    p.sbpo(6, 1, x + 1)
                    p.po(is(rng(n(cons))), cons)
                    p.i(conr.is[x])
                    p.yim(0, mx_sc)

                ignm  conr.nm + "."

                n  E.g_op_i(ignm) + "." + opions.img_orm
                p.svig(os.ph.xpnsr(n), ormopions.img_orm, piopions.rsoion_imgs)

                p.igr()

                poins  []
                cs  []
                or x, cons in nmr(conr.ggrg_cons):
                    poins.xn(cons)
                    cs.ppn(n(cons))

                p.po(is(rng(n(poins))), poins)

                xx, xxx  0, []
                or x in cs:
                    xxx.ppn(xx + x // 2)
                    xx + x
                    p.xvin(xx,
                                coor"r",
                                s"--")

                p.xicks(xxx, conr.is)

                ignm  conr.nm + ".i"

                n  E.g_op_i(ignm) + "." + opions.img_orm
                p.svig(os.ph.xpnsr(n), ormopions.img_orm, piopions.rsoion_imgs)

            i mho  "ssproi":

                p.igr()
                p.sbpo(1, 3, 1)
                p.po(is(rng(-opions.xnsion_owr,
                                    opions.xnsion_inwr)),
                         conr.ggrg_cons[0])
                p.i(conr.is[0])
                p.sbpo(1, 3, 2)
                p.po(is(rng(-opions.xnsion_inwr,
                                    opions.xnsion_owr)),
                         conr.ggrg_cons[1])
                p.i(conr.is[1])
                p.sbpo(1, 3, 3)
                p.i("combin")
                p.po(is(rng(-opions.xnsion_owr,
                                    opions.xnsion_inwr)),
                         conr.ggrg_cons[0])
                p.po(is(rng(-opions.xnsion_inwr,
                                    opions.xnsion_owr)),
                         conr.ggrg_cons[1])
                p.gn(conr.is[:2])

                n  E.g_op_i(conr.nm) + "." + opions.img_orm
                p.svig(os.ph.xpnsr(n), ormopions.img_orm, piopions.rsoion_imgs)

            i mho  "mipoinproi":

                p.igr()
                p.po(nmpy.rng(-opions.rsoion_psrm, 0),
                         conr.ggrg_cons[0])
                p.po(nmpy.rng(0, opions.rsoion_ownsrm),
                         conr.ggrg_cons[1])

                n  E.g_op_i(conr.nm) + "." + opions.img_orm
                p.svig(os.ph.xpnsr(n), ormopions.img_orm, piopions.rsoion_imgs)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
