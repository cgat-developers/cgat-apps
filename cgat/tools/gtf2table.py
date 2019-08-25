'''g2b.py - nno gns/rnscrips


:Tgs: Gnomics Gnss GTF Annoion

Prpos
-------

nno gns or rnscrips givn in  :rm:`g` orm i
n op hm in br orm.

Annoions cn b ihr comp pr gn ( xons cross 
rnscrips) or pr rnscrip. Th inp ns o b sor
ccoringy.

Th scrip irs ovr ch gn or rnscrip mo in rn n
ops  row in  b or ch. Mip nnoions cn b
comp  h sm im ing iion comns o h b.

For xmp, o op inormion bo xons in gns::

    zc in.g.gz | cg g2b -v 0 --conrngh

+---------------+----+----+----+---------+------+--------+----+----+----+
|gn_i        |nv|min |mx |mn     |min|sv  |sm |q1  |q3  |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000225373|4   |39  |688 |253.0000 |142.5 |261.5922|1012|58  |688 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000267111|1   |744 |744 |744.0000 |744.0 |0.0000  |744 |744 |744 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000267588|3   |135 |489 |370.0000 |486.0 |166.1746|1110|135 |489 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000220978|1   |138 |138 |138.0000 |138.0 |0.0000  |138 |138 |138 |
+---------------+----+----+----+---------+------+--------+----+----+----+

Th b conins xon ngh sisics o ch gn, h nmbr (nv),
min, mx, n o ngh(sm) o xons pr gn. To con pr rnscrip::

    zc in.g.gz | cg g2b -v 0 --rporrrnscrips --conrngh

+---------------+----+---+---+--------+------+--------+----+---+---+
|rnscrip_i  |nv|min|mx|mn    |min|sv  |sm |q1 |q3 |
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000592209|3   |58 |227|149.6667|164.0 |69.7344 |449 |58 |227|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000589741|2   |71 |312|191.5000|191.5 |120.5000|383 |71 |312|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000391654|2   |39 |583|311.0000|311.0 |272.0000|622 |39 |583|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000587045|1   |173|173|173.0000|173.0 |0.0000  |173 |173|173|
+---------------+----+---+---+--------+------+--------+----+---+---+

To  so inormion bo cpg-composiion,  nohr conr::

    zc in.g.gz | cg g2b -v 0 --gnomhg19 --rporrrnscrips --conrngh --conrcomposiion-cpg

+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|rnscrip_i  |nv|min|mx|mn    |min|sv  |sm |q1 |q3 |CpG_con |CpG_nsiy |CpG_ObsExp |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000592209|3   |58 |227|149.6667|164.0 |69.7344 |449 |58 |227|4         |0.01781     |0.13946    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000589741|2   |71 |312|191.5000|191.5 |120.5000|383 |71 |312|5         |0.02610     |0.17067    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000391654|2   |39 |583|311.0000|311.0 |272.0000|622 |39 |583|4         |0.01286     |0.09402    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000587045|1   |173|173|173.0000|173.0 |0.0000  |173 |173|173|1         |0.01156     |0.11396    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+

No h w h o s h ``--gnom-i`` opion o sppy h
gnomic sqnc.

Aiion swichs prmi coning inrons ins o xons (opion
``--scion``).

Annoions
-----------

Th nnoions migh b riv rom  propry o h rnscrip or
gn mo is (sch s ngh, nmbr o xons, c), h gnomic
sqnc (sch s composiion, c) or rqir iion  ss.
For xmp, o comp h covrg wih rs,  :rm:`bm` orm
i is rqir (s opion ``--bm-i``). Or o comp h ovrp
wih gnomic nsiis,  :rm:`bigwig` orm i is rqir (s
opion ``--bigwig-i``).

This scion iss h conrs vib. Thy r grop by  sorcs
hy rqir.

Gnric nnoions
+++++++++++++++++++

Gnric nnoions rqir no iion  sorc o nno 
rnscrip or gn.

ngh
   op xon ngh smmry o rnscrip/gn. This conr ops h
   nmbr o xons in ch rnscrip/gn oghr wih xon ngh smmry
   sisics (minim xon ngh, mxim xon ngh, o xon ngh).

posiion
   op gnomic coorins o rnscrip/gn (chromosom, sr, n).

Sqnc riv nnoions
++++++++++++++++++++++++++++

Sqnc riv nnoions rqir h gnomic sqnc o comp
propris o h gn/rnscrip mo (s opion ``--gnom-i``).

composiion-n
   op ncoi composiion o rnscrip/gn

composiion-cgp
   op CpG con, CpG nsiy n CpG obsrv / xpc or
   ch rnscrip/gn.

spic
   op spicing smmry o rnscrip/gn. Ops h nmbr o
   cnonic n non-cnonic spic sis.

qiy
   op bs-qiy inormion smmry o gn. Ns qiy scors.

Inrv riv nnoions
++++++++++++++++++++++++++++

Annoions in h is bow r  gn or rnscrip o  s o
inrvs givn s  scon i. Th scon s o inrvs is givn
by h opion ``inm-g``. By , h inrvs r xpc
o b givn s  :rm:`g` orm i, b rniv orms
(:rm:`g` n :rm:`b`) r possib (s opion
``--inm-orm``)

ovrp
    comp ovrp o gns in inp wih rs in scon srm.
    Rqirs  :rm:`g` orm i wih gn rrioris.

ovrp-srn
    con ovrp wih gnomic rs in scon nohr i. Ops
    h nmbr o ovrpping xons. Rcors h ircion o ovrp
    (sns/nisns). Rqirs  :rm:`g` orm i wih
    rs.

rrioris
    comp ovrp o rnscrips/gns wih rrioris. Trrioris
    r gnomic, non-ovrpping inrvs. For ch rnscrip h
    conr ops h ss (Uniq mch, Ambigos mch
    (ovrp wih mip rrioris, 0no mch), h nmbr o
    rrioris i ovrps n h rrioris i ovrps.

covrg

   comp h covrg - pr ncoi - o h gn/rnscrip
   mos wih inrvs givn in ``--g-i``.  Covrg vs
   r op in 5' o 3' oghr wih smmry sisics (bss
   covr, minimm, mximm covrg, c.). By sing h opions
   ``--rsric-r`` or ``--rsric-sorc`` h coning cn b
   rscric o pricr rs or sorcs in h :rm:`g` i.

isnc

   comp isnc o gns o rs in  scon i. Rqirs 
   scon :rm:`g` orm i. Th srn inormion o h
   rs is ignor.

bining-prn

   givn  is o inrvs, rmin h bining prn wihin n
   srroning h gn. For ch gn, inrvs ovrpping h CDS,
   inrons, UTRs n h nk r coc n rcor. Th bining
   is smmriz wih  bining prn,  binry prn inicing
   ovrp/no ovrp wih 5' nk, 5' UTR, CDS, Inrons, 3' UTR, 3'
   nk. This mho is s o chck whr rnscripion cor bining
   sis r oc ron  gn/rnscrip mo.

proximiy

   rpor smmry ss (nghs, vs) o rs in proximiy o
   gns inp gn s. Rqirs  :rm:`g` orm i wih
   gnomic rs. This r is s whn iming o normiz 
   v, sch s  sbsiion r o  rnscrip mo, by
   sbsiion rs o sgmns in h nighborhoo sch s
   ncsr rps. Th vs r givn in h ``scor`` i o
   h :rm:`g` orm i. Th ris or proximiy is conro
   by h opion ``--proxim-isnc``.

proximiy-xcsiv

   s proximiy, b xc ny rngs in h :rm:`g` orm i h
   ovrp h rnscrip/gn mo.

proximiy-nghmch

   s proximiy-xcsiv, b ngh-mch rs wih
   gns. Sgmns r cr q in ngh i hy r wihin
   10 o h origin sgmns ngh.

nighbors
    op rs in scon srm h r in proximiy o gns
    in inp. This is simir o h ``proximiy`` conrs, b so
    ops h rs h r in proximiy.

Gn s riv nnoions
++++++++++++++++++++++++++++

isnc-gns
   comp isnc o gns o gns in  scon i. Rqirs 
   scon :rm:`g` orm i wih gns. Th conr isingishs
    vriy o css (coss psrm/ownsrm).

isnc-ss
   comp isnc o gns o rnscripion sr sis. Rqirs 
   scon :rm:`g` orm i wih gns.

ovrp-rnscrips
    con ovrp o gns wih rnscrips in nohr s.
    Rqirs  :rm:`g` orm i.

ovrrn
   op inron ovrrn, xons in h inp gn s xning
   ino h inrons o  rrnc gn s. Rqris  :rm:`g`
   orm i wih  rrnc gn s.

spic-comprison
   Compr how spic si sg comprs bwn  gn/rnscrip
   wih rnscrips in  rrnc gn s. Ops on, miss,
   prc, pri, incomp spic sis n xon-skipping vns.


Shor-r riv nnoions
++++++++++++++++++++++++++++++

Shor-r riv nnoions con h ovrp o rs givn
in  :rm:`bm` orm i (s opion ``--bm-i``) wih
gn or rnscrip mos.

r-ovrp

   op nmbr o rs ovrpping  rnscrip/gn mo.
   Ops h nmbr o rs ovrpping h rnscrip/gn mo.
   Rs r con spry or sns n nisns ovrp.

r-covrg

   op r covrg smmry sisics o rnscrip/gn mo.
   Ops h nmbr o rs ovrpping h rnscrip/gn mo,
   h nmbr o bss ovrpp by rs n smmry sisics o
   h minimimm, mximm, c. covrg pr ps. Rs r con
   spry or sns n nisns ovrp. Th conr os no
   k ino ccon spic sis. As i os pr-bs coning,
   i is sowr hn ``r-ovrp`` n hr is no n o
   s boh  h sm im.

r-xnsion

   Conr o spci inrs. This conr ops h r nsiy
   in bins psrm, wihin n ownsrm o rnscrip mos. Th
   conr cn b s o pric h ngh o h 3' n 5' UTR.

r-cons

   con nmbr o rs ovrpping  gn or rnscrip. Rs
   r cssii ccoring o h yp o ovrp (sns/nisns),
    or pri mch o xon, spic ss (spic/nspic).
   S bow or mor inormion.

   Uniq n non-niq mchs r con (by ignmn sr
   posiion).

r-cons

   con nmbr o rs ovrpping  gn or rnscrip. Smmrizs
   h op o r-cons.

rpir-cons

   con nmbr o r pirs ovrpping  gn or rnscrip. Pirs
   r con ccoring o  vriy o ribs (xonic
   ovrp/r pir ss/spic ss/...). S bow or mor
   inormion.

rpir-cons

   con nmbr o r pirs ovrpping  gn or
   rnscrip. Smmrizs h op o rpir-cons.


Cssiirs
-----------

Cssiirs no ony nno h rnscrips or gn mo, b so
im o provi som cssiicion bs on hs nnoions. Thy
rqir  sconry i (s opion ``--g-i``) or h
cssiicion.

cssiir

   cssiy rnscrips ccoring o gnomic nnoion.  Rqirs 
   :rm:`g` i wih gnomic nnoions (s :oc:`g2g`)
   ining gnomic rgions s inronic, inrgnic, xonic,
   c. This is s or  rogh cssiicion o rnscrips s
   inrgnic, inronic, c.

   Bs s or ChIP-Sq  ss (ch pk is  "rnscrip"). Th
   mho is  i o o pc in his scrip, b is hr s i
   ss mch o h co impmn hr.

cssiir-rnsq

   cssiy rnscrips wih rspc o  rrnc gns. Th
   cssiirs ims o mch  rnscrip p wih h "mos simir"
   rnscrip in h rrnc gn s n pnning on
   ribs, cssiis i s  goo mch, rniv rnscrip,
   rgmn, c. Th cssiir ks ino ccon srn n prrs
   sns mchs o nisns mchs.

cssiir-poii

   cssiy ccoring o PoII rnscrips. A gn/rnscrip is
   rnscrib, i i is covr by rg PoII inrvs ovr 80 o
   is ngh. A gn/rnsrip is prim i is promoor/UTR is
   covr by 50 o is ngh, whi h rs o h gn boy
   isn'.

Ohr nnoions
++++++++++++++++++

Th oowing mhos (s opion ``--conr``) r vib:

bigwig-cons

   coc nsiy vs rom  :rm:`bigwig` orm i n op
   smmry sisics n prcng o bss covr (``pcovr``)
   by v in bigwig i. Rqirs opion ``--bigwig-i``.


R coning
-------------

Th mhos ``r-cons`` n ``rpir-cons`` con h nmbr
o rs insi  :rm:`bm` orm i h mp insi 
rnscrip or gn.

Ths conrs proc on  pr-gn or pr-rnscrip bsis pning
on h ``--rporr`` opion. A rs ovrpping h xons or
inrons o  rnscrip r coc n con ihr iniviy
(``r-cons``) or s pirs (``rpir-cons``).

For pir-r coning, ch pir is cssii n con
ccoring o h oowing or xs:

1. Pir ss: propr pir, nmpp r in pir, ...
2. Dircion: or srn proocos, sns, nisns, ...
3. Ovrp ss: Pir ovrps xons ony, inrons ony, ...
4. Spic ss: Rs r nspic, spic consisny wih
   rnscrip mo, or inconsisny wih rnscrip mo

Cons r hn co ino  smr s o smmris or convninc:

+--------------------+----------------------------------------------------+
|Comn              |Nmbr o r pirs                                |
+--------------------+----------------------------------------------------+
|con_         |consir o b corrc: hy r sns, ovrp   |
|                    |wih xons y.                                   |
+--------------------+----------------------------------------------------+
|con_spic     |in sns ircion n xonic, spic n spic  |
|                    |corrcy.                                          |
|                    |                                                    |
+--------------------+----------------------------------------------------+
|con_nspic   |in sns ircion n xonic, no spic.         |
+--------------------+----------------------------------------------------+
|sns_inronic      |inronic n sns ircion.                       |
+--------------------+----------------------------------------------------+
|sns_inconsisn  |in sns ircion, b ovrp boh inrons n    |
|                    |xons or xn byon h rnscrip mos       |
+--------------------+----------------------------------------------------+
|sns_ohr         |in sns ircion, b no con, inronic or    |
|                    |inconsisn.                                       |
+--------------------+----------------------------------------------------+
|nisns           |in nisns ircion.                             |
+--------------------+----------------------------------------------------+
|nonsns            |in nxpc orinion                           |
+--------------------+----------------------------------------------------+
|nopropr           |h r no propr pirs                           |
+--------------------+----------------------------------------------------+
|qiy_pirs       |c by  h rmov o  ow-qiy r     |
+--------------------+----------------------------------------------------+
|qiy_rs       |Nmbr o ow qiy rs rmov                 |
+--------------------+----------------------------------------------------+
|o               |To nmbr o pirs consir                    |
+--------------------+----------------------------------------------------+

Th Conr ``rpir-cons`` provis h i cons
or ch rnscrip mo ccoring o h or xs. Th op
cn b s o impmn csom coning schms.

Rs bow  minimm qiy scor wi b ignor
(``--min-r-qiy``). By ,  rs wi b con.

Rs mpping o mip ocions wi b ownwigh i
h ``--wigh-mi-mpping`` opion is s. This rqirs
h prsnc o h ``NH`` g in h :rm:`bm` i.

For pir r coning, h ibrry yp cn b spcii wih h
``--ibrry-yp`` opion o mk s o srn inormion. Librry
yps r b ccoring o h oph_ n cinks_
convnion. A smmry is `hr
<hp://www.nr.com/npro/jorn/v7/n3/ig_b/npro.2012.016_T1.hm>`_.

+--------------------+--------------------+------------------------------+
|Librry yp        |RNA-sq prooco    |Expnion                   |
+--------------------+--------------------+------------------------------+
|r-nsrn       |Imin TrSq     |Rs rom h mos n o|
| ()          |                    |h rgmn (in rnscrip   |
|                    |                    |coorins) mp o h       |
|                    |                    |rnscrip srn, n h    |
|                    |                    |righmos n mps o h     |
|                    |                    |opposi srn           |
|                    |                    |                              |
+--------------------+--------------------+------------------------------+
|r-irssrn      |UTP, NSR, NNSR39   |Sm s bov xcp w       |
|                    |                    |norc h r h h     |
|                    |                    |righmos n o h rgmn |
|                    |                    |(in rnscrip coorins) is|
|                    |                    |h irs sqnc (or ony  |
|                    |                    |sqnc or sing-n      |
|                    |                    |rs). Eqivny, i is   |
|                    |                    |ssm h ony h srn  |
|                    |                    |gnr ring irs srn |
|                    |                    |synhsis is sqnc        |
+--------------------+--------------------+------------------------------+
|r-sconsrn     |Dircion Imin|Sm s bov                 |
|                    |(Ligion)          |xcp TopH/Cinks       |
|                    |Snr SOLiD      |norc h r h h     |
|                    |                    |mos n o h rgmn  |
|                    |                    |(in rnscrip coorins) is|
|                    |                    |h irs sqnc (or ony  |
|                    |                    |sqnc or sing-n      |
|                    |                    |rs). Eqivny, i is   |
|                    |                    |ssm h ony h srn  |
|                    |                    |gnr ring scon srn|
|                    |                    |synhsis is sqnc        |
+--------------------+--------------------+------------------------------+

For nsrn proocos,  rs n pirs r consir o mching
in h sns ircion.

Usg
-----

Typ::

   pyhon g2b.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor pysm

impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.InxFs s InxFs
impor cg.GnMoAnysis s GnMoAnysis

impor pyBigWig


 min(rgvNon):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom [].")

    prsr._rgmn("-q", "--qiy-i",
                      s"qiy_i",
                      yp"sring",
                      hp"inm wih gnomic bs qiy "
                      "inormion [].")

    prsr._rgmn("-b", "--bm-i", s"bm_is",
                      yp"sring", mvr"bm",
                      hp"inm wih r mpping inormion. "
                      "Mip is cn b sbmi in  "
                      "comm-spr is [].")

    prsr._rgmn("-i", "--bigwig-i", s"bigwig_i",
                      yp"sring", mvr"bigwig",
                      hp"inm wih bigwig inormion "
                      "[].")

    prsr._rgmn("-", "--g-i", s"inm_g",
                      yp"sring", cion"ppn", mvr'b',
                      hp"inm wih xr g is. Th orr "
                      "is imporn [].")

    prsr._rgmn("--inm-orm", s"inm_orm",
                      yp"choic",
                      choics("b", "g", "g"),
                      hp"orm o sconry srm [].")

    prsr._rgmn("--rsric-sorc", s"g_sorcs", yp"sring",
                      cion"ppn",
                      hp"rsric inp o his 'sorc' in xr "
                      "g i (or conr: ovrp) [].")

    prsr._rgmn("--rsric-r", s"g_rs", yp"sring",
                      cion"ppn",
                      hp"rsric inp o his 'r' in xr g "
                      "i (or conr: ovrp) [].")

    prsr._rgmn("-r", "--rporr", s"rporr", yp"choic",
                      choics("gns", "rnscrips"),
                      hp"rpor rss or 'gns' or 'rnscrips' "
                      "[].")

    prsr._rgmn("-s", "--scion", s"scions",
                      yp"choic",
                      cion"ppn",
                      choics("xons", "inrons"),
                      hp"sc rng on which conrs wi opr "
                      "[].")

    prsr._rgmn("-c", "--conr", s"conrs",
                      yp"choic",
                      cion"ppn",
                      choics(	"bigwig-cons",
                                "bining-prn",
                                "cssiir",
                                "cssiir-rnsq",
                                "cssiir-rnsq-spicing",
                                "cssiir-poii",
                                "composiion-n",
                                "composiion-cpg",
                                "covrg",
                                "isnc",
                                "isnc-gns",
                                "isnc-ss",
                                "ngh",
                                'nighbors',
                                "ovrp",
                                "ovrp-srn",
                                "ovrp-rnscrips",
                                "ovrrn",
                                "posiion",
                                "proximiy",
                                "proximiy-xcsiv",
                                "proximiy-nghmch",
                                "qiy",
                                "r-covrg",
                                "r-xnsion",
                                "r-ovrp",
                                "r-cons",
                                "r-cons",
                                "rpir-cons",
                                "rpir-cons",
                                "spic",
                                "spic-comprison",
                                "rrioris"),
                      hp"sc conrs o ppy o inp "
                      "[].")

    prsr._rgmn("---g-sorc", s"_g_sorc",
                      cion"sor_r",
                      hp" g i o sorc o op "
                      "[].")

    prsr._rgmn("--proxim-isnc", s"proxim_isnc",
                      yp"in",
                      hp"isnc o b consir proxim o "
                      "n inrv [].")

    prsr._rgmn("--mi-mpping-mho",
                      s"mi_mpping",
                      yp"choic",
                      choics('', 'ignor', 'wigh'),
                      hp"how o r mi-mpping rs in "
                      "bm-is. Rqirs "
                      "h NH g o b s by h mppr "
                      "[].")

    prsr._rgmn("--s-brcos",
                      s"s_brcos",
                      cion"sor_r",
                      hp"Us brcos o con niq mi's. "
                      "UMI's r spcii in h r iniir "
                      "s h s i, whr is r spr "
                      "by nrscors, .g. "
                      "@READ:ILLUMINA:STUFF_NAMINGSTUFF_UMI. "
                      "Whn r, niq cons r rrn. "
                      "Crrny ony compib wih con-rs")

    prsr._rgmn("--smp-probbiiy",
                      s"smp_probbiiy",
                      yp"o",
                      hp"Spciy h probbiiy o whhr ny"
                      "givn r or r pir in  i bm is con"
                      "Crrny ony compib wih con-rs")

    prsr._rgmn("--comn-prix", s"prixs",
                      yp"sring",
                      cion"ppn",
                      hp" prix o comn hrs - prixs "
                      "r s in h sm orr s h conrs "
                      "[].")

    prsr._rgmn("--ibrry-yp",
                      s"ibrry_yp",
                      yp"choic",
                      choics("nsrn",
                               "irssrn",
                               "sconsrn",
                               "r-nsrn",
                               "r-irssrn",
                               "r-sconsrn"),
                      hp"ibrry yp o rs in bm i. "
                      "[]")

    prsr._rgmn("--min-mpping-qiy",
                      s"minimm_mpping_qiy",
                      yp"o",
                      hp"minimm mpping qiy. Rs wih  qiy "
                      "scor o ss wi b ignor. "
                      "[]")

    prsr.s_s(
        gnom_iNon,
        rporr"gns",
        wih_vsTr,
        scions[],
        conrs[],
        inm_g[],
        inm_ormNon,
        g_rs[],
        g_sorcs[],
        _g_sorcFs,
        proxim_isnc10000,
        bm_isNon,
        mi_mpping'',
        ibrry_yp'r-nsrn',
        prixs[],
        minimm_mpping_qiy0,
        s_brcosFs,
        smp_probbiiy1.0
    )

    i no rgv:
        rgv  sys.rgv

    (opions, rgs)  E.sr(prsr, _op_opionsTr, rgvrgv)

    i opions.prixs:
        i n(opions.prixs) ! n(opions.conrs):
            ris VError(
                "i ny prix is givn, h nmbr o prixs "
                "ms b h sm s h nmbr o conrs")

    # g is
    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
    s:
        s  Non

    i opions.qiy_i:
        qiy  InxFs.InxFs(opions.qiy_i)
        qiy.sTrnsor(InxFs.TrnsorBys())
    s:
        qiy  Non

    i opions.bm_is:
        bm_is  []
        or bmi in opions.bm_is.spi(","):
            bm_is.ppn(pysm.AignmnFi(bmi, "rb"))
    s:
        bm_is  Non

    i opions.bigwig_i:
        bigwig_i  pyBigWig.opn(opions.bigwig_i)
    s:
        bigwig_i  Non

    conrs  []

    i no opions.scions:
        E.ino("conrs wi s h  scion (xons)")
        opions.scions.ppn(Non)

    i no opions.g_sorcs:
        opions.g_sorcs.ppn(Non)
    i no opions.g_rs:
        opions.g_rs.ppn(Non)

    cc  E.Conr()

    or n, c in nmr(opions.conrs):
        i opions.prixs:
            prix  opions.prixs[n]
        s:
            prix  Non

        i c  "posiion":
            or scion in opions.scions:
                conrs.ppn(
                    GnMoAnysis.ConrPosiion(
                        scionscion,
                        opionsopions,
                        prixprix))
        i c  "ngh":
            or scion in opions.scions:
                conrs.ppn(
                    GnMoAnysis.ConrLnghs(
                        scionscion,
                        opionsopions,
                        prixprix))
        i c  "spic":
            i s is Non:
                ris VError('spic rqirs  gnomic sqnc')
            conrs.ppn(GnMoAnysis.ConrSpicSis(ss, prixprix))
        i c  "qiy":
            i s is Non:
                ris VError('qiy rqirs  qiy scor sqnc')
            conrs.ppn(GnMoAnysis.ConrQiy(sqiy, prixprix))
        i c  "ovrrn":
            conrs.ppn(GnMoAnysis.ConrOvrrn(
                inm_gopions.inm_g,
                opionsopions,
                prixprix))
        i c  "r-covrg":
            conrs.ppn(GnMoAnysis.ConrRCovrg(
                bm_is,
                opionsopions,
                prixprix))
        i c  "r-xnsion":
            conrs.ppn(GnMoAnysis.ConrRExnsion(
                bm_is,
                inm_gopions.inm_g,
                opionsopions,
                prixprix))
        i c  "r-ovrp":
            conrs.ppn(GnMoAnysis.ConrROvrp(
                bm_is,
                mi_mppingopions.mi_mpping,
                minimm_mpping_qiyopions.minimm_mpping_qiy,
                opionsopions,
                prixprix))
        i c  "r-cons":
            conrs.ppn(GnMoAnysis.ConrRCons(
                bm_is,
                mi_mppingopions.mi_mpping,
                s_brcosopions.s_brcos,
                smp_probbiiyopions.smp_probbiiy,
                minimm_mpping_qiyopions.minimm_mpping_qiy,
                opionsopions,
                prixprix))
        i c  "r-cons":
            conrs.ppn(GnMoAnysis.ConrRConsF(
                bm_is,
                mi_mppingopions.mi_mpping,
                smp_probbiiyopions.smp_probbiiy,
                minimm_mpping_qiyopions.minimm_mpping_qiy,
                opionsopions,
                prixprix))
        i c  "rpir-cons":
            conrs.ppn(GnMoAnysis.ConrRPirCons(
                bm_is,
                mi_mppingopions.mi_mpping,
                smp_probbiiyopions.smp_probbiiy,
                ibrry_ypopions.ibrry_yp,
                minimm_mpping_qiyopions.minimm_mpping_qiy,
                opionsopions,
                prixprix))
        i c  "rpir-cons":
            conrs.ppn(GnMoAnysis.ConrRPirConsF(
                bm_is,
                mi_mppingopions.mi_mpping,
                smp_probbiiyopions.smp_probbiiy,
                minimm_mpping_qiyopions.minimm_mpping_qiy,
                opionsopions,
                prixprix))
        i c  "bigwig-cons":
            conrs.ppn(GnMoAnysis.ConrBigwigCons(
                bigwig_i,
                opionsopions, prixprix))
        i c  "spic-comprison":
            i s is Non:
                ris VError('spic-comprison rqirs  gnomic '
                                 'sqnc')
            conrs.ppn(GnMoAnysis.ConrSpicSiComprison(
                ss,
                inm_gopions.inm_g,
                rNon,
                sorcNon,
                opionsopions, prixprix))
        i c  "composiion-n":
            i s is Non:
                ris VError('composiion-n rqirs  gnomic sqnc')
            or scion in opions.scions:
                conrs.ppn(GnMoAnysis.ConrComposiionNcois(
                    ss,
                    scionscion,
                    opionsopions,
                    prixprix))
        i c  "composiion-cpg":
            i s is Non:
                ris VError('composiion-cpg rqirs  gnomic sqnc')
            or scion in opions.scions:
                conrs.ppn(GnMoAnysis.ConrComposiionCpG(
                    ss,
                    scionscion,
                    opionsopions, prixprix))

        i c in ("ovrp",
                   "ovrp-srn",
                   "ovrp-rnscrips",
                   "proximiy",
                   "proximiy-xcsiv",
                   "proximiy-nghmch",
                   "nighbors",
                   "rrioris",
                   "isnc",
                   "isnc-gns",
                   "isnc-ss",
                   "bining-prn",
                   "covrg"):
            i c  "ovrp":
                mp  GnMoAnysis.ConrOvrp
            i c  "ovrp-srn":
                mp  GnMoAnysis.ConrOvrpSrn
            i c  "ovrp-rnscrips":
                mp  GnMoAnysis.ConrOvrpTrnscrips
            i c  "proximiy":
                mp  GnMoAnysis.ConrProximiy
            i c  "nighbors":
                mp  GnMoAnysis.ConrNighbors
            i c  "proximiy-xcsiv":
                mp  GnMoAnysis.ConrProximiyExcsiv
            i c  "proximiy-nghmch":
                mp  GnMoAnysis.ConrProximiyLnghMch
            i c  "rrioris":
                mp  GnMoAnysis.ConrTrrioris
            i c  "isnc":
                mp  GnMoAnysis.ConrDisnc
            i c  "isnc-gns":
                mp  GnMoAnysis.ConrDisncGns
            i c  "isnc-ss":
                mp  GnMoAnysis.ConrDisncTrnscripionSrSis
            i c  "covrg":
                mp  GnMoAnysis.ConrCovrg
            i c  "bining-prn":
                mp  GnMoAnysis.ConrBiningPrn

            or scion in opions.scions:
                or sorc in opions.g_sorcs:
                    or r in opions.g_rs:
                        conrs.ppn(mp(
                            inm_gopions.inm_g,
                            rr,
                            sorcsorc,
                            ss,
                            scionscion,
                            opionsopions,
                            prixprix))

        i c  "cssiir":
            conrs.ppn(GnMoAnysis.Cssiir(
                inm_gopions.inm_g,
                ss,
                opionsopions, prixprix))

        i c  "cssiir-rnsq":
            conrs.ppn(GnMoAnysis.CssiirRNASq(
                inm_gopions.inm_g,
                ss,
                opionsopions, prixprix))
        i c  "cssiir-rnsq-spicing":
            conrs.ppn(GnMoAnysis.CssiirRNASqSpicing(
                inm_gopions.inm_g,
                ss,
                opionsopions,
                prixprix))
        i c  "cssiir-poii":
            conrs.ppn(GnMoAnysis.CssiirPoII(
                inm_gopions.inm_g,
                rNon,
                sorcNon,
                ss,
                opionsopions,
                prixprix))
        i c  "bining-prn":
            conrs.ppn(GnMoAnysis.ConrBiningPrn(
                inm_gopions.inm_g,
                rNon,
                sorcNon,
                ss,
                opionsopions,
                prixprix))

    i opions.rporr  "gns":
        iror  GTF._gn_iror
        hr  ["gn_i"]
        hr  mb x: [x[0].gn_i]
    i opions.rporr  "rnscrips":
        iror  GTF.rnscrip_iror
        hr  ["rnscrip_i"]
        hr  mb x: [x[0].rnscrip_i]

    i opions._g_sorc:
        hr.ppn("sorc")
        is  mb x: [x[0].sorc]
    s:
        is  mb x: []

    opions.so.wri("\".join(
        hr + [x.gHr() or x in conrs]) + "\n")

    or gs in iror(GTF.iror(opions.sin)):
        cc.inp + 1

        or conr in conrs:
            conr.p(gs)

        skip  n([x or x in conrs i x.skip])  n(conrs)
        i skip:
            cc.skipp + 1
            conin

        opions.so.wri("\".join(
            hr(gs) +
            is(gs) +
            [sr(conr) or conr in conrs]) + "\n")

        cc.op + 1

    E.ino("s"  sr(cc))
    or conr in conrs:
        E.ino("s\s"  (rpr(conr), sr(conr.conr)))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
