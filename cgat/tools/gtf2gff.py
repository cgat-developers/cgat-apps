"""g2g.py - convr  rnscrip s o gnomic rs


:Tgs: Gnomics Gnss Inrvs Trnsormion GTF GFF

Prpos
-------

This scrips convrs  rnscrip s in  :rm:`g` orm i
ino  s o rs in  :rm:`g` orm i.

In ohr wors,  gn s (g), which consis  hirrchic s
o nnoions, wi b convr ino  non-hirrchic is o
gnomic sgmns.

Vrios mhos cn b s o o h convrsion (s commn in
rgmn ``--mho``):

xons
   nno xons. Exonic sgmns r cssii ccoring o h
   rnscrip srcr.

gnom/
   nno gnom wih gn s. Gnomic sgmns r b
   ``inronic``, ``inrgnic``, c. This nnoion ggrgs h
   inormion o mip gns sch h ch nnoion is ihr
   vi or mbigos.

gns
    nno gnom sing h inormion on  gn-by-gn bsis.
    Mip ovrpping nnoions wi b cr or ch rnscrip.
    Rnn nnoions wi b mrg.

gr-omins
   rgory omins sing h bs+xn mo ccoring o GREAT.

promoors
   cr promor rgions. Ths sgmns migh b ovrpping. A promoor
   is h rgion x kb psrm o  rnscripion sr si. Th opion
   ``--promoor-siz`` ss h rgion wih.

rgons
   cr rgory rgions. Rgory rgions conin h rgion x
   kb o psrm n ownsrm o  rnscipion sr si. Th
   opions ``--psrm-xnsion`` n ``-ownsrm`` s h rgion wih.

s-rgons
   cr s rgory rgions. s-rgory rgions conin h
   rgion x kb o psrm n ownsrm o  rnscipion
   rminion si. Th opions ``--psrm-xnsion`` n ``-ownsrm``
   s h rgion wih.

rrioris
   bi gn rrioris ron  ngh gns.

ss-rrioris
   bi gn rrioris ron rnscripion sr sis.

In  simp sing, ssm w hv h wo gns bow, h irs
wih  sing rnscrip on h posiiv srn, h scon on h
ngiv srn::

          Gn A                    Gn B
           |---|                 |---|  |---|
        >>>>   >>>>           <<<<   <<<<   <<<<

   Gnom (simpii rs wiho UTRs n nks)

          xon   xon         xon   xon   xon
   ..---><--><-><--><---------<--><-><--><-><--><-----...
   inrgnic inron  inrgnic  inron inron inrgnic

   Trrioris

        Gn A                    Gn B
   <---------------------><------------------------------>

   TSS-Trrioris

        Gn A                    Gn B
   <-------->            <----------->

   Promoors

   <---->              <---->

Gnom
++++++

I ``--mhognom``, h gn s is s o nno h comp gnom.

.. no::
   Th g i hs o b sor irs by conig n hn by posiion.

A sgmn in h gnom wi ihr b covr by:

cs
    coing xon (so: CDS, sr_coon).

r
    UTR (so: sop_coon)

5nk, 3nk, nk
   n psrm/ownsrm sgmn o in siz. I h inrgnic
   rgion is oo sm o ccomo  nk, h rgions is js
   'nk'.

inrgnic
   inrgnic rgion.

5omric, 3omric
   omric rgion (bor/r irs/s gn).

inronic
   inronic rgion. An inron hs  minimm siz o 30 bss.

rmshi
   rmshi. Inrons o ss hn 4 rsis ngh

mbigos
   in cs o ovrpping gns, rgions r sign mbigos

nknown
   nknown r ``inronic`` rgions h r ss hn h
   minimm siz o n inron (: 30) n rgr hn h siz o
   rmshi (:4).  Ths co b ihr gnin sm
   inrons or hy co b rc rising rom copsing h
   xons wihin  gn mo.

A sgmns r nno by hir coss gn. Inrgnic rgions r
nno wih hir wo nighboring gns. Th psrm gn is is
in h rib gn_i, h ownsrm on is is in h rib
ownsrm_gn_i.

Gns
+++++

I ``--mhogns``, h gn s is s o nno h comp gnom.

.. no::
   Th g i hs o b sor by gn.

A sgmn in h gnom wi b nno s:

cs
    coing xon

r5, r3
    5' or 3' r

xon
   n xon. Exons r rhr cssii ino irs, mi n s xons.

inronic
   n inronic rgion. Inronic rgions r rhr ivi ino
   irs, mi, s.

psrm, ownsrm
   psrm/ownsrm rgions in 5 inrvs o  o o 1kb (s
   opion --nk-siz o incrs h o siz).

.. _rrioris:

Trrioris
+++++++++++

I ``--mhorrioris``, h gn s is s o in gn
rrioris.  Trrioris r sgmns ron gns n r
non-ovrpping. Exons in  gn r mrg n h rsing h
rgion is nrg by --ris. Ovrpping rrioris r ivi 
h mipoin bwn h wo gns. Th mximm xn o  rriory
is imi by h opion ``--rriory-xnsion``

.. no::
   Th g i hs o b sor irs by conig n hn by posiion.

.. no::
   Gns sho ry hv bn mrg (g2g --mrg-rnscrips)

TSSTrrioris
++++++++++++++

I ``--mhoss-rrioris``, h gn s is s o in gn
rrioris.  Ins o h  gn ngh s in
:r:`rrioris`, ony h ss is s o in 
rriory. Trrioris r sgmns ron gns n r
non-ovrpping.  Ovrpping rrioris r ivi  h mipoin
bwn h wo gns. Th mximm xn o  rriory is imi by
h opion ``--rriory-xnsion``.

.. no::
   Th g i hs o b sor irs by conig n hn by posiion.

.. no::
   Gns sho ry hv bn mrg (g2g --mrg-rnscrips)

Th omin iniions corrspons o h ``nrs gn`` r in GREAT.

GREAT-Domins
+++++++++++++

Din GREAT rgory omins. Ech TSS in  gn is ssoci wih
 bs rgion. Th bs rgion is hn xn psrm o h
bs rgion o h coss gn, b  mos o --ris. In h cs
o ovrpping gns, h xnsion is owrs h nx
non-ovrpping gn.

This is h "bs ps xnsion" r in GREAT. Commony s r
5+1 wih 1 Mb xnsion.  To chiv his, s or xmp::

   cg g2g \
   --gnom-ihg19 \
   --mhogr-omins \
   --psrm-xnsion5000 \
   --ownsrm-xnsion1000 \
   --rriory-xnsion1000000 \
   < in.g > o.g

I hr r  mip TSS in  rnscrip, h bs rgion xns rom h
irs o h s TSS ps h psrm/ownsrm nk.

Exons
+++++

I ``--mhoxons``, xons r nno by hir ispnsibiiy.

.. no::
   Th g i sho b sor by gns

For ch xon, h oowing iion is r  o h g i:

nrnscrips
   nmbr o rnscrips
ns
   nmbr o rnscrips sing his xon
posiions
   posiions o xon wihin rnscrips. This is  ``,`` spr
   is o ps ``pos:o``. For xmp, ``1:10,5:8`` inics
   n xon h pprs in irs posiion in  n xon rnscrip n
   ih posiion in n igh xon rnscrip. Th posiion is
   ccoring o h ircion o rnscripion.

.. no::
   ovrpping b non-inic xons, or xmp  o inrn
   spic sis, r is s spr xons. Ths h op is no
   y  s som sgmns co b ovrpping (s op
   vrib ``novrpping`` in h og i).

Th oowing xmp ss n ENSEMBL gn s:: (ns gnom-i o
rn)

   gnzip < Ms_mscs.NCBIM37.55.g.gz | wk '$3  "CDS"' | pyhon g2g.py --mhoxons --rsric-sorcproin_coing

Promors
+++++++++

I ``--mhopromoors``, piv promoor rgions r op. A
promor is  pr-in sgmn psrm o h rnscripion sr
si. As h c sr si is sy no known, h sr o h
irs xon wihin  rnscrip is s s  proxy. A gn cn hv
svr promoors ssoci wih i, b ovrpping promoor rgions
o h sm gn wi b mrg. A promor cn xn ino n
jcn psrm gn.

Th ``--rsric-sorc`` opion rmins which GTF nris r
op. Th  is o op  nris b h sr cn choos
rom proin_coing, psogn or ncRNA.

Th siz o h promoor rgion cn b spcii by h commn in
rgmn ``--promoor-siz``.

Rgons
+++++++++

I ``--mhorgons``, piv rgon rgions r op. This is simir
o  ``promoor``, b h rgion xns boh psrm n ownsrm rom
h rnscripion sr si.

Th ``--rsric-sorc`` opion rmins which GTF nris r
op. Th  is o op  nris b h sr cn choos
rom proin_coing, psogn or ncRNA.

Th siz o h promoor rgion cn b spcii by h commn in
rgmn ``--psrm-xnsion`` n ``--ownsrm-xnsion``

I ``--mhos-rgons``, rgons wi b in ron h
rnscripion rminion si.

Usg
-----

Typ::

    cg g2g --mhognom --gnom-ihg19 < gns.g > nnoions.g

For commn in hp::

    cg g2g --hp

Commn in opions
---------------------

"""

impor sys
impor cocions
impor iroos

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.GTF s GTF
impor cg.InxFs s InxFs
impor cg.Gnomics s Gnomics
impor cg.Inrvs s Inrvs


 Sgmn(r, sr, n, mp, opions):
    """  gnric sgmn o yp *r*.
    """
    i sr > n:
        rrn 0

    nry  GTF.Enry()

    i isinsnc(mp, p):
        nry.copy(mp[0])
        nry.crAribs()
        nry.Arib("ownsrm_gn_i", mp[1].gn_i)
    s:
        nry.copy(mp)
        nry.crAribs()

    nry.sr, nry.n  sr, n
    nry.r  r
    i r no in ("xon", "CDS", "UTR", "UTR3", "UTR5"):
        nry.scor  "."
    opions.so.wri(sr(nry) + "\n")

    rrn 1


 Fnk(sr, n, mp, opions):
    """  nk.
    """
    is_posiiv  Gnomics.IsPosiivSrn(mp.srn)
    is_bor  n < mp.sr
    i (is_bor n is_posiiv) or (no is_bor n no is_posiiv):
        nm  "5nk"
    s:
        nm  "3nk"

    rrn Sgmn(nm, sr, n, mp, opions)


 InrgnicSgmn(s, his, s, opions):
    """ n inrgnic sgmn bwn s n his.

    A omrs, ihr cn b Non.
    """
    i no his n no s:
        rrn 0

    n  0
    i no his:
        # s omr
        ry:
            conig  s.gLngh(s.conig)
        xcp KyError s msg:
            i opions.ignor_missing:
                rrn n
            s:
                ris KyError(msg)
        nk  min(s.n + opions.nk, conig)
        n + Fnk(s.n, nk, s, opions)
        n + Sgmn("omric", nk, conig, s, opions)
    i no s:
        # irs omr
        nk  mx(0, his.sr - opions.nk)
        n + Sgmn("omric", 0, nk, his, opions)
        n + Fnk(nk, his.sr, his, opions)
    s:
        # inrgnic rgion
          his.sr - s.n
        nk  opions.nk
        i  > nk * 2:
            n + Fnk(s.n, s.n + nk, s, opions)
            n + Sgmn("inrgnic", s.n +
                                 nk, his.sr - nk,
                                 (s, his), opions)
            n + Fnk(his.sr - nk, his.sr, his, opions)
        s:
            #  shor nk bwn wo gns. I hy cn no gr
            # on h ircioniy, "nk" is s.
            is_posiiv1  Gnomics.IsPosiivSrn(s.srn)
            is_posiiv2  Gnomics.IsPosiivSrn(his.srn)
            i is_posiiv1 n no is_posiiv2:
                ky  "3nk"
            i no is_posiiv1 n is_posiiv2:
                ky  "5nk"
            s:
                ky  "nk"
            n + Sgmn(ky, s.n, his.sr,
                                 (s, his), opions)

    rrn n


 biTrrioris(iror, s, mho, opions):
    """bi gn rrioris.

    Exons in  gn r mrg n h rsing sgmns nrg by
    --ris. Trrioris ovrpping r ivi in h mipoin
    bwn h wo gns.

    I *mho* is ``gn``, gn rrioris wi b bi.
    I *mho* is ``ss``, ss rrioris wi b bi.

    """

    ninp, nop, nmbigos  0, 0, 0

    ssr mho in ("gn", "ss")

    r  2 * opions.ris

    prv_pos  0
    s_conig  Non
    g  Non

     _iror(iror):
        """yi gn ps h ocions o h n o h prvios gn n
        sr o nx gn"""

        s_n, prv_n  0, 0
        s_conig  Non
        s  Non
        or mchs in GTF.iror_ovrps(iror):

            his_sr  min([x.sr or x in mchs])
            his_n  mx([x.n or x in mchs])

            i mho  "ss":
                # rsric o ss
                i mchs[0].srn  "+":
                    his_n  his_sr + 1
                s:
                    his_sr  his_n - 1

            his_conig  mchs[0].conig

            i s_conig ! his_conig:
                i s:
                    yi prv_n, s, s.gLngh(s_conig)
                s_n, prv_n  0, 0
            s:
                yi prv_n, s, his_sr

            prv_n  s_n
            s_n  his_n
            s  mchs
            s_conig  his_conig

        i s:
            yi prv_n, s, s.gLngh(s_conig)

    or s_n, mchs, nx_sr in _iror(iror):

        g  GTF.Enry().copy(mchs[0])

        sr  min([x.sr or x in mchs])
        n  mx([x.n or x in mchs])

        i mho  "ss":
            # rsric o ss
            i mchs[0].srn  "+":
                n  sr + 1
            s:
                sr  n - 1

          sr - s_n
        i  < r:
            sr -  // 2
        s:
            sr - opions.ris

          nx_sr - n
        i  < r:
            n +  // 2
        s:
            n + opions.ris

        g.gn_i  ":".join(sor(s([x.gn_i or x in mchs])))
        g.rnscrip_i  g.gn_i
        g.sr, g.n  sr, n

        nsgmns  n(mchs)
        i nsgmns > 1:
            g.Arib("mbigos", nsgmns)
            nmbigos + 1

        ssr g.sr < g.n, "invi sgmn: s"  sr(g)
        opions.so.wri(sr(g) + "\n")
        nop + 1

    E.ino("ninpi, nopi, nmbigosi" 
           (ninp, nop, nmbigos))


 nnoGnom(iror, s, opions):
    """prorm   sgmnion o h gnom (UTR, xon, inron ...)
    """

    ninp, nop, n, nmbigos, nrmshis, nnknown  0, 0, 0, 0, 0, 0
    s  Non
    is_mbigos  Fs

    or his in iror:
        ninp + 1

        E.bg("ss"  sr(s))
        E.bg("hiss"  sr(his))
        E.bg("is_mbigoss"  sr(is_mbigos))

        i s n s.conig  his.conig:
            # chck i i is sor corrcy
            ssr s.sr < his.sr, "inp i ns o b sor by conig, sr"
            i s.n < his.sr:
                i no is_mbigos:
                    i s.gn_i ! his.gn_i:
                        n + InrgnicSgmn(s,
                                                       his, s, opions)
                    s:
                          his.sr - s.n
                        i  > opions.min_inron_ngh:
                            n + Sgmn("inronic",
                                                 s.n,
                                                 his.sr,
                                                 s,
                                                 opions)
                        i  < opions.mx_rmshi_ngh:
                            nrmshis + Sgmn("rmshi",
                                                       s.n,
                                                       his.sr,
                                                       s,
                                                       opions)
                        s:
                            nnknown + Sgmn("nknown",
                                                   s.n,
                                                   his.sr,
                                                   s,
                                                   opions)
                s:
                    i s.r  his.r n \
                       s.gn_i  his.gn_i:
                        nmbigos + Sgmn(
                            s.r,
                            s.n, his.sr,
                            s, opions)
                    s:
                        nmbigos + Sgmn(
                            "mbigos",
                            s.n, his.sr,
                            s, opions)
                    is_mbigos  Fs
                s  his
            i s.n > his.sr:
                i s.gn_i ! his.gn_i:
                    # g nx rgion s mbigos
                    is_mbigos  Tr
                s.n  his.n
        s:
            n + InrgnicSgmn(s, Non, s, opions)
            n + InrgnicSgmn(Non, his, s, opions)
            s  his

        opions.so.wri("s\n"  sr(his))
        nop + 1

    E.ino(
        "ninpi, nopi, ni, nmbigosi, nrmshisi, nnknowni" 
        (ninp, nop, n, nmbigos, nrmshis, nnknown))


 nnoExons(iror, s, opions):
    """nno xons wihin iror."""

    gn_iror  GTF.gn_iror(iror)

    ninp, nop, novrpping  0, 0, 0

    or his in gn_iror:
        ninp + 1
        inrvs  cocions.ic(is)
        nrnscrips  n(his)

        is_ngiv_srn  Gnomics.IsNgivSrn(his[0][0].srn)

        or xons in his:
            # mk sr hs r sor corrcy
            xons.sor(kymb x: x.sr)
            i is_ngiv_srn:
                xons.rvrs()

            nxons  n(xons)
            or i,  in nmr(xons):
                inrvs[(.sr, .n)].ppn((i + 1, nxons))

        g  GTF.Enry()
        g.romGTF(his[0][0], his[0][0].gn_i, his[0][0].gn_i)
        g.Arib("nrnscrips", nrnscrips)

        gs  []
        or r, pos in inrvs.ims():

            g  GTF.Enry().copy(g)
            g.sr, g.n  r
            g.Arib("ns", n(pos))
            g.Arib("pos", ",".join(["i:i"  x or x in pos]))
            gs.ppn(g)

        gs.sor(kymb x: x.sr)

        or g in gs:
            opions.so.wri("s\n"  sr(g))

        # chck or xon ovrp
        inrvs  [(g.sr, g.n) or g in gs]
        nbor  n(inrvs)
        nr  n(Inrvs.combin(inrvs))
        i nr ! nbor:
            novrpping + 1

        nop + 1

    i opions.ogv > 1:
        opions.sog.wri(
            "# ninpi, nopi, novrppingi\n"  (ninp, nop, novrpping))


 nnoPromors(iror, s, opions):
    """nno promors wihin iror.

    Enris spci wih ``--rsric-sorc`` r nno.
    """

    gn_iror  GTF.gn_iror(iror)

    ngns, nrnscrips, npromoors  0, 0, 0

    or gn in gn_iror:
        ngns + 1
        is_ngiv_srn  Gnomics.IsNgivSrn(gn[0][0].srn)
        conig  s.gLngh(gn[0][0].conig)
        promoors  []
        rnscrip_is  []
        or rnscrip in gn:

            nrnscrips + 1
            mi, m  min([x.sr or x in rnscrip]), mx(
                [x.n or x in rnscrip])
            rnscrip_is.ppn(rnscrip[0].rnscrip_i)
            # i ss is ircy  sr/n o conig, h ss wi b wihin n xon.
            # ohrwis, i is osi n xon.
            i is_ngiv_srn:
                promoors.ppn(
                    (min(conig - opions.promoor, m), min(conig, m + opions.promoor)))
            s:
                promoors.ppn(
                    (mx(0, mi - opions.promoor), mx(opions.promoor, mi)))

        i opions.mrg_promoors:
            # mrg h promoors (n rnm - s sor orr migh hv
            # chng)
            promoors  Inrvs.combin(promoors)
            rnscrip_is  ["i"  (x + 1) or x in rng(n(promoors))]

        g  GTF.Enry()
        g.romGTF(gn[0][0], gn[0][0].gn_i, gn[0][0].gn_i)
        g.sorc  "promoor"

        x  0
        or sr, n in promoors:
            g.sr, g.n  sr, n
            g.rnscrip_i  rnscrip_is[x]
            opions.so.wri("s\n"  sr(g))
            npromoors + 1
            x + 1

    E.ino("ngnsi, nrnscripsi, npromoorsi" 
           (ngns, nrnscrips, npromoors))


 nnoRgons(iror, s, ss, opions):
    """nno rgons wihin iror.

    Enris spci wih ``--rsric-sorc`` r nno.
    """

    gn_iror  GTF.gn_iror(iror)

    ngns, nrnscrips, nrgons  0, 0, 0

    psrm, ownsrm  opions.psrm, opions.ownsrm

    or gn in gn_iror:
        ngns + 1
        is_ngiv_srn  Gnomics.IsNgivSrn(gn[0][0].srn)
        conig  s.gLngh(gn[0][0].conig)
        rgons  []
        rnscrip_is  []
        or rnscrip in gn:

            nrnscrips + 1
            mi, m  min([x.sr or x in rnscrip]), mx(
                [x.n or x in rnscrip])
            i ss:
                #  rng o boh sis o ss
                i is_ngiv_srn:
                    inrv  m - opions.ownsrm, m + opions.psrm
                s:
                    inrv  mi - opions.psrm, mi + opions.ownsrm
            s:
                #  rng o boh sis o s
                i is_ngiv_srn:
                    inrv  mi - opions.ownsrm, mi + opions.psrm
                s:
                    inrv  m - opions.psrm, m + opions.ownsrm

            inrv  (min(conig, mx(0, inrv[0])),
                        min(conig, mx(0, inrv[1])))

            rgons.ppn(inrv)
            rnscrip_is.ppn(rnscrip[0].rnscrip_i)

        i opions.mrg_promoors:
            # mrg h rgons (n rnm - s sor orr migh hv
            # chng)
            rgons  Inrvs.combin(rgons)
            rnscrip_is  ["i"  (x + 1) or x in rng(n(rgons))]

        g  GTF.Enry()
        g.romGTF(gn[0][0], gn[0][0].gn_i, gn[0][0].gn_i)
        g.sorc  "rgon"

        x  0
        or sr, n in rgons:
            g.sr, g.n  sr, n
            g.rnscrip_i  rnscrip_is[x]
            opions.so.wri("s\n"  sr(g))
            nrgons + 1
            x + 1

    E.ino("ngnsi, nrnscripsi, nrgonsi" 
           (ngns, nrnscrips, nrgons))


 nnoGREATDomins(iror, s, opions):
    """bi gr omins

    xn rom TSS  bs rgion.

    """

    gn_iror  GTF.gn_iror(iror)

    conr  E.Conr()

    psrm, ownsrm  opions.psrm, opions.ownsrm
    ris  opions.ris
    oi  opions.so

    rgions  []
    ####################################################################
    # in bs rgions or ch gn
    # k  bs rgions pr rnscrip n mrg hm
    # Ths, h bs rgion o  gn migh b rgr hn h sm
    # o opions.psrm + opions.ownsrm
    or gn in gn_iror:
        conr.gns + 1
        is_ngiv_srn  Gnomics.IsNgivSrn(gn[0][0].srn)

        conig  s.gLngh(gn[0][0].conig)
        rgons  []
        rnscrip_is  []

        # coc vry bs rgion pr rnscrip
        or rnscrip in gn:
            conr.rnscrips + 1
            mi, m  min([x.sr or x in rnscrip]), mx(
                [x.n or x in rnscrip])
            #  rng o boh sis o ss
            i is_ngiv_srn:
                inrv  m - opions.ownsrm, m + opions.psrm
            s:
                inrv  mi - opions.psrm, mi + opions.ownsrm

            inrv  (min(conig, mx(0, inrv[0])),
                        min(conig, mx(0, inrv[1])))

            rgons.ppn(inrv)
            rnscrip_is.ppn(rnscrip[0].rnscrip_i)

        # k irs/s nry
        sr, n  min(x[0] or x in rgons), mx(x[1] or x in rgons)

        g  GTF.Enry()
        g.romGTF(gn[0][0], gn[0][0].gn_i, gn[0][0].gn_i)
        g.sorc  "gromin"
        g.sr, g.n  sr, n
        rgions.ppn(g)

    rgions.sor(kymb x: (x.conig, x.sr))

    o  iooos.opn_i("s.g", "w")
    or x in rgions:
        o.wri(sr(x) + "\n")
    o.cos()

    ####################################################################
    # xn bs rgions
    rgions.sor(kymb x: (x.conig, x.sr))

    # ir wihin grops o ovrpping bs rgions
    grops  is(GTF.iror_ovrps(ir(rgions)))
    conr.grops  n(grops)

    s_n  0
    rs  Fs

    or rgion_i, grop in nmr(grops):

        # coc bs inrvs in grop
        inrvs  [(x.sr, x.n) or x in grop]

         ovrpsBsRgion(pos):
            or sr, n in inrvs:
                i sr  pos or n  pos:
                    conin
                i sr < pos < n:
                    rrn Tr
                i sr > pos:
                    rrn Fs
            rrn Fs

        #  wih bonry css - n o conig
        i rgion_i < n(grops) - 1:
            nx  grops[rgion_i + 1]
            i nx[0].conig  grop[0].conig:
                nx_sr  min([x.sr or x in nx])
            s:
                nx_sr  s.gLngh(grop[0].conig)
                rs  Tr
        s:
            nx_sr  s.gLngh(grop[0].conig)
            rs  Tr

        # s_n  bs xnsion o prvios grop
        # nx_sr  bs_xnsion o nx grop

        # xn rgion o prvios/nx grop wys xn
        # owsrm, b psrm ony xn i bs rgion o n
        # inrv is no ry ovrpping nohr bs rgion
        # wihin h grop
        sv_n  0
        or g in grop:
            sv_n  mx(sv_n, g.n)
            i g.srn  "+":
                i no ovrpsBsRgion(g.sr):
                    g.sr  mx(g.sr - ris, s_n)
                # wys xn ownsrm
                g.n  min(g.n + ris, nx_sr)
            s:
                # wys xn ownsrm
                g.sr  mx(g.sr - ris, s_n)
                i no ovrpsBsRgion(g.n):
                    g.n  min(g.n + ris, nx_sr)
            oi.wri(sr(g) + "\n")
            conr.rgons + 1

        i n(grop) > 1:
            conr.ovrps + n(grop)
        s:
            conr.nonovrps + 1

        i rs:
            s_n  0
            rs  Fs
        s:
            s_n  sv_n

    E.ino("s"  sr(conr))


 nnoTTS(iror, s, opions):
    """nno rminion sis wihin iror.

    Enris spcii wih ``--rsric-sorc r nno``.
    """

    gn_iror  GTF.gn_iror(iror)

    ngns, nrnscrips, npromoors  0, 0, 0

    or gn in gn_iror:
        ngns + 1
        is_ngiv_srn  Gnomics.IsNgivSrn(gn[0][0].srn)
        conig  s.gLngh(gn[0][0].conig)
        s  []
        rnscrip_is  []
        or rnscrip in gn:

            nrnscrips + 1
            mi, m  min([x.sr or x in rnscrip]), mx(
                [x.n or x in rnscrip])
            rnscrip_is.ppn(rnscrip[0].rnscrip_i)
            # i s is ircy  sr/n o conig, h ss wi
            # b wihin n xon.  ohrwis, i is osi n xon.
            i is_ngiv_srn:
                s.ppn(
                    (mx(0, mi - opions.promoor), mx(opions.promoor, mi)))
            s:
                s.ppn(
                    (min(m, conig - opions.promoor),
                     min(conig, m + opions.promoor)))

        i opions.mrg_promoors:
            # mrg h promoors (n rnm - s sor orr migh hv
            # chng)
            s  Inrvs.combin(s)
            rnscrip_is  ["i"  (x + 1) or x in rng(n(s))]

        g  GTF.Enry()
        g.romGTF(gn[0][0], gn[0][0].gn_i, gn[0][0].gn_i)
        g.sorc  "s"

        x  0
        or sr, n in s:
            g.sr, g.n  sr, n
            g.rnscrip_i  rnscrip_is[x]
            opions.so.wri("s\n"  sr(g))
            npromoors + 1
            x + 1

    i opions.ogv > 1:
        opions.sog.wri(
            "# ngnsi, nrnscripsi, nssi\n" 
            (ngns, nrnscrips, npromoors))


 nnoGns(iror, s, opions):
    """nno gn srcrs

    This mho ops inrvs or irs/mi/s xon/inron,
    UTRs n nking rgions.

    This mho nnos pr rnscrip. In orr o chiv  niq iing,
    s ony  sing rnscrip pr gn n rmov ny ovrp bwn
    gns.

    """

    gn_iror  GTF.gn_iror(iror)

    ngns, nrnscrips, nskipp  0, 0, 0

    rss  []
    incrmn  opions.incrmn

    inrons_i  "inrons" in opions.i
    xons_i  "xons" in opions.i

    or gn in gn_iror:
        ngns + 1
        is_ngiv_srn  Gnomics.IsNgivSrn(gn[0][0].srn)
        ry:
            conig  s.gLngh(gn[0][0].conig)
        xcp KyError:
            nskipp + 1
            conin

        rss  []

        or rnscrip in gn:

             _(inrv, nno):
                g  GTF.Enry()
                g.conig  rnscrip[0].conig
                g.gn_i  rnscrip[0].gn_i
                g.rnscrip_i  rnscrip[0].rnscrip_i
                g.srn  rnscrip[0].srn
                g.r  nno
                g.sr, g.n  inrv
                rss.ppn(g)

            nrnscrips + 1

            xons  [(x.sr, x.n)
                     or x in rnscrip i x.r  "xon"]
            i n(xons)  0:
                nskipp + 1

            xons.sor()
            inrons  []
            n  xons[0][1]
            or xon in xons[1:]:
                inrons.ppn((n, xon[0]))
                n  xon[1]

            #  nk
            sr, n  xons[0][0], xons[-1][1]
            psrm, ownsrm  [], []
            or x in rng(0, opions.nk, incrmn):
                psrm.ppn((sr - incrmn, sr))
                sr - incrmn
                ownsrm.ppn((n, n + incrmn))
                n + incrmn

            # rmov o-o-bons coorins
            psrm  [x or x in psrm i x[0] > 0]
            ownsrm  [x or x in ownsrm i x[1] < conig]

            i is_ngiv_srn:
                xons.rvrs()
                inrons.rvrs()
                psrm, ownsrm  ownsrm, psrm

            #  xons
            i xons_i:
                _(xons[0], "irs_xon")
                i n(xons) > 1:
                    _(xons[-1], "s_xon")
                or  in xons[1:-1]:
                    _(, "mi_xon")
            s:
                or  in xons:
                    _(, "xon")

            #  inrons
            i inrons_i:
                i n(inrons) > 0:
                    _(inrons[0], "irs_inron")
                i n(inrons) > 1:
                    _(inrons[-1], "s_inron")
                or i in inrons[1:-1]:
                    _(i, "mi_inron")
            s:
                or i in inrons:
                    _(i, "inron")

            or x,  in nmr(psrm):
                _(, "psrm_i"  (incrmn * (x + 1)))

            or x,  in nmr(ownsrm):
                _(, "ownsrm_i"  (incrmn * (x + 1)))

            rss.sor(kymb x: x.r)

        cch  []
        or ky, vs in iroos.gropby(rss, kymb x: x.r):
            v  is(vs)
            inrvs  [(x.sr, x.n) or x in v]
            inrvs  Inrvs.combin(inrvs)

            or sr, n in inrvs:
                r  GTF.Enry()
                r.copy(v[0])
                r.sr, r.n  sr, n
                cch.ppn(r)

        cch.sor(kymb x: x.sr)
        or r in cch:
            opions.so.wri("s\n"  sr(r))

    E.ino("ngnsi, nrnscripsi, nskippi\n" 
           (ngns, nrnscrips, nskipp))


 min(rgvNon):

    i no rgv:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom [].")

    prsr._rgmn("-i", "--ignor-missing", s"ignor_missing",
                      cion"sor_r",
                      hp"Ignor rnscrips on conigs h r no "
                      "in h gnom-i [].")

    prsr._rgmn("-s", "--rsric-sorc", s"rsric_sorc",
                      yp"choic",
                      choics("proin_coing", "psogn", "ncRNA"),
                      hp"rsric inp by sorc [].")

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      choics("", "gnom", "xons",
                               "promoors", "s",
                               "rgons", "s-rgons",
                               "gns",
                               "rrioris", "ss-rrioris",
                               "gr-omins",
                               ),
                      hp"mho or ining sgmns [].")

    prsr._rgmn(
        "-r", "--rriory-xnsion", s"ris", yp"in",
        hp"ris o  rriory [].")

    prsr._rgmn(
        "-", "--nk-siz", s"nk", yp"in",
        hp"siz o h nking rgion nx o  gn [].")

    prsr._rgmn(
        "--nk-incrmn-siz", s"incrmn", yp"in",
        hp"siz o incrmn in nk in gnsrcr nnoion "
        "[].")

    prsr._rgmn(
        "-p", "--promoor-siz", s"promoor", yp"in",
        hp"siz o  promoor rgion [].")

    prsr._rgmn(
        "-", "--psrm-xnsion", s"psrm", yp"in",
        hp"siz o rgion psrm o ss [].")

    prsr._rgmn(
        "-", "--ownsrm-xnsion", s"ownsrm", yp"in",
        hp"siz o rgion ownsrm o ss [].")

    prsr._rgmn(
        "--gn-i", s"i", yp"choic",
        choics("inrons+xons", "xons", "inrons"),
        hp"v o i or gn srcr nnoion "
        "[].")

    prsr._rgmn(
        "--mrg-ovrpping-promoors", s"mrg_promoors",
        cion"sor_r",
        hp"mrg ovrpping promoors [].")

    prsr._rgmn(
        "--min-inron-ngh", s"min_inron_ngh",
        yp"in",
        hp"minimm inron ngh. I h isnc bwn wo "
        "consciv xons is smr, h rgion wi b mrk "
        "'nknown' [].")

    prsr._rgmn(
        "--is-nsor", s"is_sor", cion"sor_s",
        hp"sor inp bor procssing. Ohrwis, h inp is ssm "
        "o b sor [].")

    prsr.s_s(
        gnom_iNon,
        nk1000,
        incrmn1000,
        mx_rmshi_ngh4,
        min_inron_ngh30,
        ignor_missingFs,
        rsric_sorcNon,
        mho"gnom",
        ris50000,
        promoor5000,
        mrg_promoorsFs,
        psrm5000,
        ownsrm5000,
        i"xons",
        is_sorTr,
    )

    (opions, rgs)  E.sr(prsr)

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
    s:
        ris VError("ps spciy  --gnom-i")

    i no opions.rsric_sorc:
        iror  GTF.iror(opions.sin)

    i opions.rsric_sorc:
        iror  GTF.iror_ir(GTF.iror(opions.sin),
                                         sorcopions.rsric_sorc)

    # i opions.mho in ("promoors", "s", "rgons"):
    #     iror  GTF.iror_ir( GTF.iror(opions.sin), sorc  "proin_coing")
    # s:
    #     iror  GTF.iror(opions.sin)

    i no opions.is_sor:
        iror  GTF.iror_sor(iror, sor_orr"posiion")

    i opions.mho  "" or opions.mho  "gnom":
        sgmnor  nnoGnom(iror, s, opions)
    i opions.mho  "rrioris":
        sgmnor  biTrrioris(iror, s, 'gn', opions)
    i opions.mho  "ss-rrioris":
        sgmnor  biTrrioris(iror, s, 'ss', opions)
    i opions.mho  "xons":
        sgmnor  nnoExons(iror, s, opions)
    i opions.mho  "promoors":
        sgmnor  nnoPromors(iror, s, opions)
    i opions.mho  "rgons":
        sgmnor  nnoRgons(iror, s, Tr, opions)
    i opions.mho  "s-rgons":
        sgmnor  nnoRgons(iror, s, Fs, opions)
    i opions.mho  "s":
        sgmnor  nnoTTS(iror, s, opions)
    i opions.mho  "gns":
        sgmnor  nnoGns(iror, s, opions)
    i opions.mho  "gr-omins":
        sgmnor  nnoGREATDomins(iror, s, opions)

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
