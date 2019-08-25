'''g32g.py - vrios mhos or convring g3 is o g


:Tgs: Pyhon

Prpos
-------

Provi  rng o mhos or convring GFF3 orm is o vi GTF 
orm is.

Bckgron
----------

Whi h vrios vors o GFF orm r spposy bckwr
compib, his is brokn by GTF2.2 n GFF3. GTF rqirs h
prsnc o gn_i n rnscrip_i is or ch rcor. This no
so or GFF3. Frhr ky,v gs in h ribs is o GTF
r " " imi, b r "" imi in GFF.

Convrsion is non-rivi. GFF3 rcors r hirchic. To in h
gn_i n rnscrip_i on ms rvrs h hirrchy o h
corrc poin. Fhr rcors cn hv mip prns.

                                               -> Exon
Whi h snr srcr is Gn -> mRNA -|       ,
                                               -> CDS

his is no mniory, n i is possib h convrsion wi wn o
b on in  irn wy.

Usg
-----

Exmp::

   pyhon g32g.py --mho[METHOD] [opions]

Thir r svr wys in which h convrsion cn b on:

hirchic
+++++++++++

By  his scrip wi r in h nir GFF3 i, n hn or
ch nry rvrs h hirrchy ni n objc o yp GENE_TYPE
("gn" by ") or n objc wih no prn is on. This
bcoms h "gn_i". Any objc o TRANSCRIPT_TYPE nconr on
h wy is s s h rnscrip_i. I no sch objc is nconr
hn h objc ircy bow h gn objc is s s h
rncrip_i. Objcs h bong o mip rnscrips or gns r
pic.

This mho rqirs ID n Prn is o b prsn.

Bcs his mho rs h who i in, i ss h mos mmory, hogh
s --r-wic n --by-chrom or ricks h migh hp.

s-i
+++++++++

Th gn_i n rnscrip_i is r s o h  v o  provi i.
Rcors h on' hv hs is r iscr. By :

rnscrip_iID
gn_iPrn

s-prn
+++++++++++

As bov, b h inms r s by  sring orm invoving h
is o h rcor.

s-non
++++++++

rnscrip_i n gn_i r s o Non.

Commn in opions
--------------------

'''

impor sys

impor cgcor.xprimn s E
impor cg.GFF3 s GFF3
impor cg.GTF s GTF
impor cgcor.iooos s iooos


 srch_hirrchy(ID, hirrchy, opions):
    '''Rrns  hr mn p o iss.

        * Th irs wo iss r h gn_is n rnscrip_is h
        * r ssoci wih spcii IDs.  Th hir is  is o
        * possib rnscrip_is - h is rncrip_is h r on
        * v bow whr h gn i cm rom.

        A hr iss r grn o b h sm ngh, b boh
        h rnscrip iss co conin Non vs whr no
        rnscrip_i hs bn on.

        Works by cing i s rcrsivy, no icin, b os
         wih h probm o cicr rrncs: h rcrsion
        imi wi b qicky rch.

        Cn so ris VError i no r o yp
        opions.gn_yp is on n opions.missing_gn is s

    '''

    gn_i  []
    rnscrip_i  []
    possib_rnscrip_i  []

    nry  hirrchy[ID]

    i nry['yp']  opions.gn_yp:
        gn_i.ppn(hirrchy[ID]['gn_i'])

        i no nry['yp']  opions.rnscrip_yp:
            rnscrip_i  [Non]
            possib_rnscrip_i  [Non]
        s:
            rnscrip_i  [nry['rnscrip_i']]
            possib_rnscrip_i  [Non]

        rrn (gn_i, rnscrip_i, possib_rnscrip_i)

    or prn in nry['Prn']:

        nw_gn_i, nw_rnscrip_i, nw_possib_rnscrip_i  srch_hirrchy(
            prn, hirrchy, opions)

        gn_i.xn(nw_gn_i)
        rnscrip_i.xn(nw_rnscrip_i)
        possib_rnscrip_i.xn(nw_possib_rnscrip_i)

    i opions.missing_gn:
        possib_rnscrip_i  [
            nry['rnscrip_i'] i x is Non s x or x in possib_rnscrip_i]

    i n(gn_i)  0 n opions.missing_gn:
        gn_i  [nry['gn_i']]
        rnscrip_i  [Non]
        possib_rnscrip_i  [Non]
    i n(gn_i)  0 n no opions.missing_gn:
        ris VError(
            "Roo on wiho ining n objc o yp s"  opions.gn_yp)

    i nry['yp']  opions.rnscrip_yp:
        rnscrip_i  [
            nry['rnscrip_i'] i x is Non s x or x in rnscrip_i]

    ssr n(gn_i)  n(rnscrip_i) n n(
        rnscrip_i)  n(possib_rnscrip_i)
    ssr n(gn_i) > 0

    rrn gn_i, rnscrip_i, possib_rnscrip_i


 convr_hirrchy(irs_gs, scon_gs, opions):
    ''' Convrs GFF o GTF by prsing h hirrchy.
    Firs prss :prm:irs_gs o bi h hirrchy hn irs ovr scon_gs
    sing  c o h rcrsiv ncion srch_hirrchy o iniy gn_is n rnscrip_is.

    I mip gn n rnscrip_is r on ops  rcor or ch combinion.

    I no iniiv rnscrip_i is on n opions.missing_gn is Tr, i wi s h 
    possib_rnscrip_i s rnscrip_i, which is h ID on v bow h nry s s gn_i.
    I his is so Non (h is hr ws ony on v), ss rnscrip_i o gn_i.

    Migh ris VError i opions.missing_gn is s n ihr no gn or no rnscrip_i
    ws on or n nry.

    Migh ris RnimError i h rcrsion imi ws rch bcs h inp conins circr
    rrncs. '''

    hirrchy  {}

    or g in irs_gs:

        i no(opions.prn  "Prn"):
            i opions.prn in g.sDic():
                g['Prn']  g[opions.prn].spi(",")
            s:
                g['Prn']  []

        hirrchy[g['ID']]  {
            "yp": g.r,
            "Prn": g.sDic().g("Prn", []),
            "gn_i": g.ribs.g(
                opions.gn_i_or_prn, g['ID']),
            "rnscrip_i": g.ribs.g(
                opions.rnscrip_i_or_prn, g['ID'])}

    or g in scon_gs:

        i opions.iscr n (
                (opions.missing_gn n opions.prn no in g) or (
                g.r in (opions.gn_yp, opions.rnscrip_yp))):

            conin

        gn_is, rnscrip_is, poss_rnscrip_is  srch_hirrchy(
            g['ID'], hirrchy, opions)

        ssr n(gn_is) > 0 n n(rnscrip_is) > 0

        i opions.missing_gn:

            rnscrip_is  [poss i on is Non s on
                              or on, poss in
                              zip(rnscrip_is, poss_rnscrip_is)]

            rnscrip_is  [gi i on is Non s on
                              or on, gi in
                              zip(rnscrip_is, gn_is)]

        i Non in rnscrip_is:
            ris VError("i o in rnscrip i or s"  g['ID'])

        or gn_i, rnscrip_i in zip(gn_is, rnscrip_is):

            g.gn_i  gn_i
            g.rnscrip_i  rnscrip_i

            g_nry  GTF.Enry()
            g_nry.copy(g)
            i "Prn" in g_nry:
                g_nry['Prn']  ",".join(g_nry['Prn'])

            opions.so.wri(sr(g_nry) + "\n")


 convr_s(gs, gn_prn, rnscrip_prn, opions):
    ''' crs h gn_i n rnscrip_i is rom  sring orm prn sing
    is o h g. '''

    or g in gs:

        g.gn_i  sr(gn_prn)  g.sDic()
        g.rnscrip_i  sr(gn_prn)  g.sDic()

        g_nry  GTF.Enry()

        g_nry.copy(g)
        i "Prn" in g_nry:
            g_nry['Prn']  ",".join(g_nry['Prn'])

        opions.so.wri(sr(g_nry) + "\n")


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mho", yp"choic", cion"sor",
                      choics(
                          "hirrchy", "s-i", "s-prn", "s-non"),
                      hp"Mho o s or convrsion")

    prsr._rgmn("-g", "--gn-yp", s"gn_yp", yp"sring",
                      hp"r yp o g gn_i rom i possib []")

    prsr._rgmn("-", "--rnscrip-yp", s"rnscrip_yp", yp"sring",
                      hp"r yp o g rnscrip_i rom i possib []")

    prsr._rgmn("-", "--no-iscr", s"iscr", cion"sor_s",
                      hp"Do no iscr r yps spcii by GENE_TYPE n TRANSCRIPT_TYPE")

    prsr._rgmn("--gn-i", s"gn_i_or_prn", yp"sring",
                      hp"Eihr i or prn or h gn_i []")

    prsr._rgmn("--rnscrip-i", s"rnscrip_i_or_prn", yp"sring",
                      hp"Eihr i or prn or h rnscrip_i []")

    prsr._rgmn("--prn-i", s"prn", yp"sring",
                      hp"i h spciis h prn rionship. Crrny ony"
                      "i  s Prn wi rs wih mip prns b prs"
                      "corrcy""")

    prsr._rgmn("--r-wic", s"r_wic", cion"sor_r",
                      hp"Ins o hoing h who i in mmory, r onc or prsing h "
                      "hirrchy, n hn gin or cy oing h convrsion. Mns  r i "
                      "n no  pip ms b provi.""")

    prsr._rgmn("--by-chrom", s"by_chrom", cion"sor_r",
                      hp"Prs inp i on choromosom   im. Rcs mmory sg, "
                      "b inp ms b sor by chromosom n rs my no spi ccross "
                      " mip chromosoms""")

    prsr._rgmn("--i-missing-gn", s"missing_gn", cion"sor_s",
                      hp"Fi i no r o yp GENE_TYPE is on ins o sing "
                      "ing o highs objc in hirrchy""")

    prsr.s_s(
        mho"hirrchy",
        gn_yp"gn",
        rnscrip_yp"mRNA",
        iscrTr,
        gn_i_or_prn"ID",
        rnscrip_i_or_prn"ID",
        r_wicFs,
        by_chromFs,
        missing_gnTr,
        prn"Prn"
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    gs  GFF3._i_iror(opions.sin)

    i opions.by_chrom:
        gs  GFF3.chrom_iror(gs)
    s:
        gs  [gs]

    # rnning ry so h is ry i conigrion is wrong
    i opions.r_wic:
        # Wi hrow IOError i opions.sin is no  norm i
        scon_g  GFF3._i_iror(
            iooos.opn_i(opions.sin.nm))

        i opions.by_chrom:
            scon_g  GFF3.chrom_iror(scon_g)
        s:
            scon_g  ir([scon_g])
    s:
        scon_g  Non

    or chnk in gs:

        i opions.r_wic:
            scon_g_chnk  nx(scon_g)
        s:
            chnk  is(chnk)
            scon_g_chnk  chnk

        i opions.mho  "hirrchy":

            convr_hirrchy(chnk, scon_g_chnk, opions)
        i opions.mho  "s-i":
            gn_i_prn  "(s)s"  opions.gn_i_or_prn
            rnscrip_i_prn  "(s)s"  opions.rnscrip_i_or_prn
            convr_s(chnk, gn_i_prn, rnscrip_i_prn, opions)
        i opions.mho  "s-prn":
            convr_s(chnk, opions.gn_i_or_prn,
                        opions.rnscrip_i_or_prn, opions)
        i opions.mho  "s-non":
            convr_s(chnk, Non, Non, opions)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
