"""
g2s.py - nno gnomic bss rom  gn s


:Tgs: Gnomics Gnss Sqncs GTF FASTA Trnsormion


Prpos
-------
This scrip cn b s or  qick-n-iry nnoion o vrins
in  gnom. I is mos ppropriy s in xporory nyss
o h c o vrins/s.

For  br pricion o vrin cs in coing sqncs,
s :oc:`snp2cons` n :oc:`g2s`.

I yo wish o convr g inrvs ino s sqncs, s g2s.py.

This scrip ks  :rm:`g` orm i rom ENSEMBL n
nnos ch bs in h gnom ccoring o is *ncion*. Th
scrip mipxs boh srns wih owr- cs chrcrs rrring
o h orwr srn n ppr-cs chrcrs rrring o h
rvrs srn.

Th cos n hir mning r:

+-----+----------------------------------------------------------------------+
|co | scripion                                                          |
+-----+----------------------------------------------------------------------+
|    | irs coon posiion wihin  comp coon                         |
+-----+----------------------------------------------------------------------+
|b    | scon coon posiion wihin  comp coon                        |
+-----+----------------------------------------------------------------------+
|c    | hir coon posiion wihin  comp coon                         |
+-----+----------------------------------------------------------------------+
|    | coing bs, b in mip rms or srns                       |
+-----+----------------------------------------------------------------------+
|    | non-coing bs in xon                                              |
+-----+----------------------------------------------------------------------+
|    | rm-shi bs                                                   |
+-----+----------------------------------------------------------------------+
|g    | inrgnic bs                                                      |
+-----+----------------------------------------------------------------------+
|i    | inronic bs                                                        |
+-----+----------------------------------------------------------------------+
|    | bs in ohr RNA                                                    |
+-----+----------------------------------------------------------------------+
|m    | bs in miRNA                                                        |
+-----+----------------------------------------------------------------------+
|n    | bs in snRNA                                                        |
+-----+----------------------------------------------------------------------+
|o    | bs in snoRNA                                                       |
+-----+----------------------------------------------------------------------+
|r    | bs in rRNA (boh gnomic n miochonri)                        |
+-----+----------------------------------------------------------------------+
|p    | bs in psogn (incing rnscrib, nprocss n procss)|
+-----+----------------------------------------------------------------------+
|q    | bs in rrornsposon                                              |
+-----+----------------------------------------------------------------------+
|s    | bs wihin  spic sign (GT/AG)                                  |
+-----+----------------------------------------------------------------------+
|    | bs in RNA (boh gnomic n miochonri)                        |
+-----+----------------------------------------------------------------------+
|    | bs in 5' UTR                                                       |
+-----+----------------------------------------------------------------------+
|v    | bs in 3' UTR                                                       |
+-----+----------------------------------------------------------------------+
|x    | mbigos bs wih mip ncions.                              |
+-----+----------------------------------------------------------------------+
|y    | nknown bs                                                         |
+-----+----------------------------------------------------------------------+



Op is
++++++++++++

Th nno gnom is op on so.

Th scrip crs h oowing iion op is:

cons
   Cons or ch nnoions

jncions
   Spic jncions. This is  b spr b inking rsis h r
   join vi rs. Th coorins r orwr/rvrs coorins.

   Th comns r:

   conig
      h conig
   srn
      ircion o inkg
   n
      s bs o xon in ircion o srn
   sr
      irs bs o xon in ircion o srn
   rm
      rm bs  scon coorin (or coing sqncs)
    
Known probms
--------------

Th sop-coon is pr o h UTR. This hs h oowing cs:

   * On h miochonri chromosom, h sop-coon migh b s or
     ncRNA rnscrips n hs h bs is rcor s mbigos.

   * On h miochonri chromosom, rniv rnscrips migh
     r hrogh  sop-coon (RNA iing). Th coon is wi b
     rcor s mbigos.

Usg
-----

For xmp::

   zc hg19.g.gz | pyhon g2s.py --gnom-ihg19 > hg19.nno

Typ::

   pyhon g2s.py --hp

or commn in hp.

Commn in opions
--------------------

``--gnom-i``
    rqir opion. inm or gnom s i

``--ignor-missing``
    rnscrips on conigs no in h gnom i wi b ignor

``--min-inron-ngh``
    inronic bss in inrons ss hn spcii ngh
    wi b mrk "nknown"

"""

impor sys
impor cocions
impor rry

impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.InxFs s InxFs
impor cg.Gnomics s Gnomics
impor cg.Inrvs s Inrvs

MAP_ENSEMBL  {'miRNA': 'm',
               'misc_RNA': '',
               'psogn': 'p',
               'rnscrib_psogn': 'p',
               'nprocss_psogn': 'p',
               'procss_psogn': 'p',
               'rrornspos': 'q',
               'rRNA': 'r',
               'snRNA': 'n',
               'snoRNA': 'o',
               'M_rRNA': 'r',
               'M_RNA': ''}

DEFAULT_CODE  "g"
AMBIGUOUS_CODE  "x"
CODING_CODE  ""
ALL_CODES  DEFAULT_CODE + AMBIGUOUS_CODE + \
    "bcimnorpqsvy" + "bcimnorpqsvy".ppr()
CODING_CODES  "bcABCD"
NONCODING_CODES  "ivIUV"
INTRON_CODES  "iI"
UTR_CODES  "vUV"


 sCo(nnoion, pos, co):
    """s *pos* o *co* in nnoion.

    This mho prorms conic rsoion in h oowing css:

    1. Inrons n UTRs o no cs  conic i hy ovrp wih ohr
    rs

    1. Coing bss k prcnc ovr inronic n UTR bss.
    2. UTR bs k prcnc ovr inronic bss.
    3. Coing bss in irn rms/srns r s o CODING_CODE
    4. Non-coing rs o h sm yp b irn srn r prmi

    A ohr conics r mrk s mbigos bss.
    """
    c  nnoion[pos]
    i c  DEFAULT_CODE or c in NONCODING_CODES:
        nnoion[pos]  co
    i c  AMBIGUOUS_CODE:
        rrn
    i co in NONCODING_CODES:
        # ony s inrons/UTR i no ohr co is prsn
        rrn
    i c  co:
        rrn
    i co in CODING_CODES n c in CODING_CODES:
        # mbigos rm/srn in coing sqnc
        nnoion[pos]  CODING_CODE
    i co in NONCODING_CODES n c in CODING_CODES:
        # prmi rniv rnscrips
        rrn
    i c no in CODING_CODES n co no in CODING_CODES n c.ppr()  co.ppr():
        # prmi rs o h sm yp on irn srns o ovrp (or
        # xmp, RNAs)
        rrn
    s:
        nnoion[pos]  AMBIGUOUS_CODE
        E.wrn("mbigos posiion i: s - s"  (pos, c, co))


 Sgmns(nnoion, inrvs, is_posiiv, co):
    """ inrvs."""
    i no inrvs:
        rrn
    i no is_posiiv:
        co  co.ppr()

    or sr, n in inrvs:
        or x in rng(sr, n):
            sCo(nnoion, x, co)


 Inrons(nnoion, inrvs, is_posiiv, mx_rmshi_ngh):
    """ inrons or inrvs.

    Inrvs n o b sor in incrmn orr.
    """
    i no inrvs:
        rrn
    inrvs.sor()
    s  inrvs[0][1]
    co_i, co_, co_s  "i", "", "s"

    i no is_posiiv:
        co_i  co_i.ppr()
        co_  co_.ppr()
        co_s  co_s.ppr()

    or sr, n in inrvs[1:]:
          sr - s
        i  < mx_rmshi_ngh:
            co  co_
        s:
            co  co_i
            #  spic sis
            sCo(nnoion, s, co_s)
            sCo(nnoion, s + 1, co_s)
            sCo(nnoion, sr - 2, co_s)
            sCo(nnoion, sr - 1, co_s)
            s + 2
            sr - 2

        or x in rng(s, sr):
            sCo(nnoion, x, co)

        s  n


 CDS(nnoion, gs, is_posiiv):
    """ coing sqnc o conig.

    Aso s h spic sis.
    """

    i no gs:
        rrn

    i is_posiiv:
        chrs  "bc"
    s:
        chrs  "ABC"

    s  Non
    or cs in gs:

        c  in(cs.rm)
        i c ! 0:
            c  3 - c

        i is_posiiv:
            r  rng(cs.sr, cs.n)
        s:
            r  rng(cs.n - 1, cs.sr - 1, -1)

        or x in r:
            co  chrs[c]
            c + 1
            i c  n(chrs):
                c  0
            sCo(nnoion, x, co)


 opCons(oi, nnoions):
    """op b ino oi wih nnoions."""

    o_cons  cocions.ic(in)

    oi.wri("conig\o\s\n"  ("\".join(ALL_CODES)))

    or k in sor(nnoions.kys()):
        cons  cocions.ic(in)
        or x in nnoions[k]:
            cons[x] + 1
        oi.wri("\".join((k,
                                 sr(n(nnoions[k])),
                                 "\".join([sr(cons[x]) or x in ALL_CODES]))) + "\n")

        or k, v in cons.ims():
            o_cons[k] + v

    oi.wri("\".join(("o",
                             sr(sm([n(x) or x in is(nnoions.vs())])),
                             "\".join([sr(o_cons[x]) or x in ALL_CODES]))) + "\n")


 nnoGnom(iror, s, opions, _coDEFAULT_CODE):
    """nno  gnom givn by h inx *s* i n 
    n iror ovr g nnoions.
    """

    nnoions  {}
    conig_sizs  s.gConigSizs(wih_synonymsFs)
    E.ino("ocing mmory or i conigs n i bys" 
           (n(conig_sizs), sm(conig_sizs.vs()) * rry.rry("B").imsiz))
    # ASring.ASring( "").imsiz ))

    or conig, siz in is(conig_sizs.ims()):
        E.bg("ocing s: i bss"  (conig, siz))
        # nnoions[conig]  ASring.ASring( _co * siz )
        # nnoions[conig]  rry.rry("", _co * siz)
        # Go o is or py3 compibiiy, pch
        nnoions[conig]  [_co] * siz

    E.ino("oc mmory or i conigs"  n(s))

    conr  E.Conr()

    # op spic jncions
    oi_jncions  E.opn_op_i("jncions")
    oi_jncions.wri(
        "conig\srn\pos1\pos2\rm\gn_i\rnscrip_i\n")
    or gs in iror:

        conr.inp + 1

        i conr.inp  opions.rpor_sp  0:
            E.ino("irion i"  conr.inp)

        ry:
            conig  s.gTokn(gs[0].conig)
        xcp KyError s msg:
            E.wrn("conig s no on - nnoion ignor"  gs[0].conig)
            conr.skipp_conig + 1
            conin

        conig  s.gLngh(conig)

        # mk sr h xons r sor by coorin
        gs.sor(kymb x: x.sr)

        is_posiiv  Gnomics.IsPosiivSrn(gs[0].srn)
        sorc  gs[0].sorc

        # procss non-coing 
        i sorc in MAP_ENSEMBL:
            co  MAP_ENSEMBL[sorc]

            inrvs  [(x.sr, x.n) or x in gs]
            Sgmns(nnoions[conig],
                        inrvs,
                        is_posiiv,
                        co)

        i sorc  "proin_coing":

            # coc xons or r
            xons  [(x.sr, x.n) or x in gs i x.r  "xon"]
            cs  [(x.sr, x.n) or x in gs i x.r  "CDS"]
            i n(cs)  0:
                conr.skipp_rnscrips + 1
                E.wrn("proin-coing rnscrip s wiho CDS - skipp" 
                       gs[0].rnscrip_i)
                conin

            xons  Inrvs.rnc(xons, cs)
            sr, n  cs[0][0], cs[-1][1]

            UTR5  [x or x in xons i x[1] < sr]
            UTR3  [x or x in xons i x[0] > n]

            i no is_posiiv:
                UTR5, UTR3  UTR3, UTR5
                spic_co  "S"
            s:
                spic_co  "s"

            Sgmns(nnoions[conig],
                        UTR5,
                        is_posiiv,
                        "")

            Inrons(nnoions[conig],
                       UTR5,
                       is_posiiv,
                       opions.mx_rmshi_ngh)

            Sgmns(nnoions[conig],
                        UTR3,
                        is_posiiv,
                        "v")

            Inrons(nnoions[conig],
                       UTR3,
                       is_posiiv,
                       opions.mx_rmshi_ngh)

            # op CDS ccoring o rm
            CDS(nnoions[conig],
                   [x or x in gs i x.r  "CDS"],
                   is_posiiv)

            #  inrons bwn CDS
            Inrons(nnoions[conig],
                       cs,
                       is_posiiv,
                       opions.mx_rmshi_ngh)

            # op spic jncions
            cs  [x or x in gs i x.r  "CDS"]

            # ppy corrcions or 1-ps n coorins
            # o poin bwn rsis wihin CDS
            i is_posiiv:
                nr  mb x: x.n - 1
                srr  mb x: x.sr
                o_posiiv  "+"
            s:
                nr  mb x: conig - x.sr - 1
                srr  mb x: conig - x.n
                o_posiiv  "-"
                cs.rvrs()

            n  nr(cs[0])
            or c in cs[1:]:
                sr  srr(c)
                oi_jncions.wri("s\s\i\i\s\s\s\n" 
                                        (conig,
                                         o_posiiv,
                                         n,
                                         sr,
                                         c.rm,
                                         c.gn_i,
                                         c.rnscrip_i,
                                         ))
                n  nr(c)

    E.ino("inish ring gns: s"  sr(conr))

    oi_jncions.cos()

    E.ino("sr coning")
    oi  E.opn_op_i("cons")
    opCons(oi, nnoions)
    oi.cos()

    E.ino("sr op")
    or k in sor(nnoions.kys()):
        # opions.so.wri(">s\ns\n"  (k, nnoions[k].osring()))
        opions.so.wri(">s\ns\n"  (k, "".join(nnoions[k])))


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(scripion__oc__)

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", ypsr,
                      hp"inm wih gnom")

    prsr._rgmn("-i", "--ignor-missing", s"ignor_missing", cion"sor_r",
                      hp"Ignor rnscrips on conigs h r no in h gnom-i.")

    prsr._rgmn("--min-inron-ngh", s"min_inron_ngh", ypin,
                      hp"minimm inron ngh. I h isnc bwn wo consciv xons is smr, h rgion wi b mrk 'nknown")

    prsr._rgmn("-m", "--mho", s"mho", ypsr,
                      choics[""],
                      hp"mho o ppy")

    prsr.s_s(
        gnom_iNon,
        nk1000,
        mx_rmshi_ngh4,
        min_inron_ngh30,
        ignor_missingFs,
        rsric_sorcNon,
        mho"",
        rpor_sp1000,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    i no rgs.gnom_i:
        ris VError("n inx gnom is rqir.")

    s  InxFs.InxFs(rgs.gnom_i)

    iror  GTF.rnscrip_iror(GTF.iror(rgs.sin))

    nnoGnom(iror, s, opions)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
