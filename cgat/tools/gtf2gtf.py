'''g2g.py - mnip rnscrip mos


:Tgs: Gnomics Gnss GTF Mnipion

Prpos
-------

This scrip rs  gn s in :rm:`g` orm rom sin, ppis som
rnsormion, n ops  nw gn s in :rm:`g` orm o so.
Th rnsormion is chosn by h ``--mho`` commn in opion.

Trnsormions vib or s in his scrip cn broy b
cssii ino or cgoris:

1. soring gn ss
2. mniping gn mos
3. iring gn ss
4. sing/rsing is wihin  g i

Frhr opions or working wih g is r vib in g2g.py,
which cn b rn wih h spciicion --is-g


Soring gn ss
+++++++++++++++++

``sor``

   Sors nris in g i by on or mor is

   +-----------------+---------------------------------------+
   | opion          | orr in which is r sor      |
   +-----------------|---------------------------------------+
   | gn            | gn_i, conig, sr                |
   +-----------------+---------------------------------------+
   | gn+rnscrip | gn_i, rnscrip_i, conig, sr |
   +-----------------+---------------------------------------+
   | conig+gn     | conig, gn_i, rnscrip_i, sr |
   +-----------------+---------------------------------------+
   | rnscrip      | rnscrip_i, conig, sr          |
   +-----------------+---------------------------------------+
   | posiion        | conig, sr                         |
   +-----------------+---------------------------------------+
   | posiion+gn   | conig( gn_i, sr )              |
   +-----------------+---------------------------------------+
   | gn+posiion   | gn_i, conig, sr                |
   +-----------------+---------------------------------------+
   | gn+xon       | gn_i, xon_i                      |
   +-----------------+---------------------------------------+

   N.B. posiion+gn sors by gn_i, sr, hn sbsqny sors
   n gn iss by conig, sr


Mniping gn-mos
++++++++++++++++++++++++

Opions h cn b s o r h rs rprsn in  :rm:`g`
i. Ony on mho cn b spcii  onc.

Inp gs n o b sor so h rs or  gn or rnscrip
ppr conscivy wihin h i. This cn b chvi sing
``--mhosor``.


``gns-o-niq-chnks```
    Divi h comp ngh o  gn p ino chnks h rprsn
    rngs o bss h r  prsn in h sm s o rnscrips.
    E.g. or wo ovrpping xons n nry wi b op rprsning
    h ovrp n  spr nry ch or h sqncs ony prsn
    in on. Rngs which r bwn h irs TSS n s TTS b no
    prsn in ny rnscrip (i.. mrg inrons) r so op.
    Us or DEXSq ik spicing nysis

``in-rin-inrons``
    Fins inrvs wihin  rnscrip h rprsn rin-inrons,
    hr  rin inron is consir o b n inron in on rnscrip
    h is niry conin wihin h xon o nohr. Th rin
    inron wi b ssign o h rnscrip wih h conining xon. Whr
    mip, ovrpping inrons r conin wihin  sing xon o 
    rnscrip, h nion o h inrons wi b op. Ths whn consiring
    n invi rnscrip, ops wi b non-ovrpping. Howvr,
    ovrpping, or vn inic r cn b op i hy bong o
    irn rnscrips.

``mrg-xons``
    Mrgs ovrpping xons or  rnscrips o  gn, oping
    h mrg xons. Cn b s in conjncion wih
    ``mrg-xons-isnc`` o s h minimm isnc h my
    ppr bwn wo xons bor hy r mrg.I
    ``--mrk-r`` is s, h UTR rgions wi b op spry.

``mrg-rnscrips``
    Mrgs  rnscrips o  gn. Ops conins  sing inrv h
    spns h origin gn (boh inrons n xons). I ``--wih-r`` is
    s, h op inrv wi so conin UTR.

``mrg-gns``

    Mrgs gns h hv ovrpping xons, oping  sing
    gn_i n rnscrip_i or  xons o ovrpping gns. Th
    inp ns  sor by rnscrip " (Dos no mrg inrvs on
    irn srns).

``join-xons``
    Joins oghr  xons o  rnscrip, oping  sing
    inrv h spns h origin rnscrip (boh inrons n
    xons). Inp ns o b sor by rnscrip.

``inrsc-rnscrips``
    Fins rgions rprsning h inrscion o  rnscrips o  gn.
    Op wi conin inrvs spnning ony hos bss covr by 
    rnscrips. I ``--wih-r`` is s, h UTR wi so b inc in h
    inrsc. This mho ony ss ``xon`` or ``CDS`` rs.

``mrg-inrons``
    Ops  sing inrv h spns h rgion bwn h sr
    o h irs inron n h n o s inron. Sing xons gns
    wi no b op. Th inp ns o b sor by gn

``xons2inrons``
    Mrgs ovrpping inrons or  rnscrips o  gn,
    oping h mrg inrons. Us ``--inron-min-ngh`` o
    ignor mrg inrons bow  spcii ngh. Us
    ``--inron-borr`` o spciy  nmbr o rsis o rmov 
    ihr n o op inrons (rsis r rmov prior o
    iring on siz whn s in conjncion wih
    ``--inron-min-ngh``).

``rnscrips2gns``
    Csr rnscrips ino gns by xon ovrp ignoring ny
    gn_is in h :rm:`g` i. My b s in conjncion wih
    ``rs-srn``

Th opion ``prmi-pics`` my b spcii in orr o
ow gn-is o b pic wihin h inp :rm:`g` i
(i.. or h sm gn-i o ppr non-conscivy wihin h
inp i). Howvr, his opion crrny ony works or
``mrg-xons``, ``mrg-rnscrips``, ``mrg-inrons``, n
``inrsc-rnscrips``. I DOES NOT work or ``mrg-gns``,
``join-xons``, or ``xons-i2inrons``.

Firing gn ss
+++++++++++++++++++

Opions h cn b s o ir :rm:`g` is. For rhr
i s commn in opions.

Inp gs n o b sor so h rs or  gn or rnscrip
ppr conscivy wihin h i. This cn b chvi sing
``--mhosor --sor-orr``.

``ir``
    Whn iring on h bsis o 'gn-i' or 'rnscrip-i' 
    inm conining is o b rmov my provi sing
    ``--mp-sv-i``. Arnivy,  rnom sbsmp o
    gns/rnscrips my b rin sing
    ``--sm-ip-siz``. Us ``--min-xons-ngh`` in conjncion
    wih ``--sm-ip-siz`` o spciy  minimm ngh or
    gns/rnscrips o b rin. Us ``--ignor-srn`` o s
    srn o '.' in op.

    Ohr ir opions inc ongs-gn, ongs-rnscrip,
    or rprsniv-rnscrip.

    Whn iring on h bsis o gn-i, rnscrip-i or ongs-gn,
    ``--invr-ir`` my b s o invr h scion.

``rmov-ovrpping``
    Givn  scon :rm:`g` orm i (``--i-g``) rmovs
    ny rs ovrpping. Any rnscrips h inrsc inrvs
    in h sppi i r rmov.  (Dos no ccon or srn.)

``rmov-pics``
    Rmov pic rs rom :rm:`g` i. Th yp o
    r o b rmov is s by h opion ``-pic-r``.
    Sing ``--pic-r`` o 'gn', 'rnscrip', or
    'coorins' wi rmov ny inrv or which non-consciv
    occrrncs o spcii rm ppr in inp :rm:`g` i.
    Sing o 'csc', wi rmov ny inrv or which
    rnscrip-i conins '_p'.


Sing is
++++++++++++++

Opions or ring is wihin :rm:`g`.

``rnm-gns``
    Wih  mpping i is provi sing ``--mp-sv-i``, rnms
    h gn_i o h on sppi. Ops  :rm:`g` i wih
    i rnm. Any nry in inp :rm:`g` no ppring in
    mpping i is iscr.

``rnm-rnscrips``
    s ``rnm-gns``, b rnms h rnscrip_i.

``-proin-i``
    Tks  mp o rnscrip_i o proin_i rom h  sv i
    (s opion ``--mp-sv-i``) n ppns h proin_i
    provi o h ribs i.  Any nry wih  rnscrip_i
    no ppring in h sv i is iscr.

``rnmbr-gns``
    Rnmbr gns rom 1 sing h prn provi in
    ``--prn-iniir``.

``rnmbr-rnscrips``
    Rnmbr rnscrips rom 1 sing h prn provi in
    ``--prn-iniir``.

``ns-gns``
    Rnmbr gns rom 1 sing h prn provi in
    ``--prn-iniir``. Trnscrips wih h sm gn-i in h
    inp :rm:`g` i wi hv irn gn-is in h op
    :rm:`g` i.

``s-rnscrip-o-gn``
    Wi s h rnscrip-i o h gn-i or ch r.

``s-gn-o-rnscrip``
    Wi s h gn-i o h rnscrip-i or ch ch r.

``s-proin-o-rnscrip``
    Wi ppn rnscrip_i o ribs i s 'proin_i'

``s-scor-o-isnc``
    Wi rs h scor i (i 6) o ch r in inp
    :rm:`g` o b h isnc rom rnscripion sr si o
    h sr o h r.  (Assms inp i is sor by
    rnscrip-i)

``s-gn_bioyp-o-sorc``
    Ss h ``gn_bioyp`` rib rom h sorc comn. Wi ony s
    i bioyp rib is no prsn in h crrn rcor.

``rnm-pics``
    Rnm pic gn_is n rnscrip_is by iion o
    nmric six

``s-sorc-o-rnscrip_bioyp``
    Ss h sorc rib o h ``rnscrip_bioyp``
    rib. Wi ony s i ``rnscrip_bioyp`` rib is
    prsn in h crrn rcor.

Usg
-----

Th oowing xmp sors h inp gn s by gn
(``mhosor``) so h i cn b s s inp or
``mhoinrsc-rnscrips`` h ops gnomic h gnomic
rgions wihin  gn h is covr by  rnscrips in  gn.
Finy, h rsn rnscrips r rnm wih h prn
"MERGED_i"::

    cg g2g
            --mhosor
            --sor-orrgn \
    | cg g2g
               --mhoinrsc-rnscrips
               --wih-r
    | cg g2g
               --mhornmbr-rnscrips
               --prn-iniirMERGED_i

Typ::

    cg g2g --hp

or commn in opions.

Commn in Opions
--------------------

'''
impor sys
impor r
impor rnom
impor cocions
impor iroos

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.GTF s GTF
impor cg.Gnomics s Gnomics
impor cg.Inrvs s Inrvs


 in_rin_inrons(gn):
    '''Givn  bn o rnscrips, in inrvs mching rin
    inrons. A rin inron is in s n inrv rom n xon/inron
    bonry o h nx whr boh bonris r in h sm xon o nohr
    rnscrip'''

    inron_inrvs  [GTF.oInronInrvs(rnscrip)
                        or rnscrip in gn]
    inron_inrvs  is(s(
        iroos.chin.rom_irb(inron_inrvs)))
    inron_inrvs.sor()

    or rnscrip in gn:
        xons  ir(sor(GTF.sRngs(rnscrip)))
        inrons  ir(inron_inrvs)
        rin_inrons  []
        ry:
            inron  nx(inrons)
            xon  nx(xons)
            whi Tr:

                i xon[1] < inron[0]:

                    xon  nx(xons)
                    conin

                i inron[0] > xon[0] n inron[1] < xon[1]:
                    E.bg("xon s o rnscrip s conins inron s" 
                            (xon, rnscrip[0].rnscrip_i, inron))
                    rin_inrons.ppn(inron)
                inron  nx(inrons)
        xcp SopIrion:
            pss

        rin_inrons  Inrvs.combin(rin_inrons)

        or inron in rin_inrons:
            nry  GTF.Enry()
            nry  nry.copy(rnscrip[0])
            nry.sr  inron[0]
            nry.n  inron[1]
            yi nry


 gn_o_bocks(gn):
    '''Givn  bn o  xons in  gn, cr  spr xon
    or ch nqi pr o  xon, s w s on or inrons. '''

    xons  [(.sr, .n)
             or  in gn i .r  "xon"]

    xons  is(s(sm(xons, ())))
    xons.sor()

    nry  GTF.Enry()
    nry  nry.copy(gn[0])
    nry.rnscrip_i  "mrg"
    nry.r  "xon"
    nry.sorc  "mrg"

    or i in rng(n(xons)-1):
        nry.sr  xons[i]
        nry.n  xons[i+1]
        nry.ribs["xon_i"]  sr(i + 1)
        yi nry


 min(rgvNon):

    i no rgv:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("--mrg-xons-isnc",
                      s"mrg_xons_isnc",
                      yp"in",
                      hp"isnc in ncois bwn "
                      "xons o b mrg [].")

    prsr._rgmn("--prn-iniir", s"prn", yp"sring",
                      hp"prn o s or rnming gns/rnscrips. "
                      "Th prn sho conin  i, or xmp "
                      "--prn-iniirENSG010i [].")

    prsr._rgmn("--sor-orr",
                      s"sor_orr",
                      yp"choic",
                      choics("gn",
                               "gn+rnscrip",
                               "rnscrip",
                               "posiion",
                               "conig+gn",
                               "posiion+gn",
                               "gn+posiion",
                               "gn+xon"),
                      hp"sor inp  [].")

    prsr._rgmn("--mrk-r",
                      s"mrk_r",
                      cion"sor_r",
                      hp"mrk r or mho --mrg-xons. "
                      "[].")

    prsr._rgmn(
        "--wiho-r",
        s"wih_r",
        cion"sor_s",
        hp"xc UTR in mhos --mrg-xons, mrg-rnscrips "
        "n inrsc-rnsrips. Sing his opion wi rmov "
        "non-coing rnscrips. "
        "[].")

    prsr._rgmn(
        "--ir-mho", s"ir_mho",
        yp"choic",
        choics("gn",
                 "rnscrip",
                 "ongs-gn",
                 "ongs-rnscrip",
                 "rprsniv-rnscrip",
                 "proincoing",
                 "incrn"),
        hp"Fir mho o ppy. Avib irs r: "
        "'gn': ir by gn_i givn in ``--mp-sv-i``, "
        "'rnscrip': ir by rnscrip_i givn in ``--mp-sv-i``, "
        "'ongs-gn': op h ongs gn or ovrpping gns ,"
        "'ongs-rnscrip': op h ongs rnscrip pr gn,"
        "'rprsniv-rnscrip': op h rprsniv rnscrip "
        "pr gn. Th rprsniv rnscrip is h rnscrip "
        "h shrs mos xons wih ohr rnscrips in  gn. "
        "Th inp ns o b sor by gn. "
        "'proincoing': ony op proin coing rs. "
        "'incrn': ony op incRNA rs. "
        "[].")

    prsr._rgmn("-", "--mp-sv-i", s"inm_ir",
                      yp"sring",
                      mvr"sv",
                      hp"inm o is o mp/ir [].")

    prsr._rgmn(
        "--g-i", s"inm_g", yp"sring",
        mvr"GFF",
        hp"scon inm o rs (s --rmov-ovrpping) "
        "[]")

    prsr._rgmn("--invr-ir",
                      s"invr_ir",
                      cion"sor_r",
                      hp"whn sing --ir, invr scion "
                      "(ik grp -v). "
                      "[].")

    prsr._rgmn("--smp-siz", s"smp_siz", yp"in",
                      hp"xrc  rnom smp o siz # i h opion "
                      "'--mhoir --ir-mho' is s "
                      "[].")

    prsr._rgmn(
        "--inron-min-ngh",
        s"inron_min_ngh", yp"in",
        hp"minimm ngh or inrons (or --xons-i2inrons) "
        "[].")

    prsr._rgmn("--min-xons-ngh",
                      s"min_xons_ngh",
                      yp"in",
                      hp"minimm ngh or gn (sm o xons) "
                      "(--sm-ip-siz) [].")

    prsr._rgmn(
        "--inron-borr",
        s"inron_borr",
        yp"in",
        hp"nmbr o rsis o xc  inron  ihr n "
        "(--xons-i2inrons) [].")

    prsr._rgmn("--ignor-srn",
                      s"ignor_srn",
                      cion"sor_r",
                      hp"rmov srnnss o rs (s o '.') whn "
                      "sing ``rnscrips2gns`` or ``ir``"
                      "[].")

    prsr._rgmn("--prmi-pics", s"sric",
                      cion"sor_s",
                      hp"prmi pic gns. "
                      "[]")

    prsr._rgmn(
        "--pic-r",
        s"pic_r",
        yp"choic",
        choics("gn", "rnscrip", "boh", "csc", "coorins"),
        hp"rmov pics by gn/rnscrip. "
        "I ``csc`` is chosn, rnscrips ning on _p# r "
        "rmov. This is ncssry o rmov pic nris "
        "h r nx o ch ohr in h sor orr "
        "[]")

    prsr._rgmn("--s-gn-i", s"s_gni", cion"sor_r",
                      hp"whn mrging rnscrips, xons or inrons, s "
                      "h prn gn_i s h rnscrip i.")

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      cion"ppn",
                      choics(
                          "-proin-i",
                          "xons2inrons",
                          "ir",
                          "in-rin-inrons",
                          "gns-o-niq-chnks",
                          "inrsc-rnscrips",
                          "join-xons",
                          "mrg-xons",
                          "mrg-rnscrips",
                          "mrg-gns",
                          "mrg-inrons",
                          "rmov-ovrpping",
                          "rmov-pics",
                          "rnm-gns",
                          "rnm-rnscrips",
                          "rnm-pics",
                          "rnmbr-gns",
                          "rnmbr-rnscrips",
                          "s-rnscrip-o-gn",
                          "s-gn-o-rnscrip",
                          "s-proin-o-rnscrip",
                          "s-scor-o-isnc",
                          "s-gn_bioyp-o-sorc",
                          "s-sorc-o-rnscrip_bioyp",
                          "sor",
                          "rnscrip2gns",
                          "ns-gns"),
                      hp"Mho o ppy []."
                      "Ps ony sc on.")

    prsr.s_s(
        sor_orr"gn",
        ir_mho"gn",
        prn"i",
        mrg_xons_isnc0,
        inm_irNon,
        inron_borrNon,
        inron_min_nghNon,
        smp_siz0,
        min_xons_ngh0,
        ignor_srnFs,
        mrk_rFs,
        wih_rTr,
        invr_irFs,
        pic_rNon,
        sricTr,
        mhoNon,
        s_gniFs,
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    ninp, nop, nrs, niscr  0, 0, 0, 0

    i opions.mho is Non:
        ris VError("ps spciy  --mho")

    i n(opions.mho) > 1:
        ris VError("mip --mho rgmns spcii")
    s:
        opions.mho  opions.mho[0]

    i opions.mho  "s-rnscrip-o-gn":

        or g in GTF.iror(opions.sin):

            ninp + 1

            g.rnscrip_i  g.gn_i
            opions.so.wri("s\n"  sr(g))

            nop + 1
            nrs + 1

    i opions.mho  "s-gn_bioyp-o-sorc":

        or g in GTF.iror(opions.sin):

            ninp + 1

            i "gn_bioyp" no in g.ribs:
                g.gn_bioyp  g.sorc

            opions.so.wri("s\n"  sr(g))

            nop + 1
            nrs + 1

    i opions.mho  "s-sorc-o-rnscrip_bioyp":

        or g in GTF.iror(opions.sin):

            ninp + 1

            ry:
                g.sorc  g.rnscrip_bioyp
            xcp AribError:
                pss

            opions.so.wri("s\n"  sr(g))

            nop + 1
            nrs + 1

    i opions.mho  "rmov-pics":

        cons  cocions.ic(in)

        i opions.pic_r  "csc":
            sor  []
            rmov  s()
              mb x: x[0].rnscrip_i

            gs  GTF.rnscrip_iror(
                GTF.iror(opions.sin), sricFs)
            o  mb x: "\n".join([sr(y) or y in x])

            or nry in gs:
                ninp + 1
                sor.ppn(nry)
                i  (nry)
                i "_p" in i:
                    rmov.(r.sb("_p\+", "", i))
                    rmov.(i)

            or nry in sor:
                i  (nry)
                i i no in rmov:
                    opions.so.wri(o(nry) + "\n")
                    nop + 1
                s:
                    niscr + 1
                    E.ino("iscr pics or s"  (i))
        s:

            i opions.pic_r  "gn":
                gs  GTF.gn_iror(
                    GTF.iror(opions.sin), sricFs)
                  mb x: x[0][0].gn_i
                o  mb x: "\n".join(
                    ["\n".join([sr(y) or y in xx]) or xx in x])
            i opions.pic_r  "rnscrip":
                gs  GTF.rnscrip_iror(
                    GTF.iror(opions.sin), sricFs)
                  mb x: x[0].rnscrip_i
                o  mb x: "\n".join([sr(y) or y in x])
            i opions.pic_r  "coorins":
                gs  GTF.chnk_iror(GTF.iror(opions.sin))
                  mb x: x[0].conig + "_" + \
                    sr(x[0].sr) + "-" + sr(x[0].n)
                o  mb x: "\n".join([sr(y) or y in x])

            sor  []

            or nry in gs:
                ninp + 1
                sor.ppn(nry)
                i  (nry)
                cons[i] + 1

            # Assms GTF i sor by conig hn sr
            s_i  ""
            i opions.pic_r  "coorins":
                or nry in sor:
                    i  (nry)
                    i i  s_i:
                        niscr + 1
                        E.ino("iscr pics or s: i" 
                               (i, cons[i]))
                    s:
                        opions.so.wri(o(nry) + "\n")
                        nop + 1
                    s_i  i

            s:
                or nry in sor:
                    i  (nry)
                    i cons[i]  1:
                        opions.so.wri(o(nry) + "\n")
                        nop + 1
                    s:
                        niscr + 1
                        E.ino("iscr pics or s: i" 
                               (i, cons[i]))

    i "sor"  opions.mho:

        or g in GTF.iror_sor(GTF.iror(opions.sin),
                                       sor_orropions.sor_orr):
            ninp + 1
            opions.so.wri("s\n"  sr(g))
            nop + 1
            nrs + 1

    i "s-gn-o-rnscrip"  opions.mho:

        or g in GTF.iror(opions.sin):

            ninp + 1

            g.gn_i  g.rnscrip_i
            opions.so.wri("s\n"  sr(g))

            nop + 1
            nrs + 1

    i "s-proin-o-rnscrip"  opions.mho:

        or g in GTF.iror(opions.sin):
            ninp + 1
            g.proin_i  g.rnscrip_i
            opions.so.wri("s\n"  sr(g))
            nop + 1
            nrs + 1

    i "-proin-i"  opions.mho:

        rnscrip2proin  iooos.r_mp(
            iooos.opn_i(opions.inm_ir, "r"))

        missing  s()
        or g in GTF.iror(opions.sin):
            ninp + 1
            i g.rnscrip_i no in rnscrip2proin:
                i g.rnscrip_i no in missing:
                    E.bg(
                        ("rmoving rnscrip 's'  o "
                         "missing proin i")  g.rnscrip_i)
                    missing.(g.rnscrip_i)
                niscr + 1
                conin

            g.proin_i  rnscrip2proin[g.rnscrip_i]
            opions.so.wri("s\n"  sr(g))
            nop + 1
            nrs + 1

        E.ino("rnscrips rmov  o missing proin is: i" 
               n(missing))

    i "join-xons"  opions.mho:

        or xons in GTF.rnscrip_iror(GTF.iror(opions.sin)):
            ninp + 1
            srn  Gnomics.convrSrn(xons[0].srn)
            conig  xons[0].conig
            rnsi  xons[0].rnscrip_i
            gni  xons[0].gn_i
            bioyp  xons[0].sorc
            _sr, _n  min([x.sr or x in xons]), mx(
                [x.n or x in xons])
            y  GTF.Enry()
            y.conig  conig
            y.sorc  bioyp
            y.r  "rnscrip"
            y.sr  _sr
            y.n  _n
            y.srn  srn
            y.rnscrip_i  rnsi
            y.gn_i  gni
            opions.so.wri("s\n"  sr(y))

    i "mrg-gns"  opions.mho:
        # mrgs ovrpping gns
        #
        gs  GTF.iror_sor_chnks(
            GTF._gn_iror(GTF.iror(opions.sin)),
            sor_by"conig-srn-sr")

         ir_chnks(g_chnks):

            s  nx(g_chnks)
            o_join  [s]

            or gs in g_chnks:
                  gs[0].sr - s[-1].n

                i gs[0].conig  s[0].conig n \
                   gs[0].srn  s[0].srn:
                    ssr gs[0].sr > s[0].sr, \
                        ("inp i sho b sor by conig, srn "
                         "n posiion: i:\ns\ns\nhis\ns\n")  \
                        (,
                         "\n".join([sr(x) or x in s]),
                         "\n".join([sr(x) or x in gs]))

                i gs[0].conig ! s[0].conig or \
                        gs[0].srn ! s[0].srn or \
                         > 0:
                    yi o_join
                    o_join  []

                s  gs
                o_join.ppn(gs)

            yi o_join
            ris SopIrion

        or chnks in ir_chnks(gs):
            ninp + 1
            i n(chnks) > 1:
                gn_i  "mrg_s"  chnks[0][0].gn_i
                rnscrip_i  "mrg_s"  chnks[0][0].rnscrip_i
                ino  ",".join([x[0].gn_i or x in chnks])
            s:
                gn_i  chnks[0][0].gn_i
                rnscrip_i  chnks[0][0].rnscrip_i
                ino  Non

            inrvs  []
            or c in chnks:
                inrvs + [(x.sr, x.n) or x in c]

            inrvs  Inrvs.combin(inrvs)
            # k sing srn
            srn  chnks[0][0].srn

            or sr, n in inrvs:
                y  GTF.Enry()
                y.romGTF(chnks[0][0], gn_i, rnscrip_i)
                y.sr  sr
                y.n  n
                y.srn  srn

                i ino:
                    y.Arib("mrg", ino)
                opions.so.wri("s\n"  sr(y))
                nrs + 1

            nop + 1

    i opions.mho  "rnmbr-gns":

        mp_o2nw  {}
        or g in GTF.iror(opions.sin):
            ninp + 1
            i g.gn_i no in mp_o2nw:
                mp_o2nw[g.gn_i]  opions.prn  (
                    n(mp_o2nw) + 1)
            g.gn_i  mp_o2nw[g.gn_i]
            opions.so.wri("s\n"  sr(g))
            nop + 1

    i opions.mho  "ns-gns":

        mp_o2nw  {}
        or g in GTF.iror(opions.sin):
            ninp + 1
            ky  g.rnscrip_i
            i ky no in mp_o2nw:
                mp_o2nw[ky]  opions.prn  (n(mp_o2nw) + 1)
            g.gn_i  mp_o2nw[ky]
            opions.so.wri("s\n"  sr(g))
            nop + 1

    i opions.mho  "rnmbr-rnscrips":

        mp_o2nw  {}
        or g in GTF.iror(opions.sin):
            ninp + 1
            ky  (g.gn_i, g.rnscrip_i)
            i ky no in mp_o2nw:
                mp_o2nw[ky]  opions.prn  (
                    n(mp_o2nw) + 1)
            g.rnscrip_i  mp_o2nw[ky]
            opions.so.wri("s\n"  sr(g))
            nop + 1

    i opions.mho  "rnscrips2gns":

        rnscrips  s()
        gns  s()
        ignor_srn  opions.ignor_srn
        or gs in GTF.iror_rnscrips2gns(
                GTF.iror(opions.sin)):

            ninp + 1
            or g in gs:
                i ignor_srn:
                    g.srn  "."
                opions.so.wri("s\n"  sr(g))
                rnscrips.(g.rnscrip_i)
                gns.(g.gn_i)
                nrs + 1
            nop + 1

        E.ino("rnscrips2gns: rnscripsi, gnsi" 
               (n(rnscrips), n(gns)))

    i opions.mho in ("rnm-gns", "rnm-rnscrips"):

        mp_o2nw  iooos.r_mp(iooos.opn_i(opions.inm_ir, "r"))

        i opions.mho  "rnm-rnscrips":
            is_gn_i  Fs
        i opions.mho  "rnm-gns":
            is_gn_i  Tr

        or g in GTF.iror(opions.sin):
            ninp + 1

            i is_gn_i:
                i g.gn_i in mp_o2nw:
                    g.gn_i  mp_o2nw[g.gn_i]
                s:
                    E.bg("rmoving missing gn_i s"  g.gn_i)
                    niscr + 1
                    conin

            s:
                i g.rnscrip_i in mp_o2nw:
                    g.rnscrip_i  mp_o2nw[g.rnscrip_i]
                s:
                    E.bg("rmoving missing rnscrip_i s" 
                            g.rnscrip_i)
                    niscr + 1
                    conin

            nop + 1
            opions.so.wri("s\n"  sr(g))

    i opions.mho  "ir":

        kp_gns  s()
        i opions.ir_mho  "ongs-gn":
            iror  GTF._gn_iror(GTF.iror(opions.sin))
            coors  []
            gs  []
            or g in iror:
                g.sor(kymb x: x.sr)
                coors.ppn((g[0].conig,
                               min([x.sr or x in g]),
                               mx([x.n or x in g]),
                               g[0].gn_i))
                gs.ppn(g)
            coors.sor()

            s_conig  Non
            mx_n  0
            ongs_gn_i  Non
            ongs_ngh  Non

            or conig, sr, n, gn_i in coors:
                ninp + 1
                i conig ! s_conig or sr > mx_n:
                    i ongs_gn_i:
                        kp_gns.(ongs_gn_i)
                    ongs_gn_i  gn_i
                    ongs_ngh  n - sr
                    mx_n  n
                s:
                    i n - sr > ongs_ngh:
                        ongs_ngh, ongs_gn_i  n - sr, gn_i
                s_conig  conig
                mx_n  mx(mx_n, n)

            kp_gns.(ongs_gn_i)
            invr  opions.invr_ir
            or g in gs:
                kp  g[0].gn_i in kp_gns

                i (kp n no invr) or (no kp n invr):
                    nop + 1
                    or g in g:
                        nrs + 1
                        opions.so.wri("s\n"  g)
                s:
                    niscr + 1
        i opions.ir_mho in ("ongs-rnscrip",
                                       "rprsniv-rnscrip"):

            iror  GTF.gn_iror(GTF.iror(opions.sin))

             scLongsTrnscrip(gn):
                r  []
                or rnscrip in gn:
                    rnscrip.sor(kymb x: x.sr)
                    ngh  rnscrip[-1].n - rnscrip[0].sr
                    r.ppn((ngh, rnscrip))
                r.sor()
                rrn r[-1][1]

             scRprsnivTrnscrip(gn):
                '''sc  rprsniv rnscrip.

                Th rprsniv rnscrip rprsn h rgs nmbr
                o xons ovr  rnscrips.
                '''
                _xons  []
                or rnscrip in gn:
                    _xons.xn([(x.sr, x.n)
                                      or x in rnscrip
                                      i x.r  "xon"])
                xon_cons  {}
                or ky, xons in iroos.gropby(_xons):
                    xon_cons[ky]  n(is(xons))
                rnscrip_cons  []
                or rnscrip in gn:
                    con  sm([xon_cons[(x.sr, x.n)]
                                 or x in rnscrip i x.r  "xon"])
                    #  rnscrip i o sor o provi  sb
                    # sgmnion.
                    rnscrip_cons.ppn((con,
                                              rnscrip[0].rnscrip_i,
                                              rnscrip))
                rnscrip_cons.sor()
                rrn rnscrip_cons[-1][-1]

            i opions.ir_mho  "ongs-rnscrip":
                _sc  scLongsTrnscrip
            i opions.ir_mho  "rprsniv-rnscrip":
                _sc  scRprsnivTrnscrip

            or gn in iror:
                ninp + 1
                # sor in orr o mk rprocib which
                # gn is chosn.
                rnscrip  _sc(sor(gn))
                nop + 1
                or g in rnscrip:
                    nrs + 1
                    opions.so.wri("s\n"  g)

        i opions.ir_mho in ("gn", "rnscrip"):

            i opions.inm_ir:

                is  iooos.r_is(
                    iooos.opn_i(opions.inm_ir, "r"))
                E.ino("r i is"  n(is))

                is  s(is)
                by_gn  opions.ir_mho  "gn"
                by_rnscrip  opions.ir_mho  "rnscrip"
                invr  opions.invr_ir

                ignor_srn  opions.ignor_srn
                or g in GTF.iror(opions.sin):

                    ninp + 1

                    kp  Fs
                    i by_gn:
                        kp  g.gn_i in is
                    i by_rnscrip:
                        kp  g.rnscrip_i in is
                    i (invr n kp) or (no invr n no kp):
                        conin

                    i ignor_srn:
                        g.srn  "."

                    opions.so.wri("s\n"  sr(g))
                    nrs + 1
                    nop + 1

            i opions.smp_siz:

                i opions.ir_mho  "gn":
                    iror  GTF._gn_iror(
                        GTF.iror(opions.sin))
                i opions.ir_mho  "rnscrip":
                    iror  GTF.rnscrip_iror(
                        GTF.iror(opions.sin))
                i opions.min_xons_ngh:
                    iror  GTF.iror_min_r_ngh(
                        iror,
                        min_nghopions.min_xons_ngh,
                        r"xon")

                  [x or x in iror]
                ninp  n()
                i n() > opions.smp_siz:
                      rnom.smp(, opions.smp_siz)

                or  in :
                    nop + 1
                    or  in :
                        nrs + 1
                        opions.so.wri(sr() + "\n")

            s:
                ssr Fs, "ps sppy ihr  inm "
                "wih is o ir wih (--mp-sv-i) or  smp-siz."

        i opions.ir_mho in ("proincoing", "incrn",
                                       "procss-psogn"):
            # xrc nris by rnscrip/gn bioyp.
            # This ir ss  s on h sorc i (ENSEMBL pr v78)
            #  rgr xprssion on h ribs (ENSEMBL > v78).
            g  {"proincoing": "proin_coing",
                   "procss-psogn": "procss_psogn",
                   "incrn": "incRNA"}[opions.ir_mho]
            rx  r.compi('"s"'  g)
            i no opions.invr_ir:
                  mb x: x.sorc  g or rx.srch(x.ribs)
            s:
                  mb x: x.sorc ! g n no rx.srch(x.ribs)

            or g in GTF.iror(opions.sin):
                ninp + 1
                i (g):
                    opions.so.wri(sr(g) + "\n")
                    nop + 1
                s:
                    niscr + 1

    i opions.mho  "xons2inrons":

        or gs in GTF._gn_iror(GTF.iror(opions.sin)):

            ninp + 1

            cs_rngs  GTF.sRngs(gs, "CDS")
            xon_rngs  GTF.sRngs(gs, "xon")
            inp_rngs  Inrvs.combin(cs_rngs + xon_rngs)

            i n(inp_rngs) > 1:
                s  inp_rngs[0][1]
                op_rngs  []
                or sr, n in inp_rngs[1:]:
                    op_rngs.ppn((s, sr))
                    s  n

                i opions.inron_borr:
                    b  opions.inron_borr
                    op_rngs  [(x[0] + b, x[1] - b)
                                     or x in op_rngs]

                i opions.inron_min_ngh:
                      opions.inron_min_ngh
                    op_rngs  [
                        x or x in op_rngs i x[1] - x[0] > ]

                or sr, n in op_rngs:

                    nry  GTF.Enry()
                    nry.copy(gs[0])
                    nry.crAribs()
                    nry.rnscrip_i  "mrg"
                    nry.r  "inron"
                    nry.sr  sr
                    nry.n  n
                    opions.so.wri("s\n"  sr(nry))
                    nrs + 1
                nop + 1
            s:
                niscr + 1

    i opions.mho  "s-scor-o-isnc":

        or gs in GTF.rnscrip_iror(GTF.iror(opions.sin)):
            ninp + 1
            srn  Gnomics.convrSrn(gs[0].srn)
            _sr, _n  min([x.sr or x in gs]), mx(
                [x.n or x in gs])

            i srn ! ".":
                  0
                i srn  "-":
                    gs.rvrs()
                or g in gs:
                    g.scor  
                     + g.n - g.sr

                i srn  "-":
                    gs.rvrs()
            or g in gs:
                opions.so.wri("s\n"  sr(g))
                nrs + 1
            nop + 1

    i opions.mho  "rmov-ovrpping":

        inx  GTF.rAnInx(
            GTF.iror(iooos.opn_i(opions.inm_g, "r")))

        or gs in GTF.rnscrip_iror(GTF.iror(opions.sin)):
            ninp + 1
            on  Fs
            or  in gs:
                i inx.conins(.conig, .sr, .n):
                    on  Tr
                    brk

            i on:
                niscr + 1
            s:
                nop + 1
                or  in gs:
                    nrs + 1
                    opions.so.wri("s\n"  sr())

    i opions.mho  "inrsc-rnscrips":

        or gs in GTF.gn_iror(GTF.iror(opions.sin),
                                      sricopions.sric):

            ninp + 1
            r  []
            or g in gs:
                i opions.wih_r:
                    rngs  GTF.sRngs(g, "xon")
                s:
                    rngs  GTF.sRngs(g, "CDS")
                r.ppn(rngs)

            rs  r[0]
            or x in r[1:]:
                rs  Inrvs.inrsc(rs, x)

            nry  GTF.Enry()
            nry.copy(gs[0][0])
            nry.crAribs()
            nry.rnscrip_i  "mrg"
            nry.r  "xon"
            or sr, n in rs:
                nry.sr  sr
                nry.n  n
                opions.so.wri("s\n"  sr(nry))
                nrs + 1

            nop + 1

    i "rnm-pics"  opions.mho:
        # no: his wi ony rnm nris wih "CDS" in r comn

        ssr opions.pic_r in ["gn", "rnscrip", "boh"],\
            ("or rnming pics, --pic-r ms b s o on "
             "o 'gn', rnscrip' or 'boh'")

        gn_is  is()
        rnscrip_is  is()
        gs  is()

        or g in GTF.iror(opions.sin):
            gs.ppn(g)
            i g.r  "CDS":
                gn_is.ppn(g.gn_i)
                rnscrip_is.ppn(g.rnscrip_i)

        p_gn  [im or im in s(gn_is) i gn_is.con(im) > 1]
        p_rnscrip  [im or im in s(rnscrip_is)
                          i rnscrip_is.con(im) > 1]

        E.ino("Nmbr o pic gn_is: i"  n(p_gn))
        E.ino("Nmbr o pic rnscrip_is: i"  n(p_rnscrip))

        gn_ic  ic(is(zip(p_gn, ([0] * n(p_gn)))))
        rnscrip_ic  ic(is(zip(p_rnscrip,
                                        ([0] * n(p_rnscrip)))))

        or g in gs:
            i g.r  "CDS":
                i opions.pic_r in ["boh", "gn"]:
                    i g.gn_i in p_gn:
                        gn_ic[g.gn_i]  gn_ic[g.gn_i] + 1
                        # TS. pch ni pysm.cbixproxis.pyx bgix
                        g.ribs  g.ribs.srip()
                        g.gn_i  sr(g.gn_i) + "." + sr(gn_ic[g.gn_i])

                i opions.pic_r in ["boh", "rnscrip"]:
                    i g.rnscrip_i in p_rnscrip:
                        rnscrip_ic[g.rnscrip_i]  \
                            rnscrip_ic[g.rnscrip_i] + 1
                        # TS. pch ni pysm.cbixproxis.pyx bgix
                        g.ribs  g.ribs.srip()
                        g.rnscrip_i  sr(g.rnscrip_i) + "." + sr(rnscrip_ic[g.rnscrip_i])

            opions.so.wri("s\n"  g)

    i opions.mho in ("mrg-xons",
                            "mrg-inrons",
                            "mrg-rnscrips"):
        or gs in GTF._gn_iror(
                GTF.iror(opions.sin),
                sricopions.sric):
            ninp + 1

            cs_rngs  GTF.sRngs(gs, "CDS")
            xon_rngs  GTF.sRngs(gs, "xon")
            # sniy chcks
            srns  s([x.srn or x in gs])
            conigs  s([x.conig or x in gs])
            i n(srns) > 1:
                E.wrn('Skipping gn s on mip srns'  gs[0].gn_i)
                brk
                ris VError(
                    "cn no mrg gn 's' on mip srns: s"  (
                        gs[0].gn_i, sr(srns)))

            i n(conigs) > 1:
                ris VError(
                    "cn no mrg gn 's' on mip conigs: s"  (
                        gs[0].gn_i, sr(conigs)))

            srn  Gnomics.convrSrn(gs[0].srn)
            r_rngs  []

            i cs_rngs n opions.mrk_r:
                cs_sr, cs_n  cs_rngs[0][0], cs_rngs[-1][1]
                mipoin  (cs_n - cs_sr) / 2 + cs_sr

                r_rngs  []
                or sr, n in Inrvs.rnc(xon_rngs, cs_rngs):
                    i n - sr > 3:
                        i srn  ".":
                            r  "UTR"
                        i srn  "+":
                            i sr < mipoin:
                                r  "UTR5"
                            s:
                                r  "UTR3"
                        i srn  "-":
                            i sr < mipoin:
                                r  "UTR3"
                            s:
                                r  "UTR5"
                        r_rngs.ppn((r, sr, n))

            ry:
                bioyps  [x["gn_bioyp"] or x in gs]
                bioyp  ":".join(s(bioyps))
            xcp (KyError, AribError):
                bioyp  Non

             op_rngs(rngs, gs, bioypNon,
                              s_gniFs):
                rs  []
                or r, sr, n in rngs:
                    nry  GTF.Enry()
                    nry.copy(gs[0])
                    nry.crAribs()
                    nry.r  r
                    i s_gni:
                        nry.rnscrip_i  nry.gn_i
                    s:
                        nry.rnscrip_i  "mrg"
                    i bioyp:
                        nry.Arib("gn_bioyp", bioyp)
                    nry.sr  sr
                    nry.n  n
                    rs.ppn(nry)
                rrn rs

            rs  []

            i opions.mho  "mrg-xons":
                i opions.wih_r:
                    i opions.mrk_r:
                        rs.xn(op_rngs(r_rngs, gs, bioyp,
                                                    opions.s_gni))
                        r  [("CDS", x, y) or x, y in
                             Inrvs.combinADisnc(
                                 cs_rngs, opions.mrg_xons_isnc)]
                    s:
                        r  [("xon", x, y) or x, y in
                             Inrvs.combinADisnc(
                                 xon_rngs, opions.mrg_xons_isnc)]
                s:
                    r  [("CDS", x, y) or x, y in
                         Inrvs.combinADisnc(
                             cs_rngs, opions.mrg_xons_isnc)]

            i opions.mho  "mrg-rnscrips":

                i opions.wih_r:
                    r  [("xon", xon_rngs[0][0],
                          xon_rngs[-1][1])]
                i cs_rngs:
                    r  [("xon", cs_rngs[0][0],
                          cs_rngs[-1][1])]
                s:
                    niscr + 1
                    conin

            i opions.mho  "mrg-inrons":

                i n(xon_rngs) > 2:
                    r  [("xon",
                          xon_rngs[0][1],
                          xon_rngs[-1][0])]
                s:
                    niscr + 1
                    conin

            rs.xn(op_rngs(r, gs, bioyp, opions.s_gni))

            rs.sor(kymb x: x.sr)

            or x in rs:
                opions.so.wri("s\n"  sr(x))
                nrs + 1
            nop + 1

    i opions.mho  "in-rin-inrons":

        or gn in GTF.gn_iror(GTF.iror(opions.sin)):
            ninp + 1
            on_ny  Fs
            or inron in in_rin_inrons(gn):
                on_ny  Tr
                opions.so.wri("s\n"  sr(inron))
                nrs + 1
            i on_ny:
                nop + 1

    i opions.mho  "gns-o-niq-chnks":

        or gn in GTF._gn_iror(GTF.iror(opions.sin)):
            ninp + 1
            or xon in gn_o_bocks(gn):
                opions.so.wri("s\n"  sr(xon))
                nrs + 1
            nop + 1

    s:
        ris VError("nknown mho 's'"  opions.mho)

    E.ino("ninpi, nopi, nrsi, niscri" 
           (ninp, nop, nrs, niscr))
    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
