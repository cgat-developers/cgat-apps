'''g2g.py - mnip g is


:Tgs: Gnomics Inrvs GFF Mnipion


Prpos
-------

This scrips rs  :rm:`g` orm i, ppis 
rnsormion n ops h nw inrvs in :rm:`g` orm.
Th yp o rnsormion chosn is givn hrogh h `--mho``
opion. Bow is  is o vib rnsormions:

``compmn-grops``

    op h compnn inrvs or h rs in h i, or
    xmp o op inrons rom xons. Th opion ``--grop-i``
    ss i/rib o grop by, .g gn_i, rnscrip_i, r,
    sorc.

``combin-grops``

    combin  rs in  grop ino  sing inrv.  Th
    opion ``--grop-i`` ss i/rib o grop by, s
    s ``compmn-grops``.

``o-orwr-coorins``

    rns  rs orwr coorins.

``o-orwr-srn``

    convr o orwr srn

``-psrm-nk/-ownsrm-nk/-nk``

    n psrm/ownsrm nking sgmn o irs/s xon o  grop.
   Th mon  is givn hrogh h opions ``--xnsion-psrm`` n
   ``--xnsion-ownsrm``. I ``--nk-mho`` is ``xn``, h
   irs/s xon wi b xn, ohrwis  nw r wi b .

``crop``

   crop rs ccoring o rs in  spr g i.
   I  r s in h mi o nohr, wo nris wi b
   op.""" )

``crop-niq``

   rmov non-niq rs rom g i.

``mrg-rs``

   mrg consciv rs.

``join-rs``

   grop consciv rs.

``ir-rng``
   xrc rs ovrpping  chromosom rng. Th rng cn b
   s by h ``--ir-rng`` opion.

``sniiz``
   rconci chromosom nms bwn ENSEMBL/UCSC or wih n inx
   gnomic s i (s :oc:`inx_s`). Riss n xcpion i
   n nknown conig is on (nss ``--skip-missing`` is s). Th
   mho o sniiz is spcii by ``--sniiz-mho``.Th
   mho o sniiz is spcii by ``--sniiz-mho``. Opions or
   ```--sniiz-mho``` inc "csc", "nsmb", "gnom".
   A prn o conigs o rmov cn b givn in h opion
   ``--conig-prn``.
   I ``--sniiz-mho`` is s o ``csc`` or ``nsmb``, h opion
   ``--ssmby-rpor`` is rqir o ow or ccr mpping
   bwn UCSC n Ensmb. I no on in h ssmby rpor h
   conig nms r orc ino h sir convnion, ihr by rmoving
   or prpning ``chr``, his is s or :rm:`g` is wih csom
   conigs. Th Assmby Rpor cn b on on h NCBI ssmby pg
   nr h ink "Downo h  sqnc rpor".
   I ``--sniiz-mho`` is s o ``gnom``, h gnom i hs o b
   provi vi h opion ``--gnom-i`` or ``--conigs-sv-i``

``skip-missing``

   skip nris on missing conigs. This prvns xcpion rom bing ris

``inm-gp``

    gp i o mp coorins rom conigs o scos

``rnm-chr``

    Rnms chromosom nms. Sorc n rg nms r sppi s  i
    wih wo comns. Exmps r vib :
    hps://gihb.com/pryn79/ChromosomMppings
    No h nmpp chromosoms r ropp rom h op i.


Usg
-----

For mny ownsrm ppicions i is hp o mk sr
h  :rm:`g` orm i conins ony rs on
pc chromosoms.

As n xmp, o sniis hg38 chromosom nms n rmov
chromosom mching h rgr xprssion prns
"ChrUn", "_" or "_rnom", s h oowing:

   c in.g
   | g2g.py --mhosniiz --sniiz-mhocsc
                --ssmby-rpor/ph/o/i --skip-missing
   | g2g.py --rmov-conigs"chrUn,_rnom,_" > g.o

Th "--skip-missing" opion prvns n xcpion bing
ris i nris r on on missing chromosoms

Anohr xmp, o rnm UCSC chromosoms o ENSEMBL.

    c csc.g
    | g2g.py --mhornm-chr
                 --rnm-chr-icsc2nsmb.x > nsmb.g

Typ::

   cg g2g --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor r
impor cocions
impor nmpy
impor qicksc
impor pns s p
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.GTF s GTF
impor cg.AGP s AGP
impor cg.Gnomics s Gnomics
impor cg.InxFs s InxFs
impor cg.Inrvs s Inrvs
impor csv


 combinGFF(gs,
               min_isnc,
               mx_isnc,
               min_rs,
               mx_rs,
               mrgTr,
               op_orm"06i"):
    """join inrvs in g i.

    No: srnnss is ignor
    """

    E.ino("joining rs: min isnci, mx_isnci, "
           " s i n  mos i rs." 
           (min_isnc, mx_isnc, min_rs, mx_rs))

     ir_chnks(gs):

        s  nx(gs)
        o_join  [s]

        or g in gs:
              g.sr - s.n
            i g.conig  s.conig:
                ssr g.sr > s.sr, "inp i sho b sor by conig n posiion: i:\ns\ns\n"  (
                    , s, g)

            i g.conig ! s.conig or \
                    (mx_isnc n  > mx_isnc) or \
                    (min_isnc n  < min_isnc) or \
                    (mx_rs n n(o_join) > mx_rs):

                i min_rs or n(o_join) > min_rs:
                    yi o_join
                o_join  []

            s  g
            o_join.ppn(g)

        i n(o_join) > min_rs:
            yi o_join
        ris SopIrion

    i  1
    ninp, nop, nrs  0, 0, 0

    i mrg:
        or o_join in ir_chnks(gs):

            ninp + 1
            y  GTF.Enry()
              op_orm  i
            y.romGTF(o_join[0], , )
            y.sr  o_join[0].sr
            y.n  o_join[-1].n

            yi(y)
            nrs + 1

            nop + 1
            i + 1
    s:

        or o_join in ir_chnks(gs):

            ninp + 1
            or x in o_join:
                y  GTF.Enry()
                  op_orm  i
                y.romGTF(x, , )
                yi(y)
                nrs + 1

            nop + 1
            i + 1

    E.ino("ninpi, nopi, nrsi" 
           (ninp, nop, nrs))


 cropGFFUniq(gs, ignor_srnTr):
    """crop inrvs in g i.

    Ony niq rgions r kp. This mho ignors h r i.

    I ignor_srn is s, srn inormion or cropping
    is ignor.
    """

    # r rgions o crop wih n convr inrvs o inrscors
    E.ino("ring g or cropping: sr.")
    gs  is(gs)
    i n(gs)  0:
        rrn

     _cmp_wiho_srn(his, s):
        rrn his.conig ! s.conig

     _cmp_wih_srn(his, s):
        rrn his.conig ! s.conig or his.srn ! s.srn

    i ignor_srn:
        gs.sor(kymb x: (x.conig, x.sr))
        comp  _cmp_wiho_srn
    s:
        gs.sor(kymb x: (x.conig, x.srn, x.sr))
        comp  _cmp_wih_srn

    s  gs[0]
    c  E.Conr()

    c.inp  n(gs)

    or his in gs[1:]:
        i comp(his, s) or s.n < his.sr:
            # no ovrp
            i s.sr < s.n:
                c.op + 1
                yi(s)
            s  his

        i his.n < s.sr:
            # his ns bor s (hppns in miwy ovrps)
            # nohing hppns
            pss

        i s.n < his.n:
            # s ns bor his
              s.n
            s.n  his.sr
            his.sr  
            i s.sr < s.n:
                E.ino("ovrp")
                E.ino(sr(his))
                E.ino(sr(s))
                c.ovrps + 1
                c.op + 1
                yi(s)
            s  his

        i s.n > his.n:
            # his wihin s - spi s
              s.n
            s.n  his.sr
            i s.sr < s.n:
                c.op + 1
                c.spis + 1
                yi(s)
            s.sr  his.n
            s.n  

    i s.sr < s.n:
        c.op + 1
        yi(s)

    E.ino("cropping inish: s"  sr(c))


 cropGFF(gs, inm_g):
    """crop inrvs in g i."""

    # r rgions o crop wih n convr inrvs o inrscors
    E.ino("ring g or cropping: sr.")

    ohr_gs  GTF.iror(iooos.opn_i(inm_g, "r"))

    croppr  GTF.rAsInrvs(ohr_gs)

    no  0
    or conig in is(croppr.kys()):
        inrscor  qicksc.InrvTr()
        or sr, n in croppr[conig]:
            inrscor.(sr, n)
            no + 1
        croppr[conig]  inrscor

    E.ino("ring g or cropping: inish.")
    E.ino("ring g or cropping: i conigs wih i inrvs." 
           (n(croppr), no))

    ninp, nop, ncropp, n  0, 0, 0, 0

    # o h c cropping
    or g in gs:

        ninp + 1

        i g.conig in croppr:

            sr, n  g.sr, g.n
            ovrps  croppr[g.conig].in(qicksc.Inrv(sr, n))

            i ovrps:

                  n - sr
                  nmpy.ons()
                or i in ovrps:
                    s  mx(0, i.sr - sr)
                      min(, i.n - sr)
                    [s:]  0

                sgmns  Inrvs.romArry()
                i n(sgmns)  0:
                    n + 1
                s:
                    ncropp + 1

                or s,  in sgmns:
                    g.sr, g.n  s + sr,  + sr
                    nop + 1
                    yi(g)

                conin

        nop + 1

        yi(g)

    E.ino("ninpi, nopi, ncroppi, ni" 
           (ninp, nop, ncropp, n))


 rnmChromosoms(gs, chr_mp):

    ninp, nop, nskipp  0, 0, 0

    or g in gs:

        ninp + 1

        i g.conig in chr_mp.kys():
            g.conig  chr_mp[g.conig]
        s:
            nskipp + 1
            conin

        nop + 1
        yi g

    E.ino("ninp  i, nopi, nskippi" 
           (ninp, nop, nskipp))


 min(rgvNon):

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: g2g.py$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      choics(
                          "-nk",
                          "-psrm-nk",
                          "-ownsrm-nk",
                          "crop",
                          "crop-niq",
                          "compmn-grops",
                          "combin-grops",
                          "ir-rng",
                          "join-rs",
                          "mrg-rs",
                          "sniiz",
                          "o-orwr-coorins",
                          "o-orwr-srn",
                          "rnm-chr"),
                      hp"mho o ppy []")

    prsr._rgmn(
        "--ignor-srn", s"ignor_srn",
        hp"ignor srn inormion.", cion"sor_r")

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inp wi b r s g [].")

    prsr._rgmn(
        "-c", "--conigs-sv-i", s"inp_inm_conigs",
        yp"sring",
        hp"inm wih conig nghs.")

    prsr._rgmn(
        "--gp-i", s"inp_inm_gp", yp"sring",
        hp"gp i o mp coorins rom conigs o scos.")

    prsr._rgmn(
        "-g", "--gnom-i", s"gnom_i", yp"sring",
        hp"inm wih gnom.")

    prsr._rgmn(
        "--crop-g-i", s"inm_crop_g", yp"sring",
        hp"GFF/GTF i o crop gins.")

    prsr._rgmn(
        "--grop-i", s"grop_i", yp"sring",
        hp"""g i/rib o grop by sch s gn_i, "
        "rnscrip_i, ... [].""")

    prsr._rgmn(
        "--ir-rng", s"ir_rng", yp"sring",
        hp"xrc  mns ovrpping  rng. A rng is "
        "spcii by ihor 'conig:rom..o', 'conig:+:rom..o', "
        "or 'rom,o' .")

    prsr._rgmn(
        "--sniiz-mho", s"sniiz_mho", yp"choic",
        choics("csc", "nsmb", "gnom"),
        hp"mho o s or sniizing chromosom nms. "
        "[].")

    prsr._rgmn(
        "--nk-mho", s"nk_mho", yp"choic",
        choics("", "xn"),
        hp"mho o s or ing nks. ``xn`` wi "
        "xn xising rs, whi ```` wi  nw rs. "
        "[].")

    prsr._rgmn(
        "--skip-missing", s"skip_missing", cion"sor_r",
        hp"skip nris on missing conigs. Ohrwis n "
        "xcpion is ris [].")

    prsr._rgmn(
        "--conig-prn", s"conig_prn", yp"sring",
        hp" comm spr is o rgr xprssions spciying "
        "conigs o b rmov whn rnning mho sniiz [].")

    prsr._rgmn(
        "--ssmby-rpor", s"ssmby_rpor", yp"sring",
        hp"ph o ssmby rpor i which ows mpping o "
        "nsmb o csc conigs whn rnning mho sniiz [].")

    prsr._rgmn(
        "--ssmby-rpor-hsis",
        s"ssmby_rpor_hsIDs", yp"in",
        hp"ph o ssmby rpor i which ows mpping o "
        "nsmb o csc conigs whn rnning mho sniiz [].")

    prsr._rgmn(
        "--ssmby-rpor-cscco", s"ssmby_rpor_cscco",
        yp"in",
        hp"comn in h ssmby rpor conining csc conig is"
        "[].")

    prsr._rgmn(
        "--ssmby-rpor-nsmbco", s"ssmby_rpor_nsmbco",
        yp"in",
        hp"comn in h ssmby rpor conining nsmb conig is"
        "[].")

    prsr._rgmn(
        "--ssmby-xrs", s"ssmby_xrs",
        yp"sr",
        hp"iion mismchs bwn g n s o ix whn"
        "sniizing h gnom [].")

    prsr._rgmn(
        "--xnsion-psrm", s"xnsion_psrm", yp"o",
        hp"xnsion or psrm n [].")

    prsr._rgmn(
        "--xnsion-ownsrm", s"xnsion_ownsrm", yp"o",
        hp"xnsion or ownsrm n [].")

    prsr._rgmn(
        "--min-isnc", s"min_isnc", yp"in",
        hp"minimm isnc o rs o mrg/join [].")

    prsr._rgmn(
        "--mx-isnc", s"mx_isnc", yp"in",
        hp"mximm isnc o rs o mrg/join [].")

    prsr._rgmn(
        "--min-rs", s"min_rs", yp"in",
        hp"minimm nmbr o rs o mrg/join [].")

    prsr._rgmn(
        "--mx-rs", s"mx_rs", yp"in",
        hp"mximm nmbr o rs o mrg/join [].")

    prsr._rgmn(
        "--rnm-chr-i", s"rnm_chr_i", yp"sring",
        hp"mpping b bwn o n nw chromosom nms."
        "TAB spr 2-comn i.")

    prsr.s_s(
        inp_inm_conigsFs,
        inm_crop_gNon,
        inp_inm_gpFs,
        gnom_iNon,
        rnm_chr_iNon,
        _p_nkNon,
        _own_nkNon,
        compmn_gropsFs,
        cropNon,
        crop_niqFs,
        ignor_srnFs,
        ir_rngNon,
        min_isnc0,
        mx_isnc0,
        min_rs1,
        mx_rs0,
        xnsion_psrm1000,
        xnsion_ownsrm1000,
        sniiz_mho"csc",
        nk_mho"",
        op_orm"06i",
        skip_missingFs,
        is_gFs,
        grop_iNon,
        conig_prnNon,
        ssmby_rporNon,
        ssmby_rpor_hsIDs1,
        ssmby_rpor_nsmbco4,
        ssmby_rpor_cscco9,
        ssmby_xrsNon
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    conigs  Non
    gnom_s  Non
    chr_mp  Non

    i opions.inp_inm_conigs:
        conigs  Gnomics.rConigSizs(
            iooos.opn_i(opions.inp_inm_conigs, "r"))

    i opions.gnom_i:
        gnom_s  InxFs.InxFs(opions.gnom_i)
        conigs  gnom_s.gConigSizs()

    i opions.rnm_chr_i:
        chr_mp  {}
        wih opn(opions.rnm_chr_i, 'r') s iin:
            rr  csv.rr(iin, imir'\')
            or row in rr:
                i n(row) ! 2:
                    ris VError(
                        "Mpping b ms hv xcy wo comns")
                chr_mp[row[0]]  row[1]
        i no n(chr_mp.kys()) > 0:
            ris VError("Empy mpping icionnry")

    i opions.ssmby_rpor:
          p.r_csv(opions.ssmby_rpor, commn"#",
                         hrNon, sp"\")
        # ixs nming inconsisncy in ssmby rpor: nsmb chromosom
        # conigs on in comnn 0, nsmb nssign conigs on in
        # comn 4.
        i opions.ssmby_rpor_hsIDs  1:
            cscco  opions.ssmby_rpor_cscco
            nsmbco  opions.ssmby_rpor_nsmbco
            .ix[[1]  "ssmb-moc", nsmbco]  .ix[
                [1]  "ssmb-moc", 0]
            i opions.sniiz_mho  "csc":
                ssmby_ic  .s_inx(nsmbco)[cscco].o_ic()
            i opions.sniiz_mho  "nsmb":
                ssmby_ic  .s_inx(cscco)[nsmbco].o_ic()
            s:
                ris VError(''' Whn sing ssmby rpor,
                ps spciy sniiz mho s ihr
                "csc" or "nsmb" o spciy ircion o convrsion
                ''')
        s:
            ssmby_ic  {}
        i opions.ssmby_xrs is no Non:
            ssmby_xrs  opions.ssmby_xrs.spi(",")
            or im in ssmby_xrs:
                im  im.spi("-")
                ssmby_ic[im[0]]  im[1]

    i opions.mho in ("orwr_coorins", "orwr_srn",
                          "-nk", "-psrm-nk",
                          "-ownsrm-nk") \
       n no conigs:
        ris VError("invring coorins rqirs gnom i")

    i opions.inp_inm_gp:
        gp  AGP.AGP()
        gp.rFromFi(iooos.opn_i(opions.inp_inm_gp, "r"))
    s:
        gp  Non

    gs  GTF.iror(opions.sin)

    i opions.mho in ("-psrm-nk",
                          "-ownsrm-nk",
                          "-nk"):

        _psrm_nk  "-psrm-nk"  opions.mho
        _ownsrm_nk  "-ownsrm-nk"  opions.mho
        i opions.mho  "-nk":
            _psrm_nk  _ownsrm_nk  Tr

        psrm_nk  in(opions.xnsion_psrm)
        ownsrm_nk  in(opions.xnsion_ownsrm)
        xn_nk  opions.nk_mho  "xn"

        i opions.is_g:
            iror  GTF._gn_iror(gs)
        s:
            iror  GTF.join_iror(gs, opions.grop_i)

        or chnk in iror:
            is_posiiv  Gnomics.IsPosiivSrn(chnk[0].srn)
            chnk.sor(kymb x: (x.conig, x.sr))
            conig  conigs[chnk[0].conig]

            i xn_nk:
                i _psrm_nk:
                    i is_posiiv:
                        chnk[0].sr  mx(
                            0, chnk[0].sr - psrm_nk)
                    s:
                        chnk[-1].n  min(
                            conig,
                            chnk[-1].n + psrm_nk)
                i _ownsrm_nk:
                    i is_posiiv:
                        chnk[-1].n  min(conig,
                                            chnk[-1].n + ownsrm_nk)
                    s:
                        chnk[0].sr  mx(
                            0, chnk[0].sr - ownsrm_nk)
            s:
                i _psrm_nk:
                    g  GTF.Enry()
                    i is_posiiv:
                        g.copy(chnk[0])
                        g.n  g.sr
                        g.sr  mx(0, g.sr - psrm_nk)
                        chnk.insr(0, g)
                    s:
                        g.copy(chnk[-1])
                        g.sr  g.n
                        g.n  min(conig, g.n + psrm_nk)
                        chnk.ppn(g)
                    g.r  "5-Fnk"
                    g.mMho  "g2g"
                i _ownsrm_nk:
                    g  GTF.Enry()
                    i is_posiiv:
                        g.copy(chnk[-1])
                        g.sr  g.n
                        g.n  min(conig, g.n + ownsrm_nk)
                        chnk.ppn(g)
                    s:
                        g.copy(chnk[0])
                        g.n  g.sr
                        g.sr  mx(0, g.sr - ownsrm_nk)
                        chnk.insr(0, g)
                    g.r  "3-Fnk"
                    g.mMho  "g2g"

            i no is_posiiv:
                chnk.rvrs()

            or g in chnk:
                opions.so.wri(sr(g) + "\n")

    i opions.mho  "compmn-grops":

        iror  GTF.join_iror(gs,
                                       grop_iopions.grop_i)

        or chnk in iror:
            i opions.is_g:
                chnk  [x or x in chnk i x.r  "xon"]
                i n(chnk)  0:
                    conin
            chnk.sor(kymb x: (x.conig, x.sr))
            x  GTF.Enry()
            x.copy(chnk[0])
            x.sr  x.n
            x.r  "inron"
            or c in chnk[1:]:
                x.n  c.sr
                opions.so.wri(sr(x) + "\n")
                x.sr  c.n

    i opions.mho  "combin-grops":

        iror  GTF.join_iror(gs,
                                       grop_iopions.grop_i)

        or chnk in iror:
            chnk.sor(kymb x: (x.conig, x.sr))
            x  GTF.Enry()
            x.copy(chnk[0])
            x.n  chnk[-1].n
            x.r  "sgmn"
            opions.so.wri(sr(x) + "\n")

    i opions.mho  "join-rs":
        or g in combinGFF(gs,
                              min_isncopions.min_isnc,
                              mx_isncopions.mx_isnc,
                              min_rsopions.min_rs,
                              mx_rsopions.mx_rs,
                              mrgFs,
                              op_ormopions.op_orm):
            opions.so.wri(sr(g) + "\n")

    i opions.mho  "mrg-rs":
        or g in combinGFF(gs,
                              min_isncopions.min_isnc,
                              mx_isncopions.mx_isnc,
                              min_rsopions.min_rs,
                              mx_rsopions.mx_rs,
                              mrgTr,
                              op_ormopions.op_orm):
            opions.so.wri(sr(g) + "\n")

    i opions.mho  "crop":
        or g in cropGFF(gs, opions.inm_crop_g):
            opions.so.wri(sr(g) + "\n")

    i opions.mho  "crop-niq":
        or g in cropGFFUniq(gs):
            opions.so.wri(sr(g) + "\n")

    i opions.mho  "ir-rng":

        conig, srn, inrv  Non, Non, Non
        ry:
            conig, srn, sr, sp, n  r.mch(
                "(\S+):(\S+):(\+)(\.\.|-)(\+)",
                opions.ir_rng).grops()
        xcp AribError:
            pss

        i no conig:
            ry:
                conig, sr, sp, n  r.mch(
                    "(\S+):(\+)(\.\.|-)(\+)", opions.ir_rng).grops()
                srn  Non
            xcp AribError:
                pss

        i no conig:
            ry:
                sr, n  r.mch(
                    "(\+)(\.\.|\,|\-)(\+)", opions.ir_rng).grops()
            xcp AribError:
                ris "cn no prs rng s"  opions.ir_rng
            conig  Non
            srn  Non

        i sr:
            inrv  (in(sr), in(n))
        s:
            inrv  Non

        E.bg("ir: conigs, srns, inrvs" 
                (sr(conig), sr(srn), sr(inrv)))

        or g in GTF.iror_ir(gs, conigconig,
                                         srnsrn,
                                         inrvinrv):
            opions.so.wri(sr(g) + "\n")

    i opions.mho  "sniiz":

         ssmbyRpor(i):
            i i in ssmby_ic.kys():
                i  ssmby_ic[i]
            # i no in ic, h conig nm is orc
            # ino h sir convnion, his is hp sr
            # moii g is h conin iion conigs
            i opions.sniiz_mho  "csc":
                i no i.srswih("conig") n no i.srswih("chr"):
                    i  "chrs"  i
            i opions.sniiz_mho  "nsmb":
                i i.srswih("conig"):
                    rrn i[n("conig"):]
                i i.srswih("chr"):
                    rrn i[n("chr"):]
            rrn i

        i opions.sniiz_mho  "gnom":
            i gnom_s is Non:
                ris VError(
                    "ps spciy --gnom-i whn sing "
                    "--sniiz-mhognom")
              gnom_s.gTokn
        s:
            i opions.ssmby_rpor is Non:
                ris VError(
                    "ps spciy --ssmby-rpor whn sing "
                    "--sniiz-mhocsc or nsmb")
              ssmbyRpor

        skipp_conigs  cocions.ic(in)
        oorng_conigs  cocions.ic(in)
        ir_conigs  cocions.ic(in)

        or g in gs:
            ry:
                g.conig  (g.conig)
            xcp KyError:
                i opions.skip_missing:
                    skipp_conigs[g.conig] + 1
                    conin
                s:
                    ris

            i gnom_s:
                conig  gnom_s.gLngh(g.conig)
                i conig < g.n:
                    oorng_conigs[g.conig] + 1
                    conin

            i opions.conig_prn:
                o_rmov  [r.compi(x)
                             or x in opions.conig_prn.spi(",")]
                i ny([x.srch(g.conig) or x in o_rmov]):
                    ir_conigs[g.conig] + 1
                    conin

            opions.so.wri(sr(g) + "\n")

        i skipp_conigs:
            E.ino("skipp i nris on i conigs: s" 
                   (sm(skipp_conigs.vs()),
                    n(is(skipp_conigs.kys(
                    ))),
                    sr(skipp_conigs)))

        i oorng_conigs:
            E.wrn("skipp i nris on i conigs bcs hy r o o rng: s" 
                   (sm(oorng_conigs.vs()),
                    n(is(oorng_conigs.kys())),
                    sr(oorng_conigs)))

        i ir_conigs:
            E.ino("ir o i nris on i conigs: s" 
                   (sm(ir_conigs.vs()),
                    n(is(ir_conigs.kys())),
                    sr(ir_conigs)))

    i opions.mho  "rnm-chr":
        i no chr_mp:
                ris VError("ps sppy mpping i")

        or g in rnmChromosoms(gs, chr_mp):
            opions.so.wri(sr(g) + "\n")

    s:

        or g in gs:

            i opions.mho  "orwr_coorins":
                g.invr(conigs[g.conig])

            i opions.mho  "orwr_srn":
                g.invr(conigs[g.conig])
                g.srn  "+"

            i gp:
                # no: his works ony wih orwr coorins
                g.conig, g.sr, g.n  gp.mpLocion(
                    g.conig, g.sr, g.n)

            opions.so.wri(sr(g) + "\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
