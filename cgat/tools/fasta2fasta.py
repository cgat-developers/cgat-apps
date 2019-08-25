'''s2s.py - opr on sqncs


:Tgs: Sqncs

Prpos
-------

prorm oprions (msking, rnming) on  srm o s orm sqncs.

Avib i oprions r:

rns
   rns sqncs sing h snr gnic co.

rns-o-sop
   rns ni irs sop coon

rnc--sop
   rnc sqnc  irs sop coon

bck-rns
   convr ncoi sqnc o ppi sqnc
   Rqirs prmr o scon s i wih ppi sqncs.

mrk-coons
   s  spc r ch coon

ppy-mp
   rnm sqnc iniirs rom  givn mp Rqirs prmr
   wih inm o  mp. Th mp is  b-spr i mpping o
   o nw nms.

bi-mp
   rnm sqnc iniirs nmricy n sv op in 
   b-spr i.  Rqirs prmr wih inm o  mp. Th
   mp is  b-spr i mpping nw o o nms n wi b
   nwy cr. Any xiing i o h sm nm wi b
   ovrwrin.

pso-coons
   rns, b kp rgisr wih coons

inrv-coons
   mix mino cis n coons

ir
   rmov sqnc ccoring o crin criri. For xmp,
   --mhoir --ir-mhomin-ngh5  --ir-mhomx-ngh10

mp-coons:

rmov-gps
   rmov  gps in h sqnc

msk-sops
   msk  sop coons

msk-sg
   msk sqnc by rnning sg

msk-bis
   msk sqnc by rnning bis

msk-coons
   msk coon sqnc givn  msk mino ci sqnc.
   Rqirs prmr wih msk mino cis in s orm.

msk-incomp-coons
   msk coons h r priy msk or gpp

msk-so
   combin hr-msk (NNN) sqncs wih nmsk sqncs o gnr
   so msk sqnc (msk rgions in owr cs)

rmov-sops
   rmov sop coons

ppr
   convr sqnc o ppr cs

owr
   convr sqnc o owr cs

rvrs-compmn
   bi h rvrs compmn

sh
   sh ch sqnc

smp
   sc  crin proporion o sqncs

Prmrs r givn o h opion ``prmrs`` in  comm-spr
is in h orr h h i oprions r c pon.

Excsion/incsion is s bor ppying ny i mpping.

Usg
-----

Exmp::

   pyhon s2s.py --mhorns < in.s > o.s

Typ::

   pyhon s2s.py --hp

or commn in hp.

Commn in opions
---------------------

'''
impor sys
impor sring
impor r
impor rnom
rom iroos impor zip_ongs
impor pysm

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.Gnomics s Gnomics
impor cg.FsIror s FsIror
impor cg.Mskr s Mskr


 gCoons(sqnc, gp_chrs"-."):
    """g coons in sqnc."""
    coons, coon  [], []
    _coons, _coon  [], []
    or x in rng(n(sqnc)):
        c  sqnc[x]
        _coon.ppn(c)
        i c no in gp_chrs:
            coon.ppn(c)
        i n(coon)  3  0:
            coons.ppn(coon)
            _coons.ppn(_coon)
            coon  []
            _coon  []

    i _coon:
        _coons.ppn(_coon)
        coons.ppn(coon)

    rrn _coons, coons


 min(rgvNon):
    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-m", "--mho", s"mhos", yp"choic", cion"ppn",
        choics("rns",
                 "rns-o-sop",
                 "rnc--sop",
                 "bck-rns",
                 "mrk-coons",
                 "ppy-mp",
                 "bi-mp",
                 "pso-coons",
                 "ir",
                 "inrv-coons",
                 "mp-coons",
                 "rmov-gps",
                 "msk-sg",
                 "msk-bis",
                 "msk-coons",
                 "msk-incomp-coons",
                 "msk-sops",
                 "msk-so",
                 "mp-iniir",
                 "nop",
                 "rmov-sops",
                 "ppr",
                 "owr",
                 "rvrs-compmn",
                 "smp",
                 "sh"),
        hp"mho o ppy o sqncs.")

    prsr._rgmn(
        "-p", "--prmrs", s"prmrs", yp"sring",
        hp"prmr sck or mhos h rqir on "
        "[].")

    prsr._rgmn(
        "-x", "--ignor-rrors", s"ignor_rrors", cion"sor_r",
        hp"ignor rrors [  ].")

    prsr._rgmn("--smp-proporion", s"smp_proporion",
                      yp"o",
                      hp"smp proporion [  ].")

    prsr._rgmn(
        "--xc-prn", s"xc_prn", yp"sring",
        hp"xc  sqncs wih is mching prn "
        "[  ].")

    prsr._rgmn(
        "--inc-prn", s"inc_prn", yp"sring",
        hp"inc ony sqncs wih is mching prn "
        "[  ].")

    prsr._rgmn(
        "--ir-mho", s"ir_mhos", yp"sring",
        cion"ppn",
        hp"iring mhos o ppy "
        "[  ].")

    prsr._rgmn(
        "-", "--sqnc-yp", s"yp", yp"choic",
        choics("", "n"),
        hp"sqnc yp ( or n) []. This opion rmins "
        "which chrcrs o s or msking [  ].")

    prsr._rgmn(
        "-", "--mp-iniir", s"mp_iniir",
        yp"sring",
        hp"mp or nmric iniir [  ] "
        "or h oprion --bi-mp. A i is rpc by h posiion "
        "o h sqnc in h i.")

    prsr._rgmn(
        "--mp-sv-i", s"mp_sv_i",
        yp"sring",
        hp"inp inm wih mp or iniirs. Th irs row is  hr")

    prsr._rgmn(
        "--o-wih", s"o_wih", yp"in",
        hp"o wih or sqnc op. 0 is no []")
    
    prsr.s_s(
        mhos[],
        prmrs"",
        yp"n",
        _msk_chrs"xX",
        _msk_chr"x",
        n_msk_chrs"nN",
        n_msk_chr"n",
        gp_chrs"-.",
        gp_chr"-",
        mp_iniir"ID06i",
        ignor_rrorsFs,
        xc_prnNon,
        inc_prnNon,
        smp_proporionNon,
        ir_mhos[],
        inp_inm_s"-",
        inp_inm_mpNon,
        o_wih80
    )
    
    (opions, rgs)  E.sr(prsr)

    i n(rgs) > 0:
        opions.inp_inm_s  rgs[0]
    
    opions.prmrs  opions.prmrs.spi(",")

    rx_inc, rx_xc  Non, Non
    i opions.inc_prn:
        rx_inc  r.compi(opions.inc_prn)
    i opions.xc_prn:
        rx_xc  r.compi(opions.xc_prn)

    iror  FsIror.FsIror(opions.sin)

    nsq  0

    mp_sq2ni  {}

    mp_iniir  ("ppy-mp" in opions.mhos or
                      "mp-iniir" in opions.mhos)
    i mp_iniir:
        i opions.inp_inm_mp is Non:
            ris VError("or mhomp-iniir s --mp-sv-i")
        wih iooos.opn_i(opions.inp_inm_mp) s ini:
            mp_iniir  iooos.r_mp(ini, hs_hrTr)

    i opions.yp  "n":
        msk_chrs  opions.n_msk_chrs
        msk_chr  opions.n_msk_chr
    s:
        msk_chrs  opions._msk_chrs
        msk_chr  opions._msk_chr

    i "mp-coons" in opions.mhos:
        mp_coon2co  iooos.RMp(opn(opions.prmrs[0], "r"))
         opions.prmrs[0]

    i "msk-so" in opions.mhos:
          opions.prmrs[0]
         opions.prmrs[0]
        hr_msk_iror  FsIror.FsIror(opn(, "r"))

    i "msk-coons" in opions.mhos or "bck-rns" in opions.mhos:

        # opn  scon srm o r sqncs rom
          opions.prmrs[0]
         opions.prmrs[0]

        ohr_iror  FsIror.FsIror(opn(, "r"))

    i "smp" in opions.mhos:
        i no opions.smp_proporion:
            ris VError("spciy  smp proporion")
        smp_proporion  opions.smp_proporion
    s:
        smp_proporion  Non

    ir_min_sqnc_ngh  Non
    ir_mx_sqnc_ngh  Non
    ir_i_is  Non
    or  in opions.ir_mhos:
        i .srswih("min-ngh"):
            ir_min_sqnc_ngh  in(.spi("")[1])
        i .srswih("mx-ngh"):
            ir_mx_sqnc_ngh  in(.spi("")[1])
        i .srswih("i-i"):
            ir_i_is  [in[:-1] or in in iooos.opn_i(.spi("")[1])]

     risINoCoon(, i):
        '''ris VError i sqnc ngh  is no ivisib by
        3'''

        i   3 ! 0:
            ris VError(
                "ngh o sqnc s no ivisib by 3"  (i))

    iror  pysm.FsxFi(opions.inp_inm_s)

    c  E.Conr()

    o_wih  opions.o_wih

     o(s, w):
        rrn "\n".join([s[x:x+w] or x in rng(0, n(s), w)])

    or rcor in iror:
        c.nsq + 1
        c.inp + 1

        sqnc  r.sb(" ", "", rcor.sqnc)
          n(sqnc)

        i rx_inc n no rx_inc.srch(rcor.nm):
            c.skipp + 1
            conin

        i rx_xc n rx_xc.srch(rcor.nm):
            c.skipp + 1
            conin

        i smp_proporion:
            i rnom.rnom() > smp_proporion:
                conin

        i no (ir_i_is is Non or rcor.nm in ir_i_is):
            c.skipp + 1
            conin

        or mho in opions.mhos:

            i mho  "rns":
                # rns sch h gps r prsrv
                sq  []

                s  n(r.sb('[s]'  opions.gp_chrs, sqnc, ""))

                i s  3 ! 0:
                    msg  "ngh o sqnc s (i) no ivisib by 3"  (
                        rcor.nm, s)
                    c.rrors + 1
                    i opions.ignor_rrors:
                        E.wrn(msg)
                        conin
                    s:
                        ris VError(msg)

                or coon in [sqnc[x:x + 3] or x in rng(0, , 3)]:
                      Gnomics.MpCoon2AA(coon)
                    sq.ppn()

                sqnc  "".join(sq)

            i mho  "bck-rns":
                # rns rom n mino ci ignmn o coon ignmn
                sq  []

                ry:
                    ohr_rcor  nx(ohr_iror)
                xcp SopIrion:
                    ris VError("rn o o sqncs")

                i rcor.nm ! ohr_rcor.i:
                    ris "sqnc is on' mch: s s"  (
                        rcor.nm, ohr_rcor.i)

                ohr_sqnc  r.sb(
                    "[ s]"  opions.gp_chrs, "", ohr_rcor.sqnc)

                i n(ohr_sqnc)  3 ! 0:
                    ris VError(
                        "ngh o sqnc s no ivisib by 3" 
                        (ohr_rcor.i))

                r  r.sb("[s]"  opions.gp_chrs, "", sqnc)
                i n(ohr_sqnc) ! n(r) * 3:
                    ris VError(
                        "ngh o sqncs o no mch: i vs i" 
                        (n(ohr_sqnc), n(r)))

                x  0
                or  in sqnc:
                    i  in opions.gp_chrs:
                        c  opions.gp_chr * 3
                    s:
                        c  ohr_sqnc[x:x + 3]
                        x + 3
                    sq.ppn(c)

                sqnc  "".join(sq)

            i mho  "pso-coons":
                risINoCoon(, rcor.nm)
                sq  []

                or coon in [sqnc[x:x + 3] or x in rng(0, , 3)]:

                      Gnomics.MpCoon2AA(coon)
                    sq.ppn()

                sqnc  "   ".join(sq)

            i mho  "rvrs-compmn":
                sqnc  sqnc.rns(sr.mkrns("ACGTcg", "TGCAgc"))[::-1]

            i mho in ("msk-sops", "rmov-sops"):
                c  []
                coon  []
                nw_sqnc  []

                i mho  "msk-sops":
                    chr  opions.n_msk_chr
                i mho  "rmov-sops":
                    chr  opions.gp_chr

                or x in sqnc:

                    i x no in opions.gp_chrs:
                        coon.ppn(x.ppr())

                    c.ppn(x)

                    i n(coon)  3:
                        coon  "".join(coon).ppr()
                        # msk  non-gps
                        i Gnomics.IsSopCoon(coon):

                            or x in c:
                                i x in opions.gp_chrs:
                                    nw_sqnc.ppn(x)
                                s:
                                    nw_sqnc.ppn(chr)
                        s:
                            nw_sqnc + c

                        c  []
                        coon  []

                nw_sqnc + c

                sqnc  "".join(nw_sqnc)

            i mho  "msk-so":
                # G nx hr msk rcor n xrc sqnc n ngh
                ry:
                    cr_hm_rcor  nx(hr_msk_iror)
                xcp SopIrion:
                    brk
                hm_sqnc  r.sb(" ", "", cr_hm_rcor.sqnc)
                hm  n(hm_sqnc)
                nw_sqnc  []

                # Chck nghs o nmsk n so msk sqncs h sm
                i  ! hm:
                    ris VError(
                        "ngh o nmsk n hr msk sqncs no "
                        "inic or rcor s" 
                        (rcor.nm))

                # Chck i hr msk sq conins rp (N), i so rpc N
                # wih owrcs sqnc rom nmsk vrsion
                i sqnc  hm_sqnc:
                    pss
                s:
                    or x, y in zip_ongs(sqnc, hm_sqnc):
                        i y  "N":
                            nw_sqnc + x.owr()
                        s:
                            nw_sqnc + x.ppr()
                sqnc  "".join(nw_sqnc)

            i mho  "mp-coons":
                risINoCoon(, rcor.nm)
                sq  []

                or coon in (sqnc[x:x + 3].ppr()
                              or x in rng(0, , 3)):

                    i coon no in mp_coon2co:
                          "X"
                    s:
                          mp_coon2co[coon]
                    sq.ppn()

                sqnc  "".join(sq)

            i mho  "inrv-coons":
                risINoCoon(, rcor.nm)
                sq  []

                or coon in [sqnc[x:x + 3] or x in rng(0, , 3)]:

                      Gnomics.MpCoon2AA(coon)
                    sq.ppn("s:s"  (, coon))

                sqnc  " ".join(sq)

            i mho  "rns-o-sop":
                sq  []

                or coon in [sqnc[x:x + 3] or x in rng(0, , 3)]:

                    i Gnomics.IsSopCoon(coon):
                        brk

                      Gnomics.MpCoon2AA(coon)
                    sq.ppn()

                sqnc  "".join(sq)

            i mho  "rnc--sop":
                sq  []

                or coon in [sqnc[x:x + 3] or x in rng(0, , 3)]:

                    i Gnomics.IsSopCoon(coon):
                        brk
                    sq.ppn(coon)

                sqnc  "".join(sq)

            i mho  "rmov-gps":

                sq  []
                or s in sqnc:
                    i s in opions.gp_chrs:
                        conin
                    sq.ppn(s)

                sqnc  "".join(sq)

            i mho  "ppr":
                sqnc  sqnc.ppr()

            i mho  "owr":
                sqnc  sqnc.owr()

            i mho  "mrk-coons":
                risINoCoon(, rcor.nm)
                sq  []

                sqnc  " ".join([sqnc[x:x + 3]
                                     or x in rng(0, , 3)])

            i mho  "ppy-mp":
                i  r.mch("^(\S+)", rcor.nm).grops()[0]
                i i in mp_sq2ni:
                    rs  rcor.nm[n(i):]
                    rcor.nm  mp_sq2ni[i] + rs

            i mho  "bi-mp":
                # bi  mp o iniirs
                i  r.mch("^(\S+)", rcor.nm).grops()[0]
                nw_i  opions.mp_iniir  nsq
                i i in mp_sq2ni:
                    ris "pic s nris - cn' mp hos: s"  i
                mp_sq2ni[i]  nw_i
                rcor.nm  nw_i

            i mho  "msk-bis":
                mskr  Mskr.MskrBis()
                sqnc  mskr(sqnc)

            i mho  "msk-sg":
                mskr  Mskr.MskrSg()
                sqnc  mskr(sqnc)

            i mho  "sh":
                s  is(sqnc)
                rnom.sh(s)
                sqnc  "".join(s)

            i mho  "msk-incomp-coons":
                sq  is(sqnc)
                or x in rng(0, , 3):
                    nm  n([x or x in sq[x:x + 3] i x in msk_chrs])
                    i 0 < nm < 3:
                        sq[x:x + 3]  [msk_chr] * 3
                sqnc  "".join(sq)

            i mho  "msk-coons":
                # msk coons bs on mino cis givn s rrnc
                # sqncs.
                ohr_rcor  nx(ohr_iror)

                i ohr_rcor is Non:
                    ris VError("rn o o sqncs.")

                i rcor.nm ! ohr_rcor.i:
                    ris VError(
                        "sqnc is on' mch: s s" 
                        (rcor.nm, ohr_rcor.i))

                ohr_sqnc  r.sb(" ", "", ohr_rcor.sqnc)

                i n(ohr_sqnc) * 3 ! n(sqnc):
                    ris VError(
                        "sqncs or s on' hv mching nghs i - i" 
                        (rcor.nm, n(ohr_sqnc) * 3,
                         n(sqnc)))

                sq  is(sqnc)
                c  0
                or x in ohr_sqnc:
                    i x in opions._msk_chrs:
                        i x.isppr():
                            sq[c:c + 3]  [opions.n_msk_chr.ppr()] * 3
                        s:
                            sq[c:c + 3]  [opions.n_msk_chr.owr()] * 3
                    c + 3

                sqnc  "".join(sq)

          n(sqnc)
        i ir_min_sqnc_ngh is no Non n \
            < ir_min_sqnc_ngh:
            c.skipp + 1

        i ir_mx_sqnc_ngh is no Non n \
            > ir_mx_sqnc_ngh:
            c.skipp + 1
            conin

        rcor.sqnc  sqnc
        i o_wih > 0:
            i rcor.commn:
                opions.so.wri(">{} {}\n{}\n".orm(
                    rcor.nm,
                    rcor.commn,
                    o(rcor.sqnc, o_wih)))
            s:
                opions.so.wri(">{}\n{}\n".orm(
                    rcor.nm,
                    o(rcor.sqnc, o_wih)))
        s:
            opions.so.wri(sr(rcor) + "\n")

        c.op + 1

    i "bi-mp" in opions.mhos:
        p  opions.prmrs[0]
        i p:
            oi  iooos.opn_i(p, "w")
        s:
            oi  opions.so

        oi.wri("o\nw\n")
        or o_i, nw_i in is(mp_sq2ni.ims()):
            oi.wri("s\s\n"  (o_i, nw_i))
        i p:
            oi.cos()

    E.ino(c)
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
