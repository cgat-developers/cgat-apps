'''s2sq.py - sim rs rom s


:Tgs: Sqncs

Prpos
-------

Sim imin sqnc rs rom  s i. Th nmbr o
rs pr nry is rnomy sc rom h rng givn. Th primry
s cs is xpc o b h gnrion o simion RNA-Sq rs

For RNA-Sq simions, h prmrn-rcion opion ows h sr o spciy
wh rcion o h rnscrips origin rom pr-mRNA. Th sr ms so
sppy  scon s in h sm orr or h pr-mRNA
(--ini-prmrn-s). Th simion ssms  pr-mRNA r  ngh
which is no iky o b h cs or r RNA-Sq.
No: This my  o mny mor rs which ign o h mRNA hn h
pprn gron rh con. I is hror bs o kp h pr-mRNA rcion
ow (rcommn 0.01).

--cons-mho rmin whhr h --min-cons n --mx--cons mhos
r o h nmbr o rs or copis pr nry. I copis, hn h nmbr
o rs/r pirs pr rnscrip wi b rmin by h copy nmbr, h
ngh o h rnscrip n h rgmn ngh. E.g or copy nmbr  2,
rnscrip ngh  1000 bp, r ngh  100 bp, pir n, 300bp insr
siz, h nmbr o rs wi b 2 * (1000 / (( 2 * 100) + 300))  4

Opions
-------

--op-pir-n
   gnr pir-n rs (s o sing n)

--cons-mho
  sim  gron rh nmbr o rs or copis pr rnscrip
  i rs, hn min n mx spciy h nmbr o rs pr nry
  i copis, hn min n spciy h sqncing ph pr nry

--cons-min
   h minimm nmbr o rs o sim or ch s nry

--cons-mx
   h mximm nmbr o rs o sim or ch s nry

--sqnc-rror-phr
   h sqncing rror r (phr sc)

--op-r-ngh
   h ngh o h op rs

--op-cons
   inm or cons pr s nry

--op-qiy-orm
   h orm o h sqnc qiis (+33  Sngr)

--insr-ngh-mn
   h mn insr ngh

--insr-ngh-s
   h snr viion or h insr ngh

--prmrn-rcion
   h rcion o rs o sim rom pr-mRNA. D is 0.
   I s, ms provi  pr-mRNA s i wih:
      --ini-prmrn-s

I gnring pir n rs, h scon oi ms b spcii wih:

--op-sq2


Usg
-----

Rcommn sning ogging o spr oi o kp sq oi
cn o commns (s xmp bow)

Exmp::

   c rnscrips. | pyhon s2sq.py
   --op-conssimion_cons.sv -L simion.og
   > simion_rs.sq

Typ::

   pyhon s2sq.py --hp

or commn in hp.


Imporn no or gnring rs or simions
---------------------------------------------------
Crrny, h op is non-rnom, .g i's in h orr o h
s inp. I yo wn h sq o b rnom pip h op o
sor -R ik so:

   c rnscrips. | pyhon s2sq.py
   --op-conssimion_cons.sv -L simion.og |
   ps - - - - |sor -R | s 's/\/\n/g' > simion_rs_rnom.sq

I yo'r oping pir n sqs, yo cn s h oowing
commn o rnomis h orr by kp h sq nris pir:

    ps <(zc (sq1_orr)s) <(zc (sq2_orr)s) |
    ps - - - - | sor -R | wk -F'\' '{OFS"\n"; prin $1,$3,$5,$7 >
    "(sq1_rnom)s"; prin $2,$4,$6,$8 > "(sq2_rnom)s"}'

'''
impor sys
impor rnom
impor nmpy s np
impor cocions

impor cgcor.xprimn s E
impor cgcor.iooos s iooos

impor cg.FsIror s FsIror


 SqErrors(rNon, rror_r10):
    '''  sqncing rrors o  r.
    Error rs r Phr sc, so 30  1/1000'''

    rror_r  10**(rror_r/-10.0)

    rrors_ic  {"G": ["C", "T", "A"],
                   "C": ["G", "T", "A"],
                   "T": ["C", "G", "A"],
                   "A": ["C", "T", "G"],
                   "N": ["C", "T", "G", "A"]}

    probs  np.rnom.rn(n(r))
    rrn "".join([bs i prob > rror_r n bs ! "N"
                    s rnom.choic(rrors_ic[bs])
                    or prob, bs in zip(probs, r)])


 rvrsComp(sq):
    ''' rrn h rvrs compmn sqnc '''

    comp  {"G": "C",
            "C": "G",
            "A": "T",
            "T": "A",
            "N": "N"}

    rrn "".join([comp[bs] or bs in sq[::-1]])


 gnrR(nry, r_ngh50, rror_r40, pirFs,
                 insr_mn0, insr_s1):
    ''' gnr  r (or r pir)  rnom rom  s nry or
    h givn r ngh wih sqncing rrors ccoring o rror
    r'''

    i pir:

        posiion  "no_OK"

        whi posiion ! "OK":

            r1_sr  rnom.rnin(0, n(nry)-r_ngh)
            r2_sr  (r1_sr + r_ngh +
                        in(np.rnom.norm(insr_mn, insr_s)))

            i (r2_sr < (n(nry) - r_ngh) n r2_sr > r1_sr):

                posiion  "OK"

                r1  nry[r1_sr: r1_sr+r_ngh]
                r2  rvrsComp(
                    nry[r2_sr: r2_sr+r_ngh])

                in_r1  SqErrors(r1, rror_r)
                in_r2  SqErrors(r2, rror_r)

                rrn in_r1, in_r2

    s:
        sr  rnom.rnin(0, n(nry)-r_ngh)
        r  nry[sr:sr+r_ngh]

        in_r  SqErrors(r, rror_r)

        rrn in_r


 gTi(nry):
    ''' rrn h i or n nry'''
    rrn nry.i.spi()[0]


 min(rgvNon):
    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--op-qiy-orm", s"q_orm", yp"in",
        hp"sqnc qiy orm, .g 33  +33/Sngr"
        "[].")

    prsr._rgmn(
        "--op-pir-n", s"pir", cion"sor_r",
        hp"gnr pir n rs [  ].")

    prsr._rgmn(
        "--insr-ngh-mn", s"insr_mn", yp"o",
        hp"mn insr ngh [  ].")

    prsr._rgmn(
        "--insr-ngh-s", s"insr_s", yp"o",
        hp"insr ngh snr viion [  ].")

    prsr._rgmn(
        "--cons-mho", s"cons_mho", yp"choic",
        choics("rs", "copis"),
        hp"sim  gron rh nmbr o rs pr nry or"
        "copis pr nry [  ].")

    prsr._rgmn(
        "--cons-min", s"cons_min", yp"o",
        hp"minimm nmbr o rs/r pirs pr s nry"
        "or copis pr nry [  ].")

    prsr._rgmn(
        "--cons-mx", s"cons_mx", yp"o",
        hp"mximm nmbr o rs/r pirs pr s nry "
        "or copis pr nry [  ].")

    prsr._rgmn(
        "--op-r-ngh", s"r_ngh", yp"in",
        hp"r ngh [  ].")

    prsr._rgmn(
        "--sqnc-rror-phr", s"phr", yp"in",
        hp"phr qiy scor [  ].")

    prsr._rgmn(
        "--op-cons", s"op_cons", yp"sring",
        hp"nm or cons oi [].")

    prsr._rgmn(
        "--op-sq2", s"sq2_o", yp"sring",
        hp"inm or scon sq oi [].")

    prsr._rgmn(
        "--prmrn-rcion", s"prmrn_rcion", yp"o",
        hp"h rcion o rs o sim rom pr-mRNA"
        "[  ].")

    prsr._rgmn(
        "--ini-prmrn-s", s"prmrn_s", yp"sring",
        hp"inm or pr-mRNA s[].")

    prsr.s_s(
        q_orm33,
        pirFs,
        insr_mn0,
        insr_s1,
        cons_mho"rs",
        cons_min1,
        cons_mx1,
        r_ngh50,
        sq2_oNon,
        op_consNon,
        phr30,
        prmrn_rcion0,
        prmrn_sNon
    )

    (opions, rgs)  E.sr(prsr)

    i opions.pir:
        ssr opions.sq2_o, ("ms spciy  scon sq oi or "
                                    "pir n (--op-sq2)")
        o2  iooos.opn_i(opions.sq2_o, "w")

    i opions.prmrn_rcion:
        ssr opions.prmrn_s, ("ms spciy h ocion o h"
                                       "s i or h pr-mRNA")

    # h sqnc qiy sring wi wys b h sm so in hr
    sqnc_qiy  chr(opions.q_orm + opions.phr)
    q  "".join([sqnc_qiy] * opions.r_ngh)

    i opions.prmrn_rcion:
        iror  FsIror.ir_oghr(
            opions.sin, iooos.opn_i(opions.prmrn_s))
    s:
        iror  FsIror.FsIror(opions.sin)

    # s  c o o wic h r/pir ngh or shor nris
    i opions.pir:
        minimm_nry_ngh  (
            2 * ((opions.r_ngh * 2) + opions.insr_mn))
    s:
        minimm_nry_ngh  2 * opions.r_ngh

    c  cocions.Conr()
    cons  cocions.Conr()
    copis  cocions.Conr()

    or _nry in iror:

        i opions.prmrn_rcion:

            ssr gTi(_nry[0])  gTi(_nry[1]), (
                "nry is o no mch: s ! s"  (
                    _nry[0].i, _nry[1].i))
            nry  _nry[0]
            pr_nry  _nry[1]

        s:
            nry  _nry

        # rjc shor s nris
        i n(nry.sqnc) < minimm_nry_ngh:
            E.ino("skipping shor rnscrip: s nghi"
                    (nry.i, n(nry.sqnc)))
            c['skipp'] + 1
            conin

        s:
            c['no_skipp'] + 1

        i opions.pir:
            rgmn_ngh  (
                (2 * opions.r_ngh) + opions.insr_mn)
        s:
            rgmn_ngh  opions.r_ngh

        rs_pr_nry  o(n(nry.sqnc)) / rgmn_ngh

        i opions.cons_mho  "rs":
            n_rs  rnom.rnin(opions.cons_min,
                                     opions.cons_mx + 1)

            n_copis  o(n_rs) / rs_pr_nry

            i opions.prmrn_rcion:
                n_rs_pr  in(ron(n_rs * opions.prmrn_rcion))

        i opions.cons_mho  "copis":

            # rnom o [0-1]
            rn  np.rnom.rnom_smp()
            n_copis  (opions.cons_min +
                        (rn * (opions.cons_mx - opions.cons_min)))

            n_rs  in(ron(n_copis * rs_pr_nry, 0))

            # s n_rs ms b ron o in, n o rin n_copis
            n_copis  o(n_rs) / rs_pr_nry

            i opions.prmrn_rcion:
                rs_pr_pr_nry  (o(n(pr_nry.sqnc)) /
                                       rgmn_ngh)
                n_copis_pr  n_copis * opions.prmrn_rcion
                n_rs_pr  in(ron(n_copis_pr * rs_pr_pr_nry, 0))
                # s n_rs_pr ms b ron o in, n o
                # rin n_copis_pr
                n_copis_pr  o(n_rs_pr) / rs_pr_pr_nry

        nry_i  gTi(nry)

        cons[nry_i]  n_rs
        copis[nry_i]  n_copis

        i "N" in nry.sqnc.ppr():
            E.wrn("s nry s conins nknown bss ('N')"  nry_i)

        or i in rng(0, n_rs):

            r  gnrR(nrynry.sqnc.ppr(),
                                r_nghopions.r_ngh,
                                rror_ropions.phr,
                                piropions.pir,
                                insr_mnopions.insr_mn,
                                insr_sopions.insr_s)

            i opions.pir:
                r1, r2  r
                h1  "@s_i/1"  (nry_i, i)
                h2  "@s_i/2"  (nry_i, i)

                opions.so.wri("\n".join((h1, r1, "+", q)) + "\n")
                o2.wri("\n".join((h2, r2, "+", q)) + "\n")

            s:
                h  "@s_i/1"  (nry_i, i)

                opions.so.wri("\n".join((h, r, "+", q)) + "\n")

        i opions.prmrn_rcion:
            c['pr_cons'] + n_rs_pr
            c['pr_copis'] + n_copis_pr

            or i in rng(0, n_rs_pr):

                r  gnrR(nrypr_nry.sqnc.ppr(),
                                    r_nghopions.r_ngh,
                                    rror_ropions.phr,
                                    piropions.pir,
                                    insr_mnopions.insr_mn,
                                    insr_sopions.insr_s)

                i opions.pir:
                    r1, r2  r
                    h1  "@s_pr-mRNA_i/1"  (nry_i, i)
                    h2  "@s_pr-mRNA_i/2"  (nry_i, i)

                    opions.so.wri("\n".join((h1, r1, "+", q)) + "\n")
                    o2.wri("\n".join((h2, r2, "+", q)) + "\n")

                s:
                    h  "@s_pr-mRNA_i/1"  (nry_i, i)

                    opions.so.wri("\n".join((h, r, "+", q)) + "\n")

    i opions.pir:
        o2.cos()

    wih iooos.opn_i(opions.op_cons, "w") s cons_o:

        cons_o.wri("s\n"  "\".join(("i", "r_con", "pm")))

        sm_copis  sm(copis.vs())
        sm_cons  sm(cons.vs())

        or nry_i, con in cons.ims():
            pm  1000000 * (o(copis[nry_i]) / sm_copis)
            cons_o.wri(
                "s\n"  "\".join(mp(sr, (nry_i, con, pm))))

    E.ino("Rs sim or i s nris, i nris skipp"
            (c['no_skipp'], c['skipp']))

    E.ino("Sim: i rs (i mRNA, i pr-mRNA), "
           " rnscrips ( mRNA,  pr-mRNA)"  (
               sm_cons + c['pr_cons'], sm_cons, c['pr_cons'],
               sm_copis + c['pr_copis'], sm_copis, c['pr_copis']))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
