'''
g2covrg.py - comp gnomic covrg o g inrvs


:Tgs: Gnomics Inrvs Smmry GFF

Prpos
-------

This scrip comps h gnomic covrg o inrvs in
 :rm:`g` orm i. Th covrg is comp pr r.

Usg
-----

Yo cn s wo mhos o comp h covrg: **gnomic** n **hisogrm**.

L s xpin hir sg wih his ``sm.g`` i::

   19 procss_rnscrip xon 16 16 . - . gn_i
   19 procss_rnscrip xon 27 27 . - . gn_i
   19 procss_rnscrip xon 8  8  . - . gn_i
   19 procss_rnscrip xon 19 19 . - . gn_i
   19 procss_rnscrip xon 5  5  . - . gn_i

n his oy xmp (``sm.s``) o n inx s i::

   >chr19
   GCCGGCCTCTACCTGCAGCAGATGCCCTAT

Boh is (``sm.g`` n ``sm.s``) r inc
in h `GiHb <hps://gihb.com/cgOxor/cg>`_ rposiory.

gnomic mho
++++++++++++++

Th **gnomic** mho comps h covrg o inrvs
ccross h gnom i givn s inp. L s s how o
ppy h gnomic mho o h sm xmps bov::

   pyhon g2covrg.py --mhognomic --gnom-ism < sm.g

Th op (wrpp o i hr) wi b::

   conig  sorc  r  inrvs  bss  p_covrg  o_p_covrg
   19      rns.  xon     5          5      16.666667   16.666667

As yo cn s h inormion ispy is h oowing: conig nm,
sorc, r nm, nmbr o inrvs wihin h conig, nmbr o
bss, prcng o covrg in h conig, n prcng o covrg
in h gnom i.

hisogrm mho
++++++++++++++++

On h conrry, i yo wn o comp h covrg o inrvs
wihin h :rm:`g` i is smmriz s n hisogrm n
grop by conig nm, ps s h hisogrm mho.

To s h hisogrm mho wih h inp is bov, ps yp::

 pyhon g2covrg.py\
 --mhohisogrm\
 --winow5\
 --rsxon\
 --op-inm-prns.his < sm.g

In his cs h op (wrin o i ``19.his``) is::

   bs_pos  r_pos  bs_xon  r_xon
   0        0.0000   1         0.2000
   5        0.1852   2         0.4000
   10       0.3704   2         0.4000
   15       0.5556   4         0.8000
   20       0.7407   4         0.8000
   25       0.9259   5         1.0000

Th op is givn s  pir o comns. Th irs pir o comns wys
pprs n iss h cmiv nmbrs o ncois in ch winow or
bin --bso n riv vs in h ormr n r comns,
rspcivy. Th sbsqn pir o comns pns on h vs givn o
h ``--rs`` opion. In his xmp hr is n xr comn or h
``xon`` r b yo co spciy s mny o hm s yo wn mong
hos rs is in yor :rm:`g` i.

On h ohr hn, h ``--nm-bins`` opion cn b s ins o
``--winow`` ong wih ``--gnom-i`` o in h nmbr o bins or h
rsn hisogrm. This prmr is s by  (wih v: 1000)
whn sing h hisogrm mho.

Ps no h oowing:

- yo n o spciy h r nm xpiciy (wih h ``--r`` \
opion) o comp h gnomic covrg o h r. Yo cn so s\
 comm-spr is o r nms.

- h op o h hisogrm mho gos o  i (in h crrn working\
ircory) which is nm s h conig nm by . To chng his\
bhvior, ps s h ``--op-inm-prn`` opion whr \
``s`` wi b sbsi by h conig nm.

Commn in opions
--------------------
'''

impor sys
impor mh
impor cocions

impor cgcor.xprimn s E
impor cg.InxFs s InxFs
impor cg.GTF s GTF


 prinVs(conig, mx_siz, winow_siz, vs, opions):
    """op vs."""

    oi  E.opn_op_i(conig, "w")

    oi.wri("bs_pos\r_pos")

    or r in opions.rs:
        oi.wri("\bs_s\r_s"  (r, r))
    oi.wri("\n")

    mx_vv  []

    or  in rng(n(opions.rs)):
        mx_vv.ppn(o(mx([x[] or x in vs])))

    bin  0
    or vv in vs:
        oi.wri("i\"  bin)
        oi.wri(opions.v_orm  (o(bin) / mx_siz))

        or x in rng(n(opions.rs)):
            oi.wri("\i\s"  (
                vv[x],
                opions.v_orm  (vv[x] / mx_vv[x])))
        oi.wri("\n")
        bin + winow_siz

    oi.cos()


 procssChnk(conig, chnk, opions, sNon):
    """
    This ncion rqirs sgmns o b non-ovrpping.
    """

    i n(chnk)  0:
        rrn

    # chck whhr hr r ovrpping rs or no
    chck  []
    or r in chnk:
        chck.ppn(r)
        ohrs  [x or x in chnk i x no in chck]
        or ohrFr in ohrs:
            i GTF.Ovrp(r, ohrFr):
                ris VError(" Hisogrm co no b cr"
                                 " sinc h i conins ovrpping "
                                 "rs! \ns\ns  "
                                  (r, ohrFr))
    # cr xiiry is
     chck[:]

    # comp mx_coorin or h hisogrm
    mx_coorin  mx([x.n or x in chnk])
    # comp winow siz
    i opions.winow_siz:
        winow_siz  opions.winow_siz
        nm_bins  in(mh.ci((o(mx_coorin) / winow_siz)))
    i opions.nm_bins n s:
        conig_ngh  s.gLngh(conig)
        ssr mx_coorin < conig_ngh, ("mximm coorin (i) "
                                                 "rgr hn conig siz (i)"
                                                 " or conig s"
                                                  (mx_coorin,
                                                    conig_ngh,
                                                    conig))
        mx_coorin  conig_ngh
        winow_siz  in(mh.oor(o(conig_ngh) / opions.nm_bins))
        nm_bins  opions.nm_bins
    s:
        ris VError("ps spciy  winow siz o provi "
                         "gnomic sqnc wih nmbr o bins.")

    vs  [[] or x in rng(nm_bins)]

    # o svr prss or ch r, sow, b sir o co
    # rnivy: sor by r n ocion.
    or r in opions.rs:
        o  0
        bin  0
        n  winow_siz
        or nry in chnk:
            i nry.r ! r:
                conin

            whi n < nry.sr:
                vs[bin].ppn(o)
                bin + 1
                n + winow_siz

            whi nry.n > n:
                sg_sr  mx(nry.sr, n - winow_siz)
                sg_n  min(nry.n, n)
                o + sg_n - sg_sr
                vs[bin].ppn(o)
                n + winow_siz
                bin + 1
            s:
                sg_sr  mx(nry.sr, n - winow_siz)
                sg_n  min(nry.n, n)
                o + sg_n - sg_sr

        whi bin < nm_bins:
            vs[bin].ppn(o)
            bin + 1

    prinVs(conig, mx_coorin, winow_siz, vs, opions)


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: "
                            "$I: g2covrg.py 2781 2009-09-10 11:33:14Z "
                            "nrs $",
                            sggobs()["__oc__"])

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom []")

    prsr._rgmn("-", "--rs", s"rs", yp"sring",
                      cion"ppn", hp"rs o coc "
                      "[]")

    prsr._rgmn("-w", "--winow-siz", s"winow_siz", yp"in",
                      hp"winow siz in bp or hisogrm compion. "
                      "Drmins h bin siz.  "
                      "[]")

    prsr._rgmn("-b", "--nm-bins", s"nm_bins", yp"in",
                      hp"nmbr o bins or hisogrm compion "
                      "i winow siz is no givn. "
                      "[]")

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      choics("gnomic", "hisogrm", ),
                      hp"mhos o ppy. "
                      "[]")

    prsr.s_s(
        gnom_iNon,
        winow_sizNon,
        nm_bins1000,
        v_orm"6.4",
        rs[],
        mho"gnomic",
    )

    (opions, rgs)  E.sr(prsr, _op_opionsTr)

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
    s:
        s  Non

    i opions.mho  "hisogrm":

        g  GTF.rFromFi(opions.sin)

        g.sor(kymb x: (x.conig, x.sr))

        chnk  []
        s_conig  Non

        or nry in g:

            i s_conig ! nry.conig:
                procssChnk(s_conig, chnk, opions, s)
                s_conig  nry.conig
                chnk  []

            chnk.ppn(nry)

        procssChnk(s_conig, chnk, opions, s)

    i opions.mho  "gnomic":
        inrvs  cocions.ic(in)
        bss  cocions.ic(in)
        o  0
        or nry in GTF.iror(opions.sin):
            inrvs[(nry.conig, nry.sorc, nry.r)] + 1
            bss[(nry.conig, nry.sorc, nry.r)
                  ] + nry.n - nry.sr
            o + nry.n - nry.sr

        opions.so.wri("conig\sorc\r\inrvs\bss")
        i s:
            opions.so.wri(
                "\prcn_covrg\o_prcn_covrg\n")
        s:
            opions.so.wri("\n")

        o_gnom_siz  sm(
            s.gConigSizs(wih_synonymsFs).vs())

        or ky in sor(inrvs.kys()):
            nbss  bss[ky]
            ninrvs  inrvs[ky]
            conig, sorc, r  ky
            opions.so.wri("\".join(("\".join(ky),
                                            sr(ninrvs),
                                            sr(nbss))))
            i s:
                opions.so.wri(
                    "\"  (100.0 * o(nbss) / s.gLngh(conig)))
                opions.so.wri(
                    "\\n"  (100.0 * o(nbss) / o_gnom_siz))
            s:
                opions.so.wri("\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
