'''
b2s.py - g sqncs rom b i


:Tgs: Gnomics Inrvs Sqncs Convrsion BED FASTA



Prpos
-------

This scrip ops ncoi sqncs or inrvs wihin
 :rm:`b` orm i sing  corrsponing gnom i.

Usg
-----

A rqir inp o b2s.py is  cg inx gnom. To obin n
ix hmn rrnc gnom w wo yp

Exmp::
   c hg19.s | inx_s.py hg19 > hg19.og

This i wo hn srv s h --gnom-i whn w wish o xrc
sqncs rom  :rm:`b` orm i.


For xmp w co now yp::

   c in.b | pyhon b2s.py --gnom-i hg19 > o.s

Whr w k  s o gnomic inrvs (.g. rom  hmn ChIP-sq xprimn)
n op hir rspciv ncoi sqncs.


Typ::

   pyhon b2s.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor cgcor.xprimn s E
impor cg.B s B
impor cg.InxFs s InxFs
impor cg.Mskr s Mskr


 min(rgvNon):
    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnomic sqnc o rriv "
                      "sqncs rom.")

    prsr._rgmn("-m", "--mskr", s"mskr", yp"choic",
                      choics("s", "smskr", "somsk", "non"),
                      hp"ppy mskr o msk op sqncs "
                      "[].")

    prsr._rgmn("--op-mo", s"op_mo", yp"choic",
                      choics("inrvs", "righ", "sgmns"),
                      hp"wh o op. "
                      "'inrvs' gnrs  sing sqnc or "
                      "ch b inrv. 'righ' gnrs wo "
                      "sqncs, on in ch ircion, or ch b "
                      "inrv. 'sgmns' cn b s o op "
                      "sqnc rom b12 is so h sqnc ony covrs "
                      "h sgmns []")

    prsr._rgmn("--min-sqnc-ngh", s"min_ngh", yp"in",
                      hp"rqir  minimm sqnc ngh []")

    prsr._rgmn("--mx-sqnc-ngh", s"mx_ngh", yp"in",
                      hp"rqir  mximm sqnc ngh []")

    prsr._rgmn(
        "--xn-", s"xn_", yp"choic",
        choics("non", "3", "5", "boh", "3ony", "5ony"),
        hp"xn  3', 5' or boh or no ns. I 3ony or 5ony "
        "r s, ony h  sqnc is rrn []")

    prsr._rgmn(
        "--xn-by", s"xn_by", yp"in",
        hp"xn by # bss []")

    prsr._rgmn(
        "--s-srn", s"ignor_srn",
        cion"sor_s",
        hp"s srn inormion n rrn rvrs compmn "
        "on inrvs oc on h ngiv srn. "
        "[]")

    prsr.s_s(
        gnom_iNon,
        mskrNon,
        op_mo"inrvs",
        min_ngh0,
        mx_ngh0,
        xn_Non,
        xn_by100,
        ignor_srnTr,
    )

    (opions, rgs)  E.sr(prsr)

    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
        conigs  s.gConigSizs()
        s.sConvrr(InxFs.gConvrr("zro-boh-opn"))

    conr  E.Conr()
    is, sqs  [], []

    E.ino("cocing sqncs")
    or b in B.sNm(B.iror(opions.sin)):
        conr.inp + 1

        conig  s.gLngh(b.conig)

        i opions.ignor_srn:
            srn  "+"
        s:
            srn  b.srn

        i opions.op_mo  "sgmns" n b.comns  12:
            is.ppn("s s:i..i (s) s s" 
                       (b.nm, b.conig, b.sr, b.n, srn,
                        b["bockSizs"], b["bockSrs"]))
            sg_sqs  [s.gSqnc(b.conig, srn, sr, n)
                        or sr, n in b.oInrvs()]
            sqs.ppn("".join(sg_sqs))

        i (opions.op_mo  "inrvs" or
              opions.op_mo  "sgmns"):
            is.ppn("s s:i..i (s)" 
                       (b.nm, b.conig, b.sr, b.n, srn))
            sqs.ppn(
                s.gSqnc(b.conig, srn, b.sr, b.n))

        i opions.op_mo  "righ":
              b.n - b.sr

            sr, n  mx(0, b.sr - ), b.n - 
            is.ppn("s_ s:i..i (s)" 
                       (b.nm, b.conig, sr, n, srn))
            sqs.ppn(s.gSqnc(b.conig, srn, sr, n))

            sr, n  b.sr + , min(conig, b.n + )
            is.ppn("s_r s:i..i (s)" 
                       (b.nm, b.conig, sr, n, srn))
            sqs.ppn(s.gSqnc(b.conig, srn, sr, n))

    E.ino("coc i sqncs"  n(sqs))

    msk  Mskr.mskSqncs(sqs, opions.mskr)
    opions.so.wri(
        "\n".join([">s\ns"  (x, y) or x, y in zip(is, msk)]) + "\n")

    E.ino("msk i sqncs"  n(sqs))

    conr.op  n(sqs)

    E.ino("s"  conr)

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
