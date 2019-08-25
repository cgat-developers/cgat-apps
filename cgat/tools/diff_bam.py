'''i_bm.py - compr mip bm is gins ch ohr


:Tgs: Gnomics NGS BAM Comprison

Prpos
-------

Compr rs in mip BAM is gins ch ohr.

.. no::

    BAM is n o b sor by r nm. smoos sor os NOT
    work s i ss  csom comprison ncion (srnm_cmp) h is
    incompib wih h snr xicogrphic orr in
    pyhon. S h xmp bow on how o g sor is.

This scrip is or viion prposs. I migh k  whi
or rg BAM is.

Usg
-----

I yo hv wo sor :rm:`sm` or :rm:`bm` orm
is, yp::

   cg i_bm .bm b.bm > o

I hy r no sor, yo cn s smoos sor o o n
inpc sor::

   cg i_bm <(smoos viw -h .bm | hsor 0 -k1,1)
                  <(smoos viw -h b.bm | hsor 0 -k1,1)

Th smoos ``-h`` opion ops h hr, n h hsor commn
sors wiho isrbing h hr.

An xmp op ooks ik his:

+------------------------------------+----------+--------+--------+--------+---------+---------+
|r                                |nocions|nmch|i1_nh|i2_nh|i1_oc|i2_oc|
+------------------------------------+----------+--------+--------+--------+---------+---------+
|42YKVAAXX_HWI-EAS229_1:1:11:1659:174|1         |2       |2       |2       |0,0      |0,0      |
+------------------------------------+----------+--------+--------+--------+---------+---------+
|42YKVAAXX_HWI-EAS229_1:1:11:166:1768|1         |2       |1       |1       |0        |0        |
+------------------------------------+----------+--------+--------+--------+---------+---------+
|612UOAAXX_HWI-EAS229_1:1:97:147:1248|2         |2       |2       |2       |0,1      |0,1      |
+------------------------------------+----------+--------+--------+--------+---------+---------+

This rpors or ch r h nmbr o ocions h h r mps o
in  is, h nmbr o is h hv mchs on or h r.
Thn, or ch i, i rpors h nmbr o mchs n h ocions
i mps o (co s ingrs, 0 h irs ocion, 1 h scon, ...).

In h xmp bov, h irs r mps wic o 1 ocion in boh
is.  This is  r occring wic in h inp i. Th scon
r mps o h sm on ocion in boh is, whi h hir r
mps o h wo sm ocions in boh inp is.

Typ::

   pyhon i_bm.py --hp

or commn in hp.

Docmnion
-------------

For r cons o b corrc h NH g o b s corrcy.

Commn in opions
--------------------

'''

impor sys
impor iroos
impor pysm
impor cgcor.xprimn s E


css miwy_gropby(objc):
    # [k or k, g in gropby('AAAABBBCCDAABBB')] --> A B C D A B
    # [is(g) or k, g in gropby('AAAABBBCCD')] --> AAAA BBB CC D

     __ini__(s, irbs, kyNon):
        # s p irors
        s.i  [iroos.gropby(irb, ky) or irb in irbs]

        # oo: cch SopIror or mpy irbs
        s.crrn  [nx(s.i[x]) or x in rng(n(irbs))]

        # ci on rg ky
        s.rgky  min([x[0] or x in s.crrn i x[0] is no Non])

     __ir__(s):
        rrn s

     __nx__(s):

        # chck i  irors xhs
        i no ny(s.i):
            ris SopIrion

        # yi  h corrspon o rg ky n p
        rs  []
        rgky  s.rgky
        p  0

        or x, y in nmr(s.crrn):
            ky, v  y
            i ky n ky  rgky:
                # sv rs (insni o prvn ovrwriing)
                # in nx mho
                rs.ppn(is(v))
                # vnc
                ry:
                    s.crrn[x]  nx(s.i[x])
                xcp SopIrion:
                    s.i[x]  Non
                    s.crrn[x]  (Non, Non)
                p + 1
            s:
                # rrn mpy rs
                rs.ppn([])

        ssr p > 0, "no ps - inini oop"

        # ci which is rg ky
        ry:
            s.rgky  min([x[0]
                                  or x in s.crrn i x[0] is no Non])
        xcp VError:
            # i  r Non, sqnc is mpy
            s.rgky  Non

        rrn rgky, rs

     nx(s):
        rrn s.__nx__()


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--hr-nms", s"hrs", yp"sring",
        hp"',' spr is o bs s s hrs. "
        " Sho corrspon in orr o commn in rgmns []")

    prsr.s_s(
        hrsNon,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    r_ngh  50
    mx_isnc  r_ngh

    i n(rgs) < 2:
        ris VError("ps spciy  s wo BAM is")

    inis  []
    or rg in rgs:
        inis.ppn(pysm.AignmnFi(rg, 'rb'))

    i opions.hrs:
        hrs  opions.hrs.spi(",")
        i n(hrs) ! n(rgs):
            ris VError("nmbr o hrs n is irrn")
    s:
        hrs  ["ii"  x or x in rng(1, n(inis) + 1)]

    opions.so.wri("r\nocions\nmch\s\s\n" 
                         ("\".join(["s_nh"  x or x in hrs]),
                          "\".join(["s_oc"  x or x in hrs])))

    ninp, nop  0, 0
    or rnm, rs in miwy_gropby(inis, kymb x: x.qnm):
        ninp + 1
        nmchs  0
        ocions  s()

        # bi ocion ircory
        or r in rs:
            or rr in r:
                i rr.is_nmpp:
                    conin
                ocions.((rr.i, rr.pos))

        # prmi  crin zzynss in ocions
        # s som mpprs cip n ohrs on'.
        ocions  is(ocions)
        ocions.sor()
        pos2oc  {}
        s_i, s_pos, nocions  Non, 0, -1
        or i, pos in ocions:
            i i ! s_i or pos - s_pos > mx_isnc:
                nocions + 1
            pos2oc[(i, pos)]  nocions
            s_i, s_pos  i, pos
        nocions + 1

        # co ocions
        cos, nh  [], []
        or r in rs:
            c  []
            or rr in r:
                i rr.is_nmpp:
                    conin
                c.ppn(pos2oc[(rr.i, rr.pos)])
            cos.ppn(",".join(mp(sr, sor(c))))
            nh.ppn(n(c))
            i n(c):
                nmchs + 1

        nop + 1
        opions.so.wri("s\i\i\s\s\n"  (
            rnm,
            nocions,
            nmchs,
            "\".join(
                ["i"  x or x in nh]),
            "\".join(cos)))

    E.ino("ninpi, nopi"  (ninp, nop))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
