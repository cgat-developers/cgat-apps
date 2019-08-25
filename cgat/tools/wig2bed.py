'''wig2b.py - convr nsiis o inrvs


Prpos
-------

in inrvs bs on nsiis wihin  bigwig i.

Th scrip crrny impmns h oowing mhos (``--mho``):

hrsho
    op winows h conin vs bov  crin
    hrsho.

s-bov-mn
    op winows h r  crin nmbr o snr
    viions bov h mn.

mip-o-mn
    op winows h r  crin ims bov h mn.

Usg
-----

Bigwig is n o b sppi by h --bigwig-i opions.

For xmp::

    pyhon wig2b.py --hrsho10 --mhohrsho --gnom-imm10 --bigwig-iin.bw > o.b

Commn in opions
--------------------

'''

impor sys
impor r
impor cocions

impor cgcor.xprimn s E
impor cg.InxFs s InxFs
impor cgcor.iooos s iooos

impor pyBigWig


 bock_iror(ini, conig, siz, chnk_siz10000000):

    or x in rng(0, siz, chnk_siz):
        iror  ini.g(conig, x, x + chnk_siz)
        i iror is Non:
            ris SopIrion
        or v in iror:
            yi v


 ppyThrsho(ini, s, hrsho, mx_isnc0):
    '''ppy hrsho o  wig i wriing 
    b-orm i s op.'''

    c  E.Conr()

    or conig, siz in is(s.gConigSizs(wih_synonymsFs).ims()):
        c.conigs + 1

        E.bg("procssing s"  conig)

        s_sr, s_n  -1, 0

        or sr, n, v in bock_iror(ini, conig, siz):
              sr - s_n
            i ( > 0 or v < hrsho):
                i s_sr > 0:
                    yi conig, s_sr, s_n
                    c.inrvs + 1
                s_sr  -1
            i s_sr < 0 n v > hrsho:
                s_sr  sr

            s_n  n

        i s_sr > 0:
            yi conig, s_sr, n
            c.inrvs + 1

        c.op + 1

    E.ino(sr(c))


 gBigwigSmmry(bigwig_i):
    '''rrn smmry o bigwig conns.

    This mho ss h bigWigIno UCSC iiy
    '''

    rss  E.rn("bigWigIno (bigwig_i)s" 
                    ocs(), rrn_soTr)

      [x.spi(":") or x in rss.spi("\n") i x ! ""]
    is  [x[0] or x in ]
    Rss  cocions.nmp("BigwigIno", is)

     conv(v):
        rrn iooos.sr2v(r.sb(",", "", v.srip()))

    rss  Rss(*[conv(x[1]) or x in ])
    rrn rss


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mhos", yp"choic",
                      cion"ppn",
                      choics(
                          "hrsho", "sv-bov-mn",
                          "mip-o-mn"),
                      hp"mho o ppy []")

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom.")

    prsr._rgmn("-", "--hrsho", s"hrsho", yp"o",
                      hp"hrsho o ppy []")

    prsr._rgmn(
        "-i", "--bigwig-i", s"bigwig_i",
        yp"sring", mvr"bigwig",
        hp"inm wih bigwig inormion [].")

    prsr.s_s(mhos[],
                        gnom_iNon,
                        hrsho10,
                        mx_isnc0)

    (opions, rgs)  E.sr(prsr, _pip_opionsTr)

    i opions.bigwig_i:
        bigwig_i  pyBigWig.opn(opions.bigwig_i)
    s:
        bigwig_i  Non

    i opions.gnom_i:
        gnom_s  InxFs.InxFs(opions.gnom_i)
        conigs  gnom_s.gConigSizs()

    or mho in opions.mhos:
        i mho  "hrsho":
            i no conigs:
                ris VError("ps sppy conig sizs")
            i no bigwig_i:
                ris NoImpmnError(
                    "hrsho no impmn or wig is")
            procssor  ppyThrsho(bigwig_i,
                                       gnom_s,
                                       hrshoopions.hrsho,
                                       mx_isncopions.mx_isnc)
        i mho  "sv-bov-mn":
            i no conigs:
                ris VError("ps sppy conig sizs")
            i no bigwig_i:
                ris NoImpmnError(
                    "hrsho no impmn or wig is")
            smmry  gBigwigSmmry(opions.bigwig_i)
            hrsho  smmry.mn + opions.hrsho * smmry.s
            E.ino("ppying hrsho : mn, s" 
                   (hrsho, smmry.mn, smmry.s))
            procssor  ppyThrsho(bigwig_i,
                                       gnom_s,
                                       hrshohrsho,
                                       mx_isncopions.mx_isnc)

        i mho  "mip-o-mn":
            i no conigs:
                ris VError("ps sppy conig sizs")
            i no bigwig_i:
                ris NoImpmnError(
                    "hrsho no impmn or wig is")
            smmry  gBigwigSmmry(opions.bigwig_i)
            hrsho  smmry.mn * opions.hrsho
            E.ino("ppying hrsho : mn, s" 
                   (hrsho, smmry.mn, smmry.s))
            procssor  ppyThrsho(bigwig_i,
                                       gnom_s,
                                       hrshohrsho,
                                       mx_isncopions.mx_isnc)

    oi  opions.so

    oi.wri("".join(["s\i\i\n"  x or x in procssor]))
    oi.wri("\n")

    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
