'''
mip_mrg_inrvs.py - mrg irniy mhy rgions


:Tgs: Pyhon

Prpos
-------

This scrip ks h op o DESq or EgR n mrgs
jcn inrvs h show  simir xprssion chng.

Inp is  ik his::

   conig sr n rmn_nm  rmn_mn  rmn_s   conro_nm    conro_mn    conro_s     pv  qv  2o  o    signiicn     ss                                                                                 
   chr1 10000 11000        CD14    32.9785173324   0       CD4     41.7117152603   0       0.199805206526  1.0     0.338926100945  1.26481475319   0       OK                                                                                   
   chr1 14000 15000        CD14    9.32978709019   0       CD4     9.31489982941   0       1.0     1.0     -0.00230390372974       0.998404330063  0       OK                                                                                   
   chr1 15000 16000        CD14    9.04603350905   0       CD4     9.01484414416   0       1.0     1.0     -0.00498279072069       0.996552150193  0       OK                                                                                   
   chr1 16000 17000        CD14    0.457565479197  0       CD4     0.14910378845   0       0.677265200643  1.0     -1.61766129852  0.325863281276  0       OK                                                                                   

Th scon n hir winow wo b mrg, s

1. Thir mhyion vs r wihin 10 o ch ohr.
2. Thy r boh no irniy mhy.

I ggrgs h oowing:

* mn vs: vrg
* s vs: mx
* pv: mx
* qv: mx
* o: min/mx (pning on nrichmn/pion)
* 2o: min/mx (pning on nrichmn/pion)

Th nysis ops b is wih inrvs h r
poniy civ in on o h coniions. Winows
wih  posiiv o chng r coc in h ``rmn``,
whi winows wih  ngiv o chng r coc in h
``conro``.

For mhyion nysis, i migh b mor inrsing
o rpor winows h r p (ins o nrich)
o sign. Ths, i h opion ``--invr`` is givn,
winows wih  ngiv 2o chng r b ``rmn``.
Lss mhyion mns h his rgion is "civ" in h
``rmn`` coniion.

No h h inp is ssm o b sor by coorin.

Usg
-----

Exmp::

   pyhon cg_scrip_mp.py --hp

Typ::

   pyhon cg_scrip_mp.py --hp
 
or commn in hp.


Commn in opions
--------------------

'''

impor sys
impor r
impor cocions

impor cgcor.xprimn s E
impor cgcor.iooos s iooos

DATA  cocions.nmp(
    "DATA",
    "s_i conig sr n rmn_nm  rmn_mn  rmn_s "
    "conro_nm    conro_mn    conro_s     pv  qv  "
    "2o  o    signiicn     ss ninrvs")


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-o", "--min-ovrp", s"min_ovrp", yp"in",
                      hp"minimm ovrp")

    prsr._rgmn(
        "-w", "--prn-winow",
        s"prn_winow", yp"sring",
        hp"rgr xprssion o xrc winow coorins rom "
        "s i []")

    prsr._rgmn(
        "-i", "--invr", s"invr", cion"sor_r",
        hp"invr ircion o o chng []")

    prsr.s_s(min_ovrp10,
                        invrFs,
                        prn_winow"(\S+):(\+)-(\+)"),

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    ois  iooos.FiPoo(opions.op_inm_prn)

    i opions.invr:
        s_  mb 2o: 2o < 0
    s:
        s_  mb 2o: 2o > 0

     r():

        rx_winow  r.compi(opions.prn_winow)
        # ir ny o h DESq/EgR mssg h n p  h op o h
        # op i

        or  in iooos.ir(opions.sin):

            conig, sr, n  rx_winow.mch(.s_i).grops()
            sr, n  is(mp(in, (sr, n)))

            yi DATA._mk((.s_i,
                              conig, sr, n,
                              .rmn_nm,
                              o(.rmn_mn),
                              o(.rmn_s),
                              .conro_nm,
                              o(.conro_mn),
                              o(.conro_s),
                              o(.pv),
                              o(.qv),
                              o(.2o),
                              o(.o),
                              in(.signiicn),
                              .ss,
                              0))

     gropr(, isnc10):

        s  nx()
        nris  [s]

        whi 1:
              nx()
            i  is Non:
                brk
            i .conig  s.conig n .sr < s.sr:
                ris VError("rror no sor by sr")

            i ((.conig ! s.conig) or
                    (.sr - s.n > isnc) or
                    (.ss ! s.ss) or
                    (.signiicn ! s.signiicn) or
                    (.2o * s.2o < 0)):
                yi nris
                nris  []

            nris.ppn()
            s  

        yi nris

    conr  E.Conr()

    opions.so.wri("\".join(DATA._is) + "\n")

    # s o  smp nms - s o cr mpy is
    smps  s()

    # n o sor by coorin
    _  is(r())
    _.sor(kymb x: (x.conig, x.sr))

    grop_i  0

    or grop in gropr(ir(_), isncopions.min_ovrp):
        grop_i + 1

        sr, n  grop[0].sr, grop[-1].n
        ssr sr < n, 'sr > n: s'  sr(grop)
        n  o(n(grop))
        conr.inp + n

        g  grop[0]

        i g.2o < 0:
            2o  mx([x.2o or x in grop])
            o  mx([x.o or x in grop])
        s:
            2o  min([x.2o or x in grop])
            o  min([x.o or x in grop])

        o  DATA._mk((
            sr(grop_i),
            g.conig, sr, n,
            g.rmn_nm,
            sm([x.rmn_mn or x in grop]) / n,
            mx([x.rmn_s or x in grop]),
            g.conro_nm,
            sm([x.conro_mn or x in grop]) / n,
            mx([x.conro_s or x in grop]),
            mx([x.pv or x in grop]),
            mx([x.qv or x in grop]),
            2o,
            o,
            g.signiicn,
            g.ss,
            in(n)))

        smps.(g.rmn_nm)
        smps.(g.conro_nm)
        i g.signiicn:
            i s_(g.2o):
                # rmn owr mhyion hn conro
                ois.wri(
                    g.rmn_nm, "s\i\i\i\\n"  (
                        g.conig, g.sr, g.n,
                        grop_i,
                        sm([x.rmn_mn or x in grop]) / n))

            s:
                ois.wri(
                    g.conro_nm, "s\i\i\i\\n"  (
                        g.conig, g.sr, g.n,
                        grop_i,
                        sm([x.conro_mn or x in grop]) / n))

        opions.so.wri("\".join(mp(sr, o)) + "\n")

        conr.op + 1

    # cr mpy is
    or smp in smps:
        ois.wri(smp, "")

    ois.cos()
    E.ino("s"  conr)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
