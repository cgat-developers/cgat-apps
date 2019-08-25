'''sqs2sq.py - mrg rs in sq is


:Tgs: Gnomics NGS FASTQ FASTQ Mnipion

Prpos
-------

This scrip ks wo pir-n sq is n ops  sing
sq i in which rs hv mrg.

Th wo is ms b sor by r iniir.

No h his scrip is crrny  proo-o-princip impmnion
n hs no bn opimiz or sp or ncioniy.

Usg
-----

Exmp::

   pyhon sqs2sq.py myRs.1.sq.gz myRs.2.sq.gz
          --mhojoin
          > join.sq

In his xmp w k  pir o sq is, join h rs n sv
h op in :i:`join.sq`.

Typ::

   pyhon sqs2sq.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor cocions
impor copy

impor cgcor.iooos s iooos
impor cgcor.xprimn s E
impor cg.Fsq s Fsq
impor cg.Gnomics s Gnomics


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--mho", s"mho", yp"choic",
                      choics('join', ),
                      hp"mho o ppy [].")

    prsr.s_s(
        mho"join",
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs) ! 2:
        ris VError(
            "ps sppy  s wo sq is on h commnin")

    n1, n2  rgs
    c  E.Conr()
    oi  opions.so

    i opions.mho  "join":
        # mrg bs on igons in opo
        ir1  Fsq.ir(iooos.opn_i(n1))
        ir2  Fsq.ir(iooos.opn_i(n2))
        p_siz  2
        or , righ in zip(ir1, ir2):
            c.inp + 1

            # bi icionry o ps
            s1, q1  .sq, .qs
              cocions.ic(is)
            or x in rng(n(s1) - p_siz):
                [s1[x:x + p_siz]].ppn(x)

            s2, q2  righ.sq, righ.qs
            s2  Gnomics.rvrs_compmn(s2)
            q2  q2[::-1]

            # comp is o oss/igons
            oss  cocions.ic(in)
            or x in rng(n(s2) - p_siz):
                c  s2[x:x + p_siz]
                or y in [c]:
                    oss[x - y] + 1

            # in mximm igon
            sor  sor([(y, x) or x, y in is(oss.ims())])
            mx_con, mx_os  sor[-1]

            E.bg('s: mximm os  i'  (.iniir,
                                                  mx_os))

            # simp mrg sqnc
            k  n(s2) - mx_os
            mrg_sq  s1 + s2[k:]

            # simp mrg qiy scors
            mrg_qs  q1 + q2[k:]

            nw_nry  copy.copy()
            nw_nry.sq  mrg_sq
            nw_nry.qs  mrg_qs
            oi.wri(nw_nry)
            c.op + 1

    # wri oor n op bnchmrk inormion.
    E.ino("s"  sr(c))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
