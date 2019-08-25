'''
vcss_sqi.py - rorm op o vc-ss or bs oing


:Tgs: Pyhon

Prpos
-------

cr  csv spr i or oing ino  bs rom 
op o vc-ss iiy in vc-oos pckg.

Usg
-----

Exmp::

   pyhon vcss_sqi.py [is] > [oi]

Typ::

   pyhon vcss_sqi.py --hp

or commn in hp.


Commn in opions
--------------------

'''
impor os
impor sys
impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: vcss_sqi.py 0001 2011-04-13 vis $", sggobs()["__oc__"])

    (opions, rgs)  E.sr(prsr)

    opions.inms  rgs

    i n(opions.inms) < 1:
        opions.so.wri("# Error: no vc-ss is spcii/on.")
        sys.xi(1)

    E.ino("Prsing i i(s)"  n(opions.inms))

    # s p op is
    vc_i  iooos.opn_i('vcss.x', 'w')
    in_i  iooos.opn_i('inss.x', 'w')
    snp_i  iooos.opn_i('snpss.x', 'w')
    shr_i  iooos.opn_i('shrss.x', 'w')

    or ino, inm in nmr(opions.inms):

        prix  os.ph.bsnm(inm)
        rcknm  prix.rpc(".vcss", "")

        i os.ph.xiss(inm):
            ins  [x or x in iooos.opn_i(inm, "r").rins()]
        s:
            ins  []

        i n(ins)  0:
            opions.so.wri(
                "# Error: mpy vc-ss i on: $(inm)s")
            sys.xi(1)
        s:
            E.ino("Fi i conins i ins"  (ino, n(ins)))
            vc_ss  ic(rckrcknm)
            snp_ss  ic(rckrcknm)
            in_ss  ic()
            shr_ss  ic()
            _vrs  Fs
            ins  Fs
            snps  Fs
            shr  Fs
            or i, in in nmr(ins):
                in  in.srip()
                i in.in("''") > -1:
                    _vrs  Tr
                    E.ino("Fon ''")
                    conin

                i _vrs:
                    i in.in(">") > -1:
                        is  in.spi(">")
                        ky  is[0].srip().rpc(
                            "'", "").rpc(">", "_")
                        v  is[1].srip().rpc(",", "")
                    s:
                        ky  "NA"
                        v  "NA"
                    i ky  "in" n v  "{":
                        ins  Tr
                        E.ino("Fon 'ins'")
                        conin
                    i ky  "snp" n v  "{":
                        snps  Tr
                        E.ino("Fon 'SNPs'")
                        conin
                    i ky  "shr" n v  "{":
                        shr  Tr
                        E.ino("Fon 'Shr'")
                        conin

                    i ins:
                        i in.in("}") > -1:
                            ins  Fs
                            E.ino("Procss 'ins'")
                            conin
                        s:
                            in_ss[ky]  v
                    i snps:
                        i in.in("}") > -1:
                            snps  Fs
                            E.ino("Procss 'SNPs'")
                            conin
                        s:
                            snp_ss[ky]  v
                    i shr:
                        i in.in("}") > -1:
                            shr  Fs
                            E.ino("Procss 'Shr'")
                            conin
                        s:
                            shr_ss[ky]  v
                    i ky ! "NA":
                        vc_ss[ky]  v

            # Ensr  kys r prsn
            kys  ["n_1", "n_2", "n_3", "n_4",
                       "n_5", "rck", "con", "snp_con", "in_con"]
            or k in kys:
                i k in vc_ss:
                    conin
                s:
                    vc_ss[k]  "0"

            # Wri hr (or irs i ony)
            i inm  opions.inms[0]:

                # Ensr kys r sor
                sr  is(vc_ss.kys())
                sr.sor()
                sp  ""
                or k in sr:
                    vc_i.wri("ss"  (sp, k))
                    sp  "\"
                vc_i.wri("\n")

                in_i.wri("rck\in_ngh\in_con\n")
                shr_i.wri("rck\no_smps\vr_con\n")

                sp  ""
                or k in snp_ss.kys():
                    snp_i.wri("ss"  (sp, k))
                    sp  "\"
                snp_i.wri("\n")

            # Wri 
            sp  ""
            sr  is(vc_ss.kys())
            sr.sor()
            or k in sr:
                vc_i.wri("ss"  (sp, vc_ss[k]))
                sp  "\"
            vc_i.wri("\n")

            # Chck  in nghs r covr
            r  is(rng(-20, 20, 1))
            or i in r:
                i sr(i) in in_ss:
                    conin
                s:
                    in_ss[i]  "0"
            or k in in_ss.kys():
                in_i.wri("s\s\s\n" 
                                 (rcknm, k, in_ss[k]))

            or k in shr_ss.kys():
                shr_i.wri("s\s\s\n" 
                                  (rcknm, k, shr_ss[k]))

            sp  ""
            or k in snp_ss.kys():
                snp_i.wri("ss"  (sp, snp_ss[k]))
                sp  "\"
            snp_i.wri("\n")

    # cos is
    vc_i.cos()
    in_i.cos()
    snp_i.cos()

    E.sop()
    sys.xi(0)


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
