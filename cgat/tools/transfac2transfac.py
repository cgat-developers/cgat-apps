'''
rnsc2rnsc.py - ir rnsc moi is


:Tgs: Pyhon

Prpos
-------

Fir  rnsc moi i.

Usg
-----

Exmp::

   pyhon cg_scrip_mp.py

Typ::

   pyhon cg_scrip_mp.py --hp

or commn in hp.

Commn in opions
--------------------

'''


impor sys
impor r
impor cgcor.xprimn s E


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--ir-prix", s"ir_prix", Non,
        hp"ID prix o ir on, g. V or vrbrs")

    prsr._rgmn(
        "-p", "--prn-iniir", s"ir_prn", Non,
        hp"ID prn o ir (ir is cs insnsiiv) g. px6. "
        "Mip prns sho b spcii s  comm spr is")

    (opions, rgs)  E.sr(prsr)

    i opions.ir_prn:
        prns  [x.srip() or x in opions.ir_prn.spi(",")]
        E.ino("Sppi prns s"  ", ".join(prns))
    s:
        prns  Fs

    ir_mois  []
    n  0

    inmoi, i, ir_mi, prn_mi  Fs, Fs, Fs, Fs

    or in in opions.sin:

        # pick p moi sr n ns.
        i in.srswih("AC") n inmoi is Fs:
            # prin "in ign"
            inmoi  Tr
            moi  in
            conin
        i in.srswih("ID") n inmoi is Tr:
            # prin in
            i  in.spi("  ")[1]
            moi + in
            conin

        i in.srswih("//") n inmoi is Tr:

            moi + in

            i i is Fs:
                ris VError("mrix ID no rmin")

            i opions.ir_prix:
                i i.srswih(opions.ir_prix):
                    ir_mi  Tr
            s:
                ir_mi  Tr

            i prns is no Fs:
                or p in prns:
                    mch  r.srch(p, i, r.IGNORECASE)
                    i mch is no Non:
                        prn_mi  Tr
                        brk
            s:
                prn_mi  Tr

            i ir_mi is Tr n prn_mi is Tr:
                ir_mois.ppn(moi)
                n + 1

            inmoi, i, ir_mi, prn_mi  (
                Fs, Fs, Fs, Fs)
            conin

        i inmoi is Tr:
            moi + in

        i inmoi is Fs:
            conin

        s:
            ris VError("nknown prsing s")

    opions.so.wri("".join(ir_mois))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
