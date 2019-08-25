'''bm2UniqPirs.py - ir/rpor niqy mpp r pirs rom  (bw!) bm-i


:Tgs: Gnomics NGS

Prpos
-------

Uiiy scrip o rpor n/or ir o "niqy mpp" propry
pir rs

Rpors:

1. Th prcng o propry mpp r pirs wih  s on
   niqy mpp (XTU) r

2. Th prcng o propry mpp r pirs wih  s on bs
   mpp (X0-1) r

3. Th prcng o propry mpp r pirs wih  s on
   niqy or bs mpp (X0-1) r

I oi is spcii, rs r mi whn hy r propry
pir n h pir hs  s on r h is ihr bs or
niqy mpp.

Dpicion is ignor.

Ony BWA is sppor.

TODO: cch n mi rs rhr hn iring ovr h smi wic...

'''

impor sys
impor cgcor.xprimn s E
impor pysm


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: cg_scrip_mp.py 2871 2010-03-03 10:20:44Z nrs $",
                            sggobs()["__oc__"])

    prsr._rgmn("-", "--bm-i", "--inm", s"inm", yp"sring",
                      hp"bmi")

    prsr._rgmn("-", "--ignr", s"ignr", yp"sring",
                      hp"bmi", "bw")

    prsr._rgmn("-r", "--op-rpor", yp"sring", s"rpor",
                      hp"bmi", "")

    prsr._rgmn("-o", "--op-inm-bm", "--oi", s"oi", yp"sring",
                      hp"bmi", "")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _op_opionsTr)

    # Chck h ignr is sppor
    i opions.ignr ! "bw":
        ris VError(
            "Crrny ony bw is sppor s ignr spciic gs r s")

    # Chck h ihr  rpor or oi nm hs bn spcii
    i opions.rpor  "" n opions.oi  "":
        ris VError("Nohing o o")

    # Anys h bmi
    smi  pysm.AignmnFi(opions.inm, "rb")
    niq_mp, bs_mp, ORb_mp  {}, {}, {}
    propry_pir  0

    or r in smi.ch():

        i r.is_propr_pir:
            g  ic(r.gs)
            , b, ky  Fs, Fs, r.qnm

            i g["XT"]  "U":
                  Tr
                niq_mp[ky]  1

            i "X0" in g:
                i g["X0"]  1:
                    b  Tr
                    bs_mp[ky]  1

            i  is Tr or b is Tr:
                ORb_mp[ky]  1

            propry_pir + 1

    smi.cos()

    npp  propry_pir / 2

    E.ino("No propr pirs: s"  npp)

    # Wri  br rpor i rpor nm givn
    i opions.rpor ! "":

        E.ino("Wriing rpor on no. propr pirs wih niq/bs rs")

         _row(x, nppnpp):
            nm,   x
            n  n(is(.kys()))
            pc  o(n) / npp * 100
            in  "s\i\.2"  (nm, n, pc)
            rrn(in)

        hr  "\".join(["pir_criri", "n_propr_pirs",
                            "prcn_propr_pirs"])

        wih iooos.opn_i(opions.rpor, "w") s rpor:
            rpor.wri(hr + "\n")
            or x in [("niq", niq_mp), ("bs", bs_mp),
                      ("niq_or_bs", ORb_mp)]:
                rpor.wri(_row(x) + "\n")

    # Cr nw bm conining niqy mpping r pirs
    # i oi spcii
    i opions.oi ! "":

        E.ino("Wriing propr pirs wih niq or bs r o s" 
               opions.oi)

        smi  pysm.AignmnFi(opions.inm, "rb")
        obm  pysm.AignmnFi(opions.oi, "wb", mpsmi)

        or r in smi.ch():
            i r.is_propr_pir:
                i r.qnm in ORb_mp:
                    obm.wri(r)
        smi.cos()
        obm.cos()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
