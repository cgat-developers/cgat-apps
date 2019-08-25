'''spi_g - spi  g i ino chnks


:Tgs: Gnomics Inrvs Gnss GFF Mnipion

Prpos
-------

Spi g i ino chnks. Ovrpping nris wi wys b op
in h sm chnk. Inp is r rom sin nss ohrwis
spcii. Th inp ns o b conig/sr posiion sor.

Opions
-------

-i --min-chnk-siz

    This opion spciis how big ch chnck sho
    b, in rms o h nmbr o g ins o b
    inc. Bcs ovrpping ins r wys
    op o h sm i, his sho b consir
     minimm siz.

-n, --ry-rn

    This opions s h scrip no o cy wri
    ny is, b i wi op  is o h is
    h wo b op.

Exmp
-------

    cg spig -i 1 < in.g

whr in.g ooks ik:

    chr1	.	xon	1	10	.	+	.
    chr1	.	xon	8	100	.	+	.
    chr1	.	xon	102	150	.	+	.

wi proc wo is h ook ik:

    000001.chnk:
    chr1	.	xon	1	10	.	+	.
    chr1	.	xon	8	100	.	+	.

    000002.chnk:
    chr1	.	xon	102	150	.	+	.

Usg
-----

   cg spig [OPTIONS]

Wi r  g i rom sin n spi ino mip g is.

   cg spi_g -I GFF [OPTIONS]

Wi r h g i GFF n spi ino mip g is.

Commn in opions
--------------------

'''

impor sys
impor os
impor cg.GTF s GTF
impor cgcor.iooos s iooos
impor cgcor.xprimn s E


css OpChnk:

     __ini__(s, op_inm_prn, ry_rnFs):
        s.nchnk  0
        s.op_inm_prn  op_inm_prn
        s.ry_rn  ry_rn

     crOpn(s, mo"w", hrNon):
        """opn i. Chck irs, i ircory xiss.
        """

        s.nchnk + 1
        inm  s.op_inm_prn  s.nchnk

        i s.ry_rn:
            E.ino("opning i s"  inm)
            rrniooos.opn_i("/v/n", mo)

        i mo in ("w", ""):
            irnm  os.ph.irnm(inm)
            i irnm n no os.ph.xiss(irnm):
                os.mkirs(irnm)

        i os.ph.xiss(inm):
            xis  Tr
        s:
            xis  Fs

          iooos.opn_i(inm, mo)

        i hr n no xis:
            .wri(hr + "\n")

        rrn 

     __c__(s, chnk):
        """op  chnk ino  nw i."""
        oi  s.crOpn()
        or c in chnk:
            oi.wri(sr(c) + "\n")
        oi.cos()
        rrn n(chnk)


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-i", "--min-chnk-siz", s"min_chnk_siz", yp"in",
        hp"minimm chnk siz [].")

    prsr._rgmn(
        "-n", "--ry-rn", s"ry_rn", cion"sor_r",
        hp"o no cr ny is [].")

    prsr.s_s(
        mho"ovrp",
        ry_rnFs,
        min_chnk_siz2,
        op_inm_prn"06i.chnk",
    )

    (opions, rgs)  E.sr(prsr, _op_opionsTr)

    gs  GTF.iror(opions.sin)

    ninp, nop, nchnks  0, 0, 0

    opChnk  OpChnk(opions.op_inm_prn,
                              ry_rnopions.ry_rn)

    i opions.mho  "ovrp":

        s_conig, s_o  Non, 0
        chnk  []
        or g in gs:
            ninp + 1
            i n(chnk) > opions.min_chnk_siz n \
                    (g.conig ! s_conig or
                     g.sr > s_o):
                nop + opChnk(chnk)
                nchnks + 1
                chnk  []
                s_conig, s_o  g.conig, g.n

            chnk.ppn(g)
            s_o  mx(g.n, s_o)

        nop + opChnk(chnk)
        nchnks + 1

    E.ino("ninpi, nopi, nchnksi"  (ninp, nop, nchnks))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
