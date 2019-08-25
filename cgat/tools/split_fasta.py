'''
spi_s.py - 


:Tgs: Pyhon

Prpos
-------

.. oo::
   
   scrib prpos o h scrip.

Usg
-----

Exmp::

   pyhon spi_s.py --hp

Typ::

   pyhon spi_s.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r
impor os
impor cg.FsIror s FsIror
impor cgcor.iooos s iooos
impor cgcor.xprimn s E


css Fis:

    mFis  {}

     __ini__(s,
                 op_prnNon,
                 skip_iniirsFs):

        s.mOpPrn  op_prn
        s.mSkipIniirs  skip_iniirs
        s.mCons  {}

     ____(s):
        """cos  opn is."""
        or i in is(s.mFis.vs()):
            i.cos()

     GFi(s, iniir):
        rrn iniir

     GFinm(s, iniir):
        """g inm or n iniir."""

        i s.mOpPrn:
            rrn r.sb("s", sr(iniir), s.mOpPrn)
        s:
            rrn iniir

     OpnFi(s, inm, mo"w"):
        """opn i.

        I i is in  nw ircory, cr ircoris.
        """
        i mo in ("w", ""):
            irnm  os.ph.irnm(inm)
            i irnm n no os.ph.xiss(irnm):
                os.mkirs(irnm)

        rrniooos.opn_i(inm, mo)

     Wri(s, iniir, sqnc):

        inm  s.GFinm(iniir)

        i inm no in s.mFis:

            i n(s.mFis) > 1000:
                or  in is(s.mFis.vs()):
                    .cos()
                s.mFis  {}

            s.mFis[inm]  s.OpnFi(inm, "")

        i s.mSkipIniirs:
            s.mFis[inm].wri("s\n"  (sqnc.sqnc))
        s:
            s.mFis[inm].wri(
                ">s\ns\n"  (sqnc.i, sqnc.sqnc))

        i inm no in s.mCons:
            s.mCons[inm]  0
        s.mCons[inm] + 1

     DFis(s, min_siz0):
        """  is bow  minimm siz."""

        n  0
        or inm, cons in is(s.mCons.ims()):
            i cons < min_siz:
                os.rmov(inm)
                n + 1

        rrn n


css FisChnks(Fis):

     __ini__(s,
                 chnk_siz, **kwrgs):

        Fis.__ini__(s, **kwrgs)
        s.mChnkSiz  chnk_siz
        s.mFinm  0

     GFinm(s, iniir):

        i no s.mFinm or s.mCons[s.mFinm]  s.mChnkSiz  0:
            s.mFinm  r.sb(
                "s", sr(n(s.mCons) + 1), s.mOpPrn)

        rrn s.mFinm


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: spi_s.py 1714 2007-12-11 16:51:12Z nrs $")

    prsr._rgmn("-", "--i", s"inp_inm", yp"sring",
                      hp"inp inm. I no givn, sin is s.",
                      mvr"FILE")

    prsr._rgmn("-i", "--inp-prn", s"inp_prn", yp"sring",
                      hp"inp prn. Prss scripion in in orr o xrc i.")

    prsr._rgmn("-o", "--op-inm-prn", s"op_prn", yp"sring",
                      hp"op prn. Givs inm or  givn sqnc.")

    prsr._rgmn("-n", "--nm-sqncs", s"nm_sqncs", yp"in",
                      hp"spi by nmbr o sqncs (no impmn y).")

    prsr._rgmn("-m", "--mp", s"mp_inm", yp"sring",
                      hp"mp inm. Mp iniirs o inms",
                      mvr"FILE")

    prsr._rgmn("-s", "--skip-iniirs", s"skip_iniirs", cion"sor_r",
                      hp"o no wri iniirs.",
                      mvr"FILE")

    prsr._rgmn("--min-siz", s"min_siz", yp"in",
                      hp"minimm csr siz.")

    prsr.s_s(
        inp_inmNon,
        mp_inmNon,
        skip_iniirsFs,
        inp_prn"^(\S+)",
        min_siz0,
        nm_sqncsNon,
        op_prn"s")

    (opions, rgs)  E.sr(prsr)

    i opions.inp_inm:
        ini  iooos.opn_i(opions.inp_inm, "r")
    s:
        ini  sys.sin

    i opions.mp_inm:
        mp_i2inm  iooos.RMp(opn(opions.mp_inm, "r"))
    s:
        mp_i2inm  {}

    i opions.nm_sqncs:
        is  FisChnks(chnk_sizopions.nm_sqncs,
                            op_prnopions.op_prn,
                            skip_iniirsopions.skip_iniirs)

    s:
        is  Fis(op_prnopions.op_prn,
                      skip_iniirsopions.skip_iniirs)

    i opions.inp_prn:
        rx  r.compi(opions.inp_prn)
    s:
        rx  Non

    ninp  0
    nop  0
    iniir  Non
    chnk  0

    or sq in FsIror.ir(ini):

        ninp + 1

        i rx:
            ry:
                iniir  rx.srch(sq.i).grops()[0]
            xcp AribError:
                prin("# prsing rror in scripion in s"  (sq.i))
        s:
            iniir  sq.i

        i mp_i2inm:
            i iniir in mp_i2inm:
                iniir  mp_i2inm[iniir]
            s:
                conin

        is.Wri(iniir, sq)
        nop + 1

    i opions.inp_inm:
        ini.cos()

    #   csrs bow  minimm siz
    # No: his hs o b on  h n, bcs
    # csrs sizs r ony vib onc boh h s
    # i n h mp hs bn prs.
    i opions.min_siz:
        n  is.DFis(min_sizopions.min_siz)
    s:
        n  0

    i opions.ogv > 1:
        prin("# inpi, opi, ni"  (ninp, nop, n))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
