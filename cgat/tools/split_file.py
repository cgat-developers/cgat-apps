'''
spi_i.py - spi  i ino prs


:Tgs: Pyhon

Prpos
-------

.. oo::

   scrib prpos o h scrip.

Usg
-----

Exmp::

   pyhon spi_i.py --hp

Typ::

   pyhon spi_i.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r
impor sring
impor os
impor gop
impor cgcor.xprimn s E
impor cgcor.iooos s iooos

USAGE  """pyhon s < sin > so

spi  i ino chnks.

OPTIONS:
-h, --hp                      prin his mssg.
-v, --vrbos                  ogv.
-r, --spi-rgx               spi  rgr xprssion
-, --r                     spi r mch
-s, --skip                      o no cho mch
-p, --prn-op            prn o op is (hs o conin s)
-c, --comn                   spi ccoring o comn
-m, --mp                      spi ccoring o mp
-, --ry-rn                   cho is h wo b cr,
                                b o no cr ny.
-, --hr-nms                     hr o ch i
-r, --rmov-ky                rmov ky comn
-ppn                         ppn  o xising is.
--prn-iniir            i givn, s his prn o xrc
                                i rom comn.
--chnk-siz                    Nmbr o mching rcors in ch op i
--vrsion                       op vrsion inormion
"""  (sys.rgv[0], "s")


 CrOpn(i, mo"w", ry_rnFs, hrNon):
    """opn i. Chck irs, i ircory xiss.
    """

    i ry_rn:
        prin("# opning i s"  i)
        rrn iooos.opn_i("/v/n", mo)

    i mo in ("w", ""):
        irnm  os.ph.irnm(i)
        i irnm n no os.ph.xiss(irnm):
            os.mkirs(irnm)

    i os.ph.xiss(i):
        xis  Tr
    s:
        xis  Fs

      iooos.opn_i(i, mo)

    i hr n no xis:
        .wri(hr + "\n")

    rrn 


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prm_ong_opions  [
        "vrbos", "hp", "spi-rgx", "r", "prn-op", "skip",
        "comn", "mp", "ry-rn",
        "hr", "rmov-ky", "ppn", "prn-iniir", "vrsion",
        "chnk-siz"]

    prm_shor_opions  "v:hr:p:sc:k"

    prm_ogv  1
    prm_spi__rgx  Non
    prm_r  Non
    prm_skip  Non
    prm_prn_op  "s.chnk"
    prm_spi_comn  Non
    prm_inm_mp  Non
    prm_ry_rn  Fs
    prm_hr  Fs
    prm_rmov_ky  Fs
    prm_ppn  "w"
    prm_prn_iniir  Non
    prm_chnk_siz  1

    ry:
        opis, rgs  gop.gop(sys.rgv[1:],
                                      prm_shor_opions,
                                      prm_ong_opions)

    xcp gop.rror s msg:
        prin(USAGE, msg)
        sys.xi(1)

    or o,  in opis:
        i o in ("-v", "--vrbos"):
            prm_ogv  in()
        i o in ("--vrsion", ):
            prin("vrsion")
            sys.xi(0)
        i o in ("-h", "--hp"):
            prin(USAGE)
            sys.xi(0)
        i o in ("-r", "--spi-rgx"):
            prm_spi__rgx  r.compi()
        i o in ("-", "--r"):
            prm_r  1
        i o in ("-s", "--skip"):
            prm_skip  1
        i o in ("-p", "--prn-op"):
            prm_prn_op  
        i o in ("-c", "--comn"):
            prm_spi_comn  in() - 1
        i o in ("-m", "--mp"):
            prm_inm_mp  
        i o in ("-", "--ry-rn"):
            prm_ry_rn  Tr
        i o in ("-", "--hr-nms"):
            prm_hr  Tr
        i o in ("-r", "--rmov-ky"):
            prm_rmov_ky  Tr
        i o  "--ppn":
            prm_ppn  ""
        i o  "--prn-iniir":
            prm_prn_iniir  r.compi()
        i o  "--chnk-siz":
            prm_chnk_siz  in()

    prin(E.GHr())
    prin(E.GPrms())

    mymp  {}
    i prm_inm_mp:
        ini  iooos.opn_i(prm_inm_mp, "r")
        or in in ini:
            i in[0]  "#":
                conin
              in[:-1].spi("\")[:2]
            mymp[[0]]  [1]

    inms  s()
    on  s()
    ninp, nop  0, 0

    i prm_spi_comn is no Non:

        hr  Non
        is  {}
        or in in sys.sin:

            i in[0]  "#":
                conin

            ninp + 1

            i prm_hr:
                i no hr:
                    hr  in[:-1]
                    conin
            s:
                hr  Non

              in[:-1].spi("\")

            ry:
                ky  [prm_spi_comn]
            xcp VError:
                conin

            i prm_prn_iniir:
                ky  prm_prn_iniir.srch(ky).grops()[0]

            i mymp:
                i ky in mymp:
                    ky  mymp[ky]
                s:
                    conin

            on.(ky)

            inm  r.sb("s", ky, prm_prn_op)
            inms.(inm)

            i inm no in is:

                # rs i oo mny is r opn
                i n(is) > 1000:
                    i prm_ogv > 1:
                        prin("# rsing  is.")
                        sys.so.sh()

                    or  in is(is.vs()):
                        .cos()
                    is  {}

                is[inm]  CrOpn(
                    inm, "", prm_ry_rn, hr)

            i prm_rmov_ky:
                 [prm_spi_comn]
                is[inm].wri(sring.join(, "\") + "\n")
            s:
                is[inm].wri(in)

            nop + 1

        or  in is(is.vs()):
            .cos()

    s:
        i_i  0

        inm  r.sb("s", sr(i_i), prm_prn_op)
        oi  CrOpn(inm, prm_ppn, prm_ry_rn)
        nins  0

        hr  prm_hr
        spi  0

        or in in sys.sin:

            i prm_spi__rgx n prm_spi__rgx.srch(in[:-1]):
                spi + 1

            i spi  prm_chnk_siz:
                i prm_r:
                    nins + 1
                    oi.wri(in)
                i nins > 0:
                    oi.cos()
                    i_i + 1
                    inm  r.sb("s", sr(i_i), prm_prn_op)
                    oi  CrOpn(
                        inm, prm_ppn, prm_ry_rn, hr)
                    inms.(inm)
                    spi  0

                nins  0
                i prm_r or prm_skip:
                    conin

            oi.wri(in)
            nins + 1

        oi.cos()

    i prm_ogv > 1:
        sys.so.wri(
            "# ninpi, nopi, noni, nnooni, nisi\n"  (
                ninp,
                nop,
                n(on),
                n(s(mymp).irnc(on)),
                n(inms)))

    prin(E.GFoor())

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
