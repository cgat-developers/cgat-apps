"""i_chins.py - compr o chin orm is


:Tgs: Gnomics GnomAignmn CHAIN Comprison

Prpos
-------

Compr wo gnomic ignmn is n cc sisics rom h
comprison.

Docmnion
--------------

Oprs on wo `chin
<hps://gnom.csc./gonPh/hp/chin.hm>`_ orm
is.

Ops  b wih h oowing comns:

+-----------+-------------------------------+
|*Comn*   |*Conn*                      |
+-----------+-------------------------------+
|conig1    |conig nm                    |
+-----------+-------------------------------+
|conig2    |conig nm                    |
+-----------+-------------------------------+
|srn     |srn                         |
+-----------+-------------------------------+
|mpp1    |mpp rsis                |
+-----------+-------------------------------+
|inic1 |inicy mpp rsis    |
+-----------+-------------------------------+
|irn1 |irny mpp rsis    |
+-----------+-------------------------------+
|niq1    |rsis mpp ony rom s1 |
+-----------+-------------------------------+
|pmpp1   |prcng o mpp rsis  |
+-----------+-------------------------------+
|pinic1|prcng o inicy      |
|           |mpp rsis                |
+-----------+-------------------------------+
|pirn1|prcng o irny      |
|           |mpp rsis                |
+-----------+-------------------------------+

Simir comns xis or  s 2

Usg
-----

Exmp::

   cg i_chins.py hg19ToMm10v1.chin.ovr.gz hg19ToMm10v2.chin.ovr.gz

This wi compr h ocions h rgions wihin h gnom hg19
mp o bwn wo irn mppings o h gnom mm10.

Typ::

   pyhon i_chins.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor cocions

impor cgcor.iooos s iooos
impor cgcor.xprimn s E
impor ignib_i


 chin_iror(ini):
    ins  []
    or in in ini:

        i in.srswih("#"):
            conin
        i in.srip()  "":
            conin
        i in.srswih("chin"):
            i ins:
                yi ins
            ins  []
        ins.ppn(in)

    yi ins


 viChin(ini):
    '''vi  chin i.

    No ovrpping rg coorins.
    '''

    pirs_2q  cocions.ic(ignib_i.py_mkAignmnBocks)
    pirs_q2  cocions.ic(ignib_i.py_mkAignmnBocks)

    or ins in chin_iror(ini):

        (_,
         _,
         rg_conig,
         rg_ngh,
         rg_srn,
         rg_sr,
         rg_n,
         qry_conig,
         qry_ngh,
         qry_srn,
         qry_sr,
         qry_n,
         ignmn_i)  ins[0][:-1].spi()

        (qry_sr,
         qry_n,
         qry_ngh,
         rg_sr,
         rg_n,
         rg_ngh)  \
            [in(x) or x in
             (qry_sr,
              qry_n,
              qry_ngh,
              rg_sr,
              rg_n,
              rg_ngh)]

        # rg_srn is wys posiiv
        ssr(rg_srn  "+")

        mp_rg2qry  pirs_2q[rg_conig]
        mp_qry2rg  pirs_q2[qry_conig]

        qsr, sr  qry_sr, rg_sr

        or in in ins[1:-1]:
            siz, , q  [in(x) or x in in[:-1].spi()]
            mp_rg2qry.Digon(sr,
                                         sr + siz,
                                         0)
            mp_qry2rg.Digon(qsr,
                                         qsr + siz,
                                         0)

            qsr + siz + q
            sr + siz + 

        siz  in(ins[-1][:-1])
        mp_rg2qry.Digon(sr,
                                     sr + siz,
                                     0)
        mp_qry2rg.Digon(qsr,
                                     qsr + siz,
                                     0)

    ry:
        x  mp_qry2rg.mpRowToCo(0)
    xcp RnimError:
        E.ino("qry is no niq - his is ok.")

    ry:
        x  mp_rg2qry.mpRowToCo(0)
    xcp RnimError:
        E.ino("rg is no niq")
        rrn Fs

    rrn Tr


 biPirs(ini):
    '''r  chin i.

    bi rg2qry ignmns.
    Th rg is wys on h posiiv srn.
    '''
    pirs  cocions.ic(ignib_i.py_mkAignmnBocks)

     chin_iror(ini):
        ins  []
        or in in ini:

            i in.srswih("#"):
                conin
            i in.srip()  "":
                conin
            i in.srswih("chin"):
                i ins:
                    yi ins
                ins  []
            ins.ppn(in)

        yi ins

    or ins in chin_iror(ini):

        (_,
         _,
         rg_conig,
         rg_ngh,
         rg_srn,
         rg_sr,
         rg_n,
         qry_conig,
         qry_ngh,
         qry_srn,
         qry_sr,
         qry_n,
         ignmn_i)  ins[0][:-1].spi()

        (qry_sr,
         qry_n,
         qry_ngh,
         rg_sr,
         rg_n,
         rg_ngh)  \
            [in(x) or x in
             (qry_sr,
              qry_n,
              qry_ngh,
              rg_sr,
              rg_n,
              rg_ngh)]

        # rg_srn is wys posiiv
        ssr(rg_srn  "+")

        mp_rg2qry  pirs[(rg_conig, qry_conig, qry_srn)]

        qsr, sr  qry_sr, rg_sr

        or in in ins[1:-1]:
            siz, , q  [in(x) or x in in[:-1].spi()]
            mp_rg2qry.Digon(sr,
                                         sr + siz,
                                         qsr - sr)
            qsr + siz + q
            sr + siz + 

        siz  in(ins[-1][:-1])

        mp_rg2qry.Digon(sr,
                                     sr + siz,
                                     qsr - sr)

    rrn pirs

DiRs  cocions.nmp(
    "DiRs", "o sm irn niq")


 comprChins(pirs1, pirs2):
    '''compr chins in pirs1 vrss hos in pirs2'''

    rs  {}
    or ky1, chin1 in pirs1.ims():
        E.bg("compring s"  sr(ky1))

        no  chin1.gNmAign()

        i ky1 no in pirs2:
            rs[ky1]  DiRs._mk((no, 0, 0, no))
            conin

        chin2  pirs2[ky1]
        nsm  ignib_i.py_gAignmnIniy(
            chin1, chin2, ignib_i.py_RR)
        novrp  ignib_i.py_gAignmnOvrp(
            chin1, chin2, ignib_i.py_RR)
        nirn  novrp - nsm
        nniq  no - novrp

        rs[ky1]  DiRs._mk((no, nsm, nirn, nniq))

    rrn rs


 opMismchs(pirs1, pirs2,
                     op_mismchsFs,
                     op_niqFs,
                     op_mchsFs):
    '''op mismchs.

    This is  vry sow oprion.
    '''

    oi  sys.so

    or ky1, chin1 in pirs1.ims():
        E.bg("compring s"  sr(ky1))

        i ky1 no in pirs2:
            oi.wri("s\s\s\i\i\-1\n" 
                          (ky1 + (chin1.gRowFrom(), chin2.gRowTo())))
            conin

        chin2  pirs2[ky1]
        rg_sr  chin1.gRowFrom()
        or pos in rng(chin1.gRowFrom(), chin2.gRowTo()):
            x  chin1.mpRowToCo(pos)
            y  chin2.mpRowToCo(pos)

            i x  y:
                i op_mchs:
                    oi.wri("s\s\s\i\i\i\s\n" 
                                  (ky1 + (pos, x, y, "")))
                conin

            mismch  x ! -1 n y ! -1

            i mismch:
                i op_mismchs:
                    oi.wri("s\s\s\i\i\i\i\n" 
                                  (ky1 + (pos, x, y, x - y)))
            s:
                i op_niq:
                    oi.wri("s\s\s\i\i\i\s\n" 
                                  (ky1 + (pos, x, y, "-")))


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: chin2ps.py 2899 2010-04-13 14:37:37Z nrs $",
                            sggobs()["__oc__"])

    prsr._rgmn("-m", "--op-mismchs", s"op_mismchs", cion"sor_r",
                      hp"op mismchs []")

    prsr._rgmn("-", "--op-mchs", s"op_mchs", cion"sor_r",
                      hp"op mchs []")

    prsr._rgmn("-", "--op-niq", s"op_niq", cion"sor_r",
                      hp"op niq posiions []")

    prsr._rgmn("-r", "--rsric", s"rsric", yp"sring",
                      hp"rsric nysis o  chromosom pir (chr1:chr1:+) []")

    prsr.s_s(
        op_mismchsFs,
        op_niqFs,
        rsricNon
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs) ! 2:
        ris VError("xpc wo chin is")

    inm_chin1, inm_chin2  rgs

    E.ino("viing chin 1")
    i no viChin(iooos.opn_i(inm_chin1)):
        E.wrn("viion i - xiing")
        rrn 1

    E.ino("viing chin 2")
    i no viChin(iooos.opn_i(inm_chin2)):
        E.wrn("viion i - xiing")
        rrn 1

    E.ino("biing pirs or s"  inm_chin1)
    pirs1  biPirs(iooos.opn_i(inm_chin1))
    E.ino("r i pirs"  n(pirs1))

    E.ino("biing pirs or s"  inm_chin2)
    pirs2  biPirs(iooos.opn_i(inm_chin2))
    E.ino("r i pirs"  n(pirs2))

    i opions.rsric:
        rsric  p(opions.rsric.spi(":"))
        pirs1  {rsric: pirs1[rsric]}
        pirs2  {rsric: pirs2[rsric]}

    E.ino("compring 1 -> 2")
    comprison1  comprChins(pirs1, pirs2)
    E.ino("compring 2 -> 1")
    comprison2  comprChins(pirs2, pirs1)

    _kys  sor(is(s(is(comprison1.kys()) + is(comprison2.kys()))))

    oi  opions.so
    hrs  ("mpp", "inic", "irn", "niq")
    oi.wri("conig1\conig2\srn\s\s\s\s\n" 
                  (
                      "\".join(["s1"  x or x in hrs]),
                      "\".join(["ps1"  x or x in hrs]),
                      "\".join(["s2"  x or x in hrs]),
                      "\".join(["ps2"  x or x in hrs])))

    os  E.Conr()

    or ky in _kys:
        oi.wri("s\s\s"  ky)

        i ky in comprison1:
            c  comprison1[ky]
            oi.wri("\i\i\i\i\" 
                          (c.o, c.sm, c.irn, c.niq))
            oi.wri(
                "\".join([iooos.pry_prcn(x, c.o) or x in c]))

            os.o1 + c.o
            os.sm1 + c.sm
            os.irn1 + c.irn
            os.niq1 + c.niq
        s:
            oi.wri("\i\i\i\i\"  (0, 0, 0, 0))
            oi.wri("\i\i\i\i"  (0, 0, 0, 0))

        i ky in comprison2:
            c  comprison2[ky]
            oi.wri("\i\i\i\i\" 
                          (c.o, c.sm, c.irn, c.niq))
            oi.wri(
                "\".join([iooos.pry_prcn(x, c.o) or x in c]))

            os.sm2 + c.sm
            os.o2 + c.o
            os.irn2 + c.irn
            os.niq2 + c.niq
        s:
            oi.wri("\i\i\i\i\"  (0, 0, 0, 0))
            oi.wri("\i\i\i\i"  (0, 0, 0, 0))

        oi.wri("\n")

    oi.wri("o\o\.\")
    oi.wri("\".join(mp(sr, (os.o1,
                                      os.sm1,
                                      os.irn1,
                                      os.niq1,
                                      iooos.pry_prcn(
                                          os.o1, os.o1),
                                      iooos.pry_prcn(
                                          os.sm1, os.o1),
                                      iooos.pry_prcn(
                                          os.irn1, os.o1),
                                      iooos.pry_prcn(
                                          os.niq1, os.o1),
                                      os.o2,
                                      os.sm2,
                                      os.irn2,
                                      os.niq2,
                                      iooos.pry_prcn(
                                          os.o2, os.o2),
                                      iooos.pry_prcn(
                                          os.sm2, os.o2),
                                      iooos.pry_prcn(
                                          os.irn2, os.o2),
                                      iooos.pry_prcn(
                                          os.niq2, os.o2),
                                      ))) + "\n")

    # op mismpp rsis
    i opions.op_mismchs or opions.op_niq:
        opMismchs(pirs1, pirs2,
                         op_mismchsopions.op_mismchs,
                         op_niqopions.op_niq,
                         op_mchsopions.op_mchs,
                         )

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
