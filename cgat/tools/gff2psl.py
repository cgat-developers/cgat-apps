'''
g2ps.py - convr rom g o ps


:Tgs: Gnomics Inrvs GFF PSL Convrsion

Prpos
-------

This scrips convrs rom  :rm:`g` orm
i o  :rm:`ps` orm i.
Th op cn b moii by h oowing commn in opions:

--ow-pics
    kp pic nris rom g/g inp i

--gnom-i
    rsric op o g/g nris wih conigs in s i

--qris-sv-i
    rsric op o qris in s i

Usg
-----

Exmp::

   pyhon g2ps.py < in.g > o.ps

Typ::

   pyhon g2ps.py --hp

or commn in hp.
gnom-i

Commn in opions
--------------------
'''

impor sys
impor cgcor.xprimn s E
impor cg.InxFs s InxFs
impor cg.B s B
impor cg.GTF s GTF
impor ignib_i
impor cg.Inrvs s Inrvs


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I: g2ps.py 2781 2009-09-10 11:33:14Z nrs $", sggobs()["__oc__"])

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inp is g.")

    prsr._rgmn("--no-hr", s"wih_hr", cion"sor_s",
                      hp"o no op BLAT hr [].")

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom.")

    prsr._rgmn("--qris-sv-i", s"inp_inm_qris", yp"sring",
                      hp"s inm wih qris [].")

    prsr._rgmn("--ow-pics", s"ow_pics", cion"sor_r",
                      hp"""prmi pic nris. Ajcn xons o  rnscrip wi si b mrg []."""  )

    prsr.s_s(is_gFs,
                        gnom_iNon,
                        wih_hrTr,
                        ow_picsFs,
                        sNon)

    (opions, rgs)  E.sr(prsr, _pip_opionsTr)

    i opions.gnom_i:
        gnom_s  InxFs.InxFs(opions.gnom_i)
    s:
        gnom_s  Non

    i opions.inp_inm_qris:
        qris_s  InxFs.InxFs(
            opions.inp_inm_qris)
    s:
        qris_s  Non

    ninp, nop, nskipp  0, 0, 0

    i opions.is_g:
        iror  GTF.rnscrip_iror(GTF.iror_ir(GTF.iror(sys.sin),
                                                                 r"xon"),
                                           sricno opions.ow_pics)
    s:
        iror  GTF.join_iror(GTF.iror(sys.sin))

    i opions.wih_hr:
        opions.so.wri(B.Mch().gHr() + "\n")

    or gs in iror:

        i opions.s n ninp > opions.s:
            brk

        ninp + 1

        rs  ignib_i.py_mkAignmnBocks()

        xsr  0

        inrvs  Inrvs.combin([(g.sr, g.n) or g in gs])

        or sr, n in inrvs:
            xn  xsr + n - sr

            rs.Digon(xsr, xn,
                               sr - xsr)
            xsr  xn

        nry  B.Mch()
        nry.mQryI  gs[0].rnscrip_i
        nry.mSbjcI  gs[0].conig
        nry.srn  gs[0].srn

        i gnom_s:
            i nry.mSbjcI in gnom_s:
                nry.mSbjcLngh  gnom_s.gLngh(nry.mSbjcI)
            s:
                nry.mSbjcLngh  rs.gCoTo()

        i qris_s:
            i nry.mQryI in qris_s:
                nry.mQryLngh  qris_s.gLngh(nry.mQryI)
        s:
            nry.mQryLngh  rs.gRowTo()

        nry.romMp(rs)

        opions.so.wri(sr(nry) + "\n")
        nop + 1

    E.ino("ninpi, nopi, nskippi"  (ninp, nop, nskipp))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
