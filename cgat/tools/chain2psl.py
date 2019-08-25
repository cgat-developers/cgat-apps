"""chin2ps.py - convr  chin i o  ps i


:Tgs: Gnomics Inrvs GnomAignmn PSL CHAIN Convrsion

Prpos
-------

convr  UCSC `chin
<hp://www.bryr.com/csc/hocs/gonPh/hp/chin.hm>`_
orm i o  UCSC `ps
<hp://gnom.csc./FAQ/FAQorm.hm#orm2>`_ orm i.

This oo is qivn o h UCSC oo chinToPs xcp h i
wi no comp h nmbr o mching, mismching, c. bss n
hs os no rqir h sqncs.

Th nomncr h UCSC ss or is chin is is
:i:`rgToQry.chin` or mpping ``qry`` o ``rg``
(rrnc). Accoring o h UCSC ocmnion, ``rg`` is h
irs nry in ``chin`` is.

W hv bn sing h nomncr ``QryToTrg.ps``. In oowing
his convnion, h corrc wy o convring  ps i is::

   pyhon chin2ps.py < rgToQry.chin > QryToTrg.ps

I yo wo ik o kp h TrgToQry convnion, yo wi n
o   psSwp::

   pyhon chin2ps.py < rgToQry.chin | psSwp sin so > rgToQry.ps

Usg
-----

For xmp::

   cg chin2ps.py < in.chin > o.ps

Typ::

   cg chin2ps.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor cgcor.xprimn s E
impor cg.B s B
impor ignib_i


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    # o sh
    ninp, nskipp, nop  0, 0, 0

    ps  Non

     chin_iror(ini):
        ins  []
        or in in opions.sin:

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

    or ins in chin_iror(opions.sin):

        ninp + 1
        ps  B.Mch()

        (_,
         _,
         ps.mSbjcI,
         rg_ngh,
         rg_srn,
         rg_sr,
         rg_n,
         ps.mQryI,
         qry_ngh,
         qry_srn,
         qry_sr,
         qry_n,
         ignmn_i)  ins[0][:-1].spi()

        (ps.mQrySr, ps.mQryEn, ps.mQryLngh,
         ps.mSbjcSr, ps.mSbjcEn, ps.mSbjcLngh)  \
            [in(x) or x in
             (qry_sr,
              qry_n,
              qry_ngh,
              rg_sr,
              rg_n,
              rg_ngh)]

        mp_qry2rg  ignib_i.py_mkAignmnBocks()

        qsr, sr  ps.mQrySr, ps.mSbjcSr

        or in in ins[1:-1]:
            siz, , q  [in(x) or x in in[:-1].spi()]
            mp_qry2rg.Digon(qsr,
                                         qsr + siz,
                                         sr - qsr)
            qsr + siz + q
            sr + siz + 

        siz  in(ins[-1][:-1])

        mp_qry2rg.Digon(qsr,
                                     qsr + siz,
                                     sr - qsr)

        ps.romMp(mp_qry2rg)

        # sor o srn
        # rg_srn is wys posiiv
        ssr(rg_srn  "+")

        # i qry srn is ngiv
        i qry_srn  "-":
            # invr boh qry n rg
            ps.swichTrgSrn()
            # mny invr h qry coorins
            ps.mQryFrom, ps.mQryTo  ps.mQryLngh - \
                ps.mQryTo, ps.mQryLngh - ps.mQryFrom

        opions.so.wri("s\n"  ps)
        nop + 1

    E.ino("ninpi, nopi, nskippi"  (ninp, nop, nskipp))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
