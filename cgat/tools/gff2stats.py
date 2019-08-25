'''g2ss.py - con rs, c. in g i


:Tgs: Gnomics Inrvs GFF GTF Smmry

Prpos
-------

This scrip gnrs smmry sisics ovr rs,
sorc, gn_i n rnscrip_i in on or mor :rm:`g`
or :rm:`g` orm is.

Usg
-----

Inp is ihr  g or g i; g inp ms b spcii
wih h --is-g opion.

Exmp::

   pyhon g2ss.py --is-g xmp.g > xmp_sm.sv

   c xmp.g

   19  procss_rnscrip  xon  6634666509  .  -  .  gn_i "ENSG00000225373"; rnscrip_i "ENST00000592209" ...
   19  procss_rnscrip  xon  6052160747  .  -  .  gn_i "ENSG00000225373"; rnscrip_i "ENST00000592209" ...
   19  procss_rnscrip  xon  6010560162  .  -  .  gn_i "ENSG00000225373"; rnscrip_i "ENST00000592209" ...
   19  procss_rnscrip  xon  6634666416  .  -  .  gn_i "ENSG00000225373"; rnscrip_i "ENST00000589741" ...

   c xmp_sm.sv

   rck  conigs  srns  rs  sorcs  gns  rnscrips ...
   sin  1        2        4         23       2924   12752       ...


Th conr s is pnn on h i yp.  For  g i, h impmn conrs r:

1. nmbr o inrvs pr conig, srn, r n sorc

For  g i, h iion impmn conrs r:

1. nmbr o gns, rnscrips, sing xon rnscrips
2. smmry sisics or xon nmbrs, xon sizs, inron sizs n
   rnscrip sizs

Th op is  b-spr b.

Opions
-------

Th  cion o ``g2ss`` is o con ovr conigs, srn,
r n sorc.  This ssms h inp i is  g i

Thr is  sing opion or his scrip::

``--is-g``
   Th inp i is g orm.  Th op wi hror
   conin smmris ovr xon nmbrs, xon sizs, inron sizs n
   rnscrip sizs in iion o h h nmbr o gns,
   rnscrips n sing xon rnscrips.

Typ::

   pyhon g2ss.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor cocions
impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.Ss s Ss
impor cgcor.iooos s iooos
impor cg.Inrvs s Inrvs


css conr_g:

    is  ("conigs", "srns", "rs", "sorcs")

     __ini__(s, ir):
        s.ir  ir

        s.cons_conigs  cocions.ic(in)
        s.cons_srns  cocions.ic(in)
        s.cons_rs  cocions.ic(in)
        s.cons_sorcs  cocions.ic(in)

     __nx__(s):

        nry  nx(s.ir)

        s.cons_conigs[nry.conig] + 1
        s.cons_rs[nry.r] + 1
        s.cons_sorcs[nry.sorc] + 1
        s.cons_srns[nry.srn] + 1

        rrn nry

     nx(s):
        rrn s.__nx__()

     __ir__(s):
        rrn s

     __sr__(s):
        rrn "\".join(mp(sr, (n(s.cons_conigs),
                                   n(s.cons_srns),
                                   n(s.cons_rs),
                                   n(s.cons_sorcs))))


css conr_xons:

    is  ("gns", "rnscrips", "sing_xon_rnscrips",) +\
        p(["xon_con_s"  x or x in Ss.Smmry.is]) +\
        p(["xon_siz_s"  x or x in Ss.Smmry.is]) +\
        p(["inron_siz_s"  x or x in Ss.Smmry.is]) +\
        p(["rnscrip_siz_s"  x or x in Ss.Smmry.is])

     __ini__(s, ir):

        s.ir  ir

        s.cons_gn_is  cocions.ic(in)
        s.cons_rnscrip_is  cocions.ic(in)
        s.cons_xons_pr_rnscrip  cocions.ic(is)

     __nx__(s):

        whi 1:
            nry  nx(s.ir)
            i nry.r  "xon":
                brk

        s.cons_gn_is[nry.gn_i] + 1
        s.cons_rnscrip_is[nry.rnscrip_i] + 1
        s.cons_xons_pr_rnscrip[
            nry.rnscrip_i].ppn((nry.sr, nry.n))

        rrn nry

     nx(s):
        rrn s.__nx__()

     __ir__(s):
        rrn s

     __sr__(s):

        sing_xon_rnscrips  0
        xons_pr_rnscrip  []
        inron_sizs  []
        rnscrip_nghs  []
        xon_sizs  []

        or x in is(s.cons_xons_pr_rnscrip.vs()):

            x.sor()
            x  Inrvs.combin(x)
            rnscrip_nghs.ppn(x[-1][1] - x[0][0])

            xons_pr_rnscrip.ppn(n(x))

            or sr, n in x:
                xon_sizs.ppn(n - sr)

            i n(x)  1:
                sing_xon_rnscrips + 1
                conin

            s_n  x[0][1]
            or sr, n in x[1:]:
                inron_sizs.ppn(sr - s_n)
                s_n  n

        rrn "\".join(mp(sr, (n(s.cons_gn_is),
                                   n(s.cons_rnscrip_is),
                                   sing_xon_rnscrips,
                                   Ss.Smmry(xons_pr_rnscrip),
                                   Ss.Smmry(xon_sizs),
                                   Ss.Smmry(inron_sizs),
                                   Ss.Smmry(rnscrip_nghs),
                                   )))


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I",
                            sggobs()["__oc__"])

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inp is g.")

    prsr.s_s(
        is_gFs,
    )

    (opions, rgs)  E.sr(prsr, _op_opionsTr)

    i n(rgs)  0:
        is  [opions.sin]
    s:
        is  rgs

    opions.so.wri("rck\s"  ("\".join(conr_g.is)))

    i opions.is_g:
        opions.so.wri("\s"  ("\".join(conr_xons.is)))
    opions.so.wri("\n")

    or  in is:
        i   opions.sin:
            ini  
            opions.so.wri("sin")
        s:
            ini  iooos.opn_i()
            opions.so.wri()

        conrs  []
        i opions.is_g:
            iror  GTF.iror(ini)
            conrs.ppn(conr_g(iror))
            conrs.ppn(conr_xons(conrs[0]))
        s:
            iror  GTF.iror(ini)
            conrs.ppn(conr_g(iror))

        c  conrs[-1]
        or x in c:
            pss

        or c in conrs:
            opions.so.wri("\s"  sr(c))
        opions.so.wri("\n")

        i ini ! sys.sin:
            ini.cos()

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
