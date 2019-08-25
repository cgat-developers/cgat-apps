"""bm2b.py - convr bm orm i o b orm i


:Tgs: Gnomics NGS Inrvs BAM BED Convrsion

Prpos
-------

This oo convrs BAM is ino BED is sppying h inrvs
or ch r in h BAM i.  BAM is ms hv  corrsponing
inx i i. xmp.bm n xmp.bm.bi

For xmp::

   smoos viw xmp.bm

   READ1    163    1      13040   15     76M          13183   219     ...
   READ1    83     1      13183   7      76M          13040   -219    ...
   READ2    147    1      13207   0      76M          13120   -163    ...

   pyhon bm2b.py xmp.bm

   1       13039   13115   READ1     15      +
   1       13119   13195   READ2     0       +
   1       13182   13258   READ1     7       -
   1       13206   13282   READ2     0       -

By , bm2b ops ch r s  spr inrv.  Wih
h opion ``--mrg-pirs`` pir-n rs r mrg n op s
 sing inrv. Th srn is s ccoring o h irs r in 
pir.

Usg
-----

::

   cg bm2b BAMFILE [--mrg-pirs] [opions]

oprs on h i BAMFILE::

   cg bm2b [--mrg-pirs] [opions]

oprs on h sin s os::

   cg bm2b -I BAMFILE [--mrg-pirs] [opions]


To mrg pir-n rs n op rgmn inrv i. mos
mpp bs o righmos mpp bs::

   c xmp.bm | cg bm2b --mrg-pirs

   1       13119   13282   READ2     0       +
   1       13039   13258   READ1     7       +

To s mrg pirs on ony  rgion o h gnom s smoos viw::

   smoos viw -b xmp.bm 1:13000:13100 | cg bm2b --mrg-pirs

No h his wi sc rgmns wr h irs r-in-pir is in
h rgion.

Opions
-------

-m, --mrg-pirs
    Op on rgion pr rgmn rhr hn on rgion pr r,
    hs  sing rgion is cr srching rom h sr o h
    ris r in pir o h n o h scon.

    R pirs h m h oowing criri r rmov:

    * Rs whr on o h pir is nmpp
    * Rs h r no pir
    * Rs whr h pirs r mpp o irn chromosoms
    * Rs whr h h insr siz is no bwn h mx n
      min (s bow)

.. wrning::

    Mrg rgmns r wys rrn on h +v srn.
    Frgmn n poin is sim s h ignmn sr posiion
    o h scon-in-pir r + h ngh o h irs-in-pir
    r. This my  o inccrcy i yo hv n inron-wr
    ignr.

--mx-insr-siz, --min-insr-siz
    Th mximm n minimm siz o h insr h is ow whn
    sing h --mrg-pirs opion. R pirs cosr o ghr or hr
    pr hn h min n mx rpscivy r skipp.

-b, --b-orm
    Wh orm o op h rss in. Th irs n comns o h b
    i wi b op.



Typ::

   pyhon bm2b.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor sys
impor pysm
impor cgcor.xprimn s E
rom cg.BmToos.bmoos impor mrg_pirs


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$", sggobs()["__oc__"])

    prsr._rgmn("-m", "--mrg-pirs", s"mrg_pirs",
                      cion"sor_r",
                      hp"mrg pir-n rs n op inrv "
                      "or nir rgmn []. ")

    prsr._rgmn("--mx-insr-siz", s"mx_insr_siz", yp"in",
                      hp"ony mrg pir-n rs i hy r ss hn "
                      "# bss pr. "
                      " 0 rns o his ir []. ")

    prsr._rgmn("--min-insr-siz", s"min_insr_siz", yp"in",
                      hp"ony mrg pir-n rs i hy r  "
                      "s # bss pr. "
                      " 0 rns o his ir []. ")

    prsr._rgmn("--b-orm", s"b_orm", yp"choic",
                      choics('3', '4', '5', '6'),
                      hp"b orm o op. "
                      " []")

    prsr.s_s(
        rgionNon,
        c_pksNon,
        mrg_pirsNon,
        min_insr_siz0,
        mx_insr_siz0,
        b_orm'6',
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    i n(rgs)  0:
        rgs.ppn("-")

    smi  pysm.AignmnFi(rgs[0], "rb")

    opions.b_orm  in(opions.b_orm)

    i opions.mrg_pirs is no Non:
        conr  mrg_pirs(smi,
                              opions.so,
                              min_insr_sizopions.min_insr_siz,
                              mx_insr_sizopions.mx_insr_siz,
                              b_ormopions.b_orm)

        E.ino("cgory\cons\ns\n"  conr.sTb())

    s:
        # s ni_o. Fis rom sin hv no inx
        i  smi.ch(ni_oTr)

        # mor comorb cigr prsing wi
        # com wih h nx pysm rs
        BAM_CMATCH  0
        BAM_CDEL  2
        BAM_CREF_SKIP  3
        k  (BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP)
        oi  opions.so

        or r in i:
            i r.is_nmpp:
                conin

              0
            or op,  in r.cigr:
                i op in k:
                     + 

            i r.is_rvrs:
                srn  "-"
            s:
                srn  "+"
            oi.wri("s\\\s\\c\n" 
                          (r.rrnc_nm,
                           r.pos,
                           r.pos + ,
                           r.qnm,
                           r.mpq,
                           srn))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
