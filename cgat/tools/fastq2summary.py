'''sq2smmry.py - comp smmry ss or  sq i


:Tgs: Gnomics NGS Sqncs FASTQ Annoion

Prpos
-------

This scrip irs ovr  sq i n ops
smmry sisics or h comp i

Th op is  b-imi x i wih h som o oowing comns
pning on h opion spcii:


+----------------+-----------------------------------------------------------+
|*Comn*        |*Conn*                                                  |
+----------------+-----------------------------------------------------------+
|rs           |o rs in i                                        |
+----------------+-----------------------------------------------------------+
|bss           |o bss in i                                        |
+----------------+-----------------------------------------------------------+
|mn_ngh     |mn r ngh                                           |
+----------------+-----------------------------------------------------------+
|min_ngh   |min r ngh                                         |
+----------------+-----------------------------------------------------------+
|mn_qiy    |mn r qiy                                          |
+----------------+-----------------------------------------------------------+
|min_qiy  |min r qiy                                        |
+----------------+-----------------------------------------------------------+
|ni         |nmbr o bss bow qiy hrsho                    |
+----------------+-----------------------------------------------------------+


Usg
-----

Exmp::

   pyhon sq2smmry.py --gss-ormsngr < in.sq > o.sv

In his xmp w know h or  hv qiy scors orm s
sngr. Givn h imin-1.8 qiy scors r highy ovrpping
wih sngr, his opion s o sngr qiis. In  mo
h scrip my no b b o isingish highy ovrpping ss o
qiy scors.

Typ::

   pyhon sq2smmry.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor nmpy s np
impor cgcor.xprimn s E
impor cg.Fsq s Fsq


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--gss-orm", s"gss_orm", yp"choic",
        choics('sngr', 'sox', 'phr64',
                 'imin-1.8', 'ingr'),
        hp"Th  bhvior o h scrip is o gss \
        h qiy orm o h inp sq i. Th sr \
        cn spciy h qiy orm o h inp i sing \
        h --orm opion. Th scrip wi s his orm i \
        sqncs qiis r mbigos.[].")

    prsr._rgmn(
        "-", "--rg-orm", s"chng_orm",
        yp"choic", choics('sngr', 'sox', 'phr64',
                                'imin-1.8', 'ingr'),
        hp"Th scrip gsss h qiy orm o h inp \
        i n convrs qiy scors o h sinion \
        orm nss --orm is spcii [].")

    prsr.s_s(
        chng_ormNon,
        gss_ormNon,
        min_qiy10)

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.chng_orm:
        iror  Fsq.ir_convr(opions.sin,
                                         ormopions.chng_orm,
                                         gssopions.gss_orm)
    s:
        iror  Fsq.ir_gss(opions.sin,
                                       gssopions.gss_orm)

    min_qiy  opions.min_qiy
    nmbr_o_rs  0
    nmbr_o_bss  0
    r_nghs  []
    r_qiis  []
    bss_bow_min  0

    or rcor in iror:
        nmbr_o_rs + 1
        qs  rcor.oPhr()
        ngh_r  n(qs)
        nmbr_o_bss + ngh_r
        bss_bow_min + n([x or x in qs i x < min_qiy])
        r_nghs.ppn(ngh_r)
        r_qiis.ppn(np.mn(qs))

    mn_ngh  ron(np.mn(r_nghs), 2)
    min_ngh  ron(np.min(r_nghs), 2)
    mn_qiy  ron(np.mn(r_qiis), 2)
    min_qiy  ron(np.min(r_qiis), 2)

    opions.so.wri(
        "rs\bss\mn_ngh\min_ngh\mn_qiy\min_qiy\ni\n")

    opions.so.wri(
        "i\i\s\s\s\s\i\n"  (nmbr_o_rs, nmbr_o_bss,
                                          sr(mn_ngh),
                                          sr(min_ngh),
                                          sr(mn_qiy),
                                          sr(min_qiy),
                                          bss_bow_min))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
