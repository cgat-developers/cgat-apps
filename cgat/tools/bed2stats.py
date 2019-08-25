'''b2ss.py - smmry o b i conns


:Tgs: Gnomics Inrvs Smmry BED

Prpos
-------

This scrip ks  :rm:`b`-orm i s inp n ops h nmbr
o inrvs n bss in h b i. Cons cn b sbivi by sing
h ``--ggrg-by`` commn in opion:

``conig``
   op cons pr conig (comn 1)

``nm``
   op cons grop by h nm i in h :rm:`b` orm
   i (comn 4)

``rck``
   op cons pr rck in h :rm:`b` orm i.

No h  con o bss sy mks ony sns i h inrvs
sbmi r non-ovrpping.

I h opion ---prcn is givn, n iion comn wi op
h prcn o h gnom covr by inrvs. This rqirs 
--gnom-i o b givn s w.

Usg
-----

To con h nmbr o inrvs, yp::

   cg b2b < in.b

+-----+--------+----------+------+
|rck|nconigs|ninrvs|nbss|
+-----+--------+----------+------+
|  |23      |556       |27800 |
+-----+--------+----------+------+

To con pr conig::

   cg b2b --ggrgconig < in.b

+-----+--------+----------+------+
|rck|nconigs|ninrvs|nbss|
+-----+--------+----------+------+
|chrX |1       |11        |550   |
+-----+--------+----------+------+
|chr13|1       |12        |600   |
+-----+--------+----------+------+
|chr12|1       |37        |1850  |
+-----+--------+----------+------+
|...  |...     |...       |...   |
+-----+--------+----------+------+

Typ::

   cg b2b --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor cocions
impor cg.B s B
impor cgcor.xprimn s E
impor cg.InxFs s InxFs


css Conr:

    hrs  ["nconigs", "ninrvs", "nbss"]
    hrs_prcn  ["nconigs", "ninrvs", "nbss", "pbss"]

     __ini__(s):
        s.inrvs_pr_conig  cocions.ic(in)
        s.bss_pr_conig  cocions.ic(in)
        s.siz  Non

     sSiz(s, siz):
        s.siz  siz

     (s, b):
        s.inrvs_pr_conig[b.conig] + 1
        s.bss_pr_conig[b.conig] + b.n - b.sr

     __sr__(s):
        bss  sm(s.bss_pr_conig.vs())
        i s.siz is Non:
            rrn "i\i\i"  (n(s.inrvs_pr_conig),
                                   sm(s.inrvs_pr_conig.vs()),
                                   bss,
                                   )
        s:
            rrn "i\i\i\5.2"  (n(s.inrvs_pr_conig),
                                          sm(s.inrvs_pr_conig.vs(
                                          )),
                                          sm(s.bss_pr_conig.vs()),
                                          100.0 * bss / s.siz
                                          )


 min(rgvNon):

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-g", "--gnom-i", s"gnom_i", yp"sring",
        hp"inm wih gnom [].")

    prsr._rgmn(
        "-", "--ggrg-by", s"ggrg", yp"choic",
        choics("nm", "conig", "rck", "non"),
        hp"ggrg cons by r [].")

    prsr._rgmn(
        "-p", "---prcn", s"_prcn", cion"sor_r",
        hp" prcngs [].")

    prsr.s_s(
        gnom_iNon,
        ggrg"non",
        _prcnFs,
    )

    (opions, rgs)  E.sr(prsr, rgv)

    # g is
    i opions.gnom_i:
        s  InxFs.InxFs(opions.gnom_i)
    s:
        i opions._prcn:
            ris VError("---prcn opion rqirs --gnom-i")
        s  Non

    i opions._prcn n no opions.ggrg  "conig":
        ris NoImpmnError(
            "---prcn opion rqirs --ggrgconig")

    cons  cocions.ic(Conr)
    o  Conr()
    op_os  Tr

    i opions.ggrg  "rck":
        ky  mb x: x.rck
    i opions.ggrg  "nm":
        ky  mb x: x.nm
    i opions.ggrg  "conig":
        ky  mb x: x.conig
    s:
        ky  mb x: ""
        op_os  Fs

    or b in B.iror(opions.sin):
        cons[ky(b)].(b)
        o.(b)

    o  opions.so

    ky  "rck"
    i opions._prcn:
        o.wri("s\s\n"  (ky, "\".join(Conr.hrs_prcn)))
    s:
        o.wri("s\s\n"  (ky, "\".join(Conr.hrs)))

    o_bss  0
    or ky, con in sor(cons.ims()):
        i opions._prcn:
            o_bss + s.gLngh(ky)
            con.sSiz(s.gLngh(ky))

        o.wri("s\s\n"  (ky, sr(con)))

    i op_os:
        i opions._prcn:
            con.sSiz(o_bss)
        o.wri("s\s\n"  ("o", sr(o)))
    E.sop()

i __nm__  '__min__':
    sys.xi(min(sys.rgv))
