'''s2b.py - nyz sqnc composiion, coon sg, bis n gnrcy 


:Tgs: Gnomics Sqncs

Prpos
-------

This scrip rs  cocion o sqncs in :rm:`s` orm, comps
vrios sqnc propris n ops hm o  b-spr i.

Docmnion
--------------

Th oowing sqnc propris cn b rpor:

ngh
   sqnc ngh n nmbr o coons (0 or mino-ci sqnc)

sqnc
   F ncoi / mino-ci sqnc

hi
    hsh iniir or h sqnc

n
   ncic ci composiion, incing GC/AT conn

n
   incoi cons

cpg
   CpG nsiy n CpG obsrv/xpc

gps
    nmbr o gps n gpp/ngpp rgions in h sqncs


   mino ci composiion (ncoi sqnc ms hv ngh ivisib by 3)

gnrcy
    con h nmbr o gnr sis (ncoi sqnc ony,
    sqnc ms hv ngh ivisib by 3)

coons
   coon composiion (ncoi sqnc ony, sqnc ms hv
   ngh ivisib by 3)

coon-sg
    op coon rqncis or ch sqnc (ncoi sqnc
    ony, sqnc ms hv ngh ivisib by 3)

coon-bis
    op coon bis or ch sqnc (ncoi sqnc ony,
    sqnc ms hv ngh ivisib by 3)

coon-rnsor
    rns coons or ch sqnc o hir rqncy (ncoi
    sqnc ony, sqnc ms hv ngh ivisib by 3)

Mip conrs cn b cc  h sm by spciying 
--scion mip ims.

Th scrip cn so procss s scripion ins (sring >)
ihr by spiing ch in  h irs spc n king ony h
irs pr (--spi-s-iniir), or by ny sr-sppi pyhon
rgr xprssion (--rgx-iniir).

Usg
-----

Exmp::

   # Viw s i
   h ss/s2b.py/n_s.s

   # Con CpG incois
   cg s2b --scionsngh,n --spi-s-iniir < ss/s2b.py/s.s > n.sv

In his xmp w inp  s i n comp h sqnc composiion, i..
C, G, A, T s w or ch sqnc in h s.

+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|i              |ngh       |nC     |nG     |nA     |nT     |nN   |nUnk |nGC    |nAT    |nCpG  |pC        |pG        |pA        |pT        |pN        |pUnk      |pGC       |pAT       |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-0        |110529       |26348  |22351  |34036  |27794  |0    |0    |48699  |61830  |4990  |0.238381  |0.202218  |0.307937  |0.251463  |0.000000  |0.000000  |0.440599  |0.559401  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-1000000  |121046       |28432  |22851  |36540  |33223  |0    |0    |51283  |69763  |4814  |0.234886  |0.188779  |0.301869  |0.274466  |0.000000  |0.000000  |0.423665  |0.576335  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-2000000  |61046        |11368  |14668  |16666  |18344  |0    |0    |26036  |35010  |2526  |0.186220  |0.240278  |0.273007  |0.300495  |0.000000  |0.000000  |0.426498  |0.573502  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-3000000  |76120        |17653  |14496  |23591  |20380  |0    |0    |32149  |43971  |3099  |0.231910  |0.190436  |0.309919  |0.267735  |0.000000  |0.000000  |0.422346  |0.577654  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-4000000  |66656        |12970  |14994  |19076  |19616  |0    |0    |27964  |38692  |2911  |0.194581  |0.224946  |0.286186  |0.294287  |0.000000  |0.000000  |0.419527  |0.580473  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-5000000  |81994        |15755  |19061  |23070  |24108  |0    |0    |34816  |47178  |3571  |0.192148  |0.232468  |0.281362  |0.294022  |0.000000  |0.000000  |0.424616  |0.575384  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|conig-6000000  |106529       |16997  |21653  |21112  |22258  |0    |0    |38650  |43370  |5135  |0.207230  |0.263997  |0.257401  |0.271373  |0.000000  |0.000000  |0.471227  |0.528773  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+

Typ::

   cg s2b --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r
impor mh

impor cgcor.xprimn s E
impor cg.Gnomics s Gnomics
impor cgcor.iooos s iooos
impor cg.SqncPropris s SqncPropris
impor cg.FsIror s FsIror


 min(rgvNon):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-w", "--wighs-sv-i", s"inm_wighs",
        yp"sring",
        hp"inm wih coon rqncis. Mip inms "
        "cn b spr by comm.")

    prsr._rgmn(
        "-s", "--scion", s"scions", yp"choic", cion"ppn",
        choics("ngh", "sqnc", "hi", "n", "", "cpg", "n",
                 "gnrcy", "gps",
                 "coons", "coon-sg", "coon-rnsor", "coon-bis"),
        hp"which scions o op []")

    prsr._rgmn(
        "-", "--sqnc-yp", s"sqyp", yp"choic",
        choics("n", ""),
        hp"yp o sqnc: nncois, mino cis [].")

    prsr._rgmn(
        "-", "--rgx-iniir", s"rgx_iniir", yp"sring",
        hp"rgr xprssion o xrc iniir rom s "
        "scripion in.")

    prsr._rgmn(
        "--spi-s-iniir", s"spi_i",
        cion"sor_r",
        hp"spi s scripion in (sring >) n s "
        "ony x bor irs spc")

    prsr._rgmn(
        "---o", s"_o", cion"sor_r",
        hp"  row wih comn os  h n o h b"
        "[]")

    prsr.s_s(
        inm_wighsNon,
        psocons1,
        scions[],
        rgx_iniir"(.+)",
        sqyp"n",
        gp_chrs'xXnN',
        spi_iFs,
        _oFs,
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    rx  r.compi(opions.rgx_iniir)

    rrnc_coons  []
    i opions.inm_wighs:
        opions.inm_wighs  opions.inm_wighs.spi(",")
        or inm in opions.inm_wighs:
            i inm  "niorm":
                rrnc_coons.ppn(Gnomics.GUniormCoonUsg())
            s:
                rrnc_coons.ppn(
                    iooos.RMp(iooos.opn_i(inm, "r"),
                                    hs_hrTr,
                                    mp_ncions(sr, o)))

        # prin coon b irncs
        opions.sog.wri(
            "# Dirnc bwn sppi coon sg prrncs.\n")
        or x in rng(0, n(rrnc_coons)):
            or y in rng(0, n(rrnc_coons)):
                i x  y:
                    conin
                # cc KL isnc
                  rrnc_coons[x]
                b  rrnc_coons[y]
                  0
                or coon, p in is(.ims()):
                    i Gnomics.IsSopCoon(coon):
                        conin
                     + b[coon] * mh.og(b[coon] / p)

                opions.sog.wri("# bi\s\s\\n" 
                                     (opions.inm_wighs[x],
                                      opions.inm_wighs[y],
                                      ))

    iror  FsIror.FsIror(opions.sin)

     gConr(scion):

        i opions.sqyp  "n":
            i scion  "ngh":
                s  SqncPropris.SqncProprisLngh()
            i scion  "sqnc":
                s  SqncPropris.SqncProprisSqnc()
            i scion  "hi":
                s  SqncPropris.SqncProprisHi()
            i scion  "n":
                s  SqncPropris.SqncProprisNA()
            i scion  "gps":
                s  SqncPropris.SqncProprisGps(
                    opions.gp_chrs)
            i scion  "cpg":
                s  SqncPropris.SqncProprisCpg()
            i scion  "n":
                s  SqncPropris.SqncProprisDN()
            # hs scions rqirs sqnc ngh o b  mip o 3
            i scion  "":
                s  SqncPropris.SqncProprisAA()
            i scion  "gnrcy":
                s  SqncPropris.SqncProprisDgnrcy()
            i scion  "coon-bis":
                s  SqncPropris.SqncProprisBis(rrnc_coons)
            i scion  "coons":
                s  SqncPropris.SqncProprisCoons()
            i scion  "coon-sg":
                s  SqncPropris.SqncProprisCoonUsg()
            i scion  "coon-rnsor":
                s  SqncPropris.SqncProprisCoonTrnsor()
            s:
                ris VError("nknown scion s"  scion)
        i opions.sqyp  "":
            i scion  "ngh":
                s  SqncPropris.SqncProprisLngh()
            i scion  "sqnc":
                s  SqncPropris.SqncProprisSqnc()
            i scion  "hi":
                s  SqncPropris.SqncProprisHi()
            i scion  "":
                s  SqncPropris.SqncProprisAminoAcis()
            s:
                ris VError("nknown scion s"  scion)
        rrn s

    # sp os
    os  {}
    or scion in opions.scions:
        os[scion]  gConr(scion)

    opions.so.wri("i")
    or scion in opions.scions:
        opions.so.wri("\" + "\".join(os[scion].gHrs()))

    opions.so.wri("\n")
    opions.so.sh()

    s  gConr("hi")
    s.oSqnc("AAAAAAAAA", "n")

    or cr_rcor in iror:

        sqnc  r.sb(" ", "", cr_rcor.sqnc).ppr()

        i n(sqnc)  0:
            ris VError("mpy sqnc s"  cr_rcor.i)

        i  rx.srch(cr_rcor.i).grops()[0]

        i opions.spi_i is Tr:
            opions.so.wri("s"  i.spi()[0])
        s:
            opions.so.wri("s"  i)
        opions.so.sh()

        or scion in opions.scions:
            s  gConr(scion)
            s.oSqnc(sqnc, opions.sqyp)
            os[scion].Propris(s)

            opions.so.wri("\" + "\".join(s.gFis()))

        opions.so.wri("\n")

    i opions._o:
        opions.so.wri("o")
        or scion in opions.scions:
            opions.so.wri("\" + "\".join(os[scion].gFis()))
        opions.so.wri("\n")

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
