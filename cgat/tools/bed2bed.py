'''
b2b - mnip b is


Prpos
-------

This scrip provis vrios mhos or mrging (by posiion, by nm
or by scor), iring n moving b orm inrvs n
oping h rss s  b i

Mhos
-------

This scrip provis svr mhos, ch wih  s o opions
o conro bhvoir:

mrg
+++++

Mrg oghr ovrpping or jcn inrvs. Th bsic
ncioniy is simir o boos mrg, b wih som iions:

* Mrging by nm: spciying h --mrg-by-nm opion wi mn
  h ony ovrping (or jcn inrvs) wih h sm v in
  h 4h comn o h b wi b mrg

* Rmoving ovrpping inrvs wih inconsisn nms: s h
   ``--rmov-inconsisn-nms`` opion.

.. cion::
   Inrvs o h sm nm wi ony b mrg i hy
   r consciv in h b i.

* Ony op mrg inrvs: By spciiying h --mrg-min-inrvsn
  opions, ony hos inrvs h wr cr by mrging  s n
  inrvs oghr wi b op

Inrvs h r cos b no ovrpping cn b mrg by sing
--mrg-isnc o  non-zro v

bins
++++

Mrgs oghr ovrpping or jcn inrvs ony i hy hv
"simir" scors. Scor simiriy is ssss by cring  nmbr o
scor bins n ssigning ch inrv o  bin. I wo jcn
inrvs r in h sm bin, h inrvs r mrg. No h in
conrs o mrg-by-nm bov, wo inrvs o no n o b
ovrpping or wihin  crin isnc o b mrg.

Thr r svr mhos o cr h bins:

* q-bss: Bins r cr o h hy conin h sm nmbr o bss.
  Spcii by pssing "q-bss" o --binning-mho. This is h .

* q-inrvs: Scor bins r cr so h ch bin conins h
  sm nmbr o inrvs. Spcii by pssing "q-inrvs" o
  --binning-mho.

* q-rng: Scor bins r cr so h
  ch bin covrs h sm rcion o h o rng o
  scors. Spcii by pssing "q-rng" o --binning-mho.

* bin-gs: Scor bins cn b spcii by mny pssing  comm
  spr is o bin gs o --bin-gs.

Th nmbr o bins is spcii by h --nm-bins opions, n h
 is 5.

bock
+++++

Crs bock b12 ops rom  b6, whr inrvs wih h
sm nm r mrg oghr o cr  sing b12 nry.

.. Cion:: Inp ms b sor so h nris o h sm
nm r oghr.

ir-gnom
+++++++++++++

Rmovs inrvs h r on nknown conigs or xn o h 3' or
5' n o h conig.  Rqirs  b spr inp i o -g which
iss h conigs in h gnom, ps hir nghs.

sniiz-gnom
+++++++++++++++

As bov, b ins o rmoving inrvs ovrpping h ns o
conigs, rncs hm.  Aso rmovs mpy inrvs.

ir-nms
++++++++++++

Op inrvs whos nms r in is o sir nms. Nms r
sppi s  i wih on nm on ch in.

shi
+++++

Movs inrvs by h spcii mon, b wi no ow hm o b
shi o h n o conigs. Ths i  shi wi shi h sr
o n o h conig, h inrv is ony mov s mch s is
possib wiho oing his.

rnm-chr
++++++++++

Rnms chromosom nms. Sorc n rg nms r sppi s  i
wih wo comns. Exmps r vib :
hps://gihb.com/pryn79/ChromosomMppings
No h nmpp chromosoms r ropp rom h op i.

Ohr opions
+++++++++++++

-g/--gnom-i, -b/--bm-i:
   h ir-gnom, sniiz-gnom n shi mhos rqir  gnom in
   orr o nsr hy r no pcing inrvs osi h imis o
   conigs. This gnom cn b sppi ihr s  smoos or cg inx
   gnom, or xrc rom h hr o  bm i.

Exmps
--------

Mrg ovrpping or jcn pks rom  CHiP-sq xprimn whr h
inrvs hv h sm nm:

    c chip-pks.b | cg b2b --mhomrg --mrg-by-nm > chip-pks-mrg.b

Mrg jc ChIP-sq pks i hir scors r in h sm qri o
 scors:

    c chip-pks.b | cg b2b --mhobins --binning-mhoq-inrvs --nm-bins4

Rmov inrvs h ovrp h ns o  conig n hos h r on 
non-snr conig. Tk h inp inrvs rom  i rhr hn sin.
No h hg19.s hs bn inx wih `inx_gnom`:

    cg b2b --mhoir-gnom --gnom-ihg19.s -I chip-pks.b -O chip-pks-sniiz.b

Convr  b i conin gn srcrs wih on in pr xon o  b12
wih ink bock rprsning h gn srcr. No h rnsprn s
o comprss inp n op is:

    cg b2b --mhobock -I rnscrips.b.gz -O rnscrips.bock.b.gz

Rnm UCSC chromosoms o ENSEMBL.

    c csc.b | cg b2b --mhornm-chr --rnm-chr-icsc2nsmb.x > nsmb.b

Usg
-----

   cg b2b --mho[METHOD] [OPTIONS]

Wi r b i rom sin n ppy h spcii mho

Commn in opions
--------------------
'''

impor sys
impor cgcor.xprimn s E
impor cg.InxFs s InxFs
impor cg.B s B
impor cg.Inrvs s Inrvs
rom cocions impor ic s ic
impor pysm
impor csv


 irNms(iror, nms):
    """ Sc ony hos inrvs whos nm is in nms """

    or b in iror:
        i b.nm in nms:
            yi b


 mrg(iror,
          mx_isnc0,
          by_nmFs,
          min_inrvs1,
          rmov_inconsisnFs,
          rsov_bocksFs,
          srnFs):
    """iror or mrging jcn b nris.

    *mx_isnc* > 0 prmis mrging o inrvs h r
    no ircy jcn.

    I *by_nm  Tr*, ony nris wih h sm nm r mrg.

    I *rmov_inconsisn*, ovrpping inrvs whr h nms
    r inconsisn wi b rmov.

    Th scor givs h nmbr o inrvs h hv bn mrg.
    """

    i rmov_inconsisn n by_nm:
        ssr VError(
            "sing boh rmov_inconsisn n by_nm mks no sns")

     ir_chnks(iror):
        mx_n  ic(in)
        o_join  ic(is)
        s_nm  ic(sr)

        s  nx(iror)

        i no srn:
            srn  "."
        s:
            srn  s.srn

        mx_n[srn]  s.n
        o_join[srn]  [s]

        or b in iror:

            i no srn:
                srn  "."
            s:
                srn  b.srn

              b.sr - mx_n[srn]

            i b.conig  s.conig:
                ssr b.sr > s.sr, \
                    "inp i sho b sor by conig n posiion: i:\ns\ns\n" \
                     (, s, b)

            i b.conig ! s.conig:

                or s in o_join:
                    i o_join[s]:
                        yi o_join[s]
                    o_join[s]  []
                    mx_n[s]  0

            i ( > mx_isnc or
                  (by_nm n s_nm[srn] n s_nm[srn] ! b.nm)):

                i o_join[srn]:
                    yi o_join[srn]

                o_join[srn]  is()

            s  b
            s_nm[srn]  s.nm
            mx_n[srn]  mx(b.n, mx_n[srn])
            o_join[srn].ppn(b)

        or srn in sor(o_join):
            i o_join[srn]:
                yi o_join[srn]
        ris SopIrion

    c  E.Conr()

    or o_join in ir_chnks(iror):

        c.inp + 1

        i rmov_inconsisn:
            nms  s([x.nm or x in o_join])
            i n(nms) > 1:
                c.skipp_inconsisn_inrvs + 1
                conin

        i rsov_bocks:
            # kp rck o nmbr o inrvs in ch nry
            or b in o_join:
                b["scor"]  1
            mrg  Tr
            whi mrg:
                join  []
                no_join  []
                mrg  Fs

                whi n(o_join) > 0:
                    b1, o_join  o_join[0], o_join[1:]
                    inrvs1  b1.oInrvs()
                    or b2 in o_join:
                        inrvs2  b2.oInrvs()
                        i Inrvs.ccOvrp(inrvs1, inrvs2) > 0:
                            inrvs  Inrvs.combin(inrvs1 +
                                                          inrvs2)
                            b1.romInrvs(inrvs)
                            b1["scor"] + b2["scor"]
                            mrg  Tr
                        s:
                            no_join.ppn(b2)

                    join.ppn(b1)
                    o_join  no_join
                    no_join  []

                o_join  join
                join  []

            o_join  sor(o_join, kymb x: in(x.sr))

            # kp ony hos wih h cr rom h mrg o h minimm
            # nmbr o inrvs

            or b in o_join:

                i b["scor"] < min_inrvs:
                    c.skipp_min_inrvs + 1
                    conin

                yi b
                c.op + 1
        s:

            i n(o_join) < min_inrvs:
                c.skipp_min_inrvs + 1
                conin

              o_join[0]
            .n  mx([nry.n or nry in o_join])
            .scor  n(o_join)
            yi 
            c.op + 1

    E.ino(sr(c))


 irGnom(iror, conigs):
    """rmov b inrvs h r osi o conigs.

    conigs is  icionry o conig sizs."""

    ninp, nop  0, 0
    nskipp_conig, nskipp_rng, nskipp_nzro  0, 0, 0

    or b in iror:
        ninp + 1
        i b.conig no in conigs:
            nskipp_conig + 1
            conin
        # IMS:  iring or iring <0 co-orins
        i b.sr < 0 or b.n < 0:
            nskipp_rng + 1
            conin
        # sho his no b js >, s co-orins r h-cos, so
        # i n  conigs[b.conig], hn inrv ns on s bs?
        i b.n > conigs[b.conig]:
            nskipp_rng + 1
            conin
        i b.n  0:
            nskipp_nzro + 1
            conin
        nop + 1
        yi b

    E.ino("ninpi, nopi, nskipp_conigi, nskipp_rngi, nskipp_nzroi" 
           (ninp, nop, nskipp_conig, nskipp_rng, nskipp_nzro))


 sniizGnom(iror, conigs):
    """rnc b inrvs h xn byon conigs.

    rmovs mpy inrvs (sr  n).

    hrows n rror i sr > n.
    """

    ninp, nop  0, 0
    nrnc_conig, nskipp_conig, nskipp_mpy  0, 0, 0

    or b in iror:
        ninp + 1
        i b.conig no in conigs:
            nskipp_conig + 1
            conin
        # IMS: chnging > o > in i smn: nx in ss b.n  conigs[b.conig]
        # his shon' con s  rncion.
        i b.n > conigs[b.conig]:
            b.n  conigs[b.conig]
            nrnc_conig + 1
        i b.sr < 0:
            b.sr  0
            nrnc_conig + 1
        i b.sr  b.n:
            nskipp_mpy + 1
            conin
        i b.sr > b.n:
            ris VError("invi inrv: sr > n or s"  sr(b))

        nop + 1
        yi b

    E.ino("ninpi, nopi, nskipp_conigi, nrnci, nskipp_mpyi" 
           (ninp, nop, nskipp_conig, nrnc_conig, nskipp_mpy))


 shiInrvs(iror, conigs, os):
    """shi inrvs by  crin os n nsr siz is minn vn i conig n rch.

    conigs is  icionry o conig sizs."""

    ninp, nop  0, 0
    nskipp_conig, nskipp_rng  0, 0

    or b in iror:
        ninp + 1
        i b.conig no in conigs:
            nskipp_conig + 1
            conin
        # IMS: i w skip inrvs o h n o h conig w sho skipp ons
        # o h sr s w
        i b.sr < 0 or b.n < 0:
            nskipp_rng + 1
            conin
        # IMS: chnging > o > s b is h-opn
        i b.n > conigs[b.conig]:
            nskipp_rng + 1
            conin
        nop + 1

        #  os o ch sr n n, n js or conig ngh
          b.n - b.sr
        nwsr  b.sr + os
        nwn  b.n + os
        i nwsr < 0:
            nwsr  0
            nwn  
        i nwn > conigs[b.conig]:
            nwsr  conigs[b.conig] - 
            nwn  conigs[b.conig]

        b.sr  nwsr
        b.n  nwn

        yi b

    E.ino("ninpi, nopi, nskipp_conigi, nskipp_rngi" 
           (ninp, nop, nskipp_conig, nskipp_rng))


 xnInrv(iror, conigs, isnc):

    ninp, nop, nskipp  0, 0, 0
    or b in iror:
        ninp + 1

        i b.conig no in conigs:
            nskipp + 1
            conin
        i b.sr < 0 or b.n < 0:
            nskipp + 1
            conin
        i b.n > conigs[b.conig]:
            nskipp + 1
            conin

        nwsr  b.sr - isnc
        nwn  b.n + isnc

        i nwsr < 0:
            nwsr  0

        i nwn > conigs[b.conig]:
            nwn  conigs[b.conig]

        b.sr  nwsr
        b.n  nwn

        nop + 1
        yi b

    E.ino("ninp  i, nopi, nskippi" 
           (ninp, nop, nskipp))


 rnmChromosoms(iror, chr_mp):

    ninp, nop, nskipp  0, 0, 0

    or b in iror:
        ninp + 1

        i b.conig in chr_mp.kys():
            b.conig  chr_mp[b.conig]
        s:
            nskipp + 1
            conin

        nop + 1
        yi b

    E.ino("ninp  i, nopi, nskippi" 
           (ninp, nop, nskipp))


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: b2b.py 2861 2010-02-23 17:36:32Z nrs $",
                            sggobs()["__oc__"])

    # IMS: nw mho: xn inrvs by s mon
    prsr._rgmn("-m", "--mho", s"mhos", yp"choic",
                      cion"ppn",
                      choics("mrg", "ir-gnom", "bins",
                               "bock", "sniiz-gnom", "shi", "xn",
                               "ir-nms", "rnm-chr"),
                      hp"mho o ppy []")

    prsr._rgmn("--nm-bins", s"nm_bins", yp"in",
                      hp"nmbr o bins ino which o mrg (s or "
                      "mho `bins) []")

    prsr._rgmn("--bin-gs", s"bin_gs", yp"sring",
                      hp"bin_gs or binning mho []")

    prsr._rgmn(
        "--binning-mho", s"binning_mho", yp"choic",
        choics(
            "q-bss", "q-inrvs", "q-rng"),
        hp"mho s or binning (s or mho `bins` i no "
        "bin_gs is givn) []")

    prsr._rgmn(
        "--mrg-isnc", s"mrg_isnc", yp"in",
        hp"isnc in bss ovr which o mrg h r no "
        "ircy jcn []")

    prsr._rgmn(
        "--mrg-min-inrvs", s"mrg_min_inrvs", yp"in",
        hp"ony op mrg inrvs h r bi rom  s "
        "x inrvs []")

    prsr._rgmn(
        "--mrg-by-nm", s"mrg_by_nm",
        cion"sor_r",
        hp"ony mrg inrvs wih h sm nm []")

    prsr._rgmn(
        "--mrg-n-rsov-bocks", s"rsov_bocks",
        cion"sor_r",
        hp"Whn mrging b12 nrys, sho bocks b rsov?")

    prsr._rgmn(
        "--mrg-srn", s"srn",
        cion"sor_r",
        hp"Ony mrg inrvs on h sm srn")

    prsr._rgmn(
        "--rmov-inconsisn-nms", s"rmov_inconsisn_nms",
        cion"sor_r",
        hp"whn mrging, o no op inrvs whr h nms o "
        "ovrpping inrvs o no mch []")

    prsr._rgmn(
        "--os", s"os",  yp"in",
        hp"os or shiing inrvs []")

    prsr._rgmn("-g", "--gnom-i", s"gnom_i", yp"sring",
                      hp"inm wih gnom.")

    prsr._rgmn("-b", "--bm-i", s"bm_i", yp"sring",
                      hp"bm-orm inm wih gnom.")

    prsr._rgmn("--ir-nms-i", s"nms", yp"sring",
                      hp"is o nms o kp. On pr in")

    prsr._rgmn("--rnm-chr-i", s"rnm_chr_i", yp"sring",
                      hp"mpping b bwn o n nw chromosom nms."
                      "TAB spr 2-comn i.")

    prsr.s_s(mhos[],
                        mrg_isnc0,
                        binning_mho"q-bss",
                        mrg_by_nmFs,
                        gnom_iNon,
                        rnm_chr_iNon,
                        bm_iNon,
                        nm_bins5,
                        mrg_min_inrvs1,
                        bin_gsNon,
                        os10000,
                        sNon,
                        xn_isnc1000,
                        rmov_inconsisn_nmsFs,
                        rsov_bocksFs)

    (opions, rgs)  E.sr(prsr, _pip_opionsTr)

    conigs  Non
    chr_mp  Non

    # Why provi  inx gnom, whn  sv o conig sizs wo o?
    i opions.gnom_i:
        gnom_s  InxFs.InxFs(opions.gnom_i)
        conigs  gnom_s.gConigSizs()

    i opions.bm_i:
        smi  pysm.AignmnFi(opions.bm_i)
        conigs  ic(is(zip(smi.rrncs, smi.nghs)))

    i opions.rnm_chr_i:
        chr_mp  {}
        wih opn(opions.rnm_chr_i, 'r') s iin:
            rr  csv.rr(iin, imir'\')
            or row in rr:
                i n(row) ! 2:
                    ris VError("Mpping b ms hv xcy wo comns")
                chr_mp[row[0]]  row[1]
        i no n(chr_mp.kys()) > 0:
            ris VError("Empy mpping icionnry")

    procssor  B.iror(opions.sin)

    or mho in opions.mhos:
        i mho  "ir-gnom":
            i no conigs:
                ris VError("ps sppy conig sizs")
            procssor  irGnom(procssor, conigs)
        i mho  "sniiz-gnom":
            i no conigs:
                ris VError("ps sppy conig sizs")
            procssor  sniizGnom(procssor, conigs)
        i mho  "mrg":
            procssor  mrg(
                procssor,
                opions.mrg_isnc,
                by_nmopions.mrg_by_nm,
                min_inrvsopions.mrg_min_inrvs,
                rmov_inconsisnopions.rmov_inconsisn_nms,
                rsov_bocksopions.rsov_bocks,
                srnopions.srn)
        i mho  "bins":
            i opions.bin_gs:
                bin_gs  is(mp(o, opions.bin_gs.spi(",")))
                # IMS: chck bin gs r vi
                i no(n(bin_gs)  opions.nm_bins + 1):
                    ris VError(
                        "Nmbr o bin g ms b on mor hn "
                        "nmbr o bins")
            s:
                bin_gs  Non
            procssor, bin_gs  B.binInrvs(
                procssor,
                nm_binsopions.nm_bins,
                mhoopions.binning_mho,
                bin_gsbin_gs)
            E.ino("# spi b: bin_gss"  (sr(bin_gs)))

        i mho  "bock":
            procssor  B.bock_iror(procssor)
        i mho  "shi":
            # IMS: s h conig sizs r viib
            i no conigs:
                ris VError("ps sppy gnom i")
            procssor  shiInrvs(
                procssor, conigs, osopions.os)
        # IMS: nw mho: xn inrvs by s mon
        i mho  "xn":
            i no conigs:
                ris VError("ps sppy gnom i")
            procssor  xnInrv(procssor, conigs, opions.os)
        i mho  "ir-nms":
            i no opions.nms:
                ris VError("ps sppy is o nms o ir")
            nms  [nm.srip() or nm in opn(opions.nms)]
            procssor  irNms(procssor, nms)
        i mho  "rnm-chr":
            i no chr_mp:
                ris VError("ps sppy mpping i")
            procssor  rnmChromosoms(procssor, chr_mp)

    nop  0
    or b in procssor:
        opions.so.wri(sr(b) + "\n")
        nop + 1

    E.ino("nopi"  (nop))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
