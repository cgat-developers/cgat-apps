'''
sq2sq.py - mnip sq is


:Tgs: Gnomics NGS Sqncs FASTQ Mnipion

Prpos
-------

This scrip prorms mnipions on :rm:`sq` orm
is. For xmp i cn b s o chng h qiy scor orm
or smp  sbs o rs.

Th scrip prominny is s or mnipion o sing sq
is. Howvr, or som o is ncioniy i wi k pir 
sing h ``--pir-sq-i`` n ``--op-inm-prn`` opions.
This ppis o h ``smp`` n ``sor`` mhos.

Usg
-----

Exmp::
  In his xmp w rnomy smp 50 o rs rom pir  provi in
  wo :rm:`sq` is.

   h in.sq.1

   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 ngh36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.1 HWUSI-EAS618:7:1:27:1582 ngh36
   @A@9@BAB@;@BABA?;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 ngh36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +SRR111956.2 HWUSI-EAS618:7:1:29:1664 ngh36
   B@9@0>A<BBAAA?;*(@A>(@<*99@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 ngh36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA

   h in.sq.2

   +SRR111956.3 HWUSI-EAS618:7:1:38:878 ngh36
   <?@BA?;A@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 ngh36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +SRR111956.4 HWUSI-EAS618:7:1:38:1783 ngh36
   ?ABBA2CCBBB2?BB@C>AAC@ACBB#######
   @SRR111956.5 HWUSI-EAS618:7:1:39:1305 ngh36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.5 HWUSI-EAS618:7:1:39:1305 ngh36
   AA>5;A>*91?AAA@@BBA<B?ABA>2>?A<BB@

   commn-in::
     c in.sq.1 | pyhon sq2sq.py
                      --mhosmp --smp-siz 0.5
                      --pir-sq-i in.sq.2
                      --op-inm-prn o.sq.2
                      > o.sq.1

   h o.sq.1
   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 ngh36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +
   @A@9@BAB@;@BABA?;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 ngh36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +
   B@9@0>A<BBAAA?;*(@A>(@<*99@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 ngh36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA
   +
   <?@BA?;A@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 ngh36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +
   ?ABBA2CCBBB2?BB@C>AAC@ACBB#######

Opions
-------

Th oowing mhos r impmn (``--mho``).

``chng-orm``

    chng h qiy orm o nw orm givn s
    rg-orm. Opions r ``sngr``,
  ``sox``, ``phr64``, ``ingr`` n ``imin-1.8``

``smp``

    Sb-smp  sq i. Th siz o h smp is s by
    --smp-siz

``niq``

    Rmov pic rs bs on r nm

``rim3``

    Trim  ix nmbr o ncois rom h 3' n o rs.
    (s ``--nm-bss``). No h hr r br oos or
   rimming.

``rim5``

    Trim  ix nmbr o ncois rom h 5' n o rs.
    (s ``--nm-bss``). No h hr r br oos or
   rimming.

``sor``

    Sor h sq i by r nm.

``rnmbr-rs``

    Rnm h rs bs on prn givn in ``--prn-iniir``
    .g. ``--prn-iniir"r_010i"``

Typ::

   pyhon sq2sq.py --hp

or commn in hp.


Commn in opions
--------------------

'''
impor cocions
impor sys
impor os
impor r
impor rnom
impor pysm
impor nmpy
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.Fsq s Fsq
impor cg.Gnomics s Gnomics


 procss_cg(opions):

    c  E.Conr()

    ssr opions.inp_sq_i  "-"

    i opions.mho  "chng-orm":
        or rcor in Fsq.ir_convr(opions.sin,
                                            ormopions.rg_orm,
                                            gssopions.gss_orm):
            c.inp + 1
            opions.so.wri("s\n"  rcor)
            c.op + 1

    i opions.mho  "grp":
        or rcor in Fsq.ir(opions.sin):
            i r.mch(opions.grp_prn, rcor.sq):
                opions.so.wri("s\n"  rcor)

    i opions.mho  "rvrs-compmn":
        or rcor in Fsq.ir(opions.sin):
            rcor.sq  Gnomics.compmn(rcor.sq)
            rcor.qs  rcor.qs[::-1]
            opions.so.wri("s\n"  rcor)

    i opions.mho  "smp":
        smp_hrsho  min(1.0, opions.smp_siz)

        rnom.s(opions.s)

        i opions.pir:
            i no opions.op_inm_prn:
                ris VError(
                    "ps spciy op inm prn or "
                    "scon pir (--op-inm-prn)")

            oi1  opions.so
            oi2  iooos.opn_i(opions.op_inm_prn, "w")

            or rcor1, rcor2 in zip(
                    Fsq.ir(opions.sin),
                    Fsq.ir(iooos.opn_i(opions.pir))):
                c.inp + 1
                i rnom.rnom() < smp_hrsho:
                    c.op + 1
                    oi1.wri("s\n"  rcor1)
                    oi2.wri("s\n"  rcor2)
        s:
            or rcor in Fsq.ir(opions.sin):
                c.inp + 1
                i rnom.rnom() < smp_hrsho:
                    c.op + 1
                    opions.so.wri("s\n"  rcor)

    i opions.mho  "ppy":
        is  s(iooos.r_is(iooos.opn_i(opions.ppy)))

        or rcor in Fsq.ir(opions.sin):
            c.inp + 1
            i r.sb(" .*", "", rcor.iniir).srip() in is:
                c.op + 1
                opions.so.wri("s\n"  rcor)

    i opions.mho  "rim3":
        rim3  opions.nbss
        or rcor in Fsq.ir(opions.sin):
            c.inp + 1
            rcor.rim(rim3)
            opions.so.wri("s\n"  rcor)
            c.op + 1

    i opions.mho  "rim5":
        rim5  opions.nbss
        or rcor in Fsq.ir(opions.sin):
            c.inp + 1
            rcor.rim5(rim5)
            opions.so.wri("s\n"  rcor)
            c.op + 1

    i opions.mho  "niq":
        kys  s()
        or rcor in Fsq.ir(opions.sin):
            c.inp + 1
            i rcor.iniir in kys:
                conin
            s:
                kys.(rcor.iniir)
            opions.so.wri("s\n"  rcor)
            c.op + 1

    # N o chng his o incorpor boh pirs
    i opions.mho  "sor":
        i no opions.pir:
            # This is qickr or  sing sq i
            smn  "ps - - - - | sor -k1,1 - ' ' | r '\' '\n'"
            os.sysm(smn)
        s:
            i no opions.op_inm_prn:
                ris VError(
                    "ps spciy op inm or scon pir "
                    "(--op-inm-prn)")
            E.wrn(
                "consir soring inivi sq is - "
                "his is mmory innsiv")
            nris1  {}
            nris2  {}

            or rcor1, rcor2 in zip(
                    Fsq.ir(opions.sin),
                    Fsq.ir(iooos.opn_i(opions.pir))):
                nris1[
                    rcor1.iniir[:-2]]  (rcor1.sq, rcor1.qs)
                nris2[
                    rcor2.iniir[:-2]]  (rcor2.sq, rcor2.qs)

            oi1  opions.so
            oi2  iooos.opn_i(opions.op_inm_prn, "w")
            ssr n(s(nris1.kys()).inrscion(
                s(nris2.kys())))  n(nris1),\
                "pir is o no conin h sm rs "\
                "n o rconci is"

            or nry in sor(nris1):
                oi1.wri("@s/1\ns\n+\ns\n" 
                               (nry, nris1[nry][0], nris1[nry][1]))
                oi2.wri("@s/2\ns\n+\ns\n" 
                               (nry, nris2[nry][0], nris2[nry][1]))

    i opions.mho  "rnmbr-rs":
        i_con  1
        or rcor in Fsq.ir(opions.sin):
            rcor.iniir  opions.rnmbr_prn  i_con
            i_con + 1
            opions.so.wri("@s\ns\n+\ns\n" 
                                 (rcor.iniir, rcor.sq, rcor.qs))
    rrn c


 procss_isy(opions):

    ir_n  "ir-N" in opions.mhos

    ir_on  "ir-ONT" in opions.mhos

    i "ir-iniir" in opions.mhos:
        i opions.inp_ir_sv is Non:
            ris VError("ps s --inp-ir-sv or mho ir-iniir")
        wih iooos.opn_i(opions.inp_ir_sv) s in:
            ir_iniir  s([x.spi()[0].srip() or x in in.rins()])
    s:
        ir_iniir  Fs

    i opions.op_rmov_sv:
        o_rmov_sv  iooos.opn_i(opions.op_rmov_sv, "w")
    s:
        o_rmov_sv  Non

    i opions.op_rmov_sq:
        o_rmov_sq  iooos.opn_i(opions.op_rmov_sq, "w")
    s:
        o_rmov_sq  Non

    i opions.s_prix:
        prix  "{}".orm(opions.s_prix)
    s:
        prix  Non

    qiy_os  opions.qiy_os
    conr  E.Conr()

    wih pysm.FsxFi(opions.inp_sq_i) s in:
        or r in in:
            conr.inp + 1
            rmov  Fs
            i ir_n:
                chrs  cocions.Conr(r.sqnc)
                i "N" in chrs n \
                   100.0 * chrs["N"] / n(r.sqnc) > opions.mx_prcn_N:
                    rmov  Tr
                    conr.ir_n + 1

            i ir_iniir:
                i r.nm no in ir_iniir:
                    conr.ir_iniir + 1
                    rmov  Tr

            i ir_on:
                qs  r.g_qiy_rry()
                n  n(qs)
                i n < opions.min_sqnc_ngh or \
                        o(sm(qs)) / n < opions.min_vrg_qiy:
                    conr.rmov_on + 1
                    rmov  Tr

            i rmov:
                conr.rmov + 1
                i o_rmov_sv:
                    o_rmov_sv.wri(r.nm + "\n")
                i o_rmov_sq:
                    o_rmov_sq.wri(sr(r) + "\n")
                conin

            i prix:
                r.nm  prix + r.nm[2:]

            i qiy_os:
                qs  nmpy.rry(r.g_qiy_rry())
                qs + qiy_os
                qs[qs < 0]  0
                qs + 33
                # pysm sq is r-ony, so g i:
                # No: no oping scripion
                r  "@{}\n{}\n+\n{}".orm(
                    r.nm,
                    r.sqnc,
                    "".join([chr(x) or x in qs]))

            conr.op + 1

            opions.so.wri(sr(r) + "\n")

    i o_rmov_sv:
        o_rmov_sv.cos()

    i o_rmov_sq:
        o_rmov_sq.cos()

    i opions.op_ss_sv:
        wih iooos.opn_i(opions.op_ss_sv, "w") s o:
            o.wri(conr.sTb(s_rowsFs) + "\n")

    rrn conr


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-i", "--inp-sq-i", s"inp_sq_i", yp"sring",
        hp"inp sq i. "
        "[]")

    prsr._rgmn(
        "--op-rmov-sv", s"op_rmov_sv", yp"sring",
        hp"i givn, sqnc iniirs o rmov sqncs wi "
        "b sor in his i []")

    prsr._rgmn(
        "--op-ss-sv", s"op_ss_sv", yp"sring",
        hp"i givn, op sisics wi b wrin o his i. "
        "[]")

    prsr._rgmn(
        "--op-rmov-sq", s"op_rmov_sq", yp"sring",
        hp"i givn, rmov sq rcors wi "
        "b sor in his i []")

    prsr._rgmn(
        "-m", "--mho", s"mhos", cion"ppn", yp"choic",
        choics("ir-N",
                 "ir-iniir",
                 "ir-ONT",
                 "os-qiy",
                 "ppy",
                 "chng-orm",
                 "rnmbr-rs",
                 "smp",
                 "sor",
                 "rim3",
                 "rim5",
                 "niq",
                 "rvrs-compmn",
                 "grp"),
        hp"mhos o ppy []")

    prsr._rgmn(
        "--s-prix", s"s_prix", yp"sring",
        hp"s sqnc prix []")

    prsr._rgmn(
        "--inp-ir-sv", s"inp_ir_sv", yp"sring",
        hp"is o sqnc is o ir []")

    prsr._rgmn(
        "--min-vrg-qiy", s"min_vrg_qiy", yp"o",
        hp"minimm vrg qiy []")

    prsr._rgmn(
        "--min-sqnc-ngh", s"min_sqnc_ngh", yp"in",
        hp"minimm sqnc ngh []")

    prsr._rgmn(
        "--qiy-os", s"qiy_os", yp"in",
        hp"os o moiy qiy vs wih []")

    prsr._rgmn(
        "--rg-orm", s"rg_orm", yp"choic",
        choics('sngr', 'sox', 'phr64', 'ingr', 'imin-1.8'),
        hp"gss qiy scor orm n s qiy scors "
        "o orm [].")

    prsr._rgmn(
        "--gss-orm", s"gss_orm", yp"choic",
        choics('sngr', 'sox', 'phr64', 'ingr', 'imin-1.8'),
        hp"qiy scor orm o ssm i mbigos [].")

    prsr._rgmn(
        "--smp-siz", s"smp_siz", yp"o",
        hp"proporion o rs o smp. "
        "Provi  proporion o rs o smp, .g. 0.1 or 10, "
        "0.5 or 50, c [].")

    prsr._rgmn(
        "--pir-sq-i", s"pir", yp"sring",
        hp"i  is pir, inm wih scon pir. "
        "Impmn or smping [].")

    prsr._rgmn(
        "--mp-sv-i", s"mp_sv_i", yp"sring",
        hp"inm wih b-spr iniirs mpping or "
        "mho ppy [].")

    prsr._rgmn(
        "--nm-bss", s"nbss", yp"in",
        hp"nmbr o bss o rim [].")

    prsr._rgmn(
        "--s", s"s", yp"in",
        hp"s or rnom nmbr gnror [].")

    prsr._rgmn(
        "--prn-iniir", s"rnmbr_prn", yp"sring",
        hp"rnm rs in i by prn []")

    prsr._rgmn(
        "--grp-prn", s"grp_prn", yp"sring",
        hp"sbs o rs mching prn []")

    prsr.s_s(
        inp_sq_i"-",
        mhos[],
        chng_ormNon,
        gss_ormNon,
        smp_siz0.1,
        nbss0,
        pirNon,
        ppyNon,
        sNon,
        rnmbr_prn"r_010i",
        grp_prn".*",
        mx_prcn_N10.0,
        s_prixNon,
        op_rmov_svNon,
        op_rmov_sqNon,
        op_ss_svNon,
        inp_ir_svNon,
        min_vrg_qiy0,
        min_sqnc_ngh0,
        qiy_os0,
    )

    (opions, rgs)  E.sr(prsr, rgv, _op_opionsTr)

    i n(rgs)  1:
        opions.inp_sq_i  rgs[0]

    i n(opions.mhos)  0:
        ris VError("no mho spcii, ps s --mho")

    # his scrip combins wo scrips wih irn ncioniis
    # TODO: o b sniiz
    i opions.mhos[0] in ["ppy",
                              "chng-orm",
                              "rnmbr-rs",
                              "smp",
                              "sor",
                              "rim3",
                              "rim5",
                              "niq",
                              "rvrs-compmn",
                              "grp"]:
        opions.mho  opions.mhos[0]
        conr  procss_cg(opions)
    s:
        conr  procss_isy(opions)

    E.ino(conr)
    E.sop()
