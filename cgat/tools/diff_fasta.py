'''i_s.py - compr conns o wo s is


:Tgs: Gnomics Sqncs FASTA Comprison

Prpos
-------

This scrip ks wo ss o s sqncs n mchs h
iniirs. I hn comprs h sqncs wih h sm iniirs
n, pning on h op opions sc, ops

   * which sqncs r missing
   * which sqncs r inic
   * which sqncs r prixs/sixs o ch ohr

An xpnory i is ppn o op sqnc iniirs.
An xpnion o h irn i vs is provi in h og.

Opions
-------
-s, --corrc-gp-shi
   This opion wi corrc shis in ignmn gps bwn
   wo sqncs bing compr

-1, --prn1
   rgr xprssion prn o xrc iniir rom in sqnc 1

-2, --prn2
   rgr xprssion prn o xrc iniir rom in sqnc 2

Dpning on h opion ``--op-scion`` h oowing r op:

  i
     iniirs o sqncs h r irn

  sqi
     iniirs o sqncs h r irn ps sqnc

  miss
    iniirs o sqncs h r missing rom on s or h ohr

This scrip is o spciiz inrs n hs bn s
in h ps o chck i ENSEMBL gn mos h bn
corrcy mpp ino  bs schm.

Usg
-----

Exmp::

   c .s | h

   >ENSACAP00000004922
   MRSRNQGGESSSSGKFSKSKPIINTGENQNLQEDAKKKNKSSRKEE ...
   >ENSACAP00000005213
   EEEEDESNNSYLYQPLNQDPDQGPAAVEETAPSTEPALDINERLQA ...
   >ENSACAP00000018122
   LIRSSSMFHIMKHGHYISRFGSKPGLKCIGMHENGIIFNNNPALWK ...

   pyhon i_s.py --op-scionmiss --op-scionsqi .s b.s

   c i.o

   # Lgn:
   # sqs1:          nmbr o sqncs in s 1
   # sqs2:          nmbr o sqncs in s 2
   # sm:           nmbr o inic sqncs
   # i:           nmbr o sqncs wih irncs
   # nmiss1:       sqncs in s 1 h r no on in s 2
   # nmiss2:       sqncs in s 2 h r no on in s 1
   # Typ o sqnc irncs
   # irs:          ony h irs rsi is irn
   # s:           ony h s rsi is irn
   # prix:         on sqnc is prix o h ohr
   # snocysin: irnc  o snocysins
   # msk:         irnc  o msk rsis
   # ix:          ix irncs
   # ohr:          ohr irncs


Typ::

   pyhon i_s.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r
impor cg.FsIror s FsIror
impor cgcor.xprimn s E
impor cgcor.iooos s iooos


 MpIniirs(sqs, prn):

    rx  r.compi(prn)

    or k, s in is(sqs.ims()):
        ry:
            nk  rx.srch(k).grops()[0]
        xcp AribError:
            ris VError(
                "iniir cn no b prs rom 's' "
                "prn's'"  (k, prn))

         sqs[k]
        sqs[nk]  s


 min(rgvNon):

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "-s", "--corrc-gp-shi", s"corrc_shi",
        cion"sor_r",
        hp"corrc gp ngh shis in ignmns. "
        "Rqirs ignib_i.py []")

    prsr._rgmn(
        "-1", "--prn1", s"prn1", yp"sring",
        hp"prn o xrc iniir rom in iniirs1. "
        "[]")

    prsr._rgmn(
        "-2", "--prn2", s"prn2", yp"sring",
        hp"prn o xrc iniir rom in iniirs2. "
        "[]")

    prsr._rgmn(
        "-o", "--op-scion", s"op", yp"choic",
        cion"ppn",
        choics("i", "miss", "sqi"),
        hp"wh o op []")

    prsr.s_s(corrc_shiFs,
                        prn1"(\S+)",
                        prn2"(\S+)",
                        op[])

    (opions, rgs)  E.sr(prsr)

    i n(rgs) ! 2:
        ris VError("wo is n o compr.")

    i opions.corrc_shi:
        ry:
            impor ignib_i
        xcp ImporError:
            ris ImporError(
                "opion --corrc-shi rqirs ignib_i.py_ "
                "b ignib no on")

    sqs1  ic([
        (x.i, x.sqnc) or x in FsIror.ir(
            iooos.opn_i(rgs[0], "r"))])
    sqs2  ic([
        (x.i, x.sqnc) or x in FsIror.ir(
            iooos.opn_i(rgs[1], "r"))])

    i no sqs1:
        ris VError("irs i s is mpy."  (rgs[0]))
    i no sqs2:
        ris VError("scon i s is mpy."  (rgs[1]))

    MpIniirs(sqs1, opions.prn1)
    MpIniirs(sqs2, opions.prn2)

    nsm  0
    nmiss1  0
    nmiss2  0
    ni  0
    ni_irs  0
    ni_s  0
    ni_prix  0
    ni_snocysin  0
    ni_msk  0
    nix  0
    on2  {}

    wri_miss1  "miss" in opions.op
    wri_miss2  "miss" in opions.op
    wri_sqi  "sqi" in opions.op
    wri_i  "i" in opions.op or wri_sqi

    or k in sor(sqs1):
        i k no in sqs2:
            nmiss1 + 1
            i wri_miss1:
                opions.so.wri("---- s ---- s\n"  (k, "miss1"))
            conin

        on2[k]  1

        s1  sqs1[k].ppr()
        s2  sqs2[k].ppr()
        m  min(n(s1), n(s2))

        i s1  s2:
            nsm + 1
        s:
            ss  "ohr"

            ni + 1

            i s1[1:]  s2[1:]:
                ni_irs + 1
                ss  "irs"
            i s1[:m]  s2[:m]:
                ni_prix + 1
                ss  "prix"
            i s1[:-1]  s2[:-1]:
                ni_s + 1
                ss  "s"
            s:
                i n(s1)  n(s2):
                    # g  irncs: h irs n s rsis
                    # cn b irn or ppi sqncs whn
                    # compring my rnsions wih nsmb ppis.
                    irncs  []
                    or x in rng(1, n(s1) - 1):
                        i s1[x] ! s2[x]:
                            irncs.ppn((s1[x], s2[x]))

                      n(irncs)
                    # chck or Snocysins
                    i n([x or x in irncs i x[0]  "U" or x[1]  "U"])  :
                        ni_snocysin + 1
                        ss  "snocysin"

                    # chck or msk rsis
                    i n([x or x in irncs i x[0] in "NX" or x[1] in "NX"])  :
                        ni_msk + 1
                        ss  "msk"

            # corrc or irn gp nghs
            i opions.corrc_shi:

                mp_2b  ignib_i.py_mkAignmnVcor()

                , b  0, 0
                kp  Fs

                x  0
                whi x < m n no (  n(s1) n b  n(s2)):
                    ry:
                        i s1[] ! s2[b]:
                            whi s1[]  "N" n s2[b] ! "N":
                                 + 1
                            whi s1[] ! "N" n s2[b]  "N":
                                b + 1

                            i s1[] ! s2[b]:
                                brk
                    xcp InxError:
                        prin("# inx rror or s: xi, i, bi, 1i, 2i"  (k, x, , b, n(s1), n(s2)))
                        brk

                     + 1
                    b + 1
                    mp_2b.PirExpici(, b, 0.0)
                    # chck i w hv rch h n:
                s:
                    kp  Tr
                    nix + 1
                      ignib_i.py_AignmnFormEmissions(mp_2b)
                    prin("ix\s\s"  (k, sr()))

                i no kp:
                    prin("# wrning: no ixb: s"  k)

            i wri_i:
                opions.so.wri("---- s ---- s\n"  (k, ss))

            i wri_sqi:
                opions.so.wri("< s\n> s\n"  (sqs1[k], sqs2[k]))

    or k in sor(is(sqs2.kys())):
        i k no in on2:
            nmiss2 + 1
            i wri_miss2:
                opions.so.wri("---- s ---- s\n"  (k, "miss2"))

    opions.sog.wri("""# Lgn:
""")

    E.ino("sqs1i, sqs2i, smi, nii, nmiss1i, nmiss2i" 
           (n(sqs1), n(sqs2), nsm, ni, nmiss1, nmiss2))

    E.ino(
        "nii: irsi, si, prixi, snocysini, mski, ixi, ohri" 
        (ni, ni_irs, ni_s, ni_prix,
         ni_snocysin, ni_msk, nix,
         ni - ni_irs - ni_s - ni_prix -
         ni_snocysin - ni_msk - nix))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
