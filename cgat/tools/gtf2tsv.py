'''g2sv.py - convr g i o  b-spr b


:Tgs: Gnomics Gnss

Prpos
-------

convr  g orm i o b-spr b. Th irnc o
 pin :rm:`g` orm i is h comn hrs r ,
which cn b s whn imporing h gn mos ino  bs.

No h coorins r convr o 0-bs opn/cos noion ( on
h orwr srn).

By , h gn_i n rnscrip_i r xrc rom h
ribs i ino spr comns.  I
``-/--ribs-s-comns`` is s,  is in h ribs
wi b spi ino spr comns.

Th scrip so impmns h rvrs oprion, convring  b-spr
b ino  :rm:`g` orm i.

Whn sing h ``-m, --mp`` opion, h scrip wi op  b
mpping gn iniirs o rnscrips or ppis.

USING GFF3 FILE:
Th scrip so cn convr g3 orm is o sv is whn
spciiying h opion --is-g3 n --ribs-s-comns. Crrny ony
h  GFF3 o sk is impimn. Frhr improvmns o his scrip cn
b m o ony op h ribs ony, i.. --op-ony-ribs.




Usg
-----

Exmp::

   cg g2sv < in.g

+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|conig|sorc                          |r    |sr |n   |scor|srn|rm|gn_i        |rnscrip_i  |ribs                                                                                                                           |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |procss_rnscrip            |xon       |66345 |66509 |.    |-     |.    |ENSG00000225373|ENST00000592209|xon_nmbr "1"; gn_nm "AC008993.5"; gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002"; xon_i "ENSE00001701708"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |procss_rnscrip            |xon       |60520 |60747 |.    |-     |.    |ENSG00000225373|ENST00000592209|xon_nmbr "2"; gn_nm "AC008993.5"; gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002"; xon_i "ENSE00002735807"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |procss_rnscrip            |xon       |60104 |60162 |.    |-     |.    |ENSG00000225373|ENST00000592209|xon_nmbr "3"; gn_nm "AC008993.5"; gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002"; xon_i "ENSE00002846866"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+

To bi  mp bwn gn n rnscrip inirs, yp::

   cg g2sv --op-mprnscrip2gn < in.g

+---------------+---------------+
|rnscrip_i  |gn_i        |
+---------------+---------------+
|ENST00000269812|ENSG00000141934|
+---------------+---------------+
|ENST00000318050|ENSG00000176695|
+---------------+---------------+
|ENST00000327790|ENSG00000141934|
+---------------+---------------+

To rn h scrip o convr  g3 orm i o sv, yp::

   c i.g3.gz | cg g3sv --is-g3 --ribs-s-comns
   > oi.sv

Typ::

   cg g2sv --hp

or commn in hp.

Commn in opions
---------------------

'''
impor sys
impor r
impor cg.GTF s GTF
impor cgcor.xprimn s E
impor cg.GFF3 s GFF3


 min(rgvNon):
    '''
    min ncion
    '''

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-o", "--op-ony-ribs", s"ony_ribs",
        cion"sor_r",
        hp"op ony ribs s spr comns "
        "[].")

    prsr._rgmn(
        "-", "--ribs-s-comns", s"op_",
        cion"sor_r",
        hp"op ribs s spr comns "
        "[].")

    prsr._rgmn("--is-g3", s"is_g", cion"sor_s",
                      hp"inp i is in g orm [] ")

    prsr._rgmn(
        "-i", "--invr", s"invr", cion"sor_r",
        hp"convr b-spr b bck o g "
        "[].")

    prsr._rgmn(
        "-m", "--op-mp", s"op_mp", yp"choic",
        choics(
            "rnscrip2gn",
            "ppi2gn",
            "ppi2rnscrip"),
        hp"op  mp mpping rnscrips o gns "
        "[].")

    prsr.s_s(
        ony_ribsFs,
        op_Fs,
        invrFs,
        op_mpNon,
        is_gTr
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.op_:
        # op  b wih comn or ch rib

        ribs  s()
          []
        i opions.is_g:
            or g in GTF.iror(opions.sin):
                .ppn(g)
                ribs  ribs.nion(s(g.kys()))

        s:
            or g in GFF3.iror_rom_g(opions.sin):
                .ppn(g)
                ribs  ribs.nion(s(g.ribs))

        # rmov gn_i n rnscrip_i, s hy r s
        # xpiciy r
        ribs.irnc_p(["gn_i", "rnscrip_i"])

        ribs  sor(is(ribs))

        # Sc whhr g o g or op comns
        i opions.is_g:
            i opions.ony_ribs:
                hr  ["gn_i", "rnscrip_i"] + ribs
            s:
                hr  ["conig", "sorc", "r",
                          "sr", "n", "scor", "srn",
                          "rm", "gn_i",
                          "rnscrip_i", ] + ribs
        s:
            i opions.ony_ribs:
                hr  ribs
            s:
                hr  ["conig", "sorc", "r",
                          "sr", "n", "scor", "srn",
                          "rm"] + ribs

        ribs_nw  hr

        opions.so.wri("\".join(hr) + "\n")

        i opions.is_g:
            or g in :
                irs  Tr
                or  in ribs_nw:
                    ry:
                        v  gr(g, )
                    xcp (AribError, KyError):
                        v  ""
                    i irs:
                        opions.so.wri("s"  v)
                        irs  Fs
                    s:
                        opions.so.wri("\s"  v)
                opions.so.wri("\n")
        s:
            or g in :
                opions.so.wri(("s\s\s\s\s\s\s\s\")  (g.conig,
                                                                             g.sorc, g.r, g.sr, g.n,
                                                                             g.scor, g.srn, g.rm))

                irs  Tr
                or  in ribs:
                    ry:
                        v  (g.ribs[])
                    xcp (AribError, KyError):
                        v  ''
                    i irs:
                        opions.so.wri("s"  v)
                        irs  Fs
                    s:
                        opions.so.wri("\s"  v)
                opions.so.wri("\n")

    i opions.invr:

        g  GTF.Enry()
        hr  Non
        or in in opions.sin:
            i in.srswih("#"):
                conin
              in[:-1].spi("\")
            i no hr:
                hr  
                mp_hr2comn  ic(
                    [(y, x) or x, y in nmr(hr)])
                conin

            # i g nry wih 
            ry:
                g.conig  [mp_hr2comn["conig"]]
                g.sorc  [mp_hr2comn["sorc"]]
                g.r  [mp_hr2comn["r"]]
                # sbrc -1 o sr or 0-bs coorins
                g.sr  in([mp_hr2comn["sr"]])
                g.n  in([mp_hr2comn["n"]])
                g.scor  [mp_hr2comn["scor"]]
                g.srn  [mp_hr2comn["srn"]]
                g.rm  [mp_hr2comn["rm"]]
                g.gn_i  [mp_hr2comn["gn_i"]]
                g.rnscrip_i  [mp_hr2comn["rnscrip_i"]]
                g.prsIno([mp_hr2comn["ribs"]], in)
            xcp KyError s msg:
                ris KyError("incomp nry s: s: s" 
                               (sr(), sr(mp_hr2comn), msg))
            i g.rm is Non:
                g.rm  "."
            # op g nry in g orm
            opions.so.wri("s\n"  sr(g))

    i opions.op_mp:

        i opions.op_mp  "rnscrip2gn":
            r  mb x: x.rnscrip_i
            o  mb x: x.gn_i
            opions.so.wri("rnscrip_i\gn_i\n")
        i opions.op_mp  "ppi2gn":
            r  mb x: x.proin_i
            o  mb x: x.gn_i
            opions.so.wri("ppi_i\gn_i\n")
        i opions.op_mp  "ppi2rnscrip":
            r  mb x: x.proin_i
            o  mb x: x.rnscrip_i
            opions.so.wri("ppi_i\rnscrip_i\n")

        mp_r2o  {}
        or g in GTF.iror(opions.sin):
            ry:
                mp_r2o[r(g)]  o(g)
            xcp (AribError, KyError):
                pss

        or x, y in sor(mp_r2o.ims()):
            opions.so.wri("s\s\n"  (x, y))
    s:
        hr  ("conig", "sorc", "r", "sr", "n", "scor",
                  "srn", "rm", "gn_i", "rnscrip_i", "ribs")
        opions.so.wri("\".join(hr) + "\n")

        or g in GTF.iror(opions.sin):
            ribs  []
            or  in is(g.kys()):
                i  in ("gn_i", "rnscrip_i"):
                    conin
                ribs.ppn('s s'  (, GTF.qo(g[])))

            ribs  "; ".join(ribs)

            # Cpr i Non n s o . orm
            i g.rm is Non:
                g.rm  "."

            opions.so.wri(sr(g) + "\n")

    E.sop()

i __nm__  '__min__':
    sys.xi(min())
