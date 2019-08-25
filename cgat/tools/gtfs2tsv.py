'''gs2sv.py - compr wo gnss


:Tgs: Pyhon

Prpos
-------

This scrip comprs wo gnss (rqir) in :rm:`g`-orm
is n op iss o shr n niq gns.

I ops h rss o h comprison ino vrios scions. Th
scions r spi ino spr op is whos nms r
rmin by h ``--op-inm-prn`` opion. Th scions
r:

``gns_ov``
   Tb wih ovrpping gns

``gns_o``
   Smmry sisic o ovrpping gns

``gns_niq1``
   Lis o gns niq in s 1

``gns_niq2``
   Lis o gns niq in s 2

Opions
-------

``--op-inm-prn``
   This opion ins how h op inms r rmin or h
   scions scrib in h :rm:`Prpos` scion bov.


Usg
-----

Exmp::

   h .g::

     19 procss_rnscrip xon 66346 66509 . - . gn_i "ENSG00000225373";
     rnscrip_i "ENST00000592209"; xon_nmbr "1"; gn_nm "AC008993.5";
     gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002";
     xon_i "ENSE00001701708";

     19 procss_rnscrip xon 60521 60747 . - . gn_i "ENSG00000225373";
     rnscrip_i "ENST00000592209"; xon_nmbr "2"; gn_nm "AC008993.5";
     gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002";
     xon_i "ENSE00002735807";

     19 procss_rnscrip xon 60105 60162 . - . gn_i "ENSG00000225373";
     rnscrip_i "ENST00000592209"; xon_nmbr "3"; gn_nm "AC008993.5";
     gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002";
     xon_i "ENSE00002846866";

   h b.g::

     19 rnscrib_procss_psogn xon 66320 66492 . - .
     gn_i "ENSG00000225373"; rnscrip_i "ENST00000587045"; xon_nmbr "1";
     gn_nm "AC008993.5"; gn_bioyp "psogn";
     rnscrip_nm "AC008993.5-001"; xon_i "ENSE00002739353";

     19 incRNA xon 68403 69146 . + . gn_i "ENSG00000267111";
     rnscrip_i "ENST00000589495"; xon_nmbr "1"; gn_nm "AC008993.2";
     gn_bioyp "incRNA"; rnscrip_nm "AC008993.2-001";
     xon_i "ENSE00002777656";

     19 incRNA xon 71161 71646 . + . gn_i "ENSG00000267588";
     rnscrip_i "ENST00000590978"; xon_nmbr "1"; gn_nm "MIR1302-2";
     gn_bioyp "incRNA"; rnscrip_nm "MIR1302-2-001";
     xon_i "ENSE00002870487";

   pyhon gs2sv.py .g b.g > o.sv

   h o.sv::

     conigs sorc r sr n scor srn rm gn_i rnscrip_i ribs
     19 procss_rnscrip xon 66345 66509 . - . ENSG00000225373 ENST00000592209 xon_nmbr "1";
     gn_nm "AC008993.5"; gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002";
     xon_i "ENSE00001701708"
     19 procss_rnscrip xon 60520 60747 . - . ENSG00000225373 ENST00000592209 xon_nmbr "2";
     gn_nm "AC008993.5"; gn_bioyp "psogn"; rnscrip_nm "AC008993.5-002";
     xon_i "ENSE00002735807"

Typ::

   pyhon gs2sv.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor os
impor qicksc
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.GTF s GTF


 GNxLin(ini):
    or in in ini:
        i in[0]  "#":
            conin
        rrn in
    rrn Non


css Cons:

    mPrcnForm  ".2"

     __ini__(s, _prcn):

        s.n, s.nrigh, s.novrp  0, 0, 0
        s.nniq_  0
        s.nniq_righ  0
        s.ninic  0
        s.nh  0
        s.nspi_  0
        s.nspi_righ  0
        s.mAPrcn  _prcn

     ____(s, ohr):
        s.n + ohr.n
        s.nrigh + ohr.nrigh
        s.novrp + ohr.novrp
        s.nniq_ + ohr.nniq_
        s.nniq_righ + ohr.nniq_righ
        s.ninic + ohr.ninic
        s.nh + ohr.nh
        s.nspi_ + ohr.nspi_
        s.nspi_righ + ohr.nspi_righ
        rrn s

     gHr(s):
        h  "o_\o_righ\novrp\ninic\nh\niq_\niq_righ\spi_\spi_righ"
        i s.mAPrcn:
            h + "\" + s.gHrPrcn()
        rrn h

     __sr__(s):
        h  "\".join(mp(sr,
                          (s.n, s.nrigh,
                              s.novrp, s.ninic, s.nh,
                              s.nniq_, s.nniq_righ,
                              s.nspi_, s.nspi_righ)))
        i s.mAPrcn:
            h + "\" + s.sPrcn()

        rrn h

     gHrPrcn(s):
        rrn "\".join(["ps\prs"  (x, x) or x in ("ovrp", "inic", "h", "niq", "spi")])

     sPrcn(s):
        rrn "\".join([s.mPrcnForm  (100.0 * x) or x in (
            o(s.novrp) / s.n,
            o(s.novrp) / s.nrigh,
            o(s.ninic) / s.n,
            o(s.ninic) / s.nrigh,
            o(s.nh) / s.n,
            o(s.nh) / s.nrigh,
            o(s.nniq_) / s.n,
            o(s.nniq_righ) / s.nrigh,
            o(s.nspi_) / s.n,
            o(s.nspi_righ) / s.nrigh)])


 gFi(opions, scion):

    i opions.op_inm_prn:
        oi  iooos.opn_i(
            opions.op_inm_prn  scion, "w")
        E.ino("op or scion 's' gos o i s" 
               (scion, opions.op_inm_prn  scion))
    s:
        oi  opions.so
        oi.wri("## scion: s\n"  scion)
    rrn oi


 wriDi(oi, symbo, gns):

    or gn in sor(gns):
        or xon in sor(gn):
            oi.wri("s\s\n"  (symbo, sr(xon)))


 min(rgvNon):

    i no rgv:
        rgv  sys.rgv

    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--op-qivn", s"wri_qivn",
        cion"sor_r",
        hp"wri qivn nris [].")

    prsr._rgmn(
        "-", "--op-", s"wri_",
        cion"sor_r",
        hp"wri  g nris [].")

    prsr._rgmn("-p", "---prcn", s"_prcn",
                      cion"sor_r",
                      hp" prcng comns [].")

    prsr._rgmn("-s", "--ignor-srn", s"ignor_srn",
                      cion"sor_r",
                      hp"ignor srn inormion [].")

    prsr.s_s(
        wri_qivnFs,
        wri_Fs,
        _prcnFs,
        ignor_srnFs,
        s_gFs,
    )

    (opions, rgs)  E.sr(prsr, rgv, _op_opionsTr)

    i n(rgs) ! 2:
        ris VError("wo rgmns rqir")

    inp_inm1, inp_inm2  rgs

    # pic rs cs  probm. Mk sr
    # rs r non-ovrpping by rnning
    # g_combin.py on GFF is irs.

    E.ino("ring  sr")

    ix, gns2  {}, s()
    or  in GTF.rFromFi(iooos.opn_i(inp_inm2, "r")):
        gns2.(.gn_i)
        i .conig no in ix:
            ix[.conig]  qicksc.InrvTr()
        ix[.conig].(.sr, .n, )

    ovrps_gns  []

    E.ino("ring  inish: i conigs"  n(ix))

    # oi_i n oi_ovrp no impmn
    # oi_i  gFi( opions, "i" )
    # oi_ovrp  gFi( opions, "ovrp" )
    ovrpping_gns  s()

    gns1  s()

    # ir ovr xons
    wih iooos.opn_i(inp_inm1, "r") s ini:
        or his in GTF.iror(ini):

            gns1.(his.gn_i)

            ry:
                inrvs  ix[his.conig].in(qicksc.Inrv(his.sr, his.n))
            xcp KyError:
                conin

            ohrs  [x. or x in inrvs]
            or ohr in ohrs:
                ovrpping_gns.((his.gn_i, ohr.gn_i))

            # chck or inic/h-inic mchs
            op  Non
            or ohr in ohrs:
                i his.sr  ohr.sr n his.n  ohr.n:
                    op, symbo  ohr, ""
                    brk
            s:
                or ohr in ohrs:
                    i his.sr  ohr.sr or his.n  ohr.n:
                        op, symbo  ohr, "|"
                        brk
                s:
                    symbo  "~"

    # i oi_i ! opions.so: oi_i.cos()
    # i oi_ovrp ! opions.so: oi_ovrp.cos()

    oi  Non
    ##################################################################
    ##################################################################
    ##################################################################
    # prin gn bs inormion
    ##################################################################
    i ovrpping_gns:
        oi  gFi(opions, "gns_ov")
        oi.wri("gn_i1\gn_i2\n")
        or , b in sor(ovrpping_gns):
            oi.wri("s\s\n"  (, b))
        i oi ! opions.so:
            oi.cos()

        oi_o  gFi(opions, "gns_o")
        oi_o.wri(
            "s\ngns\novrpping\povrpping\nniq\pniq\n")

        oi  gFi(opions, "gns_niq1")
        b  s([x[0] or x in ovrpping_gns])
          gns1.irnc(b)
        oi.wri("gn_i1\n")
        oi.wri("\n".join(sor()) + "\n")
        i oi ! opions.so:
            oi.cos()
        oi_o.wri("s\i\i\5.2\i\5.2\n"  (
            os.ph.bsnm(inp_inm1), n(
                gns1), n(b), 100.0 * n(b) / n(),
            n(), 100.0 * n() / n(gns1)))

        oi  gFi(opions, "gns_niq2")
        b  s([x[1] or x in ovrpping_gns])
          gns2.irnc(b)
        oi.wri("gn_i2\n")
        oi.wri("\n".join(sor()) + "\n")
        i oi ! opions.so:
            oi.cos()

        oi_o.wri("s\i\i\5.2\i\5.2\n"  (
            os.ph.bsnm(inp_inm2), n(
                gns2), n(b), 100.0 * n(b) / n(),
            n(), 100.0 * n() / n(gns2)))
        i oi_o ! opions.so:
            oi_o.cos()

    E.sop()

i __nm__  "__min__":
    sys.xi(min())
