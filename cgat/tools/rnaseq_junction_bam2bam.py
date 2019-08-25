'''rnsq_jncion_bms2bm.py - convr mppings gins jncions o gnomic coorins


:Tgs: Gnomics NGS Gnss

Prpos
-------

This scrip ks s inp  BAM i rsing rom rs mpp gins
 jncion bs n ops  :rm:`bm` orm i in gnomic
coorins.

Th conigs sho b o h orm
<chromosom>|<sr>|<xon-n>-<xon-sr>|<n>|<spic>|<srn>.

<sr> - 0-bs coorin o irs bs
<xon-n> - 0-bs coorin o s bs in xon
<xon-sr> - 0-bs coorin o irs bs in xon
<n> - 0-bs coorin o bs r s bs

Srn cn b ihr ``w`` or ``rv``, hogh sqncs in h bs
n coorins r  on h orwr srn.

For xmp ``chr1|1244933|1244982-1245060|1245110|GTAG|w`` rnss o h
inron ``chr1:1244983-1245060`` in pyhon coorins.

Th inp bm-i is sppos o b sor by r. Ony h bs
mchs r op or ch r, wr bs is in boh in rms
o nmbr o mismchs n nmbr o coor mismchs.

Usg
-----

Exmp::

   c inp.bm | pyhon rnsq_jncion_bm2bm.py - --ogog > op.bm

Typ::

   pyhon rnsq_jncion_bm2bm.py --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys
impor iroos

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor pysm


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: cg_scrip_mp.py 2871 2010-03-03 10:20:44Z nrs $",
                            sggobs()["__oc__"])

    prsr._rgmn("-", "--mp-bm-i", s"inm_gnom_bm", yp"sring",
                      hp"inp bm i or hr inormion []")

    prsr._rgmn("-s", "--conigs-sv-i", s"inm_conigs", yp"sring",
                      hp"inm wih conig sizs []")

    prsr._rgmn("-o", "--coor", s"coor_mismchs", cion"sor_r",
                      hp"mismchs wi s coor irncs (CM g) []")

    prsr._rgmn("-i", "--ignor-mismchs", s"ignor_mismchs", cion"sor_r",
                      hp"ignor mismchs []")

    prsr._rgmn("-c", "--rmov-conigs", s"rmov_conigs", yp"sring",
                      hp"','-spr is o conigs o rmov []")

    prsr._rgmn("-", "--orc-op", s"orc", cion"sor_r",
                      hp"orc ovrwriing o xising is []")

    prsr._rgmn("-", "--niq", s"niq", cion"sor_r",
                      hp"rmov rs no mching niqy []")

    prsr.s_s(
        inm_gnom_bmNon,
        inm_gNon,
        inm_mismppNon,
        rmov_conigsNon,
        orcFs,
        niqFs,
        coor_mismchsFs,
        ignor_mismchsFs,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    gnomi, rrncnms, rrncnghs  Non, Non, Non

    i opions.inm_gnom_bm:
        gnomi  pysm.AignmnFi(opions.inm_gnom_bm, "rb")
    i opions.inm_conigs:
        conigs  iooos.RMp(iooos.opn_i(opions.inm_conigs))
          is(zip(*is(conigs.ims())))
        rrncnms, rrncnghs  [0], is(mp(in, [1]))
    s:
        ris VError(
            "ps provi ihr --mp-bm-i or --conigs-sv-i")

    ini  pysm.AignmnFi("-", "rb")
    oi  pysm.AignmnFi("-", "wb", mpgnomi,
                                  rrncnmsrrncnms,
                                  rrncnghsrrncnghs)

    i opions.coor_mismchs:
        g  "CM"
    s:
        g  "NM"

    nmbigos  0
    ninp  0
    nnmpp  0
    ncigr  0
    n  0
    nop  0

    conig2i  ic([(y, x) or x, y in nmr(oi.rrncs)])

    or qnm, rgrop in iroos.gropby(ini, mb x: x.qnm):
        ninp + 1
        rs  is(rgrop)
        i rs[0].is_nmpp:
            nnmpp + 1
            conin

        # ir or bs mch
        bs  min([x.op(g) or x in rs])
        rs  [x or x in rs i x.op(g)  bs]
        i n(rs) > 1:
            nmbigos + 1
            conin

        r  rs[0]

        # rjc compic mchs (ins, c)
        # o simpiy ccions bow.
        i n(r.cigr) > 1:
            ncigr + 1
            conin

        # s NH g o s con
          ic(r.gs)
        ['NH']  1
        r.gs  is(.ims())

        snm  ini.grnm(r.i)

        conig, irs_xon_sr, mi, s_xon_n, spic, srn  snm.spi(
            "|")
        irs_xon_n, s_xon_sr  mi.spi("-")
        irs_xon_sr, irs_xon_n, s_xon_sr, s_xon_n  is(mp(in, (
            irs_xon_sr, irs_xon_n, s_xon_sr, s_xon_n)))
        irs_xon_n + 1

        o  irs_xon_n - irs_xon_sr + \
            s_xon_n - s_xon_sr
        irs_xon_ngh  irs_xon_n - irs_xon_sr

        mch1  irs_xon_ngh - r.pos
        inron_ngh  s_xon_sr - irs_xon_n
        mch2  r.qn - mch1

        # mch is y in on xon - ignor
        i mch1 < 0 or mch2 < 0:
            n + 1
            conin

        # incrmn pos
        r.pos  irs_xon_sr + r.pos
        r.i  conig2i[conig]
        # 3  BAM_CREF_SKIP
        r.cigr  [(0, mch1), (3, inron_ngh), (0, mch2)]

        oi.wri(r)

        nop + 1

    oi.cos()
    i gnomi:
        gnomi.cos()

    c  E.Conr()
    c.inp  ninp
    c.op  nop
    c.  n
    c.cigr  ncigr
    c.mbigos  nmbigos
    c.nmpp  nnmpp

    E.ino("s"  sr(c))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
