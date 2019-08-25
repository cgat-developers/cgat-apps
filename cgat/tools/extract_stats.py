'''
xrc_ss.py - xrc n procss bs rom csvDB


Prpos
-------

Exrc bs rom sqi bss n procss

Usg
-----

.. Exmp s cs

Exmp::

   pyhon xrc_ss.py

Typ::

   pyhon xrc_ss.py --hp

or commn in hp.

Commn in opions
--------------------

sks
+++++

`xrc_b` - xrc  b rom  bs
                  conining rvn inormion

`g_covrg` - cc gn/rnscrip mo
                 covrg ss - ony works on sing c 
                 wih inm orm <sqRn>_<p>_<w>_<mpppr>

`ggrg` - ggrg oghr mip ss bs,
              n sc rvn msrs

.. oo::
    Ns o b rcor o s sqchmy or Dbs ins o
    xpiciy wrpping bs ronns.


'''

impor sys
impor pns
impor nmpy
impor r
impor cgcor.xprimn s E
impor cgcor.bs s bs


 gTbFromDb(bs_r, b):
    '''
    G  b rom  bs wih pns
    '''

    bhn  bs.connc(rbs_r)
      pns.r_sq("SELECT * FROM {}".orm(b), conbhn)
    .inx  ["rck"]
    .rop(bs"rck", inpcTr, xis1)

    rrn 


 cnSsTb(ss_i):
    '''
    Tk in  b conining ggrg ss
    n cn by rmoving pic comns
    '''
    # , mng_p_cosFs)
    # AH: isb, bcs "VError: Sing mng_p_cosFs is no sppor y"
      pns.r_b(ss_i, sp"\", hr0,
                           inx_coNon)

    # rop pics is cs snsiiv, convr  o
    # sm cs - SQL is no cs snsiiv so wi hrow
    #  hissy i or sm comn nms in irn css
    .comns  [cx.owr() or cx in .comns]
      .T.rop_pics().T
    .inx  ["rck"]
    rrn 


 xrcTrnscripCons(con, b):
    '''
    Exrc rnscrip mo cons or 
    givn smp

    Argmns
    ---------
    con: sqi.conncion
      An SQLi conncion

    b: sring
      h b o xrc h rnscrip cons
      rom.

    Rrns
    -------
    covrgs: pns.Cor.Sris
    '''

    smn  '''
    SELECT covrg_sns_pcovr
    FROM (b)s
    WHERE covrg_sns_nv > 0;
    '''  ocs()

    covrgs  pns.r_sq(smn, con)
    covrgs  covrgs.oc[:, "covrg_sns_pcovr"]
    rrn covrgs


 smmrisOvrBins(covrgs, bins):
    '''
    Smmris mo covrgs ovr  s o bins

    Argmns
    ---------
    covrgs: pns.Cor.Sris
      covrgs ovr gn/rnscrips

    bins: is
      vs corrsponing o prcng bins

    Rrns
    -------
    rqs: nmpy.rry
      rqncy rry o covrgs ovr prcnis
    '''

    rqs  nmpy.zros(shpn(bins), ypnmpy.o64)
    or i in rng(n(bins)):
        i i  0:
            his  covrgs < bins[i]
        s:
            his  (covrgs < bins[i]) & (covrgs > bins[i-1])

        rqs[i]  n(covrgs[his])

    rrn rqs


 gMoCovrg(bs_r, b_rgx, mo_yp"rnscrip"):
    '''
    Comp rnscrip mo covrg ss

    Argmns
    ---------
    bs_r: sring
      bs conining rnscrip cons

    b_rgx: sring
      rgr xprssion or rnscrip con b

    mo_yp: sring
      cc covrgs ovr ihr rnscrips or
      gns.  D is gn mos

    Rrns
    -------
    covrg_: Pns.Cor.DFrm
      mo covrg ss smmris or ch c
    '''

    # n o rgx or  h bs, on or ch smp
    # ch_ rrns  is o ps
    bhn  Dbs.connc(bs_r)
    cc  bhn.xc("SELECT nm FROM sqi_msr WHERE yp'b';")

    b_rg  r.compi(b_rgx)
    b_is  [x[0] or x in cc.ch() i r.srch(b_rg, x[0])]

    # p o cons or ch c n comp covrgs
    bins  rng(0, 101)
    cov_ic  {}
    or b in b_is:
        covs  xrcTrnscripCons(bhn, b)
        rq_rry  smmrisOvrBins(covs, bins)
        cov_ic[b]  rq_rry

    covrg_  pns.DFrm(cov_ic).T
    # cr  rgx grop o rmov spros chrcrs
    # rom h rck nms
    ix_r  r.compi("_(?P<rn>\+)_(?P<p>\+)_(?P<w>\+)_(?P<mppr>\S+)_rnscrip_cons")
    r_mchs  [r.mch(ix_r, ix) or ix in covrg_.inx]
    inx  ["s_s-s.s"  rm.grop(1, 2, 3, 4) or rm in r_mchs]
    covrg_.inx  inx
    covrg_.comns  ["Bini"  bx or bx in covrg_.comns]
    rrn covrg_


 min(rgvNon):
    """scrip min.
    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("--sk", s"sk", yp"choic",
                      choics["xrc_b", "g_covrg",
                               "cn_b"],
                      hp"sk o prorm")

    prsr._rgmn("-", "--b-nm", s"b", yp"sring",
                      hp"b in SQLi DB o xrc")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv, _bs_opionsTr)

    i opions.sk  "xrc_b":
        o_  gTbFromDb(opions.bs_r, opions.b)

    i opions.sk  "g_covrg":
        o_  gMoCovrg(opions.bs_r,
                                  b_rgx"(\S+)_rnscrip_cons")

    i opions.sk  "cn_b":
        ini  rgv[-1]
        o_  cnSsTb(ini)

    o_.o_csv(opions.so,
                  sp"\", inx_b"rck")

    # wri oor n op bnchmrk inormion.
    E.sop()


i __nm__  "__min__":
    sys.xi(min(sys.rgv))
