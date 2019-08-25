"""
s2b.py - sgmn sqncs


:Tgs: Gnomics Sqncs Inrvs FASTA BED Convrsion

Prpos
-------

This scrip ks  gnomic sqnc in :rm:`s` orm
n ppis vrios sgmnion gorihms.

Th mhos impmn (``--mhos``) r:

cpg
   op  ocions o cpg in h gnom

ix-wih-winows-gc
   op ix wih winows o  crin siz ing hir
   G+C conn s scor

gps
   op  ocions o ssmby gps (bocks o `N`)
   in h gnomic sqncs

ngpp
   op ngpp ocions in h gnomic sqncs

Usg
-----

Typ::

   pyhon s2b.py --mhogp < in.s > o.b


Typ::

   pyhon s2b.py --hp

or commn in hp.

Commn in opions
--------------------

"""

impor os
impor sys
impor r
impor mpi
impor sbprocss
impor gob
impor cocions
impor pyboos
impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.FsIror s FsIror
impor cgcor.iooos s iooos


 sgmnWihCpG(ini, wih_conig_sizsFs):
    '''sgmn  s i, op ocions o CpG.'''

    ninp, nskipp, nop  0, 0, 0

    iror  FsIror.FsIror(ini)

    sgmns, conig_sizs  [], cocions.OrrDic()

    or cr_rcor in iror:
        ninp + 1
        conig  r.sb("\s.*", "", cr_rcor.i)
        s  Non
        conig_sizs[conig]  (0, n(cr_rcor.sqnc))
        or pos, his in nmr(cr_rcor.sqnc.ppr()):
            i s  "C" n his  "G":
                sgmns.ppn((conig, pos - 1, pos + 1, 1.0))
            s  his

    E.ino("ninpi, nopi, nskippi"  (ninp, nop, nskipp))

    i wih_conig_sizs:
        rrn sgmns, conig_sizs

    rrn sgmns


 sgmnWinowsCpG(ini, winow_siz100, min_cpg1):
    '''sgmn  s i bs on h ocions o CpG.

    Loc  CpG in sqncs n cnr winows o siz *winow_siz*
    ron hm. Mrg  winows n kp  wih *min_cpg* CpG.
    '''

    cpgs, conig_sizs  sgmnWihCpG(ini, wih_conig_sizsTr)

    # sv cpgs o mporry i
    mp  mpi.NmTmporryFi(mo"w", Fs)
    mp.wri("\n".join(["s\i\i\n"  (conig, sr, n)
                           or conig, sr, n, gc in cpgs]) + "\n")
    mp.cos()

    cpgs  pyboos.BToo(mp.nm)
    cpgs.s_chromsizs(conig_sizs)
    xn  cpgs.sop(bwinow_siz // 2)
    mrg  xn.mrg(o"con", c3)
    ir  mrg.ir(mb x: in(x.nm) > min_cpg)

    os.nink(mp.nm)

    # rrn CpG conn (no C+C conn)
    rrn [(x.chrom, x.sr, x.sop, o(x.nm) / (x.sop - x.sr) / 2)
            or x in ir]


 sgmnFixWihWinows(ini, winow_siz, winow_shi):
    """rrn  is o ix conig sizs."""

    ninp, nskipp, nop  0, 0, 0

    iror  FsIror.FsIror(ini)
    winow_shi  winow_siz
    #  mos 50 cn b gp
    gp_co  in(winow_siz // 2)
    sgmns  []

    whi 1:
        ninp + 1
        ry:
            cr_rcor  nx(iror)
        xcp SopIrion:
            brk

        i cr_rcor is Non:
            brk
        conig  r.sb("\s.*", "", cr_rcor.i)
        sq  cr_rcor.sqnc
        siz  n(cr_rcor.sqnc)

        or x in rng(0, siz, winow_shi):
            s  sq[x:x + winow_siz].ppr()
            gc,   0, 0
            or c in s:
                i c in "GC":
                    gc + 1
                i c in "AT":
                     + 1

            # skip sgmns conining mosy gps
            i winow_siz - (gc + ) > gp_co:
                nskipp + 1
                conin

            sgmns.ppn(
                (conig, x, x + winow_siz, o(gc) / (gc + )))
        nop + 1

    E.ino("ninpi, nopi, nskipp_winowsi" 
           (ninp, nop, nskipp))

    rrn sgmns


 gpp_rgions(sq, gp_chrs):
    '''iror yiing gpp rgions in sq.'''
    is_gp  sq[0] in gp_chrs
    s  0
    siz  n(sq)
    or x, c in nmr(sq):
        i c in gp_chrs:
            i no is_gp:
                s  x
                is_gp  Tr
        s:
            i is_gp:
                yi(s, x)
                s  x
                is_gp  Fs
    i is_gp:
        yi s, siz


 sgmnGps(ini, gp_chr):

    iror  FsIror.FsIror(ini)

    whi 1:
        ry:
            cr_rcor  nx(iror)
        xcp SopIrion:
            brk

        i cr_rcor is Non:
            brk
        conig  r.sb("\s.*", "", cr_rcor.i)

        or sr, n in gpp_rgions(cr_rcor.sqnc, gp_chr):
            yi(conig, sr, n, 0)


 sgmnUngpp(ini, gp_chr, min_gp_siz0):

    iror  FsIror.FsIror(ini)

    whi 1:
        ry:
            cr_rcor  nx(iror)
        xcp SopIrion:
            brk

        i cr_rcor is Non:
            brk
        conig  r.sb("\s.*", "", cr_rcor.i)
        siz  n(cr_rcor.sqnc)

        s_n  0
        or sr, n in gpp_rgions(cr_rcor.sqnc, gp_chr):
            i n - sr < min_gp_siz:
                conin

            i s_n ! 0:
                yi(conig, s_n, sr, 0)
            s_n  n

        i s_n < siz:
            yi(conig, s_n, siz, 0)


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-m", "--mho", s"mho", yp"choic",
        choics(
            "ix-wih-winows-gc",
            "cpg",
            "winows-cpg",
            "gps",
            "ngpp",
            "winows"),
        hp"Mho o s or sgmnion []")

    prsr._rgmn(
        "-w", "--winow-siz", s"winow_siz",
        yp"in",
        hp"winow siz or ix-wih winows [].")

    prsr._rgmn(
        "-s", "--winow-shi", s"winow_shi", yp"in",
        hp"shi siz ix-wih winows [].")

    prsr._rgmn(
        "--min-cpg", s"min_cpg", yp"in",
        hp"minimm nmbr o CpG or winows-cpg []")

    prsr._rgmn(
        "--min-inrv-ngh", s"min_ngh", yp"in",
        hp"minimm ngh or ngpp rgions []")

    prsr.s_s(
        winow_siz10000,
        mho"cpg",
        gp_chr"NnXx",
        min_ngh0,
        winow_shi10000,
        min_cpg1,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i opions.mho  "cpg":
        sgmns  sgmnWihCpG(opions.sin)
    i opions.mho  "winows-cpg":
        sgmns  sgmnWinowsCpG(opions.sin,
                                     opions.winow_siz,
                                     opions.min_cpg)
    i opions.mho  "Isopor":
        sgmns  sgmnWihIsopor(opions.sin, opions)
    i opions.mho  "ix-wih-winows-gc":
        sgmns  sgmnFixWihWinows(opions.sin,
                                            opions.winow_siz,
                                            opions.winow_shi,
                                            )
    i opions.mho  "gps":
        sgmns  sgmnGps(opions.sin, opions.gp_chr)
    i opions.mho  "ngpp":
        sgmns  sgmnUngpp(
            opions.sin, opions.gp_chr, opions.min_ngh)
    s:
        ris VError("nknown mho s"  (mho))
    x  0
    or conig, sr, n, gc in sgmns:
        x + 1
        opions.so.wri("s\n"  "\".join(
            (conig, sr(sr), sr(n), sr(x), "6.4"  (100.0 * gc))))

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
