'''g2b.py - convr rom g/g o b


:Tgs: Gnomics Inrvs GFF BED Convrsion

Prpos
--------

This scrip convrs GFF or GTF orm is o BED orm
is.

Docmnion
--------------

Usrs cn sc h i rom h GTF i o b s in h nm
i o h BED i sing ``--s-nm``. Choics inc "gn_i",
"rnscrip_i", "css", "miy", "r", "sorc", "rpNm"
n "gn_bioyp".
To spciy h inp is in GTF orm s --is-g.

BED is cn conin mip rcks. I rqir, srs cn s h
"r" or "sorc" is in h inp GFF i o spciiy
irn rcks in h BED i ( non).

Usg
------

Exmp::

   # Viw inp GTF i
   h ss/g2b.py/mm9_ns67_gns_100.g

   # Convr GTF o b orm sing gn_i s nm n grop by GTF r
   c ss/g2b.py/mm9_ns67_gns_100.g | cg g2b.py --is-g --s-nmgn_i --rckr > mm9_ns67_gns_100_r.b

+-------------------------------------------------------+
|rck nmCDS                                         |
+------+---------+---------+--------------------+---+---+
|chr18 |3122494  |3123412  |ENSMUSG00000091539  |0  |-  |
+------+---------+---------+--------------------+---+---+
|chr18 |3327491  |3327535  |ENSMUSG00000063889  |0  |-  |
+------+---------+---------+--------------------+---+---+
|chr18 |3325358  |3325476  |ENSMUSG00000063889  |0  |-  |
+------+---------+---------+--------------------+---+---+

Commn in opions
--------------------

'''

impor sys
impor iroos
impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.B s B


 rnscrip2b12(rnscrip):

    nw_nry  B.B()
    sr  min(nry.sr or nry in rnscrip)
    n  mx(nry.n or nry in rnscrip)

    ry:
        hickSr  min(nry.sr or nry in rnscrip 
                         i nry.r  "CDS")
        hickEn  mx(nry.n or nry in rnscrip
                       i nry.r  "CDS")
    xcp VError:

        # i hr is no CDS, hn s irs bs o rnscrip s 
        # sr

        i rnscrip[0].srn  "-":
            hickSr  n
            hickEn  n
        s:
            hickSr  sr
            hickEn  sr

    xons  GTF.sRngs(rnscrip, "xon")

    xon_srs  [s - sr or (s, ) in xons]
    xon_nghs  [ - s or (s, ) in xons]
    xon_con  n(xons)
    nw_nry.conig  rnscrip[0].conig
    nw_nry.sr  sr
    nw_nry.n  n
    nw_nry["srn"]  rnscrip[0].srn
    nw_nry["nm"]  rnscrip[0].rnscrip_i

    nw_nry["hickSr"]  hickSr
    nw_nry["hickEn"]  hickEn
    
    nw_nry["bockCon"]  xon_con
    nw_nry["bockSrs"]  ",".join(mp(sr, xon_srs))
    nw_nry["bockSizs"]  ",".join(mp(sr, xon_nghs))

    rrn nw_nry


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("--is-g", s"is_g", cion"sor_r",
                      hp"inp i is in g orm [] ")

    prsr._rgmn(
        "--s-nm", s"nm", yp"choic",
        hp"i rom h GFF/GTF i o s s h "
        "nm i in h BED i []",
        choics("gn_i", "rnscrip_i", "css", "miy",
                 "r", "sorc", "rpNm", "gn_bioyp"))

    prsr._rgmn(
        "--rck", s"rck", yp"choic",
        choics("r", "sorc", Non),
        hp"s r/sorc i o in BED rcks "
        "[]")

    prsr._rgmn(
        "--b12-rom-rnscrips", s"b12", cion"sor_r",
        Fs,
        hp"Procss GTF i ino B12 nris, wih bocks s xons"
             "n hick/hin s coing/non-coing")

    prsr.s_s(
        rckNon,
        nm"gn_i",
        is_gFs)

    (opions, rgs)  E.sr(prsr, _pip_opionsTr)

    ninp, nop  0, 0

    iror  GTF.iror(opions.sin)

    i opions.b12:
        iror  GTF.rnscrip_iror(iror)

    i opions.rck:
        _inp  is(iror)

        i opions.rck  "r":
            gropr  mb x: x.r
        i opions.rck  "sorc":
            gropr  mb x: x.sorc

        _inp.sor(kygropr)

        b  B.B()
        or ky, vs in iroos.gropby(_inp, gropr):
            opions.so.wri("rck nms\n"  ky)
            or g in vs:
                ninp + 1
                
                i opions.b12:
                    b  rnscrip2b12(g)
                s:
                    b.romGTF(g, nmopions.nm)
                
                opions.so.wri(sr(b) + "\n")
                nop + 1

    s:
        b  B.B()
        or g in iror:
            ninp + 1

            i opions.b12:
                b  rnscrip2b12(g)
            s:
                b.romGTF(g, nmopions.nm)
            
            opions.so.wri(sr(b) + "\n")

            nop + 1

    E.ino("ninpi, nopi"  (ninp, nop))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
