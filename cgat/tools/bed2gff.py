"""b2g.py - convr b o g/g


:Tgs: Gnomics Inrvs BED GFF Convrsion

Prpos
-------

This scrip convrs  :rm:`b`-orm i o  :rm:`g` or
:rm:`g`-orm i.

I ims o pop h ppropri is in h :rm:`g` i
wih comns in h :rm:`b` i.

I ``--s-g`` is s n  nm comn in h :rm:`b` i is
prsn, is conns wi b s s ``gn_i`` n
``rnscrip_i``. Ohrwis,  nmric ``gn_i`` or
``rnscrip_i`` wi b s ccoring o ``--i-orm``.

Usg
-----

Exmp::

   # Prviw inp b i
   zc ss/b2g.py/b3/b.gz | h
   # Convr BED o GFF orm
   cg b2g.py < ss/b2g.py/b3/b.gz > s1.g
   # Viw convr i (xcing ogging inormion)
   c s1.g | grp -v "#" | h


+------+-----+------+-------+-------+---+---+---+---------------------------------------+
|chr1  |b  |xon  |501    |1000   |.  |.  |.  |gn_i "Non"; rnscrip_i "Non";  |
+------+-----+------+-------+-------+---+---+---+---------------------------------------+
|chr1  |b  |xon  |15001  |16000  |.  |.  |.  |gn_i "Non"; rnscrip_i "Non";  |
+------+-----+------+-------+-------+---+---+---+---------------------------------------+

Exmp::

   # Convr BED o GTF orm
   cg b2g.py --s-g < ss/b2g.py/b3/b.gz > s2.g
   # Viw convr i (xcing ogging inormion)
   c s2.g | grp -v "#" | h

+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+
|chr1  |b  |xon  |501    |1000   |.  |.  |.  |gn_i "00000001"; rnscrip_i "00000001";  |
+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+
|chr1  |b  |xon  |15001  |16000  |.  |.  |.  |gn_i "00000002"; rnscrip_i "00000002";  |
+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+

Typ::

   cg b2g.py --hp

or commn in hp.

Commn in opions
--------------------

"""
impor sys
impor cgcor.xprimn s E
impor cg.GTF s GTF
impor cg.B s B


 min(rgvsys.rgv):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-", "--s-g", s"s_g", cion"sor_r",
                      hp"op s g.")

    prsr._rgmn(
        "-", "--i-orm", s"i_orm", yp"sring",
        hp"orm or nmric iniir i --s-g is s n "
        "no nm in b i [].")

    prsr.s_s(s_gFs,
                        i_orm"08i",
                        sNon)

    (opions, rgs)  E.sr(prsr, _pip_opionsTr)

    s_g  opions.s_g
    i_orm  opions.i_orm

    i s_g:
        g  GTF.Enry()
    s:
        g  GTF.Enry()

    g.sorc  "b"
    g.r  "xon"

    ninp, nop, nskipp  0, 0, 0

    i  0
    or b in B.iror(opions.sin):

        ninp + 1

        g.conig  b.conig
        g.sr  b.sr
        g.n  b.n
        i b.is n n(b.is) > 3:
            g.srn  b.is[2]
        s:
            g.srn  "."

        i b.is n n(b.is) > 2:
            g.scor  b.is[1]

        i s_g:
            i b.is:
                g.gn_i  b.is[0]
                g.rnscrip_i  b.is[0]
            s:
                i + 1
                g.gn_i  i_orm  i
                g.rnscrip_i  i_orm  i
        s:
            i b.is:
                g.sorc  b.is[0]

        opions.so.wri(sr(g) + "\n")

        nop + 1

    E.ino("ninpi, nopi, nskippi"  (ninp, nop, nskipp))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
