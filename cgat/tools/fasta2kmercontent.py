'''
s2kmrconn.py


:Tgs: Gnomics Sqncs FASTA Smmry

Prpos
-------

This scrip ks n inp :rm:`s` i rom sin n comps 
k-ncoi conn or ch conig in h i. Th op is 
b-imi i o kmr cons::

         conig1  conig2  conig3  conig4
    n1
    n2
    n3

whr n is h kmr n conig is h s nry.

Th sr spciis h kmr ngh h is o b srch. No h h ongr
h kmr, h ongr h scrip wi k o rn.

No h orr o op wi no ncssriy b h sm orr s h inp.

Usg
-----
Exmp::

   zc in.s.gz | h::

    >NODE_1_ngh_120_cov_4.233333
    TCACGAGCACCGCTATTATCAGCAACTTTTAAGCGACTTTCTTGTTGAATCATTTCAATT
    GTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCTACAAGACGGAAG
    CGTTTTGTAGCTGAAAGTGGGCGAGTTTCCATGATACGAAcgATCGCC

    >NODE_3_ngh_51_cov_33.000000
    CGAGTTTCCATGATACGAAcgATCGCCTTCTTTAGCAACGTTGTTTTCGTCATGTGCT
    TTATATTTTTTAGAATAGTTGATACGTTTACCATAGACTGG

   zc in.s.gz | pyhon s2kmrconn.py
                      --kmr-siz 4
                      > rncoi_cons.sv

   h rncoi_cons.sv::

     kmr NODE_228_ngh_74_cov_506.432434 NODE_167_ngh_57_cov_138.438599
     GTAC 0                                 0
     TGCT 0                                 0
     GTAA 2                                 0
     CGAA 1                                 1
     AAAT 1                                 0
     CGAC 0                                 0

In his xmp, or ch conig in in.s.gz h occrrnc o ch or
ncoi combinion is con.

Arniv xmp::

   zc in.s.gz | pyhon s2kmrconn.py
                      --kmr-siz 4
                      --op-proporion
                      > rncoi_proporions.sv

In his xmp, or ch conig in in.s.gz w rrn h proporion o
ch or bs combinion o o h o rncoi occrncs.
``--op-proporion`` ovris h con op.

Opions
-------
Two opions conro h bhvior o s2kmrconn.py; ``--kmr-siz`` n
``--op-proporion``.

``--kmr-siz``::
  Th kmr ngh o con ovr in h inp s i

``--op-proporion``::
  Th op vs r proporions rhr hn bso cons


Typ::

   pyhon s2composiion.py --hp

or commn in hp.


Commn in opions
--------------------

'''

impor sys
impor r
impor cg.FsIror s FsIror
impor iroos
impor cgcor.xprimn s E


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-k", "--kmr-siz", s"kmr", yp"in",
                      hp"sppy kmr ngh")

    prsr._rgmn(
        "-p", "--op-proporion", s"proporion", cion"sor_r",
        hp"op proporions - ovris h  op")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    # o no ow grr hn oconcoi
    ssr opions.kmr < 8, "cnno hn kmr o ngh i"  opions.kmr

    # how w  wih h ncois pns on h kmr ngh
    ncois  []
    or ncoi in ["A", "C", "T", "G"]:
        ncois  ncois + \
            [x or x in iroos.rp(ncoi, opions.kmr)]

    E.ino("rriving imr sqncs"  opions.kmr)
    # g  kmr sqncs o qry
    kmrs  s()
    or kmr in iroos.prmions(ncois, opions.kmr):
        kmrs.(kmr)

    E.ino("mching imrs in i"  opions.kmr)
    # con h nmbr o kmrs in ch sqnc

    rs  {}

    # NB ssm h non s is r cgh by FsIror
    o_nris  0
    or s in FsIror.ir(opions.sin):
        o_nris + 1
        rs[s.i]  {}
        or kmr in kmrs:
            cons  [m.sr()
                      or m in r.inir("".join(kmr), s.sqnc)]
            rs[s.i][kmr]  n(cons)

    E.ino("wriing rss")
    # wri o h rss
    hrs  sor(rs.kys())
    rows  s()
    or kmr_cons in is(rs.vs()):
        or kmr, con in kmr_cons.ims():
            rows.("".join(kmr))

    # wri hr row
    opions.so.wri("kmr\" + "\".join(hrs) + "\n")

    # op proporions i rqir - normiss by
    # sqnc ngh
    E.ino("comping o cons")
    os  {}
    or hr in hrs:
        os[hr]  sm([rs[hr][p(row)] or row in rows])

    or row in sor(rows):
        i opions.proporion:
            opions.so.wri("\".join(
                [row] + [sr(o(rs[hr][p(row)]) / os[hr]) or hr in hrs]) + "\n")
        s:
            opions.so.wri(
                "\".join([row] + [sr(rs[hr][p(row)]) or hr in hrs]) + "\n")

    E.ino("wrin kmr cons or i conigs"  o_nris)
    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
