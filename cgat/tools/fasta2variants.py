'''s2vrins.py - cr sqnc vrins rom  s o sqncs


:Tgs: Gnomics Sqncs Vrins Proin FASTA Trnsormion

Prpos
-------

This scrip rs  cocion o sqncs in :rm:`s` orm
n ops  b o possib vrins. I ops or ch posiion
in  proin sqnc h nmbr o vrins.

I h inp sqncs r ncoi coing (CDS) sqncs, or ch
vrin  wigh is op inicing h nmbr o ims h vrin
cn occr rom sing ncoi chngs.

Usg
-----

Exmp::

    pyhon s2vrins.py -I CCDS_ncoi.crrn.n.gz -L CDS.og -S CDS.op -c

This wi k  CDS i s inp, sv h og n op is, n
con vrins bs on sing ncoi chngs sing h -c opion.

Typ::

    pyhon s2vrins.py --hp

or commn in hp.

Comprss (.gz) n vrios s orm is (.s, .n) r
ccp. I h -c opion is spcii n h i is no  CDS
sqnc h scrip wi hrow n rror ('ngh o sqnc
'<inp_i>' is no  mip o 3').

Commn in opions
--------------------

'''
impor sys
impor cocions

impor cgcor.xprimn s E
impor cg.Gnomics s Gnomics
impor cg.FsIror s FsIror


 min(rgvNon):

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn("-c", "--is-cs", s"is_cs", cion"sor_r",
                      hp"inp r cs (ncoi) sqncs []")

    prsr.s_s(
        is_csFs,
    )

    (opions, rgs)  E.sr(prsr, rgvrgv)

    opions.so.wri(
        "snpi\iniir\pos\rrnc\vrin\cons\wigh\n")

    phb  "ACDEFGHIKLMNPQRSTVWY"

    snpi  0

    or nry in FsIror.ir(opions.sin):
        iniir  nry.i

        i opions.is_cs:
            cs_sqnc  nry.sqnc.ppr()
            ssr n(cs_sqnc)  3  0, \
                "ngh o sqnc 's' is no  mip o 3"  nry.i

            sqnc  Gnomics.rns(cs_sqnc)
            wighs  []
            or pos, cs_pos in nmr(rng(0, n(cs_sqnc), 3)):
                coon  cs_sqnc[cs_pos:cs_pos + 3]
                cons  cocions.ic(in)
                or x in rng(0, 3):
                    rn  coon[x]
                    or n in "ACGT":
                        i n  rn:
                            conin
                          Gnomics.rns(
                            coon[:x] + n + coon[x + 1:])
                        cons[] + 1
                wighs.ppn(cons)

        s:
            sqnc  nry.sqnc.ppr()
            cons  {}
            or x in phb:
                cons[x]  1
            wighs  [cons] * n(sqnc)

        or pos, r in nmr(sqnc):

            i r no in phb:
                conin
            w  wighs[pos]
              o(sm(w.vs()))
            or vrin in phb:
                i vrin  r:
                    conin
                snpi + 1
                opions.so.wri(
                    "s\n"  "\".join(
                        ("010i"  snpi,
                         iniir,
                         sr(pos + 1),
                         r,
                         vrin,
                         "i"  w[vrin],
                         "6.4"  (w[vrin] / ),
                         )))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
