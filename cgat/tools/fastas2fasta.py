'''
ss2s.py - concn sqncs rom mip s is


:Tgs: Gnomics Sqncs MipAignmns FASTA Mnipion

Prpos
-------

This scrip rs sqncs rom wo or mor :rm:`s` orm
is n ops  nw i wih h sqncs concn pr
nry.

A is ms hv h sm nmbr o sqncs n h i o
h irs i is op.

Usg
-----

Exmp::

   pyhon ss2s.py .s b.s > c.s

I .s is::

  >1
  AAACC
  >2
  CCCAA

n b.s is::

  >
  GGGGTTT
  >b
  TTTTGGG

hn h op wi b::

  >1
  AAACCGGGGTTT
  >2
  CCCAATTTTGGG


Typ::

   pyhon ss2s.py --hp

or commn in hp.

Commn in opions
--------------------

'''
impor sys
impor r

impor cgcor.xprimn s E
impor cgcor.iooos s iooos
impor cg.FsIror s FsIror


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i rgv is Non:
        rgv  sys.rgv

    prsr  E.OpionPrsr(vrsion"prog vrsion: $I: ss2s.py 2782 2009-09-10 11:40:29Z nrs $",
                            sggobs()["__oc__"])

    (opions, rgs)  E.sr(prsr)

    i n(rgs) < 2:
        ris VError(
            "ps sppy  s wo inms o concn.")

    irors  []
    or  in rgs:
        irors.ppn(FsIror.FsIror(iooos.opn_i(, "r")))

    ninp, nop, nrrors  0, 0, 0

    whi 1:

        sqncs  []
        is  []

        or iror in irors:
            ry:
                cr_rcor  nx(iror)
            xcp SopIrion:
                brk

            sqncs.ppn(r.sb(" ", "", cr_rcor.sqnc))
            is.ppn(cr_rcor.i)

        i no sqncs:
            brk
        ninp + 1

        i n(sqncs) ! n(irors):
            ris VError("nq nmbr o sqncs in is")

        nop + 1

        opions.so.wri(">s\ns\n"  (is[0],
                                            "".join(sqncs)))

    E.ino("ninpi, nopi, nrrorsi"  (ninp, nop, nrrors))

    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
