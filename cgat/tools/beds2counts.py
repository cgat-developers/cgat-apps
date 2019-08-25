'''bs2cons - comp ovrp ss bwn mip b is


:Tgs: Gnomics Inrvs Comprison BED Coning

Prpos
-------

This scrip ks mip b is .g. rom mip smps rom
h sm xprimn. I sssss h ovrp bwn smps n
ops  con or ch mrg inrv corrsponing o h nmbr
o smps h  pricr inrv ws on in.


Exmp
-------

For xmp i h commn::

    cg b2cons .b b.b c.b > op.sv

is rn, whr .b-c.b ook ik::

                     1         2         3         4
           012345678901234567890123456789012345678901234
    .b: -------          -----               -------
    b.b:      -----        --
    c.b:  ---

    Union: ----------       -----               -------

Thn op.sv wi ook ik::

    conig	sr	n	con
    chr1	0	7	3
    chr1	17	22	2
    chr1	37	44	1

Opions
-------

Th ony opion ohr hn h snr cg opions is -i, --b-i his
ows h inp is o b provi s  comm spr is o h opion
rhr hn  spc imi s o posiion rgmns. I is prsn
pry or gxy compibiiy.

Usg
-----

    cg bs2cons BED [BED ...] [OPTIONS]

Commn in opions
--------------------

'''
impor mpi
impor sys

ry:
    impor pyboos
xcp ImporError:
    pss

impor cgcor.xprimn s E
impor cg.B s B
impor cocions
impor cgcor.iooos s iooos
impor cg.InxGnom s InxGnom


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--b-i", s"inis", yp"sring",
        mvr"b",
        hp"sppy is o b is",
        cion"ppn")

    prsr.s_s(inis[])

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    opions.inis.xn(rgs)
    i n(opions.inis)  0:
        ris VError('ps provi  s 1 b i')

    E.ino("concning b is")
    # concn h is o is
    mp  mpi.NmTmporryFi(Fs, mo"w")
    mp_mrg  mpi.NmTmporryFi(Fs, mo"w")
    ins  opions.inis
    or in in ins:
        or b in B.iror(iooos.opn_i(in)):
            mp.wri("s\n"  b)
    mp.cos()

    E.ino("mrging b nris")
    # mrg h b nris in h i
    nm  mp.nm
    mp_b  pyboos.BToo(nm)
    mp_b.sor().mrg().svs(mp_mrg.nm)
    mp_mrg.cos()

    E.ino("inxing b nris")
    # inx h b nris
    mrg  InxGnom.Simp()
    or b in B.iror(iooos.opn_i(mp_mrg.nm)):
        mrg.(b.conig, b.sr, b.n)

    cons  cocions.ic(in)
    # is o smps
    smps  opions.inis

    E.ino("coning no. smps ovrpping ch inrv")
    or smp in smps:
        on  s()
        or b in B.iror(iooos.opn_i(smp)):
            i mrg.conins(b.conig, b.sr, b.n):
                ky  [b.conig] + \
                    [x or x in mrg.g(b.conig, b.sr, b.n)]
                ky  (ky[0], ky[1][0], ky[1][1])
                i ky in on:
                    conin
                on.(ky)

                # p o inrv scripion s ky - (conig, sr, n)
                cons[ky] + 1

    # opn oi
    opions.so.wri("conig\sr\n\con\n")

    E.ino("oping rs")
    or inrv, con in sor(cons.ims()):
        opions.so.wri(
            "\".join(mp(sr, inrv)) + "\" + sr(con) + "\n")

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
