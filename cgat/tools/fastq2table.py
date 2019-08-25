'''sq2b.py - comp ss on rs in sq is


:Tgs: Gnomics NGS Sqncs FASTQ Annoion

Prpos
-------

This scrip irs ovr  sq i n ops
smmry sisics or ch r.

Th op is  b-imi x i wih h oowing comns:

+----------+----------------------------------------+
|*Comn*  |*Conn*                               |
+----------+----------------------------------------+
|r      |r iniir prsn in inp sq  |
|          |i                                    |
+----------+----------------------------------------+
|ni   |nmbr o rs h  bow Q10     |
+----------+----------------------------------------+
|nN        |nmbr o mbigos bs cs (N)      |
+----------+----------------------------------------+
|nv      |nmbr o bss in h r             |
+----------+----------------------------------------+
|min       |minimm bs qiy scor or h r |
+----------+----------------------------------------+
|mx       |mximm bs qiy or h r       |
+----------+----------------------------------------+
|mn      |mn bs qiy or h r          |
+----------+----------------------------------------+
|min    |min bs qiy or h r        |
+----------+----------------------------------------+
|sv    |snr viion o qiy scors   |
|          |or h r                            |
+----------+----------------------------------------+
|sm       |sm o qiy scors or h r      |
+----------+----------------------------------------+
|q1        |25h prcni o qiy scors or   |
|          |h r                                |
+----------+----------------------------------------+
|q3        |25h prcni o qiy scors or   |
|          |h r                                |
+----------+----------------------------------------+

Usg
-----

Exmp::

   cg sq2b --gss-ormsngr < in.sq > o.sv

In his xmp w know h or  hv qiy scors orm s
sngr. Givn h imin-1.8 qiy scors r highy ovrpping
wih sngr, his opion s o sngr qiis. In  mo
h scrip my no b b o isingish highy ovrpping ss o
qiy scors.

I w provi wo rs o h scrip::

   @DHKW5DQ1:308:D28FGACXX:5:2211:8051:4398
   ACAATGTCCTGATGTGAATGCCCCTACTATTCAGATCGCTTAGGGCATGC
   +
   B1?DFDDHHFFHIJJIJGGIJGFIEE9CHIIFEGGIIJGIGIGIIDGHI
   @DHKW5DQ1:308:D28FGACXX:5:1315:15039:83265
   GAATGCCCCTACTATTCAGATCGCTTAGGGCATGCGTCGCATGTGAGTAA
   +
   @@@FDFFFHGHHHJIIIJIGHIJJIGHGHC9FBFBGHIIEGHIGC>F@FA

w g h oowing b s op:

+-----------------------------------------+-------+--+----+-------+-------+-------+-------+------+---------+-------+-------+
|r                                     |ni|nN|nv|min    |mx    |mn   |min |sv|sm      |q1     |q3     |
+-----------------------------------------+-------+--+----+-------+-------+-------+-------+------+---------+-------+-------+
|DHKW5DQ1:308:D28FGACXX:5:2211:8051:4398  |0      |0 |50  |16.0000|41.0000|37.2000|38.0000|4.4900|1860.0000|36.0000|40.0000|
+-----------------------------------------+-------+--+----+-------+-------+-------+-------+------+---------+-------+-------+
|DHKW5DQ1:308:D28FGACXX:5:1315:15039:83265|0      |0 |50  |24.0000|41.0000|37.0200|38.0000|3.5916|1851.0000|36.0000|40.0000|
+-----------------------------------------+-------+--+----+-------+-------+-------+-------+------+---------+-------+-------+


Typ::

   cg sq2b --hp

or commn in hp.

Commn in opions
--------------------

'''

impor sys

impor cgcor.xprimn s E
impor cg.Ss s Ss
impor cg.Fsq s Fsq


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(vrsion"prog vrsion: $I$",
                            sggobs()["__oc__"])

    prsr._rgmn(
        "--gss-orm", s"gss_orm", yp"choic",
        choics(
            'sngr', 'sox', 'phr64', 'imin-1.8', 'ingr'),
        hp"Th  bhvior o h scrip is o gss h qiy "
        "orm o h inp sq i. Th sr cn spciy h "
        "qiy orm o h inp i sing h --gss-orm opion. "
        "Th scrip wi s his orm i h "
        "sqnc qiis r mbigos.[].")

    prsr._rgmn(
        "--rg-orm", s"rg_orm", yp"choic",
        choics(
            'sngr', 'sox', 'phr64', 'imin-1.8', 'ingr'),
        hp"Th scrip wi convr qiy scors o h sinion "
        "orm nss [].")

    prsr.s_s(
        rg_ormNon,
        gss_ormNon,
        min_qiy10,
    )

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    c  E.Conr()

    i opions.rg_orm:
        iror  Fsq.ir_convr(opions.sin,
                                         ormopions.rg_orm,
                                         gssopions.gss_orm)
    s:
        iror  Fsq.ir_gss(opions.sin,
                                       gssopions.gss_orm)

    opions.so.wri("r\ni\nN\s\n" 
                         ("\".join(Ss.Smmry().gHrs())))

    min_qiy  opions.min_qiy

    or rcor in iror:
        c.inp + 1
        qs  rcor.oPhr()
        ni  n([x or x in qs i x < min_qiy])
        nns  rcor.sq.con("N") + rcor.sq.con(".")
        opions.so.wri("s\i\i\s\n"  (rcor.iniir,
                                                   ni,
                                                   nns,
                                                   sr(Ss.Smmry(qs))
                                                   ))
        c.op + 1

    # wri oor n op bnchmrk inormion.
    E.ino("s"  sr(c))
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
