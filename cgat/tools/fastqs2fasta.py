'''
sqs2s.py - inrv wo sq is


:Tgs: Gnomics NGS FASTQ FASTA Convrsion

Prpos
-------

This scrip is s o inrv wo :rm:`sq`-orm is
(pir ) ino  sing :rm:`s`-orm i. R1 is
oow by r2 in h rsn i.

:rm:`sq` is MUST b sor by r iniir.

Usg
-----

For xmp::

   cg sqs2s \
         --irs-sq-iin.sq.1.gz \
         --scon-sq-iin.sq.2.gz > o.s

I :i:`in.sq.1.gz` ooks ik his::

    @r1_rom_gi|387760314|r|NC_017594.1|_Srpococcs_siv_#0/1
    TTCTTGTTGAATCATTTCAATTGTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCT
    +
    ??A???ABBDDBDDEDGGFGAFHHCHHIIIDIHGIFIHHFICIHDHIHIFIFIIIIIIHFHIFHIHHHH
    @r3_rom_gi|315441696|r|NC_014814.1|_Mycobcrim_givm_#0/1
    ATGAACGCGGCCGAGCAACACCGCCACCACGTGAATCGGTGGTTCTACGACTGCCCGTCGGCCTTCCACC
    +

n :i:`in.sq.2.gz` ooks ik his::

    A??A?B??BDBDDDBDGGFA>CFCFIIIIIIF;HFIGHCIGHIHHEHHHIIHHFDHH-HD-IDHHHGIHG
    @r1_rom_gi|387760314|r|NC_017594.1|_Srpococcs_siv_#0/2
    ACCTTCGTTTCCAAGGTGCAGCAGGTCAACTTGATCAAACTGCCCCTTTGAACGAAGTGAAAAAACAAAT
    +
    A????@BBDBDDADABGFGFFEHHHIEHHII@IIHIHHIDHCCIHIIIHHIEI5HIHFHIEHIHCHHC)
    @r3_rom_gi|315441696|r|NC_014814.1|_Mycobcrim_givm_#0/2
    GGGAGCCTGCAGCGCCGCCGCGACTGCATCGCCGCGGCCGGCATCGTGGGATGGACGGTGCGTCAGACGC
    +
    ???A?9BBDDD5@DDDGFFGFFHIIIHHIHBFHIIHIIHHH>HEIHHFI>FFHGIIHHHDHCCFIHFIHD

hn h op wi b::

  >r1_rom_gi|387760314|r|NC_017594.1|_Srpococcs_siv_#0/1
  TTCTTGTTGAATCATTTCAATTGTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCT
  >r1_rom_gi|387760314|r|NC_017594.1|_Srpococcs_siv_#0/2
  ACCTTCGTTTCCAAGGTGCAGCAGGTCAACTTGATCAAACTGCCCCTTTGAACGAAGTGAAAAAACAAAT
  >r3_rom_gi|315441696|r|NC_014814.1|_Mycobcrim_givm_#0/1
  ATGAACGCGGCCGAGCAACACCGCCACCACGTGAATCGGTGGTTCTACGACTGCCCGTCGGCCTTCCACC
  >r3_rom_gi|315441696|r|NC_014814.1|_Mycobcrim_givm_#0/2
  GGGAGCCTGCAGCGCCGCCGCGACTGCATCGCCGCGGCCGGCATCGTGGGATGGACGGTGCGTCAGACGC
  >r4_rom_gi|53711291|r|NC_006347.1|_Bcrois_rgiis_#0/1
  GAGGGATCAGCCTGTTATCCCCGGAGTACCTTTTATCCTTTGAGcgGTCCCTTCCATACGGAAACACC
  >r4_rom_gi|53711291|r|NC_006347.1|_Bcrois_rgiis_#0/2
  CAACCGTGAGCTCAGTGAAATTGTAGTATCGGTGAAGATGCcgTACCCGcgGGGACGAAAAGACCC
  >r5_rom_gi|325297172|r|NC_015164.1|_Bcrois_snir_#0/1
  TGCGGCGAAATACCAGCCCATGCCCCGTCCCCAGAATTCCTTGGAGCAGCCTTTGTGAGGTTCGGCTTTG
  >r5_rom_gi|325297172|r|NC_015164.1|_Bcrois_snir_#0/2
  AACGGCACGCACAATGCCGACCGCTACAAAAAGGCTGCCGACTGGCTCCGCAATTACCTGGTGAACGACT


Typ::

   cg sqs2s --hp

or commn in hp.


Commn in opions
--------------------

'''

impor sys
rom iroos impor zip_ongs

impor cgcor.iooos s iooos
impor cg.Fsq s Fsq
impor cgcor.xprimn s E


css PirRError(Excpion):

    '''
    xcpion ris whn rs rn' pir -
    co b no sor or is o irn nghs
    '''


 min(rgvNon):
    """scrip min.

    prss commn in opions in sys.rgv, nss *rgv* is givn.
    """

    i no rgv:
        rgv  sys.rgv

    # sp commn in prsr
    prsr  E.OpionPrsr(
        vrsion"prog vrsion: $I$",
        sggobs()["__oc__"])

    prsr._rgmn(
        "-", "--irs-sq-i", s"sq1", yp"sring",
        hp"sppy r1 sq i")
    prsr._rgmn(
        "-b", "--scon-sq-i", s"sq2", yp"sring",
        hp"sppy r2 sq i")

    #  common opions (-h/--hp, ...) n prs commn in
    (opions, rgs)  E.sr(prsr, rgvrgv)

    i rgs n n(rgs)  2:
        opions.sq1, opions.sq2  rgs

    sq1  iooos.opn_i(opions.sq1)
    sq2  iooos.opn_i(opions.sq2)

    E.ino("iring ovr sq is")
    1_con  0
    or 1, 2 in zip_ongs(Fsq.ir(sq1),
                              Fsq.ir(sq2)):
        i no (1 n 2) or (no 2 n 1):
            ry:
                ris PirRError(
                    "npir rs c. Ar is sor? r "
                    "is o q ngh?")
            xcp PirRError s :
                ris PirRError().wih_rcbck(sys.xc_ino()[2])
        s:
            ssr 1.iniir.nswih("/1") n \
                2.iniir.nswih("/2"), \
                "Rs in i 1 ms n wih /1 n rs in i 2 wih /2"
            opions.so.wri(
                ">s\ns\n>s\ns\n" 
                (1.iniir, 1.sq, 2.iniir, 2.sq))
            1_con + 1

    E.ino("op: i pirs"  1_con)

    # wri oor n op bnchmrk inormion.
    E.sop()

i __nm__  "__min__":
    sys.xi(min(sys.rgv))
