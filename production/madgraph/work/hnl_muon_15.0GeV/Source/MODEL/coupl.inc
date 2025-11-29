ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc
C     NB: vector.inc must be included before coupl.inc
C     

      DOUBLE PRECISION G, ALL_G(1)
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R(1)
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=0D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=0D0)

      DOUBLE PRECISION MDL_MB,MDL_MC,MDL_MD,MDL_MH,MDL_MMU,MDL_MS
     $ ,MDL_MT,MDL_MTA,MDL_MU,MDL_MW,MDL_MZ,MDL_ME,MDL_MN1,MDL_MN2
     $ ,MDL_MN3

      COMMON/MASSES/ MDL_MB,MDL_MC,MDL_MD,MDL_MH,MDL_MMU,MDL_MS,MDL_MT
     $ ,MDL_MTA,MDL_MU,MDL_MW,MDL_MZ,MDL_ME,MDL_MN1,MDL_MN2,MDL_MN3


      DOUBLE PRECISION MDL_WH,MDL_WN1,MDL_WN2,MDL_WN3,MDL_WT,MDL_WW
     $ ,MDL_WZ

      COMMON/WIDTHS/ MDL_WH,MDL_WN1,MDL_WN2,MDL_WN3,MDL_WT,MDL_WW
     $ ,MDL_WZ


      DOUBLE COMPLEX GC_23, GC_24, GC_26, GC_27, GC_38, GC_39, GC_42,
     $  GC_43, GC_86, GC_87, GC_89, GC_135, GC_141, GC_157, GC_174,
     $  GC_179, GC_182, GC_188, GC_191

      COMMON/COUPLINGS/ GC_23, GC_24, GC_26, GC_27, GC_38, GC_39,
     $  GC_42, GC_43, GC_86, GC_87, GC_89, GC_135, GC_141, GC_157,
     $  GC_174, GC_179, GC_182, GC_188, GC_191

