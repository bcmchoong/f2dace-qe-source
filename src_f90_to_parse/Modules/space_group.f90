!# 1 "space_group.f90"
!
! Copyright (C) 2014 Federico Zadra
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE space_group
  !! Space groups, symmetries, Bravais lattices.
  USE kinds, ONLY: DP
  IMPLICIT NONE
!# 13 "space_group.f90"
  REAL(DP), PARAMETER :: unterz=(1.0_DP)/(3.0_DP)
  REAL(DP), PARAMETER :: duterz=(2.0_DP)/(3.0_DP)
  REAL(DP), PARAMETER :: unsest=(1.0_DP)/(6.0_DP)
  REAL(DP), PARAMETER :: cisest=(5.0_DP)/(6.0_DP)
!# 18 "space_group.f90"
   SAVE
   PRIVATE
   PUBLIC sym_brav, find_equivalent_tau
   
 CONTAINS
   
   SUBROUTINE sym_brav(space_group_number,sym_n,ibrav)
   
   !! Sym\_brav:  
   !! - input: spacegroup number;  
   !! - output: symmetries number (\(\text{sym_n}\)), Bravais
   !! lattice number (\(\text{ibrav}\))
!# 31 "space_group.f90"
      INTEGER, INTENT(IN) :: space_group_number
      INTEGER, INTENT(OUT) :: sym_n,ibrav
      
      simmetria: SELECT CASE (space_group_number)
      !Triclinic 1-2
      CASE (1)
         sym_n=1
         ibrav=14
      CASE (2)
         sym_n=2
         ibrav=14
      !Monoclinic 3-15
      CASE (3) !P2
         sym_n=2
         ibrav=12
      CASE (4) !P2(1)
         sym_n=2
         ibrav=12
      CASE (5) !C2
         sym_n=2
         ibrav=13
      CASE (6) !PM
         sym_n=2
         ibrav=12
      CASE (7) !Pc
         sym_n=2
         ibrav=12
      CASE (8) !Cm
         sym_n=2
         ibrav=13
      CASE (9) !Cc
         sym_n=2
         ibrav=13
      CASE (10) !P2/m
         sym_n=4
         ibrav=12
      CASE (11) !P2(1)/m
         sym_n=4
         ibrav=12
      CASE (12) !C2/m
         sym_n=4
         ibrav=13
      CASE (13) !P2/c
         sym_n=4
         ibrav=12
      CASE (14) !P2(1)/c
         sym_n=4
         ibrav=12
      CASE (15) !C2/c
         sym_n=4
         ibrav=13
      !Orthorhombic
      CASE (16) !P222
         sym_n=4
         ibrav=8
      CASE (17) !P222(1)
         sym_n=4
         ibrav=8
      CASE (18) !P2(1)2(1)2
         sym_n=4
         ibrav=8
      CASE (19) !P2(1)2(1)2(1)
         sym_n=4
         ibrav=8
      CASE (20) !C222(1)
         sym_n=4
         ibrav=9
      CASE (21) !C222
         sym_n=4
         ibrav=9
      CASE (22) !F222
         sym_n=4
         ibrav=10
      CASE (23) !I222
         sym_n=4
         ibrav=11
      CASE (24) !I2(1)2(1)2(1)
         sym_n=4
         ibrav=11
      CASE (25) !Pmm2
         sym_n=4
         ibrav=8
      CASE (26) !Pmc2(1)
         sym_n=4
         ibrav=8
      CASE (27) !Pcc2
         sym_n=4
         ibrav=8
      CASE (28) !Pma2
         sym_n=4
         ibrav=8
      CASE (29) !Pca2(1)
         sym_n=4
         ibrav=8
      CASE (30) !Pnc2
         sym_n=4
         ibrav=8
      CASE (31) !Pmn2(1)
         sym_n=4
         ibrav=8
      CASE (32) !Pba2
         sym_n=4
         ibrav=8
      CASE (33) !Pna2(1)
         sym_n=4
         ibrav=8
      CASE (34) !Pnn2
         sym_n=4
         ibrav=8
      CASE (35) !Cmm2
         sym_n=4
         ibrav=9
      CASE (36) !Cmc2(1)
         sym_n=4
         ibrav=9
      CASE (37) !Ccc2
         sym_n=4
         ibrav=9
      CASE (38) !Amm2
         sym_n=4
         ibrav=91
      CASE (39) !Abm2
         sym_n=4
         ibrav=91
      CASE (40) !Ama2
         sym_n=4
         ibrav=91
      CASE (41) !Aba2
         sym_n=4
         ibrav=91
      CASE (42) !Fmm2
         sym_n=4
         ibrav=10
      CASE (43) !Fdd2
         sym_n=4
         ibrav=10
      CASE (44) !Imm2
         sym_n=4
         ibrav=11
      CASE (45) !Iba2
         sym_n=4
         ibrav=11
      CASE (46) !Ima2
         sym_n=4
         ibrav=11
      CASE (47) !Pmmm
         sym_n=8
         ibrav=8
      CASE (48) !Pnnn
         sym_n=8
         ibrav=8
      CASE (49) !Pccm
         sym_n=8
         ibrav=8
      CASE (50) !Pban
         sym_n=8
         ibrav=8
      CASE (51) !Pmma
         sym_n=8
         ibrav=8
      CASE (52) !Pnna
         sym_n=8
         ibrav=8
      CASE (53) !Pmna
         sym_n=8
         ibrav=8
      CASE (54) !Pcca
         sym_n=8
         ibrav=8
      CASE (55) !Pbam
         sym_n=8
         ibrav=8
      CASE (56) !Pccn
         sym_n=8
         ibrav=8
      CASE (57) !Pbcm
         sym_n=8
         ibrav=8
      CASE (58) !Pnnm
         sym_n=8
         ibrav=8
      CASE (59) !Pmmn
         sym_n=8
         ibrav=8
      CASE (60) !Pbcn
         sym_n=8
         ibrav=8
      CASE (61) !Pbca
         sym_n=8
         ibrav=8
      CASE (62) !Pnma
         sym_n=8
         ibrav=8
      CASE (63) !Cmcm
         sym_n=8
         ibrav=9
      CASE (64) !Cmca
         sym_n=8
         ibrav=9
      CASE (65) !Cmmm
         sym_n=8
         ibrav=9
      CASE (66) !Cccm
         sym_n=8
         ibrav=9
      CASE (67) !Cmma
         sym_n=8
         ibrav=9
      CASE (68) !Ccca
         sym_n=8
         ibrav=9
      CASE (69) !Fmmm
         sym_n=8
         ibrav=10
      CASE (70) !Fddd
         sym_n=8
         ibrav=10
      CASE (71) !Immm
         sym_n=8
         ibrav=11
      CASE (72) !Ibam
         sym_n=8
         ibrav=11
      CASE (73) !Ibca
         sym_n=8
         ibrav=11
      CASE (74) !Imma
         sym_n=8
         ibrav=11
      !Tetragonal
      CASE (75) !P4
         sym_n=4
         ibrav=6
      CASE (76) !P4(1)
         sym_n=4
         ibrav=6
      CASE (77) !P4(2)
         sym_n=4
         ibrav=6
      CASE (78) !P4(3)
         sym_n=4
         ibrav=6
      CASE (79) !I4
         sym_n=4
         ibrav=7
      CASE (80) !I4(1)
         sym_n=4
         ibrav=7
      CASE (81) !P-4
         sym_n=4
         ibrav=6
      CASE (82) !I-4
         sym_n=4
         ibrav=7
      CASE (83) !P4/m
         sym_n=8
         ibrav=6
      CASE (84) !P4(2)/m
         sym_n=8
         ibrav=6
      CASE (85) !P4/n
         sym_n=8
         ibrav=6
      CASE (86) !P4(2)/n
         sym_n=8
         ibrav=6
      CASE (87) !I4/m
         sym_n=8
         ibrav=7
      CASE (88) !I4(1)/a
         sym_n=8
         ibrav=7
      CASE (89) !P422
         sym_n=8
         ibrav=6
      CASE (90) !P42(1)2
         sym_n=8
         ibrav=6
      CASE (91) !P4(1)22
         sym_n=8
         ibrav=6
      CASE (92) !P4(1)2(1)2
         sym_n=8
         ibrav=6
      CASE (93) !P4(2)22
         sym_n=8
         ibrav=6
      CASE (94) !P4(2)2(1)2
         sym_n=8
         ibrav=6
      CASE (95) !P4(3)22
         sym_n=8
         ibrav=6
      CASE (96) !P4(3)2(1)2
         sym_n=8
         ibrav=6
      CASE (97) !I422
         sym_n=8
         ibrav=7
      CASE (98) !I4(1)22
         sym_n=8
         ibrav=7
      CASE (99) !P4mm
         sym_n=8
         ibrav=6
      CASE (100) !P4bm
         sym_n=8
         ibrav=6
      CASE (101) !P4(2)cm
         sym_n=8
         ibrav=6
      CASE (102) !P4(2)nm
         sym_n=8
         ibrav=6
      CASE (103) !P4cc
         sym_n=8
         ibrav=6
      CASE (104) !P4nc
         sym_n=8
         ibrav=6
      CASE (105) !P4(2)mc
         sym_n=8
         ibrav=6
      CASE (106) !P4(2)bc
         sym_n=8
         ibrav=6
      CASE (107) !I4mm
         sym_n=8
         ibrav=7
      CASE (108) !I4cm
         sym_n=8
         ibrav=7
      CASE (109) !I4(!)md
         sym_n=8
         ibrav=7
      CASE (110) !I4(1)cd
         sym_n=8
         ibrav=7
      CASE (111) !P-42m
         sym_n=8
         ibrav=6
      CASE (112) !P-42c
         sym_n=8
         ibrav=6
      CASE (113) !P-42(1)m
         sym_n=8
         ibrav=6
      CASE (114) !P-42(1)c
         sym_n=8
         ibrav=6
      CASE (115) !P-4m2
         sym_n=8
         ibrav=6
      CASE (116) !P-4c2
         sym_n=8
         ibrav=6
      CASE (117) !P-4b2
         sym_n=8
         ibrav=6
      CASE (118) !P-4n2
         sym_n=8
         ibrav=6
      CASE (119) !I-4m2
         sym_n=8
         ibrav=7
      CASE (120) !I-4c2
         sym_n=8
         ibrav=7
      CASE (121) !I-42m
         sym_n=8
         ibrav=7
      CASE (122) !I-42d
         sym_n=8
         ibrav=7
      CASE (123) !P4/mmm
         sym_n=16
         ibrav=6
      CASE (124) !P4/mcc
         sym_n=16
         ibrav=6
      CASE (125) !P4/nbm
         sym_n=16
         ibrav=6
      CASE (126) !P4/nnc
         sym_n=16
         ibrav=6
      CASE (127) !P4/mbm
         sym_n=16
         ibrav=6
      CASE (128) !P4/mnc
         sym_n=16
         ibrav=6
      CASE (129) !P4/nmm
         sym_n=16
         ibrav=6
      CASE (130) !P4/ncc
         sym_n=16
         ibrav=6
      CASE (131) !P4(2)/mmc
         sym_n=16
         ibrav=6
      CASE (132) !P4(2)/mcm
         sym_n=16
         ibrav=6
      CASE (133) !P4(2)nbc
         sym_n=16
         ibrav=6
      CASE (134) !P4(2)/nnm
         sym_n=16
         ibrav=6
      CASE (135) !P4(2)/mbc
         sym_n=16
         ibrav=6
      CASE (136) !P4(2)/mnm
         sym_n=16
         ibrav=6
      CASE (137) !P4(2)/nmc
         sym_n=16
         ibrav=6
      CASE (138) !P4(2)/ncm
         sym_n=16
         ibrav=6
      CASE (139) !I4/mmm
         sym_n=16
         ibrav=7
      CASE (140) !I4/mcm
         sym_n=16
         ibrav=7
      CASE (141) !I4(1)/amd
         sym_n=16
         ibrav=7
      CASE (142) !I4(1)/acd
         sym_n=16
         ibrav=7
      ! Trigonal
      CASE (143) !P3
         sym_n=3
         ibrav=4
      CASE (144)
         sym_n=3
         ibrav=4
      CASE (145)
         sym_n=3
         ibrav=4
      CASE (146) !R3
         sym_n=3
         ibrav=5
      CASE (147)
         sym_n=6
         ibrav=4
      CASE (148) !R-3
         sym_n=6
         ibrav=5
      CASE (149)
         sym_n=6
         ibrav=4
      CASE (150) 
         sym_n=6
         ibrav=4
      CASE (151)
         sym_n=6
         ibrav=4
      CASE (152)
         sym_n=6
         ibrav=4
      CASE (153)
         sym_n=6
         ibrav=4
      CASE (154)
         sym_n=6
         ibrav=4
      CASE (155) !R32
         sym_n=6
         ibrav=5
      CASE (156)
         sym_n=6
         ibrav=4
      CASE (157) 
         sym_n=6
         ibrav=4
      CASE (158)
         sym_n=6
         ibrav=4
      CASE (159)
         sym_n=6
         ibrav=4
      CASE (160) !R3m
         sym_n=6
         ibrav=5
      CASE (161) !R3c
         sym_n=6
         ibrav=5
      CASE (162)
         sym_n=12
         ibrav=4
      CASE (163)
         sym_n=12
         ibrav=4
      CASE (164)
         sym_n=12
         ibrav=4
      CASE (165)
         sym_n=12
         ibrav=4
      CASE (166) !R-3m
         sym_n=12
         ibrav=5
      CASE (167) !R-3c
         sym_n=12
         ibrav=5
      ! Exagonal
      CASE (168)
         sym_n=6
         ibrav=4
      CASE (169)
         sym_n=6
         ibrav=4
      CASE (170)
         sym_n=6
         ibrav=4
      CASE (171)
         sym_n=6
         ibrav=4
      CASE (172)
         sym_n=6
         ibrav=4
      CASE (173)
         sym_n=6
         ibrav=4
      CASE (174)
         sym_n=6
         ibrav=4
      CASE (175)
         sym_n=12
         ibrav=4
      CASE (176)
         sym_n=12
         ibrav=4
      CASE (177)
         sym_n=12
         ibrav=4
      CASE (178)
         sym_n=12
         ibrav=4
      CASE (179)
         sym_n=12
         ibrav=4
      CASE (180)
         sym_n=12
         ibrav=4
      CASE (181)
         sym_n=12
         ibrav=4
      CASE (182)
         sym_n=12
         ibrav=4
      CASE (183)
         sym_n=12
         ibrav=4
      CASE (184)
         sym_n=12
         ibrav=4
      CASE (185)
         sym_n=12
         ibrav=4
      CASE (186)
         sym_n=12
         ibrav=4
      CASE (187)
         sym_n=12
         ibrav=4
      CASE (188)
         sym_n=12
         ibrav=4
      CASE (189)
         sym_n=12
         ibrav=4
      CASE (190)
         sym_n=12
         ibrav=4
      CASE (191)
         sym_n=24
         ibrav=4
      CASE (192)
         sym_n=24
         ibrav=4
      CASE (193)
         sym_n=24
         ibrav=4
      CASE (194)
         sym_n=24
         ibrav=4
      !Cubic
      CASE (195)
         sym_n=12
         ibrav=1
      CASE (196)
         sym_n=12
         ibrav=2
      CASE (197)
         sym_n=12
         ibrav=3
      CASE (198)
         sym_n=12
         ibrav=1
      CASE (199)
         sym_n=12
         ibrav=3
      CASE (200)
         sym_n=24
         ibrav=1
      CASE (201)
         sym_n=24
         ibrav=1
      CASE (202)
         sym_n=24
         ibrav=2
      CASE (203)
         sym_n=24
         ibrav=2
      CASE (204)
         sym_n=24
         ibrav=3
      CASE (205)
         sym_n=24
         ibrav=1
      CASE (206)
         sym_n=24
         ibrav=3
      CASE (207)
         sym_n=24
         ibrav=1
      CASE (208)
         sym_n=24
         ibrav=1
      CASE (209)
         sym_n=24
         ibrav=2
      CASE (210)
         sym_n=24
         ibrav=2
      CASE (211)
         sym_n=24
         ibrav=3
      CASE (212)
         sym_n=24
         ibrav=1
      CASE (213)
         sym_n=24
         ibrav=1
      CASE (214)
         sym_n=24
         ibrav=3
      CASE (215)
         sym_n=24
         ibrav=1
      CASE (216)
         sym_n=24
         ibrav=2
      CASE (217)
         sym_n=24
         ibrav=3
      CASE (218)
         sym_n=24
         ibrav=1
      CASE (219)
         sym_n=24
         ibrav=2
      CASE (220)
         sym_n=24
         ibrav=3
      CASE (221)
         sym_n=48
         ibrav=1
      CASE (222)
         sym_n=48
         ibrav=1
      CASE (223)
         sym_n=48
         ibrav=1
      CASE (224)
         sym_n=48
         ibrav=1
      CASE (225)
         sym_n=48
         ibrav=2
      CASE (226)
         sym_n=48
         ibrav=2
      CASE (227)
         sym_n=48
         ibrav=2
      CASE (228)
         sym_n=48
         ibrav=2
      CASE (229)
         sym_n=48
         ibrav=3
      CASE (230)
         sym_n=48
         ibrav=3
      END SELECT simmetria
      RETURN
   END SUBROUTINE sym_brav
!# 736 "space_group.f90"
   SUBROUTINE find_equivalent_tau(space_group_number,inco,outco,i,unique)
!# 738 "space_group.f90"
   !sel_grup ->   input   space_group_number
   !         inco coordinate
   !         i element index
   !      output outco coordinates
!# 743 "space_group.f90"
      INTEGER, INTENT(IN) :: space_group_number,i
      REAL(DP),dimension(:,:), INTENT(IN) :: inco
      REAL(DP),dimension(:,:,:), INTENT(OUT) :: outco
      character(LEN=1), INTENT(IN) :: unique
!# 748 "space_group.f90"
      INTEGER :: k,j
      simmetria: SELECT CASE (space_group_number)
      !*****************************************
      !Triclinic 1-2
      CASE (1)
             CALL find_equiv_1  ( i, inco, outco )
      CASE (2)
             CALL find_equiv_2  ( i, inco, outco )
      CASE (3)
             CALL find_equiv_3  ( i, inco, unique, outco )
      CASE (4)
             CALL find_equiv_4  ( i, inco, unique, outco )
      CASE (5)
             CALL find_equiv_5  ( i, inco, unique, outco )
      CASE (6)
             CALL find_equiv_6  ( i, inco, unique, outco )
      CASE (7)
             CALL find_equiv_7  ( i, inco, unique, outco )
      CASE (8)
             CALL find_equiv_8  ( i, inco, unique, outco )
      CASE (9)
             CALL find_equiv_9  ( i, inco, unique, outco )
      CASE (10)
             CALL find_equiv_10 ( i, inco, unique, outco )
      CASE (11)
             CALL find_equiv_11 ( i, inco, unique, outco )
      CASE (12)
             CALL find_equiv_12 ( i, inco, unique, outco )
      CASE (13)
             CALL find_equiv_13 ( i, inco, unique, outco )
      CASE (14)
             CALL find_equiv_14 ( i, inco, unique, outco )
      CASE (15)
             CALL find_equiv_15 ( i, inco, unique, outco )
      CASE (16) !P222
             CALL find_equiv_16 ( i, inco, outco )
      CASE (17) !P222(1)
             CALL find_equiv_17 ( i, inco, outco )
      CASE (18) !P2(1)2(1)2
             CALL find_equiv_18 ( i, inco, outco )
      CASE (19) !P2(1)2(1)2(1)
             CALL find_equiv_19 ( i, inco, outco )
      CASE (20) !C222(1)
             CALL find_equiv_20 ( i, inco, outco )
      CASE (21) !C222
             CALL find_equiv_21 ( i, inco, outco )
      CASE (22) !F222
             CALL find_equiv_22 ( i, inco, outco )
      CASE (23) !I222
             CALL find_equiv_23 ( i, inco, outco )
      CASE (24) !I2(1)2(1)2(1)
             CALL find_equiv_24 ( i, inco, outco )
      CASE (25) !Pmm2
             CALL find_equiv_25 ( i, inco, outco )
      CASE (26) !Pmc2(1)
             CALL find_equiv_26 ( i, inco, outco )
      CASE (27) !Pcc2
             CALL find_equiv_27 ( i, inco, outco )
      CASE (28) !Pma2
             CALL find_equiv_28 ( i, inco, outco )
      CASE (29) !Pca2(1)
             CALL find_equiv_29 ( i, inco, outco )
      CASE (30) !Pnc2
             CALL find_equiv_30 ( i, inco, outco )
      CASE (31) !Pmn2(1)
             CALL find_equiv_31 ( i, inco, outco )
      CASE (32) !Pba2
             CALL find_equiv_32 ( i, inco, outco )
      CASE (33) !Pna2(1)
             CALL find_equiv_33 ( i, inco, outco )
      CASE (34) !Pnn2
             CALL find_equiv_34 ( i, inco, outco )
      CASE (35) !Cmm2
             CALL find_equiv_35 ( i, inco, outco )
      CASE (36) !Cmc2(1)
             CALL find_equiv_36 ( i, inco, outco )
      CASE (37) !Ccc2
             CALL find_equiv_37 ( i, inco, outco )
      CASE (38) !Amm2
             CALL find_equiv_38 ( i, inco, outco )
      CASE (39) !Abm2
             CALL find_equiv_39 ( i, inco, outco )
      CASE (40) !Ama2
             CALL find_equiv_40 ( i, inco, outco )
      CASE (41) !Aba2
             CALL find_equiv_41 ( i, inco, outco )
      CASE (42) !Fmm2
             CALL find_equiv_42 ( i, inco, outco )
      CASE (43) !Fdd2
             CALL find_equiv_43 ( i, inco, outco )
      CASE (44) !Imm2
             CALL find_equiv_44 ( i, inco, outco )
      CASE (45) !Iba2
             CALL find_equiv_45 ( i, inco, outco )
      CASE (46) !Ima2
             CALL find_equiv_46 ( i, inco, outco )
      CASE (47) !Pmmm
             CALL find_equiv_47 ( i, inco, outco )
      CASE (48) !Pnnn
             CALL find_equiv_48 ( i, inco, unique, outco )
      CASE (49) !Pccm
             CALL find_equiv_49 ( i, inco, outco )
      CASE (50) !Pban
             CALL find_equiv_50 ( i, inco, unique, outco )
      CASE (51) !Pmma
             CALL find_equiv_51 ( i, inco, outco )
      CASE (52) !Pnna
             CALL find_equiv_52 ( i, inco, outco )
      CASE (53) !Pmna
             CALL find_equiv_53 ( i, inco, outco )
      CASE (54) !Pcca
             CALL find_equiv_54 ( i, inco, outco )
      CASE (55) !Pbam
             CALL find_equiv_55 ( i, inco, outco )
      CASE (56) !Pccn
             CALL find_equiv_56 ( i, inco, outco )
      CASE (57) !Pbcm
             CALL find_equiv_57 ( i, inco, outco )
      CASE (58) !Pnnm
             CALL find_equiv_58 ( i, inco, outco )
      CASE (59) !Pmmn
             CALL find_equiv_59 ( i, inco, unique, outco )
      CASE (60) !Pbcn
             CALL find_equiv_60 ( i, inco, outco )
      CASE (61) !Pbca
             CALL find_equiv_61 ( i, inco, outco )
      CASE (62) !Pnma
             CALL find_equiv_62 ( i, inco, outco )
      CASE (63) !Cmcm
             CALL find_equiv_63 ( i, inco, outco )
      CASE (64) !Cmca
             CALL find_equiv_64 ( i, inco, outco )
      CASE (65) !Cmmm
             CALL find_equiv_65 ( i, inco, outco )
      CASE (66) !Cccm
             CALL find_equiv_66 ( i, inco, outco )
      CASE (67) !Cmma
             CALL find_equiv_67 ( i, inco, outco )
      CASE (68) !Ccca
             CALL find_equiv_68 ( i, inco, unique, outco )
      CASE (69) !Fmmm
             CALL find_equiv_69 ( i, inco, outco )
      CASE (70) !Fddd
             CALL find_equiv_70 ( i, inco, unique, outco )
      CASE (71) !Immm
             CALL find_equiv_71 ( i, inco, outco )
      CASE (72) !Ibam
             CALL find_equiv_72 ( i, inco, outco )
      CASE (73) !Ibca
             CALL find_equiv_73 ( i, inco, outco )
      CASE (74) !Imma
             CALL find_equiv_74 ( i, inco, outco )
      CASE (75) !P4
             CALL find_equiv_75 ( i, inco, outco )
      CASE (76) !P4(1)
             CALL find_equiv_76 ( i, inco, outco )
      CASE (77) !P4(2)
             CALL find_equiv_77 ( i, inco, outco )
      CASE (78) !P4(3)
             CALL find_equiv_78 ( i, inco, outco )
      CASE (79) !I4
             CALL find_equiv_79 ( i, inco, outco )
      CASE (80) !I4(1)
             CALL find_equiv_80 ( i, inco, outco )
      CASE (81) !P-4
             CALL find_equiv_81 ( i, inco, outco )
      CASE (82) !I-4
             CALL find_equiv_82 ( i, inco, outco )
      CASE (83) !P4/m
             CALL find_equiv_83 ( i, inco, outco )
      CASE (84) !P(2)/m
             CALL find_equiv_84 ( i, inco, outco )
      CASE (85) !P4/n
             CALL find_equiv_85 ( i, inco, unique, outco )
      CASE (86) !P4(2)/n
             CALL find_equiv_86 ( i, inco, unique, outco )
      CASE (87) !I4/m
             CALL find_equiv_87 ( i, inco, outco )
      CASE (88) !I4(1)/a
             CALL find_equiv_88 ( i, inco, unique, outco )
      CASE (89) !P422
             CALL find_equiv_89 ( i, inco, outco )
      CASE (90) !P42(1)2
             CALL find_equiv_90 ( i, inco, outco )
      CASE (91) !P4(1)22
             CALL find_equiv_91 ( i, inco, outco )
      CASE (92) !P4(1)2(1)2
             CALL find_equiv_92 ( i, inco, outco )
      CASE (93) !P4(2)22
             CALL find_equiv_93 ( i, inco, outco )
      CASE (94) !P4(2)2(1)2
             CALL find_equiv_94 ( i, inco, outco )
      CASE (95) !P4(3)22
             CALL find_equiv_95 ( i, inco, outco )
      CASE (96) !P4(3)2(1)2
             CALL find_equiv_96 ( i, inco, outco )
      CASE (97) !I422
             CALL find_equiv_97 ( i, inco, outco )
      CASE (98) !I4(1)22
             CALL find_equiv_98 ( i, inco, outco )
      CASE (99) !P4mm
             CALL find_equiv_99 ( i, inco, outco )
      CASE (100) !P4bm
             CALL find_equiv_100( i, inco, outco )
      CASE (101) !P4(2)cm
             CALL find_equiv_101( i, inco, outco )
      CASE (102) !P4(2)nm
             CALL find_equiv_102( i, inco, outco )
      CASE (103) !P4cc
             CALL find_equiv_103( i, inco, outco )
      CASE (104) !P4nc
             CALL find_equiv_104( i, inco, outco )
      CASE (105) !P4(2)mc
             CALL find_equiv_105( i, inco, outco )
      CASE (106) !P4(2)bc
             CALL find_equiv_106( i, inco, outco )
      CASE (107) !I4mm
             CALL find_equiv_107( i, inco, outco )
      CASE (108) !I4cm
             CALL find_equiv_108( i, inco, outco )
      CASE (109) !I4(1)md
             CALL find_equiv_109( i, inco, outco )
      CASE (110) !I4(1)cd
             CALL find_equiv_110( i, inco, outco )
      CASE (111) !P-42m
             CALL find_equiv_111( i, inco, outco )
      CASE (112) !P-42c
             CALL find_equiv_112( i, inco, outco )
      CASE (113) !P-42(1)m
             CALL find_equiv_113( i, inco, outco )
      CASE (114) !P-42(1)c
             CALL find_equiv_114( i, inco, outco )
      CASE (115) !P-4m2
             CALL find_equiv_115( i, inco, outco )
      CASE (116) !P-4c2
             CALL find_equiv_116( i, inco, outco )
      CASE (117) !P-4b2
             CALL find_equiv_117( i, inco, outco )
      CASE (118) !P-4n2
             CALL find_equiv_118( i, inco, outco )
      CASE (119) !I-4m2
             CALL find_equiv_119( i, inco, outco )
      CASE (120) !I-4c2
             CALL find_equiv_120( i, inco, outco )
      CASE (121) !I-42m
             CALL find_equiv_121( i, inco, outco )
      CASE (122) !I-42d
             CALL find_equiv_122( i, inco, outco )
      CASE (123) !P4/mmm
             CALL find_equiv_123( i, inco, outco )
      CASE (124) !P4/mcc
             CALL find_equiv_124( i, inco, outco )
      CASE (125) !P4/nbm
             CALL find_equiv_125( i, inco, unique, outco )
      CASE (126) !P4/nnc
             CALL find_equiv_126( i, inco, unique, outco )
      CASE (127) !P4/mbm
             CALL find_equiv_127( i, inco, outco )
      CASE (128) !P4/mnc
             CALL find_equiv_128( i, inco, outco )
      CASE (129)
             CALL find_equiv_129( i, inco, unique, outco )
      CASE (130) !P4/ncc
             CALL find_equiv_130( i, inco, unique, outco )
      CASE (131) !P4(2)/mmc
             CALL find_equiv_131( i, inco, outco )
      CASE (132) !P4(2)mcm
             CALL find_equiv_132( i, inco, outco )
      CASE (133) !P4(2)/nbc
             CALL find_equiv_133( i, inco, unique, outco )
      CASE (134) !P4(2)/nnm
             CALL find_equiv_134( i, inco, unique, outco )
      CASE (135) !P4(2)/mbc
             CALL find_equiv_135( i, inco, outco )
      CASE (136) !P4(2)mnm
             CALL find_equiv_136( i, inco, outco )
      CASE (137) !P4(2)/nmc
             CALL find_equiv_137( i, inco, unique, outco )
      CASE (138) !P4(2)/ncm
             CALL find_equiv_138( i, inco, unique, outco )
      CASE (139) !I4/mmm
             CALL find_equiv_139( i, inco, outco )
      CASE (140) !I4/mcm
             CALL find_equiv_140( i, inco, outco )
      CASE (141) !I4(1)amd
             CALL find_equiv_141( i, inco, unique, outco )
      CASE (142) !I4(1)/acd
             CALL find_equiv_142( i, inco, unique, outco )
      CASE (143) !P3
             CALL find_equiv_143( i, inco, outco )
      CASE (144) !P3(1)
             CALL find_equiv_144( i, inco, outco )
      CASE (145) !P3(2)
             CALL find_equiv_145( i, inco, outco )
      CASE (146) !R3
             CALL find_equiv_146( i, inco, unique, outco )
      CASE (147) !P-3
             CALL find_equiv_147( i, inco, outco )
      CASE (148) !R-3
             CALL find_equiv_148( i, inco, unique, outco )
      CASE (149) !P312
             CALL find_equiv_149( i, inco, outco )
      CASE (150) !P321
             CALL find_equiv_150( i, inco, outco )
      CASE (151) !P3(1)12
             CALL find_equiv_151( i, inco, outco )
      CASE (152) !P3(1)21
             CALL find_equiv_152( i, inco, outco )
      CASE (153) !P3(2)12
             CALL find_equiv_153( i, inco, outco )
      CASE (154) !P3(2)21
             CALL find_equiv_154( i, inco, outco )
      CASE (155) !R32
             CALL find_equiv_155( i, inco, unique, outco )
      CASE (156) !P3m1
             CALL find_equiv_156( i, inco, outco )
      CASE (157) !P31m
             CALL find_equiv_157( i, inco, outco )
      CASE (158) !P3c1
             CALL find_equiv_158( i, inco, outco )
      CASE (159) !P31c
             CALL find_equiv_159( i, inco, outco )
      CASE (160) !R3m
             CALL find_equiv_160( i, inco, unique, outco )
      CASE (161) !R3c
             CALL find_equiv_161( i, inco, unique, outco )
      CASE (162) !P-31m
             CALL find_equiv_162( i, inco, outco )
      CASE (163) !P-31c
             CALL find_equiv_163( i, inco, outco )
      CASE (164) !P-3m1
             CALL find_equiv_164( i, inco, outco )
      CASE (165) !P-3c1
             CALL find_equiv_165( i, inco, outco )
      CASE (166) !R-3m
             CALL find_equiv_166( i, inco, unique, outco )
      CASE (167) !R-3c
             CALL find_equiv_167( i, inco, unique, outco )
      CASE (168) !P6
             CALL find_equiv_168( i, inco, outco )
      CASE (169) !P6(1)
             CALL find_equiv_169( i, inco, outco )
      CASE (170) !P6(5)
             CALL find_equiv_170( i, inco, outco )
      CASE (171) !P6(2)
             CALL find_equiv_171( i, inco, outco )
      CASE (172) !P6(4)
             CALL find_equiv_172( i, inco, outco )
      CASE (173) !P6(3)
             CALL find_equiv_173( i, inco, outco )
      CASE (174) !P-6
             CALL find_equiv_174( i, inco, outco )
      CASE (175) !P6/m
             CALL find_equiv_175( i, inco, outco )
      CASE (176) !P6(3)/m
             CALL find_equiv_176( i, inco, outco )
      CASE (177) !P622
             CALL find_equiv_177( i, inco, outco )
      CASE (178) !P(1)22
             CALL find_equiv_178( i, inco, outco )
      CASE (179) !P6(5)22
             CALL find_equiv_179( i, inco, outco )
      CASE (180) !P6(2)22
             CALL find_equiv_180( i, inco, outco )
      CASE (181) !P6(4)22
             CALL find_equiv_181( i, inco, outco )
      CASE (182) !6(3)22
             CALL find_equiv_182( i, inco, outco )
      CASE (183) !P6mm
             CALL find_equiv_183( i, inco, outco )
      CASE (184) !P6cc
             CALL find_equiv_184( i, inco, outco )
      CASE (185) !P6(3)cm
             CALL find_equiv_185( i, inco, outco )
      CASE (186) !P(3)mc
             CALL find_equiv_186( i, inco, outco )
      CASE (187) !P-6m2
             CALL find_equiv_187( i, inco, outco )
      CASE (188) !P-6c2
             CALL find_equiv_188( i, inco, outco )
      CASE (189) !P-62m
             CALL find_equiv_189( i, inco, outco )
      CASE (190) !P-62c
             CALL find_equiv_190( i, inco, outco )
      CASE (191) !P6/mmm
             CALL find_equiv_191( i, inco, outco )
      CASE (192) !P6/mmc
             CALL find_equiv_192( i, inco, outco )
      CASE (193) !P6(3)/mcm
             CALL find_equiv_193( i, inco, outco )
      CASE (194)
             CALL find_equiv_194( i, inco, outco )
      CASE (195) !P23
             CALL find_equiv_195( i, inco, outco )
      CASE (196) !F23
             CALL find_equiv_196( i, inco, outco )
      CASE (197) !I23
             CALL find_equiv_197( i, inco, outco )
      CASE (198) !P2(1)3
             CALL find_equiv_198( i, inco, outco )
      CASE (199) !I2(1)3
             CALL find_equiv_199( i, inco, outco )
      CASE (200) !Pm-3
             CALL find_equiv_200( i, inco, outco )
      CASE(201) !Pn-3
             CALL find_equiv_201( i, inco, unique, outco )
      CASE (202) !Fm-3
             CALL find_equiv_202( i, inco, outco )
      CASE (203) !Fd-3
             CALL find_equiv_203( i, inco, unique, outco )
      CASE (204) !Im-3
             CALL find_equiv_204( i, inco, outco )
      CASE (205) !Pa-3
             CALL find_equiv_205( i, inco, outco )
      CASE (206) !Ia-3
             CALL find_equiv_206( i, inco, outco )
      CASE (207) !P432
             CALL find_equiv_207( i, inco, outco )
      CASE (208) !P4(2)32
             CALL find_equiv_208( i, inco, outco )
      CASE (209) !F432
             CALL find_equiv_209( i, inco, outco )
      CASE (210) !F4(1)32
             CALL find_equiv_210( i, inco, outco )
      CASE (211) !I432
             CALL find_equiv_211( i, inco, outco )
      CASE (212) !P4(3)32
             CALL find_equiv_212( i, inco, outco )
      CASE (213) !P4(1)32
             CALL find_equiv_213( i, inco, outco )
      CASE (214) !I4(1)32
             CALL find_equiv_214( i, inco, outco )
      CASE (215) !P-43m
             CALL find_equiv_215( i, inco, outco )
      CASE (216) !F-43m
             CALL find_equiv_216( i, inco, outco )
      CASE (217) !I-43m
             CALL find_equiv_217( i, inco, outco )
      CASE (218) !P-43n
             CALL find_equiv_218( i, inco, outco )
      CASE (219) !F-43c
             CALL find_equiv_219( i, inco, outco )
      CASE (220) !I-43d
             CALL find_equiv_220( i, inco, outco )
      CASE (221) !Pm-3m
             CALL find_equiv_221( i, inco, outco )
      CASE (222) !Pn-3n
             CALL find_equiv_222( i, inco, unique, outco )
      CASE (223) !Pm-3n
             CALL find_equiv_223( i, inco, outco )
      CASE (224) !Pn-3m
             CALL find_equiv_224( i, inco, unique, outco )
      CASE (225) !Fm-3m
             CALL find_equiv_225( i, inco, outco )
      CASE (226) !Fm-3c
             CALL find_equiv_226( i, inco, outco )
      CASE (227) !Fd-3m
             CALL find_equiv_227( i, inco, unique, outco )
      CASE (228) !Fd-3c
             CALL find_equiv_228( i, inco, unique, outco )
      CASE (229)
             CALL find_equiv_229( i, inco, outco )
      CASE (230)
             CALL find_equiv_230( i, inco, outco )
      END SELECT simmetria
!# 1214 "space_group.f90"
    END SUBROUTINE find_equivalent_tau
!# 1216 "space_group.f90"
SUBROUTINE find_equiv_1  ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1222 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
END SUBROUTINE find_equiv_1  
!# 1227 "space_group.f90"
SUBROUTINE find_equiv_2  ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1233 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         outco(k,2,i)=-inco(k,i)
         END DO
      !*****************************************
      !Monoclinic 3-15
END SUBROUTINE find_equiv_2  
!# 1241 "space_group.f90"
SUBROUTINE find_equiv_3  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1248 "space_group.f90"
         !x,y,z
         !-x,y,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1254 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
!# 1260 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
END SUBROUTINE find_equiv_3  
!# 1267 "space_group.f90"
SUBROUTINE find_equiv_4  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1274 "space_group.f90"
         !x,y,z
         !-X,Y+1/2,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1280 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)+0.5_DP
         outco(3,2,i)=-inco(3,i)
         END IF
!# 1286 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         END IF
!# 1292 "space_group.f90"
END SUBROUTINE find_equiv_4  
!# 1294 "space_group.f90"
SUBROUTINE find_equiv_5  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1301 "space_group.f90"
         !X,Y,Z identita
         !-X,Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
!# 1313 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
END SUBROUTINE find_equiv_5  
!# 1320 "space_group.f90"
SUBROUTINE find_equiv_6  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1327 "space_group.f90"
         !ID
         !x,-y,z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1333 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
!# 1339 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_6  
!# 1346 "space_group.f90"
SUBROUTINE find_equiv_7  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1353 "space_group.f90"
         !ID
         !x,-y,1/2+z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1359 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         END IF
!# 1365 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=0.5_DP+inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_7  
!# 1372 "space_group.f90"
SUBROUTINE find_equiv_8  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1379 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
!# 1386 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
!# 1392 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_8  
!# 1399 "space_group.f90"
SUBROUTINE find_equiv_9  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1406 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= X,-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
!# 1413 "space_group.f90"
         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         END IF
!# 1419 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)+0.5_DP
         outco(3,2,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_9  
!# 1426 "space_group.f90"
SUBROUTINE find_equiv_10 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1433 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= -X,-Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1441 "space_group.f90"
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
!# 1456 "space_group.f90"
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_10 
!# 1472 "space_group.f90"
SUBROUTINE find_equiv_11 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1479 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,1/2+Y,-Z
         !symmetry= -X,-Y,-Z
         !symmetry= X,1/2-Y,Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP+inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=inco(3,i)
         END IF
!# 1501 "space_group.f90"
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=inco(2,i)
         outco(3,4,i)=0.5_DP-inco(3,i)
         END IF
END SUBROUTINE find_equiv_11 
!# 1517 "space_group.f90"
SUBROUTINE find_equiv_12 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1524 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= -X,-Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1532 "space_group.f90"
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
!# 1547 "space_group.f90"
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_12 
!# 1562 "space_group.f90"
SUBROUTINE find_equiv_13 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1569 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= -X,-Y,-Z
         !symmetry= X,-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1577 "space_group.f90"
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=0.5_DP-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF
!# 1592 "space_group.f90"
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_13 
!# 1608 "space_group.f90"
SUBROUTINE find_equiv_14 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1615 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,-Z
         !symmetry= -X,1/2+Y,1/2-Z
         !symmetry= X,1/2-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1623 "space_group.f90"
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF
!# 1638 "space_group.f90"
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP-inco(2,i)
         outco(3,3,i)=0.5_DP+inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=0.5_DP-inco(3,i)
         END IF
END SUBROUTINE find_equiv_14 
!# 1654 "space_group.f90"
SUBROUTINE find_equiv_15 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1661 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= -X,-Y,-Z
         !symmetry= X,-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
!# 1669 "space_group.f90"
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=0.5_DP-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=3
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF
!# 1684 "space_group.f90"
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=3
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
!# 1699 "space_group.f90"
      !*****************************************
      !Orthorhombic 16-74
END SUBROUTINE find_equiv_15 
!# 1703 "space_group.f90"
SUBROUTINE find_equiv_16 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1709 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
END SUBROUTINE find_equiv_16 
!# 1730 "space_group.f90"
SUBROUTINE find_equiv_17 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1736 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,1/2+Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= X,-Y,-Z
!# 1741 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
END SUBROUTINE find_equiv_17 
!# 1758 "space_group.f90"
SUBROUTINE find_equiv_18 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1764 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= 1/2-X,1/2+Y,-Z
         !symmetry= 1/2+X,1/2-Y,-Z
!# 1769 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=0.5_DP-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=0.5_DP+inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=-inco(3,i)
!# 1785 "space_group.f90"
END SUBROUTINE find_equiv_18 
!# 1787 "space_group.f90"
SUBROUTINE find_equiv_19 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1793 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= 1/2-X,-Y,1/2+Z
         !symmetry= -X,1/2+Y,1/2-Z
         !symmetry= 1/2+X,1/2-Y,-Z
!# 1798 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=0.5_DP-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=0.5_DP+inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=-inco(3,i)
END SUBROUTINE find_equiv_19 
!# 1815 "space_group.f90"
SUBROUTINE find_equiv_20 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1822 "space_group.f90"
         ! symmetry= X,Y,Z
         !symmetry= -X,-Y,1/2+Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= X,-Y,-Z
!# 1827 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
!# 1843 "space_group.f90"
END SUBROUTINE find_equiv_20 
!# 1845 "space_group.f90"
SUBROUTINE find_equiv_21 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1851 "space_group.f90"
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z
!# 1856 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
!# 1873 "space_group.f90"
END SUBROUTINE find_equiv_21 
!# 1875 "space_group.f90"
SUBROUTINE find_equiv_22 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1881 "space_group.f90"
         ! symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z
!# 1886 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
END SUBROUTINE find_equiv_22 
!# 1903 "space_group.f90"
SUBROUTINE find_equiv_23 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1909 "space_group.f90"
         !id
         !-x,-y,z
         !x,,y,-z
         !x,-y,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
END SUBROUTINE find_equiv_23 
!# 1930 "space_group.f90"
SUBROUTINE find_equiv_24 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1936 "space_group.f90"
         !id
         !-x+1/2,-y,z+1/2
         !-x,1/2+y,1/2-z
         !x+1/2,-y+1/2,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
!# 1956 "space_group.f90"
END SUBROUTINE find_equiv_24 
!# 1958 "space_group.f90"
SUBROUTINE find_equiv_25 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1964 "space_group.f90"
         !id
         !-x,-y,z
         !+x,-y,+z
         !-x,y,z
!# 1969 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 1985 "space_group.f90"
END SUBROUTINE find_equiv_25 
!# 1987 "space_group.f90"
SUBROUTINE find_equiv_26 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 1993 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !+x,-y,+z+1/2
         !-x,y,z
!# 1998 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2014 "space_group.f90"
END SUBROUTINE find_equiv_26 
!# 2016 "space_group.f90"
SUBROUTINE find_equiv_27 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2022 "space_group.f90"
         !id
         !-x,-y,z
         !+x,-y,+z+1/2
         !-x,y,z+1/2
!# 2027 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2043 "space_group.f90"
END SUBROUTINE find_equiv_27 
!# 2045 "space_group.f90"
SUBROUTINE find_equiv_28 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2051 "space_group.f90"
         !id
         !-x,-y,z
         !1/2+x,-y,z
         !1/2-x,y,z
!# 2056 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2072 "space_group.f90"
END SUBROUTINE find_equiv_28 
!# 2074 "space_group.f90"
SUBROUTINE find_equiv_29 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2080 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !1/2+x,-y,z
         !1/2-x,y,z+1/2
!# 2085 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2101 "space_group.f90"
END SUBROUTINE find_equiv_29 
!# 2103 "space_group.f90"
SUBROUTINE find_equiv_30 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2109 "space_group.f90"
         !id
         !-x,-y,z
         !+x,1/2-y,z+1/2
         !-x,y+1/2,z+1/2
!# 2114 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2130 "space_group.f90"
END SUBROUTINE find_equiv_30 
!# 2132 "space_group.f90"
SUBROUTINE find_equiv_31 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2138 "space_group.f90"
         !id
         !1/2-x,-y,z+1/2
         !1/2+x,-y,z+1/2
         !-x,y,z
!# 2143 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2159 "space_group.f90"
END SUBROUTINE find_equiv_31 
!# 2161 "space_group.f90"
SUBROUTINE find_equiv_32 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2167 "space_group.f90"
         !id
         !-x,-y,z
         !1/2+x,1/2-y,z
         !1/2-x,1/2+y,z
!# 2172 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
!# 2188 "space_group.f90"
END SUBROUTINE find_equiv_32 
!# 2190 "space_group.f90"
SUBROUTINE find_equiv_33 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2196 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !1/2+x,1/2-y,z
         !1/2-x,1/2+y,z+1/2
!# 2201 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2217 "space_group.f90"
END SUBROUTINE find_equiv_33 
!# 2219 "space_group.f90"
SUBROUTINE find_equiv_34 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2225 "space_group.f90"
         !id
         !-x,-y,z
         !1/2+x,1/2-y,1/2+z
         !1/2-x,1/2+y,1/2+z
!# 2230 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2246 "space_group.f90"
END SUBROUTINE find_equiv_34 
!# 2248 "space_group.f90"
SUBROUTINE find_equiv_35 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2254 "space_group.f90"
         !id
         !-x,-y,z
         !+x,-y,z
         !-x,+y,z
!# 2259 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2275 "space_group.f90"
END SUBROUTINE find_equiv_35 
!# 2277 "space_group.f90"
SUBROUTINE find_equiv_36 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2283 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !+x,-y,z+1/2
         !-x,+y,z
!# 2288 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2304 "space_group.f90"
END SUBROUTINE find_equiv_36 
!# 2306 "space_group.f90"
SUBROUTINE find_equiv_37 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2312 "space_group.f90"
         !id
         !-x,-y,z
         !+x,-y,z+1/2
         !-x,+y,z+1/2
!# 2317 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 2333 "space_group.f90"
END SUBROUTINE find_equiv_37 
!# 2335 "space_group.f90"
SUBROUTINE find_equiv_38 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2341 "space_group.f90"
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z
!# 2346 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2362 "space_group.f90"
END SUBROUTINE find_equiv_38 
!# 2364 "space_group.f90"
SUBROUTINE find_equiv_39 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2370 "space_group.f90"
         !id
         !-x,-y,z
         !x,-y+1/2,z
         !-x,y+1/2,z
!# 2375 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
!# 2391 "space_group.f90"
END SUBROUTINE find_equiv_39 
!# 2393 "space_group.f90"
SUBROUTINE find_equiv_40 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2399 "space_group.f90"
         !id
         !-x,-y,z
         !x+1/2,-y,z
         !-x+1/2,y,z
!# 2404 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2420 "space_group.f90"
END SUBROUTINE find_equiv_40 
!# 2422 "space_group.f90"
SUBROUTINE find_equiv_41 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2428 "space_group.f90"
         !id
         !-x,-y,z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
!# 2433 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
!# 2449 "space_group.f90"
END SUBROUTINE find_equiv_41 
!# 2451 "space_group.f90"
SUBROUTINE find_equiv_42 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2457 "space_group.f90"
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z
!# 2462 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2478 "space_group.f90"
END SUBROUTINE find_equiv_42 
!# 2480 "space_group.f90"
SUBROUTINE find_equiv_43 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2486 "space_group.f90"
         !id
         !-x,-y,z
         !x+1/4,-y+1/4,z+1/4
         !-x+1/4,y+1/4,z+1/4
!# 2491 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.25_DP
         outco(2,3,i)=-inco(2,i)+0.25_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=-inco(1,i)+0.25_DP
         outco(2,4,i)=+inco(2,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.25_DP
!# 2507 "space_group.f90"
END SUBROUTINE find_equiv_43 
!# 2509 "space_group.f90"
SUBROUTINE find_equiv_44 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2515 "space_group.f90"
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z
!# 2520 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2536 "space_group.f90"
END SUBROUTINE find_equiv_44 
!# 2538 "space_group.f90"
SUBROUTINE find_equiv_45 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2544 "space_group.f90"
         !id
         !-x,-y,z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
!# 2549 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
!# 2565 "space_group.f90"
END SUBROUTINE find_equiv_45 
!# 2567 "space_group.f90"
SUBROUTINE find_equiv_46 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2573 "space_group.f90"
         !id
         !-x,-y,z
         !x+1/2,-y,z
         !-x+1/2,y,z
!# 2578 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
!# 2594 "space_group.f90"
END SUBROUTINE find_equiv_46 
!# 2596 "space_group.f90"
SUBROUTINE find_equiv_47 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2602 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z
!# 2611 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 2643 "space_group.f90"
END SUBROUTINE find_equiv_47 
!# 2645 "space_group.f90"
SUBROUTINE find_equiv_48 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2653 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !x+1/2,-y+/2,z+1/2
         !-x+1/2,y+1/2,z+1/2
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF
!# 2695 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x+1/2,+y,-z+1/2
         !+x,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y,z+1/2
         !-x,y+1/2,z+1/2
!# 2705 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF
!# 2738 "space_group.f90"
END SUBROUTINE find_equiv_48 
!# 2740 "space_group.f90"
SUBROUTINE find_equiv_49 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2746 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x,-y,z+1/2
         !-x,y,z+1/2
!# 2755 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 2787 "space_group.f90"
END SUBROUTINE find_equiv_49 
!# 2789 "space_group.f90"
SUBROUTINE find_equiv_50 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2797 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
!# 2807 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
         END IF
!# 2840 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x+1/2,+y,-z
         !+x,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y,z
         !-x,y+1/2,z
!# 2850 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
         END IF
!# 2883 "space_group.f90"
END SUBROUTINE find_equiv_50 
!# 2885 "space_group.f90"
SUBROUTINE find_equiv_51 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2891 "space_group.f90"
         !id
         !-x+1/2,-y,z
         !-x,+y,-z
         !+x+1/2,-y,-z
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y,z
         !-x+1/2,y,z
!# 2900 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 2932 "space_group.f90"
END SUBROUTINE find_equiv_51 
!# 2934 "space_group.f90"
SUBROUTINE find_equiv_52 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2940 "space_group.f90"
         !id
         !-x+1/2,-y,z
         !-x+1/2,+y+1/2,-z+1/2
         !+x,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x+1/2,-y+1/2,z+1/2
         !-x,y+1/2,z+1/2
!# 2949 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 2981 "space_group.f90"
END SUBROUTINE find_equiv_52 
!# 2983 "space_group.f90"
SUBROUTINE find_equiv_53 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 2989 "space_group.f90"
         !id
         !-x+1/2,-y,z+1/2
         !-x+1/2,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x+1/2,-y,z+1/2
         !-x,y,z
!# 2998 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 3030 "space_group.f90"
END SUBROUTINE find_equiv_53 
!# 3032 "space_group.f90"
SUBROUTINE find_equiv_54 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3038 "space_group.f90"
         !id
         !-x+1/2,-y,z
         !-x,+y,-z+1/2
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y,z+1/2
         !-x+1/2,y,z+1/2
!# 3047 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 3079 "space_group.f90"
END SUBROUTINE find_equiv_54 
!# 3081 "space_group.f90"
SUBROUTINE find_equiv_55 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3087 "space_group.f90"
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
!# 3096 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 3128 "space_group.f90"
END SUBROUTINE find_equiv_55 
!# 3130 "space_group.f90"
SUBROUTINE find_equiv_56 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3136 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x,-y+1/2,z+1/2
         !-x+1/2,y,z+1/2
!# 3145 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 3177 "space_group.f90"
END SUBROUTINE find_equiv_56 
!# 3179 "space_group.f90"
SUBROUTINE find_equiv_57 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3185 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x,y+1/2,z
!# 3194 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 3226 "space_group.f90"
END SUBROUTINE find_equiv_57 
!# 3228 "space_group.f90"
SUBROUTINE find_equiv_58 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3234 "space_group.f90"
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z+1/2
!# 3243 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 3275 "space_group.f90"
END SUBROUTINE find_equiv_58 
!# 3277 "space_group.f90"
SUBROUTINE find_equiv_59 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3285 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !x,-y,z
         !-x,y,z
!# 3295 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
         END IF
!# 3328 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y+1/2,-z
         !+x+1/2,-y,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x,-y+1/2,z
         !-x+1/2,y,z
!# 3338 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
         END IF
!# 3371 "space_group.f90"
END SUBROUTINE find_equiv_59 
!# 3373 "space_group.f90"
SUBROUTINE find_equiv_60 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3379 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-x,+y,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z+1/2
         !x,-y,z+1/2
         !-x+1/2,y+1/2,z
!# 3388 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 3420 "space_group.f90"
END SUBROUTINE find_equiv_60 
!# 3422 "space_group.f90"
SUBROUTINE find_equiv_61 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3428 "space_group.f90"
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z
!# 3437 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 3469 "space_group.f90"
END SUBROUTINE find_equiv_61 
!# 3471 "space_group.f90"
SUBROUTINE find_equiv_62 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3477 "space_group.f90"
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z
         !+x+1/2,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z
         !-x+1/2,y+1/2,z+1/2
!# 3486 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 3518 "space_group.f90"
END SUBROUTINE find_equiv_62 
!# 3520 "space_group.f90"
SUBROUTINE find_equiv_63 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3526 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-x,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z+1/2
         !x,-y,z+1/2
         !-x,y,z
!# 3535 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 3567 "space_group.f90"
END SUBROUTINE find_equiv_63 
!# 3569 "space_group.f90"
SUBROUTINE find_equiv_64 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3575 "space_group.f90"
         !id
         !-x,-y+1/2,z+1/2
         !-x,+y,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z+1/2
         !x,-y,z+1/2
         !-x+1/2,y+1/2,z
!# 3584 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 3616 "space_group.f90"
END SUBROUTINE find_equiv_64 
!# 3618 "space_group.f90"
SUBROUTINE find_equiv_65 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3624 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z
!# 3633 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 3665 "space_group.f90"
END SUBROUTINE find_equiv_65 
!# 3667 "space_group.f90"
SUBROUTINE find_equiv_66 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3673 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x,-y,z+1/2
         !-x,y,z+1/2
!# 3682 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 3714 "space_group.f90"
END SUBROUTINE find_equiv_66 
!# 3716 "space_group.f90"
SUBROUTINE find_equiv_67 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3722 "space_group.f90"
         !id
         !-x,-y+1/2,z
         !-x,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x,y+1/2,-z
         !x,-y+1/2,z
         !-x,y,z
!# 3731 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 3763 "space_group.f90"
END SUBROUTINE find_equiv_67 
!# 3765 "space_group.f90"
SUBROUTINE find_equiv_68 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3773 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y+1/2,-z+1/2
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y,z+1/2
!# 3783 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF
!# 3816 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y,-z
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y+,z+1/2
         !-x+1/2,y,z+1/2
!# 3826 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF
!# 3859 "space_group.f90"
END SUBROUTINE find_equiv_68 
!# 3861 "space_group.f90"
SUBROUTINE find_equiv_69 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3867 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z
!# 3876 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 3908 "space_group.f90"
END SUBROUTINE find_equiv_69 
!# 3910 "space_group.f90"
SUBROUTINE find_equiv_70 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 3918 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/4,-y+1/4,-z+1/4
         !x+1/4,y+1/4,-z+1/4
         !x+1/4,-y+1/4,z+1/4
         !-x+1/4,y+1/4,z+1/4
!# 3928 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.25_DP
         outco(2,5,i)=-inco(2,i)+0.25_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.25_DP
         outco(2,6,i)=+inco(2,i)+0.25_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.25_DP
         outco(2,7,i)=-inco(2,i)+0.25_DP
         outco(3,7,i)=+inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.25_DP
         outco(2,8,i)=+inco(2,i)+0.25_DP
         outco(3,8,i)=+inco(3,i)+0.25_DP
         END IF
!# 3961 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+3/4,-y+3/4,z
         !-x+3/4,+y,-z+3/4
         !+x,-y+3/4,-z+3/4
         !-x,-y,-z
         !x+3/4,y+3/4,-z
         !x+3/4,-y,z+3/4
         !-x,y+3/4,z+3/4
!# 3971 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.75_DP
         outco(2,6,i)=+inco(2,i)+0.75_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.75_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.75_DP
         outco(3,8,i)=+inco(3,i)+0.75_DP
         END IF
!# 4004 "space_group.f90"
END SUBROUTINE find_equiv_70 
!# 4006 "space_group.f90"
SUBROUTINE find_equiv_71 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4012 "space_group.f90"
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z
!# 4021 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 4053 "space_group.f90"
END SUBROUTINE find_equiv_71 
!# 4055 "space_group.f90"
SUBROUTINE find_equiv_72 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4061 "space_group.f90"
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
!# 4070 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 4102 "space_group.f90"
END SUBROUTINE find_equiv_72 
!# 4104 "space_group.f90"
SUBROUTINE find_equiv_73 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4110 "space_group.f90"
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z
!# 4119 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 4151 "space_group.f90"
END SUBROUTINE find_equiv_73 
!# 4153 "space_group.f90"
SUBROUTINE find_equiv_74 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4159 "space_group.f90"
         !id
         !-x,-y+1/2,z
         !-x,+y+1/2,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y+1/2,-z
         !x,-y+1/2,z
         !-x,y,z
!# 4168 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
!# 4200 "space_group.f90"
      !*****************************************
      !Tetragonal 75-142
!# 4203 "space_group.f90"
END SUBROUTINE find_equiv_74 
!# 4205 "space_group.f90"
SUBROUTINE find_equiv_75 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4211 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
!# 4216 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
!# 4232 "space_group.f90"
END SUBROUTINE find_equiv_75 
!# 4234 "space_group.f90"
SUBROUTINE find_equiv_76 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4240 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
!# 4245 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
!# 4261 "space_group.f90"
END SUBROUTINE find_equiv_76 
!# 4263 "space_group.f90"
SUBROUTINE find_equiv_77 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4269 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
!# 4274 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
!# 4290 "space_group.f90"
END SUBROUTINE find_equiv_77 
!# 4292 "space_group.f90"
SUBROUTINE find_equiv_78 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4298 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y,x,z+3/4
         !y,-x,z+1/4
!# 4303 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.25_DP
!# 4319 "space_group.f90"
END SUBROUTINE find_equiv_78 
!# 4321 "space_group.f90"
SUBROUTINE find_equiv_79 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4327 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
!# 4332 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
!# 4348 "space_group.f90"
END SUBROUTINE find_equiv_79 
!# 4350 "space_group.f90"
SUBROUTINE find_equiv_80 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4356 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
!# 4361 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
!# 4377 "space_group.f90"
END SUBROUTINE find_equiv_80 
!# 4379 "space_group.f90"
SUBROUTINE find_equiv_81 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4385 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,x,-z
!# 4390 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
!# 4406 "space_group.f90"
END SUBROUTINE find_equiv_81 
!# 4408 "space_group.f90"
SUBROUTINE find_equiv_82 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4414 "space_group.f90"
         !id
         !-x,-y,z
         !+y,-x,-z
         !-y,x,-z
!# 4419 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
!# 4435 "space_group.f90"
END SUBROUTINE find_equiv_82 
!# 4437 "space_group.f90"
SUBROUTINE find_equiv_83 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4443 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x-z
!# 4452 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 4484 "space_group.f90"
END SUBROUTINE find_equiv_83 
!# 4486 "space_group.f90"
SUBROUTINE find_equiv_84 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4492 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z+1/2
         !-y,x-z+1/2
!# 4501 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 4533 "space_group.f90"
END SUBROUTINE find_equiv_84 
!# 4535 "space_group.f90"
SUBROUTINE find_equiv_85 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4543 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z
         !y+1/2,-x+1/2,z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y,-x,-z
         !-y,x-z
!# 4553 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
         END IF
!# 4586 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y+1/2,x,z
         !y,-x+1/2,z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x,-z
         !-y,x+1/2,-z
!# 4596 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         END IF
!# 4629 "space_group.f90"
END SUBROUTINE find_equiv_85 
!# 4631 "space_group.f90"
SUBROUTINE find_equiv_86 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4638 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !y,-x,-z
         !-y,x-z
!# 4648 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
         END IF
!# 4681 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y,x+1/2,z+1/2
         !y+1/2,-x,z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y,-x+1/2,-z+1/2
         !-y+1/2,x,-z+1/2
!# 4691 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         END IF
!# 4724 "space_group.f90"
END SUBROUTINE find_equiv_86 
!# 4726 "space_group.f90"
SUBROUTINE find_equiv_87 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4732 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x-z
!# 4741 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 4773 "space_group.f90"
END SUBROUTINE find_equiv_87 
!# 4775 "space_group.f90"
SUBROUTINE find_equiv_88 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4782 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !-x,-y+1/2,-z+1/4
         !x+1/2,y,-z+3/4
         !y,-x,-z
         !-y+1/2,x+1/2,-z+1/2
!# 4792 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         END IF
!# 4825 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y,z+1/2
         !-y+3/4,x+1/4,z+1/4
         !y+3/4,-x+3/4,z+3/4
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !y+1/4,-x+3/4,-z+3/4
         !-y+1/4,x+1/4,-z+1/4
!# 4835 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.75_DP
         outco(2,3,i)=+inco(1,i)+0.25_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.75_DP
         outco(2,4,i)=-inco(1,i)+0.75_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=+inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.25_DP
         END IF
!# 4868 "space_group.f90"
END SUBROUTINE find_equiv_88 
!# 4870 "space_group.f90"
SUBROUTINE find_equiv_89 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4876 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z
         !-y,-x-z
!# 4885 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 4917 "space_group.f90"
END SUBROUTINE find_equiv_89 
!# 4919 "space_group.f90"
SUBROUTINE find_equiv_90 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4925 "space_group.f90"
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z
         !y+1/2,-x+1/2,z
         !-x+1/2,+y+1/2,-z
         !x+1/2,-y+1/2,-z
         !y,x,-z
         !-y,-x-z
!# 4934 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 4966 "space_group.f90"
END SUBROUTINE find_equiv_90 
!# 4968 "space_group.f90"
SUBROUTINE find_equiv_91 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 4974 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
!# 4983 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.25_DP
!# 5015 "space_group.f90"
END SUBROUTINE find_equiv_91 
!# 5017 "space_group.f90"
SUBROUTINE find_equiv_92 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5023 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y+1/2,x+1/2,z+1/4
         !y+1/2,-x+1/2,z+3/4
         !-x+1/2,+y+1/2,-z+1/4
         !x+1/2,-y+1/2,-z+3/4
         !y,x,-z
         !-y,-x-z+1/2
!# 5032 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 5064 "space_group.f90"
END SUBROUTINE find_equiv_92 
!# 5066 "space_group.f90"
SUBROUTINE find_equiv_93 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5072 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z+1/2
         !-y,-x-z+1/2
!# 5081 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 5113 "space_group.f90"
END SUBROUTINE find_equiv_93 
!# 5115 "space_group.f90"
SUBROUTINE find_equiv_94 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5121 "space_group.f90"
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !-x+1/2,+y+1/2,-z+1/2
         !x+1/2,-y+1/2,-z+1/2
         !y,x,-z
         !-y,-x,-z
!# 5130 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 5162 "space_group.f90"
END SUBROUTINE find_equiv_94 
!# 5164 "space_group.f90"
SUBROUTINE find_equiv_95 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5170 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y,x,z+3/4
         !y,-x,z+1/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+1/4
         !-y,-x-z+3/4
!# 5179 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.75_DP
!# 5211 "space_group.f90"
END SUBROUTINE find_equiv_95 
!# 5213 "space_group.f90"
SUBROUTINE find_equiv_96 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5219 "space_group.f90"
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
!# 5228 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 5260 "space_group.f90"
END SUBROUTINE find_equiv_96 
!# 5262 "space_group.f90"
SUBROUTINE find_equiv_97 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5268 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z
         !-y,-x-z
!# 5277 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 5309 "space_group.f90"
END SUBROUTINE find_equiv_97 
!# 5311 "space_group.f90"
SUBROUTINE find_equiv_98 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5317 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
!# 5326 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 5358 "space_group.f90"
END SUBROUTINE find_equiv_98 
!# 5360 "space_group.f90"
SUBROUTINE find_equiv_99 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5366 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z
         !y,x,z
!# 5375 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 5407 "space_group.f90"
END SUBROUTINE find_equiv_99 
!# 5409 "space_group.f90"
SUBROUTINE find_equiv_100( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5415 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
!# 5424 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 5456 "space_group.f90"
END SUBROUTINE find_equiv_100
!# 5458 "space_group.f90"
SUBROUTINE find_equiv_101( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5464 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !-y,-x,+z
         !y,x,z
!# 5473 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 5505 "space_group.f90"
END SUBROUTINE find_equiv_101
!# 5507 "space_group.f90"
SUBROUTINE find_equiv_102( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5513 "space_group.f90"
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !-y,-x,+z
         !y,x,z
!# 5522 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 5554 "space_group.f90"
END SUBROUTINE find_equiv_102
!# 5556 "space_group.f90"
SUBROUTINE find_equiv_103( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5562 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !-y,-x,+z+1/2
         !y,x,z+1/2
!# 5571 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 5603 "space_group.f90"
END SUBROUTINE find_equiv_103
!# 5605 "space_group.f90"
SUBROUTINE find_equiv_104( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5611 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
!# 5620 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 5652 "space_group.f90"
END SUBROUTINE find_equiv_104
!# 5654 "space_group.f90"
SUBROUTINE find_equiv_105( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5660 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z+1/2
         !y,x,z+1/2
!# 5669 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 5701 "space_group.f90"
END SUBROUTINE find_equiv_105
!# 5703 "space_group.f90"
SUBROUTINE find_equiv_106( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5709 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
!# 5718 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 5750 "space_group.f90"
END SUBROUTINE find_equiv_106
!# 5752 "space_group.f90"
SUBROUTINE find_equiv_107( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5758 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z
         !y,x,z
!# 5767 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 5799 "space_group.f90"
END SUBROUTINE find_equiv_107
!# 5801 "space_group.f90"
SUBROUTINE find_equiv_108( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5807 "space_group.f90"
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
!# 5816 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 5848 "space_group.f90"
END SUBROUTINE find_equiv_108
!# 5850 "space_group.f90"
SUBROUTINE find_equiv_109( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5856 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !+x,-y,+z
         !-x+1/2,+y+1/2,+z+1/2
         !-y,-x+1/2,+z+1/2
         !y+1/2,x,z+1/2
!# 5865 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.75_DP
!# 5897 "space_group.f90"
END SUBROUTINE find_equiv_109
!# 5899 "space_group.f90"
SUBROUTINE find_equiv_110( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5905 "space_group.f90"
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !+x,-y,+z+1/2
         !-x+1/2,+y+1/2,+z
         !-y,-x+1/2,+z+3/4
         !y+1/2,x,z+1/4
!# 5914 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.25_DP
!# 5946 "space_group.f90"
END SUBROUTINE find_equiv_110
!# 5948 "space_group.f90"
SUBROUTINE find_equiv_111( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 5954 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x,+y,-z
         !+x,-y,-z
         !-y,-x,+z
         !y,x,z
!# 5963 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 5995 "space_group.f90"
END SUBROUTINE find_equiv_111
!# 5997 "space_group.f90"
SUBROUTINE find_equiv_112( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6003 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-y,-x,+z+1/2
         !y,x,z+1/2
!# 6012 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 6044 "space_group.f90"
END SUBROUTINE find_equiv_112
!# 6046 "space_group.f90"
SUBROUTINE find_equiv_113( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6052 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
!# 6061 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
!# 6093 "space_group.f90"
END SUBROUTINE find_equiv_113
!# 6095 "space_group.f90"
SUBROUTINE find_equiv_114( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6101 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
!# 6110 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
!# 6142 "space_group.f90"
END SUBROUTINE find_equiv_114
!# 6144 "space_group.f90"
SUBROUTINE find_equiv_115( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6150 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
!# 6159 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 6191 "space_group.f90"
END SUBROUTINE find_equiv_115
!# 6193 "space_group.f90"
SUBROUTINE find_equiv_116( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6199 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
!# 6208 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 6240 "space_group.f90"
END SUBROUTINE find_equiv_116
!# 6242 "space_group.f90"
SUBROUTINE find_equiv_117( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6248 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !+y+1/2,+x+1/2,-z
         !-y+1/2,-x+1/2,-z
!# 6257 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
!# 6289 "space_group.f90"
END SUBROUTINE find_equiv_117
!# 6291 "space_group.f90"
SUBROUTINE find_equiv_118( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6297 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !+y+1/2,+x+1/2,-z+1/2
         !-y+1/2,-x+1/2,-z+1/2
!# 6306 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 6338 "space_group.f90"
END SUBROUTINE find_equiv_118
!# 6340 "space_group.f90"
SUBROUTINE find_equiv_119( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6346 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
!# 6355 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
!# 6387 "space_group.f90"
END SUBROUTINE find_equiv_119
!# 6389 "space_group.f90"
SUBROUTINE find_equiv_120( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6395 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
!# 6404 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
!# 6436 "space_group.f90"
END SUBROUTINE find_equiv_120
!# 6438 "space_group.f90"
SUBROUTINE find_equiv_121( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6444 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
!# 6453 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)
!# 6485 "space_group.f90"
END SUBROUTINE find_equiv_121
!# 6487 "space_group.f90"
SUBROUTINE find_equiv_122( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6493 "space_group.f90"
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y,-z+3/4
         !+x+1/2,-y,-z+3/4
         !-y+1/2,-x,+z+3/4
         !+y+1/2,+x,+z+3/4
!# 6502 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.75_DP
!# 6534 "space_group.f90"
END SUBROUTINE find_equiv_122
!# 6536 "space_group.f90"
SUBROUTINE find_equiv_123( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6542 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 6559 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
!# 6623 "space_group.f90"
END SUBROUTINE find_equiv_123
!# 6625 "space_group.f90"
SUBROUTINE find_equiv_124( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6631 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z+1/2
         !-x,y,z+1/2
         !-y,-x,z+1/2
         !y,x,z+1/2
!# 6648 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
!# 6712 "space_group.f90"
END SUBROUTINE find_equiv_124
!# 6714 "space_group.f90"
SUBROUTINE find_equiv_125( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6721 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x+1/2,-z
         !-y+1/2,x+1/2,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
         !-y+1/2,-x+1/2,z
         !y+1/2,x+1/2,z
!# 6739 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
!# 6804 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y+1/2,+x,+z
         !+y,-x+1/2,+z
         !-x+1/2,+y,-z
         !+x,-y+1/2,-z
         !+y,+x,-z
         !-y+1/2,-x+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x,-z
         !-y,x+1/2,-z
         !x+1/2,-y,z
         !-x,y+1/2,z
         !-y,-x,z
         !y+1/2,x+1/2,z
!# 6822 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
!# 6887 "space_group.f90"
END SUBROUTINE find_equiv_125
!# 6889 "space_group.f90"
SUBROUTINE find_equiv_126( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 6896 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !y+1/2,-x+1/2,-z+1/2
         !-y+1/2,x+1/2,-z+1/2
         !x+1/2,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z+1/2
         !-y+1/2,-x+1/2,z+1/2
         !y+1/2,x+1/2,z+1/2
!# 6914 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 6979 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 6997 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 7062 "space_group.f90"
END SUBROUTINE find_equiv_126
!# 7064 "space_group.f90"
SUBROUTINE find_equiv_127( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7070 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7087 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
!# 7151 "space_group.f90"
END SUBROUTINE find_equiv_127
!# 7153 "space_group.f90"
SUBROUTINE find_equiv_128( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7159 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7176 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
!# 7240 "space_group.f90"
END SUBROUTINE find_equiv_128
!# 7242 "space_group.f90"
SUBROUTINE find_equiv_129( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7249 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,+x+1/2,+z
         !+y+1/2,-x+1/2,+z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y+1/2,-x+1/2,z
         !y+1/2,x+1/2,z
!# 7267 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
!# 7332 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7350 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF
!# 7415 "space_group.f90"
END SUBROUTINE find_equiv_129
!# 7417 "space_group.f90"
SUBROUTINE find_equiv_130( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7424 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7442 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 7507 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7525 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 7590 "space_group.f90"
END SUBROUTINE find_equiv_130
!# 7592 "space_group.f90"
SUBROUTINE find_equiv_131( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7598 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7615 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
!# 7679 "space_group.f90"
END SUBROUTINE find_equiv_131
!# 7681 "space_group.f90"
SUBROUTINE find_equiv_132( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7687 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7704 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
!# 7768 "space_group.f90"
END SUBROUTINE find_equiv_132
!# 7770 "space_group.f90"
SUBROUTINE find_equiv_133( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7777 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7795 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 7860 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7878 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 7943 "space_group.f90"
END SUBROUTINE find_equiv_133
!# 7945 "space_group.f90"
SUBROUTINE find_equiv_134( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 7952 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 7970 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF
!# 8035 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8053 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
!# 8118 "space_group.f90"
END SUBROUTINE find_equiv_134
!# 8120 "space_group.f90"
SUBROUTINE find_equiv_135( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8126 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8143 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
!# 8207 "space_group.f90"
END SUBROUTINE find_equiv_135
!# 8209 "space_group.f90"
SUBROUTINE find_equiv_136( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8215 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8232 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
!# 8296 "space_group.f90"
END SUBROUTINE find_equiv_136
!# 8298 "space_group.f90"
SUBROUTINE find_equiv_137( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8305 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8323 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 8388 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8406 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
!# 8471 "space_group.f90"
END SUBROUTINE find_equiv_137
!# 8473 "space_group.f90"
SUBROUTINE find_equiv_138( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8480 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8498 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
!# 8563 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8581 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF
!# 8646 "space_group.f90"
END SUBROUTINE find_equiv_138
!# 8648 "space_group.f90"
SUBROUTINE find_equiv_139( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8654 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8671 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
!# 8735 "space_group.f90"
END SUBROUTINE find_equiv_139
!# 8737 "space_group.f90"
SUBROUTINE find_equiv_140( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8743 "space_group.f90"
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8760 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
!# 8824 "space_group.f90"
END SUBROUTINE find_equiv_140
!# 8826 "space_group.f90"
SUBROUTINE find_equiv_141( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 8833 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8851 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.25_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.75_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         END IF
!# 8916 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 8934 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.25_DP
         outco(2,3,i)=+inco(1,i)+0.75_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.25_DP
         outco(2,4,i)=-inco(1,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=+inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=-inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.75_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.75_DP
         outco(2,11,i)=-inco(1,i)+0.25_DP
         outco(3,11,i)=-inco(3,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=+inco(1,i)+0.75_DP
         outco(3,12,i)=-inco(3,i)+0.25_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=+inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         END IF
!# 8999 "space_group.f90"
END SUBROUTINE find_equiv_141
!# 9001 "space_group.f90"
SUBROUTINE find_equiv_142( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9008 "space_group.f90"
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 9026 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.25_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.75_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         END IF
!# 9091 "space_group.f90"
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z
!# 9109 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.25_DP
         outco(2,3,i)=+inco(1,i)+0.75_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.25_DP
         outco(2,4,i)=-inco(1,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=+inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=-inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.25_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.75_DP
         outco(2,11,i)=-inco(1,i)+0.25_DP
         outco(3,11,i)=-inco(3,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=+inco(1,i)+0.75_DP
         outco(3,12,i)=-inco(3,i)+0.25_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=+inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         END IF
!# 9174 "space_group.f90"
      !*****************************************
      !Trigonal 143-167
!# 9177 "space_group.f90"
END SUBROUTINE find_equiv_142
!# 9179 "space_group.f90"
SUBROUTINE find_equiv_143( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9185 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
!# 9197 "space_group.f90"
END SUBROUTINE find_equiv_143
!# 9199 "space_group.f90"
SUBROUTINE find_equiv_144( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9205 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
!# 9217 "space_group.f90"
END SUBROUTINE find_equiv_144
!# 9219 "space_group.f90"
SUBROUTINE find_equiv_145( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9225 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
!# 9237 "space_group.f90"
END SUBROUTINE find_equiv_145
!# 9239 "space_group.f90"
SUBROUTINE find_equiv_146( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9246 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=+inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=+inco(1,i)
         END IF
!# 9260 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         END IF
!# 9274 "space_group.f90"
END SUBROUTINE find_equiv_146
!# 9276 "space_group.f90"
SUBROUTINE find_equiv_147( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9282 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=-inco(3,i)
END SUBROUTINE find_equiv_147
!# 9307 "space_group.f90"
SUBROUTINE find_equiv_148( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9314 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=+inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(3,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(2,i)
         !S=6
         outco(1,6,i)=-inco(2,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(1,i)
         END IF
!# 9340 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=-inco(3,i)
         END IF
END SUBROUTINE find_equiv_148
!# 9367 "space_group.f90"
SUBROUTINE find_equiv_149( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9373 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
!# 9397 "space_group.f90"
END SUBROUTINE find_equiv_149
!# 9399 "space_group.f90"
SUBROUTINE find_equiv_150( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9405 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
!# 9429 "space_group.f90"
END SUBROUTINE find_equiv_150
!# 9431 "space_group.f90"
SUBROUTINE find_equiv_151( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9437 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+duterz
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+unterz
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
!# 9461 "space_group.f90"
END SUBROUTINE find_equiv_151
!# 9463 "space_group.f90"
SUBROUTINE find_equiv_152( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9469 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+duterz
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+unterz
!# 9493 "space_group.f90"
END SUBROUTINE find_equiv_152
!# 9495 "space_group.f90"
SUBROUTINE find_equiv_153( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9501 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+unterz
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+duterz
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
!# 9525 "space_group.f90"
END SUBROUTINE find_equiv_153
!# 9527 "space_group.f90"
SUBROUTINE find_equiv_154( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9533 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+unterz
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+duterz
!# 9557 "space_group.f90"
END SUBROUTINE find_equiv_154
!# 9559 "space_group.f90"
SUBROUTINE find_equiv_155( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9566 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(1,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(2,i)
         END IF
!# 9592 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         END IF
!# 9618 "space_group.f90"
END SUBROUTINE find_equiv_155
!# 9620 "space_group.f90"
SUBROUTINE find_equiv_156( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9626 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !s=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)
!# 9650 "space_group.f90"
END SUBROUTINE find_equiv_156
!# 9652 "space_group.f90"
SUBROUTINE find_equiv_157( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9658 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=+inco(3,i)
!# 9682 "space_group.f90"
END SUBROUTINE find_equiv_157
!# 9684 "space_group.f90"
SUBROUTINE find_equiv_158( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9690 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !s=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
!# 9714 "space_group.f90"
END SUBROUTINE find_equiv_158
!# 9716 "space_group.f90"
SUBROUTINE find_equiv_159( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9722 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
!# 9746 "space_group.f90"
END SUBROUTINE find_equiv_159
!# 9748 "space_group.f90"
SUBROUTINE find_equiv_160( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9755 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=inco(3,i)
         outco(2,4,i)=inco(2,i)
         outco(3,4,i)=inco(1,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(3,i)
         outco(3,6,i)=inco(2,i)
         END IF
!# 9781 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)
         END IF
!# 9807 "space_group.f90"
END SUBROUTINE find_equiv_160
!# 9809 "space_group.f90"
SUBROUTINE find_equiv_161( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9816 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=inco(3,i)+0.5_DP
         outco(2,4,i)=inco(2,i)+0.5_DP
         outco(3,4,i)=inco(1,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(2,i)+0.5_DP
         outco(2,5,i)=inco(1,i)+0.5_DP
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)+0.5_DP
         outco(2,6,i)=inco(3,i)+0.5_DP
         outco(3,6,i)=inco(2,i)+0.5_DP
         END IF
!# 9842 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         END IF
!# 9868 "space_group.f90"
END SUBROUTINE find_equiv_161
!# 9870 "space_group.f90"
SUBROUTINE find_equiv_162( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9876 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(2,i)
         outco(2,10,i)=+inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)
!# 9924 "space_group.f90"
END SUBROUTINE find_equiv_162
!# 9926 "space_group.f90"
SUBROUTINE find_equiv_163( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9932 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(2,i)
         outco(2,10,i)=+inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
!# 9980 "space_group.f90"
END SUBROUTINE find_equiv_163
!# 9982 "space_group.f90"
SUBROUTINE find_equiv_164( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 9988 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)
!# 10036 "space_group.f90"
END SUBROUTINE find_equiv_164
!# 10038 "space_group.f90"
SUBROUTINE find_equiv_165( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10044 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
!# 10092 "space_group.f90"
END SUBROUTINE find_equiv_165
!# 10094 "space_group.f90"
SUBROUTINE find_equiv_166( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10101 "space_group.f90"
         IF (unique=='1') THEN
         !Rhombohedral
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(1,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=-inco(2,i)
         outco(2,9,i)=-inco(3,i)
         outco(3,9,i)=-inco(1,i)
         !S=10
         outco(1,10,i)=inco(3,i)
         outco(2,10,i)=inco(2,i)
         outco(3,10,i)=inco(1,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=+inco(1,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(3,i)
         outco(3,12,i)=+inco(2,i)
         END IF
!# 10152 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)
         END IF
!# 10202 "space_group.f90"
END SUBROUTINE find_equiv_166
!# 10204 "space_group.f90"
SUBROUTINE find_equiv_167( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10211 "space_group.f90"
         IF (unique=='1') THEN
         !Rhombohedral
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(1,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)+0.5_DP
         outco(2,5,i)=-inco(1,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(3,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=-inco(2,i)
         outco(2,9,i)=-inco(3,i)
         outco(3,9,i)=-inco(1,i)
         !S=10
         outco(1,10,i)=inco(3,i)+0.5_DP
         outco(2,10,i)=inco(2,i)+0.5_DP
         outco(3,10,i)=inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=+inco(1,i)+0.5_DP
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)+0.5_DP
         outco(2,12,i)=+inco(3,i)+0.5_DP
         outco(3,12,i)=+inco(2,i)+0.5_DP
         END IF
!# 10262 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
         END IF
!# 10312 "space_group.f90"
      !*****************************************
      !Exagonal 168-194
END SUBROUTINE find_equiv_167
!# 10316 "space_group.f90"
SUBROUTINE find_equiv_168( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10322 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
!# 10346 "space_group.f90"
END SUBROUTINE find_equiv_168
!# 10348 "space_group.f90"
SUBROUTINE find_equiv_169( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10354 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+cisest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unsest
!# 10378 "space_group.f90"
END SUBROUTINE find_equiv_169
!# 10380 "space_group.f90"
SUBROUTINE find_equiv_170( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10386 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unsest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+cisest
!# 10410 "space_group.f90"
END SUBROUTINE find_equiv_170
!# 10412 "space_group.f90"
SUBROUTINE find_equiv_171( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10418 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+duterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unterz
!# 10442 "space_group.f90"
END SUBROUTINE find_equiv_171
!# 10444 "space_group.f90"
SUBROUTINE find_equiv_172( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10450 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+duterz
!# 10474 "space_group.f90"
END SUBROUTINE find_equiv_172
!# 10476 "space_group.f90"
SUBROUTINE find_equiv_173( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10482 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
!# 10506 "space_group.f90"
END SUBROUTINE find_equiv_173
!# 10508 "space_group.f90"
SUBROUTINE find_equiv_174( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10514 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)
!# 10538 "space_group.f90"
END SUBROUTINE find_equiv_174
!# 10540 "space_group.f90"
SUBROUTINE find_equiv_175( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10546 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(2,i)
         outco(2,11,i)=+inco(1,i)-inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)+inco(2,i)
         outco(2,12,i)=-inco(1,i)
         outco(3,12,i)=-inco(3,i)
!# 10594 "space_group.f90"
END SUBROUTINE find_equiv_175
!# 10596 "space_group.f90"
SUBROUTINE find_equiv_176( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10602 "space_group.f90"
                  DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(2,i)
         outco(2,11,i)=+inco(1,i)-inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)+inco(2,i)
         outco(2,12,i)=-inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
!# 10650 "space_group.f90"
END SUBROUTINE find_equiv_176
!# 10652 "space_group.f90"
SUBROUTINE find_equiv_177( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10658 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
!# 10706 "space_group.f90"
END SUBROUTINE find_equiv_177
!# 10708 "space_group.f90"
SUBROUTINE find_equiv_178( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10714 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+cisest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unsest
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+unterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+duterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+cisest
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+unsest
!# 10762 "space_group.f90"
END SUBROUTINE find_equiv_178
!# 10764 "space_group.f90"
SUBROUTINE find_equiv_179( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10770 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unsest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+cisest
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+duterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+unterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+unsest
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+cisest
!# 10818 "space_group.f90"
END SUBROUTINE find_equiv_179
!# 10820 "space_group.f90"
SUBROUTINE find_equiv_180( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10826 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+duterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unterz
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+duterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+unterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+duterz
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+unterz
!# 10874 "space_group.f90"
END SUBROUTINE find_equiv_180
!# 10876 "space_group.f90"
SUBROUTINE find_equiv_181( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10882 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+duterz
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+unterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+duterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+unterz
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+duterz
!# 10930 "space_group.f90"
END SUBROUTINE find_equiv_181
!# 10932 "space_group.f90"
SUBROUTINE find_equiv_182( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10938 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
!# 10986 "space_group.f90"
END SUBROUTINE find_equiv_182
!# 10988 "space_group.f90"
SUBROUTINE find_equiv_183( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 10994 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)
!# 11042 "space_group.f90"
END SUBROUTINE find_equiv_183
!# 11044 "space_group.f90"
SUBROUTINE find_equiv_184( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11050 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
!# 11098 "space_group.f90"
END SUBROUTINE find_equiv_184
!# 11100 "space_group.f90"
SUBROUTINE find_equiv_185( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11106 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)
!# 11154 "space_group.f90"
END SUBROUTINE find_equiv_185
!# 11156 "space_group.f90"
SUBROUTINE find_equiv_186( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11162 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
!# 11210 "space_group.f90"
END SUBROUTINE find_equiv_186
!# 11212 "space_group.f90"
SUBROUTINE find_equiv_187( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11218 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
!# 11266 "space_group.f90"
END SUBROUTINE find_equiv_187
!# 11268 "space_group.f90"
SUBROUTINE find_equiv_188( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11274 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
!# 11322 "space_group.f90"
END SUBROUTINE find_equiv_188
!# 11324 "space_group.f90"
SUBROUTINE find_equiv_189( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11330 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=inco(3,i)
!# 11378 "space_group.f90"
END SUBROUTINE find_equiv_189
!# 11380 "space_group.f90"
SUBROUTINE find_equiv_190( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11386 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=inco(3,i)+0.5_DP
!# 11434 "space_group.f90"
END SUBROUTINE find_equiv_190
!# 11436 "space_group.f90"
SUBROUTINE find_equiv_191( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11442 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)
!# 11538 "space_group.f90"
END SUBROUTINE find_equiv_191
!# 11540 "space_group.f90"
SUBROUTINE find_equiv_192( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11546 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)+0.5_DP
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)+0.5_DP
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)+0.5_DP
!# 11642 "space_group.f90"
END SUBROUTINE find_equiv_192
!# 11644 "space_group.f90"
SUBROUTINE find_equiv_193( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11650 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)+0.5_DP
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)
!# 11746 "space_group.f90"
END SUBROUTINE find_equiv_193
!# 11748 "space_group.f90"
SUBROUTINE find_equiv_194( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11754 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)+0.5_DP
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)+0.5_DP
      !*****************************************
      !Cubic 195-230
END SUBROUTINE find_equiv_194
!# 11853 "space_group.f90"
SUBROUTINE find_equiv_195( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11859 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
!# 11907 "space_group.f90"
END SUBROUTINE find_equiv_195
!# 11909 "space_group.f90"
SUBROUTINE find_equiv_196( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11915 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
!# 11963 "space_group.f90"
END SUBROUTINE find_equiv_196
!# 11965 "space_group.f90"
SUBROUTINE find_equiv_197( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 11971 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
!# 12019 "space_group.f90"
END SUBROUTINE find_equiv_197
!# 12021 "space_group.f90"
SUBROUTINE find_equiv_198( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12027 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
!# 12075 "space_group.f90"
END SUBROUTINE find_equiv_198
!# 12077 "space_group.f90"
SUBROUTINE find_equiv_199( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12083 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
!# 12131 "space_group.f90"
END SUBROUTINE find_equiv_199
!# 12133 "space_group.f90"
SUBROUTINE find_equiv_200( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12139 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)
!# 12235 "space_group.f90"
END SUBROUTINE find_equiv_200
!# 12237 "space_group.f90"
SUBROUTINE find_equiv_201( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12244 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)+0.5_DP
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(3,i)+0.5_DP
         outco(2,17,i)=-inco(1,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)+0.5_DP
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)+0.5_DP
         outco(2,21,i)=-inco(3,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)+0.5_DP
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         END IF
!# 12342 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)+0.5_DP
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)+0.5_DP
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)+0.5_DP
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)
         END IF
!# 12440 "space_group.f90"
END SUBROUTINE find_equiv_201
!# 12442 "space_group.f90"
SUBROUTINE find_equiv_202( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12448 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)
!# 12544 "space_group.f90"
END SUBROUTINE find_equiv_202
!# 12546 "space_group.f90"
SUBROUTINE find_equiv_203( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12553 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)+0.25_DP
         outco(2,13,i)=-inco(2,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=+inco(1,i)+0.25_DP
         outco(2,14,i)=+inco(2,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(1,i)+0.25_DP
         outco(2,15,i)=-inco(2,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.25_DP
         outco(2,16,i)=+inco(2,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=-inco(3,i)+0.25_DP
         outco(2,17,i)=-inco(1,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(3,i)+0.25_DP
         outco(2,18,i)=+inco(1,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.25_DP
         outco(2,19,i)=+inco(1,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(3,i)+0.25_DP
         outco(2,20,i)=-inco(1,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=-inco(2,i)+0.25_DP
         outco(2,21,i)=-inco(3,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(2,i)+0.25_DP
         outco(2,22,i)=-inco(3,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.25_DP
         outco(2,23,i)=inco(3,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.25_DP
         outco(2,24,i)=+inco(3,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         END IF
!# 12651 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.75_DP
         outco(3,6,i)=-inco(2,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.75_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.75_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.75_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.75_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.75_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.75_DP
         outco(3,11,i)=-inco(1,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=-inco(3,i)+0.75_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.25_DP
         outco(2,14,i)=+inco(2,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)+0.25_DP
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.25_DP
         outco(2,19,i)=+inco(1,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)+0.25_DP
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)+0.25_DP
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.25_DP
         outco(2,24,i)=+inco(3,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)
         END IF
!# 12749 "space_group.f90"
END SUBROUTINE find_equiv_203
!# 12751 "space_group.f90"
SUBROUTINE find_equiv_204( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12757 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)
!# 12853 "space_group.f90"
END SUBROUTINE find_equiv_204
!# 12855 "space_group.f90"
SUBROUTINE find_equiv_205( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12861 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)+0.5_DP
!# 12957 "space_group.f90"
END SUBROUTINE find_equiv_205
!# 12959 "space_group.f90"
SUBROUTINE find_equiv_206( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 12965 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)+0.5_DP
!# 13061 "space_group.f90"
END SUBROUTINE find_equiv_206
!# 13063 "space_group.f90"
SUBROUTINE find_equiv_207( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13069 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
!# 13165 "space_group.f90"
END SUBROUTINE find_equiv_207
!# 13167 "space_group.f90"
SUBROUTINE find_equiv_208( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13173 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
!# 13269 "space_group.f90"
END SUBROUTINE find_equiv_208
!# 13271 "space_group.f90"
SUBROUTINE find_equiv_209( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13277 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
!# 13373 "space_group.f90"
END SUBROUTINE find_equiv_209
!# 13375 "space_group.f90"
SUBROUTINE find_equiv_210( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13381 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
!# 13477 "space_group.f90"
END SUBROUTINE find_equiv_210
!# 13479 "space_group.f90"
SUBROUTINE find_equiv_211( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13485 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
!# 13581 "space_group.f90"
END SUBROUTINE find_equiv_211
!# 13583 "space_group.f90"
SUBROUTINE find_equiv_212( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13589 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.25_DP
         outco(2,13,i)=inco(1,i)+0.75_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.25_DP
         outco(2,17,i)=+inco(3,i)+0.75_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.75_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.75_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.25_DP
         outco(2,21,i)=inco(2,i)+0.75_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.75_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.75_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
!# 13685 "space_group.f90"
END SUBROUTINE find_equiv_212
!# 13687 "space_group.f90"
SUBROUTINE find_equiv_213( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13693 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP
!# 13789 "space_group.f90"
END SUBROUTINE find_equiv_213
!# 13791 "space_group.f90"
SUBROUTINE find_equiv_214( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13797 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP
!# 13893 "space_group.f90"
END SUBROUTINE find_equiv_214
!# 13895 "space_group.f90"
SUBROUTINE find_equiv_215( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 13901 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)
!# 13997 "space_group.f90"
END SUBROUTINE find_equiv_215
!# 13999 "space_group.f90"
SUBROUTINE find_equiv_216( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14005 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)
!# 14101 "space_group.f90"
END SUBROUTINE find_equiv_216
!# 14103 "space_group.f90"
SUBROUTINE find_equiv_217( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14109 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)
!# 14205 "space_group.f90"
END SUBROUTINE find_equiv_217
!# 14207 "space_group.f90"
SUBROUTINE find_equiv_218( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14213 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=-inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=-inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=-inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=-inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=-inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=inco(1,i)+0.5_DP
!# 14309 "space_group.f90"
END SUBROUTINE find_equiv_218
!# 14311 "space_group.f90"
SUBROUTINE find_equiv_219( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14317 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=-inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=-inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=-inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=-inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=-inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=inco(1,i)+0.5_DP
!# 14413 "space_group.f90"
END SUBROUTINE find_equiv_219
!# 14415 "space_group.f90"
SUBROUTINE find_equiv_220( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14421 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.25_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=-inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=-inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.25_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=-inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.75_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=-inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.25_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.75_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=-inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=-inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=inco(1,i)+0.75_DP
!# 14517 "space_group.f90"
END SUBROUTINE find_equiv_220
!# 14519 "space_group.f90"
SUBROUTINE find_equiv_221( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14525 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
!# 14717 "space_group.f90"
END SUBROUTINE find_equiv_221
!# 14719 "space_group.f90"
SUBROUTINE find_equiv_222( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 14726 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)+0.5_DP
         outco(2,25,i)=-inco(2,i)+0.5_DP
         outco(3,25,i)=-inco(3,i)+0.5_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.5_DP
         outco(2,29,i)=-inco(1,i)+0.5_DP
         outco(3,29,i)=-inco(2,i)+0.5_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.5_DP
         outco(2,33,i)=-inco(3,i)+0.5_DP
         outco(3,33,i)=-inco(1,i)+0.5_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF
!# 14920 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF
!# 15114 "space_group.f90"
END SUBROUTINE find_equiv_222
!# 15116 "space_group.f90"
SUBROUTINE find_equiv_223( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 15122 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
!# 15314 "space_group.f90"
END SUBROUTINE find_equiv_223
!# 15316 "space_group.f90"
SUBROUTINE find_equiv_224( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 15323 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.5_DP
         outco(2,25,i)=-inco(2,i)+0.5_DP
         outco(3,25,i)=-inco(3,i)+0.5_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.5_DP
         outco(2,29,i)=-inco(1,i)+0.5_DP
         outco(3,29,i)=-inco(2,i)+0.5_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.5_DP
         outco(2,33,i)=-inco(3,i)+0.5_DP
         outco(3,33,i)=-inco(1,i)+0.5_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF
!# 15517 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF
!# 15711 "space_group.f90"
END SUBROUTINE find_equiv_224
!# 15713 "space_group.f90"
SUBROUTINE find_equiv_225( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 15719 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
!# 15911 "space_group.f90"
END SUBROUTINE find_equiv_225
!# 15913 "space_group.f90"
SUBROUTINE find_equiv_226( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 15919 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
!# 16111 "space_group.f90"
END SUBROUTINE find_equiv_226
!# 16113 "space_group.f90"
SUBROUTINE find_equiv_227( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 16120 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.25_DP
         outco(2,25,i)=-inco(2,i)+0.25_DP
         outco(3,25,i)=-inco(3,i)+0.25_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.25_DP
         outco(2,26,i)=inco(2,i)+0.75_DP
         outco(3,26,i)=-inco(3,i)+0.75_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.75_DP
         outco(2,27,i)=-inco(2,i)+0.75_DP
         outco(3,27,i)=inco(3,i)+0.25_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.75_DP
         outco(2,28,i)=inco(2,i)+0.25_DP
         outco(3,28,i)=inco(3,i)+0.75_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.25_DP
         outco(2,29,i)=-inco(1,i)+0.25_DP
         outco(3,29,i)=-inco(2,i)+0.25_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.75_DP
         outco(2,30,i)=inco(1,i)+0.25_DP
         outco(3,30,i)=inco(2,i)+0.75_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.25_DP
         outco(2,31,i)=inco(1,i)+0.75_DP
         outco(3,31,i)=-inco(2,i)+0.75_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.75_DP
         outco(2,32,i)=-inco(1,i)+0.75_DP
         outco(3,32,i)=inco(2,i)+0.25_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.25_DP
         outco(2,33,i)=-inco(3,i)+0.25_DP
         outco(3,33,i)=-inco(1,i)+0.25_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.75_DP
         outco(2,34,i)=-inco(3,i)+0.75_DP
         outco(3,34,i)=inco(1,i)+0.25_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.75_DP
         outco(2,35,i)=inco(3,i)+0.25_DP
         outco(3,35,i)=inco(1,i)+0.75_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.25_DP
         outco(2,36,i)=inco(3,i)+0.75_DP
         outco(3,36,i)=-inco(1,i)+0.75_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF
!# 16314 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.25_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.25_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.75_DP
         outco(3,6,i)=-inco(2,i)+0.25_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.75_DP
         outco(2,7,i)=-inco(1,i)+0.25_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.25_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.75_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.25_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.75_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.75_DP
         outco(3,11,i)=-inco(1,i)+0.25_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=-inco(3,i)+0.25_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.25_DP
         outco(2,26,i)=inco(2,i)+0.75_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.75_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.25_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.25_DP
         outco(3,28,i)=inco(3,i)+0.75_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.25_DP
         outco(3,30,i)=inco(2,i)+0.75_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.25_DP
         outco(2,31,i)=inco(1,i)+0.75_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.75_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.25_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.75_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.25_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.25_DP
         outco(3,35,i)=inco(1,i)+0.75_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.25_DP
         outco(2,36,i)=inco(3,i)+0.75_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF
!# 16508 "space_group.f90"
END SUBROUTINE find_equiv_227
!# 16510 "space_group.f90"
SUBROUTINE find_equiv_228( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 16517 "space_group.f90"
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.75_DP
         outco(2,25,i)=-inco(2,i)+0.75_DP
         outco(3,25,i)=-inco(3,i)+0.75_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.75_DP
         outco(2,26,i)=inco(2,i)+0.25_DP
         outco(3,26,i)=-inco(3,i)+0.25_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.25_DP
         outco(2,27,i)=-inco(2,i)+0.25_DP
         outco(3,27,i)=inco(3,i)+0.75_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.25_DP
         outco(2,28,i)=inco(2,i)+0.75_DP
         outco(3,28,i)=inco(3,i)+0.25_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.75_DP
         outco(2,29,i)=-inco(1,i)+0.75_DP
         outco(3,29,i)=-inco(2,i)+0.75_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.25_DP
         outco(2,30,i)=inco(1,i)+0.75_DP
         outco(3,30,i)=inco(2,i)+0.25_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.75_DP
         outco(2,31,i)=inco(1,i)+0.25_DP
         outco(3,31,i)=-inco(2,i)+0.25_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.25_DP
         outco(2,32,i)=-inco(1,i)+0.25_DP
         outco(3,32,i)=inco(2,i)+0.75_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.75_DP
         outco(2,33,i)=-inco(3,i)+0.75_DP
         outco(3,33,i)=-inco(1,i)+0.75_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.25_DP
         outco(2,34,i)=-inco(3,i)+0.25_DP
         outco(3,34,i)=inco(1,i)+0.75_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.25_DP
         outco(2,35,i)=inco(3,i)+0.75_DP
         outco(3,35,i)=inco(1,i)+0.25_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.75_DP
         outco(2,36,i)=inco(3,i)+0.25_DP
         outco(3,36,i)=-inco(1,i)+0.25_DP
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF
!# 16711 "space_group.f90"
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.25_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.25_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.25_DP
         outco(3,6,i)=-inco(2,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.25_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.75_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.25_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.75_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.25_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.25_DP
         outco(3,11,i)=-inco(1,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.25_DP
         outco(2,12,i)=-inco(3,i)+0.75_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.75_DP
         outco(2,26,i)=inco(2,i)+0.25_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.25_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.75_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.75_DP
         outco(3,28,i)=inco(3,i)+0.25_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.75_DP
         outco(3,30,i)=inco(2,i)+0.25_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.75_DP
         outco(2,31,i)=inco(1,i)+0.25_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.25_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.75_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.25_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.75_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.75_DP
         outco(3,35,i)=inco(1,i)+0.25_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.75_DP
         outco(2,36,i)=inco(3,i)+0.25_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF
!# 16905 "space_group.f90"
END SUBROUTINE find_equiv_228
!# 16907 "space_group.f90"
SUBROUTINE find_equiv_229( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 16913 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
!# 17105 "space_group.f90"
END SUBROUTINE find_equiv_229
!# 17107 "space_group.f90"
SUBROUTINE find_equiv_230( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k
!# 17113 "space_group.f90"
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)+0.75_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.25_DP
         outco(2,38,i)=inco(1,i)+0.25_DP
         outco(3,38,i)=inco(3,i)+0.25_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)+0.75_DP
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.75_DP
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)+0.75_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.75_DP
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.25_DP
         outco(2,43,i)=inco(3,i)+0.25_DP
         outco(3,43,i)=inco(2,i)+0.25_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)+0.75_DP
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)+0.75_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)+0.75_DP
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.75_DP
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.25_DP
         outco(2,48,i)=inco(2,i)+0.25_DP
         outco(3,48,i)=inco(1,i)+0.25_DP
       END SUBROUTINE FIND_EQUIV_230
!# 17306 "space_group.f90"
     END MODULE space_group
