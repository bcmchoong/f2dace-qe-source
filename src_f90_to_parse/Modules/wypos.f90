!# 1 "wypos.f90"
!
! Copyright (C) 2014 Federico Zadra
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wy_pos
!# 10 "wypos.f90"
!! Main subroutine: \(\texttt{wypos}\). Converts atomic positions
!! given in Wyckoff convention: multiplicity-letter + parameter(s),
!! to crystal positions.
!# 14 "wypos.f90"
USE kinds,  ONLY : DP
IMPLICIT NONE
!# 17 "wypos.f90"
SAVE
PRIVATE
!# 20 "wypos.f90"
PUBLIC wypos
!# 22 "wypos.f90"
CONTAINS
   SUBROUTINE wypos(tau,wp,inp,space_group_number,uniqueb,&
                                      rhombohedral,origin_choice)
   !-----------------------------------------------------------
   !! Convert atomic positions given in Wyckoff convention:
   !! multiplicity-letter + parameter(s), to crystal positions.
   !-----------------------------------------------------------
!# 30 "wypos.f90"
      REAL(DP), DIMENSION(3), INTENT(OUT) :: tau
      REAL(DP), INTENT(IN) :: inp(3)
      !! parameter(s) (if needed)
      CHARACTER(LEN=*), INTENT (IN) :: wp
      !! Wyckoff label (e.g. 8c)
      INTEGER, INTENT(IN) :: space_group_number
      LOGICAL, INTENT(IN) :: uniqueb, rhombohedral
      INTEGER, INTENT(IN) :: origin_choice
!# 39 "wypos.f90"
      tau=1.d5
!# 41 "wypos.f90"
      SELECT CASE (space_group_number)
         CASE (1, 4, 7, 9, 19, 29, 33, 76, 78, 144, 145, 169, 170)
            ! groups having only one set of positions "Na x y z": do nothing
         CASE (2) !P-1
             CALL wypos_2  ( wp, tau )
         CASE (3) !P2
             CALL wypos_3  ( wp, inp, uniqueb, tau )
         CASE (5) !C2
             CALL wypos_5  ( wp, inp, uniqueb, tau )
         CASE (6) !Pm
             CALL wypos_6  ( wp, inp, uniqueb, tau )
         CASE (8) !Cm
             CALL wypos_8  ( wp, inp, uniqueb, tau )
         CASE (10) !P2/m
             CALL wypos_10 ( wp, inp, uniqueb, tau )
         CASE (11) !P2(1)/m
             CALL wypos_11 ( wp, inp, uniqueb, tau )
         CASE (12) !C2/m
             CALL wypos_12 ( wp, inp, uniqueb, tau )
         CASE (13) !P2/c
             CALL wypos_13 ( wp, inp, uniqueb, tau )
         CASE (14) !-P2(1)/c
             CALL wypos_14 ( wp, inp, uniqueb, tau )
         CASE (15) !C2/c
             CALL wypos_15 ( wp, inp, uniqueb, tau )
         CASE (16) !P222
             CALL wypos_16 ( wp, inp, tau )
         CASE (17) !P222(1)
             CALL wypos_17 ( wp, inp, tau )
         CASE (18) !P2(1)2(1)2
             CALL wypos_18 ( wp, inp, tau )
         CASE (20) !C222(1)
             CALL wypos_20 ( wp, inp, tau )
         CASE (21) !C222
             CALL wypos_21 ( wp, inp, tau )
         CASE (22) !F222
             CALL wypos_22 ( wp, inp, tau )
         CASE (23) !I222
             CALL wypos_23 ( wp, inp, tau )
         CASE (24) !I2(1)2(1)2(1)
             CALL wypos_24 ( wp, inp, tau )
         CASE (25) !Pmm2
             CALL wypos_25 ( wp, inp, tau )
         CASE (26) !Pmc2(1)
             CALL wypos_26 ( wp, inp, tau )
         CASE (27) !Pcc2
             CALL wypos_27 ( wp, inp, tau )
         CASE (28) !Pma2
             CALL wypos_28 ( wp, inp, tau )
         CASE (30) !Pca2(1)
             CALL wypos_30 ( wp, inp, tau )
         CASE (31) !Pmn2(1)
             CALL wypos_31 ( wp, inp, tau )
         CASE (32) !Pba2
             CALL wypos_32 ( wp, inp, tau )
         CASE (34) !Pnn2
             CALL wypos_34 ( wp, inp, tau )
         CASE (35) !Cmm2
             CALL wypos_35 ( wp, inp, tau )
         CASE (36) !Cmc2(1)
             CALL wypos_36 ( wp, inp, tau )
         CASE (37) !Ccc2
             CALL wypos_37 ( wp, inp, tau )
         CASE (38) !Amm2
             CALL wypos_38 ( wp, inp, tau )
         CASE (39) !Aem2
             CALL wypos_39 ( wp, inp, tau )
         CASE (40) !Ama2
             CALL wypos_40 ( wp, inp, tau )
         CASE (41) !Aea2
             CALL wypos_41 ( wp, inp, tau )
         CASE (42) !Fmm2
             CALL wypos_42 ( wp, inp, tau )
         CASE (43) !Fdd2
             CALL wypos_43 ( wp, inp, tau )
         CASE (44) !Imm2
             CALL wypos_44 ( wp, inp, tau )
         CASE (45) !Iba2
             CALL wypos_45 ( wp, inp, tau )
         CASE (46) !Ima2
             CALL wypos_46 ( wp, inp, tau )
         CASE (47) !Pmmm
             CALL wypos_47 ( wp, inp, tau )
         CASE (48) !Pnnn
             CALL wypos_48 ( wp, inp, origin_choice, tau )
         CASE (49) !Pccm
             CALL wypos_49 ( wp, inp, tau )
         CASE (50) !Pban
             CALL wypos_50 ( wp, inp, origin_choice, tau )
         CASE (51) !Pmma
             CALL wypos_51 ( wp, inp, tau )
         CASE (52) !Pnna
             CALL wypos_52 ( wp, inp, tau )
         CASE (53) !Pmna
             CALL wypos_53 ( wp, inp, tau )
         CASE (54) !Pcca
             CALL wypos_54 ( wp, inp, tau )
         CASE (55) !Pbam
             CALL wypos_55 ( wp, inp, tau )
         CASE (56) !Pccn
             CALL wypos_56 ( wp, inp, tau )
         CASE (57) !Pbcm
             CALL wypos_57 ( wp, inp, tau )
         CASE (58) !Pnnm
             CALL wypos_58 ( wp, inp, tau )
         CASE (59) !Pmmn
             CALL wypos_59 ( wp, inp, origin_choice, tau )
         CASE (60) !Pbcn
             CALL wypos_60 ( wp, inp, tau )
         CASE (61) !Pbca
             CALL wypos_61 ( wp, inp, tau )
         CASE (62) !Pnma
             CALL wypos_62 ( wp, inp, tau )
         CASE (63) !Cmcm
             CALL wypos_63 ( wp, inp, tau )
         CASE (64) !Cmce
             CALL wypos_64 ( wp, inp, tau )
         CASE (65) !Cmmm
             CALL wypos_65 ( wp, inp, tau )
         CASE (66) !Cccm
             CALL wypos_66 ( wp, inp, tau )
         CASE (67) !Cmma
             CALL wypos_67 ( wp, inp, tau )
         CASE (68) !Ccce
             CALL wypos_68 ( wp, inp, origin_choice, tau )
         CASE (69) !Fmmm
             CALL wypos_69 ( wp, inp, tau )
         CASE (70) !Fddd
             CALL wypos_70 ( wp, inp, origin_choice, tau )
         CASE (71) !Immm
             CALL wypos_71 ( wp, inp, tau )
         CASE (72) !Ibam
             CALL wypos_72 ( wp, inp, tau )
         CASE (73) !Ibca
             CALL wypos_73 ( wp, inp, tau )
         CASE (74) !Imma
             CALL wypos_74 ( wp, inp, tau )
         CASE (75) !P4
             CALL wypos_75 ( wp, inp, tau )
         CASE (77) !P4(2)
             CALL wypos_77 ( wp, inp, tau )
         CASE (79) !I4(2)
             CALL wypos_79 ( wp, inp, tau )
         CASE (80) !I4(1)
             CALL wypos_80 ( wp, inp, tau )
         CASE (81) !P-4
             CALL wypos_81 ( wp, inp, tau )
         CASE (82) !I-4
             CALL wypos_82 ( wp, inp, tau )
         CASE (83) !P4/m
             CALL wypos_83 ( wp, inp, tau )
         CASE (84)
             CALL wypos_84 ( wp, inp, tau )
         CASE (85)
             CALL wypos_85 ( wp, inp, origin_choice, tau )
         CASE (86)
             CALL wypos_86 ( wp, inp, origin_choice, tau )
         CASE (87) !I4/m
             CALL wypos_87 ( wp, inp, tau )
         CASE (88) !I4(1)/a
             CALL wypos_88 ( wp, inp, origin_choice, tau )
         CASE (89) !P422
             CALL wypos_89 ( wp, inp, tau )
         CASE (90) !P42(1)2
             CALL wypos_90 ( wp, inp, tau )
         CASE (91) !P4(1)22
             CALL wypos_91 ( wp, inp, tau )
         CASE (92) !P4(1)2(1)2
             CALL wypos_92 ( wp, inp, tau )
         CASE (93) !P4(2)22
             CALL wypos_93 ( wp, inp, tau )
         CASE (94) !P4(2)2(1)2
             CALL wypos_94 ( wp, inp, tau )
         CASE (95) !P4(3)22
             CALL wypos_95 ( wp, inp, tau )
         CASE (96) !P4(2)2(1)2
             CALL wypos_96 ( wp, inp, tau )
         CASE (97) !I422
             CALL wypos_97 ( wp, inp, tau )
         CASE (98) !I4(1)22
             CALL wypos_98 ( wp, inp, tau )
         CASE (99) !P4mm
             CALL wypos_99 ( wp, inp, tau )
         CASE (100) !P4bm
             CALL wypos_100( wp, inp, tau )
         CASE (101) !P4(2)cm
             CALL wypos_101( wp, inp, tau )
         CASE (102) !P4(2)nm
             CALL wypos_102( wp, inp, tau )
         CASE (103) !P4cc
             CALL wypos_103( wp, inp, tau )
         CASE (104) !P4nc
             CALL wypos_104( wp, inp, tau )
         CASE (105) !P4(2)mc
             CALL wypos_105( wp, inp, tau )
         CASE (106) !P4(2)bc
             CALL wypos_106( wp, inp, tau )
         CASE (107) !I4mm
             CALL wypos_107( wp, inp, tau )
         CASE (108) !I4cm
             CALL wypos_108( wp, inp, tau )
         CASE (109) !I4(1)md
             CALL wypos_109( wp, inp, tau )
         CASE (110) !I4(1)cd
             CALL wypos_110( wp, inp, tau )
         CASE (111) !P-42m
             CALL wypos_111( wp, inp, tau )
         CASE (112) !P-42c
             CALL wypos_112( wp, inp, tau )
         CASE (113) !P-42(1)m
             CALL wypos_113( wp, inp, tau )
         CASE (114) !P-42(1)c
             CALL wypos_114( wp, inp, tau )
         CASE (115) !P-4m2
             CALL wypos_115( wp, inp, tau )
         CASE (116) !P4c2
             CALL wypos_116( wp, inp, tau )
         CASE (117) !P-4b2
             CALL wypos_117( wp, inp, tau )
         CASE (118) !P-4n2
             CALL wypos_118( wp, inp, tau )
         CASE (119) !I-4m2
             CALL wypos_119( wp, inp, tau )
         CASE (120) !I-4c2
             CALL wypos_120( wp, inp, tau )
         CASE (121) !I-42m
             CALL wypos_121( wp, inp, tau )
         CASE (122) !I-42d
             CALL wypos_122( wp, inp, tau )
         CASE (123) !P4/mmm
             CALL wypos_123( wp, inp, tau )
         CASE (124) !P4/mmc
             CALL wypos_124( wp, inp, tau )
         CASE (125) !P/nbm
             CALL wypos_125( wp, inp, origin_choice, tau )
         CASE (126)
             CALL wypos_126( wp, inp, origin_choice, tau )
         CASE (127) !P4/mbm
             CALL wypos_127( wp, inp, tau )
         CASE (128) !P4/mnc
             CALL wypos_128( wp, inp, tau )
         CASE (129) !P4/nmm
             CALL wypos_129( wp, inp, origin_choice, tau )
         CASE (130) !P4/ncc
             CALL wypos_130( wp, inp, origin_choice, tau )
         CASE (131) !P4(2)/mmc
             CALL wypos_131( wp, inp, tau )
         CASE (132) !P4(2)mcm
             CALL wypos_132( wp, inp, tau )
         CASE (133) !P4(2)/nbc
             CALL wypos_133( wp, inp, origin_choice, tau )
         CASE (134) !P4(2)/nnm
             CALL wypos_134( wp, inp, origin_choice, tau )
         CASE (135) !P3(2)/mbc
             CALL wypos_135( wp, inp, tau )
         CASE (136) !P4(2)/mnm
             CALL wypos_136( wp, inp, tau )
         CASE (137) !P4(2)/nmc
             CALL wypos_137( wp, inp, origin_choice, tau )
         CASE (138) !P4(2)/ncm
             CALL wypos_138( wp, inp, origin_choice, tau )
         CASE (139) !I4/mmm
             CALL wypos_139( wp, inp, tau )
         CASE (140) !I4/mcm
             CALL wypos_140( wp, inp, tau )
         CASE (141) !I4(1)/amd
             CALL wypos_141( wp, inp, origin_choice, tau )
         CASE (142) !I4(1)/acd
             CALL wypos_142( wp, inp, origin_choice, tau )
         CASE (143) !P3
             CALL wypos_143( wp, inp, tau )
         CASE (146) !R3
             CALL wypos_146( wp, inp, rhombohedral, tau )
         CASE (147) !P-3
             CALL wypos_147( wp, inp, tau )
         CASE (148) !R-3
             CALL wypos_148( wp, inp, rhombohedral, tau )
         CASE (149) !P312
             CALL wypos_149( wp, inp, tau )
         CASE (150) !P321
             CALL wypos_150( wp, inp, tau )
         CASE (151) !P3(1)12
             CALL wypos_151( wp, inp, tau )
         CASE (152) !P3(1)21
             CALL wypos_152( wp, inp, tau )
         CASE (153) !P3(2)12
             CALL wypos_153( wp, inp, tau )
         CASE (154) !3(2)21
             CALL wypos_154( wp, inp, tau )
         CASE (155) !R32
             CALL wypos_155( wp, inp, rhombohedral, tau )
         CASE (156) !P-3m1
             CALL wypos_156( wp, inp, tau )
         CASE (157) !P31m
             CALL wypos_157( wp, inp, tau )
         CASE (158) !P3c1
             CALL wypos_158( wp, inp, tau )
         CASE (159) !P31c
             CALL wypos_159( wp, inp, tau )
         CASE (160) !R3m
             CALL wypos_160( wp, inp, rhombohedral, tau )
         CASE (161) !R3c
             CALL wypos_161( wp, inp, rhombohedral, tau )
         CASE (162) !P-31m
             CALL wypos_162( wp, inp, tau )
         CASE (163) !P-31c
             CALL wypos_163( wp, inp, tau )
         CASE (164) !P-3m1
             CALL wypos_164( wp, inp, tau )
         CASE (165) !P-3c1
             CALL wypos_165( wp, inp, tau )
         CASE (166) !R-3m
             CALL wypos_166( wp, inp, rhombohedral, tau )
         CASE (167) !R-3c
             CALL wypos_167( wp, inp, rhombohedral, tau )
         CASE (168) !P6
             CALL wypos_168( wp, inp, tau )
         CASE (171) !P6/m
             CALL wypos_171( wp, inp, tau )
         CASE (172) !P6(4)
             CALL wypos_172( wp, inp, tau )
         CASE (173) !P6(3)
             CALL wypos_173( wp, inp, tau )
         CASE (174) !P-6
             CALL wypos_174( wp, inp, tau )
         CASE (175) !P6/m
             CALL wypos_175( wp, inp, tau )
         CASE (176) !P6(3)/m
             CALL wypos_176( wp, inp, tau )
         CASE (177) !P622
             CALL wypos_177( wp, inp, tau )
         CASE (178) !P6(1)22
             CALL wypos_178( wp, inp, tau )
         CASE (179) !P6(5)22
             CALL wypos_179( wp, inp, tau )
         CASE (180) !P6(2)22
             CALL wypos_180( wp, inp, tau )
         CASE (181) !P6(4)22
             CALL wypos_181( wp, inp, tau )
         CASE (182) !P6(3)22
             CALL wypos_182( wp, inp, tau )
         CASE (183) !P6mm
             CALL wypos_183( wp, inp, tau )
         CASE (184) !P6cc
             CALL wypos_184( wp, inp, tau )
         CASE (185) !P6(3)cm
             CALL wypos_185( wp, inp, tau )
         CASE (186) !P6(3)mc
             CALL wypos_186( wp, inp, tau )
         CASE (187) !P-6m2
             CALL wypos_187( wp, inp, tau )
         CASE (188) !P-6c2
             CALL wypos_188( wp, inp, tau )
         CASE (189) !P-62m
             CALL wypos_189( wp, inp, tau )
         CASE (190) !P-62c
             CALL wypos_190( wp, inp, tau )
         CASE (191) !P6/mmm
             CALL wypos_191( wp, inp, tau )
         CASE (192) !P6/mcc
             CALL wypos_192( wp, inp, tau )
         CASE (193) !P6(3)/mcm
             CALL wypos_193( wp, inp, tau )
         CASE (194) !P6(3)mmc
             CALL wypos_194( wp, inp, tau )
         CASE (195) !P23
             CALL wypos_195( wp, inp, tau )
         CASE (196) !F23
             CALL wypos_196( wp, inp, tau )
         CASE (197) !I23
             CALL wypos_197( wp, inp, tau )
         CASE (198) !P2(1)3
             CALL wypos_198( wp, inp, tau )
         CASE (199) !I2(1)3
             CALL wypos_199( wp, inp, tau )
         CASE (200) !Pm-3
             CALL wypos_200( wp, inp, tau )
         CASE (201) !Pn-3
             CALL wypos_201( wp, inp, origin_choice, tau )
         CASE (202) !Fm-3
             CALL wypos_202( wp, inp, tau )
         CASE (203) !Fd-3
             CALL wypos_203( wp, inp, origin_choice, tau )
         CASE (204) ! Im-3
             CALL wypos_204( wp, inp, tau )
         CASE (205) !Pa-3
             CALL wypos_205( wp, inp, tau )
         CASE (206) !Ia-3
             CALL wypos_206( wp, inp, tau )
         CASE (207) !P432
             CALL wypos_207( wp, inp, tau )
         CASE (208) !P4(2)32
             CALL wypos_208( wp, inp, tau )
         CASE (209) !F432
             CALL wypos_209( wp, inp, tau )
         CASE (210) !F4(1)32
             CALL wypos_210( wp, inp, tau )
         CASE (211) !I432
             CALL wypos_211( wp, inp, tau )
         CASE (212) !P4(3)32
             CALL wypos_212( wp, inp, tau )
         CASE (213) !P4(1)32
             CALL wypos_213( wp, inp, tau )
         CASE (214) !I4(I)32
             CALL wypos_214( wp, inp, tau )
         CASE (215) !P-43m
             CALL wypos_215( wp, inp, tau )
         CASE (216) !F-43m
             CALL wypos_216( wp, inp, tau )
         CASE (217) !I-43m
             CALL wypos_217( wp, inp, tau )
         CASE (218) !P-43n
             CALL wypos_218( wp, inp, tau )
         CASE (219) !F-43c
             CALL wypos_219( wp, inp, tau )
         CASE (220) !I-43d
             CALL wypos_220( wp, inp, tau )
         CASE (221) !Pm-3m
             CALL wypos_221( wp, inp, tau )
         CASE (222) !Pn-3n
             CALL wypos_222( wp, inp, origin_choice, tau )
         CASE (223) !Pm-3n
             CALL wypos_223( wp, inp, tau )
         CASE (224) !Pn-3m
             CALL wypos_224( wp, inp, origin_choice, tau )
         CASE (225) !Fm-3m
             CALL wypos_225( wp, inp, tau )
         CASE (226) !Fm-3c
             CALL wypos_226( wp, inp, tau )
         CASE (227) !Fd-3m
             CALL wypos_227( wp, inp, origin_choice, tau )
         CASE (228) !Fd-3c
             CALL wypos_228( wp, inp, origin_choice, tau )
         CASE (229) !Im-3m
             CALL wypos_229( wp, inp, tau )
         CASE (230) !Ia-3d
             CALL wypos_230( wp, inp, tau )
          CASE DEFAULT
            CALL errore('wypos','group not recognized',1)
          END SELECT
!# 482 "wypos.f90"
         IF (tau(1)==1.d5.OR.tau(2)==1.d5.OR.tau(3)==1.d5) THEN
            IF (inp(1)==1.d5.OR.inp(2)==1.d5.OR.inp(3)==1.d5) THEN
               CALL errore('wypos','wyckoff position not found',1)
            ELSE
               CALL infomsg('wypos','wyckoff position not found, assuming x y z')
               tau(:)=inp(:)
            END IF
         END IF
!# 491 "wypos.f90"
      END SUBROUTINE wypos
!# 493 "wypos.f90"
SUBROUTINE wypos_2  ( wp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF
!# 531 "wypos.f90"
          END SUBROUTINE wypos_2
!# 533 "wypos.f90"
SUBROUTINE wypos_3  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
!# 539 "wypos.f90"
   IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF
!# 558 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 578 "wypos.f90"
          END SUBROUTINE wypos_3
!# 580 "wypos.f90"
SUBROUTINE wypos_5  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF
!# 597 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 609 "wypos.f90"
          END SUBROUTINE wypos_5
!# 611 "wypos.f90"
SUBROUTINE wypos_6  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF
!# 628 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF
!# 640 "wypos.f90"
          END SUBROUTINE wypos_6
!# 642 "wypos.f90"
SUBROUTINE wypos_8  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
!# 655 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF
            ENDIF
!# 663 "wypos.f90"
          END SUBROUTINE wypos_8
!# 665 "wypos.f90"
SUBROUTINE wypos_10 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF
!# 730 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF
!# 790 "wypos.f90"
          END SUBROUTINE wypos_10
!# 792 "wypos.f90"
SUBROUTINE wypos_11 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)
               ENDIF
!# 821 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 845 "wypos.f90"
          END SUBROUTINE wypos_11
!# 847 "wypos.f90"
SUBROUTINE wypos_12 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
!# 892 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
!# 894 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF
            ENDIF
!# 933 "wypos.f90"
          END SUBROUTINE wypos_12
!# 935 "wypos.f90"
SUBROUTINE wypos_13 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
!# 943 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
!# 969 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
!# 971 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
!# 997 "wypos.f90"
            ENDIF
!# 999 "wypos.f90"
          END SUBROUTINE wypos_13
!# 1001 "wypos.f90"
SUBROUTINE wypos_14 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
!# 1009 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
!# 1027 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
!# 1029 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF
!# 1047 "wypos.f90"
            ENDIF
!# 1049 "wypos.f90"
          END SUBROUTINE wypos_14
!# 1051 "wypos.f90"
SUBROUTINE wypos_15 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
!# 1059 "wypos.f90"
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25
               ENDIF
!# 1081 "wypos.f90"
            ELSEIF (.NOT.uniqueb) THEN
!# 1083 "wypos.f90"
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
!# 1105 "wypos.f90"
            ENDIF
!# 1107 "wypos.f90"
          END SUBROUTINE wypos_15
!# 1109 "wypos.f90"
SUBROUTINE wypos_16 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1197 "wypos.f90"
          END SUBROUTINE wypos_16
!# 1199 "wypos.f90"
SUBROUTINE wypos_17 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 1222 "wypos.f90"
          END SUBROUTINE wypos_17
!# 1224 "wypos.f90"
SUBROUTINE wypos_18 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1239 "wypos.f90"
          END SUBROUTINE wypos_18
!# 1241 "wypos.f90"
SUBROUTINE wypos_20 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 1256 "wypos.f90"
          END SUBROUTINE wypos_20
!# 1258 "wypos.f90"
SUBROUTINE wypos_21 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
!# 1263 "wypos.f90"
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
!# 1309 "wypos.f90"
          END SUBROUTINE wypos_21
!# 1311 "wypos.f90"
SUBROUTINE wypos_22 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
!# 1358 "wypos.f90"
          END SUBROUTINE wypos_22
!# 1360 "wypos.f90"
SUBROUTINE wypos_23 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1407 "wypos.f90"
          END SUBROUTINE wypos_23
!# 1409 "wypos.f90"
SUBROUTINE wypos_24 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
!# 1428 "wypos.f90"
          END SUBROUTINE wypos_24
!# 1430 "wypos.f90"
SUBROUTINE wypos_25 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1469 "wypos.f90"
          END SUBROUTINE wypos_25
!# 1471 "wypos.f90"
SUBROUTINE wypos_26 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1486 "wypos.f90"
          END SUBROUTINE wypos_26
!# 1488 "wypos.f90"
SUBROUTINE wypos_27 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1511 "wypos.f90"
          END SUBROUTINE wypos_27
!# 1513 "wypos.f90"
SUBROUTINE wypos_28 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1532 "wypos.f90"
          END SUBROUTINE wypos_28
!# 1534 "wypos.f90"
SUBROUTINE wypos_30 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 1549 "wypos.f90"
          END SUBROUTINE wypos_30
!# 1551 "wypos.f90"
SUBROUTINE wypos_31 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
          END SUBROUTINE wypos_31
!# 1563 "wypos.f90"
SUBROUTINE wypos_32 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1578 "wypos.f90"
          END SUBROUTINE wypos_32
!# 1580 "wypos.f90"
SUBROUTINE wypos_34 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1595 "wypos.f90"
          END SUBROUTINE wypos_34
!# 1597 "wypos.f90"
SUBROUTINE wypos_35 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1624 "wypos.f90"
          END SUBROUTINE wypos_35
!# 1626 "wypos.f90"
SUBROUTINE wypos_36 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1637 "wypos.f90"
          END SUBROUTINE wypos_36
!# 1639 "wypos.f90"
SUBROUTINE wypos_37 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
!# 1658 "wypos.f90"
END SUBROUTINE wypos_37 
!# 1660 "wypos.f90"
SUBROUTINE wypos_38 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1687 "wypos.f90"
END SUBROUTINE wypos_38 
!# 1689 "wypos.f90"
SUBROUTINE wypos_39 ( wp, inp, tau )
  CHARACTER(LEN=*), INTENT(in)  :: wp
  REAL(dp), INTENT(in) :: inp(3)
  REAL(dp), INTENT(out) :: tau (3)
  
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
!# 1708 "wypos.f90"
END SUBROUTINE wypos_39 
!# 1710 "wypos.f90"
SUBROUTINE wypos_40 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1725 "wypos.f90"
END SUBROUTINE wypos_40 
!# 1727 "wypos.f90"
SUBROUTINE wypos_41 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 1738 "wypos.f90"
END SUBROUTINE wypos_41 
!# 1740 "wypos.f90"
SUBROUTINE wypos_42 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 1763 "wypos.f90"
END SUBROUTINE wypos_42 
!# 1765 "wypos.f90"
SUBROUTINE wypos_43 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 1776 "wypos.f90"
END SUBROUTINE wypos_43 
!# 1778 "wypos.f90"
SUBROUTINE wypos_44 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1801 "wypos.f90"
END SUBROUTINE wypos_44 
!# 1803 "wypos.f90"
SUBROUTINE wypos_45 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 1818 "wypos.f90"
END SUBROUTINE wypos_45 
!# 1820 "wypos.f90"
SUBROUTINE wypos_46 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 1835 "wypos.f90"
END SUBROUTINE wypos_46 
!# 1837 "wypos.f90"
SUBROUTINE wypos_47 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4u') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4v') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4w') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4x') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4y') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4z') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 1948 "wypos.f90"
END SUBROUTINE wypos_47 
!# 1950 "wypos.f90"
SUBROUTINE wypos_48 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
!# 2007 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.75_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 2059 "wypos.f90"
END SUBROUTINE wypos_48 
!# 2061 "wypos.f90"
SUBROUTINE wypos_49 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4p') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 2136 "wypos.f90"
END SUBROUTINE wypos_49 
!# 2138 "wypos.f90"
SUBROUTINE wypos_50 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
!# 2195 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 2247 "wypos.f90"
END SUBROUTINE wypos_50 
!# 2249 "wypos.f90"
SUBROUTINE wypos_51 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 2300 "wypos.f90"
END SUBROUTINE wypos_51 
!# 2302 "wypos.f90"
SUBROUTINE wypos_52 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
!# 2324 "wypos.f90"
END SUBROUTINE wypos_52 
!# 2326 "wypos.f90"
SUBROUTINE wypos_53 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 2365 "wypos.f90"
END SUBROUTINE wypos_53 
!# 2367 "wypos.f90"
SUBROUTINE wypos_54 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 2393 "wypos.f90"
END SUBROUTINE wypos_54 
!# 2395 "wypos.f90"
SUBROUTINE wypos_55 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 2434 "wypos.f90"
END SUBROUTINE wypos_55 
!# 2436 "wypos.f90"
SUBROUTINE wypos_56 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=inp(1)
            ENDIF
!# 2459 "wypos.f90"
END SUBROUTINE wypos_56 
!# 2461 "wypos.f90"
SUBROUTINE wypos_57 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
!# 2484 "wypos.f90"
END SUBROUTINE wypos_57 
!# 2486 "wypos.f90"
SUBROUTINE wypos_58 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 2523 "wypos.f90"
END SUBROUTINE wypos_58 
!# 2525 "wypos.f90"
SUBROUTINE wypos_59 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF
!# 2558 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 2586 "wypos.f90"
END SUBROUTINE wypos_59 
!# 2588 "wypos.f90"
SUBROUTINE wypos_60 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 2608 "wypos.f90"
END SUBROUTINE wypos_60 
!# 2610 "wypos.f90"
SUBROUTINE wypos_61 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
!# 2626 "wypos.f90"
END SUBROUTINE wypos_61 
!# 2628 "wypos.f90"
SUBROUTINE wypos_62 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
!# 2648 "wypos.f90"
END SUBROUTINE wypos_62 
!# 2650 "wypos.f90"
SUBROUTINE wypos_63 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
!# 2685 "wypos.f90"
END SUBROUTINE wypos_63 
!# 2687 "wypos.f90"
SUBROUTINE wypos_64 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 2718 "wypos.f90"
END SUBROUTINE wypos_64 
!# 2720 "wypos.f90"
SUBROUTINE wypos_65 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 2796 "wypos.f90"
END SUBROUTINE wypos_65 
!# 2798 "wypos.f90"
SUBROUTINE wypos_66 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 2853 "wypos.f90"
END SUBROUTINE wypos_66 
!# 2855 "wypos.f90"
SUBROUTINE wypos_67 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
!# 2918 "wypos.f90"
END SUBROUTINE wypos_67 
!# 2920 "wypos.f90"
SUBROUTINE wypos_68 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
!# 2961 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 2997 "wypos.f90"
END SUBROUTINE wypos_68 
!# 2999 "wypos.f90"
SUBROUTINE wypos_69 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16o') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 3066 "wypos.f90"
END SUBROUTINE wypos_69 
!# 3068 "wypos.f90"
SUBROUTINE wypos_70 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
!# 3106 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 3138 "wypos.f90"
END SUBROUTINE wypos_70 
!# 3140 "wypos.f90"
SUBROUTINE wypos_71 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 3203 "wypos.f90"
END SUBROUTINE wypos_71 
!# 3205 "wypos.f90"
SUBROUTINE wypos_72 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 3252 "wypos.f90"
END SUBROUTINE wypos_72 
!# 3254 "wypos.f90"
SUBROUTINE wypos_73 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF
!# 3282 "wypos.f90"
END SUBROUTINE wypos_73 
!# 3284 "wypos.f90"
SUBROUTINE wypos_74 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF
!# 3327 "wypos.f90"
END SUBROUTINE wypos_74 
!# 3329 "wypos.f90"
SUBROUTINE wypos_75 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 3348 "wypos.f90"
END SUBROUTINE wypos_75 
!# 3350 "wypos.f90"
SUBROUTINE wypos_77 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 3369 "wypos.f90"
END SUBROUTINE wypos_77 
!# 3371 "wypos.f90"
SUBROUTINE wypos_79 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 3386 "wypos.f90"
END SUBROUTINE wypos_79 
!# 3388 "wypos.f90"
SUBROUTINE wypos_80 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 3399 "wypos.f90"
END SUBROUTINE wypos_80 
!# 3401 "wypos.f90"
SUBROUTINE wypos_81 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 3436 "wypos.f90"
END SUBROUTINE wypos_81 
!# 3438 "wypos.f90"
SUBROUTINE wypos_82 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 3469 "wypos.f90"
END SUBROUTINE wypos_82 
!# 3471 "wypos.f90"
SUBROUTINE wypos_83 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 3522 "wypos.f90"
END SUBROUTINE wypos_83 
!# 3524 "wypos.f90"
SUBROUTINE wypos_84 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 3571 "wypos.f90"
END SUBROUTINE wypos_84 
!# 3573 "wypos.f90"
SUBROUTINE wypos_85 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
!# 3606 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 3634 "wypos.f90"
END SUBROUTINE wypos_85 
!# 3636 "wypos.f90"
SUBROUTINE wypos_86 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
!# 3669 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 3697 "wypos.f90"
END SUBROUTINE wypos_86 
!# 3699 "wypos.f90"
SUBROUTINE wypos_87 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 3739 "wypos.f90"
END SUBROUTINE wypos_87 
!# 3741 "wypos.f90"
SUBROUTINE wypos_88 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
!# 3770 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 3794 "wypos.f90"
END SUBROUTINE wypos_88 
!# 3796 "wypos.f90"
SUBROUTINE wypos_89 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF
!# 3863 "wypos.f90"
END SUBROUTINE wypos_89 
!# 3865 "wypos.f90"
SUBROUTINE wypos_90 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 3896 "wypos.f90"
END SUBROUTINE wypos_90 
!# 3898 "wypos.f90"
SUBROUTINE wypos_91 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.375_DP
            ENDIF
!# 3917 "wypos.f90"
END SUBROUTINE wypos_91 
!# 3919 "wypos.f90"
SUBROUTINE wypos_92 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF
!# 3930 "wypos.f90"
END SUBROUTINE wypos_92 
!# 3932 "wypos.f90"
SUBROUTINE wypos_93 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ENDIF
!# 3999 "wypos.f90"
END SUBROUTINE wypos_93 
!# 4001 "wypos.f90"
SUBROUTINE wypos_94 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 4032 "wypos.f90"
END SUBROUTINE wypos_94 
!# 4034 "wypos.f90"
SUBROUTINE wypos_95 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.625_DP
            ENDIF
!# 4053 "wypos.f90"
END SUBROUTINE wypos_95 
!# 4055 "wypos.f90"
SUBROUTINE wypos_96 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF
!# 4066 "wypos.f90"
END SUBROUTINE wypos_96 
!# 4068 "wypos.f90"
SUBROUTINE wypos_97 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ENDIF
!# 4115 "wypos.f90"
END SUBROUTINE wypos_97 
!# 4117 "wypos.f90"
SUBROUTINE wypos_98 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=-inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF
!# 4148 "wypos.f90"
END SUBROUTINE wypos_98 
!# 4150 "wypos.f90"
SUBROUTINE wypos_99 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4181 "wypos.f90"
END SUBROUTINE wypos_99 
!# 4183 "wypos.f90"
SUBROUTINE wypos_100( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4202 "wypos.f90"
END SUBROUTINE wypos_100
!# 4204 "wypos.f90"
SUBROUTINE wypos_101( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 4227 "wypos.f90"
END SUBROUTINE wypos_101
!# 4229 "wypos.f90"
SUBROUTINE wypos_102( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 4248 "wypos.f90"
END SUBROUTINE wypos_102
!# 4250 "wypos.f90"
SUBROUTINE wypos_103( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4269 "wypos.f90"
END SUBROUTINE wypos_103
!# 4271 "wypos.f90"
SUBROUTINE wypos_104( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4286 "wypos.f90"
END SUBROUTINE wypos_104
!# 4288 "wypos.f90"
SUBROUTINE wypos_105( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4315 "wypos.f90"
END SUBROUTINE wypos_105
!# 4317 "wypos.f90"
SUBROUTINE wypos_106( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4332 "wypos.f90"
END SUBROUTINE wypos_106
!# 4334 "wypos.f90"
SUBROUTINE wypos_107( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 4357 "wypos.f90"
END SUBROUTINE wypos_107
!# 4359 "wypos.f90"
SUBROUTINE wypos_108( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4378 "wypos.f90"
END SUBROUTINE wypos_108
!# 4380 "wypos.f90"
SUBROUTINE wypos_109( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 4395 "wypos.f90"
END SUBROUTINE wypos_109
!# 4397 "wypos.f90"
SUBROUTINE wypos_110( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 4408 "wypos.f90"
END SUBROUTINE wypos_110
!# 4410 "wypos.f90"
SUBROUTINE wypos_111( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 4473 "wypos.f90"
END SUBROUTINE wypos_111
!# 4475 "wypos.f90"
SUBROUTINE wypos_112( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4535 "wypos.f90"
END SUBROUTINE wypos_112
!# 4537 "wypos.f90"
SUBROUTINE wypos_113( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4564 "wypos.f90"
END SUBROUTINE wypos_113
!# 4566 "wypos.f90"
SUBROUTINE wypos_114( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4590 "wypos.f90"
END SUBROUTINE wypos_114
!# 4592 "wypos.f90"
SUBROUTINE wypos_115( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 4644 "wypos.f90"
END SUBROUTINE wypos_115
!# 4646 "wypos.f90"
SUBROUTINE wypos_116( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4690 "wypos.f90"
END SUBROUTINE wypos_116
!# 4692 "wypos.f90"
SUBROUTINE wypos_117( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ENDIF
!# 4732 "wypos.f90"
END SUBROUTINE wypos_117
!# 4734 "wypos.f90"
SUBROUTINE wypos_118( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 4773 "wypos.f90"
END SUBROUTINE wypos_118
!# 4775 "wypos.f90"
SUBROUTINE wypos_119( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 4819 "wypos.f90"
END SUBROUTINE wypos_119
!# 4821 "wypos.f90"
SUBROUTINE wypos_120( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ENDIF
!# 4860 "wypos.f90"
END SUBROUTINE wypos_120
!# 4862 "wypos.f90"
SUBROUTINE wypos_121( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 4905 "wypos.f90"
END SUBROUTINE wypos_121
!# 4907 "wypos.f90"
SUBROUTINE wypos_122( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF
!# 4931 "wypos.f90"
END SUBROUTINE wypos_122
!# 4933 "wypos.f90"
SUBROUTINE wypos_123( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
!# 4938 "wypos.f90"
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8r') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8s') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8t') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 5021 "wypos.f90"
END SUBROUTINE wypos_123
!# 5023 "wypos.f90"
SUBROUTINE wypos_124( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 5083 "wypos.f90"
END SUBROUTINE wypos_124
!# 5085 "wypos.f90"
SUBROUTINE wypos_125( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
!# 5146 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 5202 "wypos.f90"
END SUBROUTINE wypos_125
!# 5204 "wypos.f90"
SUBROUTINE wypos_126( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
!# 5210 "wypos.f90"
   IF (origin_choice==1) THEN
!# 5212 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
!# 5254 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 5298 "wypos.f90"
END SUBROUTINE wypos_126
!# 5300 "wypos.f90"
SUBROUTINE wypos_127( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 5352 "wypos.f90"
END SUBROUTINE wypos_127
!# 5354 "wypos.f90"
SUBROUTINE wypos_128( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 5393 "wypos.f90"
END SUBROUTINE wypos_128
!# 5395 "wypos.f90"
SUBROUTINE wypos_129( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
!# 5445 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
!# 5448 "wypos.f90"
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 5491 "wypos.f90"
END SUBROUTINE wypos_129
!# 5493 "wypos.f90"
SUBROUTINE wypos_130( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
!# 5526 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 5554 "wypos.f90"
END SUBROUTINE wypos_130
!# 5556 "wypos.f90"
SUBROUTINE wypos_131( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 5632 "wypos.f90"
END SUBROUTINE wypos_131
!# 5634 "wypos.f90"
SUBROUTINE wypos_132( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 5701 "wypos.f90"
END SUBROUTINE wypos_132
!# 5703 "wypos.f90"
SUBROUTINE wypos_133( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.0_DP
               ENDIF
!# 5752 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.00_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 5796 "wypos.f90"
END SUBROUTINE wypos_133
!# 5798 "wypos.f90"
SUBROUTINE wypos_134( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 5859 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 5914 "wypos.f90"
            ENDIF
!# 5916 "wypos.f90"
END SUBROUTINE wypos_134
!# 5918 "wypos.f90"
SUBROUTINE wypos_135( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 5957 "wypos.f90"
END SUBROUTINE wypos_135
!# 5959 "wypos.f90"
SUBROUTINE wypos_136( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 6006 "wypos.f90"
END SUBROUTINE wypos_136
!# 6008 "wypos.f90"
SUBROUTINE wypos_137( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 6045 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 6077 "wypos.f90"
END SUBROUTINE wypos_137
!# 6079 "wypos.f90"
SUBROUTINE wypos_138( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF
!# 6124 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 6164 "wypos.f90"
END SUBROUTINE wypos_138
!# 6166 "wypos.f90"
SUBROUTINE wypos_139( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 6229 "wypos.f90"
END SUBROUTINE wypos_139
!# 6231 "wypos.f90"
SUBROUTINE wypos_140( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF
!# 6282 "wypos.f90"
END SUBROUTINE wypos_140
!# 6284 "wypos.f90"
SUBROUTINE wypos_141( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 6325 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.875_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 6360 "wypos.f90"
            ENDIF
!# 6362 "wypos.f90"
END SUBROUTINE wypos_141
!# 6364 "wypos.f90"
SUBROUTINE wypos_142( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
!# 6397 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.125_DP
               ENDIF
!# 6424 "wypos.f90"
            ENDIF
!# 6426 "wypos.f90"
END SUBROUTINE wypos_142
!# 6428 "wypos.f90"
SUBROUTINE wypos_143( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF
!# 6447 "wypos.f90"
END SUBROUTINE wypos_143
!# 6449 "wypos.f90"
SUBROUTINE wypos_146( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ELSE !If HEXAGONAL
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 6469 "wypos.f90"
END SUBROUTINE wypos_146
!# 6471 "wypos.f90"
SUBROUTINE wypos_147( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
!# 6503 "wypos.f90"
END SUBROUTINE wypos_147
!# 6505 "wypos.f90"
SUBROUTINE wypos_148( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF
!# 6534 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
            ENDIF
!# 6558 "wypos.f90"
END SUBROUTINE wypos_148
!# 6560 "wypos.f90"
SUBROUTINE wypos_149( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 6611 "wypos.f90"
END SUBROUTINE wypos_149
!# 6613 "wypos.f90"
SUBROUTINE wypos_150( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
!# 6644 "wypos.f90"
END SUBROUTINE wypos_150
!# 6646 "wypos.f90"
SUBROUTINE wypos_151( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=5.0_DP/6.0_DP
            ENDIF
!# 6661 "wypos.f90"
END SUBROUTINE wypos_151
!# 6663 "wypos.f90"
SUBROUTINE wypos_152( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=5.0_DP/6.0_DP
            ENDIF
!# 6678 "wypos.f90"
END SUBROUTINE wypos_152
!# 6680 "wypos.f90"
SUBROUTINE wypos_153( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/6.0_DP
            ENDIF
!# 6695 "wypos.f90"
END SUBROUTINE wypos_153
!# 6697 "wypos.f90"
SUBROUTINE wypos_154( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/6.0_DP
            ENDIF
!# 6712 "wypos.f90"
END SUBROUTINE wypos_154
!# 6714 "wypos.f90"
SUBROUTINE wypos_155( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
!# 6743 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
            ENDIF
!# 6767 "wypos.f90"
END SUBROUTINE wypos_155
!# 6769 "wypos.f90"
SUBROUTINE wypos_156( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
!# 6792 "wypos.f90"
END SUBROUTINE wypos_156
!# 6794 "wypos.f90"
SUBROUTINE wypos_157( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 6813 "wypos.f90"
END SUBROUTINE wypos_157
!# 6815 "wypos.f90"
SUBROUTINE wypos_158( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF
!# 6834 "wypos.f90"
END SUBROUTINE wypos_158
!# 6836 "wypos.f90"
SUBROUTINE wypos_159( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF
!# 6851 "wypos.f90"
END SUBROUTINE wypos_159
!# 6853 "wypos.f90"
SUBROUTINE wypos_160( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 6870 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9b') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 6882 "wypos.f90"
END SUBROUTINE wypos_160
!# 6884 "wypos.f90"
SUBROUTINE wypos_161( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
!# 6897 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 6905 "wypos.f90"
END SUBROUTINE wypos_161
!# 6907 "wypos.f90"
SUBROUTINE wypos_162( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 6958 "wypos.f90"
END SUBROUTINE wypos_162
!# 6960 "wypos.f90"
SUBROUTINE wypos_163( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 6999 "wypos.f90"
END SUBROUTINE wypos_163
!# 7001 "wypos.f90"
SUBROUTINE wypos_164( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
!# 7044 "wypos.f90"
END SUBROUTINE wypos_164
!# 7046 "wypos.f90"
SUBROUTINE wypos_165( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
!# 7077 "wypos.f90"
END SUBROUTINE wypos_165
!# 7079 "wypos.f90"
SUBROUTINE wypos_166( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 7120 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='18h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 7156 "wypos.f90"
END SUBROUTINE wypos_166
!# 7158 "wypos.f90"
SUBROUTINE wypos_167( wp, inp, rhombohedral,  tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6e') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ENDIF
!# 7187 "wypos.f90"
            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='18d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 7211 "wypos.f90"
END SUBROUTINE wypos_167
!# 7213 "wypos.f90"
SUBROUTINE wypos_168( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 7232 "wypos.f90"
END SUBROUTINE wypos_168
!# 7234 "wypos.f90"
SUBROUTINE wypos_171( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 7249 "wypos.f90"
END SUBROUTINE wypos_171
!# 7251 "wypos.f90"
SUBROUTINE wypos_172( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF
!# 7266 "wypos.f90"
END SUBROUTINE wypos_172
!# 7268 "wypos.f90"
SUBROUTINE wypos_173( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF
!# 7283 "wypos.f90"
END SUBROUTINE wypos_173
!# 7285 "wypos.f90"
SUBROUTINE wypos_174( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 7336 "wypos.f90"
END SUBROUTINE wypos_174
!# 7338 "wypos.f90"
SUBROUTINE wypos_175( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 7389 "wypos.f90"
END SUBROUTINE wypos_175
!# 7391 "wypos.f90"
SUBROUTINE wypos_176( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
!# 7430 "wypos.f90"
END SUBROUTINE wypos_176
!# 7432 "wypos.f90"
SUBROUTINE wypos_177( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 7491 "wypos.f90"
END SUBROUTINE wypos_177
!# 7493 "wypos.f90"
SUBROUTINE wypos_178( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 7508 "wypos.f90"
END SUBROUTINE wypos_178
!# 7510 "wypos.f90"
SUBROUTINE wypos_179( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.75_DP
            ENDIF
!# 7524 "wypos.f90"
END SUBROUTINE wypos_179
!# 7526 "wypos.f90"
SUBROUTINE wypos_180( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 7573 "wypos.f90"
END SUBROUTINE wypos_180
!# 7575 "wypos.f90"
SUBROUTINE wypos_181( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF
!# 7622 "wypos.f90"
END SUBROUTINE wypos_181
!# 7624 "wypos.f90"
SUBROUTINE wypos_182( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF
!# 7663 "wypos.f90"
END SUBROUTINE wypos_182
!# 7665 "wypos.f90"
SUBROUTINE wypos_183( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
!# 7692 "wypos.f90"
END SUBROUTINE wypos_183
!# 7694 "wypos.f90"
SUBROUTINE wypos_184( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF
!# 7713 "wypos.f90"
END SUBROUTINE wypos_184
!# 7715 "wypos.f90"
SUBROUTINE wypos_185( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 7734 "wypos.f90"
END SUBROUTINE wypos_185
!# 7736 "wypos.f90"
SUBROUTINE wypos_186( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
!# 7755 "wypos.f90"
END SUBROUTINE wypos_186
!# 7757 "wypos.f90"
SUBROUTINE wypos_187( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6n') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF
!# 7820 "wypos.f90"
END SUBROUTINE wypos_187
!# 7822 "wypos.f90"
SUBROUTINE wypos_188( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
!# 7873 "wypos.f90"
END SUBROUTINE wypos_188
!# 7875 "wypos.f90"
SUBROUTINE wypos_189( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 7926 "wypos.f90"
END SUBROUTINE wypos_189
!# 7928 "wypos.f90"
SUBROUTINE wypos_190( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF
!# 7968 "wypos.f90"
END SUBROUTINE wypos_190
!# 7970 "wypos.f90"
SUBROUTINE wypos_191( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12o') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF
!# 8045 "wypos.f90"
END SUBROUTINE wypos_191
!# 8047 "wypos.f90"
SUBROUTINE wypos_192( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF
!# 8102 "wypos.f90"
END SUBROUTINE wypos_192
!# 8104 "wypos.f90"
SUBROUTINE wypos_193( wp, inp, tau )
   
   REAL(dp), INTENT(in) :: inp(3)
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF
!# 8157 "wypos.f90"
END SUBROUTINE wypos_193
!# 8159 "wypos.f90"
SUBROUTINE wypos_194( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ENDIF
!# 8210 "wypos.f90"
END SUBROUTINE wypos_194
!# 8212 "wypos.f90"
SUBROUTINE wypos_195( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF
!# 8255 "wypos.f90"
END SUBROUTINE wypos_195
!# 8257 "wypos.f90"
SUBROUTINE wypos_196( wp, inp, tau )
  CHARACTER(LEN=*), INTENT(in)  :: wp
  REAL(dp), INTENT(in) :: inp(3)
  REAL(dp), INTENT(out) :: tau (3)
  
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
!# 8292 "wypos.f90"
END SUBROUTINE wypos_196
!# 8294 "wypos.f90"
SUBROUTINE wypos_197( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF
!# 8321 "wypos.f90"
END SUBROUTINE wypos_197
!# 8323 "wypos.f90"
SUBROUTINE wypos_198( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF
!# 8334 "wypos.f90"
END SUBROUTINE wypos_198
!# 8336 "wypos.f90"
SUBROUTINE wypos_199( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
!# 8351 "wypos.f90"
END SUBROUTINE wypos_199
!# 8353 "wypos.f90"
SUBROUTINE wypos_200( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 8404 "wypos.f90"
END SUBROUTINE wypos_200
!# 8406 "wypos.f90"
SUBROUTINE wypos_201( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF
!# 8443 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF
!# 8475 "wypos.f90"
END SUBROUTINE wypos_201
!# 8477 "wypos.f90"
SUBROUTINE wypos_202( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 8516 "wypos.f90"
END SUBROUTINE wypos_202
!# 8518 "wypos.f90"
SUBROUTINE wypos_203( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
!# 8550 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ENDIF
            ENDIF
!# 8578 "wypos.f90"
END SUBROUTINE wypos_203
!# 8580 "wypos.f90"
SUBROUTINE wypos_204( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 8615 "wypos.f90"
END SUBROUTINE wypos_204
!# 8617 "wypos.f90"
SUBROUTINE wypos_205( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF
!# 8636 "wypos.f90"
END SUBROUTINE wypos_205
!# 8638 "wypos.f90"
SUBROUTINE wypos_206( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
!# 8661 "wypos.f90"
END SUBROUTINE wypos_206
!# 8663 "wypos.f90"
SUBROUTINE wypos_207( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF
!# 8711 "wypos.f90"
END SUBROUTINE wypos_207
!# 8713 "wypos.f90"
SUBROUTINE wypos_208( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ENDIF
!# 8768 "wypos.f90"
END SUBROUTINE wypos_208
!# 8770 "wypos.f90"
SUBROUTINE wypos_209( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
!# 8813 "wypos.f90"
END SUBROUTINE wypos_209
!# 8815 "wypos.f90"
SUBROUTINE wypos_210( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='16d') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
!# 8850 "wypos.f90"
END SUBROUTINE wypos_210
!# 8852 "wypos.f90"
SUBROUTINE wypos_211( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ENDIF
!# 8895 "wypos.f90"
END SUBROUTINE wypos_211
!# 8897 "wypos.f90"
SUBROUTINE wypos_212( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
!# 8920 "wypos.f90"
END SUBROUTINE wypos_212
!# 8922 "wypos.f90"
SUBROUTINE wypos_213( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.375_DP
               tau(2)=0.375_DP
               tau(3)=0.375_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ENDIF
!# 8945 "wypos.f90"
END SUBROUTINE wypos_213
!# 8947 "wypos.f90"
SUBROUTINE wypos_214( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='12c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.625_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
!# 8987 "wypos.f90"
END SUBROUTINE wypos_214
!# 8989 "wypos.f90"
SUBROUTINE wypos_215( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9032 "wypos.f90"
END SUBROUTINE wypos_215
!# 9034 "wypos.f90"
SUBROUTINE wypos_216( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9073 "wypos.f90"
END SUBROUTINE wypos_216
!# 9075 "wypos.f90"
SUBROUTINE wypos_217( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9110 "wypos.f90"
END SUBROUTINE wypos_217
!# 9112 "wypos.f90"
SUBROUTINE wypos_218( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF
!# 9151 "wypos.f90"
END SUBROUTINE wypos_218
!# 9153 "wypos.f90"
SUBROUTINE wypos_219( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF
!# 9188 "wypos.f90"
END SUBROUTINE wypos_219
!# 9190 "wypos.f90"
SUBROUTINE wypos_220( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='12a') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=0.875_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF
!# 9213 "wypos.f90"
END SUBROUTINE wypos_220
!# 9215 "wypos.f90"
SUBROUTINE wypos_221( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9275 "wypos.f90"
END SUBROUTINE wypos_221
!# 9277 "wypos.f90"
SUBROUTINE wypos_222( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
!# 9318 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ENDIF
!# 9354 "wypos.f90"
END SUBROUTINE wypos_222
!# 9356 "wypos.f90"
SUBROUTINE wypos_223( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9407 "wypos.f90"
END SUBROUTINE wypos_223
!# 9409 "wypos.f90"
SUBROUTINE wypos_224( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
!# 9462 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF
!# 9510 "wypos.f90"
END SUBROUTINE wypos_224
!# 9512 "wypos.f90"
SUBROUTINE wypos_225( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='96k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9563 "wypos.f90"
END SUBROUTINE wypos_225
!# 9565 "wypos.f90"
SUBROUTINE wypos_226( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp 
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='64g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96h') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9608 "wypos.f90"
END SUBROUTINE wypos_226
!# 9610 "wypos.f90"
SUBROUTINE wypos_227( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF
!# 9651 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF
!# 9687 "wypos.f90"
END SUBROUTINE wypos_227
!# 9689 "wypos.f90"
SUBROUTINE wypos_228( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF
!# 9726 "wypos.f90"
            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.875_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF
!# 9758 "wypos.f90"
END SUBROUTINE wypos_228
!# 9760 "wypos.f90"
SUBROUTINE wypos_229( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='48j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='48k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
!# 9811 "wypos.f90"
END SUBROUTINE wypos_229
!# 9813 "wypos.f90"
SUBROUTINE wypos_230( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='16a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16b') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
          END SUBROUTINE wypos_230
END MODULE
