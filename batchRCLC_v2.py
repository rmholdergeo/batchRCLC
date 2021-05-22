# last edited by Robert Holder on 21 May 2021

# ##### Background on the develompent of this code ######
#     The version of the RCLC algorithm used in this code is a modified version of Tom Chacko's 1996
#     original. The modifications have been made by Chris McFarlane (04-97), David Pattison (05-97 and 10-01),
#     and Tom Chacko (05-02).
#
#     This python code was written by Robert Holder in May 2021 as a modfication and expansion of
#     python code and a web application (https://rclc.facsci.ualberta.ca/) developed by
#     Justin Widney at the University of Alberta, as part of a summer research internship
#     under the supervision of Tom Chacko in 2018. The reason Holder developed this code
#     was to run to be able to run large batches of RCLC calculations with a single file input
#     (a list of mineral compositions).
#
#     If you use this code, please cite:
#         The original development of this idea:
#             Chacko T, Lamb M, Farquhar J (1996) Ultra-high temperature metamorphism in the Kerala Khondalite Belt. In: Santosh M &
#               Yoshida M (eds) Gondwana Research Group Memoir 3: The Archean and Proterozoic Terrains in Southern India within
#               Gondwana, pp. 157–165. Osaka, Japan: Field Research Publishers.
#         Expansion of the idea and the BASIC version of the algorithm:
#             Pattison DRM, Chacko T, Farquhar J, McFarlane CRM (2003). Temperatures of granulite facies metamorphism: constraints
#               from experimental phase equilibria and thermobarometry corrected for late Fe–Mg exchange. Journal of Petrology 44, 867–900.
#         The first python version of the algorithm and the accompanying web application:
#             Pattison DRM, Chacko T, Farquhar J, McFarlane CRM, Widney J (2019) Program “RCLC”: Garnet–Orthopyroxene Thermobarometry
#               Corrected for Late Fe–Mg Exchange. Journal of Petrology 60, 1107–1108. https://doi.org/10.1093/petrology/egz018
#         This code (allowing for batch processing of large quantities of data):
#             Horton F, Holder RM, Swindle C (2021) An extensive record of orogenesis recorded in a Madagascar granulite. Journal
#               of Metamorphic Geology. [In revision at the time this document was written.]
#
# ##### What this code does #####
#     RCLC calculates P and T of metamorphism from garnet-othropyroxene-plagioclase-quartz
#     equilibria (the GOPS barometer and Al-in-opx thermometer), corrected for Fe-Mg exchange
#     among minerals (garnet and orthopyroxene as well as cordierite and biotite if they
#     are relevant for a particular sample). In other words, this is Al-in-orthopyroxene
#     thermobarometry, corrected for retrograde Fe-Mg diffusion.
#
#     The Fe-Mg exchange correction is done by adjusting the Fe-Mg ratios of the mafic
#     minerals according to their modes (mass balance) and Fe-Mg exchange equilibria:
#     garnet-orthopyroxene, garnet-cordierite, and garnet-biotite
#
#     The Thermodynamic data (end-members and a-x relationships for solid solutions) are
#     from TWQ202b (files BA96A.DAT and BA96A.SLN), except for BT, which uses the solution
#     model of McMullin (1991 CanMin). For references, see Berman (1991 CanMin), Berman
#     and Aranovich (1996 CMP), and Aranovich and Berman (1997 AmMin).
#
#     This code has two 'runmodes' in which any number of desired calculations can be done:
#
#       1) PT is calculated sequentially from input mineral compositions.
#           For example, if 4 garnet, 4 orthopyroxene, and 4 plagioclase compositions were
#           entered into the input file, the code would calculate PT for:
#               gar1-opx1-pl1, gar2-opx2-pl2, gar3-opx3-pl3, gar4-opx4-pl4.
#           (biotite and cordierite could also be included)
#           This runmode is only available if an equal number of compositions were
#           input for each mineral.
#
#       2) PT is calculated for EVERY POSSIBLE combination of input mineral compositions.
#           For example if 4 garnet, 3 orthopyroxene, and 2 plagioclase compositions were
#           entered into the input file, the code would calculate PT for:
#               gar1-opx1-pl1, gar1-opx1-pl2, gar1-opx2-pl1, gar1-opx2-pl2, gar1-opx3-pl1, gar1-opx3-pl2
#               gar2-opx1-pl1, gar2-opx1-pl2, gar2-opx2-pl1, gar2-opx2-pl2, gar1-opx3-pl1, gar2-opx3-pl2
#               gar3-opx1-pl1, gar3-opx1-pl2, gar3-opx2-pl1, gar3-opx2-pl2, gar3-opx3-pl1, gar1-opx3-pl2
#               gar4-opx1-pl1, gar4-opx1-pl2, gar4-opx2-pl1, gar4-opx2-pl2, gar4-opx3-pl1, gar4-opx3-pl2
#           (biotite and cordierite could also be included)
#           If an unequal number of compositions were input for each mineral, the code will run this
#           mode by default.
#           A great advantage of this runmode is to determine the senstivity of the calculated PT to the choice
#           of input mineral compositions.
#
# ##### What you need to run this code ######
#     The required input data are mineral cations for garnet (normalized to 12 O),
#     orthpyroxene (normalized to 6 O), and plagioclase (normalized to 8 O) as well as
#     modes for each of these minerals. Optional input data are cations in cordierite
#     (normalized to 18 O) and biotite (normalized to 11 O [10 O and 2 (OH,Cl,F]); if
#     cordierite or biotite are included, their modes must also be included.
#
#     Mineral compositions input files are tab-delmited text documents: 'gar.txt',
#     'opx.txt', 'pl.txt', 'bt.txt', and 'crd.txt'. Any number of mineral compositions
#     can be included in the input files. The order of cations in the input files is
#     the same for all minerals: Si, Ti, Al, Cr, Fe3, Fe2, Mn, Mg, Ca, Na, K. If any
#     of these elements was not measured / was not present, a value of 0 should be entered.
#     Example input of five orthopyroxene measurements (Fe3 calculated using Tim
#     Holland's AX2 program):
#
#       Si	Ti	Al	Cr	Fe3	Fe2	Mn	Mg	Ca	Na	K
#       1.798	0.006	0.356	0.001	0.035	0.56	0.003	1.238	0.004	0	0
#       1.791	0.006	0.364	0.002	0.039	0.556	0.002	1.236	0.004	0	0
#       1.803	0.005	0.36	0.002	0.021	0.571	0.002	1.23	0.004	0.001	0
#       1.787	0.005	0.366	0.001	0.049	0.546	0.002	1.238	0.005	0.001	0
#       1.797	0.009	0.366	0.001	0.022	0.569	0.002	1.229	0.005	0.001	0
#
#     Ferric iron: Ferric iron in orthopyroxene, in particular, can moderately influence the calculated
#     Al-in-opx temperature (although typically, the effect is < 50C). The input files allow for entering
#     Fe3 and Fe2 contents of minerals separately. If Fe3 is not being considered, enter zeros for that
#     column of the input file and put all Fe as Fe2.
#
# ##### Uncertainty ######
#     quantifying and reporting uncertainty in phase-equilibrium calculations is
#     very difficult. This code does not output an uncertainty. The commonly quoted
#     rule-of-thumb +/- 50 C and +/- 1 kbar is probably appropriate
#     (Pattison and Chacko, pers comm); however, due to the attempt to correct for
#     diffusional re-equilibration mineral compositions during cooling, this
#     thermobarometer is likely more accurate than most other methods of thermobarometry
#     applied to granulite-facies rocks.
#
#     Al-in-orthpyroxene thermobarometry (as used in this code) has been
#     evaluated in comparison with other methods of thermobarometry on ultrahigh-temperature
#     rocks. It agrees well with estimates of peak metamorphism from petrographic interpretations
#     of peak mineral assemblages in equilibrium-assemblage diagrams (EADs: aka pseudosections), Al-in-opx
#     isopleths in EADs, as well as the highest Zr concentrations measured in rutile.
#     Although a straight-forward, easily quotable uncertainty cannot be give, this
#     method of thermobarometry appears to be one of the most reliable for constraining
#     the peak PT of granulite-facies metamorphism.
#
#     There is additional uncertainty in the choice of which mineral compositions to use.
#     Runmode 2 of this code can be used to evaluate the sensitivty of the calculation to
#     such choices by calculating PT for all possible combinations of input mineral data.
#     e.g. Horton, Holder, and Swindle (2021 JMG)


########################################################
################# IMPORTING functions ##################
########################################################
from re import search as rsearch
from numpy import array as nparray
from math import exp
from math import log
from csv import writer as csvwriter
########################################################
################# END IMPORTING LIBRARIES ##############
########################################################


########################################################
########### DEFINE FUNCTIONS FOR THE PROGRAM ###########
########################################################
def ALOPX1(): # Al site occupancy model 1
    global aXFEOPX, aXMGOPX, aXAL_M1
    aXFEOPX, aXMGOPX, aXAL_M1  = aFE2OPX/2, aMGOPX/2, (aALOPX-(2-aSIOPX))/2

def ALOPX2(): # Al site occupancy model 2
    global aXFEOPX, aXMGOPX, aXAL_M1
    aXFEOPX, aXMGOPX, aXAL_M1 = aFE2OPX/2, aMGOPX/2, (aALOPX/2)/2

def ALOPX3(): # Al site occupancy model 3
    global aXFEOPX, aXMGOPX, aXAL_M1
    aXFEOPX = aFE2OPX / (aFE2OPX + aMGOPX + aMNOPX + aCAOPX + (aALOPX / 2))
    aXMGOPX = aMGOPX / (aFE2OPX + aMGOPX + aMNOPX + aCAOPX + (aALOPX / 2))
    aXAL_M1 = (aALOPX / 2) / (aFE2OPX + aMGOPX + aMNOPX + aCAOPX + (aALOPX / 2))

def ALOPX4(): # Al site occupancy model 4
    global aXFEOPX, aXMGOPX, aXAL_M1
    aXFEOPX, aXMGOPX, aXAL_M1 = aFE2OPX/2, aMGOPX/2, (aALOPX-aFE3OPX-aCROPX-(2*aTIOPX))/4

def CP(): # calculates H and S of minerals at T
    global DATASET, HALM, HPY, HGR, HAN, HBQ, HEN, HFS, HALOPX, HPHL, HANN, HCRD, HFECRD
    global SALM, SPY, SGR, SAN, SBQ, SEN, SFS, SALOPX, SPHL, SANN, SCRD, SFECRD, TK, P
    # STANDARD STATE ENTHALPIES AND HEAT CAPACITY EXPRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN
    HALM = DATASET[0][0] + ( (DATASET[0][2] * (TK - 298.15)) + ((2 * DATASET[0][ 3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[0][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[0][5] * ((TK**-2) - (298.15**-2))))
    HPY = DATASET[1][0] + ((DATASET[1][2] * (TK - 298.15)) + ((2 * DATASET[1][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[1][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[1][5] * ((TK**-2) - (298.15**-2))))
    HGR = DATASET[2][0] + ((DATASET[2][2] * (TK - 298.15)) + ((2 * DATASET[2][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[2][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[2][5] * ((TK**-2) - (298.15**-2))))
    HAN = DATASET[3][0] + ((DATASET[3][2] * (TK - 298.15)) + ((2 * DATASET[3][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[3][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[3][5] * ((TK**-2) - (298.15**-2))))
    HBQ = DATASET[4][ 0] + ((DATASET[4][2] * (TK - 298.15)) + ((2 * DATASET[4][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[4][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[4][5] * ((TK**-2) - (298.15**-2))))
    HEN = DATASET[5][ 0] + ((DATASET[5][ 2] * (TK - 298.15)) + ((2 * DATASET[5][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[5][ 4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[5][5] * ((TK**-2) - (298.15**-2))))
    HFS = DATASET[6][ 0] + ((DATASET[6][ 2] * (TK - 298.15)) + ((2 * DATASET[6][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[6][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[6][ 5] * ((TK**-2) - (298.15**-2))))
    HALOPX = DATASET[7][0] + ((DATASET[7][2] * (TK - 298.15)) + ((2 * DATASET[7][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[7][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[7][5] * ((TK**-2) - (298.15**-2))))
    HPHL = DATASET[8][ 0] + ((DATASET[8][ 2] * (TK - 298.15)) + ((2 * DATASET[8][ 3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[8][ 4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[8][ 5] * ((TK**-2) - (298.15**-2))))
    HANN = DATASET[9][0] + ((DATASET[9][2] * (TK - 298.15)) + ((2 * DATASET[9][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[9][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[9][5] * ((TK**-2) - (298.15**-2))))
    HCRD = DATASET[10][0] + ((DATASET[10][2] * (TK - 298.15)) + ((2 * DATASET[10][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[10][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[10][5] * ((TK**-2) - (298.15**-2))))
    HFECRD = DATASET[11][ 0] + ((DATASET[11][ 2] * (TK - 298.15)) + ((2 * DATASET[11][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[11][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[11][5] * ((TK**-2) - (298.15**-2))))
    # STANDARD STATE ENTROPIES AND HEAT CAPACITY EXPRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN
    SALM = DATASET[0][1] + ((DATASET[0][ 2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[0][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[0][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[0][5]) * ((TK ** -3) - (298.15 ** -3))))
    SPY = DATASET[1][1] + ((DATASET[1][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[1][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[1][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[1][5]) * ((TK ** -3) - (298.15 ** -3))))
    SGR = DATASET[2][1] + ((DATASET[2][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[2][3]) * ((TK ** -.5) - (298.15 ** -.5))) -  ((.5 * DATASET[2][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[2][5]) * ((TK ** -3) - (298.15 ** -3))))
    SAN = DATASET[3][1] + ((DATASET[3][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[3][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[3][4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[3][5]) * ((TK ** -3) - (298.15 ** -3))))
    SBQ = DATASET[4][1] + ((DATASET[4][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[4][ 3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[4][ 4]) * ((TK**-2) - (298.15**-2))) -    (((1 / 3) * DATASET[4][5]) * ((TK ** -3) - (298.15 ** -3))))
    SEN = DATASET[5][1] + ((DATASET[5][2] * ((log(TK)) - (log(298.15)))) -     ((2 * DATASET[5][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[5][ 4]) * ((TK**-2) - (298.15**-2))) -    (((1 / 3) * DATASET[5][ 5]) * ((TK ** -3) - (298.15 ** -3))))
    SFS = DATASET[6][1] + ((DATASET[6][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[6][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[6][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[6][5]) * ((TK ** -3) - (298.15 ** -3))))
    SALOPX = DATASET[7][1] + ((DATASET[7][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[7][3]) * ((TK ** -.5) - (298.15 ** -.5))) -  ((.5 * DATASET[7][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[7][5]) * ((TK ** -3) - (298.15 ** -3))))
    SPHL = DATASET[8][1] + ((DATASET[8][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[8][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[8][ 4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[8][5]) * ((TK ** -3) - (298.15 ** -3))))
    SANN = DATASET[9][1] + ((DATASET[9][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[9][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[9][4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[9][5]) * ((TK ** -3) - (298.15 ** -3))))
    SCRD = DATASET[10][ 1] + ((DATASET[10][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[10][ 3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[10][ 4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[10][5]) * ((TK ** -3) - (298.15 ** -3))))
    SFECRD = DATASET[11][1] + ((DATASET[11][2] * ((log(TK)) - (log(298.15)))) - ((2 * DATASET[11][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[11][ 4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[11][5]) * ((TK ** -3) - (298.15 ** -3))))

def VOLUMEPT(): # Calculates V of minerals at P and T
    global PBARS, VALM, VPY, VGR, VAN, VBQ, VEN, VFS, VALOPX, VPHL, VANN, VCRD, VFECRD, TK, P
    # STANDARD STATE VOLUMES AND EXPANSION AND COMPRESSIBILITY EXPRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN
    VALM = 11.524 * (1 + (.0000185989054 * (TK - 298)) + (7.4711E-09 * ((TK - 298) ** 2)) + (-.000000570324 * PBARS) + (4.344E-13 * (PBARS ** 2)))
    VPY = 11.311 * (1 + (0.0000225186544 * (TK - 298)) + (3.7044E-09 * ((TK - 298) ** 2)) + (-.000000576209 * PBARS) + (4.42E-13 * (PBARS ** 2)))
    VGR = 12.538 * (1 + (0.0000189942017 * (TK - 298)) + (7.9756E-09 * ((TK - 298) ** 2)) + (-.0000006539136 * PBARS) + (1.635E-12 * (PBARS ** 2)))
    VAN = 10.075 * (1 + (.0000109181141 * (TK - 298)) + (4.1985E-09 * ((TK - 298) ** 2)) + (-.0000012724268 * PBARS) + (3.1762E-12 * (PBARS ** 2)))
    VBQ = 2.37 * (1 + (0 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000012382672 * PBARS) + (7.0871E-12 * (PBARS ** 2)))
    VEN = 3.133 * (1 + (.0000246558172 * (TK - 298)) + (7.467E-09 * ((TK - 298) ** 2)) + (-.0000007493458 * PBARS) + (4.467E-13 * (PBARS ** 2)))
    VFS = 3.295 * (1 + (.0000314064017 * (TK - 298)) + (8.04E-09 * ((TK - 298) ** 2)) + (-.0000009111044 * PBARS) + (3.034E-13 * (PBARS ** 2)))
    VALOPX = 3.093 * (1 + (.0000246558172 * (TK - 298)) + (7.467E-09 * ((TK - 298) ** 2)) + (-.0000007493458 * PBARS) + (4.467E-13 * (PBARS ** 2)))
    VPHL = 14.971 * (1 + (.0000344473262 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016969784 * PBARS) + (0 * (PBARS ** 2)))
    VANN = 15.487 * (1 + (.0000344473262 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016969784 * PBARS) + (0 * (PBARS ** 2)))
    VCRD = 23.311 * (1 + (.0000030028742 * (TK - 298)) + (1.8017E-09 * ((TK - 298) ** 2)) + (-.0000011582515 * PBARS) + (0 * (PBARS ** 2)))
    VFECRD = 23.706 * (1 + (.0000042647431 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016998228 * PBARS) + (0 * (PBARS ** 2)))

def GARNET():
    global AGR, APY, AAL, GAMMAGAR, PBARS, TK, P, XCAGAR, XMGGAR, XFEGAR, XMNGAR
    # GARNET ACTIVITIES FOR CA-FE-MG-MN GARNET WITH THE MODEL IN TWQ202B - BA96a.SLN OF BERMAN
    X1, X2, X3, X4 = XCAGAR, XMGGAR, XFEGAR, XMNGAR
    W112 = (85529) - (TK * 18.79) + (PBARS * .21)
    W122 = 50874.9 - (TK * 18.79) + (PBARS * .02)
    W113 = 24025.5 - (TK * 9.43) + (PBARS * .17)
    W133 = 9876.2 - (TK * 9.43) + (PBARS * .09)
    W223 = 1307.4 + (PBARS * .01)
    W233 = 2092.4 + (PBARS * .06)
    W123 = 86852.8 - (TK * 28.22) + (PBARS * .28)
    W124 = 82759.9 - (TK * 28.79) + (PBARS * .1)
    W134 = 7053.9 + (TK * 30.01) + (PBARS * .13)
    W234 = (6361) + (TK * 29.44) + (PBARS * .04)
    W224 = (14558) - (TK * (10)) + (PBARS * .04)
    W244 = (14558) - (TK * (10)) + (PBARS * .04)
    W344 = (158) + (TK * 35.1) + (PBARS * .04)
    W334 = (-(19952)) + (TK * 43.78) + (PBARS * .04)

    # This "fixme" comment was left in Widney's code, although it isn't clear why.
    # Things seem to work corrrectly. I have left it here in case someone else notices anything wrong
    ### FIXME: ###
    W114 = 0
    W144 = 0

    TERM1GR = (W112 * ((2 * X1 * X2) - (2 * (X1 ** 2) * X2))) + (W122 * ((X2 ** 2) - (2 * X1 * (X2 ** 2))))
    TERM2GR = (W113 * ((2 * X1 * X3) - (2 * (X1 ** 2) * X3))) + (W133 * ((X3 ** 2) - (2 * X1 * (X3 ** 2))))
    TERM3GR = (W114 * ((2 * X1 * X4) - (2 * (X1 ** 2) * X4))) + (W144 * ((X4 ** 2) - (2 * X1 * (X4 ** 2))))
    TERM4GR = (W223 * ((-2) * (X2 ** 2) * X3)) + (W233 * ((-2) * X2 * (X3 ** 2)))
    TERM5GR = (W224 * ((-2) * (X2 ** 2) * X4)) + (W244 * ((-2) * X2 * (X4 ** 2)))
    TERM6GR = (W334 * ((-2) * (X3 ** 2) * X4)) + (W344 * ((-2) * X3 * (X4 ** 2)))
    TERM7GR = (W123 * ((X2 * X3) - (2 * X1 * X2 * X3))) + (W124 * ((X2 * X4) - (2 * X1 * X2 * X4)))
    TERM8GR = (W134 * ((X3 * X4) - (2 * X1 * X3 * X4))) + (W234 * ((-2) * X2 * X3 * X4))
    SUMGR = TERM1GR + TERM2GR + TERM3GR + TERM4GR + TERM5GR + TERM6GR + TERM7GR + TERM8GR
    GAMMAGR = exp(SUMGR / (3 * 8.314 * TK))
    AGR = ((X1 * GAMMAGR) ** 3)

    TERM1PY = (W112 * ((X1 ** 2) - (2 * X2 * (X1 ** 2)))) + (W122 * ((2 * X1 * X2) - (2 * (X2 ** 2) * X1)))
    TERM2PY = (W113 * ((-2) * (X1 ** 2) * X3)) + (W133 * ((-2) * X1 * (X3 ** 2)))
    TERM3PY = (W114 * ((-2) * (X1 ** 2) * X4)) + (W144 * ((-2) * X1 * (X4 ** 2)))
    TERM4PY = (W223 * ((2 * X2 * X3) - (2 * (X2 ** 2) * X3))) + (W233 * ((X3 ** 2) - (2 * X2 * (X3 ** 2))))
    TERM5PY = (W224 * ((2 * X2 * X4) - (2 * (X2 ** 2) * X4))) + (W244 * ((X4 ** 2) - (2 * X2 * (X4 ** 2))))
    TERM6PY = (W334 * ((-2) * (X3 ** 2) * X4)) + (W344 * ((-2) * X3 * (X4 ** 2)))
    TERM7PY = (W123 * ((X1 * X3) - (2 * X1 * X2 * X3))) + (W124 * ((X1 * X4) - (2 * X1 * X2 * X4)))
    TERM8PY = (W134 * ((-2) * X1 * X3 * X4)) + (W234 * ((X3 * X4) - (2 * X2 * X3 * X4)))
    SUMPY = TERM1PY + TERM2PY + TERM3PY + TERM4PY + TERM5PY + TERM6PY + TERM7PY + TERM8PY
    GAMMAPY = exp(SUMPY / (3 * 8.314 * TK))
    APY = ((X2 * GAMMAPY) ** 3)

    TERM1AL = (W112 * ((-2) * (X1 ** 2) * X2)) + (W122 * ((-2) * X1 * (X2 ** 2)))
    TERM2AL = (W113 * ((X1 ** 2) - (2 * X3 * (X1 ** 2)))) + (W133 * ((2 * X1 * X3) - (2 * (X3 ** 2) * X1)))
    TERM3AL = (W114 * ((-2) * (X1 ** 2) * X4)) + (W144 * ((-2) * X1 * (X4 ** 2)))
    TERM4AL = (W223 * ((X2 ** 2) - (2 * X3 * (X2 ** 2)))) + (W233 * ((2 * X2 * X3) - (2 * (X3 ** 2) * X2)))
    TERM5AL = (W224 * ((-2) * (X2 ** 2) * X4)) + (W244 * ((-2) * X2 * (X4 ** 2)))
    TERM6AL = (W334 * ((2 * X3 * X4) - (2 * (X3 ** 2) * X4))) + (W344 * ((X4 ** 2) - (2 * X3 * (X4 ** 2))))
    TERM7AL = (W123 * ((X1 * X2) - (2 * X1 * X2 * X3))) + (W124 * ((-2) * X1 * X2 * X4))
    TERM8AL = (W134 * ((X1 * X4) - (2 * X1 * X3 * X4))) + (W234 * ((X2 * X4) - (2 * X2 * X3 * X4)))
    SUMAL = TERM1AL + TERM2AL + TERM3AL + TERM4AL + TERM5AL + TERM6AL + TERM7AL + TERM8AL
    GAMMAAL = exp(SUMAL / (3 * 8.314 * TK))
    AAL = ((X3 * GAMMAAL) ** 3)
    GAMMAGAR = GAMMAAL / GAMMAPY

def PLAGIOCLASE():
    global AAN, XAB, XSAN, XAN, TK, P
    # THIS SUBROUTINE CALCULATES PLAGIOCLASE ACTIVITIES WITH THE MODEL
    # IN TWQ202B - BA96a.SLN OF BERMAN. IT IS THE MODEL OF FUHRMAN AND LINDSLEY (1988)
    # WITH WORABAN MODIFIED BY BERMAN FOR TWQ202B.
    WABOR = 18.81 - (TK * .0103) + (P * .39)
    WORAB = 27.32 - (TK * .0103) + (P * .39)
    WABAN = 28.226
    WANAB = 8.471
    WANOR = 52.468 - (P * .12)
    WORAN = 47.396
    WORABAN = 100.0455 - (TK * .0103) - (P * .76)
    FIRSTTERM = WORAB * (XAB * XSAN * (.5 - XAN - (2 * XAB)))
    SECONDTERM = WABOR * (XAB * XSAN * (.5 - XAN - (2 * XSAN)))
    THIRDTERM = WORAN * ((2 * XSAN * XAN * (1 - XAN)) + (XAB * XSAN * (.5 - XAN)))
    FOURTHTERM = WANOR * (((XSAN ** 2) * (1 - (2 * XAN))) + (XAB * XSAN * (.5 - XAN)))
    FifTHTERM = WABAN * ((2 * XAB * XAN * (1 - XAN)) + (XAB * XSAN * (.5 - XAN)))
    SIXTHTERM = WANAB * (((XAB ** 2) * (1 - (2 * XAN))) + (XAB * XSAN * (.5 - XAN)))
    SEVENTHTERM = WORABAN * (XSAN * XAB * (1 - (2 * XAN)))
    AAN = exp((FIRSTTERM + SECONDTERM + THIRDTERM + FOURTHTERM + FifTHTERM + SIXTHTERM + SEVENTHTERM) / (.008314 * TK))
    AAN = (XAN * (((1 + XAN) ** 2) / 4)) * AAN

def BIOTITE():
    global GAMMABT, TK, XFEBT, XALBT, XTIBT, XMGBT
    # THIS SUBROUTINE CALCULATES ANNITE AND PHLOGOPITE ACTIVITIES WITH THE MODEL
    # OF MCMULLIN (1991). IT IS THUS DIFFERENT FROM THE MODEL IN TWQ202B - BA96a.SLN OF BERMAN.
    WMGFE = 0
    WMGTI = 58.865
    WMGAL = 75
    WFETI = 30.921
    WFEAL = 63.721
    WTIAL = 0
    RTGAMMAMGBT = ((XFEBT ** 2) * WMGFE) + ((XTIBT ** 2) * WMGTI) + ((XALBT ** 2) * WMGAL) + (XFEBT * XTIBT * (WMGFE + WMGTI - WFETI)) + (XFEBT * XALBT * (WMGFE + WMGAL - WFEAL)) + (XTIBT * XALBT * (WMGTI + WMGAL - WTIAL))
    RTGAMMAFEBT = ((XMGBT ** 2) * WMGFE) + ((XTIBT ** 2) * WFETI) + ((XALBT ** 2) * WFEAL) + (XMGBT * XTIBT * (WMGFE + WFETI - WMGTI)) + (XMGBT * XALBT * (WMGFE + WFEAL - WMGAL)) + (XTIBT * XALBT * (WFETI + WFEAL - WTIAL))
    GAMMAMGBT = exp(RTGAMMAMGBT / (.008314 * TK))
    GAMMAFEBT = exp(RTGAMMAFEBT / (.008314 * TK))
    GAMMABT = GAMMAMGBT / GAMMAFEBT

def CORDIERITE():
    global GAMMACRD, TK, MGRATIOCRD
    # THIS SUBROUTINE CALCULATES MGCRD AND FECRD ACTIVITIES WITH THE MODEL
    # IN TWQ202B - BA96a.SLN OF BERMAN
    WCRD = -1754.7
    RTGAMMAMGCRD = WCRD * ((1 - MGRATIOCRD) ** 2)
    RTGAMMAFECRD = WCRD * (MGRATIOCRD ** 2)
    GAMMAMGCRD = exp(RTGAMMAMGCRD / (8.314 * TK))
    GAMMAFECRD = exp(RTGAMMAFECRD / (8.314 * TK))
    GAMMACRD = GAMMAMGCRD / GAMMAFECRD

def ORTHOPYROXENE():
    global AEN, AFS, AALOPX, GAMMAOPX, PBARS, TK, FERATIOOPX, XAL_M1, XMGOPX, XFEOPX
    # THIS SUBROUTINE CALCULATES OPX ACTIVITIES WITH THE MODEL
    # IN TWQ202B - BA96a.SLN OF BERMAN. IT IS BASED ON THE MODEL OF ARANOVICH AND BERMAN (1997)
    # NOTE THAT THE ALOPX ACTIVITY MODEL INCLUDES A DARKEN CORRECTION OF THE FORM
    # RTLN(GAMMA)ALOPX =RTLN(GAMMA)ALOPX + FE/(FE+MG)*(DH-T*DS)
    W12 = -4543.8 + (TK * 3.36)
    W23 = -32213.3 - (PBARS * .69)
    W13 = -26944.5 - (PBARS * .58)
    RTGAMMAMGOPX = (W12 * (XFEOPX - (XFEOPX * XMGOPX))) - (W23 * XFEOPX * XAL_M1) + (W13 * (XAL_M1 - (XMGOPX * XAL_M1)))
    RTGAMMAFEOPX = (W12 * (XMGOPX - (XFEOPX * XMGOPX))) + (W23 * (XAL_M1 - (XFEOPX * XAL_M1))) - (W13 * XMGOPX * XAL_M1)
    RTGAMMAALOPX = (-1 * (W12 * XFEOPX * XMGOPX)) + (W23 * (XFEOPX - (XFEOPX * XAL_M1))) + (W13 * (XMGOPX - (XMGOPX * XAL_M1))) + (FERATIOOPX * ((24307) - (TK * 14.404) + (PBARS * .185)))
    GAMMAMGOPX = exp(RTGAMMAMGOPX / (8.314 * TK))
    GAMMAFEOPX = exp(RTGAMMAFEOPX / (8.314 * TK))
    GAMMAALOPX = exp( RTGAMMAALOPX / (8.314 * TK) )
    AEN = XMGOPX * GAMMAMGOPX
    AFS = XFEOPX * GAMMAFEOPX
    AALOPX = XAL_M1 * GAMMAALOPX
    GAMMAOPX = GAMMAMGOPX / GAMMAFEOPX

def outputfunc():
    # this function builds the output variables from the output of each iteration of the main program
    global TFEALIout, TFEALI, PFEALIout, PFEALI, TGAROPXIout, TGAROPXI, PGAROPXIout
    global PGAROPXI, TGARBTIout, TGARBTI, PGARBTIout, PGARBTI, TGARCRDIout, TGARCRDI
    global PGARCRDIout, PGARCRDI, TCout, TC, Pout, P, TGAROPXout, TGAROPX, TGARBTout, TGARBT, TGARCRDout, TGARCRD
    TFEALIout.append(TFEALI)
    PFEALIout.append(PFEALI)
    TGAROPXIout.append(TGAROPXI)
    PGAROPXIout.append(PGAROPXI)
    TGARBTIout.append(TGARBTI)
    PGARBTIout.append(PGARBTI)
    TGARCRDIout.append(TGARCRDI)
    PGARCRDIout.append(PGARCRDI)
    TCout.append(TC)
    Pout.append(P)
    TGAROPXout.append(TGAROPX)
    TGARBTout.append(TGARBT)
    TGARCRDout.append(TGARCRD)

########################################################
######## END DEFINING FUNCTIONS FOR THE PROGRAM ########
########################################################

########################################################
############### DEFINE THE MAIN PROGRAM ################
########################################################

def RCLCfunction():
    #SOME OF THESE MIGHT NOT BE USED. VESTIGES OF WIDNEY'S CODE THAT WAS MODIFIED TO MAKE THIS
    global XMNGAR, XFEGAR, XCAGAR, XMGGAR, XFEGARI, MGRATIOGA, MGRATIOGARI, MODXCAGAR
    global XFEOPX, XMGOPX, XAL_M1
    global FEGAR, MNGAR, MGGAR, CAGAR
    global SIOPX, TIOPX, ALOPX, CROPX, FE3OPX, FE2OPX, MNOPX, MGOPX, CAOPX
    global FECRD, MNCRD, MGCRD
    global CAPL, NAPL, KPL
    global SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT
    global skip_bt, skip_crd
    global DATASET, HALM, HPY, HGR, HAN, HBQ, HEN, HFS, HALOPX, HPHL, HANN, HCRD, HFECRD
    global SALM, SPY, SGR, SAN, SBQ, SEN, SFS, SALOPX, SPHL, SANN, SCRD, SFECRD
    global VALM, VPY, VGR, VAN, VBQ, VEN, VFS, VALOPX, VPHL, VANN, VCRD, VFECRD
    global AGR, APY, AAL, GAMMAGAR
    global MGRATIOCRD
    global TK, XFEBT, XALBT, XTIBT, XMGBT
    global TFEALI, PFEALI, TGAROPXI, PGAROPXI, TGARBTI, PGARBTI, TGARCRDI, PGARCRDI
    global TC, P, TGAROPX, TGARBT, TGARCRD
    global FERATIOOPX
    global AAN, XAB, XSAN, XAN
    global PBARS, TK, P

    # ORTHOPYROXENE mole graction calculations
    XFEOPXI = XFEOPX
    MGRATIOOPX = MGOPX / (MGOPX + FE2OPX)
    MGRATIOOPXI = MGRATIOOPX
    FERATIOOPX = 1 - MGRATIOOPX
    TOTOPX = XFEOPX + XMGOPX + XAL_M1+(TIOPX/2)+(FE3OPX/2)+(CROPX/2)+(MNOPX/2)+(CAOPX/2)

    # Garnet mole fraction calculations
    XMGGAR = MGGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XFEGAR = FEGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XCAGAR = CAGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XMNGAR = MNGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XFEGARI = XFEGAR
    MGRATIOGAR = MGGAR / (MGGAR + FEGAR)
    MGRATIOGARI = MGRATIOGAR
    MODXCAGAR = (CAGAR + MNGAR) / (CAGAR + MNGAR + FEGAR + MGGAR)

    # BIOTITE MOLE FRACTION CALCULATIONS
    if  (minmodes['bt'] < 0.01) or skip_bt:
            XFEBT, XMGBT, XTIBT, XALBT, MGRATIOBT = 0,0,0,0,0
    else:
            ALIVBT = 4.0 - SIBT
            ALVIBT = ALBT - ALIVBT
            XFEBT = FEBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XMGBT = MGBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XALBT = ALVIBT/ (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XTIBT = TIBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            MGRATIOBT = MGBT / (MGBT + FEBT)
    XFEBTI, MGRATIOBTI = XFEBT, MGRATIOBT

    # CORDIERITE MOLE FRACTION CALCULATIONS
    if  (minmodes['crd'] < 0.01) or skip_crd:
            XFECRD, XMGCRD, XMNCRD, MGRATIOCRD = 0,0,0,0
    else:
            XFECRD = FECRD / (FECRD + MGCRD + MNCRD)
            XMGCRD = MGCRD / (FECRD + MGCRD + MNCRD)
            XMNCRD = MNCRD / (FECRD + MGCRD + MNCRD)
            MGRATIOCRD = MGCRD / (MGCRD + FECRD)
    XFECRDI, MGRATIOCRDI = XFECRD, MGRATIOCRD

    # PLAGIOCLASE MOLE FRACTIONS
    XAN, XAB, XSAN = CAPL/(CAPL+NAPL+KPL), NAPL/(NAPL+CAPL+KPL), KPL/(NAPL+CAPL+KPL)

    #  VOLUME FRACTIONS OF FE-MG MINERALS FROM MODE
    VFGAR = minmodes['gar'] / (minmodes['gar'] + minmodes['opx'] + minmodes['crd'] + minmodes['bt'])
    VFOPX = minmodes['opx'] / (minmodes['gar'] + minmodes['opx'] + minmodes['crd'] + minmodes['bt'])
    VFCRD = minmodes['crd'] / (minmodes['gar'] + minmodes['opx'] + minmodes['crd'] + minmodes['bt'])
    VFBT =  minmodes['bt']  / (minmodes['gar'] + minmodes['opx'] + minmodes['crd'] + minmodes['bt'])

    # CONVERT VOLUME FRACTION MINERALS TO MOLE FRACTION MINERALS
    #  DENSITIES
    DENSGAR = (DENSFEGAR * XFEGAR) + (DENSMGGAR * XMGGAR) + (DENSCAGAR * XCAGAR) + (DENSMNGAR * XMNGAR)
    DENSOPX = (DENSFEOPX * (1 - MGRATIOOPX)) + (DENSMGOPX * MGRATIOOPX)
    DENSCRD = (DENSFECRD * (1 - MGRATIOCRD)) + (DENSMGCRD * MGRATIOCRD)
    DENSBT =  (DENSFEBT  * (1 - MGRATIOBT))  + (DENSMGBT  * MGRATIOBT)
    #  MOLECULAR WEIGHTS
    MWOPX = (SIOPX * 28.1) + (TIOPX*47.9) + (ALOPX * 26.1) + (CROPX*(52)) + (FE3OPX*55.8) + (FE2OPX * 55.8) + (MGOPX * 24.3) + (MNOPX * 54.9) + (CAOPX * 40.1) + (6 * 16)
    MWGAR = (3.00  * 28.1) + (2.00  * 26.1) + (FEGAR * 55.8) + (MGGAR * 24.3) + (MNGAR * 54.9) + (CAGAR * 40.1) + (12 * 16)
    if (minmodes['crd'] > 0.01) and not skip_crd:
        MWCRD = (5.00  * 28.1) + (4.00  * 26.1) + (FECRD * 55.8) + (MGCRD * 24.3) + (MNCRD * 54.9) + (18 * 16)
    if (minmodes['bt']>0.01) and not skip_bt:
        MWBT =  (SIBT * 28.1)  + (TIBT * 47.9)  + (ALBT * 26.1)  + (FEBT * 55.8)  + (MNBT * 54.9)  + (MGBT * 24.3)  + (NABT * 23) + (KBT * 39.1) + (11 * 16) + 2
    #   MOLES OF MINERALS
    MOLEGAR = (VFGAR * DENSGAR) / MWGAR
    MOLEOPX = (VFOPX * DENSOPX) / MWOPX
    if (minmodes['crd'] < 0.01) or skip_crd:
        MOLECRD = 0
    else:
        MOLECRD = (VFCRD * DENSCRD) / MWCRD
    if (minmodes['bt'] < 0.01) or skip_bt:
        MOLEBT = 0
    else:
        MOLEBT = (VFBT * DENSBT) / MWBT
    # MOLES OF FE-MG COMPONENTS OF MINERALS
    MOLEFEMGGAR = MOLEGAR * (FEGAR + MGGAR)
    MOLEFEMGOPX = MOLEOPX * (FE2OPX + MGOPX)
    if (minmodes['crd'] > 0.01) and not skip_crd:
        MOLEFEMGCRD = MOLECRD * (FECRD + MGCRD)
    else:
        MOLEFEMGCRD = 0
    if (minmodes['bt'] > 0.01) and not skip_bt:
        MOLEFEMGBT =  MOLEBT  * (FEBT  + MGBT )
    else:
        MOLEFEMGBT = 0
    # MOLE FRACTION OF FE-MG COMPONENTS OF MINERALS
    MFGAR = MOLEFEMGGAR / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFOPX = MOLEFEMGOPX / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFCRD = MOLEFEMGCRD / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFBT =  MOLEFEMGBT  / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    #  CALCULATE XMG ROCK (THIS APPEARS TO HAVE BEEN USED AS A TEST IN WIDNEY'S ORIGINAL CODE. NOT NEEDED ANYMORE BUT HAVEN'T DELETED YET WHILE I CHECK THE REST OF THE CODE)
    XMGROCK = (MGRATIOGAR * MFGAR) + (MGRATIOOPX * MFOPX) + (MGRATIOCRD * MFCRD) + (MGRATIOBT * MFBT)

    #  CALCULATE GRT-OPX FE-MG  -  GRT-OPX-PL-QTZ (FE-END MEMBER)INTERSECTION
    TK, P, PBARS = 1123.85, 6, 6000 #INITIAL GUESSES 850 C and 6 kbar
    CP() #CALCULATE H AND S AT STARTING GUESSES
    for J in range(10): # SHOULD CONVERGE IN < 10 ITERATIONS
        # CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE
            # the following comment was left in by Widney, although the problem seems to have been fixed, whatever it was
            # FIXME: Causing an error on numbers,
        GARNET()
        PLAGIOCLASE()
        ORTHOPYROXENE()
        VOLUMEPT()
        DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
        DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
        DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
        KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
        P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (log(KGAPES)))) / DELTAVGAPES
        PBARS = P * 1000
        # CALCULATE GRT-OPX FE-MG EXCHANGE TEMP AT THIS PRESSURE
        DELTAHFEMGOPX = (((1 * HEN) + ((1 / 3) * HALM)) - ((1 * HFS) + ((1 / 3) * HPY))) / 1000
        DELTASFEMGOPX = (((1 * SEN) + ((1 / 3) * SALM)) - ((1 * SFS) + ((1 / 3) * SPY))) / 1000
        DELTAVFEMGOPX = ((1 * VEN) + ((1 / 3) * VALM)) - ((1 * VFS) + ((1 / 3) * VPY))
        GAMMAFEMGOPX = GAMMAGAR * GAMMAOPX
        KDGAROPX = (XFEGAR * XMGOPX) / (XMGGAR * XFEOPX)
        TGAROPX = (DELTAHFEMGOPX + (P * DELTAVFEMGOPX)) / (DELTASFEMGOPX - (.008314 * log(KDGAROPX)) - (.008314 * log(GAMMAFEMGOPX)))
        TK = TGAROPX
        TCGAROPX = TGAROPX - 273
        CP() # UPDATE H AND S AND REITERATE
    TGAROPXI = TCGAROPX
    PGAROPXI = P

    #  CALCULATE GRT-CRD FE-MG  -  GRT-OPX-PL-QTZ (FE END MEMBER) INTERSECTION IF CORDIERITE IS BEING CONSIDERED
    if (minmodes['crd'] > 0.01) and not skip_crd:
        TK, P, PBARS = 1123.85, 6, 6000 #INITIAL GUESSES 850 C and 6 kbar
        CP() #CALCULATE H AND S AT STARTING GUESSES
        for J in range(10): # SHOULD CONVERGE IN < 10 ITERATIONS
            # CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE
            GARNET()
            PLAGIOCLASE()
            ORTHOPYROXENE()
            VOLUMEPT()
            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000
            # CALCULATE GRT-CRD FE-MG EXCHANGE TEMP AT THIS PRESSURE
            CORDIERITE()
            DELTAHFEMGCRD = (((.5 * HCRD) + ((1 / 3) * HALM)) - ((.5 * HFECRD) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGCRD = (((.5 * SCRD) + ((1 / 3) * SALM)) - ((.5 * SFECRD) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGCRD = ((.5 * VCRD) + ((1 / 3) * VALM)) - ((.5 * VFECRD) + ((1 / 3) * VPY))
            GAMMAFEMGCRD = GAMMAGAR * GAMMACRD
            KDGARCRD = (XFEGAR * XMGCRD) / (XMGGAR * XFECRD)
            TKGARCRD = (DELTAHFEMGCRD + (P * DELTAVFEMGCRD)) / (DELTASFEMGCRD - (.008314 * log(KDGARCRD)) - (.008314 * log(GAMMAFEMGCRD)))
            TK = TKGARCRD
            TGARCRD = TKGARCRD - 273
            CP() # UPDATE H AND S AND REITERATE
        TGARCRDI = TGARCRD
        PGARCRDI = P
    else:
        TGARCRD = 0
        TGARCRDI = 0
        PGARCRD = 0
        PGARCRDI = 0

    #  CALCULATE GRT-BT FE-MG  -  GRT-OPX-PL-QTZ (FE END MEMBER) INTERSECTION IF BIOTITE IS BEING CONSIDERED
    if (minmodes['bt'] > 0.01) and not skip_bt:
        TK, P, PBARS = 1123.85, 600, 6000 #INITIAL GUESSES 850 C and 6 kbar
        CP() #CALCULATE H AND S AT STARTING GUESSES
        for J in range (10): #SHOULD CONVERGE IN LESS THAN 10 ITERATIONS
            # CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE
            GARNET()
            PLAGIOCLASE()
            ORTHOPYROXENE()
            VOLUMEPT()
            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000
            # CALCULATE GRT-BT FE-MG EXCHANGE TEMP AT THIS PRESSURE
            BIOTITE()
            DELTAHFEMGBT = ((((1 / 3) * HPHL) + ((1 / 3) * HALM)) - (((1 / 3) * HANN) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGBT = ((((1 / 3) * SPHL) + ((1 / 3) * SALM)) - (((1 / 3) * SANN) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGBT = (((1 / 3) * VPHL) + ((1 / 3) * VALM)) - (((1 / 3) * VANN) + ((1 / 3) * VPY))
            GAMMAGARBT = GAMMABT * GAMMAGAR
            KDGARBT = (XFEGAR * XMGBT) / (XMGGAR * XFEBT)
            TKGARBT = (DELTAHFEMGBT + (P * DELTAVFEMGBT)) / (DELTASFEMGBT - (.008314 * log(KDGARBT)) - (.008314 * log(GAMMAGARBT)))
            TK = TKGARBT
            TGARBT = TKGARBT - 273
            CP() # UPDATE H AND S AND REITERATE
        TGARBTI = TGARBT
        PGARBTI = P
    else:
        TGARBT = 0
        TGARBTI = 0
        PGARBT = 0
        PGARBTI = 0

    # CALCULATE CONVERGED INTERSECTION OF GRT-OPX AL-SOLUBILITY AND GRT-OPX-PL-QTZ USING
    # FE-END MEMBER EXPRESSIONS.
    # CONVERGENCE APPROACH - 1. CALCULATE INITIAL INTERSECTION OF GRT-OPX AL-SOLUB AND GRT-OPX-PL-QTZ.
    # 2. CHANGE KD GRT-OPX (AND if  APPLICABLE KD GRT-CRD AND KD GRT-BT) SO COINCIDES WITH 1.
    # 3. ADJUST FE/MG RATIOS OF FE-MG MINERALS TO SATISFY KD'S.
    # 4. REPEAT 10 TIMES (I = 1 TO 10) TO GET CONVERGENCE.
    # ASSUME INITIAL TEMP TO BEGIN
    TK, P, PBARS = 1123.85, 600, 6000 #INITIAL GUESSES 850 C and 6 kbar

    for I in range(10): #should converge in <10 iterations
        CP() # calculate H and S for starting PT guesses
        # CALCULATE INTERSECTION OF FE-AL-OPX AND GRT-OPX-PL-QTZ IN 10 ITERATIONS (J = 1 TO 10)
        for J in range (10): #should converge in < 10 iterations
            # CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE
            GARNET()
            PLAGIOCLASE()
            ORTHOPYROXENE()
            VOLUMEPT()
            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000
            # CALCULATE FE-ALOPX TEMPERATURE AT THIS PRESSURE
            DELTAHFEAL = ((HALOPX + (3 * HFS)) - HALM) / 1000
            DELTASFEAL = ((SALOPX + (3 * SFS)) - SALM) / 1000
            DELTAVFEAL = (VALOPX + (3 * VFS)) - VALM
            KFEAL = ((AFS ** 3) * AALOPX) / AAL
            TFEAL = (DELTAHFEAL + (P * DELTAVFEAL)) / (DELTASFEAL - (.008314 * log(KFEAL)))
            TK = TFEAL
            TC = TK - 273
            CP() #update H and S for next iteration
        if  (I==0): #for the first iteration, define the initial T and P calculated
            TFEALI = TC
            PFEALI = P

        GARNET()
        ORTHOPYROXENE()
        if  (minmodes['crd'] > 0.01) and not skip_crd:
            CORDIERITE()
        if  (minmodes['bt'] > 0.01) and not skip_bt:
            BIOTITE()
        VOLUMEPT()

        # CALCULATES A CORRECTED KD(GRT-OPX(FE-MG))
        DELTAHFEMGOPX = (((1 * HEN) + ((1 / 3) * HALM)) - ((1 * HFS) + ((1 / 3) * HPY))) / 1000
        DELTASFEMGOPX = (((1 * SEN) + ((1 / 3) * SALM)) - ((1 * SFS) + ((1 / 3) * SPY))) / 1000
        DELTAVFEMGOPX = ((1 * VEN) + ((1 / 3) * VALM)) - ((1 * VFS) + ((1 / 3) * VPY))
        GAMMAFEMGOPX = GAMMAGAR * GAMMAOPX
        KDGAROPX = ((TK * DELTASFEMGOPX) - DELTAHFEMGOPX - (P * DELTAVFEMGOPX) - (.008314 * TK * (log(GAMMAFEMGOPX)))) / (.008314 * TK)
        KDGAROPX = exp(KDGAROPX)
        # CALCULATES A CORRECTED KD(GRT-CRD) if cordierite is  being considered
        if  (minmodes['crd'] > 0.01) and not skip_crd:
            DELTAHFEMGCRD = (((.5 * HCRD) + ((1 / 3) * HALM)) - ((.5 * HFECRD) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGCRD = (((.5 * SCRD) + ((1 / 3) * SALM)) - ((.5 * SFECRD) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGCRD = ((.5 * VCRD) + ((1 / 3) * VALM)) - ((.5 * VFECRD) + ((1 / 3) * VPY))
            GAMMAFEMGCRD = GAMMAGAR * GAMMACRD
            KDGARCRD = ((TK * DELTASFEMGCRD) - DELTAHFEMGCRD - (P * DELTAVFEMGCRD) - (.008314 * TK * (log(GAMMAFEMGCRD)))) / (.008314 * TK)
            KDGARCRD = exp(KDGARCRD)
        # CALCULATES A CORRECTED KD(GRT-BT) if biotite is being considered
        if  (minmodes['bt'] > 0.01) and not skip_bt:
            DELTAHFEMGBT = ((((1 / 3) * HPHL) + ((1 / 3) * HALM)) - (((1 / 3) * HANN) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGBT = ((((1 / 3) * SPHL) + ((1 / 3) * SALM)) - (((1 / 3) * SANN) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGBT = (((1 / 3) * VPHL) + ((1 / 3) * VALM)) - (((1 / 3) * VANN) + ((1 / 3) * VPY))
            GAMMAGARBT = GAMMABT * GAMMAGAR
            KDGARBT = ((TK * DELTASFEMGBT) - DELTAHFEMGBT - (P * DELTAVFEMGBT) - (.008314 * TK * (log(GAMMAGARBT)))) / (.008314 * TK)
            KDGARBT = exp(KDGARBT)

        # QUADRATIC SOLUTION TO CORRECTED MG-RATIOS OF MINERALS
        for L in range(10): #should converge in < 10 iterations
            # QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND OPX
            A = MFOPX - (KDGAROPX * MFOPX)
            B = (MFOPX * KDGAROPX) + MFGAR + (XMGROCK * KDGAROPX) - XMGROCK - (KDGAROPX * MGRATIOCRD * MFCRD) - (KDGAROPX * MGRATIOBT * MFBT) + (MFCRD * MGRATIOCRD) + (MFBT * MGRATIOBT)
            C = (MGRATIOCRD * MFCRD * KDGAROPX) + (MGRATIOBT * MFBT * KDGAROPX) - (XMGROCK * KDGAROPX)
            MGRATIOOPX = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
            MGRATIOGAR = (XMGROCK - (MGRATIOOPX * MFOPX) - (MGRATIOCRD * MFCRD) - (MGRATIOBT * MFBT)) / MFGAR
            # QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND BIOTITE if biotite is being considered
            if  (minmodes['bt'] > 0.01) and not skip_bt:
                A = MFBT - (KDGARBT * MFBT)
                B = (MFBT * KDGARBT) + MFGAR + (XMGROCK * KDGARBT) - XMGROCK - (KDGARBT * MGRATIOOPX * MFOPX) - (KDGARBT * MGRATIOCRD * MFCRD) + (MFOPX * MGRATIOOPX) + (MFCRD * MGRATIOCRD)
                C = (MGRATIOOPX * MFOPX * KDGARBT) + (MGRATIOCRD * MFCRD * KDGARBT) - (XMGROCK * KDGARBT)
                MGRATIOBT = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
                MGRATIOGAR = (XMGROCK - (MGRATIOBT * MFBT) - (MGRATIOOPX * MFOPX) - (MGRATIOCRD * MFCRD)) / MFGAR
            # QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND CORDIERITE if cordierite is being considered
            if  (minmodes['crd'] > 0.01) and not skip_crd:
                A = MFCRD - (KDGARCRD * MFCRD)
                B = (MFCRD * KDGARCRD) + MFGAR + (XMGROCK * KDGARCRD) - XMGROCK - (KDGARCRD * MGRATIOOPX * MFOPX) - (KDGARCRD * MGRATIOBT * MFBT) + (MFOPX * MGRATIOOPX) + (MFBT * MGRATIOBT)
                C = (MGRATIOOPX * MFOPX * KDGARCRD) + (MGRATIOBT * MFBT * KDGARCRD) - (XMGROCK * KDGARCRD)
                MGRATIOCRD = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
                MGRATIOGAR = (XMGROCK - (MGRATIOCRD * MFCRD) - (MGRATIOOPX * MFOPX) - (MGRATIOBT * MFBT)) / MFGAR

        FERATIOOPX = 1 - MGRATIOOPX
        XMGOPX = (MGRATIOOPX) * ((FE2OPX + MGOPX) / 2)
        XFEOPX = (1 - MGRATIOOPX) * ((FE2OPX + MGOPX) / 2)
        XMGGAR = (MGRATIOGAR) * ((FEGAR + MGGAR) / (FEGAR + MGGAR + CAGAR + MNGAR))
        XFEGAR = (1 - MGRATIOGAR) * ((FEGAR + MGGAR) / (FEGAR + MGGAR + CAGAR + MNGAR))
        XMGCRD = MGRATIOCRD * (XFECRD + XMGCRD)
        XFECRD = (1 - MGRATIOCRD) * (XFECRD + XMGCRD)
        XMGBT =  MGRATIOBT * (XFEBT + XMGBT)
        XFEBT =  (1 - MGRATIOBT) * (XFEBT + XMGBT)
    #NEXT I

    # CALCULATE GRT-OPX FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T
    KDGAROPX = (XFEGAR * XMGOPX) / (XMGGAR * XFEOPX)
    TGAROPX = (DELTAHFEMGOPX + (P * DELTAVFEMGOPX)) / (DELTASFEMGOPX - (.008314 * log(KDGAROPX)) - (.008314 * log(GAMMAFEMGOPX)))-273
    # CALCULATE GRT-CRD FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T, if cordierite is being considered
    if  (minmodes['crd'] > 0.01) and not skip_crd:
        KDGARCRD = (XFEGAR * XMGCRD) / (XMGGAR * XFECRD)
        TGARCRD = (DELTAHFEMGCRD + (P * DELTAVFEMGCRD)) / (DELTASFEMGCRD - (.008314 * log(KDGARCRD)) - (.008314 * log(GAMMAFEMGCRD)))-273
    # CALCULATE GRT-BT FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T, if biotite is being considered
    if  (minmodes['bt'] > 0.01) and not skip_bt:
        KDGARBT = (XFEGAR * XMGBT) / (XMGGAR * XFEBT)
        TGARBT = (DELTAHFEMGBT + (P * DELTAVFEMGBT)) / (DELTASFEMGBT - (.008314 * log(KDGARBT)) - (.008314 * log(GAMMAGARBT)))-273
    # CHECK if  RECALCULATED XMGROCK IS THE SAME AS THE INITIAL XMGROCKI (this TEST is not currently being used, but has not been deleted while I look through the rest of the code)
    TEST = (MGRATIOGAR * MFGAR) + (MGRATIOOPX * MFOPX) + (MGRATIOCRD * MFCRD) + (MGRATIOBT * MFBT)

########################################################
########## END DEFINING THE MAIN PROGRAM ###############
########################################################

########################################################
######### IMPORTING COMPOSITIONAL DATA & MODES #########
########################################################

#import opx formula normalized to 6 oxygen
aSIOPX, aTIOPX, aALOPX, aCROPX, aFE3OPX, aFE2OPX, aMNOPX, aMGOPX, aCAOPX = [],[],[],[],[],[],[],[],[]
with open('opx.txt') as opxdata:
    for line in opxdata:
        line = line.rstrip()
        if rsearch('^[A-z]', line) == None:
            line = line.split()
            aSIOPX.append(float(line[0]))
            aTIOPX.append(float(line[1]))
            aALOPX.append(float(line[2]))
            aCROPX.append(float(line[3]))
            aFE3OPX.append(float(line[4]))
            aFE2OPX.append(float(line[5]))
            aMNOPX.append(float(line[6]))
            aMGOPX.append(float(line[7]))
            aCAOPX.append(float(line[8]))
#converted to numpy arrays for ease of use in ALOPX FUNCTIONS
aSIOPX, aTIOPX, aALOPX = nparray(aSIOPX), nparray(aTIOPX), nparray(aALOPX)
aCROPX, aFE3OPX, aFE2OPX = nparray(aCROPX), nparray(aFE3OPX), nparray(aFE2OPX)
aMNOPX, aMGOPX, aCAOPX = nparray(aMNOPX), nparray(aMGOPX), nparray(aCAOPX)

#import garnet formula normalized to 12 oxygen
aFEGAR, aMNGAR, aMGGAR, aCAGAR = [],[],[],[]
with open('gar.txt') as gardata:
    for line in gardata:
        line = line.rstrip()
        if rsearch('^[A-z]', line) == None:
            line = line.split()
            aFEGAR.append(float(line[5]))
            aMNGAR.append(float(line[6]))
            aMGGAR.append(float(line[7]))
            aCAGAR.append(float(line[8]))
#converted to numpy arrays for consistency with OPX data
aFEGAR, aMNGAR, aMGGAR, aCAGAR = nparray(aFEGAR), nparray(aMNGAR), nparray(aMGGAR), nparray(aCAGAR)

#import plagioclase formula normalized to 8 oxygen
aCAPL, aNAPL, aKPL = [],[],[]
with open('pl.txt') as pldata:
    for line in pldata:
        line = line.rstrip()
        if rsearch('^[A-z]', line) == None:
            line = line.split()
            aCAPL.append(float(line[8]))
            aNAPL.append(float(line[9]))
            aKPL.append(float(line[10]))
#converted to numpy arrays for consistency with OPX data
aCAPL, aNAPL, aKPL = nparray(aCAPL), nparray(aNAPL), nparray(aKPL)

#import cordierite data normalized to 18 oxygen
aFECRD, aMNCRD, aMGCRD = [],[],[]
try:
    with open('crd.txt') as crddata:
        for line in crddata:
            line = line.rstrip()
            if rsearch('^[A-z]', line) == None:
                line = line.split()
                aFECRD.append(float(line[5]))
                aMNCRD.append(float(line[6]))
                aMGCRD.append(float(line[7]))
    #converted to numpy arrays for consistency with OPX data
    aFECRD, aMNCRD, aMGCRD = nparray(aFECRD), nparray(aMNCRD), nparray(aMGCRD)
    skip_crd = False
except:
    print('no cordierite compositional file found in directory\ncordierite will not be considered in the calculation\n')
    skip_crd = True

#import biotite formula normalized to 11 oxygen (10 O + 2 OH), optional
aSIBT, aTIBT, aALBT, aFEBT, aMNBT, aMGBT, aNABT, aKBT = [],[],[],[],[],[],[],[]
try:
    with open('bt.txt') as btdata:
        for line in btdata:
            line = line.rstrip()
            if rsearch('^[A-z]', line) == None:
                line = line.split()
                aSIBT.append(float(line[0]))
                aTIBT.append(float(line[1]))
                aALBT.append(float(line[2]))
                aFEBT.append(float(line[5]))
                aMNBT.append(float(line[6]))
                aMGBT.append(float(line[7]))
                aNABT.append(float(line[9]))
                aKBT.append(float(line[10]))
    #converted to numpy arrays for consistency with OPX data
    aSIBT, aALBT, aFEBT = nparray(aSIBT), nparray(aALBT), nparray(aFEBT)
    aMNBT, aMGBT, aNABT, aKBT = nparray(aMNBT), nparray(aMGBT), nparray(aNABT), nparray(aKBT)
    skip_bt = False
except:
    print('no biotite compositional file found in directory\nbiotite will not be considered in the calculation\n')
    skip_bt = True

#import mineral modes
minmodes = dict()
with open('modes.txt') as modes:
    for line in modes:
        line = line.rstrip().split()
        minmodes[line[0]]=float(line[1])

########################################################
####### END IMPORTING COMPOSITIONAL DATA & MODES #######
########################################################

########################################################
############# CHOOSE AL IN OPX SITE MODEL ##############
########################################################

print ("YOU HAVE A CHOICE FOR CALCULATING XALM IN OPX.\nTHE FOLLOWING FORMULAE ASSUME A 6-OXYGEN OPX FORMULA.")
print ("1: XAL_M1 = Al - (2 - Si)")
print ("2: XAL_M1 = Al/2")
print ("3: XAL_M1 = (Al/2) / (Fe2+ + Mg + Mn + Ca + (Al/2) )")
print ("4: XAL_M1 = (Al - Fe3+ - Cr - (2*Ti) ) / 2 \n")
num = int( input("Please enter 1,2,3,4: "))
if  num == 1:
    ALOPX1()
if  num == 2:
    ALOPX2()
if  num == 3:
    ALOPX3()
if  num == 4:
    ALOPX4()

########################################################
######### END CHOOSING AL IN OPX SITE MODEL ############
########################################################

########################################################
################# Thermodynamic data ###################
########################################################

DATASET = [[-5265317.1, 341.5824, 621.4269, -3287.931, -15081040, 2211865100], \
[-6284733.7, 268.8, 590.9042, -2826.956, -13320810, 1260328500], \
[-6632861, 255.15, 573.43042, -2039.405, -18887168, 2319311872], \
[-4228730, 200.1861, 439.36938, -3734.149, 0, -317023232], \
[-908626.8, 44.2068, 80.01199, -240.276, -3546684, 491568384], \
[-1546037.1, 66.18, 166.5795, -1200.588, -2270560, 279150300], \
[-1192860, 96.5587, 174.2024, -1392.959, -454390, -37711400], \
[ -1631665.9, 35.375, 119.38, 774.808, -6509130, 422877600], \
[-6216676.7, 325.9239, 610.37988, -2083.781, -21533008, 2841040896], \
[-5155234.4, 405.01, 727.208, -4775.04, -13831900, 2119060000], \
[-9161425.7, 416.2714, 954.3865, -7962.274, -2317258, -370214090], \
[-8429860.2, 482.8282, 983.479, -8403.659, -1870290, -85683500]]
DENSFEGAR, DENSMGGAR, DENSCAGAR, DENSMNGAR, DENSFEOPX = 4.33, 3.54, 3.56, 4.19, 3.96
DENSMGOPX, DENSMGCRD, DENSFECRD, DENSFEBT, DENSMGBT = 3.21, 2.53, 2.78, 3.3, 2.7

########################################################
############## end Thermodynamic data ##################
########################################################

########################################################
##### run calcs for various compositional combos #######
########################################################
# this section works, but is not particularly elegant
# redundant code could be consolidated if desired

#output variables
TFEALIout = ['Fe-Al T init']
PFEALIout = ['Fe-Al P init']
TGAROPXIout = ['gar-opx Fe-Mg T init']
PGAROPXIout = ['gar-opx Fe-Mg P init']
TGARBTIout = ['gar-bt Fe-Mg T init']
PGARBTIout = ['gar-bt Fe-Mg P init']
TGARCRDIout = ['gar-crd Fe-Mg T init']
PGARCRDIout = ['gar-crd Fe-Mg P init']
TCout = ['Fe-Al T final']
Pout = ['Fe-Al P final']
TGAROPXout = ['gar-opx Fe-Mg T final']
TGARBTout = ['gar-bt Fe-Mg T final']
TGARCRDout = ['gar-crd Fe-Mg T init']
calctracker = ['analyses used'] #will be used to track which mineral combos were used for each calculation

# Determine run mode from user. either
# run calculations in sequence (gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN)
# or run every possible combination of the input mineral analyses
if len(aMGGAR) == len(aMGOPX) == len(aCAPL):
    if not skip_bt and not skip_crd:
        if len(aMGBT) == len(aMGGAR) == len(aMGCRD):
            print('You entered an equal number of GAR, OPX, PL, CRD, and BT analyses.\nWould you like to:')
            print('1: Run them in sequence (gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN) or')
            print('2: Run every possible combination of the input mineral analyses?\n')
            runmode = int(input('Enter 1 or 2: '))
        else:
            print('You did not enter an equal number of analyses for each mineral.')
            print('Runmode 2 will be used: run every possible combination of the input mineral analyses.\n')
            print('If you would rather use runmode 1 (run in sequence: gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN),\n include the same number of analyses for each mineral in your input files.')
            runmode = 2
    elif not skip_bt:
        if len(aMGBT) == len(aMGGAR):
            print('You entered an equal number of GAR, OPX, PL, and BT analyses.\nWould you like to:')
            print('1: Run them in sequence (gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN) or')
            print('2: Run every possible combination of the input mineral analyses?\n')
            runmode = int(input('Enter 1 or 2: '))
        else:
            print('You did not enter an equal number of analyses for each mineral.')
            print('Runmode 2 will be used: run every possible combination of the input mineral analyses.\n')
            print('If you would rather use runmode 1 (run in sequence: gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN),\n include the same number of analyses for each mineral in your input files.')
            runmode = 2
    elif not skip_crd:
        if len(aMGCRD) == len(aMGGAR):
            print('You entered an equal number of GAR, OPX, PL, and CRD analyses.\nWould you like to:')
            print('1: Run them in sequence (gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN) or')
            print('2: Run every possible combination of the input mineral analyses?\n')
            runmode = int(input('Enter 1 or 2: '))
        else:
            print('You did not enter an equal number of analyses for each mineral.')
            print('Runmode 2 will be used: run every possible combination of the input mineral analyses.\n')
            print('If you would rather use runmode 1 (run in sequence: gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN),\n include the same number of analyses for each mineral in your input files.')
            runmode = 2
    else:
        print('You entered an equal number of GAR, OPX, and PL analyses.\nWould you like to:')
        print('1: Run them in sequence (gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN) or')
        print('2: Run every possible combination of the input mineral analyses?\n')
        runmode = int(input('Enter 1 or 2: '))
else:
    print('You did not enter an equal number of analyses for each mineral.')
    print('Runmode 2 will be used: run every possible combination of the input mineral analyses.\n')
    print('If you would rather use runmode 1 (run in sequence: gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN),\n include the same number of analyses for each mineral in your input files.')
    runmode = 2

# run input mineral data in sequence if user selected runmode 1: gar1-opx1-pl1, gar2-opx2-pl2... garN-opxN-plN
if runmode == 1:
    for i in range(len(aSIOPX)): #iterate through mineral compositions
        SIOPX, TIOPX, ALOPX, CROPX, FE3OPX, FE2OPX, MNOPX, MGOPX, CAOPX = aSIOPX[i], aTIOPX[i], aALOPX[i], aCROPX[i], aFE3OPX[i], aFE2OPX[i], aMNOPX[i], aMGOPX[i], aCAOPX[i]
        XFEOPX, XMGOPX, XAL_M1 = aXFEOPX[i], aXMGOPX[i], aXAL_M1[i]
        FEGAR, MNGAR, MGGAR, CAGAR = aFEGAR[i], aMNGAR[i], aMGGAR[i], aCAGAR[i]
        CAPL, NAPL, KPL = aCAPL[i], aNAPL[i], aKPL[i]
        if not skip_bt: # if using biotite
            SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT = aSIBT[i], aTIBT[i], aALBT[i], aFEBT[i], aMNBT[i], aMGBT[i], aNABT[i], aKBT[i]
        if not skip_crd: # if using cordierite
            FECRD, MNCRD, MGCRD = aFECRD[i], aMNCRD[i], aMGCRD[i]
        RCLCfunction()
        outputfunc()
        calctracker.append('calculation'+str(i+1))
# Run every possible combination of input mineral analyses if use selected runmode 2
elif runmode == 2:
    for i in range(len(aSIOPX)): #iterate through each opx
        SIOPX, TIOPX, ALOPX, CROPX, FE3OPX, FE2OPX, MNOPX, MGOPX, CAOPX = aSIOPX[i], aTIOPX[i], aALOPX[i], aCROPX[i], aFE3OPX[i], aFE2OPX[i], aMNOPX[i], aMGOPX[i], aCAOPX[i]
        XFEOPX, XMGOPX, XAL_M1 = aXFEOPX[i], aXMGOPX[i], aXAL_M1[i]
        for ii in range(len(aFEGAR)): #iterate through each garnet
            FEGAR, MNGAR, MGGAR, CAGAR = aFEGAR[ii], aMNGAR[ii], aMGGAR[ii], aCAGAR[ii]
            for iii in range(len(aCAPL)): #iterate through each plagioclase cations per 8 oxgyen
                CAPL, NAPL, KPL = aCAPL[iii], aNAPL[iii], aKPL[iii]
                if not skip_bt and not skip_crd: # if including biotite and cordierite
                    for iv in range(len(aFECRD)): #iterate through each cordierite
                        FECRD, MNCRD, MGCRD = aFECRD[iv], aMNCRD[iv], aMGCRD[iv]
                        for v in range(len(aSIBT)): #biotite cations per 11 oxygen (10 O + 2(OH))
                            SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT = aSIBT[v], aTIBT[v], aALBT[v], aFEBT[v], aMNBT[v], aMGBT[v], aNABT[v], aKBT[v]
                            RCLCfunction()
                            outputfunc()
                            calctracker.append('opx'+str(i+1)+' gar'+str(ii+1)+' pl'+str(iii+1)+' bt'+str(v+1)+' crd'+str(iv+1))
                elif not skip_bt: # if including biotite, but not cordierite
                    for v in range(len(aSIBT)): #iterate through each biotite composition
                        SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT = aSIBT[v], aTIBT[v], aALBT[v], aFEBT[v], aMNBT[v], aMGBT[v], aNABT[v], aKBT[v]
                        RCLCfunction()
                        outputfunc()
                        calctracker.append('opx'+str(i+1)+' gar'+str(ii+1)+' pl'+str(iii+1)+' bt'+str(v+1))
                elif not skip_crd: # if including cordierite but not biotite
                    for iv in range(len(aFECRD)): # iterate through each cordierite composition
                        FECRD, MNCRD, MGCRD = aFECRD[iv], aMNCRD[iv], aMGCRD[iv]
                        RCLCfunction()
                        outputfunc()
                        calctracker.append('opx'+str(i+1)+' gar'+str(ii+1)+' pl'+str(iii+1)+' crd'+str(iv+1))
                else: # if not using biotite or cordierite
                    RCLCfunction()
                    outputfunc()
                    calctracker.append('opx'+str(i+1)+' gar'+str(ii+1)+' pl'+str(iii+1))
print('\ndone with calculations\n')
########################################################
## Done running calcs for various compositional combos #
########################################################

########################################################
################ outputting results ####################
########################################################

results = [calctracker,TCout,Pout,TGAROPXout,TGARBTout,TGARCRDout,TFEALIout,PFEALIout,TGAROPXIout,PGAROPXIout,TGARBTIout,PGARBTIout,TGARCRDIout,PGARCRDIout]
with open('outputfile.csv', 'w', newline='') as f:
    w = csvwriter(f)
    w.writerows(results)
print('calculation results saved to outputfile.csv\n')
########################################################
############## Done outputting results #################
########################################################
