# batchRCLC
python code, input files and instructions for running RCLC Al-in-opx thermobarometry

##### Background on the develompent of this code ######
    The version of the RCLC algorithm used in this code is a modified version of Tom Chacko's 1996
    original. The modifications have been made by Chris McFarlane (04-97), David Pattison (05-97 and 10-01),
    and Tom Chacko (05-02).
    
    This python code was written by Robert Holder in May 2021 as a modfication and expansion of
    python code and a web application (https://rclc.facsci.ualberta.ca/) developed by
    Justin Widney at the University of Alberta, as part of a summer research internship
    under the supervision of Tom Chacko in 2018. The reason Holder developed this code
    was to run to be able to run large batches of RCLC calculations with a single file input
    (a list of mineral compositions).
    
    If you use this code, please cite:
        The original development of this idea:
            Chacko T, Lamb M, Farquhar J (1996) Ultra-high temperature metamorphism in the Kerala Khondalite Belt. In: Santosh M &
              Yoshida M (eds) Gondwana Research Group Memoir 3: The Archean and Proterozoic Terrains in Southern India within
              Gondwana, pp. 157–165. Osaka, Japan: Field Research Publishers.
        Expansion of the idea and the BASIC version of the algorithm:
            Pattison DRM, Chacko T, Farquhar J, McFarlane CRM (2003). Temperatures of granulite facies metamorphism: constraints
              from experimental phase equilibria and thermobarometry corrected for late Fe–Mg exchange. Journal of Petrology 44, 867–900.
        The first python version of the algorithm and the accompanying web application:
            Pattison DRM, Chacko T, Farquhar J, McFarlane CRM, Widney J (2019) Program “RCLC”: Garnet–Orthopyroxene Thermobarometry
              Corrected for Late Fe–Mg Exchange. Journal of Petrology 60, 1107–1108. https://doi.org/10.1093/petrology/egz018
        This code (allowing for batch processing of large quantities of data):
            Horton F, Holder RM, Swindle C (2021) An extensive record of orogenesis recorded in a Madagascar granulite. Journal
            of Metamorphic Geology. [In revision at the time this document was written.]
    
##### What this code does #####
    RCLC calculates P and T of metamorphism from garnet-othropyroxene-plagioclase-quartz
    equilibria (the GOPS barometer and Al-in-opx thermometer), corrected for Fe-Mg exchange
    among minerals (garnet and orthopyroxene as well as cordierite and biotite if they
    are relevant for a particular sample). In other words, this is Al-in-orthopyroxene
    thermobarometry, corrected for retrograde Fe-Mg diffusion.
    
    The Fe-Mg exchange correction is done by adjusting the Fe-Mg ratios of the mafic
    minerals according to their modes (mass balance) and Fe-Mg exchange equilibria:
    garnet-orthopyroxene, garnet-cordierite, and garnet-biotite
    
    The Thermodynamic data (end-members and a-x relationships for solid solutions) are
    from TWQ202b (files BA96A.DAT and BA96A.SLN), except for BT, which uses the solution
    model of McMullin (1991 CanMin). For references, see Berman (1991 CanMin), Berman
    and Aranovich (1996 CMP), and Aranovich and Berman (1997 AmMin).
    
    This code has two 'runmodes' in which any number of desired calculations can be done:
    
      1) PT is calculated sequentially from input mineral compositions.
          For example, if 4 garnet, 4 orthopyroxene, and 4 plagioclase compositions were
          entered into the input file, the code would calculate PT for:
              gar1-opx1-pl1, gar2-opx2-pl2, gar3-opx3-pl3, gar4-opx4-pl4.
          (biotite and cordierite could also be included)
          This runmode is only available if an equal number of compositions were
          input for each mineral.
    
      2) PT is calculated for EVERY POSSIBLE combination of input mineral compositions.
          For example if 4 garnet, 3 orthopyroxene, and 2 plagioclase compositions were
          entered into the input file, the code would calculate PT for:
              gar1-opx1-pl1, gar1-opx1-pl2, gar1-opx2-pl1, gar1-opx2-pl2, gar1-opx3-pl1, gar1-opx3-pl2
              gar2-opx1-pl1, gar2-opx1-pl2, gar2-opx2-pl1, gar2-opx2-pl2, gar1-opx3-pl1, gar2-opx3-pl2
              gar3-opx1-pl1, gar3-opx1-pl2, gar3-opx2-pl1, gar3-opx2-pl2, gar3-opx3-pl1, gar1-opx3-pl2
              gar4-opx1-pl1, gar4-opx1-pl2, gar4-opx2-pl1, gar4-opx2-pl2, gar4-opx3-pl1, gar4-opx3-pl2
          (biotite and cordierite could also be included)
          If an unequal number of compositions were input for each mineral, the code will run this
          mode by default.
          A great advantage of this runmode is to determine the senstivity of the calculated PT to the choice
          of input mineral compositions.
    
##### What you need to run this code ######
    The required input data are mineral cations for garnet (normalized to 12 O),
    orthpyroxene (normalized to 6 O), and plagioclase (normalized to 8 O) as well as
    modes for each of these minerals. Optional input data are cations in cordierite
    (normalized to 18 O) and biotite (normalized to 11 O [10 O and 2 (OH,Cl,F]); if
    cordierite or biotite are included, their modes must also be included.
    
    Mineral compositions input files are tab-delmited text documents: 'gar.txt',
    'opx.txt', 'pl.txt', 'bt.txt', and 'crd.txt'. Any number of mineral compositions
    can be included in the input files. The order of cations in the input files is
    the same for all minerals: Si, Ti, Al, Cr, Fe3, Fe2, Mn, Mg, Ca, Na, K. If any
    of these elements was not measured / was not present, a value of 0 should be entered.
    Example input of five orthopyroxene measurements (Fe3 calculated using Tim
    Holland's AX2 program):
    
      Si	Ti	Al	Cr	Fe3	Fe2	Mn	Mg	Ca	Na	K
      1.798	0.006	0.356	0.001	0.035	0.56	0.003	1.238	0.004	0	0
      1.791	0.006	0.364	0.002	0.039	0.556	0.002	1.236	0.004	0	0
      1.803	0.005	0.36	0.002	0.021	0.571	0.002	1.23	0.004	0.001	0
      1.787	0.005	0.366	0.001	0.049	0.546	0.002	1.238	0.005	0.001	0
      1.797	0.009	0.366	0.001	0.022	0.569	0.002	1.229	0.005	0.001	0
    
    Ferric iron: Ferric iron in orthopyroxene, in particular, can moderately influence the calculated
    Al-in-opx temperature (although typically, the effect is < 50C). The input files allow for entering
    Fe3 and Fe2 contents of minerals separately. If Fe3 is not being considered, enter zeros for that
    column of the input file and put all Fe as Fe2.
    
##### Uncertainty ######
    quantifying and reporting uncertainty in phase-equilibrium calculations is
    very difficult. This code does not output an uncertainty. The commonly quoted
    rule-of-thumb +/- 50 C and +/- 1 kbar is probably appropriate
    (Pattison and Chacko, pers comm); however, due to the attempt to correct for
    diffusional re-equilibration mineral compositions during cooling, this
    thermobarometer is likely more accurate than most other methods of thermobarometry
    applied to granulite-facies rocks.
    
    Al-in-orthpyroxene thermobarometry (as used in this code) has been
    evaluated in comparison with other methods of thermobarometry on ultrahigh-temperature
    rocks. It agrees well with estimates of peak metamorphism from petrographic interpretations
    of peak mineral assemblages in equilibrium-assemblage diagrams (EADs: aka pseudosections), Al-in-opx
    isopleths in EADs, as well as the highest Zr concentrations measured in rutile.
    Although a straight-forward, easily quotable uncertainty cannot be give, this
    method of thermobarometry appears to be one of the most reliable for constraining
    the peak PT of granulite-facies metamorphism.
    
    There is additional uncertainty in the choice of which mineral compositions to use.
    Runmode 2 of this code can be used to evaluate the sensitivty of the calculation to
    such choices by calculating PT for all possible combinations of input mineral data.
    e.g. Horton, Holder, and Swindle (2021 JMG)
