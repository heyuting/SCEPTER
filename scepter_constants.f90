module scepter_constants
    implicit none
!-----------------------------
#ifdef mod_basalt_cmp
#include <basalt_defines.h>
#endif 
        public
        
        integer,parameter :: nt = 50000000

        real(kind=8),parameter :: rg = 8.3d-3   !  kJ mol^-1 K^-1
        real(kind=8),parameter :: rg2 = 8.2d-2  ! L mol^-1 atm K^-1

        real(kind=8),parameter :: tempk_0 = 273d0
        real(kind=8),parameter :: sec2yr = 60d0*60d0*24d0*365d0

        real(kind=8),parameter :: n2c_g1 = 0.1d0 ! N to C ratio for OM-G1; Could be related to reactivity cf. Janssen 1996
        real(kind=8),parameter :: n2c_g2 = 0.1d0
        real(kind=8),parameter :: n2c_g3 = 0.1d0

        real(kind=8),parameter :: oxa2c_g2 = 0.01d0 ! amount of oxlate (C2O4=) released (DEF)
        ! real(kind=8),parameter :: oxa2c_g2 = 0.5d0 ! amount of oxlate (C2O4=) released  
        ! real(kind=8),parameter :: oxa2c_g2 = 0.2d0 ! amount of oxlate (C2O4=) released 
        ! real(kind=8),parameter :: oxa2c_g2 = 0.1d0 ! amount of oxlate (C2O4=) released  
        ! real(kind=8),parameter :: oxa2c_g2 = 0.0d0 ! amount of oxlate (C2O4=) released 

        real(kind=8),parameter :: ac2c_g2 = 0.01d0 ! amount of acetate (CH3COO-) released

        real(kind=8),parameter :: fr_an_ab = 0.0d0 ! Anorthite fraction for albite (Beerling et al., 2020); 0.0 - 0.1
        real(kind=8),parameter :: fr_an_olg = 0.2d0 ! Anorthite fraction for oligoclase (Beerling et al., 2020); 0.1 - 0.3
        real(kind=8),parameter :: fr_an_and = 0.4d0 ! Anorthite fraction for andesine (Beerling et al., 2020); 0.3 - 0.5
        real(kind=8),parameter :: fr_an_la = 0.6d0 ! Anorthite fraction for labradorite (Beerling et al., 2020); 0.5 - 0.7
        real(kind=8),parameter :: fr_an_by = 0.8d0 ! Anorthite fraction for bytownite (Beerling et al., 2020); 0.7 - 0.9
        real(kind=8),parameter :: fr_an_an = 1.0d0 ! Anorthite fraction for anorthite (Beerling et al., 2020); 0.9 - 1.0

        real(kind=8),parameter :: fr_hb_cpx = 0.5d0 ! Hedenbergite fraction for clinopyroxene; 0.0 - 1.0
        real(kind=8),parameter :: fr_fer_opx = 0.5d0 ! Ferrosilite fraction for orthopyroxene; 0.0 - 1.0
        real(kind=8),parameter :: fr_fer_agt = 0.0d0 ! Ferrosilite (and Hedenbergite; or Fe/(Fe+Mg)) fraction for Augite; 0.0 - 1.0; Beerling et al 2020 
        real(kind=8),parameter :: fr_opx_agt = 0.0d0 ! OPX (or 1 - Ca:(Fe+Mg)) fraction for Augite; 0.0 - 1.0; from Beerling et al 2020
        real(kind=8),parameter :: fr_napx_agt = 0.1d0 ! Na fraction for Augite (or Na/(Ca+Fe+Mg)); 0.0 - 1.0; from Beerling et al 2020

#ifndef mod_basalt_cmp
        real(kind=8),parameter :: fr_si_gbas = 1d0 ! Si fraction of glass basalt; Pollyea and Rimstidt 2017 (referring to basalt used by Oelkers and Gislason (2001) and Gundbrandsson et al. (2011)
        real(kind=8),parameter :: fr_al_gbas = 0.358d0 ! Al fraction of glass basalt
        real(kind=8),parameter :: fr_na_gbas = 0.079d0 ! Na fraction of glass basalt
        real(kind=8),parameter :: fr_k_gbas = 0.08d0 ! K fraction of glass basalt
        real(kind=8),parameter :: fr_mg_gbas = 0.281d0 ! Mg fraction of glass basalt
        real(kind=8),parameter :: fr_ca_gbas = 0.264d0 ! Ca fraction of glass basalt
        real(kind=8),parameter :: fr_fe2_gbas = 0.190d0 ! Fe2 fraction of glass basalt

        real(kind=8),parameter :: fr_si_cbas = 1d0 ! Si fraction of clystaline basalt; Pollyea and Rimstidt 2017 (referring to basalt used by Oelkers and Gislason (2001) and Gundbrandsson et al. (2011)
        real(kind=8),parameter :: fr_al_cbas = 0.358d0 ! Al fraction of clystaline basalt
        real(kind=8),parameter :: fr_na_cbas = 0.079d0 ! Na fraction of clystaline basalt
        real(kind=8),parameter :: fr_k_cbas = 0.08d0 ! K fraction of clystaline basalt
        real(kind=8),parameter :: fr_mg_cbas = 0.281d0 ! Mg fraction of clystaline basalt
        real(kind=8),parameter :: fr_ca_cbas = 0.264d0 ! Ca fraction of clystaline basalt
        real(kind=8),parameter :: fr_fe2_cbas = 0.190d0 ! Fe2 fraction of clystaline basalt
#else
        real(kind=8),parameter :: fr_si_gbas  = def_bas_si_fr
        real(kind=8),parameter :: fr_al_gbas  = def_bas_al_fr
        real(kind=8),parameter :: fr_na_gbas  = def_bas_na_fr
        real(kind=8),parameter :: fr_k_gbas   = def_bas_k_fr
        real(kind=8),parameter :: fr_mg_gbas  = def_bas_mg_fr
        real(kind=8),parameter :: fr_ca_gbas  = def_bas_ca_fr
        real(kind=8),parameter :: fr_fe2_gbas = def_bas_fe2_fr
        real(kind=8),parameter :: fr_si_cbas  = def_bas_si_fr
        real(kind=8),parameter :: fr_al_cbas  = def_bas_al_fr
        real(kind=8),parameter :: fr_na_cbas  = def_bas_na_fr
        real(kind=8),parameter :: fr_k_cbas   = def_bas_k_fr
        real(kind=8),parameter :: fr_mg_cbas  = def_bas_mg_fr
        real(kind=8),parameter :: fr_ca_cbas  = def_bas_ca_fr
        real(kind=8),parameter :: fr_fe2_cbas = def_bas_fe2_fr
#endif

        real(kind=8),parameter :: mvka = 99.52d0 ! cm3/mol; molar volume of kaolinite; Robie et al. 1978
        real(kind=8),parameter :: mvfo = 43.79d0 ! cm3/mol; molar volume of Fo; Robie et al. 1978
        real(kind=8),parameter :: mvab_0 = 100.07d0 ! cm3/mol; molar volume of Ab(NaAlSi3O8); Robie et al. 1978 
        real(kind=8),parameter :: mvan_0 = 100.79d0 ! cm3/mol; molar volume of An (CaAl2Si2O8); Robie et al. 1978
        !===========================================
        ! Mineral Molar Volumes (cm3/mol)
        !===========================================
        ! Primary Silicates
        real(kind=8),parameter :: mvab = fr_an_ab*mvan_0 + (1d0-fr_an_ab)*mvab_0 ! cm3/mol; molar volume of albite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
        real(kind=8),parameter :: mvan = fr_an_an*mvan_0 + (1d0-fr_an_an)*mvab_0 ! cm3/mol; molar volume of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
        real(kind=8),parameter :: mvby = fr_an_by*mvan_0 + (1d0-fr_an_by)*mvab_0 ! cm3/mol; molar volume of bytownite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
        real(kind=8),parameter :: mvla = fr_an_la*mvan_0 + (1d0-fr_an_la)*mvab_0 ! cm3/mol; molar volume of labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mvand = fr_an_and*mvan_0 + (1d0-fr_an_and)*mvab_0 ! cm3/mol; molar volume of andesine (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mvolg = fr_an_olg*mvan_0 + (1d0-fr_an_olg)*mvab_0 ! cm3/mol; molar volume of oligoclase (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mvcc = 36.934d0 ! cm3/mol; molar volume of Cc (CaCO3); Robie et al. 1978
        real(kind=8),parameter :: mvpy = 23.94d0 ! cm3/mol; molar volume of Pyrite (FeS2); Robie et al. 1978
        real(kind=8),parameter :: mvamal = 31.956d0 ! cm3/mol; assuming amorphous Al has molar volume of Gibsite (Al(OH)3); Robie et al. 1978
        real(kind=8),parameter :: mvgb = 31.956d0 ! cm3/mol; molar volume of Gibsite (Al(OH)3); Robie et al. 1978
        ! Carbonates and Sulfides   
        real(kind=8),parameter :: mvct = 108.5d0 ! cm3/mol; molar volume of Chrysotile (Mg3Si2O5(OH)4); Robie et al. 1978
        real(kind=8),parameter :: mvfa = 46.39d0 ! cm3/mol; molar volume of Fayalite (Fe2SiO4); Robie et al. 1978
        real(kind=8),parameter :: mvdlm = 64.34d0 ! cm3/mol; molar volume of Dolomite (CaMg(CO3)2); Robie et al. 1978

        ! Amorphous and Secondary Mineralsfa = 46.39d0 ! cm3/mol; molar volume of Fayalite (Fe2SiO4); Robie et al. 1978
        real(kind=8),parameter :: mvamfe3 = 20.82d0 ! cm3/mol; assuming amorphous Fe(OH)3 has molar volume of Goethite (FeO(OH)); Robie et al. 1978
        real(kind=8),parameter :: mvgt = 20.82d0 ! cm3/mol; molar volume of Goethite (FeO(OH)); Robie et al. 1978
        real(kind=8),parameter :: mvcabd = 129.77d0 ! cm3/mol; molar volume of Ca-beidellite (Ca(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvkbd = 134.15d0 ! cm3/mol; molar volume of K-beidellite (K(1/3)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvnabd = 130.73d0 ! cm3/mol; molar volume of Na-beidellite (Na(1/3)Al(7/3)S

        ! Clay Mineralsi(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvmgbd = 128.73d0 ! cm3/mol; molar volume of Mg-beidellite (Mg(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvcasp = 134.359d0 ! cm3/mol; molar volume of Ca-saponite (Ca(1/6)Mg3Al(1/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvksp = 138.745d0 ! cm3/mol; molar volume of K-saponite (K(1/3)Mg3Al(1/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvnasp = 135.320d0 ! cm3/mol; molar volume of Na-saponite (Na(1/3)Mg3Al(1/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvmgsp = 132.602d0 ! cm3/mol; molar volume of Mg-saponite (Mg(1/6)Mg3Al(1/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvdp = 66.09d0 ! cm3/mol; molar volume of Diopside (MgCaSi2O6);  Robie et al. 1978
        real(kind=8),parameter :: mvhb = 248.09d0/3.55d0 ! cm3/mol; molar volume of Hedenbergite (FeCaSi2O6); from a webpage
        real(kind=8),parameter :: mvcpx = fr_hb_cpx*mvhb + (1d0-fr_hb_cpx)*mvdp  ! cm3/mol; molar volume of clinopyroxene (FexMg(1-x)CaSi2O6); assuming simple ('ideal'?) mixing
    
        ! Pyroxenes and Related Minerals
        real(kind=8),parameter :: mvkfs = 108.72d0 ! cm3/mol; molar volume of K-feldspar (KAlSi3O8); Robie et al. 1978
        real(kind=8),parameter :: mvom = 30d0/1.5d0 ! cm3/mol; molar volume of OM (CH2O); calculated assuming 30 g/mol of molar weight and 1.2 g/cm3 of density (Mayer et al., 2004; Ruhlmann et al.,2006)
        real(kind=8),parameter :: mvomb = 30d0/1.5d0 ! cm3/mol; assumed to be same as mvom
        real(kind=8),parameter :: mvg1 = 30d0/1.5d0 ! cm3/mol; assumed to be same as mvom
        real(kind=8),parameter :: mvg2 = 30d0/1.5d0 ! cm3/mol; assumed to be same as mvom
        real(kind=8),parameter :: mvg3 = 30d0/1.5d0 ! cm3/mol; assumed to be same as mvom
        real(kind=8),parameter :: mvamsi = 25.739d0 ! cm3/mol; molar volume of amorphous silica taken as cristobalite (SiO2); Robie et al. 1978
        real(kind=8),parameter :: mvphsi = 25.739d0 ! cm3/mol; molar volume of phytolith silica taken as cristobalite (SiO2); Robie et al. 1978
        real(kind=8),parameter :: mvarg = 34.15d0 ! cm3/mol; molar volume of aragonite; Robie et al. 1978
        real(kind=8),parameter :: mvhm = 30.274d0 ! cm3/mol; molar volume of hematite; Robie et al. 1978
        real(kind=8),parameter :: mvill = 139.35d0 ! cm3/mol; molar volume of illite (K0.6Mg0.25Al2        real(kind=8),parameter :: mvqtz = 22.688d0 ! cm3/mol; molar volume of quartz (SiO2); Robie et al. 1978
        real(kind=8),parameter :: mvanl = 97.49d0 ! cm3/mol; molar volume of analcime (NaAlSi2O6*H2O); Robie et al. 1978
        real(kind=8),parameter :: mvnph = 54.16d0 ! cm3/mol; molar volume of nepheline (NaAlSiO4); Robie et al. 1978

        ! Other Minerals.16d0 ! cm3/mol; molar volume of nepheline (NaAlSiO4); Robie et al. 1978
        real(kind=8),parameter :: mvqtz = 22.688d0 ! cm3/mol; molar volume of quartz (SiO2); Robie et al. 1978
        real(kind=8),parameter :: mvgps = 74.69d0 ! cm3/mol; molar volume of gypsum (CaSO4*2H2O); Robie et al. 1978
        real(kind=8),parameter :: mvtm = 272.92d0 ! cm3/mol; molar volume of tremolite (Ca2Mg5(Si8O22)(OH)2); Robie et al. 1978
        real(kind=8),parameter :: mven = 31.31d0 ! cm3/mol; molar volume of enstatite (MgSiO3); te (FeSiO3); Robie and Hemingway 1995
        real(kind=8),parameter :: mvfer = 33.00d0 ! cm3/mol; molar volume of ferrosilite (FeSiO3); Robie and Hemingway 1995
        real(kind=8),parameter :: mvopx = fr_fer_opx*mvfer +(1d0-fr_fer_opx)*mven !  cm3/mol; molar volume of clinopyroxene (FexMg(1-x)SiO3); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mvmscv = 140.71d0 ! cm3/mol; molar volume of muscovite (KAl2(AlSi3O10)(OH)2); Robie et al. 1978
        real(kind=8),parameter :: mvplgp = 149.91d0 ! cm3/mol; molar volume of phlogopite (KMg3(AlSi3O10)(OH)2); Robie et al. 1978
        real(kind=8),parameter :: mvantp = 274.00d0 ! cm3/mol; molar volume of anthophyllite (Mg7Si8O22(OH)2); Robie and Bethke 1962
        real(kind=8),parameter :: mvsplt = 285.600d0 ! cm3/mol; molar volume of sepiolite (Mg4Si6O15(OH)2:6H2O); Wolery and Jove-Colon 2004
        real(kind=8),parameter :: mvjd = 60.4d0 ! cm3/mol; molar volume of jadeite (NaAlSi2O6); Robie et al. 1978
        real(kind=8),parameter :: mvwls = 39.93d0 ! cm3/mol; molar volume of wollastonite (CaSiO3); Robie et al. 1978
        real(kind=8),parameter :: mvagt = ( & 
                                        & (fr_fer_agt*mvfer +(1d0-fr_fer_agt)*mven)*2d0*fr_opx_agt &! (Fe2xyMg2(1-x)ySi2yO6y)
                                        & + (fr_fer_agt*mvhb + (1d0-fr_fer_agt)*mvdp)*(1d0-fr_opx_agt) &! (Fex(1-y)Mg(1-x)(1-y)Ca(1-y)Si2(1-y)O6(1-y))
                                        & ) * (1d0-fr_napx_agt) &! non-Na pyroxene 
                                        & + ( &
                                        & mvjd &! NaAlSi2O6  
                                        & ) *  fr_napx_agt ! fraction of Na pyroxene
                                        !  cm3/mol; molar volume of augite 
                                        ! (Fe(2xy+x(1-y))Mg(2y-2xy+1+xy-x-y)Ca(1-y)Si2O6 = Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6) ! non-Na pyroxene
                                        ! Fe(xy+x)(1-z)Mg(y-xy+1-x)(1-z)Ca(1-y)(1-z)Si2(1-z)O6(1-z) + NazAlzSi2zO6z
                                        ! = Fe(xy+x)(1-z)Mg(y-xy+1-x)(1-z)Ca(1-y)(1-z)NazAlzSi2O6
                                        ! ; assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mvamnt = 46.40173913043478d0 ! cm3/mol; molar volume of ammonium nitrate (NH4NO3); density 1.725 g/cm3 (at 20C) from wikipedea 
        real(kind=8),parameter :: mvfe2o = 12d0 ! cm3/mol; molar volume of ferrous oxide; Robie et al. 1978
        real(kind=8),parameter :: mvmgo = 11.248d0 ! cm3/mol; molar volume of periclase; Robie et al. 1978
        real(kind=8),parameter :: mvk2o = 40.38d0 ! cm3/mol; molar volume of dipotasium monoxide; Robie et al. 1978
        real(kind=8),parameter :: mvcao = 16.764d0 ! cm3/mol; molar volume of calcium monoxide; Robie et al. 1978
        real(kind=8),parameter :: mvna2o = 25.88d0 ! cm3/mol; molar volume of disodium monoxide; Robie et al. 1978
        real(kind=8),parameter :: mval2o3 = 25.575d0 ! cm3/mol; molar volume of corundum; Robie et al. 1978
        real(kind=8),parameter :: mvsio2 = mvqtz ! cm3/mol; molar volume of SiO2; fake material that easily dissolves with otherwise quartz property
        ! real(kind=8),parameter :: mvgbas = ( &
                                        ! & fr_si_gbas*mvamsi + fr_al_gbas/2d0*mval2o3 + fr_na_gbas/2d0*mvna2o &
                                        ! & + fr_k_gbas/2d0*mvk2o + fr_ca_gbas*mvcao + fr_mg_gbas*mvmgo + fr_fe2_gbas*mvfe2o &
                                        ! & ) ! assuming simply mixing molar volume?
        ! real(kind=8),parameter :: mvcbas = ( &
                                        ! & fr_si_gbas*mvqtz + fr_al_cbas/2d0*mval2o3 + fr_na_cbas/2d0*mvna2o &
                                        ! & + fr_k_cbas/2d0*mvk2o + fr_ca_cbas*mvcao + fr_mg_cbas*mvmgo + fr_fe2_cbas*mvfe2o &
                                        ! & ) ! assuming simply mixing molar volume? 
        real(kind=8),parameter :: mvep = 139.10d0 ! cm3/mol; molar volume of epidote; Ca2FeAl2Si3O12OH; Gottschalk 2004 originally from Holland and Powell 1998
        real(kind=8),parameter :: mvclch = 211.470d0 ! cm3/mol; molar volume of clinochlore; Mg5Al2Si3O10(OH)8; Roots 1994 Eur. J. Mineral.
        real(kind=8),parameter :: mvsdn = 109.05d0 ! cm3/mol; molar volume of sanidine; Robie et al. 1978
        real(kind=8),parameter :: mvcdr = 584.95d0/2.65d0  ! cm3/mol; molar volume of cordierite Mg2Al4Si5O18 from molar weight and density from http://www.webmineral.com/data/Cordierite.shtml#.YZtITrqIaUk
        real(kind=8),parameter :: mvleu =  88.39d0 ! cm3/mol; molar volume of leucite KAlSi206; Robie et al. 1978
        real(kind=8),parameter :: mvkcl =  37.524d0 ! cm3/mol; molar volume of sylvite; from Robie et al. 1978
        real(kind=8),parameter :: mvgac =  60.052d0/1.27d0 ! cm3/mol; from CH3COOH Molar mass 60.052 g·mol−1 and Density 1.27 g/cm3 from Wikipedea https://en.wikipedia.org/wiki/Acetic_acid
        real(kind=8),parameter :: mvmesmh =  213.25d0/0.56d0 ! cm3/mol; from MES monohydrate Molar mass 213.25 g·mol−1 and Density 0.56 g/cm3 from https://www.emdmillipore.com/US/en/product/2-Morpholinoethanesulfonic-acid-monohydrate,MDA_CHEM-106126
        real(kind=8),parameter :: mvims =  68.077d0/1.23d0 ! cm3/mol; from MES monohydrate Molar mass 68.077 g·mol−1 and Density 1.23 g/cm3 from https://en.wikipedia.org/wiki/Imidazole
        real(kind=8),parameter :: mvteas = 149.190d0/1.124d0 ! cm3/mol; from MES monohydrate Molar mass 149.190 g·mol−1 and Density 1.124 g/cm3 from https://en.wikipedia.org/wiki/Triethanolamine
        real(kind=8),parameter :: mvnaoh = 39.9971d0/2.13d0 ! cm3/mol; from NaOH Molar mass 39.9971 g·mol−1 and Density 2.13 g/cm3 from https://en.wikipedia.org/wiki/Sodium_hydroxide
        real(kind=8),parameter :: mvnaglp = 216.036d0/1.5d0 ! cm3/mol; from sodium glycerophosphate Molar mass 216.036 g·mol−1 and assuming Density 1.5 g/cm3 
        real(kind=8),parameter :: mvcacl2 = 50.75d0 ! cm3/mol; molar volume of hydrophilite; from Robie et al. 1978
        real(kind=8),parameter :: mvnacl = 27.015d0 ! cm3/mol; molar volume of halite; from Robie et al. 1978
        real(kind=8),parameter :: mvcaso4 = 45.94d0 ! cm3/mol; molar volume of anhydrite; from Robie et al. 1978
        real(kind=8),parameter :: mvinrt = mvka ! cm3/mol; molar volume of kaolinite; Robie et al. 1978
        ! real(kind=8),parameter :: mvinrt = mvqtz ! cm3/mol; molar volume of quartz; Robie et al. 1978

                                        
        real(kind=8),parameter :: mwtka = 258.162d0 ! g/mol; formula weight of Ka; Robie et al. 1978
        real(kind=8),parameter :: mwtfo = 140.694d0 ! g/mol; formula weight of Fo; Robie et al. 1978
        real(kind=8),parameter :: mwtab_0 = 262.225d0 ! g/mol; formula weight of Ab; Robie et al. 1978
        real(kind=8),parameter :: mwtan_0 = 278.311d0 ! g/mol; formula weight of An; Robie et al. 1978
        real(kind=8),parameter :: mwtab = fr_an_ab*mwtan_0 + (1d0-fr_an_ab)*mwtab_0 ! g/mol; formula weight of albte (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtan = fr_an_an*mwtan_0 + (1d0-fr_an_an)*mwtab_0 ! g/mol; formula weight of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtby = fr_an_by*mwtan_0 + (1d0-fr_an_by)*mwtab_0 ! g/mol; formula weight of bytownite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtla = fr_an_la*mwtan_0 + (1d0-fr_an_la)*mwtab_0 ! g/mol; formula weight of labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtand = fr_an_and*mwtan_0 + (1d0-fr_an_and)*mwtab_0 ! g/mol; formula weight of andesine (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtolg = fr_an_olg*mwtan_0 + (1d0-fr_an_olg)*mwtab_0 ! g/mol; formula weight of oligoclase (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtcc = 100.089d0 ! g/mol; formula weight of Cc; Robie et al. 1978
        real(kind=8),parameter :: mwtpy = 119.967d0 ! g/mol; formula weight of Py; Robie et al. 1978
        real(kind=8),parameter :: mwtamal = 78.004d0 ! g/mol; assuming amorphous Al(OH)3 has weight of Gb; Robie et al. 1978
        real(kind=8),parameter :: mwtgb = 78.004d0 ! g/mol; formula weight of Gb; Robie et al. 1978
        real(kind=8),parameter :: mwtct = 277.113d0 ! g/mol; formula weight of Ct; Robie et al. 1978
        real(kind=8),parameter :: mwtfa = 203.778d0 ! g/mol; formula weight of Fa; Robie et al. 1978
        real(kind=8),parameter :: mwtamfe3 = 88.854d0 ! g/mol; assuming amorphous Fe(OH)3 has formula weight of Gt; Robie et al. 1978
        real(kind=8),parameter :: mwtgt = 88.854d0 ! g/mol; formula weight of Gt; Robie et al. 1978
        real(kind=8),parameter :: mwtcabd = 366.6252667d0 ! g/mol; formula weight of Cabd calculated from atmoic weight
        real(kind=8),parameter :: mwtkbd = 372.9783667d0 ! g/mol; formula weight of Kbd calculated from atmoic weight
        real(kind=8),parameter :: mwtnabd = 367.6088333d0 ! g/mol; formula weight of Nabd calculated from atmoic weight
        real(kind=8),parameter :: mwtmgbd = 363.9964333d0 ! g/mol; formula weight of Mgbd calculated from atmoic weight
        real(kind=8),parameter :: mwtcasp = 385.5777533d0 ! g/mol; formula weight of Ca-saponite calculated from atmoic weight
        real(kind=8),parameter :: mwtksp = 391.93052d0 ! g/mol; formula weight of K-saponite calculated from atmoic weight
        real(kind=8),parameter :: mwtnasp = 386.56101d0 ! g/mol; formula weight of Na-saponite calculated from atmoic weight
        real(kind=8),parameter :: mwtmgsp = 382.9485867d0 ! g/mol; formula weight of Mg-saponite calculated from atmoic weight
        real(kind=8),parameter :: mwtdp = 216.553d0 ! g/mol;  Robie et al. 1978
        real(kind=8),parameter :: mwthb = 248.09d0 ! g/mol; from a webpage
        real(kind=8),parameter :: mwtcpx = fr_hb_cpx*mwthb + (1d0-fr_hb_cpx)*mwtdp ! g/mol; formula weight of clinopyroxene (FexMg(1-x)CaSi2O6); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtkfs = 278.33d0 ! g/mol; formula weight of Kfs; Robie et al. 1978
        real(kind=8),parameter :: mwtom = 30d0 ! g/mol; formula weight of CH2O
        real(kind=8),parameter :: mwtomb = 30d0 ! g/mol; formula weight of CH2O
        real(kind=8),parameter :: mwtg1 = 30d0 ! g/mol; formula weight of CH2O
        real(kind=8),parameter :: mwtg2 = 30d0 ! g/mol; formula weight of CH2O
        real(kind=8),parameter :: mwtg3 = 30d0 ! g/mol; formula weight of CH2O
        real(kind=8),parameter :: mwtamsi = 60.085d0 ! g/mol; formula weight of amorphous silica
        real(kind=8),parameter :: mwtphsi = 60.085d0 ! g/mol; formula weight of phytolith silica
        real(kind=8),parameter :: mwtarg = 100.089d0 ! g/mol; formula weight of aragonite
        real(kind=8),parameter :: mwtdlm = 184.403d0 ! g/mol; formula weight of dolomite
        real(kind=8),parameter :: mwthm = 159.692d0 ! g/mol; formula weight of hematite
        real(kind=8),parameter :: mwtill = 383.90053d0 ! g/mol; formula weight of Ill calculated from atmoic weight
        real(kind=8),parameter :: mwtanl = 220.155d0 ! g/mol; formula weight of analcime
        real(kind=8),parameter :: mwtnph = 142.055d0 ! g/mol; formula weight of nepheline
        real(kind=8),parameter :: mwtqtz = 60.085d0 ! g/mol; formula weight of quartz
        real(kind=8),parameter :: mwtgps = 172.168d0 ! g/mol; formula weight of gypsum
        real(kind=8),parameter :: mwttm = 812.374d0 ! g/mol; formula weight of tremolite
        real(kind=8),parameter :: mwten = 100.389d0 ! g/mol; formula weight of enstatite
        real(kind=8),parameter :: mwtfer = 131.931d0 ! g/mol; formula weight of ferrosilite
        real(kind=8),parameter :: mwtopx = fr_fer_opx*mwtfer + (1d0 -fr_fer_opx)*mwten ! g/mol; formula weight of clinopyroxene (FexMg(1-x)SiO3); assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtmscv = 398.311d0 ! g/mol; formula weight of muscovite
        real(kind=8),parameter :: mwtplgp = 417.262d0 ! g/mol; formula weight of phlogopite
        real(kind=8),parameter :: mwtantp = 780.976d0 ! g/mol; formula weight of anthophyllite
        real(kind=8),parameter :: mwtsplt = 647.8304d0 ! g/mol; formula weight of sepiolite
        real(kind=8),parameter :: mwtjd = 202.140d0 ! g/mol; formula weight of jadeite; Robie et al. 1978
        real(kind=8),parameter :: mwtwls = 116.164d0 ! g/mol; formula weight of wollastonite; Robie et al. 1978
        real(kind=8),parameter :: mwtagt = ( &
                                        & (fr_fer_agt*mwtfer +(1d0-fr_fer_agt)*mwten)*2d0*fr_opx_agt &! (Fe2xyMg2(1-x)ySi2yO6y) or y(Fe2xMg2(1-x)Si2O6)
                                        & + (fr_fer_agt*mwthb + (1d0-fr_fer_agt)*mwtdp)*(1d0-fr_opx_agt) &! (Fex(1-y)Mg(1-x)(1-y)Ca(1-y)Si2(1-y)O6(1-y)) or (1-y)(FexMg(1-x)CaSi2O6)
                                        & ) * (1d0 - fr_napx_agt) &! (1 - z) (Fex(1-y+2y)Mg(1-x)(1-y+2y)Ca(1-y)Si2O6)
                                        & + ( &
                                        & mwtjd &!  NaAlSi2O6
                                        & ) * fr_napx_agt  ! z
                                        !  g/mol; formula weight of augite 
                                        ! (Fe(2xy+x(1-y))Mg(2y-2xy+1+xy-x-y)Ca(1-y)Si2O6 = Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6) ! non-Na pyroxene
                                        ! Fe(xy+x)(1-z)Mg(y-xy+1-x)(1-z)Ca(1-y)(1-z)Si2(1-z)O6(1-z) + NazAlzSi2zO6z
                                        ! = Fe(xy+x)(1-z)Mg(y-xy+1-x)(1-z)Ca(1-y)(1-z)NazAlzSi2O6
                                        ! ; assuming simple ('ideal'?) mixing
        real(kind=8),parameter :: mwtamnt = 80.043d0 ! g/mol; formula weight of ammonium nitrate 
        real(kind=8),parameter :: mwtfe2o = 71.846d0 ! g/mol; molar weight of ferrous oxide; Robie et al. 1978 
        real(kind=8),parameter :: mwtmgo = 40.304d0 ! g/mol; molar weight of periclase; Robie et al. 1978 
        real(kind=8),parameter :: mwtk2o = 94.195d0 ! g/mol; molar weight of dipotasium monoxide; Robie et al. 1978 
        real(kind=8),parameter :: mwtcao = 56.079d0 ! g/mol; molar weight of calcium oxide; Robie et al. 1978 
        real(kind=8),parameter :: mwtna2o = 61.979d0 ! g/mol; molar weight of disodium monoxide; Robie et al. 1978 
        real(kind=8),parameter :: mwtal2o3 = 101.962d0 ! g/mol; molar weight of corundum; Robie et al. 1978 
        real(kind=8),parameter :: mwtsio2 = mwtqtz ! g/mol; molar weight of SiO2
        real(kind=8),parameter :: mwtep = 483.22675d0 ! g/mol; molar weight of epidote; Ca2FeAl2Si3O12OH; calculated from elemental molar weight
        real(kind=8),parameter :: mwtclch = 555.79754d0 ! g/mol; molar weight of clinochlore; Mg5Al2Si3O10(OH)8; calculated from elemental molar weight
        real(kind=8),parameter :: mwtsdn = 278.333d0 ! cm3/mol; molar weight of sanidine; Robie et al. 1978
        real(kind=8),parameter :: mwtcdr = 584.95d0 ! cm3/mol; from molar weight of cordierite Mg2Al4Si5O18 from http://www.webmineral.com/data/Cordierite.shtml#.YZtITrqIaUk
        real(kind=8),parameter :: mwtleu =  218.248d0 ! cm3/mol; molar volume of leucite KAlSi206; Robie et al. 1978
        real(kind=8),parameter :: mwtkcl =  74.551d0 ! g·mol−1; molar weight of sylvite; Robie et al. 1978
        real(kind=8),parameter :: mwtgac =  60.052d0 ! g·mol−1; from CH3COOH Molar mass 60.052 g·mol−1 and Density 1.27 g/cm3 from Wikipedea https://en.wikipedia.org/wiki/Acetic_acid
        real(kind=8),parameter :: mwtmesmh =  213.25d0 ! g·mol−1; from MES monohydrate Molar mass 213.25 g·mol−1 and Density 0.56 g/cm3 from https://www.emdmillipore.com/US/en/product/2-Morpholinoethanesulfonic-acid-monohydrate,MDA_CHEM-106126
        real(kind=8),parameter :: mwtims =  68.077d0 ! g/mol; from Imidazole monohydrate Molar mass 68.077 g·mol−1 and Density 1.23 g/cm3 from https://en.wikipedia.org/wiki/Imidazole
        real(kind=8),parameter :: mwtteas = 149.190d0 ! g/mol; from TEA Molar mass 149.190 g·mol−1 and Density 1.124 g/cm3 from https://en.wikipedia.org/wiki/Triethanolamine
        real(kind=8),parameter :: mwtnaoh = 39.9971d0 ! g/mol; from NaOH Molar mass 39.9971 g·mol−1 and Density 2.13 g/cm3 from https://en.wikipedia.org/wiki/Sodium_hydroxide
        real(kind=8),parameter :: mwtnaglp = 216.036d0 ! g/mol; from NaOH Molar mass 216.036 g·mol−1 from https://en.wikipedia.org/wiki/Sodium_glycerophosphate
        real(kind=8),parameter :: mwtcacl2 = 110.986d0 ! g/mol; molar weight of hydrophilite; Robie et al. 1978
        real(kind=8),parameter :: mwtnacl = 58.443d0 ! g/mol; molar weight of halite; Robie et al. 1978 
        real(kind=8),parameter :: mwtcaso4 = 136.138d0 ! g/mol; molar weight of halite; Robie et al. 1978 
        real(kind=8),parameter :: mwtinrt = mwtka ! g/mol; formula weight of Ka; Robie et al. 1978
        ! real(kind=8),parameter :: mwtinrt = mwtqtz ! g/mol; formula weight of quartz; Robie et al. 1978
 
        real(kind=8),parameter :: mvgbas = ( &
                                        & fr_si_gbas*mwtamsi + fr_al_gbas/2d0*mwtal2o3 + fr_na_gbas/2d0*mwtna2o &
                                        & + fr_k_gbas/2d0*mwtk2o + fr_ca_gbas*mwtcao + fr_mg_gbas*mwtmgo + fr_fe2_gbas*mwtfe2o &
                                        & )/3.0d0  ! assuming 3.0 g/cm3 particle density
        real(kind=8),parameter :: mvcbas = ( &
                                        & fr_si_cbas*mwtamsi + fr_al_cbas/2d0*mwtal2o3 + fr_na_cbas/2d0*mwtna2o &
                                        & + fr_k_cbas/2d0*mwtk2o + fr_ca_cbas*mwtcao + fr_mg_cbas*mwtmgo + fr_fe2_cbas*mwtfe2o &
                                        & )/3.0d0 ! assuming 3.0 g/cm3 particle density

        real(kind=8),parameter :: mwtgbas = ( &
                                        & fr_si_gbas*mwtamsi + fr_al_gbas/2d0*mwtal2o3 + fr_na_gbas/2d0*mwtna2o &
                                        & + fr_k_gbas/2d0*mwtk2o + fr_ca_gbas*mwtcao + fr_mg_gbas*mwtmgo + fr_fe2_gbas*mwtfe2o &
                                        & ) ! assuming simply mixing molar weight?
        real(kind=8),parameter :: mwtcbas = ( &
                                        & fr_si_cbas*mwtamsi + fr_al_cbas/2d0*mwtal2o3 + fr_na_cbas/2d0*mwtna2o &
                                        & + fr_k_cbas/2d0*mwtk2o + fr_ca_cbas*mwtcao + fr_mg_cbas*mwtmgo + fr_fe2_cbas*mwtfe2o &
                                        & ) ! assuming simply mixing molar weight?
        
        ! cation molar weight from PHREEQC.DAT
        real(kind=8),parameter :: mwtaqna   = 22.9898d0  
        real(kind=8),parameter :: mwtaqk    = 39.102d0  
        real(kind=8),parameter :: mwtaqca   = 40.08d0  
        real(kind=8),parameter :: mwtaqmg   = 24.312d0  
        real(kind=8),parameter :: mwtaqal   = 26.9815d0  
        real(kind=8),parameter :: mwtaqsi   = 28.0843d0  
        real(kind=8),parameter :: mwtaqfe2  = 55.847d0  
        real(kind=8),parameter :: mwtaqfe3  = 55.847d0  
        ! anions 
        real(kind=8),parameter :: mwtaqno3  = 14.0067d0 + 3*16.0d0
        real(kind=8),parameter :: mwtaqso4  = 32.064d0 + 4*16.0d0
        real(kind=8),parameter :: mwtaqcl   = 35.453d0  
        real(kind=8),parameter :: mwtaqoxa  = 1.008d0 + 2*12.0111d0 + 4*16.0d0  ! tracer of oxalate: OxaH- = HC2O4-
        real(kind=8),parameter :: mwtaqac   = 3*1.008d0 + 2*12.0111d0 + 3*16.0d0  ! tracer of acetate: AcO- = C₂H₃O⁻₂
        real(kind=8),parameter :: mwtaqmes  = 12*1.008d0 + 6*12.0111d0 + 14.0067d0 + 4*16.0d0 + 32.064d0  ! tracer of MES: MESH- = C6H12NO4S-
        real(kind=8),parameter :: mwtaqim   = 4*1.008d0 + 3*12.0111d0 + 2*14.0067d0   ! tracer of Imidazole: C3H4N2
        real(kind=8),parameter :: mwtaqtea  = 15*1.008d0 + 6*12.0111d0 + 14.0067d0 + 3*16.0d0   ! tracer of Triethanolamine: C6H15NO3
        real(kind=8),parameter :: mwtaqglp  = 8*1.008d0 + 3*12.0111d0 + 14.0067d0 + 6*16.0d0 + 30.9738d0 !  tracer of glycerophosphate: HL- = C3H8NO6P-
        
        real(kind=8),parameter :: mvblk = mvka ! for bulk soil assumed to be equal to kaolinite
        real(kind=8),parameter :: mwtblk = mwtka
       
        ! real(kind=8),parameter :: w = 1.0d-4 ! m yr^-1, uplift rate
      
        real(kind=8),parameter :: disp_FULL = 1d0

        ! integer,parameter :: nrec_prof = 22
#ifndef nrec_prof_in
        integer,parameter :: nrec_prof = 20
        logical :: linear_rectime = .false. 
#else
        integer,parameter :: nrec_prof = nrec_prof_in
        logical :: linear_rectime = .true. 
#endif 
        integer,parameter :: nrec_flx = 60

        ! type of uplift vs porosity relationship
        ! #ifndef iwtypein 
        ! #define iwtypein  0
        ! #endif 
        ! parameter(iwtype = iwtypein)
        integer,parameter :: iwtype_cnst = 0
        integer,parameter :: iwtype_pwcnst = 1
        integer,parameter :: iwtype_spwcnst = 2
        integer,parameter :: iwtype_flex = 3
        
#ifndef imixtype_background_in 
        integer,parameter :: imixtype_background_in = 1
#endif 
#ifndef imixtype_OM_in 
        integer,parameter :: imixtype_OM_in = 1
#endif 

        integer,parameter :: imixtype_nobio = 0
        integer,parameter :: imixtype_fick = 1
        integer,parameter :: imixtype_turbo2 = 2
        integer,parameter :: imixtype_till = 3
        integer,parameter :: imixtype_labs = 4

        integer,parameter :: iroughtype_smooth = 0
        integer,parameter :: iroughtype_NSB07 = 1
        integer,parameter :: iroughtype_BM00 = 2
        integer,parameter :: iroughtype_Letal21 = 3

        integer,parameter::nsp_sld_all = 80
        integer,parameter::nsp_aq_ph = 17
        integer,parameter::nsp_aq_all = 17
        integer,parameter::nsp_gas_ph = 2
        integer,parameter::nsp_gas_all = 4
        integer,parameter::nrxn_ext_all = 13

        ! an attempt to record psd
        integer,parameter :: nps = 50 ! bins for particle size 
        ! real(kind=8),parameter :: ps_min = 0.1d-6 ! min particle size (0.1 um)
        real(kind=8),parameter :: ps_min = 10d-9 ! min particle size (10 nm)
        ! real(kind=8),parameter :: ps_min = 100d-9 ! min particle size (100 nm)
        real(kind=8),parameter :: ps_max = 10d-3 ! max particle size (10 mm)
        real(kind=8),parameter :: pi = 4d0*atan(1d0) ! 
        ! real(kind=8),parameter :: psd_th_0 = 1d-12 ! 
        ! real(kind=8),parameter :: psd_th_0 = 1d-3 ! 
        real(kind=8),parameter :: psd_th_0 = 1d0 ! 
        real(kind=8),parameter :: tol_dvd = 1d-4 ! 
        ! integer,parameter :: nps_rain_char = 4
        integer,parameter :: nflx_psd = 6
#ifndef full_flux_report
        integer,parameter ::nsp_saveall = 1
#endif 
        integer,parameter ::idust = 15
        integer,parameter :: nph = 101

end module scepter_constants