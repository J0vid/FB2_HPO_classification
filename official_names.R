#official names script####
#key between fb2 names and hpo names####
# load("Classification_demo_legacy/demo_objects.Rdata")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# save(atlas, d.meta, front.face, PC.eigenvectors, PC.scores, synd.mshape, phenotype.df, hpo, hpo.pos, hdrda.mod, hdrda.df, official.names, file = "data.Rdata")
load("data.Rdata")

official.names <- levels(hdrda.df$synd)

# hdrda.df$synd[hdrda.df$synd == "Cleft_Lip_Palate"] <- "CLEFT LIP/PALATE"
#no 18p tetrasomy, goldenhar, klinefelter, xxyy, xxx, Trisomy 13,  in hpo?
#OMIM names####
official.names[official.names == "1p36 Deletion"] <- "CHROMOSOME 1P36 DELETION SYNDROME"
official.names[official.names == "22q 11.2 Deletion Syndrome"] <- "CHROMOSOME 22q11.2 DELETION SYNDROME, DISTAL"
official.names[official.names == "5p Deletion Cri du Chat"] <- "CRI-DU-CHAT SYNDROME"
official.names[official.names == "Achondroplasia"] <- "#100800 ACHONDROPLASIA; ACH"
official.names[official.names == "Angelman Syndrome"] <- "ANGELMAN SYNDROME"
official.names[official.names == "Apert Syndrome"] <- "#101200 APERT SYNDROME;;ACROCEPHALOSYNDACTYLY, TYPE I; ACS1;;ACS IAPERT-CROUZON DISEASE, INCLUDED;;ACROCEPHALOSYNDACTYLY, TYPE II, INCLUDED;;ACS II, INCLUDED;;VOGT CEPHALODACTYLY, INCLUDED"
official.names[official.names == "Cardiofaciocutaneous Syndrome"] <- "CARDIOFACIOCUTANEOUS SYNDROME 1; CFC1"
official.names[official.names == "CHARGE Syndrome"] <- "CHARGE SYNDROME"
official.names[official.names == "Cleft Lip/Palate"] <- "CLEFT LIP/PALATE WITH ABNORMAL THUMBS AND MICROCEPHALY"
official.names[official.names == "Cockayne Syndrome"] <- "COCKAYNE SYNDROME A; CSA"
official.names[official.names == "Coffin-Siris Syndrome"] <- "COFFIN-SIRIS SYNDROME 1; CSS1"
official.names[official.names == "Cohen Syndrome"] <- "COHEN SYNDROME"
official.names[official.names == "Cornelia de Lange Syndrome"] <- "CORNELIA DE LANGE SYNDROME 1"
official.names[official.names == "Costello Syndrome"] <- "COSTELLO SYNDROME; CSTLO"
official.names[official.names == "Crouzon Syndrome"] <- "CROUZON SYNDROME"
official.names[official.names == "Down Syndrome"] <- "DOWN SYNDROMETRISOMY 21"
official.names[official.names == "Ehlers Danlos Syndrome"] <- "EHLERS-DANLOS SYNDROME, TYPE I"
official.names[official.names == "Fragile X"] <- "FRAGILE X MENTAL RETARDATION SYNDROME"
official.names[official.names == "Jacobsen Syndrome"] <- "JACOBSEN SYNDROME"
official.names[official.names == "Joubert Syndrome"] <- "JOUBERT SYNDROME 1"
official.names[official.names == "Kabuki Syndrome"] <- "KABUKI SYNDROME 1"
official.names[official.names == "Loeys-Dietz Syndrome"] <- "LOEYS-DIETZ SYNDROME, TYPE 1A LOEYS-DIETZ AORTIC ANEURYSM SYNDROME"
official.names[official.names == "Marfan Syndrome"] <- "MARFAN SYNDROME; MFS"
official.names[official.names == "Moebius Syndrome"] <- "MOEBIUS SYNDROME; MBS"
official.names[official.names == "Mucopolysaccharidosis"] <- "MUCOPOLYSACCHARIDOSIS TYPE IIIA"
official.names[official.names == "Nager Syndrome"] <- "ACROFACIAL DYSOSTOSIS 1, NAGER TYPE; AFD1"
official.names[official.names == "Neurofibromatosis"] <- "NEUROFIBROMATOSIS, TYPE I"
official.names[official.names == "Noonan Syndrome"] <- "#163950 NOONAN SYNDROME 1; NS1;;NOONAN SYNDROME;;MALE TURNER SYNDROME;;FEMALE PSEUDO-TURNER SYNDROME;;TURNER PHENOTYPE WITH NORMAL KARYOTYPEPTERYGIUM COLLI SYNDROME, INCLUDED"
official.names[official.names == "Osteogenesis Imperfecta" | official.names == "Osteogenesis imperfecta"] <- "OSTEOGENESIS IMPERFECTA, TYPE I"
official.names[official.names == "Phelan McDermid Syndrome"] <- "PHELAN-MCDERMID SYNDROME; PHMDS"
official.names[official.names == "Pierre Robin Sequence"] <- "311895 PIERRE ROBIN SEQUENCE WITH FACIAL AND DIGITAL ANOMALIES"
official.names[official.names == "Pitt-Hopkins Syndrome"] <- "PITT-HOPKINS SYNDROME"
official.names[official.names == "Pseudoachondroplasia"] <- "PSEUDOACHONDROPLASIA"
official.names[official.names == "Rett Syndrome" | official.names == "Rett Syndrome_Other" | official.names == "Rett Syndrome_CDKL5"] <- "RETT SYNDROME; RTT"
official.names[official.names == "Rhizomelic Chondrodysplasia Punctata"] <- "#215100 RHIZOMELIC CHONDRODYSPLASIA PUNCTATA, TYPE 1; RCDP1;;PEROXISOME BIOGENESIS DISORDER 9; PBD9;;CHONDRODYSPLASIA PUNCTATA, RHIZOMELIC FORM; CDPR;;CHONDRODYSTROPHIA CALCIFICANS PUNCTATA"
official.names[official.names == "Rubinstein-Taybi Syndrome"] <- "RUBINSTEIN-TAYBI SYNDROME 1; RSTS1"
official.names[official.names == "Russell Silver Syndrome"] <- "SILVER-RUSSELL SYNDROME; SRS"
official.names[official.names == "Smith-Lemli-Opitz Syndrome"] <- "SMITH-LEMLI-OPITZ SYNDROME; SLOS"
official.names[official.names == "Smith-Magenis Syndrome"] <- "SMITH-MAGENIS SYNDROME; SMS"
official.names[official.names == "Sotos Syndrome"] <- "#117550 SOTOS SYNDROME 1; SOTOS1;;SOTOS SYNDROME;;CEREBRAL GIGANTISM;;CHROMOSOME 5q35 DELETION SYNDROME"
official.names[official.names == "Spondyloepiphyseal Dysplasia"] <- "SPONDYLOEPIPHYSEAL DYSPLASIA WITH CONGENITAL JOINT DISLOCATIONS"
official.names[official.names == "Stickler Syndrome"] <- "STICKLER SYNDROME, TYPE I"
official.names[official.names == "Treacher Collins Syndrome"] <- "TREACHER COLLINS-FRANCESCHETTI SYNDROME"
official.names[official.names == "Trisomy 18"] <- "TRISOMY 18-LIKE SYNDROME"
official.names[official.names == "Turner Syndrome"] <- "MENTAL RETARDATION, X-LINKED, SYNDROMIC, TURNER TYPE; MRXST"
official.names[official.names == "Van der Woude Syndrome"] <- "VAN DER WOUDE SYNDROME"
official.names[official.names == "Williams-Beuren Syndrome"] <- "WILLIAMS-BEUREN SYNDROME; WBS"
official.names[official.names == "X-Linked Hypohidrotic Ectodermal Dysplasia"] <- "ECTODERMAL DYSPLASIA 1, HYPOHIDROTIC, X-LINKED; XHED"
official.names[official.names == "Craniofrontonasal Dysplasia"] <- "CRANIOFRONTONASAL SYNDROME; CFNS"
official.names[official.names == "Muenke Syndrome"] <- "Muenke syndrome"
official.names[official.names == "Craniosynostosis"] <- "CRANIOSYNOSTOSIS 1"
official.names[official.names == "Pfeiffer Syndrome"] <- "#101600 PFEIFFER SYNDROME;;ACROCEPHALOSYNDACTYLY, TYPE V; ACS5;;ACS V;;NOACK SYNDROMECRANIOFACIAL-SKELETAL-DERMATOLOGIC DYSPLASIA, INCLUDED"
official.names[official.names == "Ectodermal Dysplasia"] <- "FACIAL ECTODERMAL DYSPLASIA"
official.names[official.names == "Autism Spectrum Disorder"] <- "INTELLECTUAL DEVELOPMENTAL DISORDER WITH AUTISM AND SPEECH DELAY; IDDAS"
official.names[official.names == "Wolf-Hirschhorn syndrome"] <- "#194190 WOLF-HIRSCHHORN SYNDROME; WHS;;CHROMOSOME 4p16.3 DELETION SYNDROME;;PITT-ROGERS-DANKS SYNDROME; PRDS;;PITT SYNDROME"
official.names[official.names == "Simpson-Golabi-Behmel Syndrome"] <- "SIMPSON-GOLABI-BEHMEL SYNDROME, TYPE 1"
official.names[official.names == "Microtia"] <- "MICROTIA-ANOTIA"
official.names[official.names == "Gaucher Disease"] <- "GAUCHER DISEASE, TYPE I"
official.names[official.names == "Prader-Willi Syndrome"] <- "PRADER-WILLI SYNDROME; PWS"
official.names[official.names == "Polycystic Kidney Disease"] <- "POLYCYSTIC KIDNEY DISEASE 4 WITH OR WITHOUT POLYCYSTIC LIVER DISEASE; PKD4"
official.names[official.names == "Saethre-Chotzen Syndrome"] <- "#101400 SAETHRE-CHOTZEN SYNDROME; SCS;;ACROCEPHALOSYNDACTYLY, TYPE III; ACS3;;ACS III;;CHOTZEN SYNDROME;;ACROCEPHALY, SKULL ASYMMETRY, AND MILD SYNDACTYLYSAETHRE-CHOTZEN SYNDROME WITH EYELID ANOMALIES, INCLUDED;;BLEPHAROPHIMOSIS, EPICANTHUS INVERSUS, AND PTOSIS 3, FORMERLY, INCLUDED;BPES3, FORMERLY, INCLUDED"
official.names[official.names == "Branchio-Oto-Renal Syndrome"] <- "BRANCHIOOTORENAL SYNDROME 1"
official.names[official.names == "Beckwith-Weidemann Syndrome"] <- "BECKWITH-WIEDEMANN SYNDROME"
official.names[official.names == "15q26.3 Deletion"] <- "#612626 CHROMOSOME 15q26-qter DELETION SYNDROME;;DRAYER SYNDROME"
official.names[official.names == "Diastrophic Dysplasia"] <- "#222600 DIASTROPHIC DYSPLASIA;;DTD;;DDDIASTROPHIC DYSPLASIA, BROAD BONE-PLATYSPONDYLIC VARIANT, INCLUDED"
official.names[official.names == "Ectrodactyly-Ectodermal Dysplasia-Cleft Lip/Palate"] <- "%129900 ECTRODACTYLY, ECTODERMAL DYSPLASIA, AND CLEFT LIP/PALATE SYNDROME1; EEC1;;EEC;;EEC SYNDROME 1"
official.names[official.names == "Goltz Syndrome"] <- "FOCAL DERMAL HYPOPLASIA; FDH"
official.names[official.names == "18p Deletion"] <- "CHROMOSOME 18P DELETION SYNDROME"
official.names[official.names == "15q Duplication"] <- "CHROMOSOME 15q11-q13 DUPLICATION SYNDROME"
official.names[official.names == "Opitz GBBB Syndrome"] <- "OPITZ GBBB SYNDROME, TYPE I; GBBB1"
official.names[official.names == "Hemifacial Microsomia"] <- "HEMIFACIAL MICROSOMIA"
official.names[official.names == "Epileptic Encephalopathy Early Infantile Type 2"] <- "#300672 EPILEPTIC ENCEPHALOPATHY, EARLY INFANTILE, 2; EIEE2;;INFANTILE SPASM SYNDROME, X-LINKED 2; ISSX2"
official.names[official.names == "Mowat-Wilson Syndrome"] <- "MOWAT-WILSON SYNDROME; MOWS"
official.names[official.names == "Zellweger Syndrome"] <- "PEROXISOME BIOGENESIS DISORDER 1A (ZELLWEGER)"
official.names[official.names == "Coffin-Lowry Syndrome"] <- "COFFIN-LOWRY SYNDROME; CLS"
official.names[official.names == "Bardet-Biedl Syndrome"] <- "BARDET-BIEDL SYNDROME 1; BBS1"
official.names[official.names == "Pallister-Killian Syndrome"] <- "PALLISTER-KILLIAN SYNDROME"


#orpha names####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

official.names <- levels(hdrda.df$synd)

#match to new HPO database to get better frequency info
phenotype_2022 <- readr::read_delim("phenotype-2022.hpoa", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 4)

#cut down to only physical abnormalities
phenotype_2022 <- phenotype_2022[phenotype_2022$Aspect == "P", ]

View(phenotype_2022[,c(1,2,4,8)])

#no goldenhar, klinefelter, in hpo?
official.names[official.names == "18p Tetrasomy"] <- "Tetrasomy 18p"
official.names[official.names == "XXYY"] <- "48,XXYY syndrome"
official.names[official.names == "XXX"] <- "Trisomy X" 
official.names[official.names == "Trisomy 13"] <- "Trisomy 13" 
official.names[official.names == "1p36 Deletion"] <- "Chromosome 1p36 deletion syndrome"
official.names[official.names == "22q 11.2 Deletion Syndrome"] <- "Chromosome 22q11.2 deletion syndrome, distal"
official.names[official.names == "5p Deletion Cri du Chat"] <- "Cri du Chat Syndrome (5p deletion)"
official.names[official.names == "Achondroplasia"] <- "Achondroplasia"
official.names[official.names == "Angelman Syndrome"] <- "Angelman syndrome"
official.names[official.names == "Apert Syndrome"] <- "Apert syndrome"
official.names[official.names == "Cardiofaciocutaneous Syndrome"] <- "Cardiofaciocutaneous syndrome"
official.names[official.names == "CHARGE Syndrome"] <- "CHARGE syndrome"
official.names[official.names == "Cleft Lip/Palate"] <- "Cleft lip/palate"
official.names[official.names == "Cockayne Syndrome"] <- "Cockayne syndrome"
official.names[official.names == "Coffin-Siris Syndrome"] <- "Coffin-Siris syndrome"
official.names[official.names == "Cohen Syndrome"] <- "Cohen syndrome"
official.names[official.names == "Cornelia de Lange Syndrome"] <- "Cornelia de Lange syndrome"
official.names[official.names == "Costello Syndrome"] <- "Costello syndrome"
official.names[official.names == "Crouzon Syndrome"] <- "Crouzon syndrome"
official.names[official.names == "Down Syndrome"] <- "Down syndrome"
official.names[official.names == "Ehlers Danlos Syndrome"] <- "Classical Ehlers-Danlos syndrome"
official.names[official.names == "Fragile X"] <- "Fragile X syndrome"
official.names[official.names == "Jacobsen Syndrome"] <- "Jacobsen syndrome"
official.names[official.names == "Joubert Syndrome"] <- "Joubert syndrome"
official.names[official.names == "Kabuki Syndrome"] <- "Kabuki syndrome"
official.names[official.names == "Loeys-Dietz Syndrome"] <- "Loeys-Dietz syndrome"
official.names[official.names == "Marfan Syndrome"] <- "Marfan syndrome"
official.names[official.names == "Moebius Syndrome"] <- "Moebius syndrome"
official.names[official.names == "Mucopolysaccharidosis"] <- "Mucopolysaccharidosis type 1"
official.names[official.names == "Nager Syndrome"] <- "Nager syndrome"
official.names[official.names == "Neurofibromatosis"] <- "Neurofibromatosis type 1"
official.names[official.names == "Noonan Syndrome"] <- "Noonan syndrome"
official.names[official.names == "Osteogenesis Imperfecta" | official.names == "Osteogenesis imperfecta"] <- "Osteogenesis imperfecta"
official.names[official.names == "Phelan McDermid Syndrome"] <- "Phelan-Mcdermid syndrome"
official.names[official.names == "Pierre Robin Sequence"] <- "Contractures-developmental delay-Pierre Robin syndrome"
official.names[official.names == "Pitt-Hopkins Syndrome"] <- "Pitt-Hopkins syndrome"
official.names[official.names == "Pseudoachondroplasia"] <- "Pseudoachondroplasia"
official.names[official.names == "Rett Syndrome" | official.names == "Rett Syndrome_Other" | official.names == "Rett Syndrome_CDKL5"] <- "Rett syndrome"
official.names[official.names == "Rhizomelic Chondrodysplasia Punctata"] <- "Rhizomelic chondrodysplasia punctata"
official.names[official.names == "Rubinstein-Taybi Syndrome"] <- "Rubinstein-Taybi syndrome"
official.names[official.names == "Russell Silver Syndrome"] <- "Silver-Russell syndrome"
official.names[official.names == "Smith-Lemli-Opitz Syndrome"] <- "Smith-Lemli-Opitz syndrome"
official.names[official.names == "Smith-Magenis Syndrome"] <- "Smith-Magenis syndrome"
official.names[official.names == "Sotos Syndrome"] <- "Sotos syndrome"
official.names[official.names == "Spondyloepiphyseal Dysplasia"] <- "Spondyloepiphyseal dysplasia congenita"
official.names[official.names == "Stickler Syndrome"] <- "Stickler syndrome"
official.names[official.names == "Treacher Collins Syndrome"] <- "Treacher-Collins syndrome"
official.names[official.names == "Trisomy 18"] <- "Trisomy 18"
official.names[official.names == "Turner Syndrome"] <- "Turner syndrome"
official.names[official.names == "Van der Woude Syndrome"] <- "Van der Woude syndrome"
official.names[official.names == "Williams-Beuren Syndrome"] <- "Williams syndrome"
official.names[official.names == "X-Linked Hypohidrotic Ectodermal Dysplasia"] <- "X-linked hypohidrotic ectodermal dysplasia"
official.names[official.names == "Craniofrontonasal Dysplasia"] <- "Craniofrontonasal dysplasia"
official.names[official.names == "Muenke Syndrome"] <- "Muenke syndrome"
official.names[official.names == "Craniosynostosis"] <- "Craniosynostosis and dental anomalies"
official.names[official.names == "Pfeiffer Syndrome"] <- "Pfeiffer syndrome"
official.names[official.names == "Ectodermal Dysplasia"] <- "Cranioectodermal dysplasia"
official.names[official.names == "Autism Spectrum Disorder"] <- "Autism spectrum disorder due to AUTS2 deficiency"
official.names[official.names == "Wolf-Hirschhorn syndrome"] <- "Wolf-Hirschhorn syndrome"
official.names[official.names == "Simpson-Golabi-Behmel Syndrome"] <- "Simpson-Golabi-Behmel syndrome"
official.names[official.names == "Microtia"] <- "Microtia"
official.names[official.names == "Gaucher Disease"] <- "Gaucher disease"
official.names[official.names == "Prader-Willi Syndrome"] <- "Prader-Willi syndrome"
official.names[official.names == "Polycystic Kidney Disease"] <- "Autosomal recessive polycystic kidney disease"
official.names[official.names == "Saethre-Chotzen Syndrome"] <- "Saethre-Chotzen syndrome"
official.names[official.names == "Branchio-Oto-Renal Syndrome"] <- "Branchiootorenal syndrome 1"
official.names[official.names == "Beckwith-Weidemann Syndrome"] <- "Beckwith-Wiedemann syndrome"
official.names[official.names == "15q26.3 Deletion"] <- "Chromosome 15q26-qter deletion syndrome"
official.names[official.names == "Diastrophic Dysplasia"] <- "Diastrophic dysplasia"
official.names[official.names == "Ectrodactyly-Ectodermal Dysplasia-Cleft Lip/Palate"] <- "Ectrodactyly, ectodermal dysplasia, and cleft lip/palate syndrome1"
official.names[official.names == "Goltz Syndrome"] <- "Focal dermal hypoplasia"
official.names[official.names == "18p Deletion"] <- "Chromosome 18P deletion syndrome"
official.names[official.names == "15q Duplication"] <- "Chromosome 15q11-q13 duplication syndrome"
official.names[official.names == "Opitz GBBB Syndrome"] <- "Opitz GBBB syndrome"
official.names[official.names == "Hemifacial Microsomia"] <- "Hemifacial microsomia"
official.names[official.names == "Epileptic Encephalopathy Early Infantile Type 2"] <- "Early infantile epileptic encephalopathy"
official.names[official.names == "Mowat-Wilson Syndrome"] <- "Mowat-Wilson syndrome"
official.names[official.names == "Zellweger Syndrome"] <- "Zellweger syndrome"
official.names[official.names == "Coffin-Lowry Syndrome"] <- "Coffin-Lowry syndrome"
official.names[official.names == "Bardet-Biedl Syndrome"] <- "Bardet-Biedl syndrome"
official.names[official.names == "Pallister-Killian Syndrome"] <- "Pallister-Killian syndrome"













