#####prior in SPSS:
# Descriptive statistics for Liking_1, Liking_2, Meta_Liking_1, Meta_Liking_2, SM_Global rater1, rater2 and rater3
# Conduct Cohens Kappa, Cronbachs Alpha & ICC Index for interraterreliability 


library (foreign)
library (TripleR)
# Load neccessary packages

RRSM <- read.spss ("Raw_data_Ratings_Social_Mimicry.sav", to.data.frame=TRUE)
# Define an object with data set


####################################################SRM-ANalysen#################################################

#Univariate SRM-ANalysen (Table 1 Main Text; Table A2 Supplemenatry Analyses)

SMfactorcompounds.latent <- RR(Zsocial_mimicry_global_rater1_rater2_rater3/Zsocial_mimicry_actions_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)
SMfactorcompounds.latent
#Number of Actions & Global Social Mimicry rating building one "SM-Faktor" -> Results of Table 1

SMZ_Mean_GlobalandActions.compounds<- RR(ZSM_Mean_GlobalandActions ~ subject * target | group, data=RRSM)
SMZ_Mean_GlobalandActions.compounds
# Z-standardized Social Mimicry Mean Global and Number of Actions -> das gleiche manifest zum Auslesen der Werte

SMcompounds<- RR(social_mimicry_global_rater1_rater2_rater3 ~ subject * target | group, data=RRSM)
SMcompounds
# SM Global rating rater1, rater2 and rater3 -> Results of Table A2 "Global Rating"

SMactioncompounds.rater1_rater2_rater3 <- RR(social_mimicry_actions_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)
SMactioncompounds.rater1_rater2_rater3
# Number of SM Actions rater1, rater2 and rater3 -> Results of Table A2 "Number of actions"

Liking_1compounds <- RR(liking_1 ~ subject * target | group, data=RRSM)
Liking_1compounds
# Liking pre dyadic interaction general -> Results of Table A2 "Liking Pre-Conversation" (Achtung Tipp-Fehler in der Tabelle bei der Perceiver Variance)

Liking_2compounds <- RR(liking_2 ~ subject * target | group, data=RRSM)
Liking_2compounds
# Liking post dyadic interaction global -> Results of Table A2 "Liking Post-Conversation"

Metaliking_1compounds<- RR(meta_liking_1 ~ subject * target | group, data=RRSM)
Metaliking_1compounds
# Metaliking pre dyadic interaction global -> Results of Table A2 "Meta-Liking Pre-Conversation"

Metaliking_2compounds <- RR(meta_liking_2 ~ subject * target | group, data=RRSM)
Metaliking_2compounds
# Metaliking post dyadic interaction global -> Results of Table A2 "Meta-Liking Post-Conversation"




####################################################Export Process für das Pfadmodell#################################################

### Export relationship effects  (saved as .dat data):

write.csv(SMcompounds$effectsRel, "Relationshipeffect data/Relationsship_SMglobal_rater123.dat")
#Social Mimicry global rating rater1/rater2/rater3

write.csv(SMactioncompounds.rater1_rater2_rater3$effectsRel, "Relationshipeffect data/Relationsship_SMActions_rater123.dat")
# Number of SM Actions rater1/rater2/rater3

write.csv(SMZ_Mean_GlobalandActions.compounds$effectsRel, "Relationshipeffect data/Relationsship_SMZ_Mean_GlobalandActions.compounds.dat")
#...

write.csv(Liking_1compounds$effectsRel, "Relationshipeffect data/Relationsship_L1_lik_1.dat")
# Liking pre dyadic interaction

write.csv(Liking_2compounds$effectsRel, "Relationshipeffect data/Relationsship_L2_lik_2.dat")
# Liking post dyadic interaction

write.csv(Metaliking_1compounds$effectsRel, "Relationshipeffect data/Relationsship_ML1_metalik1.dat")
# Metaliking pre dyadic interaction

write.csv(Metaliking_2compounds$effectsRel, "Relationshipeffect data/Relationsship_ML2_metalik2.dat")
# Metaliking post dyadic interaction

#Wir haben damit  dann ein Mixed Model in SPSS gerechnet nach Vorlage aus dem Kenny-Buch: 
#
#********Liking Modell with control for Liking pre interaction* -> entspricht Figure 2 linke Seite.
#
#Hinweis: PersonA = Intrapersonal("doing mimicry") und PersonB = Interpersonal("being mimicked"), die Variablen wurden entsprechend "gekreuzt" im Datensatz kodiert
#
#MIXED
#Zrelationship_liking2_personA with Zrelationship_ZSM_Mean_GlobalandActions_personA Zrelationship_ZSM_Mean_GlobalandActions_personB Zrelationship_liking1_personA by dyad
#/fixed = Zrelationship_ZSM_Mean_GlobalandActions_personA Zrelationship_ZSM_Mean_GlobalandActions_personB Zrelationship_liking1_personA
#/print = solution testcov
#/random = dyad.

#und dann haben wir noch entsprechend ein paralleles Modell für Meta-Liking gerechnet (Figure 1, rechte Seite)


####################################################Export Process für die Korrelation mit Persönlichkeit auf individueller Ebene#################################################

### Export actor effects (saved as .dat data):
# -> export of SM Actor effects for correlations with Big5

write.csv(SMcompounds$effects, "Actoreffect data/Actor_SMglobal.dat")
write.csv(SMactioncompounds.rater1_rater2_rater3$effects, "Actoreffect data/Actor_SMactions.dat")
write.csv(SMZ_Mean_GlobalandActions.compounds$effects, "Actor_SMZ_Mean_GA.dat")
# es wurden dann in SPSS Partialkorrelationen mit den Big Five (Mittel aus Selbst und Fremd) unter Kontrolle von dummy-kodierten Gruppenvariablen gerechnet
#-> diese Ergebnisse finden sich dann in Table A3

#*calculate mean variables for self and third party ratings:.
#compute extra = mean (self_extra, F_extra_5B).
#EXECUTE.

#*Correlation SM(Global + Actions) and extraversion.
#PARTIAL CORR
#/VARIABLES=ZSM_Mean_GlobalandActions_actor extra BY d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17 d18 d19 d20 d21 d22 d23 d24 d25
#/SIGNIFICANCE=TWOTAIL
#/MISSING=LISTWISE.

