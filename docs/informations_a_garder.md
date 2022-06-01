# interprétation des résultats
## 05/05/2022
CGHcall: la colonne endpos est construite à partir de startpos+20. garder ça en tête lors de l'interprétation des résultats.

## 06/05/2022
rCGH: peut paralléliser les tâches, mais pas sur windows -> potentiellement plus rapide. revoir quelles tâches.

oncoscanR transforme les données: nombre de segments -> pourcentage d'altération par bras. Il réduira donc le nombre de segments. Ce n'est pas forcément inintéressant, car c'est un filtrage qui pourrait avoir du sens -> à tester.

## 09/05/2022
- La segmentation DNAcopy (CGHcall, rCGH) ne tient pas compte du BAF; ASPCF, si. On peut donc rater des segments de CN-NE avec DNAcopy, pas ASPCF.
- CGHcall: Les données brutes des sondes comprennent 360 sondes sans valeur (NaN) sur les 200 000 et quelques. Elles sont supprimées par la fonction preprocess() de CGHcall, qui effectue d'autres altérations des données, notamment combler les trous laissés par les données manquantes en estimant leur valeur. Tous les échantillons testés présentent un NaN pour ces sondes, le problème vient donc de la construction des données.
- CGHcall: pour le calling, on peut spécifier un nombre de CPUs (installer le package snowfall est alors requis).
- CGHcall: pour le calling, on peut calculer les prior probabilities des segments sur tout le génome ou sur chaque bras chromosomique.
- [X] ~~*Les informations des chromosomes sexuels ne sont pas utilisées par oncoscanR, mais bien par les autres. c'est à savoir pour comparer les GIs.*~~ [2022-05-23] Les infos des chromosomes sexuels sont surtout à retirer des autres outils car on ne doit pas calculer le GI avec.
## 10/05/2022
- CGHcall: les valeurs de call sont en ``log ratio``, donc les plots produits ne sont pas à interpréter en copy number mais bien en LRR. Le Gi est calculé en prenant ça en compte, mais ce serait mieux de normaliser l'output pour qu'il ait le même que rCGH. Il suffit de faire +2 à toutes les valeurs.'

## 16/05/2022
- Parmi les outils à considérer au-delà de ceux sur lesquels j'ai travaillé, il y a CGHclassify (cité dans l'article de CGHcall) et GISTIC (cité dans la revue comparative de 6 outils) `edit: le lien menant au code de CGHclassify est cassé.`
- CGHcall permet d'utiliser non pas les segments mais les bras chromosomiques dans le modèle de mélange. à tester pour calculer le GI! -> et voir la différence de GI obtenu.

## 18/05/2022
* [ ] rCGH: pour voir si l'estimation du nombre de copies est pertinente, comparer le plot des segments bruts et des segments callés. Si ça remanie trop, réfléchir à filtrer à la main les altérations pertinentes. Trouver les cas où l'estimation n'est pas pertinente et en parler dans le rapport.

## 19/05/2022
Un Hidden Markov Model peut être dérivé d'un modèle de mélange (cf. `https://en.wikipedia.org/wiki/Mixture_model#Topics_in_a_document`).

## 20/05/2022
La citation du site d'Agilent Sureprint a un problème visuel de parenthèses quand elle apparaît dans la biblio à la fin du rapport.
* [X] ~~*ajouter l'info que le calcul du GI sur Agilent est actuellement fait manuellement.*~~ [2022-05-24]
* [ ] Résultats / Discussion: cf cahier bleu au 19/05. plusieurs choses à tester.
* [ ] rapport: questions à poser , choses à discuter: cf. cahier bleu au 19/05.
* [X] ~~*rapport: oncoscanR ne fait pas de lissage, pourtant c'est ce que la figure pipeline dit. supprimer ça de la figure.*~~ [2022-05-23] SI, oncoscanR fait un lissage (smoothing) en fusionnant les segments distants de moins de 300kbp.
* [X] ~~*noter dans fonctionnement_des_outils.md la licence de chaque package: peut-on en faire ce qu'on veut? -> ASCAT & EaCoN*~~ [2022-05-30] Oui, on ne peut juste pas modifier le code de rCGH.
* [X] ~~*rCGH: changer le pipeline de rCGH*~~ [2022-05-23] Pourquoi? comment?
* [ ] Etat de l'art, pipeline CGHcall: ajouter "CBS" à la segmentation pour être raccord avec le pipeline rCGH

## 23/05/2022
* [ ] oncoscanR: faire un plot qui montre les segments d'altération sans les superposer -> on peut voir si un LOH est supérieur à un segment de perte, et surtout voir l'impact du filtrage sur les segments. cf. cahier
* [ ] $vérifier que le calcul du GI d'oncoscanR est correct: $
    * [X] ~~*1. on doit calculer si le bras est en gain: somme des segments gain et amp*~~ [2022-05-30]
    * [X] ~~*2. on doit calculer si le bras est en perte: somme des segments loss et LOH*~~ [2022-05-30]
    * [X] ~~*on a donc un vecteur de bras en gain, et un vecteur de bras en perte*~~ [2022-05-30]
    * [ ] 3. on itère ces 2 vecteurs: pour chaque bras, si il subit une altération (qu'elle soit un gain ou une perte), on ajoute 1 dans le compte des altérations et on ajoute son chromosome dans la liste des chromosomes altérés.
    * [ ] enfin on calcule le GI à partir du compte des altérations et des chromosomes.
    bonus: faire un plot des segments de chaque altération: afficher |---| pour chaque segment
* [ ] le score TDplus d'oncoscanR est implémenté dans oncoscanR comme le nombre de segments en *gain* dont la longueur est comprise entre 1 et 10 Mb. Mais on n'a pas de garantie que ces gains sont bien issus de duplications en tandem.
* [ ] Exploiter le tableau comparatif pour dire ce que chaque outil ne fait *pas*. -> discussion?
* [ ] mettre les 4 organigrammes dans une image? et faire des symboles pour chaque étape / montrer un plot des donneés obtenues à chaque étape
* [ ] Si le temps le permet, regarder quel package est utilisé dans quel cancer -> regarder où l'article de tel outil est cité. pour faire une ouverture, à la fin du rapport.
* [ ] sortir les profils WGV de chaque échantillon avec ChAS pour pouvoir comparer avec les plots de segtables des 4 outils
* [X] ~~*`résultats: ` utiliser le cahier bleu au 10/05*~~ [2022-05-31]
* [ ] `résultats: ` pour tous les outils, montrer l'effet des nettoyages en faisant tourner les données avec et sans.
* [ ] oncoscanR: faire un plot qui montre le WGV d'un échantillon et indique les bras déclarés comme altérés et les bras non déclarés comme altérés. ça permet de visualiser l'effet du seuil et de *voir* ce qui est fait par cet outil. Voir powerpoint oncoscanR, diapo 13
* [ ] rCGH: regarder comment le nombre de copies est estimé. voir rCGH_dev.R .
* [ ] etat de l'art: rCGH: dans le texte, préciser que la normalisation et l'estimation du nombre de copies sont 2 étapes différentes et les séparer également dans la fig pipeline. et aussi regarder pourquoi l'estimation du nb copies se fait pdt la seg *et* pdt la normalisation.
* [ ] Mat Met , rCGH, seg CBS: détailler plus (p-value, distribution de référence...)
* [ ] seg DNA copy dans matmet: revoir comment la taille de la fenetre coulissante est utilisée.
* [ ] matmet rCGH normalisation: changer le plot pour qu'il corresponde à l'échantillon 5 ! je fais le reste des plots avec.
* [ ] matmet cghcall preprocess: utiliser 5LD
* [ ] matmet cghcall normalisation 2: utiliser 5LD.
* [ ] matmet cghcall preprocess: voir l'estimation des valeurs manquantes par impute.knn . mais plus tard.
## 24/05/2022
* [ ] résultats ASCAT: le critère d'optimization O est O = sum(gof_baf)\*(1-w) + sum(gof_lrr)\*w + penalty\*Q où Q est le nombre de segments. Par défaut penalty=50 (dans EaCoN), mais c'est 70 dans ASCAT par défaut. tester la valeur de 70, voir si ça change des choses notamment pour l'échantillon 17 qui avait un GI de 400+.
*Pendant la réunion d'aujourd'hui, j'ai pris ces notes:*
* [X] ~~*GI Agilent: l'ajouter aux plots de corrélation des GIs.*~~ [2022-05-30]
* [X] ~~*faire le plot comme indiqué par Claire. cf cahier bleu.*~~ [2022-05-30]
* [ ] Dans le rapport, afficher tous les pipelines horizontalement, et éventuellement en annexe placer les images des pipelines à la verticale + les plots de données qui correspondent à chaque étape.
* [ ] Rapport, Res & Discussion, perspectives: parler du fait qu'on peut changer le seuil de 300 kbp d'oncoscanR pour jouer sur le calcul du GI.
## 30/05/2022
* [X] ~~*ASCAT: renseigner la cellularité et voir si on obtient de meilleurs résultats*~~ [2022-05-30] Non. Mais comparer cette cellularité avec la cellularité estimée par lame HES.
* [ ] ASCAT.R: pourquoi y a-t-il des problèmes dans la sauvegarde des segTables?
## 31/05/2022
* [ ] Parler de ça dans la discussion: Laetitia m'indique que certaines altérations font parfois une dizaine de sondes. C'est à garder en tête pour discuter d'oncoscanR qui estime les altérations sur un bras entier, ou tous les outils qui lissent avant de déterminer les altérations.
`dans les scripts: `
* [ ] relancer tous calculs GI avec 1-RV.
* [ ] 3 grp -> 2 grp. on enlève l'intermédiaire.
* [ ] refaire plot distri avec 2 grp
* [ ] refaire plot cellularité: mêmes échelles en X et Y; vérifier les valeurs.
* [ ] vérifier qu'ASCAT donne bien le % de cellules tumorales
* [ ] courbes ROC: regarder packages R, pROC
* [ ] plot distri: relier tous les points. Ensuite, voir si c'est plus pertinent de le faire seulement avec certains points.
* [ ] changer, dans le tableur excel, les GI agilent que j'utilise pour les E 3 et 18. cf cahier.
* [X] ~~*Demander à Gaetan une photo de lame HES cerclée ou non mais scannée sur laquelle on voit facilement les cellules tumorales ou non tumorales.*~~ [2022-05-31] Claire a des photos de ce genre.
* [ ] illustrer le plot de cellularité par une photo de lame HES cerclée; voir avec Claire.
* [ ] retirer E 17 des plots et résultats
* [ ] Mat Met ASCAT: changer le plot de distance entre les segments et les entiers non nuls les plus proches. cf cahier au 23/05
* [ ] Mat Met ASCAT: eventuellement ajouter un plot qui résume comment on peut déterminer la ploïdie (ou le CN) et la cellularité à partir de log RAtio et BAF. cf cahier au 23/05, ou slides ASCAT 45-47 (donc dans `C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\docs\docs_I_made\images`) ou à partir de la vidéo `https://www.biodiscovery.com/videos/ascat-algorithm`
* [ ] résultats: utiliser les bons GI d'Agilent pour les échantillons qui en ont 2. cf cahier.
* [ ] ajouter un label à chaque point du plot distri.
* [X] ~~*CGHcall ## make it so a list containing only 1 cghCall result can't exist. do this at the end of pipelineCGHcall(). cf getPrbLvSegments()*~~ [2022-06-01]