# questions de bioinfo -> Elodie, Jennifer, Julie
~~* [ ] voir avec elodie la première formule de DNAcopy.~~
* [X] ~~*ASCAT : voir w/ elodie pour les équations du calling? pas forcément important pour nous. cf cahier 26/04*~~ [2022-05-04]
* [ ] Pour Rihab: comment appliquer le même style à plusieurs éléments du même panel? -> CSS 
* [ ] Pour élodie: comment dire à un panel de prendre la place restante de la fenêtre ? -> parameters de CGHcall.



# questions de bio -> Laetitia, Claire
- pour laetitia, cf cahier :
    - [X] ~~*peut-on avoir START et STOP pour des probeset.txt?*~~ [2022-04-06] a priori non, il faudra demander à Tony
    - [X] ~~*pourquoi a-t-on plus de segments (100 exactement) si on exporte segments.txt à partir de Workflow analysis par rapport à la vue Segments de ChAS?*~~ [2022-04-06] beaucoup de ces segments sont à 0 de log ratio, donc négligés par ChAS.
    * [X] ~~*pourquoi, dans ChAS, a-t-on parfois 24 chromosomes? est-ce le chromosome Y?*~~ [2022-04-27] chr24 = chrX et chr25 = chrY
    * [X] ~~*pourquoi ChAS n'exporte-t-il pas les segments des chromosomes sexuels?*~~ [2022-05-02] pour des raisons de confidentialité.
    * [X] ~~*Je ne peux pas générer l'OSCHP de 1-RV à partir des fichiers CEL, il y a sûrement un problème avec le fichier AT . Demander à laetitia si c'est possible de me le renvoyer*~~ [2022-06-28]
    * [X] ~~*Ce qu'on appelle la récurrence, c'est l'apparition de métastases?*~~ [2022-06-28] c'est la rechute
    * [X] ~~*L'estimation de la cellularité est-elle faite en utilisant une lame différente de celle utilisée pour l'analyse CGH?*~~ [2022-05-12] la lame est issue du bloc, donc les deux informations viennent du même endroit
    * [X] ~~*le logiciel propriétaire d'affymétrix est-il bien `OncoScan Console 1.3` ?*~~ [2022-05-23] Non, c'est ChAS.
    * [X] ~~*qu'est-ce qu'une perte hétérozygote/homozygote? C'est différent d'une LOH.*~~ [2022-05-31] 
# questions sur OncoscanCNV/ChAS/... -> Tony
* [X] ~~*peut-on exporter start et end en colonnes de la vue graph de ChAS? sinon, a-t-on un fichier qui indique la couverture de chaque sonde? c'est-à-dire les positions de début et de fin. car CGHcall a besoin de ces 2 infos en input, et ChAS ne les exporte pas a priori.*~~ [2022-05-13] a priori non mais cf cahier à la date d'aujourd'hui
* [X] ~~*Ou alors, a-t-on un fichier qui indique la couverture de chaque sonde, ou la taille de chaque sonde?*~~ [2022-05-13] ce sont des positions de SNP. cf cahier aussi.
* [X] ~~*peut-on visualiser la couverture des sondes par rapport au génome?*~~ [2022-05-13] -> oui car on a les positions de SNP pour chaque sonde.
* [ ] On peut exporter les données des segments à partir de la vue interne de ChAS, mais c'est fastidieux car il faut traiter chaque échantillon individuellement. Dans analysis setup, on peut générer en batch un fichier qui contient tous les segments des OSCHP qu'on lui donne -> plus rapide. J'aimerais bien pouvoir choisir les colonnes que je veux à ce niveau, est-ce possible? (je veux choisir les colonnes car pour OncoscanR, on a besoin de la colonne Full Location: c'est primordial.)
* [X] ~~*comment ChAS estime le copy number? modèle de mélange de gaussiennes?*~~ [2022-05-12]
* [ ] -> voir dans ChASRUO.pdf pour voir comment les segments sont estimés aussi.
* [ ] de combien est la compaction des log Ratios?
* [ ] est-ce que ChAS utilise le programme en ligne de commande Array Power Tools (APT, aussi appelé Affymetrix Power Tools ou Analysis Power Tools) pour déterminer le nombre de copies?

