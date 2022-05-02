# questions de bioinfo -> Elodie, Jennifer, Julie
* [ ] voir avec elodie la première formule de DNAcopy.
* [ ] ASCAT : voir w/ elodie pour les équations du calling? pas forcément important pour nous. cf cahier 26/04
# questions de bio -> Laetitia, Claire
- pour laetitia, cf cahier :
    - [X] ~~*peut-on avoir START et STOP pour des probeset.txt?*~~ [2022-04-06] a priori non, il faudra demander à Tony
    - [X] ~~*pourquoi a-t-on plus de segments (100 exactement) si on exporte segments.txt à partir de Workflow analysis par rapport à la vue Segments de ChAS?*~~ [2022-04-06] beaucoup de ces segments sont à 0 de log ratio, donc négligés par ChAS.
    * [X] ~~*pourquoi, dans ChAS, a-t-on parfois 24 chromosomes? est-ce le chromosome Y?*~~ [2022-04-27] chr24 = chrX et chr25 = chrY
    * [ ] pourquoi ChAS n'exporte-t-il pas les segments des chromosomes sexuels?
# questions sur OncoscanCNV/ChAS/... -> Tony
* [ ] peut-on exporter start et end en colonnes de la vue graph de ChAS? sinon, a-t-on un fichier qui indique la couverture de chaque sonde? c'est-à-dire les positions de début et de fin. car CGHcall a besoin de ces 2 infos en input, et ChAS ne les exporte pas a priori.
* [ ] peut-on visualiser la couverture des sondes par rapport au génome?
* [ ] On peut exporter les données des segments à partir de la vue interne de ChAS, mais c'est fastidieux car il faut traiter chaque échantillon individuellement. Dans l'analysis setup, on peut générer en batch un fichier qui contient tous les segments des OSCHP qu'on lui donne -> plus rapide. J'aimerais bien pouvoir choisir les colonnes que je veux à ce niveau, est-ce possible? (je veux choisir les colonnes car pour OncoscanR, on a besoin de la colonne Full Location: c'est primordial.)
* [ ] comment ChAS estime le copy number? modèle de mélange de gaussiennes?