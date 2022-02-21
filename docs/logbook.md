# <span style="color:#999900"> 14/02/2022
tomorrow: meeting on the planning of the internship

I received pdfs from claire:
- Projet_GIRONDE_synopsis23052019.pdf
- 2021_Christinat_HRD_Oncoscan.pdf
- modpathol20153.pdf
- Nanocind_signature_S._CROCE.pdf

See notes_on_articles.md .\
About the synopsis:
## Projet_GIRONDE_synopsis23052019.pdf
le pronostic des lésions musculaires tumorales peut passer par des techniques comme l'index génomique (GI) = degré de complexité moléculaire & d'instabilité génomique. bon prédicteur de l'agressivité d'une tumeur.
cela est fait à partir de la technologie Agilent, mais l'équipe veut étendre ça à la technologie Oncoscan, plus récente. pour cela, il est important de comparer ces deux technologies sur cette même technique.
Cela est pertinent :
- en terme de diagnostique, car certaines tumeurs ne peuvent pas être traitées par les méthodes morphologiques.
- en terme de pronostique et de thérapeutique, cette méthode pourrait être appliquée à différentes tumeurs, même si on travaillera dans un premier temps sur les tumeurs stromales gastro-intestinales (GIST).

obj: transposer le calcul du GI et détermination du seuil de classification des tumeurs from Agilent to Affymetrix.
Oncoscan a une résolution et une couverture génomique plus large. cela devrait permettre une plus grande précision dans le diagnostique.


## Suivre la manip cette semaine
Remise en contexte: qu'est-ce que cette manip cherche à faire?
-> produire un index génomique
    -> comment est calculé l'index génomique?
    GI = A^2/C
    où A = total number of alterations (segmental gains and losses)
    et C = the number of involved chromosomes.
    -> C est déterminé par l'expérience
    -> A : voir les résultats.

15:50
J'ai suivi la manip avec Laetitia. Il semble qu'on fasse du SNP genotyping, et pas du Loci Capture.
Cette image illustre la manipulation dans son ensemble:
![](./oncoscan__MIP_probe.png)
D'où viennent les échantillons liquides d'ADN?
Laetitia reçoit des lames colorées Hématoxyline Éosine Safran (HES) + un FFPE: bloc de cellules tumorales extrait chirurgicalement. à l'aide de la lame (issue du bloc), elle sait où récupérer les cellules les plus tumorales du bloc. l'ADN est ensuite extrait de ces cellules sous forme liquide.
cet ADN est appelé l'ADNg dans le protocole. Après avoir ajouté des sondes MIP, on le chauffe à 95°C pour passer de double à simple brin, puis on descend à 58°C pendant 2h pour laisser les sondes MIP s'associer aux brins d'ADN. Cette association est appelée Annealing. La structure en anneau ainsi formée est centrée sur un nucléotide du brin d'ADN, le seul qui n'est pas couvert par la sonde.
l'étape suivante est le gap filling. on sépare le résultat de l'annealing en 2: un tube recevra du AT, l'autre du CG. les gaps seront complétés de manière complémentaire. ultimement, cela nous apprendra... qu'est-ce que ça nous apprend?
Bref, une exonucléase est ensuite ajoutée pour dégrader les brins d'ADN libres (cela inclut les sondes non accrochées) et ne garder ainsi que les anneaux. ces derniers sont ensuite clivés, ce qui donne la structure suivante:

%%%%###########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Où:\
%%%%% = Site 1 de clivage\
\##### = tag correspondant à la séquence d'ADN\
@@@@@ = région homologue à la séquence d'ADN visée.

Deux PCR sont ensuite effectuées pour multiplier les sondes, puis une digestion Hae III est appliquée pour cliver la séquence ADN du reste de la sonde:  

%%%%###########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

L'hybridation sur puce affymetrix a ensuite lieu. Le tag effectue cette hybridation. Ce dernier est spécifique de la région d'ADNg. Ainsi, le nombre de copies de ce tag permet de connaître le nombre de copies de la région à laquelle il correspond. Le statut allélique de ces copies pourra également être déterminé en fonction de si ces copies sont présentes sur la puce AT ou la puce GC (car on utilise une puce par paire de bases.)


## TODO:
Comment fonctionne Oncoscan CNV? Je veux savoir:
- sous quelle forme apparaissent les résultats
- quel est le but de la méthode?
- qu'est-ce que le gap filling nous apprend?







Questions:
- l'annealing laisse la place à un nucléotide. pourtant plus tard on vient boucher les trous avec des AT ou des CG, spécifiquement. quid?
- pourquoi ne faire du gap filling que sur AT et GC quand on pourrait le faire sur les 4 nucléotides? parce qu'ils sont complémentaires, mais ça ne répond que partiellement à la question. En fait si: un SNP consiste en un A/T qui se transform en C/G. peu importe lequel de la paire, c'est comme ça il semblerait.



aujourd'hui, j'ai fait 9h30-18h00 avec 30 minutes de pause le midi.
Donc 8h de travail.

# <span style="color:#999900"> 15/02/2022
Je continue de suivre la manip avec Laetitia. ce matin, on s'est arrêtés après avoir fini la PCR 1.
Pour comprendre le fonctionnement global de la manip, voir https://en.wikipedia.org/wiki/DNA_microarray , la vidéo explique très bien le fonctionnement.

PCR: chercher "youtube PCR" donne une très bonne vidéo.
1. on chauffe fort pour séparer les 2 brins d'ADN
2. à une température plus basse, on laisse les primers s'hybrider sur les séquences d'intérêt
3. on chauffe un peu pour que les Taq polymérases se fixent aux primers. Ces dernières répliquent les brins d'ADN aux régions concernées.
les étapes 2 et 3 sont répétées pour plusieurs cycles, ce qui double le nombre de copies de à chaque étape, résultant en une augmentation exponentielle.

Utiliser la partie précédente dans l'écriture du rapport de stage.

Je commence à bien connaître le protocole. Le but final est de déterminer l'index génomique et je cite `/home/waren/Desktop/stage_M2/sent_by_claire/sujet_stage_projet_GIRONDE__copie_from_administration.pdf` : "bien appréhender les critères utilisés [en utilisant les puces Agilent] pour la détection des variants afin de pouvoir proposer et développer une approche automatisée [...]\[à partir des] données Affymetrix".
je cite également `Projet_GIRONDE_synopsis23052019.pdf`: "L’objectif de cette étude est de transposer le calcul de l’index génomique et la détermination du seuil de classification des tumeurs de la technologie Agilent à la technologie Affymetrix/Oncoscan".

─> comment l'index génomique est-il déterminé pour Agilent?
    GI = A^2/C
    où A = total number of alterations (segmental gains and losses)
    et C = the number of involved chromosomes.
    -> C est déterminé par le protocole. on travaille sur l'humain donc 23 paires de chromosomes.
    -> A est calculé par l'expérience.

─> comment la détermination du seuil de classification des tumeurs est-elle faite pour Agilent? parle-t-on de la valeur de 10 pour le GI?

TODO: lire la prise de notes que j'ai fait avec claire et élodie. noter les infos ici, et les choses à faire.
arrivée à 8h30 ----> départ à 17h. 30 min pause midi.
Donc 8h de travail.



# <span style="color:#999900"> 16/02/2022
J'ai suivi la fin de la manip avec Laetitia. J'ai vu comment les résultats sortaient d'affymétrix (.CEL, .ARR, .DAT), et le logiciel qui est utilisé pour les traiter. Cependant, je vais utiliser autre chose pour traiter ces fichiers, certainement des packages R.

TODO:
- lire la prise de notes que j'ai fait avec claire et élodie. noter les infos ici, et les choses à faire.
- voir les infos que contiennent les .CEL, .DAT et .ARR
- chercher les logiciels/packages qui lisent ces fichiers
- faire un bon récap de la manip et éclaircir les points obscurs:
    - PCR double brin
    - SNP?
- ~~organiser les dossiers comme élodie l'a indiqué~~
- bloquer le Jeudi 3 mars à 11h: https://u-bordeaux-fr.zoom.us/j/82532998606?pwd=a0Q3aWZ3ZjZMdC9udXcxem85clJPUT09
- noter les mardis de 12h à 14h: faire un point avec Élodie et Claire (et Sabrina!)
- noter qu'au début de chaque semaine, je dois aller voir la team CGH pour savoir qui analyse les résultats de la semaine et quand.

J'ai obtenu les codes pour me connecter à un ordi fixe de l'institut. la question actuelle est: ai-je un mail \@bordeaux.unicancer.fr ?

Un repository git est créé et a été cloné sur le PC fixe de bergonié et sur le mien. la version la plus avancé est cassebriques.

Prochaine étape: récupérer des .CEL , .DAT , .ARR, et les explorer avec un logiciel/package R que j'aurai trouvé. sinon relire la todo list.
arrivée à 9h10 ----> départ à 17h10. 30 min pause midi.
Donc 7h30 de travail



# <span style="color:#999900"> 17/02/2022

TODO:
- lire la prise de notes que j'ai fait avec claire et élodie. noter les infos ici, et les choses à faire.
- voir les infos que contiennent les .CEL, .DAT et .ARR
- chercher les logiciels/packages qui lisent ces fichiers --> le package R de l'article.
- faire un bon récap de la manip et éclaircir les points obscurs:
    - ~~PCR double brin~~\
    les 2 brins sont formés séparément come le montre la vidéo, et se lient l'un à l'autre pour former le produit final double brin.
    - SNP?
- ~~organiser les dossiers comme élodie l'a indiqué~~
- bloquer le Jeudi 3 mars à 11h: https://u-bordeaux-fr.zoom.us/j/82532998606?pwd=a0Q3aWZ3ZjZMdC9udXcxem85clJPUT09
- noter les mardis de 12h à 14h: faire un point avec Élodie et Claire (et Sabrina!)
- noter qu'au début de chaque semaine, je dois aller voir la team CGH pour savoir qui analyse les résultats de la semaine et quand.
- ~~lire en profondeur l'article arm-level.~~ CCL & discussion: "notre méthode est aussi efficace que les experts humains et bien plus rapide."
    - savoir expliquer cette méthode et comment l'utiliser.
    - prendre des données auprès de l'équipe technique. apporter une clé USB demain! récupérer les .ARR, .DAT et surtout .CEL .
    - l'utiliser sur nos données.

Ai installé R et Rstudio sur Bergonié, ai téléchargé le package R de l'article arm-level à https://github.com/yannchristinat/oncoscanR-public.
faire un push propre sur cass et (tenter de) le pull sur bergonié.
le push c'est bon mais bergonié ne peut même pas faire de commits. la commande `git add` à elle seule trigger ce message:
```
C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie>git add
warning: unable to access 'P://.gitconfig': Permission denied
warning: unable to access 'P://.gitconfig': Permission denied
warning: unable to access 'P://.gitconfig': Permission denied
fatal: unknown error occurred while reading the configuration files
```
voir la réponse de jennifer pour ça. potentiellement demander un linux

arrivée à 10h -> départ à 16h50 = 6h20 de travail
total cumulé sur la semaine: 29h50.
pour faire 35h, reste 5h10.



# <span style="color:#999900"> 18/02/2022
ai résolu le pb qui m'empêchait de pull en supprimant la variable d'envt HOMEPATH: https://stackoverflow.com/questions/14774159/git-warning-unable-to-access-p-gitconfig-invalid-argument

Ai fait un push clean. je travaille maintenant essentiellement sur bergo.

TODO:
- lire la prise de notes que j'ai fait avec claire et élodie. noter les infos ici, et les choses à faire.
- faire un bon récap de la manip et éclaircir les points obscurs:
    - ~~PCR double brin~~\
    les 2 brins sont formés séparément come le montre la vidéo, et se lient l'un à l'autre pour former le produit final double brin.
    - SNP?
- ~~organiser les dossiers comme élodie l'a indiqué~~
- bloquer le Jeudi 3 mars à 11h: https://u-bordeaux-fr.zoom.us/j/82532998606?pwd=a0Q3aWZ3ZjZMdC9udXcxem85clJPUT09
- noter les mardis de 12h à 14h: faire un point avec Élodie et Claire (et Sabrina!)
- noter qu'au début de chaque semaine, je dois aller voir la team CGH pour savoir qui analyse les résultats de la semaine et quand.
- ~~voir les infos que contiennent les .CEL, .DAT et .ARR~~
- ~~chercher les logiciels/packages qui lisent ces fichiers~~ --> le package R de l'article.
- ~~lire en profondeur l'article arm-level.~~ CCL & discussion: "notre méthode est aussi efficace que les experts humains et bien plus rapide."
    - savoir expliquer cette méthode et comment l'utiliser.
    - ~~prendre des données auprès de l'équipe technique. apporter une clé USB demain! récupérer les .ARR, .DAT et surtout .CEL .~~
    - l'utiliser sur nos données.

.ARR = du xml

Laetitia m'a passé les données anonymisées d'affymetrix. Je fais passer le premier échantillon dans le package R.

path to R.exe: "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\R.exe"
path to Oncoscan-R script: "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\run_oncoscan_workflow.R"
path to first txt file: "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data\1-RV.OSCHP.segments.txt"

to set an environment variable:
> setx r_exe "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\R.exe"

to use it:
> %r_exe%

output:
```
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R est un logiciel libre livré sans AUCUNE GARANTIE.
Vous pouvez le redistribuer sous certaines conditions.
Tapez 'license()' ou 'licence()' pour plus de détails.

R est un projet collaboratif avec de nombreux contributeurs.
Tapez 'contributors()' pour plus d'information et
'citation()' pour la façon de le citer dans les publications.

Tapez 'demo()' pour des démonstrations, 'help()' pour l'aide
en ligne ou 'help.start()' pour obtenir l'aide au format HTML.
Tapez 'q()' pour quitter R.
>
```

commands entered to create envt variables:
> setx r_exe "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\Rscript.exe"\
> setx oncos-r "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\oncoscan-workflow.R"\
> setx data_folder "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data"\

command to run workflow: 
> r_exe 


Unrelated but the theme Shades of Purple is cool but for markdown editing I prefer built-in theme Monokai Dimmed.
Also on Bergo I save my keybindings.json file in C:\Users\e.bordron\Documents .

I created 2 folders for data: raw_data, a backup folder, and working_data, which I will be working on. Before doing that, I unintentionaly lost the file segments.txt for the first sample by sending the output of a command into it.

I did this:

```
> "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\Rscript.exe" "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\run_oncoscan_workflow.R"  "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data\1-RV.OSCHP.segments.txt" M

Erreur dans UseMethod("collector_value") :
  pas de méthode pour 'collector_value' applicable pour un objet de classe "c('collector_skip', 'collector')"
Appels : workflow_oncoscan.run ... load_chas -> read_tsv -> <Anonymous> -> collector_value
De plus : Message d'avis :
The following named parsers don't match the column names: CN State, Type, Full Location
Exécution arrêtée
```

Prochaine chose à faire: modifier les colonnes du .txt pour qu'il y ait seulement le s3 dont le package a besoin.

nouveau problème quand je lance cette ligne de commande:
> "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\Rscript.exe" "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\run_oncoscan_workflow.R"  "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data\woring_data\2-AD\2-ADREC.RC.OSCHP.segments.txt" F

> Accès refusé.

et une pop-up apparaît: 
```
Cette application ne peut pas s'exécuter sur votre PC.
```
même quand je vais dans C:\Users\e.bordron\Documents\R\R-4.1.2\bin et que je fais:
> Rscript.exe

le même problème survient.
essayer de redémarrer.

Je viens de redémarrer. je me rends compte sur `https://helpdeskgeek.com/windows-10/how-to-fix-this-app-cant-run-on-your-pc-in-windows-10/` que Rscript.exe est en 32 bit.
edit: j'ai les 2 Rscript: le 32bit et le 64bit.
le chemin du 32bit: `C:\Users\e.bordron\Documents\R\R-4.1.2\bin`
le chemin du 64bit: `C:\Users\e.bordron\Documents\R\R-4.1.2\bin\x64`

J'ai le problème suivant, en tapant la bonne ligne de commande:
> "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\x64\Rscript.exe" "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\run_oncoscan_workflow.R"  "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data\working_data\2-AD\2-ADREC.RC.OSCHP.segments.txt" F
```
Erreur dans load_chas(chas.fn, oncoscan.cov) : Parsing ChAS file failed.
Appels : workflow_oncoscan.run -> load_chas
De plus : Message d'avis :
The following named parsers don't match the column names: Full Location
Exécution arrêtée
```

on en revient à cette solution: Prochaine chose à faire: modifier les colonnes du .txt pour qu'il y ait seulement les 3 dont le package a besoin:  `Type, CN State and Full Location`
Je viens d'essayer: les colonnes Type et CN state sont bien présentes, mais pas la colonne Full Location. J'aurais du lui dire plus tôt, mais je vais demander à Laetita si il est possible d'avoir cette 3eme colonne. Si ce n'est pas possible, peut-être qu'il est possible de récupérer cette information à partir d'autres colonnes, auquel cas je regarderai d'autres moyens d'obtenir des données de ces .CEL. parser moi-même ce fichier est aussi une option


arrivée à 10h10 -> 17:30 = 6h50 de travail et 30 min de pause le midi.
1h40 est déjà faite pour la semaine prochaine.



# <span style="color:#999900"> 21/02/2022

Je regarde vite fait si l'une des colonnes contient par hasard l'information "Full location", sinon je demande à Laetitia.
La colonne `Full location` semble être un arrondi au 100 des valeurs de position contenues dans `Microarray Nomenclature`. Je fais une colonne Full location à partir de ça. en R:
- import CSV as dataframe
- get value between parenthesis: `(754,192-145,095,477)` from column `Microarray Nomenclature`
- get chromosome from column `Chromosome`
- create column `Full location` that contains such values: `chr7:129199300-129813700`

fait. je lance oncoscan-R dessus:
> "C:\Users\e.bordron\Documents\R\R-4.1.2\bin\x64\Rscript.exe" "C:\Users\e.bordron\Documents\R\R-4.1.2\library\oncoscanR\bin\run_oncoscan_workflow.R"  "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\data\working_data\2-AD\2-ADREC.RC.OSCHP.segments_FULL_LOCATION.txt" F

Cela me donne une erreur: 
```
Erreur dans load_chas(chas.fn, oncoscan.cov) : Parsing ChAS file failed.
Appels : workflow_oncoscan.run -> load_chas
De plus : Message d'avis :
The following named parsers don't match the column names: Full Location
Exécution arrêtée
```

je le lance aussi dans R. voir csv_formatting.R. le message d'erreur est:
```
Error in load_chas(chas.fn, oncoscan.cov) : Parsing ChAS file failed.                                                                                       
In addition: Warning message:
The following named parsers don't match the column names: Full Location 
```
Or, la colonne Full location est bien écrite.

résolu. maintenant:
> workflow_oncoscan.run("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD/2-ADREC.RC.OSCHP.segments_FULL_LOCA ..." ... [TRUNCATED] 
```
Error in if (length(parm) == 0 || seg_start > end(parm)) {:

    missing value where TRUE/FALSE needed
In addition: Warning messages:
1: In load_chas(chas.fn, oncoscan.cov) : NAs introduced by coercion
2: In load_chas(chas.fn, oncoscan.cov) : NAs introduced by coercion
```

arrivée à 10h25, départ à 19:25
