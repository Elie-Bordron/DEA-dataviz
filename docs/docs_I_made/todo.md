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

28/02/2022
- faire un plan pour intro 

01/03/2022
- ~~changer le script csv_formatting.R pour qu'il enlève les virgules de la colonne Full Location~~
- les différents indices pouvant êre calculés par oncoscan-R sont-ils intéressants? de quoi témoignent-ils?
    -> voir surtout ceux qui sont cités dans un article.
- dans notes_on_articles.md, voir la dernière ligne de l'article Nanocind. lire l'article en question.
- dans cet article Nanocind (qui est la thèse de Sabrina Croce), on conclut que la signature Nanocind est un meilleur indice de prédiction que le GI. voir si on ne pourrait pas l'appliquer à nos données? Sabrina m'en avait parlé mais je ne vois pas de mention de ça dans mes notes. 
    edit: page 4, table 1, je vois que CINSARC et Nanocind sont peut-être la même chose.  
    
02/03/2022
- ~~refaire un ticket? ... je ne sais pas a priori non, sur l'ent le ticket a toujours l'air en cours de gestion. garder un oeil là-dessus. J'ai renvoyé un mail sans avoir de réponse aujourd'hui.~~
- ~~si pas de réponse demain midi au sujet du ticket, en refaire un.~~
    ~~je refais un ticket.~~
- Lire `LRR_and_BAF.pdf` pour mieux comprendre ces 2 indices.
- ~~chercher si on peut exporter des cnchp à partir d'oncoscan.~~
    -> non, selon le manuel de ChAS.
- Chercher aussi la différence entre oncoscan et Affymetrix SNP 6.0 comme elodie l'indique
- voir si les infos requises par rCGH, qui proviennent d'un CNCHP, peuvent être trouvées dans mes OSCHP à l'aide de HDFView .
- à https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/index.html, regarder les différents formats de fichier.
- ~~voir si je peux accéer à la base de données ennov. en ligne~~
    -> oui
- ~~lire l'article qui décrit ASCAT: ASCAT.pdf~~  

10/03/2022
- lire les articles des packages. je commence par CGHcall. noter dans fonctionnement_des_outils.md et notes_on_articles.md ce que je trouve.
- ~~répondre à la question: ça fait quoi quand on lance CGHcall sur plusieurs échantillons en même temps? est-ce que ça change les résultats par rapport à lancer les analyses échantillon par échantillon? Si oui, quels calculs sont faits?~~
pour répondre à ça, comprendre le fonctionnement du le sa segmentation couplée au mixture model serait un bon point de départ. je lis donc Picard, voir notes_on_articles.md.
Réponse à ça: lancer CGHcall sur un groupe d'E ou sur les E un à un change bien les résultats. En effet, les sondes sont classées par segment (ou par bras si le mm alternatif est utilisé), et pour cela, la distribution de *toutes* les sondes est utilisée.  

23/03/2022
- ~~finir de noter pour ASCAT, il reste 2 phrases.~~

24/03/2022
- ~~lire ASCAT article pour comprendre EaCoN~~
* [ ] lire DNAcopy article pour comprendre ~~segmentation CGHcall~~ et ~~format CBS~~ (utilisé dans plusieurs outils)
le format CBS (output DNAcopy) est décrit dans fonctionnement_des_outils.md. 
* [ ] écrire ce que j'ai compris de la segmentation par DNAcopy. C'est dans le cahier.
* [ ] décrire le pipeline pour chaque outil dans fonctionnement_des_outils.md.
* [ ] connaître précisément les spécificités de chaque outil: qu'apporte-t-il par rapport aux autres? quelles sont ses limites par rapport aux autres?
* [ ] mettre à jour le tableau comparatif
* [ ] indiquer à quoi sert quel score HRD donné par OncoscanR (cf cahier)
* [ ] lire Mat Met ascat pour savoir comment ASCAT procède pour avoir les LOH et les CN-neutral events alors que la CGH ordinaire ne permet pas d'avoir ça. cela entre dans le cadre de "connaître les spécificités de chaque outil".
* [ ] me remettre en mémoire rCGH, faire le pipeline.
* [ ] EaCoN, lors de l'estimation du nombre de copies, utilise des paramètres "gamma", et recommande de lire les pages d'aide du package R ASCAT pour savoir de quoi il s'agit. regarder ça.
* [ ] transférer les infos sur CGHcall du cahier vers fonctionnement_des_outils.md.
* [ ] comment ASCAT estime la cellularité?
~~comment ASCAT estime le copy number (ASCN)?~~
En se basant sur les valeurs de log ratio et de BAF données par les puces Oncoscan.
* [ ] Un algorithme dit "ASPCF" traite les données BAF et log ratio. voir ce que c'est et où il est utilisé.
* [X] ~~*ASCAT Mat et Met from cahier to md*~~ [2022-04-04]

J'ai le choix:
* [ ] pipeline pour les 4 packages (manque oncoscanR et rCGH) + comparaison tableau
* [ ] slides pour 2 packages (CGHcall et EaCoN) mais EaCoN produit surtout des graphes -> parler de ces derniers, et de ce qu'ASCAT promet.
* [ ] présenter tous les packages sans les faire tourner? -> que faire pour ça?
* [ ] faire un tableau-pipeline pour montrer ce qu'ils font de *similaire*, ça permet de montrer aussi ce qu'ils font de *différent*. pour cela, les 4 pipelines sont nécessaires.  
***je fais ça.*** d'abord, les 4 pipelines.

28/02/2022
- ~~envoyer un rappel pour la réunion de demain.~~  
29/02/2022  
* [X] ~~*remplir le CR de réunion*~~ [2022-03-29]  
* [ ] se renseigner sur Sequenza  
* [X] ~~*voir définition BAF et allele difference avec laetitia. voir cahier, 29/03/22*~~ [2022-03-30]
* [X] ~~*partager les articles de CGHcall et EaCoN (donc ASCAT).*~~ [2022-04-04]
* [X] ~~*faire le tri des choses à faire qui sont notées dans le cahier à la date d'aujourd'hui, et les noter ici.*~~ [2022-03-30]
* [ ] définir et pouvoir expliquer les Copy Number-Neutral Events (sera certainement une diapo pour EaCoN)
* [ ] se renseigner sur normalisation d'EaCoN par APT
* [ ] Centralisation et calling d'EaCoN dans la segmentation: qu'est-ce?
* [ ] normaliser les organigrammes pour qu'ils soient faciles à comparer
* [ ] ajouter l'output de chaque étape dans les organigrammes
* [ ] quand claire m'aura partagé le dossier GIRONDE , j'y ajouterai les infos sur chaque package. un dossier pour chaque.
* [ ] pour la formation ChAS: demander à Tony si ChAS utilise un modèle de mixture gaussien
* [ ] puis-je voir la distribution gaussienne de toutes les données? oui, voir logratios_for_hist dans CGHcall.R . ça plot une unique courbe de gauss; peut-être qu'avec de vraies données ce serait plus parlant. voir sur slack, Elodie m'avait fait un retour.
31/03/2022
* [ ] lire l'article sur la normalisation de rCGH
01/04/2022
* [ ] se renseigner sur sequenza.
* [ ] à demander à Tony Sierra: a-t-on un moyen d'exporter en batch des fichiers segments.txt ?
* [ ] comparer la validation des techniques: les auteurs ont validé les résultats d'OncoscanR en les comparant avec l'annotation d'un expert sur 25 cas. quid de CGHcall, EaCoN et rCGH?
* [ ] parmi les articles que m'a envoyé claire, en chercher un qui contient des infos sur le système AT/GC des puces oncoscan. voir cahier
* [X] ~~*noter les infos sur l'export de segments.txt dans logbook*~~ [2022-04-01]
* [ ] noter ce que j'ai surligné dans l'article de rCGH
04/04/2022
* [ ] ASPCF: cahier to md


