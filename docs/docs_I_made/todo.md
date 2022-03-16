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
- ~~refaire un ticket? ... je ne sais pas a priori non, sur l'ent le ticket a toujours l'air en cours de gestion. garder un oeil là-dessus.~~ J'ai renvoyé un mail sans avoir de réponse aujourd'hui.
- ~~si pas de réponse demain midi au sujet du ticket, en refaire un.~~
    je refais un ticket.
- Lire `LRR_and_BAF.pdf` pour mieux comprendre ces 2 indices.
- ~~chercher si on peut exporter des cnchp à partir d'oncoscan.~~
    -> non, selon le manuel de ChAS.
- Chercher aussi la différence entre oncoscan et Affymetrix SNP 6.0 comme elodie l'indique
- voir si les infos requises par rCGH, qui proviennent d'un CNCHP, peuvent être trouvées dans mes OSCHP à l'aide de HDFView .
- à https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/index.html, regarder les différents formats de fichier.
- ~~voir si je peux accéer à la base de données ennov. en ligne~~
    -> oui
- lire l'article qui décrit ASCAT: ASCAT.pdf
10/03/2022
- lire les articles des packages. je commence par CGHcall. noter dans fonctionnement_des_outils.md et notes_on_articles.md ce que je trouve.
- répondre à la question: ça fait quoi quand on lance CGHcall sur plusieurs échantillons en même temps? est-ce que ça change les résultats par rapport à lancer les analyses échantillon par échantillon? Si oui, quels calculs sont faits?
pour répondre à ça, comprendre le fonctionnement du le sa segmentation couplée au mixture model serait un bon point de départ. je lis donc Picard, voir notes_on_articles.md.
Réponse à ça: lancer CGHcall sur un groupe d'E ou sur les E un à un change bien les résultats. En effet, les sondes sont classées par segment (ou par bras si le mm alternatif est utilisé), et pour cela, la distribution de *toutes* les sondes est utilisée.

