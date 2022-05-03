# Intro
Contexte: 
Dans le traitement du cancer, la prise de décision thérapeutique passe par des analyses de l'ADN.  `pourquoi? aussi par la protéomique, ...`
Les cancers sont souvent caractérisés, pour certaines régions chromosomiques, par un nombre de copies altéré.  `quels cancers? quelles régions? pourquoi? -> parce que les gènes du pathway HR empêchent le dvp des tumeurs, et leur perte est nécessaire à certains cancers.`
La technologie array-CGH permet d'étudier le nombre de copies de l'ADN d'un individu. `Quelles autres technologies permettent la`
Un score, appelé Index Génomique, permet d'évaluer le degré d'altération d'un profil.  `parler d'autres scores: Signature Nanocind...`
Ce score est utilisé en routine à l'institut Bergonié avec la technologie Agilent pour classifier les tumeurs (bénigne / maligne).  
La technologie offrant une meilleure résolution, OncoScan CNV, est utilisée depuis peu, mais le calcul du GI n'a pas été caractérisé pour cet outil plus précis.   ``"Dans la CGH, on retrouve une technologie récente, OncoScan, ayant une meilleure résolution qu' Agilent."``
Le but de ce stage est de transposer le calcul du GI à cette nouvelle technologie, et de déterminer le nouveau seuil permettant de classifier les tumeurs.   
D'autre part, cette procédure (pas la technique, le calcul du score) nécessite aujourd'hui des manipulations humaines.
L'objectif est également d'automatiser le processus de détermination du GI dans un pipeline d'analyse de routine.  

Dit plus concrètement, l'objectif du stage est de mettre en place un outil qui lira les données OncoScan CNV et calculera l'index génomique automatiquement, et d'utiliser ses résultats pour définir un nouveau seuil



# Etat de l'art
définir le HRD, donc le homologous recombination pathway  
définir la CGH, donc le principe de comparer le signal avec une référence pour avoir une valeur de log raito
définir Agilent 

Comment le GI est déterminé ?

Parler des outils; segmentation -> DNAcopy, calling...


## parler du GI: 
La signature Nanocind est un meilleur indicateur que le GI, mais non applicable dans notre cas.



# Méthodes
```
En arrivant à cette partie, on doit savoir:
- pourquoi un outil qui apporte le calcul d'autres scores HRD est intéressant
    - qu'est-ce que les scores HRD
        - qu'est-ce que le HR pathway et quel est son lien avec le cancer
- quels sont ces outils
    - que font-ils
    - pourquoi on les compare
    - comment ils ont été choisis
```
Quatre outils permettant d'intégrer les données OncoScan sont comparés.  

Comparaisons: comment elles ont été menées. Prise en compte des scores HRD en plus; importance du nombre d'étapes manuelles que chaque outil ajoute, des remaniements des données qui sont à éviter... 
Aussi, parler du fait que le GI va garder la même proportionnalité d'un échantillon à l'autre, a priori, mais les valeurs en elles-même vont changer (une plus grande résolution impliquant un plus grand nombre d'altérations), d'où le besoin de redéfinir le seuil.


# Resultats
courbes ROC, tests de performance... Quel outil est le plus adapté.  
Comparer ce qui est comparable: combien de temps pour obtenir le GI? -> temps par échantillon. Un échantillon plus lourd fait-il *1.5 sur le temps de calcul d'un outil, et *10 pour un autre? enlève-t-on les préprocess, filtrages et smoothings (oui.)?
La proportionnalité avec le GI obtenu manuellement est-elle conservée? quels échantillons

# Discussion
Les + et les - de tous les outils. tableau?
quel outil est le plus performant -> précision; vitesse.

OncoScan apporte, par rapport à Agilent, plus qu'une meilleure résolution: La CGH indique la notion de statut allélique, ce qui peut être un paramètre à considérer pour un nouveau score HRD. Prendre en compte la perte d'allèles de gènes précis liés à la HRD serait intéressant pour caractériser au mieux les échantillons tumoraux.
