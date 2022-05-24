# <span style="color:#0099ff"> Introduction
- § Contexte sur la cancérologie: 
    - L'institut Bergonié
    - Le diagnostic nécessite d'identifier la tumeur
    - comment classer les tumeurs:
        - critères morphologiques/anatomiques
        - signatures et scores moléculaires
            - GI
- § La technologie array-CGH: 
    - permet de calculer des scores moléculaires
    - Son fonctionnement
    - Existe en parallèle du Whole Genome / Whole Exome
    
- § à l'Institut Bergonié, ces technologies sont utilisées:
    - Agilent (utilisée pour calculer le GI)
    - CytoScan (simplement le citer)
    - OncoScan
    - OncoScan > Agilent (résolution)
    - on veut donc étendre le calcul du GI à oncoScan


- § Transfert de connaissances
    - échange Biologie-Bioinformatique
    - J'ai suivi la manip
    - Les biologistes ont un aperçu plus en profondeur des outils
    - vulgarisation -> slides


- § En quoi OncoScan est mieux qu'Agilent: cf projet GIRONDE.
    - donc comparaison nécessaire sur les mêmes échantillons (parler surtout des + d'oncoscan par rapport à Agilent, détailer le mode opératoire+ tard)
    - Préciser qu'on va comparer des outils: identifier un outil qui permet de calculer le GI et de mettre en place cette comparaison
    - Aussi: question de l'analyse en routine même si il y aura toujours un regard
    - Question: Peut-on établir la correspondance du calcul de l’index génomique et la détermination du seuil de classification des tumeurs de la technologie Agilent à la technologie Affymetrix pour pouvoir en faire bénéficier les patients analysés au sein du laboratoire dans le cadre diagnostic ou thérapeutique?



<div style="page-break-after: always"></div>

# <span style="color:#0099ff"> Etat de l'art
- § Calcul du GI par Agilent
    CF les liens du projet GIRONDE

- § ChAS:
    - En quoi est-ce utile
    - En quoi est-ce limité: 
        - boîte noire
        - pas de possibilité de calculer ce score ni d'autres
        - On n'a jamais validé les résultats de ce logiciel à l'aide d'une autre technologie donc on ne connaît pas ses biais

- § OncoScanTM Console 1.3 (Logiciel propriétaire d'affymetrix):
    - En quoi est-ce utile
    - En quoi est-ce limité

- § On n'utilise pas ces deux outils, on en cherche un autre qui sera plus adapté
    - ce que l'outil devra faire:
        - calculer le GI
            - donc prendre en input les données OncoScan
        - être ~~utilisable en routine~~ automatisable  même si il y aura toujours un regard
        - être open source, c'est plus reproductible

- § oncoscanR:
    pipeline

- § rCGH:
    pipeline
    
- § CGHcall:
    pipeline 

- § ASCAT:
    pipeline

<div style="page-break-after: always"></div>

# <span style="color:#0099ff"> Matériel et Méthodes
## Matériel
- Données utilisées
    - préciser que je travaille sur des données brutes/traitées, quels fichiers.
    - ADN tumoral, FFPE, tumeurs GIST (biopsie), patients... quelques phrases pour tout ça

- Logiciels utilisés
    ` tableau comparatif par étape: l'outil fait cete étape ou non.`
    - ChAS version X
    - Rstudio version X
    - R version X
    - packages R: 
        - oncoscanR version X -> article...
        - CGHcall version X -> article...
        - ASCAT version X
        - rCGH version X
    
## Méthodes
- Nativement, chaque outil a des + et des -. -> détailler exaustivement, sélectionner les figures les plus pertinentes de chaque ppt, puis-je en fusionner 2 ensembles, ...

- Comment compare-t-on les outils?
    - § GI
        - comment il est calculé: à partir des segments d'altération
        - ce qui a dû être adapté pour que les 4 outils fassent ce calcul
        - Y a-t-il une proportionnalité similaire entre Agilent et OncoScan? -> graphe de corrélation (*ici, on annonce que cette question est posée dans la partie résultats*)
        - GGally
    - § performance
        - vitesse
        - précision
            - courbes ROC



<div style="page-break-after: always"></div>

# <span style="color:#0099ff"> Resultats & Discussion

- Sachant les plus et moins de chaque outil, on fait une comparaison des outils:
    - GI
        - Y a-t-il une proportionnalité similaire entre Agilent et OncoScan? -> graphe de corrélation (*ici, on répond à cette question*)
        - Y a-t-il des outils qui n'en expriment pas?
    - performance
        - vitesse: quel est l'outil le plus rapide
        - précision: quel est l'outil le plus précis? (On prend aussi en compte la spécificité, et d'autres paramètres... )
            - courbes ROC
        - Sur cette comparaison, quel est le meilleur compromis vitesse/précision?
- Au vu de ce compromis, du GI obtenu et des points + et - de chaque outil vus dans l'état de l'art, quel est l'outil retenu?

- Conclusion sur la question
