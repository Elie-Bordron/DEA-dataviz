
# <span style="color:#ff9999">Automated Detection of Arm-Level Alterations for Individual Cancer Patients in the Clinical Setting
2021_Christinat_HRD_Oncoscan.pdf

Copy number alterations, a genetic event, promotes tumor development. these events are used as predictive biomarkers in clinical care. they are roughly classified as arm-level or focal. genome-wide techniques exist to classify arm-level ones, but challenges exist:
- How to define an arm-level alteration? there is no consensus on it.
- there is a lack of tools to compute them for individual patients.
To answer this, using OncoScan, clinical samples were analyzed. The results indicate respectively that:
- sum of altered segments was a better indicator than longest segment to define an arm-level alteration, BUT Some of the discordances ultimately were attributed to human error.
- a new software has been made publicly available in routine analyses (https://doi.org/10.1016/j.jmoldx.2021.08.003)

Les altérations arm-level sont plus longues que les focales: elles contiennent des centaines de gènes. Leur effet sur les phénotypes de cancer est connu, mais il est compliqué, quand ce n'est pas impossible, de lier cet effet à des gènes individuels.
Cependant, les avancées technologiques récentes en génomique offrent des indications sur la valeur prédictive des arm-level alterations: on associe par exemple la perte de tels bras de tels chromosomes à un critère principal de diagnostic pour les oligodendrogliomes. Il y a plusieurs autres exemples de ça, également pour le pronostic. D'autre part,  l'aneuploidie (le fait pour une cellule de ne pas être euploide, donc d'avoir un nombre anormal de chromosomes) a un intérêt dans le pronostic mais n'est pas toujours utilisé dans la prise de décision thérapeutique.
Certaines méthodes utilisent du FFPE mais ces manips sont sujetes à dégrader l'ADN. la technique Oncoscan est conçue et optimisée pour le FFPE. D'autre part, GISTIC est une méthode largement utilisée pour estimer les Copy Number Alterations (CNA).
Enfin, la définition d'arm-level varie dans les exemples présents dans la littérature. Le seuil d'altération d'un bras chromosomique à partir duquel on ne considère plus l'altération comme focale mais comme arm-level varie au cas par cas.
Cette étude vise à définir précisément les arm-level alterations qui peuvent être appliquées dans un contexte clinique.

En utilisant la technique Oncoscan sur des échantillons FFPE, des .CEL ont été produits, puis convertis en .OSCHP (OncoScan array data) qui ont été analysés à l'aide du logiciel Chromosome Analysis Suite (ChAS) avec le génome de référence: hg19. Les résultats ont été évalués manuellement. Précision: les segments de moins de 50 marqueurs ou moins de 50 Kbp ont été ignorés.



## Estimation of Percentage Arm Alteration
les données de nombre de copies ont été utilisées pour  détecter des arm-level alterations (ALA), puis ont été processées par un script R


# <span style="color:#ff9999">Uterine smooth muscle tumor analysis by comparative genomic hybridization: a useful diagnostic tool in challenging lesions
modpathol20153.pdf

the diagnosis of STUMP tumors is often challenging. The authors proposed to test the hypothesis that GI could split those tumors in two groups: benign and malignant. A study was conducted, using a threshold of 10 to classify the STUMP.
In conclusion, array-CGH is a useful technique to classify these tumors.


# <span style="color:#ff9999">The Nanocind Signature Is an Independent Prognosticator of Recurrence and Death in Uterine Leiomyosarcomas

Nanocind_signature_S._CROCE.pdf

Uterine leiomyosarcoma is an aggressive tumor responsible for a signiﬁcant proportion of uterine cancer–related deaths. Plus, using the FIGO staging system, it is currently impossible to predict the clinical outcome of stage I leiomyosarcomas. However, the authors published in 2010 a transcriptomic signature (67 genes related to chromosome biogenesis, mitosis control, and chromosome segregation), which has proven since its predicting efficiency over different cancer types. Plus, it has been successfully used with NanoCind (Nanostring) technology, which makes it usable routinely.
Uterine leiomyosarcoma were analyzed with the Nanocind signature. The process split the group in two groups. This result was validated.
In conclusion, the NanoCind signature is a powerful prognostic indicator that outperforms FIGO staging and the genomic index. Plus, GI is platform-dependent. 




