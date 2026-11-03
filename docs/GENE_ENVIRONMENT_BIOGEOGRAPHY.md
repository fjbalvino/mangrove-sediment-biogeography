# Biogeografía de genes funcionales en sedimentos de manglar

Protocolo analítico y revisión de herramientas para los 191 metagenomas funcionales

Fecha de revisión 29 de octubre de 2026

## Decisión y alcance

El objetivo es explicar cómo cambia el repertorio genético de las transformaciones biogeoquímicas entre manglares y con qué gradientes ambientales se asocia. Proponemos integrar tres resultados: mapas de ciclos con valores por gen, modelos de asociación gen ambiente y una comparación de la estructura geográfica de funciones y taxonomía. El esquema bioquímico será el contexto de interpretación de resultados cuantitativos, no una red de correlaciones interpretada como intercambio metabólico.

Trabajaremos con los 191 metagenomas aceptados por Fran. La matriz permite comenzar una biogeografía de familias funcionales representadas por KOs; no conserva por sí sola la identidad de cada variante de secuencia. El salto posterior a genes individuales, variantes y portadores requerirá recuperar el vínculo entre secuencias, anotaciones, contigs y, cuando existan, MAGs. Esa ampliación no impide desarrollar ahora el análisis comunitario.

El análisis tendrá dos niveles complementarios. El nivel amplio considerará los 17.482 KOs detectados al menos una vez, con filtros de estimabilidad fijados antes de buscar efectos. El nivel dirigido utilizará el panel de 418 KOs consultados, las 97 reglas y las 65 categorías de procesos. Un KO puede participar en varios procesos; no duplicaremos sus observaciones ni contaremos esas pertenencias como replicaciones.

En esta entrega se ejecutó la auditoría de entradas, disponibilidad ambiental y rango del diseño. La revisión de herramientas incluye documentación y código seleccionado. No se instalaron ni se compararon experimentalmente todos los programas, y todavía no se estimaron asociaciones gen ambiente. Las hipótesis y ejemplos de historias ecológicas que siguen son propuestas, no resultados del conjunto de datos.

## Antecedentes que se pueden adaptar

Un KO es una categoría de ortólogos de KEGG: agrupa genes con una función de referencia, aunque algunas familias necesitan evidencia adicional para distinguir sustratos o direcciones. En este protocolo, abundancia significa representación relativa de asignaciones a esa categoría. El ADN informa de capacidad genética; medir cuánto aumenta una enzima requeriría cuantificar su proteína, y medir su actividad requeriría otro tipo de ensayo. Una asociación ambiental indica que dos variables cambian conjuntamente bajo el modelo, sin establecer por sí sola causalidad.

Ocean Gene Atlas 2.0 es el precedente visual más cercano: integra distribución geográfica de genes, abundancias y contexto ambiental, con filtros interactivos. Su aplicación pública está construida sobre catálogos oceánicos; no admite nuestra matriz de manglar como sustitución automática. Utilizaremos su organización conceptual, no sus resultados ni sus valores de referencia. Su versión actual documenta código en GitLab; la versión inicial también publicó código en GitHub [1].

El estudio Global biogeography of microbes driving ocean ecological status under climate change integra perfiles taxonómicos, marcadores biogeoquímicos, ambiente y predicción, con código en ZhenyanZhang/Ecological_status [2]. Es útil como precedente de integración. Su cuantificación mediante lecturas y equivalentes genómicos difiere de nuestros conteos históricos, de modo que no trasladaremos directamente sus normalizaciones ni sus proyecciones climáticas.

Faure y colaboradores relacionaron familias de proteínas con gradientes ambientales del océano y publicaron el código en EmileFaure/MAGsProteinFunctionalClusters [3]. Su trabajo respalda combinar una búsqueda funcional amplia con interpretación de rutas. Nuestra validación reservará estudios o conjuntos de sitios, y mantendrá la selección de variables dentro del entrenamiento, en lugar de asumir que cualquier estrategia publicada de selección y validación es transferible.

La búsqueda se dirigió a herramientas con entradas WGS o matrices funcionales, reconstrucción de ciclos, modelos multivariables y validación espacial. Incluyó búsquedas públicas en GitHub y literatura primaria, lectura de README, documentación y funciones seleccionadas. No constituye un inventario exhaustivo de todos los repositorios de GitHub. Se priorizaron compatibilidad, trazabilidad, entrada requerida y límites de interpretación, no número de estrellas.

## Herramientas revisadas y decisión de uso

| Herramienta o repositorio | Aporte y entrada | Decisión para este proyecto |
| --- | --- | --- |
| CNPS.cycle — yuezhengfu/CNPS.cycle | Resúmenes y diagramas de C N P S desde anotaciones y perfiles; portadores mediante tablas de genes y taxonomía [4] | Referencia visual y contraste de curación. Adaptar fórmulas y manejo de datos no disponibles. |
| DiTing — xuechunxu/DiTing y SilentGene/DiTing | Lecturas y ensamblajes para perfiles y esquemas de ciclos; el README original remite a una continuación 2.0.2 con Snakemake [5] | Evaluar cuando se recupere cuantificación por lecturas. La visualización puede adaptarse sin afirmar que se ejecutó el flujo completo. |
| METABOLIC — AnantharamanLab/METABOLIC | Anotación de genomas y perfiles biogeoquímicos; la modalidad comunitaria utiliza lecturas [6] | Segunda fase para secuencias y portadores. No convertir nuestra matriz en una ejecución ficticia de METABOLIC. |
| bigecyhmm — ArnaudBelcour/bigecyhmm | Búsqueda HMM en proteínas, reglas y gráficos de ciclos personalizables [7] | Candidato ligero para verificar familias cuando tengamos secuencias. Su abundancia ponderada por organismos no equivale a abundancia de genes. |
| Tabigecy — ArnaudBelcour/tabigecy | Predicción de funciones a partir de afiliaciones taxonómicas [7] | No usar como evidencia principal de genes observados por WGS. |
| DRAM — BortonWrightonLabs/DRAM | Anotación y resúmenes metabólicos de ensamblajes o genes [8] | Complemento futuro. La rama revisada anuncia DRAM2 en desarrollo; fijar versión y validar antes de incorporarla. |
| ggkegg — noriakis/ggkegg | Gráficos de KEGG con tidygraph ggraph y ggplot2 [9] | Útil para correspondencias y figuras estáticas. La geometría final de los ciclos será curada. |
| Escher — opencobra/escher | Editor y visor de mapas metabólicos con datos sobre genes y reacciones [10] | Candidato visual interactivo; no usar una superposición de valores como estimación de flujo. |
| MaAsLin 3 — biobakery/maaslin3 | Asociaciones multivariables de representación y detección, con covariables y repetición [11] | Motor principal candidato para modelos por KO, con parámetros explícitos y comprobación por gen. |
| ANCOMBC — FrederickHuangLin/ANCOMBC | Análisis composicional con ajustes y repetición [12] | Sensibilidad condicionada a supuestos y tipo de cuantificación. No principal para asumir abundancia absoluta a partir de proteínas recuperadas. |
| vegan y adespatial | Ordenación, modelos multivariables, partición y estructura espacial [13] | Núcleo del análisis de perfiles completos y comparación entre ambiente y espacio. |
| mgcv | Modelos aditivos y mixtos para respuestas no lineales [14] | Curvas de unos pocos candidatos con soporte suficiente; no ajuste indiscriminado de miles de curvas. |
| blockCV — rvalavi/blockCV | Particiones geográficas, ambientales y por grupos para validación [15] | Usar en predicción y evaluación de transferencia, respetando además el estudio de origen. |
| microeco — ChiLiubio/microeco | Organización y exploración ecológica de perfiles y ambiente [16] | Capa opcional; no delegar la elección del modelo ni usar inferencias taxonómicas como sustituto de genes medidos. |
| FeGenie y PCycDB | Anotación especializada de hierro y fósforo a partir de secuencias [17,18] | Ampliación posterior de familias difíciles de resolver con KOs generales. |
| KEGGDecoder — bjtully/BioData/KEGGDecoder | Completitud de rutas desde anotaciones de genomas [19] | Contraste de definiciones; completitud no es abundancia y un perfil comunitario no es un genoma. |

La combinación recomendada es curación propia versionada, gráficos con ggraph o SVG, mapas con sf y Leaflet, modelos con MaAsLin 3 y vegan, y validación agrupada con blockCV. Escher es una alternativa de representación, no una dependencia imprescindible. El método no requiere instalar simultáneamente todos los paquetes de la tabla.

La inspección de R/Ncycleabundance.R de CNPS.cycle detectó que rellena KOs no encontrados con cero y resume algunos conjuntos mediante promedios o sumas. En nuestro conjunto, eso transformaría KOs de anammox no disponibles en ausencia aparente. También conserva asignaciones históricas Amo/Pmo. Reutilizar sus gráficos exige conservar nuestros estados de información y revisar la identidad funcional. La función revisada y su hash están registrados en config/gene_environment/repository_review.json [20].

## Lo que realmente contiene el conjunto actual

La auditoría conserva run_accession como identificador y no modifica la matriz ni los metadatos. Se comprobó la selección de muestras y se calcularon los siguientes valores desde las entradas actuales.

| Elemento | Resultado auditado | Consecuencia analítica |
| --- | --- | --- |
| Perfiles funcionales | 191 | Base fija de trabajo. No se reabre la curación del hábitat. |
| KOs de la matriz | 17.508 filas; 17.482 detectados y 26 con cero en los 191 | La exploración amplia puede exceder el panel de ciclos. |
| Estudios | 44 | Los perfiles no son 191 observaciones intercambiables. |
| Coordenadas | 184 perfiles; 57 pares de coordenadas redondeados | Son agrupaciones cartográficas provisionales, no 57 unidades biológicas verificadas. |
| BIO1 | 178 perfiles; 34 valores distintos; variación dentro de 7 estudios | El efecto climático ajustado por estudio depende de una parte limitada del conjunto. |
| Carbono orgánico cartográfico | 172 perfiles; 44 valores; variación dentro de 8 estudios | Distinguir asociación entre estudios de contraste ambiental dentro de estudio. |
| pH cartográfico | 172 perfiles; 21 valores; variación dentro de 8 estudios | Revisar codificación antes de interpretar unidades. |
| BIO1 carbono y pH conjuntamente | 172 perfiles; 75 combinaciones estudio coordenada; 53 grupos de coordenadas | La disponibilidad conjunta no equivale a 172 réplicas ambientales. |
| BioSamples repetidos | Dos identificadores asociados a tres runs cada uno | Resolver bibliotecas o muestreos relacionados; modelar dependencia sin borrar los registros originales. |
| Rendimiento de anotación | 1.212 a 2.019.349 asignaciones KO por perfil | No comparar conteos brutos como si fueran intensidad biológica. |

El diseño descriptivo a nivel estudio coordenada para BIO1 carbono y pH tiene rango 39 con estudio y 42 al añadir esas tres variables, en 75 filas. Hay información algebraica adicional para estimarlas, pero esto no garantiza potencia, independencia ni estimabilidad para cada KO. Cada gen tendrá su propio número de valores positivos y su propio diagnóstico de rango.

Todos los estudios están contenidos en un solo realm. Añadir realm a todos los efectos fijos de estudio incrementa el rango en cero. Ese efecto condicional no es identificable. Las comparaciones regionales entre estudios siguen siendo posibles bajo supuestos de comparabilidad explícitos; el Pacífico oriental tropical está representado por un solo estudio y tendrá interpretación descriptiva restringida.

## Variables ambientales y significado ecológico

Los scripts del repositorio documentan WorldClim 2.1 y SoilGrids para covariables cartográficas. La versión histórica de la tabla no conserva en cada columna todos los detalles de extracción necesarios. Por tanto, la procedencia general está orientada por código, pero versión, resolución, profundidad efectiva, desplazamiento de píxel y conversiones deben reconciliarse con los archivos originales de esa tabla.

BIO1 representa temperatura media anual de una climatología; WorldClim 2.1 corresponde a 1970–2000. No es la temperatura del sedimento al muestrear [21]. SoilGrids estima propiedades del suelo a partir de modelos y covariables; incluye carbono orgánico y pH, con incertidumbre espacial. No mide carbono disuelto, lábil o disponible para los microorganismos del sedimento [22]. El uso de clima dentro de los modelos de SoilGrids también impide interpretar sus productos y WorldClim como evidencias mecanísticas totalmente independientes.

La tabla histórica contiene phh2o_mean entre 47 y 79 y soc_mean entre 47 y 1.215. Son compatibles con pH multiplicado por diez y carbono en dg/kg, las codificaciones documentadas por ISRIC [23]. Una conversión confirmada daría pH 4,7–7,9 y carbono 4,7–121,5 g/kg. Son conversiones candidatas, no cambios aplicados silenciosamente. BIO1 ya presenta 21,96–28,47, compatible con grados Celsius; no debe dividirse nuevamente por diez de manera automática.

Human Footprint presenta valores distintos en seis grupos de coordenadas redondeadas. Se revisarán origen, resolución y extracción antes de resumirlo como propiedad única del sitio. Los índices de GBIF describen registros de ocurrencia y están condicionados por esfuerzo de observación; no equivalen a un inventario vegetal tomado junto al sedimento.

En la tabla examinada no hay columnas estructuradas de salinidad, carbono disuelto o disponible, potencial redox, oxígeno, temperatura del sedimento, fecha de muestreo ni profundidad. Esto no demuestra que nunca se midieran: podrían recuperarse de estudios y suplementos. El flujo admite incorporarlas con unidad, método, profundidad, fecha y fuente. Las asociaciones con carbono disponible quedan condicionadas a esa recuperación; las asociaciones con carbono orgánico cartográfico pueden explorarse ahora con su nombre correcto.

Se priorizarán BIO1, carbono orgánico cartográfico y pH una vez verificadas sus unidades. Precipitación y estacionalidad se evaluarán como familias secundarias con predictores representativos o ejes de PCA interpretados por sus cargas. No se introducirán los 19 BIO y todas las propiedades del suelo en una sola fórmula. La selección inicial se basará en significado, redundancia y soporte, no en el menor valor P.

## Preguntas e hipótesis de trabajo

H1 plantea que un repertorio funcional ampliamente compartido coexiste con variación geográfica en genes de ramas concretas de los ciclos. Evaluaremos prevalencia observada, representación relativa y variabilidad entre sitios. La detección ubicua no implica esencialidad; los conteos históricos no permiten afirmar redundancia de portadores.

H2 plantea que parte de esa variación se asocia con gradientes ambientales, incluso dentro de estudios con varios sitios. Un ejemplo dirigido es contrastar genes de formación de amonio y pasos de desnitrificación frente a carbono orgánico y pH. El sentido esperado se dejará como hipótesis, porque disponibilidad de nitrato, oxígeno y otros factores no medidos pueden cambiarlo.

H3 compara el recambio taxonómico y funcional en las mismas muestras. Podrían existir comunidades taxonómicamente diferentes con perfiles funcionales parecidos, o diferencias en ambos niveles. Esa comparación describirá concordancia o desacoplamiento entre perfiles; no asignará un gen a un taxón ni demostrará resiliencia. Para portadores se necesitará vínculo de secuencia.

H4 examina si los resultados se transfieren a un estudio o sitio no usado para ajustar el modelo. Un patrón visual que desaparece al retirar un estudio se reportará como dependiente de ese estudio, aunque tenga un valor P pequeño en el conjunto completo. No se fijará una historia de captura de carbono, emisiones o recuperación antes de observar resultados compatibles.

## Construcción de las respuestas funcionales

La primera respuesta será la representación relativa de cada KO: asignaciones a ese KO divididas por el total de asignaciones KO del perfil. El denominador se calculará sobre la matriz completa, antes de seleccionar ciclos o aplicar filtros por gen. Para la figura puede expresarse por millón; no se denominará TPM ni abundancia celular. La proporción es condicional al conjunto anotado y no incluye secuencias no anotadas que no están representadas en esta tabla.

Se mantendrán por separado el valor observado, el estado de detección y la disponibilidad de la anotación. Cero observado, fila no disponible y asignación ambigua tendrán estados distintos. Una sensibilidad con umbral de dos registros evaluará dependencia de detecciones únicas; no será un límite de detección validado. Los valores positivos se conservarán sin pseudoconteo para el modelo de representación condicional.

El catálogo de reacciones incluirá reaction_id, ciclo, sustrato, producto, gen, KO, regla, alternativa enzimática, alcance de dirección, referencia y evidencia requerida. La flecha describe la reacción de referencia; el color describe el gen. Una familia ambigua no generará automáticamente una flecha doble. En NifHDK mostraremos las tres casillas y el estado de la regla; no sumaremos las subunidades como medida de flujo.

Los procesos con varias arquitecturas se conservarán separados. Un índice de conjunto solo se añadirá si su fórmula se justifica y se mantiene la información por miembro. Los cocientes entre grupos, por ejemplo marcadores de formación de amonio frente a pasos de desnitrificación, serán contrastes de representación genética. No se interpretarán como tasas, eficiencia, estequiometría ni balance de gases. No se sumarán KOs compartidos varias veces en un total global de ciclo.

Para la comparación multivariada amplia se emplearán log cocientes con un conjunto de referencia documentado y manejo explícito de ceros, más una sensibilidad Hellinger o presencia ausencia. La selección de referencia y el reemplazo de ceros se fijarán antes de probar los gradientes; se evaluará cuánto cambian los resultados. Un CLR de cada ciclo responde una pregunta interna a ese ciclo, distinta del CLR calculado sobre un fondo funcional común.

## Modelos de asociación entre genes y ambiente

La estrategia principal candidata es un modelo en dos partes: probabilidad de detectar el KO y representación relativa entre perfiles donde se detecta. MaAsLin 3 permite ambas respuestas [11]. La primera se describirá como detección observada, no ocupación biológica verdadera; la segunda como cambio relativo del repertorio recuperado, no aumento de concentración de enzima.

El estimando ambiental prioritario es el cambio en representación relativa de un KO asociado con un gradiente entre sitios comparables, ajustado por estudio. La fórmula inicial contendrá los predictores ambientales verificados y efectos fijos de estudio. Se añadirá dependencia de sitio o BioSample solo cuando la estructura pueda estimarse; si esa parametrización es singular se usará una unidad agregada validada y se informará el cambio de estimando. No se eliminarán registros del inventario de 191.

El modelo con estudio fijo identifica la asociación a partir de variación dentro de estudios. Se mostrará junto al patrón entre sitios y estudios, que tiene una interpretación más descriptiva. Cuando haya suficiente información se estimarán pendientes por estudio y heterogeneidad; un efecto aleatorio de estudio no elimina por sí solo confusión entre ambiente y métodos. No se incluirá realm junto a todos los estudios fijos para reclamar un efecto regional separado.

Las variables ambientales se centrarán y escalarán para comparar efectos, conservando las unidades originales en las curvas. Se examinará la dependencia espacial residual. Si persiste, se contrastará una representación espacial de baja dimensión o incertidumbre por bloques, siempre que haya suficientes unidades y rango. Si ambiente, espacio y estudio no pueden separarse, se reportará esa asociación compartida sin atribuirla exclusivamente al ambiente.

Configuración propuesta de MaAsLin 3: entrada relativa calculada con el denominador completo, normalization = NONE, transform = LOG, median_comparison_abundance = FALSE, median_comparison_prevalence = FALSE y subtract_median = FALSE. La intención es contrastar asociaciones relativas frente a cero, no estimar abundancias absolutas bajo la hipótesis de una mediana de efectos nula. Los argumentos se confirmarán contra la versión fijada y se ensayarán en un conjunto sintético antes de ejecutar [24].

Para la detección se evaluará el total de asignaciones como indicador imperfecto de recuperación, además del estudio. Para la representación relativa se usará como diagnóstico y sensibilidad, porque ya participa en el denominador y su ajuste mecánico puede inducir asociaciones. Si se recuperan longitud ensamblada, cobertura o fracción anotada, se preferirán controles de medición con una justificación explícita. Ninguna normalización elimina todo el sesgo de ensamblaje.

Un filtro operativo inicial exigirá para estimar representación al menos diez grupos de sitio provisionales con detección y al menos tres estudios; para detección, suficientes grupos con presencia y con cero. Estos son pisos de admisión que deberán aumentar si la fórmula exige más información, no garantías de validez. Cada ajuste verificará rango, soporte del predictor, grados de libertad, convergencia, separación, observaciones influyentes y residuos. Los genes no estimables permanecerán en el atlas descriptivo.

La familia primaria de pruebas comprenderá KOs únicos del panel curado por predictores ambientales prioritarios y las dos respuestas, con corrección BH conjunta de los valores P válidos. La exploración del resto del repertorio tendrá una familia secundaria declarada. Se reportarán también los fallos y los tamaños muestrales de cada ajuste. No se seleccionarán ciclos después de observar resultados para reducir artificialmente la corrección.

Para unos pocos candidatos se comparará una tendencia lineal con una curva GAM de baja complejidad, ajustada por REML y validada fuera de sitios o estudios [14]. No se declararán umbrales ecológicos porque una curva flexible dibuje una inflexión. Se presentará el cambio entre valores ambientales observados, su intervalo de incertidumbre y los datos que sostienen esa zona de la curva.

## Biogeografía del repertorio completo y comparación con taxonomía

Se describirán el perfil funcional global y los perfiles de cada ciclo con ordenaciones reproducibles. RDA o dbRDA examinarán la relación con conjuntos ambientales reducidos, espacio y estudio, según escala de la respuesta. Se reportará varianza ajustada y las fracciones compartidas; una fracción ambiental no demuestra un mecanismo causal [13]. El número de ejes espaciales se limitará por el número de sitios y por rango, no por el número de perfiles.

Las coordenadas mundiales no se introducirán como distancias euclidianas en grados. Se usarán distancias geodésicas y una representación espacial compatible, con sensibilidad regional cuando las separaciones oceánicas hagan poco interpretable una única estructura continua. La distancia geográfica será contexto de distribución, no una medida directa de dispersión por corrientes o conectividad de manglar.

Las permutaciones deberán respetar la unidad del contraste. Permutar dentro de un sitio no prueba una variable ambiental constante en ese sitio. Para covariables de sitio se intercambiarán unidades completas dentro de estudios solo cuando exista intercambiabilidad; con muy pocas permutaciones distintas se priorizarán efectos e incertidumbre. No se tratarán todas las parejas de distancias como réplicas independientes ni se usará Mantel parcial como control universal del espacio.

Se cotejarán perfiles taxonómicos y funcionales en el mismo conjunto de IDs mediante ordenaciones emparejadas, Procrustes y modelos ambientales comparables. La bandera de disponibilidad taxonómica es positiva en los 191 metadatos, pero las matrices taxonómicas concretas deberán verificarse antes de este paso. Diferencias de escalas, filtrado y cuantificación impedirán comparar índices de diversidad como si midieran lo mismo.

Los modelos predictivos serán secundarios: primero un modelo ambiental sencillo y un modelo basal; solo después un modelo flexible si mejora el error fuera de estudio. El modelo de asociación con efectos fijos de estudio no se usará sin cambios para predecir un estudio desconocido: el modelo predictivo deberá producir una predicción sin disponer de su coeficiente de estudio. Toda imputación, transformación, selección y ajuste se hará dentro del entrenamiento. Se dejarán fuera estudios completos y se añadirá validación espacial o por componentes que mantenga juntos estudios y sitios compartidos [15]. Si ese cruce deja muy pocos bloques útiles, se declarará la transferencia no estimable. No se producirán superficies globales fuera del dominio ambiental muestreado.

## Atlas y figuras que cuentan la historia

El atlas mantendrá una geometría fija por ciclo. El mapa seleccionará localidad, estudio y muestra; las casillas junto a las reacciones mostrarán cada gen con una escala común entre localidades. Se conservará una vista de valores individuales y otra de resúmenes por unidad de muestreo validada. La variación, el número de muestras y los estudios de origen serán visibles. Los siete perfiles sin coordenadas continuarán en tablas y diagramas.

La vista observada mostrará representación relativa y detección. La vista ambiental coloreará los mismos genes por el efecto ajustado de carbono, temperatura u otra variable. Una tercera vista comparará predicciones en dos valores ambientales dentro del soporte observado, manteniendo las demás covariables bajo una definición común. Nunca se mezclarán en una misma leyenda abundancia, coeficiente de asociación y certeza funcional.

El manuscrito podría organizarse en cuatro figuras: distribución y soporte del muestreo; recambio funcional global y concordancia con taxonomía; ciclos con casillas por gen y contraste entre localidades; efectos ambientales con curvas, intervalos y estabilidad al retirar estudios. Las tablas completas quedarían como resultados asociados, conservando también genes sin asociación y ajustes fallidos.

La historia de cada candidato seguirá cinco preguntas. Qué gen o conjunto cambió. Qué reacción representa con certeza suficiente. Dónde se observa. Con qué gradiente se asocia tras el ajuste. Qué explicación alternativa sigue siendo compatible. Por ejemplo, una asociación de nrfA con carbono orgánico cartográfico podría sugerir un cambio en el potencial de formación de amonio; no demostraría disponibilidad de carbono ni una tasa mayor de DNRA sin las mediciones pertinentes.

La red principal será el mapa bioquímico curado. Las redes de covariación, si se conservan, serán exploratorias y secundarias: no se usarán para afirmar cooperación, transferencia de metabolitos ni estabilidad ecológica. La coherencia de varias piezas de una ruta será evidencia convergente de potencial, sin asumir que pertenecen al mismo organismo.

## Robustez y condiciones para interpretar

Se contrastarán resultados con detección de uno y dos registros, referencias composicionales alternativas justificadas, tratamiento de ceros, inclusión de perfiles de menor recuperación y retirada de estudios completos. Se comprobará si el signo y el intervalo del efecto se sostienen y si existe solapamiento ambiental entre grupos. No se eliminarán muestras por hacer desaparecer o aparecer una asociación.

Antes de considerar un candidato como hallazgo prioritario, deberá tener identidad funcional suficientemente específica, soporte ambiental no extrapolado, modelo válido, efecto e incertidumbre informados y estabilidad razonable frente a las sensibilidades pertinentes. La asociación entre estudios sin réplica interna se etiquetará como exploratoria. No se exigirá una narrativa positiva: un repertorio conservado o un efecto no identificable también son resultados informativos.

Se mantendrán por separado las incertidumbres de anotación, medición, muestreo y modelo. Disponer de más metagenomas puede mejorar el soporte, pero no corrige automáticamente una anotación ambigua, una escala ambiental incorrecta o la ausencia de mediciones de proceso. La actualización de Dani entrará como una nueva versión de datos con el mismo contrato; si cambia de conteos de proteínas a cobertura de genes, habrá que revisar también el modelo de medición.

## Implementación y criterios de avance

Etapa 1. Congelar entradas y reconciliar unidades, fuentes ambientales, BioSamples y grupos de muestreo. Mantener los 191 registros y producir un informe de disponibilidad por objetivo. Esta entrega ya incluye el auditor 12A, ejecutado sobre la matriz actual, y cuatro pruebas sintéticas de integridad.

Etapa 2. Completar el catálogo de reacciones direccionales a partir de las reglas existentes. Hacer un piloto del nitrógeno que muestre valores por gen, incertidumbre de anotación y diferencias entre observación y efecto ambiental. El piloto incluirá sitios de distintos estudios y niveles de recuperación, elegidos por diseño, no por resultados significativos.

Etapa 3. Verificar el motor estadístico con datos sintéticos que reproduzcan estructura de estudios, dependencia y ceros. Evaluar errores de tipo I y sensibilidad a confusión. Ejecutar después el panel dirigido y, con los mismos contratos, la exploración amplia. Registrar versiones de paquetes y bases, semillas, fórmulas y diagnósticos.

Etapa 4. Extender las figuras a carbono, metano, azufre, fósforo y los procesos de hierro que sean evaluables. Integrar modelos globales, comparación taxonómica y validación geográfica. La prioridad será explicar conjuntos coherentes de reacciones y mostrar efectos transferibles, no maximizar el número de genes con significación.

Etapa 5. Repetir con la matriz actualizada, manteniendo un manifiesto de muestras comunes y nuevas. Comparar efectos en el subconjunto compartido para separar cambios de cuantificación o anotación de cambios en cobertura geográfica. La repetición sobre las mismas muestras no se presentará como validación independiente.

La reconstrucción de portadores y resolución de familias ambiguas utilizará secuencias de proteínas, genes nucleotídicos y contigs; MAGs cuando estén disponibles. Para cuantificación más sólida se necesitarán longitudes génicas y mapeo de lecturas a una referencia común con manejo explícito de multimapeo y fragmentación. Para afirmar actividad o balances biogeoquímicos se requerirán mediciones adicionales; el análisis actual conserva su utilidad como potencial genético.

## Árbol de decisiones aplicado

| Nodo | Decisión actual | Estado |
| --- | --- | --- |
| D0 Datos | WGS de ADN resumido en conteos de proteínas asignadas a KO; 191 perfiles | Viable |
| D1 Objetivo | Describir geografía y asociar repertorio con ambiente; predicción secundaria | Viable con alcance definido |
| D2 Diseño | 44 estudios y agrupaciones geográficas provisionales; dependencia por resolver | Condicional para inferencia |
| D3 Medición | Normalización relativa explícita y auditoría de codificaciones ambientales | Condicional antes del ajuste |
| D4 Ruta | Matriz para primera fase; secuencias y lecturas para resolución y cuantificación posterior | Viable ahora |
| D5 Función | Genes por reacción con evidencia direccional y estados de información separados | Viable o condicional por reacción |
| D6 Estadística | Modelos por KO ajustados por estudio y modelos multivariados con rango verificado | Condicional por gen y contraste |
| D7 Robustez | Sensibilidades de medición y retirada de estudios; validación agrupada | Diseñada pendiente de ejecución |
| D8 Entrega | Código y definiciones en GitHub; matrices figuras y documentos fuera del repositorio | Aplicado a esta entrega |

No estimable con la matriz actual sola: concentración o actividad de enzimas, flujo neto, portadores taxonómicos, variantes de genes y asociación con carbono disponible no medido. Condicional: dirección de Mcr, Dsr y otras familias ambiguas; efectos ambientales con escaso contraste interno; predicción a nuevos dominios. Viable: atlas descriptivo, perfiles relativos y diseño reproducible para contrastes que superen los controles.

## Reproducción de la auditoría

Desde la raíz del repositorio, ejecutar con las entradas externas originales y una carpeta de resultados nueva:

```bash
python scripts/12_gene_environment/12A_AUDIT_GENE_ENVIRONMENT_DESIGN.py \
  --metadata /ruta/metadata_functional.csv \
  --matrix /ruta/functional_KO_counts_original.tsv.gz \
  --samples config/samples/functional_accessions.txt \
  --unit predicted_protein_records \
  --outdir /ruta/nueva/auditoria_ambiente
```

La salida incluye soporte por variable, variación dentro de estudio, conflictos dentro de coordenadas, BioSamples repetidos, rangos de diseño y hashes de entradas y código. No corrige datos, no imputa valores ni estima asociaciones. Los cuatro controles sintéticos se ejecutan con unittest sobre tests/test_gene_environment_audit.py. Las definiciones analíticas están en config/gene_environment/analysis_plan.json. Las huellas de archivos revisados de herramientas están en repository_review.json; son hashes de archivos concretos, no versiones de todo el programa.

## Referencias y documentación verificada

[1] Vernette C y colaboradores. The Ocean Gene Atlas v2.0 online exploration of the biogeography and phylogeny of plankton genes. Nucleic Acids Research. 2022;50:W516–W526. https://doi.org/10.1093/nar/gkac420. Código https://gitlab.osupytheas.fr/ocean_atlas/oga.

[2] Zhang Z y colaboradores. Global biogeography of microbes driving ocean ecological status under climate change. Nature Communications. 2024. https://doi.org/10.1038/s41467-024-49124-0. Código https://github.com/ZhenyanZhang/Ecological_status.

[3] Faure E y colaboradores. Towards omics based predictions of planktonic functional composition from environmental data. Nature Communications. 2021. https://doi.org/10.1038/s41467-021-24547-1. Código https://github.com/EmileFaure/MAGsProteinFunctionalClusters.

[4] Yue Z y colaboradores. CNPS.cycle streamlining shotgun metagenomic data analysis for biogeochemical cycles. mSystems. 2025;10:e01021-25. https://doi.org/10.1128/msystems.01021-25. Código https://github.com/yuezhengfu/CNPS.cycle.

[5] Xue CX y colaboradores. DiTing A Pipeline to Infer and Compare Biogeochemical Pathways From Metagenomic and Metatranscriptomic Data. Frontiers in Microbiology. 2021;12:698286. https://doi.org/10.3389/fmicb.2021.698286. Original https://github.com/xuechunxu/DiTing. Continuación señalada por el original https://github.com/SilentGene/DiTing.

[6] Zhou Z y colaboradores. METABOLIC high throughput profiling of microbial genomes for functional traits metabolism biogeochemistry and community scale functional networks. Microbiome. 2022;10:33. https://doi.org/10.1186/s40168-021-01213-8. Código https://github.com/AnantharamanLab/METABOLIC.

[7] Documentación de bigecyhmm y Tabigecy, consultada el 29 de octubre de 2026. https://github.com/ArnaudBelcour/bigecyhmm y https://github.com/ArnaudBelcour/tabigecy. Publicación vinculada por sus autores https://doi.org/10.1093/bioinformatics/btaf230.

[8] Shaffer M y colaboradores. DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Research. 2020. https://pmc.ncbi.nlm.nih.gov/articles/PMC7498326/. Desarrollo revisado https://github.com/BortonWrightonLabs/DRAM.

[9] ggkegg. Documentación y ejemplos oficiales. https://github.com/noriakis/ggkegg y https://noriakis.github.io/software/ggkegg. Consulta 29 de octubre de 2026.

[10] King ZA y colaboradores. Escher A Web Application for Building Sharing and Embedding Data Rich Visualizations of Biological Pathways. PLoS Computational Biology. 2015;11:e1004321. https://doi.org/10.1371/journal.pcbi.1004321. Desarrollo actual https://github.com/opencobra/escher.

[11] Nickols WA Kuntz T Shen J y colaboradores. MaAsLin 3 refining and extending generalized multivariable linear models for meta omic association discovery. Nature Methods. 2026;23:554–564. https://doi.org/10.1038/s41592-025-02923-9. Código https://github.com/biobakery/maaslin3.

[12] Lin H y Peddada SD. Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods. Publicación en línea 2023. https://doi.org/10.1038/s41592-023-02092-7. Código https://github.com/FrederickHuangLin/ANCOMBC.

[13] Documentación oficial de vegan y adespatial. https://github.com/vegandevs/vegan y https://github.com/adeverse/adespatial. Consulta 29 de octubre de 2026. Funciones y supuestos deben fijarse con la versión instalada.

[14] Wood SN. mgcv documentación del paquete y del modelo aditivo generalizado. https://stat.ethz.ch/CRAN/web/packages/mgcv/index.html. Código espejo de CRAN https://github.com/cran/mgcv. Consulta 29 de octubre de 2026.

[15] Valavi R Elith J Lahoz Monfort JJ y Guillera Arroita G. blockCV An R package for generating spatially or environmentally separated folds for k fold cross validation of species distribution models. Methods in Ecology and Evolution. 2019;10:225–232. https://doi.org/10.1111/2041-210X.13107. Código https://github.com/rvalavi/blockCV.

[16] microeco. Documentación del proyecto y tutorial oficial. https://github.com/ChiLiubio/microeco y https://chiliubio.github.io/microeco_tutorial/. Consulta 29 de octubre de 2026.

[17] Garber AI y colaboradores. FeGenie A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Frontiers in Microbiology. 2020;11:37. https://doi.org/10.3389/fmicb.2020.00037. Código https://github.com/Arkadiy-Garber/FeGenie.

[18] PCycDB a comprehensive and accurate database for fast analysis of phosphorus cycling genes. Microbiome. 2022. https://doi.org/10.1186/s40168-022-01292-1. Código https://github.com/ZengJiaxiong/Phosphorus-cycling-database.

[19] KEGGDecoder documentación y código original. https://github.com/bjtully/BioData/tree/master/KEGGDecoder. Consulta 29 de octubre de 2026.

[20] CNPS.cycle función Ncyc.abundance. https://github.com/yuezhengfu/CNPS.cycle/blob/main/R/Ncycleabundance.R. Archivo revisado con blob SHA bc5605c939030c851d39035344198144c0a0a851. Consulta 29 de octubre de 2026.

[21] WorldClim. Historical climate data versión 2.1 y definición de variables bioclimáticas. https://worldclim.org/data/worldclim21.html y https://worldclim.org/data/bioclim.html. Consulta 29 de octubre de 2026.

[22] ISRIC. SoilGrids Global gridded soil information. https://isric.org/explore/soilgrids. Consulta 29 de octubre de 2026.

[23] ISRIC. SoilGrids layers y factores de conversión. https://dev-docs.isric.containers.wur.nl/globaldata/soilgrids/SoilGrids_faqs_01.html. Consulta 29 de octubre de 2026.

[24] MaAsLin 3 manual de usuario y documentación de parámetros. https://github.com/biobakery/maaslin3 y https://www.bioconductor.org/packages/release/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html. Consulta 29 de octubre de 2026. El DESCRIPTION revisado indica 1.5.8 y R mayor o igual a 4.4; no se afirma que esa versión haya sido instalada aquí.
