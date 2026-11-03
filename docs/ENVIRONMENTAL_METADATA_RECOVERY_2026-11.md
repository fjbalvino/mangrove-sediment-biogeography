# Recuperación de campos ambientales existentes — 2–3 noviembre de 2026

Se revisó el universo de 401 muestras. Hay 399 con coordenadas; SRR31131347 y SRR31131364 se conservan en la base y se excluyen únicamente de esta extracción. El subconjunto funcional conserva sus 191 muestras. No se aplicó el manifiesto taxonómico de 341 ni se revisó la curación ecológica del usuario.

Se evaluaron todas las coordenadas con algún campo ambiental faltante: 62 muestras en 31 coordenadas. La actualización rellena celdas vacías de las 34 variables ambientales que ya existían. Conserva los nombres de archivo, 401 × 86 y 191 × 83 columnas/filas, el orden, las coordenadas, los identificadores, las clasificaciones y todos los valores previos. Los diagnósticos viven fuera de los CSV canónicos y fuera de git.

## Resultado de la ejecución

| Archivo | Clima: celdas | Huella: celdas | Suelo: celdas | Total | Muestras modificadas |
|---|---:|---:|---:|---:|---:|
| Metadata_reinicio_401_reconciliada.csv | 304 | 16 | 152 | 472 | 19 |
| metadata_functional.csv | 133 | 7 | 56 | 196 | 7 |

La segunda tabla es un subconjunto de la primera: no sumar sus celdas como observaciones independientes. Además de las 16 muestras georreferenciadas recientemente, se recuperaron los ocho campos de suelo para tres muestras que ya tenían coordenadas. No se recuperó GBIF: las 20 consultas por coordenada respondieron HTTP 503. Un error de servicio nunca se codifica como cero.

## Compatibilidad histórica y decisiones

1. **WorldClim v2.1, 1970–2000, resolución 30 segundos de arco.** Los controles identifican esta resolución, aunque un script antiguo declara 2,5 minutos. Se usan BIO1–BIO19 sin aplicar un factor adicional de 0,1 a las temperaturas. BIO4 mantiene su definición nativa. Se descartó el piloto de 2,5 minutos.
2. **Huella humana 2020: HII HF3-1, sin OSM.** Se compararon ambas variantes oficiales: HF3-1 coincide exactamente en 14/20 controles, frente a 1/20 para HF3-2. Se preserva la codificación numérica histórica. Las discrepancias restantes se conservan en el diagnóstico, no se sobrescriben.
3. **SoilGrids 2.0: media predictiva de la capa 0–5 cm, enteros nativos.** Coinciden exactamente los ocho atributos en 15/20 controles; pH en 16/20. El control de profundidades distingue 0–5 cm de un promedio 0–30 cm. Por tanto, estas columnas NO deben describirse como un promedio de 0–30 cm. Por ejemplo, pH 61 representa 6,1; el CSV conserva 61. Se usa WCS porque la API REST está pausada. Los recortes WCS devuelven ceros sin etiqueta NoData en el mar: se enmascaran mediante pH <= 0 antes de buscar una celda válida. No se incorporan ceros de relleno.
4. **GBIF:** se consultan también las coordenadas que ya existían. Los índices vacíos asociados a 12 registros históricos con cero observaciones son indefinidos, no un valor cero de diversidad. No se mezclan nuevos índices con un bloque histórico que conserva cero observaciones. La extracción actual de GBIF del repositorio usa rarefacción y adelgazamiento espacial; los scripts antiguos usan otras reglas, incluso radios crecientes. Hace falta reconciliar el método para toda la cohorte antes de mezclar resultados. El adaptador conserva las respuestas/errores y no incorpora GBIF en el parche de esta corrida.

Se utiliza la celda que contiene el punto cuando tiene dato. Si está vacía, se busca la celda válida más cercana a un máximo de 2.000 m, registrando centro y distancia. El límite no se amplía para forzar datos. Para suelo, todas las propiedades se extraen en el mismo anclaje. La búsqueda por centro de celda sustituye la selección antigua del primer píxel válido de una ventana. Los controles que difieren y las coordenadas aproximadas requieren análisis de sensibilidad posterior. Completar campos no convierte un proxy espacial en una medición del sedimento.

## Scripts revisados y usados

- 02A se ejecutó con el manifiesto de 399 accesiones, generado desde la base completa.
- Para 02B y 02C, los adaptadores `recover_rasters.py` y `recover_worldclim30.py` usan `metadata_pipeline_common.py`, preservan las unidades comprobadas y extraen las capas oficiales identificadas arriba.
- Para 02D, `recover_gbif.py` usa las constantes, limitador de consultas y ecuaciones del repositorio; registra errores sin rellenar la base.
- Para 02E, `recover_soil.py` sustituye únicamente el acceso REST por WCS y conserva la profundidad/codificación histórica comprobada. La opción 0–30 cm del script general no se aplicó a estas columnas.
- 02F se revisó: REALM, PROVINCE y ECOREGION estaban completos para las 399 muestras. No se reasignaron categorías curadas.
- 02G pasó con 401 y 191 filas, respectivamente. Su comprobación de completitud de 14 predictores es un diagnóstico técnico; no autoriza por sí misma un GAM ni resuelve los problemas de diseño.
- `merge_verified.py` exige campos existentes, celdas vacías, datos finitos, radio permitido, coherencia de textura, ausencia de errores y concordancia numérica entre la base completa y la funcional. También impone un bloqueo técnico ante pérdida general de coincidencia con controles históricos; ese umbral no es una prueba estadística de equivalencia.

## Reejecución

Dependencias: Python 3, numpy, pandas, rasterio, pyproj, requests; tqdm es opcional para los módulos originales. Mantener los datos y productos fuera del repositorio.

```bash
export RECOVERY_WORKDIR=/ruta/privada/recuperacion_ambiental
export METADATA_MASTER=/ruta/Metadata_reinicio_401_reconciliada.csv
export METADATA_FUNCTIONAL=/ruta/metadata_functional.csv
bash workflow/02_RECOVER_MISSING_ENVIRONMENTAL_FIELDS.sh
```

Antes de ejecutar, colocar en `$RECOVERY_WORKDIR/cache/` `hii_no_osm_2020-01-01.tif` y los mosaicos oficiales `tile_<ID>_wc2.1_30s_bio.tif` que cubran las coordenadas. Esta corrida usa 32, 33, 34 y 41; el último aporta controles. `download_asset.py URL DESTINO` descarga con HTTPS, registra URL, cabeceras y SHA-256, y reutiliza una descarga existente. Un archivo existente debe cotejarse con el manifiesto de procedencia conservado. No reemplazarlo por otra versión con el mismo nombre.

El flujo escribe copias con los nombres canónicos en `outputs/`; no modifica los archivos de entrada. La reposición de los archivos maestros y del CSV dentro del paquete de análisis ocurre solamente después de validar el parche. El paquete mantiene los resultados anteriores; los modelos y figuras dependientes de metadatos deberán volver a ejecutarse.

Validación de esta entrega: preservación de todas las celdas ajenas al parche, esquemas, orden, exclusiones y concordancia entre cohortes; repetición desde caché; rechazo de la variante incompatible de HII; 02G para ambos tamaños. Las tablas de cambios, pendientes, controles, errores de API y recortes originales se conservan en el paquete de evidencia privado. Este commit contiene exclusivamente código y documentación.

## Fuentes oficiales

- WorldClim v2.1: https://www.worldclim.org/data/worldclim21.html
- Distribución oficial de mosaicos, función `worldclim_tile`: https://github.com/rspatial/geodata/blob/master/R/worldclim.R
- HII y variantes HF3-1/HF3-2: https://wcshumanfootprint.org/data-access
- SoilGrids: https://docs.isric.org/globaldata/soilgrids/SoilGrids_faqs_02.html
- SoilGrids WCS: https://docs.isric.org/globaldata/soilgrids/WCS_from_Python.html
- Índice y VRT oficiales: https://files.isric.org/soilgrids/latest/data/
- GBIF API: https://techdocs.gbif.org/en/openapi/v1/occurrence
