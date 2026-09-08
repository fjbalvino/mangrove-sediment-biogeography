# Guía breve: GitHub, Zenodo y publicación

## Etapa 1 — Construcción privada

1. Crear el repositorio vacío en la cuenta `fjbalvino`.
2. Mantenerlo privado mientras se completa el bloque funcional 07.
3. Subir esta estructura a la rama `main`.
4. Confirmar que la validación automática quede en verde.
5. Trabajar mediante ramas cortas y pull requests; no editar directamente una
   versión ya publicada.

## Etapa 2 — Preparar los datos

1. Construir un paquete Zenodo sin lecturas FASTQ.
2. Incluir tablas de acceso, metadata curada, matrices, objetos canónicos,
   árbol GTDB, resultados principales y un manifiesto SHA-256.
3. Documentar licencia, versión, fecha, creador y procedencia de cada activo.
4. Probar la descarga y los checksums en un directorio limpio.

## Etapa 3 — Primera versión citable

1. Conectar el repositorio con Zenodo.
2. Cambiar el repositorio a público cuando autores y datos estén aprobados.
3. Crear una versión GitHub `v1.0.0` con notas claras.
4. Zenodo archivará esa versión y asignará un DOI.
5. Añadir el DOI a `README.md`, `CITATION.cff`, `DATA_AVAILABILITY.md` y al
   manuscrito.
6. Usar el DOI de versión para reproducir exactamente el trabajo y el DOI de
   concepto para referirse al software a través de futuras versiones.

## Regla práctica

- GitHub responde: **¿cómo se hizo?**
- Zenodo responde: **¿con qué archivos exactos?**
- SRA responde: **¿dónde están las lecturas originales?**

