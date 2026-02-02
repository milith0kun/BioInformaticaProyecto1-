# Informe Técnico: Sistema de Análisis Genómico de E. coli K-12 MG1655

## 1. Arquitectura del Sistema

### 1.1 Estructura de Directorios
```
proyecto-ecoli/
├── ncbi_dataset.zip                    # Archivo descargado de NCBI (raíz)
├── backend/
│   ├── app/
│   │   ├── main.py
│   │   ├── api/
│   │   │   └── routes/
│   │   │       ├── files.py
│   │   │       ├── analysis.py
│   │   │       └── export.py
│   │   ├── core/
│   │   │   ├── file_detector.py
│   │   │   ├── genome_parser.py
│   │   │   └── analyzers/
│   │   │       ├── codon_analyzer.py
│   │   │       ├── gene_analyzer.py
│   │   │       └── validator.py
│   │   └── models/
│   │       └── schemas.py
│   ├── cache/
│   └── requirements.txt
├── frontend/
│   ├── src/
│   │   ├── components/
│   │   │   ├── FileManager.jsx
│   │   │   ├── AnalysisDashboard.jsx
│   │   │   ├── CodonVisualization.jsx
│   │   │   ├── GeneStatistics.jsx
│   │   │   └── DataExport.jsx
│   │   ├── services/
│   │   │   └── api.js
│   │   └── App.jsx
│   ├── package.json
│   └── vite.config.js
└── README.md
```

## 2. Stack Tecnológico

### 2.1 Backend
- **Framework**: FastAPI 0.109.0
- **Servidor**: Uvicorn con soporte estándar
- **Bioinformática**: BioPython 1.83
- **Procesamiento datos**: Pandas 2.2.0, NumPy 1.26.3
- **Visualización**: Matplotlib 3.8.2, Seaborn 0.13.1
- **Validación**: Pydantic 2.5.3
- **CORS**: fastapi.middleware.cors

### 2.2 Frontend
- **Framework**: React 18.2.0
- **Bundler**: Vite 5.0.8
- **Cliente HTTP**: Axios 1.6.5
- **Gráficos**: Recharts 2.10.3, D3.js 7.8.5
- **Tablas**: AG-Grid React 31.0.1
- **UI**: Tailwind CSS 3.4.1, HeadlessUI 1.7.18
- **Iconos**: Heroicons 2.1.1
- **Notificaciones**: React Hot Toast 2.4.1
- **Exportación**: jsPDF 2.5.1, XLSX 0.18.5
- **Utilidades**: Lodash 4.17.21, PapaParse 5.4.1

### 2.3 Librerías de Bioinformática Específicas

#### BioPython 1.83 - Módulos Utilizados

**Módulo Bio.SeqIO**
- Parsing archivo GenBank (.gbff)
- Lectura de secuencias en formato FASTA
- Función: `SeqIO.parse()` y `SeqIO.read()`
- Ubicación: `backend/app/core/genome_parser.py`

**Módulo Bio.Seq**
- Representación de secuencias como objetos Seq
- Búsqueda de patrones (codones ATG, TAA, TAG, TGA)
- Función: `seq.count()`, `seq.find()`
- Ubicación: `backend/app/core/analyzers/codon_analyzer.py`

**Módulo Bio.SeqUtils**
- Cálculo automático de contenido GC
- Función: `gc_fraction(sequence)`
- Ubicación: `backend/app/core/genome_parser.py`

**Módulo Bio.SeqFeature**
- Extracción de features (genes, CDS)
- Parsing de anotaciones GenBank
- Obtención de qualifiers (locus_tag, product)
- Función: `feature.extract()`, `feature.qualifiers`
- Ubicación: `backend/app/core/genome_parser.py`

**Módulo Bio.SeqRecord**
- Gestión completa de información genómica
- Acceso a anotaciones del archivo GenBank
- Ubicación: `backend/app/core/genome_parser.py`

**Módulo Bio.Entrez (opcional)**
- Validación de metadata del genoma
- Verificación de organismo correcto
- Ubicación: `backend/app/core/analyzers/validator.py`

## 3. Funcionalidades del Sistema

### 3.1 Detección y Carga de Archivos
El backend debe automáticamente:
- Extraer `ncbi_dataset.zip` ubicado en la raíz del proyecto
- Identificar archivos por extensión: `.gbff`, `.fna`, `.gff`, `.gtf`, `.fasta`, `.jsonl`
- Priorizar archivo GenBank (`.gbff`) como fuente principal
- Cachear el parsing para evitar reprocesamiento
- Exponer endpoint que liste archivos detectados con metadata

### 3.2 Análisis de Codones
Implementar:
- Búsqueda exhaustiva de codón ATG en secuencia completa
- Conteo de codones de terminación: TAA, TAG, TGA
- Cálculo de densidad por kilobase (ATG/kb)
- Proporción relativa de cada codón de terminación
- Comparación con cantidad de genes anotados
- Generación de JSON y CSV con resultados

### 3.3 Análisis de Genes
Extraer del archivo GenBank:
- Número total de genes anotados
- Número de CDS (coding sequences)
- Longitud de cada gen
- Contenido GC global y por gen
- Posiciones genómicas (start, end, strand)
- Densidad génica (genes/Mb)
- Distribución estadística de tamaños

### 3.4 Validación Científica
Verificar contra valores de referencia:
- Total de genes: ~4,300
- Pares de bases: ~4.6 millones
- Contenido GC: ~50.8%
- Calcular desviaciones porcentuales
- Generar alertas si desviación > 5%

### 3.5 Visualizaciones Interactivas
React debe renderizar:
- **Gráfico de barras**: Distribución de codones ATG vs terminación
- **Histograma**: Distribución de tamaños de genes
- **Gráfico de dispersión**: Contenido GC por gen
- **Tabla filtrable**: Lista completa de genes con búsqueda
- **Panel estadístico**: Métricas clave con tarjetas
- **Visor comparativo**: Valores calculados vs referencia

### 3.6 Exportación de Datos
Permitir descargar:
- Resultados completos en JSON
- Tablas en CSV y Excel
- Informe en PDF con gráficos embebidos
- Gráficos individuales en PNG

## 4. API REST Endpoints

### 4.1 Gestión de Archivos
```
GET  /api/files/detect          # Lista archivos disponibles
POST /api/files/extract          # Extrae ZIP si no está extraído
GET  /api/files/{filename}/info  # Metadata de archivo específico
```

### 4.2 Análisis
```
POST /api/analysis/codons        # Ejecuta análisis de codones
POST /api/analysis/genes         # Ejecuta análisis de genes
GET  /api/analysis/validate      # Compara con valores referencia
GET  /api/analysis/complete      # Ejecuta todos los análisis
GET  /api/analysis/status        # Estado de análisis en curso
```

### 4.3 Resultados
```
GET /api/results/codons          # Resultados análisis codones
GET /api/results/genes           # Resultados análisis genes
GET /api/results/statistics      # Estadísticas generales
```

### 4.4 Exportación
```
GET /api/export/json             # Exportar todo en JSON
GET /api/export/csv/{type}       # CSV por tipo de análisis
GET /api/export/pdf              # Generar informe PDF
```

## 5. Flujo de Ejecución

### 5.1 Inicio del Sistema
1. Backend inicia en puerto 8000
2. Verifica existencia de `ncbi_dataset.zip` en raíz
3. Si existe, extrae contenido a `backend/extracted/`
4. Escanea archivos extraídos
5. Identifica archivo GenBank principal
6. Frontend inicia en puerto 5173
7. React conecta con backend vía Axios

### 5.2 Proceso de Análisis
1. Usuario accede a interfaz React
2. Frontend solicita lista de archivos disponibles
3. Backend retorna archivos con metadata
4. Usuario selecciona archivo GenBank o usa default
5. Frontend inicia análisis completo
6. Backend ejecuta en paralelo:
   - Parsing del GenBank con BioPython
   - Análisis de codones ATG y terminación
   - Extracción de información génica
   - Cálculos estadísticos
   - Validación contra referencia
7. Resultados se cachean en `backend/cache/`
8. Frontend recibe datos y renderiza visualizaciones
9. Usuario interactúa con gráficos y tablas
10. Usuario exporta resultados en formato deseado

## 6. Requisitos de Implementación

### 6.1 Backend
- Manejo de archivos ZIP con `zipfile`
- Parser GenBank usando `Bio.SeqIO`
- Búsqueda de patrones con expresiones regulares
- Procesamiento de DataFrames con Pandas
- Generación de gráficos con Matplotlib (para PDF)
- Cache de resultados con hash de archivo
- Endpoints asíncronos para análisis pesados
- CORS configurado para localhost:5173

### 6.2 Frontend
- Componentes React funcionales con hooks
- Estado global con Context API o hooks personalizados
- Llamadas API con Axios interceptors
- Visualizaciones responsivas con Recharts
- Tablas con filtrado y ordenamiento (AG-Grid)
- Loading states y error handling
- Notificaciones de progreso
- Diseño responsive con Tailwind

## 7. Algoritmos Clave con BioPython

### 7.1 Conteo de Codones (con BioPython)

```python
from Bio import SeqIO
from Bio.Seq import Seq

# Cargar genoma usando Bio.SeqIO
record = SeqIO.read("genome.gbff", "genbank")
sequence = record.seq  # Objeto Bio.Seq

# Buscar ATG usando Bio.Seq.count()
atg_count = sequence.count("ATG")

# Buscar codones stop
taa_count = sequence.count("TAA")
tag_count = sequence.count("TAG")
tga_count = sequence.count("TGA")

# Calcular densidad (por kilobase)
genome_length = len(sequence)
atg_density = (atg_count / genome_length) * 1000
```

### 7.2 Extracción de Genes (con BioPython)

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Cargar genoma
record = SeqIO.read("genome.gbff", "genbank")
sequence = record.seq

genes_data = []
for feature in record.features:
    if feature.type == "gene" or feature.type == "CDS":
        # Extraer secuencia del gen usando BioPython
        gene_seq = feature.extract(sequence)
        
        gene_info = {
            "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "length": len(gene_seq),
            "strand": feature.location.strand,
            # Calcular GC usando Bio.SeqUtils.gc_fraction
            "gc_content": gc_fraction(gene_seq) * 100,
            "product": feature.qualifiers.get("product", [None])[0]
        }
        genes_data.append(gene_info)
```

### 7.3 Validación
- Comparar valores calculados con constantes de referencia
- Calcular: desviación = abs(calculado - referencia) / referencia * 100
- Clasificar: PASS (<5%), WARNING (5-10%), FAIL (>10%)

## 8. Formato de Datos

### 8.1 Respuesta Análisis Codones
```json
{
  "genome_length": 4641652,
  "atg_count": 4498,
  "atg_density": 0.969,
  "stop_codons": {
    "TAA": {"count": 2850, "percentage": 63.2},
    "TAG": {"count": 326, "percentage": 7.2},
    "TGA": {"count": 1332, "percentage": 29.6}
  },
  "gene_comparison": {
    "annotated_genes": 4321,
    "atg_found": 4498,
    "difference": 177
  }
}
```

### 8.2 Respuesta Análisis Genes
```json
{
  "total_genes": 4321,
  "total_cds": 4140,
  "genome_length": 4641652,
  "gc_content": 50.79,
  "gene_density": 931.2,
  "size_statistics": {
    "mean": 950.3,
    "median": 846,
    "min": 51,
    "max": 7320
  },
  "genes": [
    {
      "locus_tag": "b0001",
      "start": 190,
      "end": 255,
      "length": 66,
      "strand": 1,
      "product": "thrL",
      "gc_content": 48.5
    }
  ]
}
```

## 9. Consideraciones de Rendimiento

- Cachear resultados de análisis con timestamp
- Lazy loading de tabla de genes
- Virtualización de listas largas
- Compresión de respuestas JSON
- Paginación en endpoints de genes
- Web Workers para procesamiento pesado en frontend
- Streaming de análisis de archivos grandes

## 10. Entregables del Proyecto

### 10.1 Repositorio GitHub
- Código fuente completo
- README con instrucciones de instalación
- requirements.txt y package.json
- .gitignore configurado
- Documentación de API

### 10.2 Informe Científico
- Introducción: E. coli como modelo, importancia del análisis
- Metodología: descripción técnica completa
- Resultados: tablas, gráficos, estadísticas
- Discusión: interpretación biológica, validación
- Conclusiones: hallazgos principales

### 10.3 Aplicación Funcional
- Backend ejecutable con un comando
- Frontend ejecutable con un comando
- Interfaz completamente funcional
- Todas las visualizaciones operativas
- Exportación funcionando

Este informe proporciona la especificación completa para que una IA pueda estructurar e implementar el sistema completamente funcional sin información irrelevante.