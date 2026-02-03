# 🧬 Sistema de Análisis Genómico de E. coli K-12 MG1655

> **Plataforma web completa para análisis bioinformático** con búsqueda en NCBI, visualizaciones interactivas y validación con IA.

![Python](https://img.shields.io/badge/Python-3.12-blue?logo=python)
![FastAPI](https://img.shields.io/badge/FastAPI-0.128-green?logo=fastapi)
![React](https://img.shields.io/badge/React-18.2-61DAFB?logo=react)
![Claude AI](https://img.shields.io/badge/Claude_AI-3.5_Haiku-purple?logo=anthropic)

## 📋 Tabla de Contenidos

1. [Descripción](#-descripción)
2. [Características](#-características)
3. [Arquitectura](#-arquitectura)
4. [Instalación](#-instalación)
5. [Uso](#-uso)
6. [API Endpoints](#-api-endpoints)
7. [Stack Tecnológico](#-stack-tecnológico)
8. [Tests](#-tests)
9. [Estructura del Proyecto](#-estructura-del-proyecto)
10. [Validación Científica con IA](#-validación-científica-con-ia)
11. [Contribución](#-contribución)

---

## 📖 Descripción

Sistema de análisis genómico diseñado para estudiar el genoma de **Escherichia coli K-12 MG1655** (RefSeq: NC_000913.3), el organismo modelo más importante en biología molecular.

### ¿Qué hace este sistema?

1. **Busca y descarga genomas** directamente desde NCBI Datasets API v2
2. **Analiza codones** (ATG, TAA, TAG, TGA) con estadísticas de densidad
3. **Extrae información génica** (4,651 genes, ~4,318 CDS)
4. **Valida resultados científicamente** usando Claude AI (Anthropic)
5. **Visualiza datos** con gráficos interactivos
6. **Exporta resultados** en JSON, CSV y PDF

---

## ✨ Características

### 🔬 Análisis Bioinformático
- **Análisis de Codones**: Conteo exhaustivo de codones de inicio (ATG) y terminación (TAA, TAG, TGA)
- **Análisis de Genes**: Extracción de información génica, estadísticas de tamaño, contenido GC
- **Densidad Génica**: Cálculo de genes/Mb y cobertura del genoma

### 🌐 Integración con NCBI
- Búsqueda de genomas por nombre de organismo
- Descarga automática desde NCBI Datasets API v2
- Soporte para múltiples genomas (E. coli, Bacillus, Saccharomyces, etc.)
- Lista de genomas populares pre-configurados

### 🤖 Validación con IA (Claude 3.5 Haiku)
- Validación científica de resultados
- Contexto biológico interpretativo
- Detección de discrepancias
- Recomendaciones para análisis adicionales

### 📊 Visualizaciones
- Gráficos de barras para distribución de codones
- Histogramas de tamaños de genes
- Gráficos de dispersión GC% vs longitud
- Tablas interactivas con AG Grid (búsqueda, paginación, ordenamiento)

### 📤 Exportación
- JSON completo
- CSV por tipo de análisis
- PDF con gráficos embebidos

---

## 🏗️ Arquitectura

```
┌─────────────────────────────────────────────────────────────┐
│                        FRONTEND                              │
│                    React + Vite + Tailwind                   │
│  ┌─────────────┐ ┌───────────────┐ ┌──────────────────────┐ │
│  │ GenomeSelector │ │AnalysisDashboard│ │ CodonVisualization │ │
│  └─────────────┘ └───────────────┘ └──────────────────────┘ │
│  ┌─────────────┐ ┌───────────────┐ ┌──────────────────────┐ │
│  │ GeneStatistics │ │ AIValidation  │ │   DataExport       │ │
│  └─────────────┘ └───────────────┘ └──────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼ Axios + Proxy
┌─────────────────────────────────────────────────────────────┐
│                        BACKEND API                           │
│                    FastAPI + Uvicorn                         │
│  ┌─────────────┐ ┌───────────────┐ ┌──────────────────────┐ │
│  │ /api/genome │ │ /api/analysis │ │    /api/export       │ │
│  └─────────────┘ └───────────────┘ └──────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
         │                   │                    │
         ▼                   ▼                    ▼
┌─────────────┐     ┌─────────────┐      ┌─────────────┐
│ NCBI API v2 │     │  BioPython  │      │  Claude AI  │
│  (Datasets) │     │  (Parsing)  │      │ (Anthropic) │
└─────────────┘     └─────────────┘      └─────────────┘
```

---

## 🚀 Instalación

### Prerrequisitos

- Python 3.10+
- Node.js 18+
- npm o yarn

### 1. Clonar repositorio

```bash
git clone https://github.com/tu-usuario/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis
```

### 2. Backend

```bash
cd backend

# Crear entorno virtual
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# o: venv\Scripts\activate  # Windows

# Instalar dependencias
pip install -r requirements.txt

# Configurar variables de entorno
cp .env.example .env
# Editar .env y agregar CLAUDE_API_KEY
```

### 3. Frontend

```bash
cd frontend
npm install
```

### 4. Iniciar servidores

**Terminal 1 - Backend:**
```bash
cd backend
source venv/bin/activate
uvicorn app.main:app --reload --port 8000
```

**Terminal 2 - Frontend:**
```bash
cd frontend
npm run dev
```

### 5. Acceder a la aplicación

Abrir navegador en: **http://localhost:5173**

---

## 📱 Uso

### Flujo de trabajo guiado (4 pasos)

1. **🔍 Buscar Genoma**
   - Escribe "Escherichia coli" en el buscador
   - O pega un accession number (ej: `GCF_000005845.2`)
   - O usa un genoma popular de la lista

2. **📥 Seleccionar**
   - Descarga el genoma desde NCBI
   - Selecciona de la lista de genomas descargados
   - Activa el genoma para análisis

3. **⚡ Analizar**
   - Haz clic en "Ejecutar Análisis Completo"
   - El sistema analiza codones y genes
   - Valida contra valores de referencia

4. **📊 Resultados**
   - Dashboard con métricas principales
   - Visualización de codones
   - Estadísticas de genes
   - Validación con IA
   - Exportar datos

---

## 🔌 API Endpoints

### Gestión de Genomas

| Método | Endpoint | Descripción |
|--------|----------|-------------|
| GET | `/api/genome/search/{query}` | Buscar genomas en NCBI |
| GET | `/api/genome/info/{accession}` | Info de genoma específico |
| POST | `/api/genome/download` | Descargar genoma |
| GET | `/api/genome/downloaded` | Listar descargados |
| POST | `/api/genome/activate/{accession}` | Activar para análisis |
| GET | `/api/genome/popular` | Genomas populares |

### Análisis

| Método | Endpoint | Descripción |
|--------|----------|-------------|
| POST | `/api/analysis/codons` | Analizar codones |
| POST | `/api/analysis/genes` | Analizar genes |
| GET | `/api/analysis/complete` | Análisis completo |
| GET | `/api/analysis/validate` | Validar vs referencia |
| POST | `/api/analysis/ai-validation` | Validar con Claude AI |
| GET | `/api/analysis/status` | Estado del análisis |

### Exportación

| Método | Endpoint | Descripción |
|--------|----------|-------------|
| GET | `/api/export/json` | Exportar JSON |
| GET | `/api/export/csv/{type}` | Exportar CSV |
| GET | `/api/export/pdf` | Generar PDF |

---

## 🔧 Stack Tecnológico

### Backend

| Tecnología | Versión | Uso |
|------------|---------|-----|
| FastAPI | 0.128.0 | Framework API REST |
| Uvicorn | 0.27.0 | Servidor ASGI |
| BioPython | 1.83 | Parsing genómico |
| Pandas | 2.2.0 | Procesamiento datos |
| Anthropic | 0.45.0 | Claude AI API |
| Pydantic | 2.5.3 | Validación datos |

### Frontend

| Tecnología | Versión | Uso |
|------------|---------|-----|
| React | 18.2.0 | UI Framework |
| Vite | 5.0.8 | Build tool |
| Tailwind CSS | 3.4.1 | Estilos |
| Recharts | 2.10.3 | Gráficos |
| AG Grid | 31.0.1 | Tablas |
| Axios | 1.6.5 | HTTP Client |

### APIs Externas

| Servicio | Uso |
|----------|-----|
| NCBI Datasets API v2 | Búsqueda y descarga de genomas |
| Claude AI (Anthropic) | Validación científica |

---

## 🧪 Tests

### Ejecutar tests del backend

```bash
cd backend
source venv/bin/activate
pip install pytest pytest-asyncio httpx
python -m pytest tests/ -v
```

### Cobertura de tests

```
✓ Health endpoints (2 tests)
✓ Genome endpoints (4 tests)
✓ Analysis endpoints (2 tests)
✓ Export endpoints (2 tests)
✓ Codon Analyzer (3 tests)
✓ Gene Analyzer (1 test)
✓ Validator (2 tests)
✓ NCBI Downloader (1 test)

Total: 17 tests
```

---

## 📁 Estructura del Proyecto

```
proyecto/
├── backend/
│   ├── app/
│   │   ├── main.py                 # FastAPI app
│   │   ├── api/
│   │   │   └── routes/
│   │   │       ├── files.py        # Gestión de archivos
│   │   │       ├── analysis.py     # Análisis genómico
│   │   │       ├── genome.py       # NCBI integration
│   │   │       └── export.py       # Exportación
│   │   ├── core/
│   │   │   ├── file_detector.py    # Detección de archivos
│   │   │   ├── genome_parser.py    # Parser GenBank
│   │   │   ├── ncbi_downloader.py  # Descarga NCBI
│   │   │   └── analyzers/
│   │   │       ├── codon_analyzer.py    # Análisis codones
│   │   │       ├── gene_analyzer.py     # Análisis genes
│   │   │       ├── validator.py         # Validación básica
│   │   │       └── ai_validator.py      # Validación IA
│   │   └── models/
│   │       └── schemas.py          # Pydantic schemas
│   ├── tests/
│   │   └── test_api.py            # Tests unitarios
│   ├── genomes/                   # Genomas descargados
│   ├── cache/                     # Cache de análisis
│   ├── requirements.txt
│   └── .env
│
├── frontend/
│   ├── src/
│   │   ├── App.jsx                # Componente principal
│   │   ├── components/
│   │   │   ├── GenomeSelector.jsx      # Búsqueda NCBI
│   │   │   ├── AnalysisDashboard.jsx   # Dashboard
│   │   │   ├── CodonVisualization.jsx  # Gráficos codones
│   │   │   ├── GeneStatistics.jsx      # Estadísticas genes
│   │   │   ├── AIValidation.jsx        # Validación IA
│   │   │   ├── DataExport.jsx          # Exportación
│   │   │   └── FileManager.jsx         # Gestión archivos
│   │   ├── services/
│   │   │   └── api.js             # Axios service
│   │   └── index.css              # Estilos Tailwind
│   ├── package.json
│   └── vite.config.js
│
└── README.md
```

---

## 🤖 Validación Científica con IA

### Configuración

1. Obtener API key de [Anthropic Console](https://console.anthropic.com/)
2. Agregar al archivo `.env`:
   ```
   CLAUDE_API_KEY=sk-ant-api03-xxxxx
   ```

### Qué valida la IA

- **Análisis de Codones**: Proporciones TAA/TAG/TGA, densidad ATG
- **Análisis de Genes**: Desviación vs referencia NC_000913.3
- **Validación Global**: Consistencia integral del análisis

### Campos de respuesta

| Campo | Descripción |
|-------|-------------|
| `is_valid` | Resultado de validación (true/false) |
| `confidence` | Nivel de confianza (0-100%) |
| `key_findings` | Hallazgos clave |
| `scientific_context` | Contexto biológico |
| `discrepancies` | Discrepancias detectadas |
| `recommendations` | Acciones recomendadas |

---

## 📊 Valores de Referencia

### E. coli K-12 MG1655 (NC_000913.3)

| Métrica | Valor |
|---------|-------|
| Tamaño del genoma | 4,641,652 bp |
| Total de genes | 4,651 |
| Total CDS | 4,318 |
| Contenido GC | 50.79% |
| Densidad génica | ~1,000 genes/Mb |

### Proporción típica de codones de terminación

| Codón | Porcentaje esperado |
|-------|---------------------|
| TAA | ~64% |
| TGA | ~30% |
| TAG | ~6% |

---

## 📝 Algoritmos BioPython

### Conteo de Codones

```python
from Bio import SeqIO

record = SeqIO.read("genome.gbff", "genbank")
sequence = record.seq

atg_count = sequence.count("ATG")
taa_count = sequence.count("TAA")
tag_count = sequence.count("TAG")
tga_count = sequence.count("TGA")
```

### Extracción de Genes

```python
from Bio.SeqUtils import gc_fraction

for feature in record.features:
    if feature.type == "gene":
        gene_seq = feature.extract(sequence)
        gc_content = gc_fraction(gene_seq) * 100
```

---

## 🤝 Contribución

1. Fork el repositorio
2. Crear rama feature (`git checkout -b feature/nueva-funcionalidad`)
3. Commit cambios (`git commit -m 'Add: nueva funcionalidad'`)
4. Push a la rama (`git push origin feature/nueva-funcionalidad`)
5. Abrir Pull Request

---

## 📜 Licencia

MIT License - ver [LICENSE](LICENSE)

---

## 👨‍💻 Autor

**Proyecto de Bioinformática**  
Universidad 2026

---

## 🙏 Agradecimientos

- [NCBI](https://www.ncbi.nlm.nih.gov/) por la API Datasets
- [BioPython](https://biopython.org/) por las herramientas genómicas
- [Anthropic](https://anthropic.com/) por Claude AI
