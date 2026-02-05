# 🧬 Sistema de Análisis Genómico de E. coli

> Plataforma web para análisis bioinformático con búsqueda en NCBI, visualizaciones interactivas y validación con IA.

![Python](https://img.shields.io/badge/Python-3.12-blue?logo=python)
![FastAPI](https://img.shields.io/badge/FastAPI-0.128-green?logo=fastapi)
![React](https://img.shields.io/badge/React-18.2-61DAFB?logo=react)
![Claude AI](https://img.shields.io/badge/Claude_AI-3.5_Haiku-purple?logo=anthropic)

## ✨ Características

### 🔍 Búsqueda y Descarga de Genomas
- **NCBI Datasets API v2** - Búsqueda en tiempo real de genomas
- Descarga automática de archivos GenBank (.gbff)
- Soporte para múltiples genomas: E. coli, Bacillus, Saccharomyces, etc.
- Detección de duplicados (GCA/GCF del mismo genoma)

### 🧬 Análisis Bioinformático
- **Codones**: ATG (inicio), TAA/TAG/TGA (terminación)
- **Genes**: Extracción con BioPython, estadísticas de tamaño y GC%
- **Validación dinámica**: Compara contra promedio del grupo o rangos bacterianos
- **Comparación multi-genoma**: Análisis simultáneo de múltiples cepas

### 📊 Visualizaciones
- **Gráficos de barras**: Distribución de codones de terminación
- **Histogramas**: Tamaños de genes y CDS
- **Scatter plots**: GC% vs longitud génica
- **Tablas interactivas**: AG Grid con búsqueda, ordenamiento y paginación

### 🤖 Validación con IA (Claude 3.5 Haiku)
- Análisis contextual de resultados genómicos
- Validación científica automática
- Detección de discrepancias vs valores esperados
- Recomendaciones para análisis adicionales

### 📤 Exportación de Datos
- **JSON**: Análisis completo estructurado
- **CSV**: Por tipo (codones, genes, validación)
- **PDF**: Informe con gráficos embebidos

## 🚀 Instalación

### Backend
```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
cp .env.example .env  # Añadir CLAUDE_API_KEY
uvicorn app.main:app --reload
```

### Frontend
```bash
cd frontend
npm install
npm run dev
```

Abrir: **http://localhost:5173**

## 📱 Uso

1. **Buscar** - Escribe "Escherichia coli" o pega un accession (ej: \`GCF_000005845.2\`)
2. **Seleccionar** - Descarga y activa el genoma
3. **Analizar** - Ejecuta análisis completo
4. **Resultados** - Dashboard con métricas, visualizaciones y exportación

## 🔧 Stack

| Backend | Frontend | APIs Externas |
|---------|----------|---------------|
| FastAPI + Uvicorn | React + Vite | NCBI Datasets API v2 |
| BioPython (parsing) | Tailwind CSS | Claude AI (Anthropic) |
| Pandas (datos) | Recharts (gráficos) | |
| Pydantic (validación) | AG Grid (tablas) | |

## 🧬 Algoritmos y Métodos

### Análisis de Codones
```python
# BioPython - Conteo de secuencias
from Bio import SeqIO
record = SeqIO.read("genomic.gbff", "genbank")
atg_count = record.seq.count("ATG")  # Codón de inicio
taa_count = record.seq.count("TAA")  # Stop codon
```

### Extracción de Genes
```python
# Iteración sobre features del GenBank
for feature in record.features:
    if feature.type == "gene":
        gene_seq = feature.extract(record.seq)
        gc_content = gc_fraction(gene_seq) * 100
```

### Validación Dinámica
- **Multi-genoma**: Compara contra promedio µ ± 2σ del grupo
- **Single-genoma**: Valida contra rangos típicos bacterianos
- **Tolerancias**: GC% ±5%, longitud ±10%, genes ±15%

### IA - Claude 3.5 Haiku
```python
# Validación contextual con IA
response = client.messages.create(
    model="claude-3-5-haiku-20241022",
    messages=[{
        "role": "user",
        "content": f"Valida: {genome_data}"
    }]
)
```

**Output de IA incluye:**
- ✅ `is_valid`: true/false
- 📊 `confidence`: 0-100%
- 🔍 `key_findings`: Hallazgos clave
- 🧬 `scientific_context`: Contexto biológico
- ⚠️ `discrepancies`: Anomalías detectadas
- 💡 `recommendations`: Análisis sugeridos

## 📁 Estructura

```
├── backend/
│   ├── app/
│   │   ├── api/routes/     # Endpoints API
│   │   ├── core/           # Lógica de negocio
│   │   │   └── analyzers/  # Analizadores
│   │   └── models/         # Schemas Pydantic
│   └── genomes/            # Genomas descargados
│
└── frontend/
    └── src/
        ├── components/     # Componentes React
        └── services/       # API client
```

## 📊 Valores de Referencia (E. coli K-12 MG1655)

| Métrica | Valor |
|---------|-------|
| Tamaño genoma | 4.64 Mb |
| Genes totales | 4,651 |
| Contenido GC | 50.79% |

---

**Proyecto de Bioinformática** | Universidad 2026
