# 🧬 Sistema de Análisis Genómico de E. coli

> Plataforma web para análisis bioinformático con búsqueda en NCBI, visualizaciones interactivas y validación con IA.

![Python](https://img.shields.io/badge/Python-3.12-blue?logo=python)
![FastAPI](https://img.shields.io/badge/FastAPI-0.128-green?logo=fastapi)
![React](https://img.shields.io/badge/React-18.2-61DAFB?logo=react)
![Claude AI](https://img.shields.io/badge/Claude_AI-3.5_Haiku-purple?logo=anthropic)

## ✨ Características

- 🔍 **Búsqueda NCBI** - Descarga genomas directamente desde NCBI Datasets API
- 🧬 **Análisis de codones** - ATG, TAA, TAG, TGA con estadísticas de densidad
- 📊 **Visualizaciones** - Gráficos interactivos con Recharts y tablas AG Grid
- 🤖 **Validación IA** - Claude AI para validación científica de resultados
- 📤 **Exportación** - JSON, CSV y PDF con gráficos embebidos
- ⚖️ **Comparación multi-genoma** - Compara múltiples cepas de E. coli

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

| Backend | Frontend |
|---------|----------|
| FastAPI + Uvicorn | React + Vite |
| BioPython | Tailwind CSS |
| Anthropic (Claude) | Recharts + AG Grid |

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
