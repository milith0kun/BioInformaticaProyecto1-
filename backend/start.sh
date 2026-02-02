#!/bin/bash

# E. coli K-12 MG1655 Genomic Analysis System
# Script de inicio del backend

echo "🧬 Iniciando Backend - E. coli Analysis API"
echo "============================================"

# Cambiar al directorio del backend
cd "$(dirname "$0")"

# Verificar si existe el entorno virtual
if [ ! -d "venv" ]; then
    echo "📦 Creando entorno virtual..."
    python3 -m venv venv
fi

# Activar entorno virtual
echo "🔧 Activando entorno virtual..."
source venv/bin/activate

# Instalar dependencias
echo "📥 Instalando dependencias..."
pip install -r requirements.txt -q

# Crear directorio de cache si no existe
mkdir -p cache

# Iniciar servidor
echo ""
echo "🚀 Iniciando servidor en http://localhost:8000"
echo "📚 Documentación disponible en http://localhost:8000/docs"
echo ""
cd app
python -m uvicorn main:app --host 0.0.0.0 --port 8000 --reload
