#!/bin/bash

# E. coli K-12 MG1655 Genomic Analysis System
# Script de inicio del frontend

echo "🎨 Iniciando Frontend - E. coli Analysis Dashboard"
echo "==================================================="

# Cambiar al directorio del frontend
cd "$(dirname "$0")"

# Verificar si node_modules existe
if [ ! -d "node_modules" ]; then
    echo "📦 Instalando dependencias..."
    npm install
fi

# Iniciar servidor de desarrollo
echo ""
echo "🚀 Iniciando servidor en http://localhost:5173"
echo ""
npm run dev
