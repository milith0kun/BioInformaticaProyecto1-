#!/bin/bash

# E. coli K-12 MG1655 Genomic Analysis System
# Script principal de inicio

echo "========================================================"
echo "  🧬 E. coli K-12 MG1655 - Sistema de Análisis Genómico"
echo "========================================================"
echo ""

PROJECT_DIR="$(dirname "$0")"
cd "$PROJECT_DIR"

# Función para iniciar backend
start_backend() {
    echo "🔧 Iniciando Backend..."
    cd backend
    
    # Crear y activar entorno virtual si no existe
    if [ ! -d "venv" ]; then
        python3 -m venv venv
    fi
    source venv/bin/activate
    pip install -r requirements.txt -q
    mkdir -p cache
    
    cd app
    python -m uvicorn main:app --host 0.0.0.0 --port 8000 --reload &
    BACKEND_PID=$!
    cd ../..
    echo "✅ Backend iniciado (PID: $BACKEND_PID)"
}

# Función para iniciar frontend
start_frontend() {
    echo "🎨 Iniciando Frontend..."
    cd frontend
    
    if [ ! -d "node_modules" ]; then
        npm install
    fi
    
    npm run dev &
    FRONTEND_PID=$!
    cd ..
    echo "✅ Frontend iniciado (PID: $FRONTEND_PID)"
}

# Trap para limpiar procesos al salir
cleanup() {
    echo ""
    echo "🛑 Deteniendo servicios..."
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    echo "👋 ¡Hasta luego!"
    exit 0
}
trap cleanup SIGINT SIGTERM

# Iniciar servicios
start_backend
sleep 2
start_frontend

echo ""
echo "========================================================"
echo "  🌐 Backend:  http://localhost:8000"
echo "  🌐 Frontend: http://localhost:5173"
echo "  📚 API Docs: http://localhost:8000/docs"
echo "========================================================"
echo ""
echo "Presiona Ctrl+C para detener todos los servicios"
echo ""

# Mantener script corriendo
wait
