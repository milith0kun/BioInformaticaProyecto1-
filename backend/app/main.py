"""
FastAPI main application for E. coli K-12 MG1655 Genomic Analysis System
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import os
import sys
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.api.routes import files, analysis, export, genome
from app.api.routes import ncbi_enhanced, chat
from app.core.file_detector import FileDetector


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown events"""
    # Startup: Initialize file detector and check for ZIP
    print("🧬 Iniciando Sistema de Análisis Genómico E. coli K-12...")
    
    # Get project root (two levels up from app/)
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    app.state.project_root = project_root
    
    # Initialize file detector
    app.state.file_detector = FileDetector(project_root)
    
    # Check and extract ZIP if exists
    zip_path = os.path.join(project_root, "ncbi_dataset.zip")
    if os.path.exists(zip_path):
        print(f"📦 Encontrado ncbi_dataset.zip, extrayendo...")
        app.state.file_detector.extract_zip(zip_path)
    
    # Also check for already extracted ncbi_dataset folder
    ncbi_folder = os.path.join(project_root, "ncbi_dataset")
    if os.path.exists(ncbi_folder):
        print(f"📁 Encontrada carpeta ncbi_dataset, escaneando archivos...")
        app.state.file_detector.scan_directory(ncbi_folder)

    # NEW: Check for genomes managed by NCBIDownloader (genomes/ folder)
    genomes_dir = os.path.join(project_root, "genomes")
    if os.path.exists(genomes_dir):
        print(f"📁 Encontrada estructura de genomas múltiples en: {genomes_dir}")
        try:
            from app.core.ncbi_downloader import NCBIDownloader
            api_key = os.environ.get("NCBI_API_KEY")
            downloader = NCBIDownloader(genomes_dir, api_key)
            app.state.ncbi_downloader = downloader
            
            # Auto-activate the first available genome found
            available = downloader.get_available_genomes_dir()
            if available:
                # Pick the first one (or ideally the last modified, but list is usually enough)
                latest = available[0]
                accession = latest['accession']
                print(f"🔄 Auto-activando genoma persistente: {accession}")
                
                # Load metadata
                downloader.get_genome_info(accession)
                
                # Update file detector to point to this genome's extracted files
                if latest.get('extracted_dir') and os.path.exists(latest['extracted_dir']):
                    app.state.file_detector.scan_directory(latest['extracted_dir'])
                    print(f"   ↳ Archivos escaneados en FileDetector: {len(app.state.file_detector.detected_files)}")
        except Exception as e:
            print(f"⚠️ Error al inicializar genomas: {e}")
    
    # Create cache directory
    cache_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache")
    os.makedirs(cache_dir, exist_ok=True)
    app.state.cache_dir = cache_dir
    
    print("✅ Sistema listo")
    
    yield
    
    # Shutdown
    print("👋 Cerrando sistema...")


app = FastAPI(
    title="E. coli K-12 MG1655 Genomic Analysis API",
    description="API para análisis genómico de E. coli K-12 MG1655",
    version="1.0.0",
    lifespan=lifespan
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for development
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(files.router, prefix="/api/files", tags=["Files"])
app.include_router(analysis.router, prefix="/api/analysis", tags=["Analysis"])
app.include_router(export.router, prefix="/api/export", tags=["Export"])
app.include_router(genome.router, prefix="/api/genome", tags=["Genome"])
app.include_router(ncbi_enhanced.router, prefix="/api/ncbi", tags=["NCBI Enhanced"])
app.include_router(chat.router, prefix="/api/chat", tags=["AI Chat"])


@app.get("/")
async def root():
    return {
        "message": "E. coli K-12 MG1655 Genomic Analysis API",
        "version": "1.0.0",
        "docs": "/docs"
    }


@app.get("/health")
async def health_check():
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
