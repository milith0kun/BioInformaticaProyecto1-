    # Plan de Mejoras - Sistema de Análisis Genómico

## Resumen
Mejorar el sistema usando la API de NCBI correctamente, agregar nuevas funcionalidades bioinformáticas, mejorar el chat IA, y redeseñar el dashboard.

## Cambios Backend

### 1. Nuevos endpoints NCBI (ncbi_service.py)
- GET /api/ncbi/proteins/{accession} - Obtener proteínas del genoma
- GET /api/ncbi/gene/{gene_id} - Detalles de un gen específico (posición, función, etc.)
- GET /api/ncbi/sequence/{accession} - Secuencia completa del genoma
- GET /api/ncbi/codon-usage/{accession} - Tabla de uso de codones completa (64 codones)
- GET /api/ncbi/gene-location/{accession}/{position} - Qué gen está en una posición dada

### 2. Chat IA Interactivo (chat_routes.py)  
- POST /api/chat/message - Chat interactivo con biólogo experto
- El chat tendrá contexto del genoma actual + datos de NCBI
- Citará referencias de GenBank/PubMed
- Usará APIs de NCBI para búsqueda de literatura

### 3. Análisis avanzado de codones (mejorar codon_analyzer.py)
- Tabla RSCU (Relative Synonymous Codon Usage)
- Codon Adaptation Index (CAI)
- Análisis por aminoácido
- Uso de codones por marco de lectura

### 4. Parser mejorado (genome_parser.py)
- Extraer proteínas (traducciones CDS)
- Extraer secuencia completa legible 
- Mapeo posición → gen

## Cambios Frontend

### 5. Nuevo componente: GenomeViewer (mapa circular + secuencia)
- Mapa circular del genoma con genes posicionados
- Vista de secuencia con navegación por posición
- Indicador de gen en posición actual

### 6. Nuevo componente: ProteinViewer
- Lista de proteínas del genoma
- Información de cada proteína

### 7. Nuevo componente: CodonUsageTable
- Tabla completa de 64 codones
- RSCU, CAI
- Gráficos de uso por aminoácido

### 8. Componente: InteractiveChat
- Chat tipo ChatGPT con el biólogo experto
- Historial de mensajes
- Referencias clickeables a NCBI/PubMed

### 9. Dashboard mejorado
- Rediseño con layout más informativo
- Gráficos adicionales (distribución GC por ventana, etc.)
- Más interactividad

## Prioridad de implementación
1. Backend: endpoints NCBI nuevos + análisis codones mejorado
2. Backend: Chat IA interactivo  
3. Frontend: Nuevos componentes (GenomeViewer, ProteinViewer, CodonUsage)
4. Frontend: Chat interactivo
5. Frontend: Rediseño dashboard
