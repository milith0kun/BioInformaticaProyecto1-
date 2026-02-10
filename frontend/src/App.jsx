import { useState, useEffect } from 'react'
import toast, { Toaster } from 'react-hot-toast'
import FileManager from './components/FileManager'
import AnalysisDashboard from './components/AnalysisDashboard'
import CodonVisualization from './components/CodonVisualization'
import GeneStatistics from './components/GeneStatistics'
import DataExport from './components/DataExport'
import AIValidation from './components/AIValidation'
import GenomeSelector from './components/GenomeSelector'
import GenomeComparison from './components/GenomeComparison'
import GeneFilter from './components/GeneFilter'
import MultiGenomeAnalyzer from './components/MultiGenomeAnalyzer'
import GenomeMultiSelector from './components/GenomeMultiSelector'
import ComparisonResults from './components/ComparisonResults'
import InteractiveChat from './components/InteractiveChat'
import GenomeViewer from './components/GenomeViewer'
import ProteinViewer from './components/ProteinViewer'
import CodonUsageTable from './components/CodonUsageTable'
import CentralDogma from './components/CentralDogma'
import GCWindowViewer from './components/GCWindowViewer'
import RNAAnalysis from './components/RNAAnalysis'
import BLASTSearch from './components/BLASTSearch'
import PhylogeneticTree from './components/PhylogeneticTree'
import FunctionalCategories from './components/FunctionalCategories'
import CAIAnalysis from './components/CAIAnalysis'
import { api } from './services/api'

function App() {
  // Estado principal del workflow
  const [currentStep, setCurrentStep] = useState(1)
  const [files, setFiles] = useState([])
  const [analysisData, setAnalysisData] = useState(null)
  const [aiValidation, setAiValidation] = useState(null)
  const [isLoading, setIsLoading] = useState(false)
  const [isValidatingAI, setIsValidatingAI] = useState(false)
  const [currentGenome, setCurrentGenome] = useState(null)
  const [downloadedGenomes, setDownloadedGenomes] = useState([])
  const [selectedGenomes, setSelectedGenomes] = useState([]) // Múltiples genomas seleccionados
  const [comparisonResult, setComparisonResult] = useState(null) // Resultado de comparación

  // Vista activa para resultados (después del paso 3)
  const [activeView, setActiveView] = useState('dashboard')

  // Estado para sidebar móvil y preferencias
  const [sidebarOpen, setSidebarOpen] = useState(false)
  const [sidebarSearch, setSidebarSearch] = useState('')

  // Cargar preferencias de sidebar desde localStorage
  useEffect(() => {
    const saved = localStorage.getItem('sidebarPreferences')
    if (saved) {
      try {
        const prefs = JSON.parse(saved)
        if (prefs.defaultOpen !== undefined && window.innerWidth >= 1024) {
          setSidebarOpen(prefs.defaultOpen)
        }
      } catch (e) {
        console.error('Error loading sidebar preferences:', e)
      }
    }
  }, [])

  // Guardar preferencias cuando cambia el estado del sidebar
  useEffect(() => {
    if (window.innerWidth >= 1024) {
      localStorage.setItem('sidebarPreferences', JSON.stringify({ defaultOpen: sidebarOpen }))
    }
  }, [sidebarOpen])

  useEffect(() => {
    loadFiles()
    loadDownloadedGenomes()
  }, [])

  const loadFiles = async () => {
    try {
      const response = await api.detectFiles()
      setFiles(response.files)
    } catch (error) {
      console.error('Error loading files:', error)
    }
  }

  const loadDownloadedGenomes = async () => {
    try {
      const result = await api.getDownloadedGenomes()
      setDownloadedGenomes(result.genomes || [])
    } catch (error) {
      console.error('Error:', error)
    }
  }

  const handleGenomeDownloaded = () => {
    loadDownloadedGenomes()
    loadFiles()
  }

  const handleGenomeActivated = (result) => {
    setCurrentGenome(result.genome_info)
    setSelectedGenomes([result.genome_info.accession])
    loadFiles()
    setAnalysisData(null)
    setAiValidation(null)
    setComparisonResult(null)
    setCurrentStep(3) // Avanzar al paso de análisis
    toast.success(`Genoma ${result.genome_info.accession} seleccionado`)
  }

  // Manejar selección múltiple
  const handleMultipleSelection = (accessions) => {
    setSelectedGenomes(accessions)
    // Si hay al menos 1, tomar el primero como "actual" para compatibilidad
    if (accessions.length > 0) {
      const genome = downloadedGenomes.find(g => g.accession === accessions[0])
      if (genome) {
        setCurrentGenome(genome)
      }
    }
  }

  // Ejecutar análisis (1 o más genomas)
  const runAnalysis = async () => {
    if (selectedGenomes.length === 0) {
      toast.error('Seleccione al menos un genoma')
      return
    }

    setIsLoading(true)
    const failedGenomes = []
    const successfulGenomes = []

    try {
      // Flujo unificado: descargar genomas y comparar (funciona para 1 o más genomas)
      toast.loading(`Descargando ${selectedGenomes.length} genoma${selectedGenomes.length > 1 ? 's' : ''}...`, { id: 'analysis' })

      // Descargar todos los genomas seleccionados
      for (let i = 0; i < selectedGenomes.length; i++) {
        const accession = selectedGenomes[i]

        // Verificar si ya está descargado
        const isAlreadyDownloaded = downloadedGenomes.some(g => g.accession === accession)

        if (isAlreadyDownloaded) {
          toast.loading(`Genoma ${i + 1}/${selectedGenomes.length} ya descargado ✓`, { id: 'analysis', duration: 1000 })
          successfulGenomes.push(accession)
          continue // Saltar al siguiente
        }

        toast.loading(`Descargando genoma ${i + 1}/${selectedGenomes.length}: ${accession}...`, { id: 'analysis' })

        try {
          // Iniciar descarga
          await api.downloadGenome({
            accession: accession,
            include_gbff: true,
            include_gff: true,
            include_fasta: true
          })

          // Esperar a que la descarga termine
          let downloadComplete = false
          let attempts = 0
          const maxAttempts = 60 // 60 segundos máximo por genoma

          while (!downloadComplete && attempts < maxAttempts) {
            await new Promise(resolve => setTimeout(resolve, 1000)) // Esperar 1 segundo
            const status = await api.getGenomeDownloadStatus(accession)
            if (status.status === 'completed') {
              downloadComplete = true
              successfulGenomes.push(accession)
              toast.loading(`Genoma ${i + 1}/${selectedGenomes.length} descargado ✓`, { id: 'analysis' })
            } else if (status.status === 'error') {
              throw new Error(status.message || 'Error en descarga')
            }
            attempts++
          }

          if (!downloadComplete) {
            throw new Error(`Timeout descargando ${accession}`)
          }
        } catch (downloadError) {
          console.warn(`Error descargando ${accession}:`, downloadError)
          failedGenomes.push({ accession, error: downloadError.message })
          toast.error(`${accession} falló: ${downloadError.message.substring(0, 50)}...`, { duration: 4000 })
        }
      }

      toast.loading('Ejecutando análisis comparativo...', { id: 'analysis' })

      // Esperar 2 segundos adicionales para asegurar que los archivos estén listos
      await new Promise(resolve => setTimeout(resolve, 2000))

      // Comparar solo los genomas exitosos
      const comparison = await api.compareGenomes(successfulGenomes)
      setComparisonResult(comparison)

      // Siempre activar el último genoma exitoso para análisis detallado (Dashboard, Codones, Genes, etc.)
      toast.loading('Obteniendo análisis detallado con codones...', { id: 'analysis' })
      try {
        const lastGenome = successfulGenomes[successfulGenomes.length - 1]
        await api.activateGenome(lastGenome)
        await loadFiles()
        const detailedAnalysis = await api.runCompleteAnalysis()
        setAnalysisData(detailedAnalysis)
      } catch (detailError) {
        console.warn('No se pudo obtener análisis detallado:', detailError)
        setAnalysisData(null)
      }

      setCurrentStep(3)
      setActiveView('comparison')

      // Actualizar lista de genomas descargados
      await loadDownloadedGenomes()

      // Mensaje final con resumen
      if (failedGenomes.length > 0) {
        toast.success(
          `✓ ${successfulGenomes.length} genoma${successfulGenomes.length > 1 ? 's' : ''} • ✗ ${failedGenomes.length} falló${failedGenomes.length > 1 ? 'fallaron' : ''}`,
          { id: 'analysis', duration: 5000 }
        )
      } else {
        toast.success(`${successfulGenomes.length} genoma${successfulGenomes.length > 1 ? 's' : ''} analizado${successfulGenomes.length > 1 ? 's' : ''}`, { id: 'analysis' })
      }
    } catch (error) {
      console.error('Analysis error:', error)
      toast.error('Error en el análisis: ' + (error.response?.data?.detail || error.message), { id: 'analysis' })
    } finally {
      setIsLoading(false)
    }
  }

  const runAIValidation = async () => {
    if (!analysisData) {
      toast.error('Ejecute el análisis primero')
      return
    }

    setIsValidatingAI(true)
    try {
      toast.loading('Validando con IA...', { id: 'ai-validation' })
      const result = await api.validateWithAI()
      setAiValidation(result)
      setActiveView('ai')
      toast.success('Validación completada', { id: 'ai-validation' })
    } catch (error) {
      console.error('AI Validation error:', error)
      toast.error('Error en validación IA', { id: 'ai-validation' })
    } finally {
      setIsValidatingAI(false)
    }
  }

  // Definir los pasos del workflow
  const steps = [
    { id: 1, name: 'Seleccionar', description: 'Elegir genomas' },
    { id: 2, name: 'Analizar', description: 'Procesar datos' },
    { id: 3, name: 'Resultados', description: 'Ver análisis' },
  ]

  // Organización mejorada de vistas por categorías
  const navigationSections = [
    {
      title: 'Inicio',
      items: [
        { id: 'dashboard', name: 'Dashboard', description: 'Vista general del análisis' },
      ]
    },
    {
      title: 'Análisis Comparativo',
      items: [
        { id: 'comparison', name: 'Comparación de Genomas', description: 'Comparar múltiples genomas' },
      ]
    },
    {
      title: 'Análisis Genómico',
      items: [
        { id: 'genome-map', name: 'Mapa del Genoma', description: 'Visualización del genoma completo' },
        { id: 'codons', name: 'Análisis de Codones', description: 'Frecuencia de codones' },
        { id: 'codon-usage', name: 'Tabla de Uso de Codones', description: 'Tabla detallada de uso' },
        { id: 'genes', name: 'Estadísticas de Genes', description: 'Análisis de genes' },
        { id: 'proteins', name: 'Análisis de Proteínas', description: 'Visualización de proteínas' },
      ]
    },
    {
      title: 'Análisis de RNA',
      items: [
        { id: 'rna', name: 'tRNA y rRNA', description: 'Análisis de RNA' },
      ]
    },
    {
      title: 'Herramientas',
      items: [
        { id: 'dogma', name: 'Dogma Central', description: 'DNA → RNA → Proteína' },
        { id: 'gc-window', name: 'Ventana GC', description: 'Contenido GC por ventana' },
        { id: 'blast', name: 'BLAST Search', description: 'Búsqueda de similitud' },
        { id: 'filter', name: 'Filtro de Genes', description: 'Filtrado avanzado' },
      ]
    },
    {
      title: 'Análisis Avanzado',
      items: [
        { id: 'phylo', name: 'Árbol Filogenético', description: 'Análisis evolutivo' },
        { id: 'cog', name: 'Categorías COG', description: 'Clasificación funcional' },
        { id: 'cai', name: 'Índice CAI', description: 'Adaptación de codones' },
      ]
    },
    {
      title: 'IA y Utilidades',
      items: [
        { id: 'chat', name: 'Chat IA', description: 'Asistente inteligente' },
        { id: 'ai', name: 'Validación IA', description: 'Validación con IA' },
        { id: 'files', name: 'Archivos', description: 'Gestión de archivos' },
        { id: 'export', name: 'Exportar Datos', description: 'Exportar resultados' },
      ]
    }
  ]

  return (
    <div className="h-screen bg-[#f8fafc] flex overflow-hidden font-sans">
      <Toaster position="top-right" />

      {/* Sidebar Navigation */}
      <aside className={`fixed lg:static inset-y-0 left-0 z-50 w-80 bg-white border-r border-slate-200 transform transition-transform duration-500 ease-in-out ${sidebarOpen ? 'translate-x-0' : '-translate-x-full lg:translate-x-0'} flex flex-col h-screen shadow-2xl lg:shadow-none`}>

        {/* Sidebar Header */}
        <div className="p-8 border-b border-slate-100 bg-white">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-4">
              <div className="w-12 h-12 bg-blue-600 rounded-2xl flex items-center justify-center shadow-lg shadow-blue-200 transition-transform hover:rotate-12">
                <svg className="w-7 h-7 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
              </div>
              <div>
                <h1 className="text-sm font-black text-slate-900 uppercase tracking-tighter">Genómica <span className="text-blue-600">Lab</span></h1>
                <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest">Análisis Molecular</p>
              </div>
            </div>
            <button
              onClick={() => setSidebarOpen(false)}
              className="lg:hidden p-2 rounded-xl hover:bg-slate-100 transition-colors"
            >
              <svg className="w-5 h-5 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>
        </div>

        {/* Workflow Progress */}
        {currentStep < 3 && (
          <div className="p-6 bg-slate-50/50 border-b border-slate-100 flex-shrink-0">
            <div className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] mb-4 px-2">Pipeline Progress</div>
            <div className="space-y-3">
              {steps.map((step) => (
                <button
                  key={step.id}
                  onClick={() => {
                    if (step.id <= currentStep || (step.id === 2 && downloadedGenomes.length > 0)) {
                      setCurrentStep(step.id)
                    }
                  }}
                  disabled={step.id > currentStep && !(step.id === 2 && downloadedGenomes.length > 0)}
                  className={`w-full flex items-center gap-4 px-4 py-3 rounded-2xl text-sm transition-all duration-300 ${
                    step.id === currentStep
                      ? 'bg-blue-600 text-white shadow-xl shadow-blue-200'
                      : step.id < currentStep
                      ? 'bg-blue-50 text-blue-700 hover:bg-blue-100'
                      : 'bg-slate-100/50 text-slate-400 cursor-not-allowed opacity-50'
                  }`}
                >
                  <div className={`w-7 h-7 rounded-xl flex items-center justify-center text-xs font-black ${
                    step.id === currentStep
                      ? 'bg-white/20'
                      : step.id < currentStep
                      ? 'bg-blue-100'
                      : 'bg-slate-200'
                  }`}>
                    {step.id < currentStep ? (
                      <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                        <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                      </svg>
                    ) : step.id}
                  </div>
                  <div className="text-left flex-1">
                    <div className="font-black uppercase tracking-tighter">{step.name}</div>
                    <div className="text-[10px] opacity-70 font-bold">{step.description}</div>
                  </div>
                </button>
              ))}
            </div>
          </div>
        )}

        {/* Navigation Sections */}
        {currentStep === 3 && (
          <>
            <div className="p-6 border-b border-slate-50 flex-shrink-0">
              <div className="relative group">
                <input
                  type="text"
                  value={sidebarSearch}
                  onChange={(e) => setSidebarSearch(e.target.value)}
                  placeholder="Filtrar herramientas..."
                  className="w-full pl-10 pr-4 py-3 text-xs font-bold border border-slate-100 bg-slate-50 rounded-2xl focus:outline-none focus:ring-4 focus:ring-blue-500/5 focus:border-blue-500/30 transition-all placeholder-slate-400"
                />
                <svg className="w-4 h-4 text-slate-400 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-blue-500 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                </svg>
              </div>
            </div>

            <nav className="flex-1 overflow-y-auto p-6 sidebar-scroll min-h-0 space-y-8">
              {navigationSections.map((section, idx) => {
                const filteredItems = section.items.filter(item =>
                  item.name.toLowerCase().includes(sidebarSearch.toLowerCase()) ||
                  item.description.toLowerCase().includes(sidebarSearch.toLowerCase())
                )

                if (filteredItems.length === 0 && sidebarSearch) return null

                return (
                  <div key={idx} className="space-y-2">
                    <h2 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-3">
                      {section.title}
                    </h2>
                    <div className="space-y-1">
                      {filteredItems.map((item) => (
                        <button
                          key={item.id}
                          onClick={() => {
                            setActiveView(item.id)
                            setSidebarOpen(false)
                            setSidebarSearch('')
                          }}
                          className={`w-full flex flex-col items-start px-4 py-3 rounded-2xl text-xs transition-all duration-300 ${
                            activeView === item.id
                              ? 'bg-blue-50 text-blue-700 shadow-sm border border-blue-100/50'
                              : 'text-slate-600 hover:bg-slate-50 hover:text-slate-900'
                          }`}
                        >
                          <span className="font-black uppercase tracking-tight">{item.name}</span>
                          <span className="text-[9px] text-slate-400 font-bold mt-0.5 uppercase tracking-wider">{item.description}</span>
                        </button>
                      ))}
                    </div>
                  </div>
                )
              })}
            </nav>
          </>
        )}

        {/* Active Genome HUD */}
        {currentGenome && currentStep === 3 && (
          <div className="p-6 border-t border-slate-50 bg-slate-50/30 flex-shrink-0">
            <div className="bg-white rounded-2xl p-4 border border-slate-100 shadow-sm">
              <div className="flex items-center gap-3 mb-2">
                <div className="w-2 h-2 bg-emerald-500 rounded-full animate-pulse shadow-[0_0_8px_rgba(16,185,129,0.5)]"></div>
                <span className="font-mono text-[10px] font-black text-blue-600 tracking-widest">{currentGenome.accession}</span>
              </div>
              <p className="text-[10px] font-bold text-slate-500 uppercase truncate tracking-tight">{currentGenome.organism_name}</p>
            </div>
          </div>
        )}
      </aside>

      {/* Main Content Area */}
      <div className="flex-1 flex flex-col min-w-0 h-screen overflow-hidden">

        {/* Top Bar */}
        <header className="bg-white/80 backdrop-blur-2xl border-b border-slate-100 px-8 py-6 flex-shrink-0 z-30">
          <div className="flex flex-col gap-4">
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-6">
                <button
                  onClick={() => setSidebarOpen(true)}
                  className="lg:hidden p-3 bg-slate-50 rounded-2xl hover:bg-slate-100 transition-colors border border-slate-100"
                >
                  <svg className="w-6 h-6 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" />
                  </svg>
                </button>

                <div>
                  <h2 className="text-2xl font-black text-slate-900 tracking-tighter uppercase italic">
                    {currentStep === 1 && 'Seleccionar Genomas'}
                    {currentStep === 2 && 'Ejecutar Análisis'}
                    {currentStep === 3 && navigationSections.flatMap(s => s.items).find(item => item.id === activeView)?.name}
                  </h2>
                  {currentStep === 3 && (
                    <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] mt-1">
                      {navigationSections.flatMap(s => s.items).find(item => item.id === activeView)?.description}
                    </p>
                  )}
                </div>
              </div>

              <div className="flex items-center gap-3">
                {currentStep === 3 && (
                  <button
                    onClick={() => setCurrentStep(2)}
                    className="flex items-center gap-3 px-5 py-2.5 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-xl shadow-slate-200 transition-all hover:bg-blue-600 hover:-translate-y-1 active:scale-95"
                  >
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M12 4v16m8-8H4" />
                    </svg>
                    Nueva Sesión
                  </button>
                )}
              </div>
            </div>
          </div>
        </header>

        {/* Main Workspace */}
        <main className="flex-1 overflow-y-auto p-8 lg:p-12 min-h-0 custom-scrollbar">
          <div className="max-w-[1600px] mx-auto">
            {/* Step Content */}
            {currentStep === 1 && (
              <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-700">
                <div className="max-w-3xl">
                  <h3 className="text-3xl font-black text-slate-900 tracking-tight mb-4 uppercase italic">Configuración de Experimento</h3>
                  <p className="text-slate-500 font-medium leading-relaxed">
                    Selecciona los genomas de interés. El sistema sincronizará los datos directamente desde las bases de datos de NCBI.
                  </p>
                </div>

                <GenomeMultiSelector
                  downloadedGenomes={downloadedGenomes}
                  selectedGenomes={selectedGenomes}
                  onSelectionChange={handleMultipleSelection}
                  onRefresh={loadDownloadedGenomes}
                />

                <div className="flex justify-end p-8 bg-white rounded-[3rem] border border-slate-100 shadow-xl relative overflow-hidden group">
                  <div className="absolute inset-0 bg-blue-600 translate-y-full group-hover:translate-y-0 transition-transform duration-700"></div>
                  <button
                    onClick={() => {
                      if (selectedGenomes.length === 0) {
                        toast.error('Selecciona al menos un genoma')
                        return
                      }
                      setCurrentStep(2)
                    }}
                    disabled={selectedGenomes.length === 0}
                    className="relative flex items-center gap-4 text-slate-900 group-hover:text-white transition-colors font-black uppercase tracking-widest text-sm"
                  >
                    Siguiente Fase: {selectedGenomes.length} Genomas
                    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M17 8l4 4m0 0l-4 4m4-4H3" />
                    </svg>
                  </button>
                </div>
              </div>
            )}

            {currentStep === 2 && (
              <div className="max-w-4xl mx-auto space-y-10 animate-in fade-in zoom-in-95 duration-700">
                <div className="text-center space-y-4">
                  <div className="inline-flex p-4 bg-blue-50 rounded-[2rem] mb-4">
                    <svg className="w-10 h-10 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                    </svg>
                  </div>
                  <h3 className="text-4xl font-black text-slate-900 uppercase italic tracking-tighter">Motor de Análisis</h3>
                  <p className="text-slate-500 font-bold uppercase text-[10px] tracking-[0.3em]">Preparando entorno de cómputo</p>
                </div>

                <div className="bg-white rounded-[3rem] border border-slate-100 p-12 shadow-2xl space-y-10">
                  <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
                    {selectedGenomes.map((accession, i) => (
                      <div key={accession} className="flex items-center gap-4 p-5 bg-slate-50 rounded-3xl border border-slate-100 transition-all hover:bg-white hover:shadow-xl group/item">
                        <span className="w-10 h-10 bg-white shadow-sm rounded-2xl flex items-center justify-center text-xs font-black text-blue-600 border border-slate-100 group-hover/item:bg-blue-600 group-hover/item:text-white transition-all">
                          {i + 1}
                        </span>
                        <div className="flex-1 min-w-0 font-mono text-xs font-black text-slate-700 tracking-tighter uppercase">{accession}</div>
                      </div>
                    ))}
                  </div>

                  <button
                    onClick={runAnalysis}
                    disabled={isLoading}
                    className="w-full py-6 bg-slate-900 text-white font-black uppercase tracking-[0.2em] rounded-[2rem] transition-all hover:bg-blue-600 hover:-translate-y-2 shadow-2xl shadow-blue-200 disabled:opacity-50 flex items-center justify-center gap-6 group"
                  >
                    {isLoading ? (
                      <>
                        <div className="w-6 h-6 border-4 border-white/20 border-t-white rounded-full animate-spin"></div>
                        <span className="tracking-widest">Sincronizando...</span>
                      </>
                    ) : (
                      <>
                        <span>Iniciar Procesamiento Genómico</span>
                        <svg className="w-6 h-6 group-hover:translate-x-2 transition-transform" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M13 10V3L4 14h7v7l9-11h-7z" />
                        </svg>
                      </>
                    )}
                  </button>
                </div>
              </div>
            )}

            {currentStep === 3 && (
              <div className="space-y-10 animate-in fade-in duration-1000">
                {selectedGenomes.length > 1 && (
                  <div className="flex items-center justify-between p-6 bg-white rounded-3xl border border-slate-100 shadow-sm relative overflow-hidden">
                    <div className="absolute left-0 top-0 bottom-0 w-1.5 bg-blue-600"></div>
                    <div className="pl-4">
                      <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-1">Active Comparison</p>
                      <p className="text-sm font-black text-slate-900 uppercase italic tracking-tighter">{selectedGenomes.length} Genomas en entorno de trabajo</p>
                    </div>
                    <button onClick={() => setCurrentStep(1)} className="px-6 py-2.5 bg-slate-50 hover:bg-slate-100 text-[10px] font-black uppercase tracking-widest text-slate-600 rounded-2xl transition-all border border-slate-100">
                      Reconfigurar
                    </button>
                  </div>
                )}

                <div className="bg-white rounded-[3rem] border border-slate-100 p-10 shadow-sm min-h-[800px]">
                  {activeView === 'dashboard' && <AnalysisDashboard analysisData={analysisData} isLoading={isLoading} status={analysisData ? 'completed' : 'idle'} />}
                  {activeView === 'codons' && <CodonVisualization codonData={analysisData?.codons} />}
                  {activeView === 'comparison' && <ComparisonResults comparisonResult={comparisonResult} selectedGenomes={selectedGenomes} />}
                  {activeView === 'genes' && <GeneStatistics geneData={analysisData?.genes} />}
                  {activeView === 'filter' && <GeneFilter hasAnalysis={!!analysisData} />}
                  {activeView === 'files' && <FileManager files={files} onRefresh={loadFiles} selectedGenomes={selectedGenomes} />}
                  {activeView === 'genome-map' && <GenomeViewer />}
                  {activeView === 'dogma' && <CentralDogma />}
                  {activeView === 'proteins' && <ProteinViewer />}
                  {activeView === 'codon-usage' && <CodonUsageTable />}
                  {activeView === 'gc-window' && <GCWindowViewer />}
                  {activeView === 'rna' && <RNAAnalysis />}
                  {activeView === 'blast' && <BLASTSearch />}
                  {activeView === 'phylo' && <PhylogeneticTree />}
                  {activeView === 'cog' && <FunctionalCategories />}
                  {activeView === 'cai' && <CAIAnalysis />}
                  {activeView === 'chat' && <InteractiveChat hasAnalysis={!!analysisData} currentGenome={currentGenome} />}
                  {activeView === 'ai' && <AIValidation validationData={aiValidation} isValidating={isValidatingAI} onValidate={runAIValidation} hasAnalysis={!!analysisData} comparisonData={comparisonResult} selectedGenomes={selectedGenomes} currentGenome={currentGenome} />}
                  {activeView === 'export' && <DataExport hasData={analysisData !== null} comparisonData={comparisonResult} currentGenome={currentGenome} selectedGenomes={selectedGenomes} />}
                </div>
              </div>
            )}
          </div>
        </main>

        {/* Bottom Bar */}
        <footer className="bg-white border-t border-slate-100 py-4 px-10 flex-shrink-0 flex items-center justify-between">
          <div className="flex items-center gap-4">
            <div className="w-2 h-2 bg-blue-600 rounded-full animate-pulse shadow-[0_0_8px_rgba(37,99,235,0.5)]"></div>
            <span className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em]">LABORATORY SYSTEM ONLINE</span>
          </div>
          <div className="text-[9px] font-bold text-slate-300 uppercase tracking-widest">
            VITE + FASTAPI CORE
          </div>
        </footer>
      </div>
    </div>
  )
}

export default App
