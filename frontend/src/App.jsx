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
    <div className="h-screen bg-slate-50 flex overflow-hidden">
      <Toaster position="top-right" />

      {/* Sidebar Navigation */}
      <aside className={`fixed lg:static inset-y-0 left-0 z-50 w-72 bg-white border-r border-slate-200 transform transition-transform duration-300 ease-in-out ${sidebarOpen ? 'translate-x-0' : '-translate-x-full lg:translate-x-0'} flex flex-col h-screen`}>

        {/* Sidebar Header */}
        <div className="p-6 border-b border-slate-200">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 bg-gradient-to-br from-teal-600 to-emerald-600 rounded-lg flex items-center justify-center">
                <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
              </div>
              <div>
                <h1 className="text-sm font-bold text-slate-800">Análisis Genómico</h1>
                <p className="text-xs text-slate-500">NCBI Datasets v2</p>
              </div>
            </div>
            <button
              onClick={() => setSidebarOpen(false)}
              className="lg:hidden p-1 rounded hover:bg-slate-100"
            >
              <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>
        </div>

        {/* Workflow Progress - Solo visible en pasos 1 y 2 */}
        {currentStep < 3 && (
          <div className="p-4 bg-slate-50 border-b border-slate-200 flex-shrink-0">
            <div className="text-xs font-medium text-slate-600 mb-2">Progreso del Flujo</div>
            <div className="space-y-2">
              {steps.map((step) => (
                <button
                  key={step.id}
                  onClick={() => {
                    if (step.id <= currentStep || (step.id === 2 && downloadedGenomes.length > 0)) {
                      setCurrentStep(step.id)
                    }
                  }}
                  disabled={step.id > currentStep && !(step.id === 2 && downloadedGenomes.length > 0)}
                  className={`w-full flex items-center gap-2 px-3 py-2 rounded-lg text-sm transition-all ${
                    step.id === currentStep
                      ? 'bg-teal-600 text-white'
                      : step.id < currentStep
                      ? 'bg-teal-50 text-teal-700 hover:bg-teal-100'
                      : 'bg-slate-100 text-slate-400 cursor-not-allowed'
                  }`}
                >
                  <div className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-medium ${
                    step.id === currentStep
                      ? 'bg-white/20'
                      : step.id < currentStep
                      ? 'bg-teal-100'
                      : 'bg-slate-200'
                  }`}>
                    {step.id < currentStep ? (
                      <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
                        <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                      </svg>
                    ) : step.id}
                  </div>
                  <div className="text-left flex-1">
                    <div className="font-medium">{step.name}</div>
                    <div className="text-xs opacity-80">{step.description}</div>
                  </div>
                </button>
              ))}
            </div>
          </div>
        )}

        {/* Navigation Sections - Solo visible en paso 3 */}
        {currentStep === 3 && (
          <>
            {/* Búsqueda rápida en sidebar */}
            <div className="p-4 border-b border-slate-200 flex-shrink-0">
              <div className="relative">
                <input
                  type="text"
                  value={sidebarSearch}
                  onChange={(e) => setSidebarSearch(e.target.value)}
                  placeholder="Buscar herramienta..."
                  className="w-full pl-9 pr-4 py-2 text-sm border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-teal-400 focus:border-transparent"
                />
                <svg className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                </svg>
                {sidebarSearch && (
                  <button
                    onClick={() => setSidebarSearch('')}
                    className="absolute right-3 top-1/2 -translate-y-1/2 text-slate-400 hover:text-slate-600"
                  >
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                    </svg>
                  </button>
                )}
              </div>
            </div>

            <nav className="flex-1 overflow-y-auto p-4 sidebar-scroll min-h-0">
              {navigationSections.map((section, idx) => {
                // Filtrar items basado en búsqueda
                const filteredItems = section.items.filter(item =>
                  item.name.toLowerCase().includes(sidebarSearch.toLowerCase()) ||
                  item.description.toLowerCase().includes(sidebarSearch.toLowerCase())
                )

                // Si no hay items que coincidan, no mostrar la sección
                if (filteredItems.length === 0 && sidebarSearch) return null

                return (
                  <div key={idx} className="mb-6">
                    <h2 className="text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2 px-3">
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
                          className={`w-full flex flex-col items-start px-3 py-2 rounded-lg text-sm transition-all sidebar-nav-item ${
                            activeView === item.id
                              ? 'bg-teal-50 text-teal-700 font-medium active'
                              : 'text-slate-700 hover:bg-slate-100'
                          }`}
                        >
                          <span>{item.name}</span>
                          <span className="text-xs text-slate-500 mt-0.5">{item.description}</span>
                        </button>
                      ))}
                    </div>
                  </div>
                )
              })}

              {/* Mensaje si no hay resultados */}
              {sidebarSearch && navigationSections.every(section =>
                section.items.every(item =>
                  !item.name.toLowerCase().includes(sidebarSearch.toLowerCase()) &&
                  !item.description.toLowerCase().includes(sidebarSearch.toLowerCase())
                )
              ) && (
                <div className="text-center py-8">
                  <svg className="w-12 h-12 text-slate-300 mx-auto mb-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                  </svg>
                  <p className="text-sm text-slate-500">No se encontraron herramientas</p>
                  <button
                    onClick={() => setSidebarSearch('')}
                    className="text-xs text-teal-600 hover:text-teal-700 mt-2"
                  >
                    Limpiar búsqueda
                  </button>
                </div>
              )}
            </nav>
          </>
        )}

        {/* Genome Info - Solo visible cuando hay genoma y estamos en paso 3 */}
        {currentGenome && currentStep === 3 && (
          <div className="p-4 border-t border-slate-200 bg-slate-50 flex-shrink-0">
            <div className="text-xs font-medium text-slate-600 mb-2">Genoma Activo</div>
            <div className="bg-white rounded-lg p-3 border border-slate-200">
              <div className="flex items-center gap-2 mb-1">
                <span className="w-2 h-2 bg-emerald-500 rounded-full"></span>
                <span className="font-mono text-xs font-semibold text-teal-700">{currentGenome.accession}</span>
              </div>
              <p className="text-xs text-slate-600 truncate">{currentGenome.organism_name}</p>
              {selectedGenomes.length > 1 && (
                <p className="text-xs text-slate-500 mt-1">
                  Comparando {selectedGenomes.length} genomas
                </p>
              )}
            </div>
          </div>
        )}
      </aside>

      {/* Overlay para cerrar sidebar en móvil */}
      {sidebarOpen && (
        <div
          onClick={() => setSidebarOpen(false)}
          className="fixed inset-0 bg-black/50 z-40 lg:hidden"
        />
      )}

      {/* Main Content Area */}
      <div className="flex-1 flex flex-col min-w-0 h-screen overflow-hidden">

        {/* Top Bar */}
        <header className="bg-white border-b border-slate-200 px-4 py-3 lg:px-6 lg:py-4 flex-shrink-0 z-30">
          <div className="flex flex-col gap-2">
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-4">
                {/* Hamburger Menu */}
                <button
                  onClick={() => setSidebarOpen(true)}
                  className="lg:hidden p-2 rounded-lg hover:bg-slate-100"
                >
                  <svg className="w-6 h-6 text-slate-700" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" />
                  </svg>
                </button>

                {/* Page Title */}
                <div>
                  <h2 className="text-lg font-semibold text-slate-800">
                    {currentStep === 1 && 'Seleccionar Genomas'}
                    {currentStep === 2 && 'Ejecutar Análisis'}
                    {currentStep === 3 && navigationSections.flatMap(s => s.items).find(item => item.id === activeView)?.name}
                  </h2>
                  {currentStep === 3 && (
                    <p className="text-xs text-slate-500">
                      {navigationSections.flatMap(s => s.items).find(item => item.id === activeView)?.description}
                    </p>
                  )}
                </div>
              </div>

              {/* Actions */}
              <div className="flex items-center gap-2">
                {currentStep === 3 && (
                  <button
                    onClick={() => setCurrentStep(2)}
                    className="hidden sm:flex items-center gap-2 px-3 py-2 text-sm text-slate-600 hover:text-slate-800 hover:bg-slate-100 rounded-lg transition-all"
                  >
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                    </svg>
                    Nuevo análisis
                  </button>
                )}
              </div>
            </div>

            {/* Breadcrumbs */}
            <nav className="flex items-center gap-2 text-sm overflow-x-auto">
              <button
                onClick={() => setCurrentStep(1)}
                className="text-slate-500 hover:text-slate-700 transition-colors whitespace-nowrap"
              >
                Inicio
              </button>
              {currentStep >= 2 && (
                <>
                  <svg className="w-4 h-4 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                  </svg>
                  <button
                    onClick={() => setCurrentStep(2)}
                    className={`transition-colors whitespace-nowrap ${
                      currentStep === 2 ? 'text-slate-800 font-medium' : 'text-slate-500 hover:text-slate-700'
                    }`}
                  >
                    Análisis
                  </button>
                </>
              )}
              {currentStep === 3 && (
                <>
                  <svg className="w-4 h-4 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                  </svg>
                  <span className="text-slate-800 font-medium whitespace-nowrap">
                    {navigationSections.flatMap(s => s.items).find(item => item.id === activeView)?.name || 'Resultados'}
                  </span>
                </>
              )}
            </nav>
          </div>
        </header>

        {/* Main Content */}
        <main className="flex-1 overflow-y-auto p-4 lg:p-6 min-h-0">

          {/* PASO 1: Seleccionar Múltiples Genomas */}
          {currentStep === 1 && (
            <div className="max-w-6xl mx-auto">
              <div className="mb-6">
                <p className="text-slate-600">
                  Selecciona uno o más genomas para analizar. Los genomas se descargarán automáticamente desde NCBI.
                </p>
              </div>

              <GenomeMultiSelector
                downloadedGenomes={downloadedGenomes}
                selectedGenomes={selectedGenomes}
                onSelectionChange={handleMultipleSelection}
                onRefresh={loadDownloadedGenomes}
              />

              {/* Botón de continuar */}
              <div className="flex justify-end items-center pt-4 bg-white rounded-xl border border-slate-200 p-4 mt-6">
                <button
                  onClick={() => {
                    if (selectedGenomes.length === 0) {
                      toast.error('Selecciona al menos un genoma')
                      return
                    }
                    setCurrentStep(2)
                  }}
                  disabled={selectedGenomes.length === 0}
                  className="flex items-center gap-2 px-6 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-lg font-medium hover:shadow-lg disabled:opacity-50 transition-all"
                >
                  Continuar con {selectedGenomes.length} genoma{selectedGenomes.length !== 1 ? 's' : ''}
                  <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7l5 5m0 0l-5 5m5-5H6" />
                  </svg>
                </button>
              </div>
            </div>
          )}

          {/* PASO 2: Descargar y Analizar */}
          {currentStep === 2 && (
            <div className="max-w-4xl mx-auto">
              <div className="mb-6">
                <p className="text-slate-600">
                  {selectedGenomes.length === 1
                    ? 'Ejecuta el análisis completo del genoma seleccionado'
                    : `Ejecuta el análisis comparativo de ${selectedGenomes.length} genomas`}
                </p>
              </div>

              {/* Resumen de genomas seleccionados */}
              {selectedGenomes.length > 0 ? (
                <div className="bg-white rounded-xl border border-slate-200 p-6">
                  <div className="mb-6">
                    <div className="inline-flex items-center gap-2 bg-teal-50 px-4 py-2 rounded-lg">
                      <span className="font-semibold text-teal-700">
                        {selectedGenomes.length} genoma{selectedGenomes.length !== 1 ? 's' : ''} seleccionado{selectedGenomes.length !== 1 ? 's' : ''}
                      </span>
                    </div>
                    <p className="text-slate-600 text-sm mt-2">
                      {selectedGenomes.length === 1
                        ? 'Análisis completo de codones, genes, estadísticas y validación con IA'
                        : 'Análisis comparativo de tamaño, genes, GC%, densidad génica y genes extremos'}
                    </p>
                  </div>

                  {/* Lista de genomas */}
                  <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3 mb-6">
                    {selectedGenomes.map((accession, i) => {
                      const genome = downloadedGenomes.find(g => g.accession === accession)
                      return (
                        <div key={accession} className="flex items-center gap-3 p-3 bg-slate-50 rounded-lg border border-slate-200">
                          <span className="w-8 h-8 bg-teal-600 text-white rounded-full flex items-center justify-center text-sm font-bold flex-shrink-0">
                            {i + 1}
                          </span>
                          <div className="flex-1 min-w-0">
                            <div className="font-mono text-sm text-teal-700 truncate font-medium">{accession}</div>
                            {genome?.organism_name && (
                              <div className="text-xs text-slate-500 truncate">{genome.organism_name}</div>
                            )}
                          </div>
                        </div>
                      )
                    })}
                  </div>

                  {/* Botón de análisis */}
                  <button
                    onClick={runAnalysis}
                    disabled={isLoading}
                    className="w-full py-4 bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700 text-white font-semibold rounded-lg transition-all shadow hover:shadow-lg disabled:opacity-50 flex items-center justify-center gap-3"
                  >
                    {isLoading ? (
                      <>
                        <svg className="animate-spin h-5 w-5" viewBox="0 0 24 24">
                          <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                          <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                        </svg>
                        <span>Analizando genomas...</span>
                      </>
                    ) : (
                      <>
                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4" />
                        </svg>
                        <span>
                          {selectedGenomes.length === 1
                            ? 'Ejecutar Análisis Completo'
                            : 'Ejecutar Análisis Comparativo'}
                        </span>
                      </>
                    )}
                  </button>
                </div>
              ) : (
                <div className="bg-white rounded-xl border border-slate-200 p-12 text-center">
                  <div className="w-16 h-16 mx-auto bg-amber-50 rounded-xl flex items-center justify-center mb-4">
                    <svg className="w-8 h-8 text-amber-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                    </svg>
                  </div>
                  <h3 className="font-semibold text-slate-700 mb-2">No hay genomas seleccionados</h3>
                  <p className="text-slate-500 text-sm mb-4">Selecciona al menos un genoma para analizar</p>
                  <button
                    onClick={() => setCurrentStep(1)}
                    className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 transition-all text-sm font-medium"
                  >
                    Ir a selección
                  </button>
                </div>
              )}

              <div className="flex justify-start pt-4">
                <button
                  onClick={() => setCurrentStep(1)}
                  className="flex items-center gap-2 px-4 py-2 text-slate-600 hover:text-slate-800 text-sm hover:bg-slate-100 rounded-lg transition-all"
                >
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
                  </svg>
                  Volver a selección
                </button>
              </div>
            </div>
          )}

          {/* PASO 3: Resultados */}
          {currentStep === 3 && (
            <div className="max-w-full">
              {/* Genome Info Bar */}
              {currentGenome && selectedGenomes.length > 1 && (
                <div className="bg-teal-50 border border-teal-100 rounded-lg p-4 mb-6">
                  <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
                    <div>
                      <div className="text-sm text-slate-700 mb-1">
                        <span className="font-medium">Comparación:</span> {selectedGenomes.length} genomas
                      </div>
                      <p className="text-xs text-slate-600">
                        Las vistas de Dashboard, Codones, Genes y Validación IA muestran análisis detallado de {currentGenome.accession}
                      </p>
                    </div>
                    <button
                      onClick={() => setCurrentStep(1)}
                      className="text-sm text-teal-600 hover:text-teal-700 font-medium whitespace-nowrap"
                    >
                      Cambiar selección
                    </button>
                  </div>
                </div>
              )}

              {/* Content */}
              <div className="bg-white rounded-xl border border-slate-200 p-6">
                {activeView === 'dashboard' && (
                  <AnalysisDashboard
                    analysisData={analysisData}
                    isLoading={isLoading}
                    status={analysisData ? 'completed' : 'idle'}
                  />
                )}
                {activeView === 'codons' && (
                  <CodonVisualization codonData={analysisData?.codons} />
                )}
                {activeView === 'comparison' && (
                  <ComparisonResults
                    comparisonResult={comparisonResult}
                    selectedGenomes={selectedGenomes}
                  />
                )}
                {activeView === 'genes' && (
                  <GeneStatistics geneData={analysisData?.genes} />
                )}
                {activeView === 'filter' && (
                  <GeneFilter hasAnalysis={!!analysisData} />
                )}
                {activeView === 'files' && (
                  <FileManager
                    files={files}
                    onRefresh={loadFiles}
                    selectedGenomes={selectedGenomes}
                  />
                )}
                {activeView === 'genome-map' && (
                  <GenomeViewer />
                )}
                {activeView === 'dogma' && (
                  <CentralDogma />
                )}
                {activeView === 'proteins' && (
                  <ProteinViewer />
                )}
                {activeView === 'codon-usage' && (
                  <CodonUsageTable />
                )}
                {activeView === 'gc-window' && (
                  <GCWindowViewer />
                )}
                {activeView === 'rna' && (
                  <RNAAnalysis />
                )}
                {activeView === 'blast' && (
                  <BLASTSearch />
                )}
                {activeView === 'phylo' && (
                  <PhylogeneticTree />
                )}
                {activeView === 'cog' && (
                  <FunctionalCategories />
                )}
                {activeView === 'cai' && (
                  <CAIAnalysis />
                )}
                {activeView === 'chat' && (
                  <InteractiveChat
                    hasAnalysis={!!analysisData}
                    currentGenome={currentGenome}
                  />
                )}
                {activeView === 'ai' && (
                  <AIValidation
                    validationData={aiValidation}
                    isValidating={isValidatingAI}
                    onValidate={runAIValidation}
                    hasAnalysis={!!analysisData}
                    comparisonData={comparisonResult}
                    selectedGenomes={selectedGenomes}
                    currentGenome={currentGenome}
                  />
                )}
                {activeView === 'export' && (
                  <DataExport
                    hasData={analysisData !== null}
                    comparisonData={comparisonResult}
                    currentGenome={currentGenome}
                    selectedGenomes={selectedGenomes}
                  />
                )}
              </div>
            </div>
          )}
        </main>

        {/* Footer */}
        <footer className="bg-slate-800 text-slate-400 py-3 px-6 border-t border-slate-700 flex-shrink-0">
          <div className="flex flex-col sm:flex-row justify-between items-center gap-2 text-sm">
            <div className="flex items-center gap-2">
              <svg className="w-4 h-4 text-teal-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
              <span className="text-xs">Sistema de Análisis Genómico</span>
            </div>
            <div className="flex items-center gap-3 text-xs">
              <span>FastAPI + BioPython</span>
              <span>•</span>
              <span>React + Vite</span>
            </div>
          </div>
        </footer>
      </div>
    </div>
  )
}

export default App
