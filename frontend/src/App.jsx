import { useState, useEffect } from 'react'
import toast, { Toaster } from 'react-hot-toast'
import FileManager from './components/FileManager'
import AnalysisDashboard from './components/AnalysisDashboard'
import CodonVisualization from './components/CodonVisualization'
import GeneStatistics from './components/GeneStatistics'
import DataExport from './components/DataExport'
import AIValidation from './components/AIValidation'
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
import ConceptMap from './components/ConceptMap'
import NetworkMap from './components/NetworkMap'
import FloatingMapButton from './components/FloatingMapButton'
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
  const [globalSearch, setGlobalSearch] = useState('')

  // Vista activa para resultados (después del paso 3)
  const [activeView, setActiveView] = useState('dashboard')

  // Efecto para hacer scroll al inicio cuando cambia la vista
  useEffect(() => {
    const mainArea = document.querySelector('main')
    if (mainArea) {
      mainArea.scrollTo({ top: 0, behavior: 'smooth' })
    }
  }, [activeView])

  // Estado para sidebar móvil y preferencias
  const [sidebarOpen, setSidebarOpen] = useState(false)
  const [miniSidebar, setMiniSidebar] = useState(false) // Estado para sidebar minimizado
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
      localStorage.setItem('sidebarPreferences', JSON.stringify({ defaultOpen: sidebarOpen, mini: miniSidebar }))
    }
  }, [sidebarOpen, miniSidebar])

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

    console.log('🚀 [ANALYSIS] Iniciando para:', selectedGenomes)
    setIsLoading(true)
    const failedGenomes = []
    const successfulGenomes = []

    try {
      toast.loading(`Sincronizando ${selectedGenomes.length} genoma${selectedGenomes.length > 1 ? 's' : ''}...`, { id: 'analysis' })

      for (let i = 0; i < selectedGenomes.length; i++) {
        const accession = selectedGenomes[i]
        console.log(`📂 [${i + 1}/${selectedGenomes.length}] Procesando: ${accession}`)

        const isAlreadyDownloaded = downloadedGenomes.some(g => g.accession === accession)

        if (isAlreadyDownloaded) {
          console.log(`✅ [${accession}] Ya disponible localmente.`)
          successfulGenomes.push(accession)
          toast.loading(`Genoma ${i + 1}/${selectedGenomes.length} listo ✓`, { id: 'analysis' })
          continue
        }

        toast.loading(`Descargando ${accession} desde NCBI...`, { id: 'analysis' })

        try {
          console.log(`📡 [${accession}] Solicitando descarga...`)
          const initResponse = await api.downloadGenome({
            accession: accession,
            include_gbff: true,
            include_gff: true,
            include_fasta: true
          })
          console.log(`📥 [${accession}] Respuesta de inicio:`, initResponse)

          let downloadComplete = false
          let attempts = 0
          const maxAttempts = 120 // Aumentamos a 2 minutos para genomas grandes

          console.log(`⏳ [${accession}] Iniciando polling de estado...`)
          while (!downloadComplete && attempts < maxAttempts) {
            await new Promise(resolve => setTimeout(resolve, 1500)) // Esperar 1.5s entre encuestas

            const status = await api.getGenomeDownloadStatus(accession)
            console.log(`📊 [${accession}] Intento ${attempts + 1}: ${status.status} - ${status.message}`)

            if (status.status === 'completed') {
              downloadComplete = true
              successfulGenomes.push(accession)
              toast.loading(`${accession} descargado con éxito ✓`, { id: 'analysis' })
            } else if (status.status === 'error') {
              throw new Error(status.message || 'Error en servidor NCBI')
            }
            attempts++
          }

          if (!downloadComplete) {
            throw new Error(`Tiempo de espera agotado para ${accession}`)
          }
        } catch (downloadError) {
          console.error(`❌ [${accession}] Falló la descarga:`, downloadError)
          failedGenomes.push({ accession, error: downloadError.message })
          toast.error(`${accession} falló: ${downloadError.message.substring(0, 40)}...`, { duration: 4000 })
        }
      }

      if (successfulGenomes.length === 0) {
        throw new Error('No se pudo descargar ningún genoma para el análisis.')
      }

      console.log('🧪 [ANALYSIS] Genomas exitosos:', successfulGenomes)
      toast.loading('Ejecutando motores de comparación...', { id: 'analysis' })

      // Comparar genomas
      const comparison = await api.compareGenomes(successfulGenomes)
      setComparisonResult(comparison)
      console.log('📊 [ANALYSIS] Comparación completada.')

      // Activar el último para análisis detallado
      const lastGenome = successfulGenomes[successfulGenomes.length - 1]
      console.log(`🎯 [ANALYSIS] Activando genoma principal: ${lastGenome}`)

      await api.activateGenome(lastGenome)
      await loadFiles()

      toast.loading('Calculando métricas moleculares finales...', { id: 'analysis' })
      const detailedAnalysis = await api.runCompleteAnalysis()
      setAnalysisData(detailedAnalysis)

      setCurrentStep(3)
      setActiveView('comparison')
      await loadDownloadedGenomes()

      toast.success(
        failedGenomes.length > 0
          ? `Análisis parcial: ${successfulGenomes.length} OK, ${failedGenomes.length} Error`
          : `Análisis completado: ${successfulGenomes.length} genomas procesados`,
        { id: 'analysis' }
      )

    } catch (error) {
      console.error('🔥 [FATAL ERROR]:', error)
      toast.error('Error crítico: ' + (error.response?.data?.detail || error.message), { id: 'analysis' })
    } finally {
      setIsLoading(false)
      console.log('🏁 [ANALYSIS] Flujo terminado.')
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

  // Iconos para la navegación (SVG paths)
  const icons = {
    dashboard: "M4 5a1 1 0 011-1h6a1 1 0 011 1v6a1 1 0 01-1 1H5a1 1 0 01-1-1V5zM14 5a1 1 0 011-1h6a1 1 0 011 1v6a1 1 0 01-1 1h-6a1 1 0 01-1-1V5zM4 14a1 1 0 011-1h6a1 1 0 011 1v6a1 1 0 01-1 1H5a1 1 0 01-1-1v-6zM14 14a1 1 0 011-1h6a1 1 0 011 1v6a1 1 0 01-1 1h-6a1 1 0 01-1-1v-6z",
    comparison: "M8 7h12m0 0l-4-4m4 4l-4 4m0 6H4m0 0l4 4m-4-4l4-4",
    'genome-map': "M9 20l-5.447-2.724A1 1 0 013 16.382V5.618a1 1 0 011.447-.894L9 7m0 13l6-3m-6 3V7m6 10l4.553 2.276A1 1 0 0021 18.382V7.618a1 1 0 00-.553-.894L15 4m0 13V4m0 0L9 7",
    codons: "M3 10h18M3 14h18m-9-4v8m-7-4h14M7 6h10",
    'codon-usage': "M9 17v-2m3 2v-4m3 4v-6m2 10H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z",
    genes: "M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 002 2h2a2 2 0 002-2",
    proteins: "M20 7l-8-4-8 4m16 0l-8 4m8-4v10l-8 4m0-10L4 7m8 4v10M4 7v10l8 4",
    rna: "M13 10V3L4 14h7v7l9-11h-7z", // Using generic thunder here as placeholder
    dogma: "M13 5l7 7-7 7M5 5l7 7-7 7",
    'gc-window': "M7 12l3-3 3 3 4-4M8 21l4-4 4 4M3 4h18M4 4h16v12a1 1 0 01-1 1H5a1 1 0 01-1-1V4z",
    blast: "M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z",
    'concept-map': "M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z",
    'network-map': "M13 10V3L4 14h7v7l9-11h-7z", // Flash
    phylo: "M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253",
    cog: "M11 3.055A9.001 9.001 0 1020.945 13H11V3.055z",
    cai: "M9 7h6m0 10v-3m-3 3h.01M9 17h.01M9 14h.01M12 14h.01M15 11h.01M12 11h.01M9 11h.01M7 21h10a2 2 0 002-2V5a2 2 0 00-2-2H7a2 2 0 00-2 2v14a2 2 0 002 2z",
    chat: "M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z",
    ai: "M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z",
    files: "M5 19a2 2 0 01-2-2V7a2 2 0 012-2h4l2 2h4a2 2 0 012 2v8a2 2 0 01-2 2H5z",
    export: "M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"
  }

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
        { id: 'concept-map', name: 'Mapa Conceptual', description: 'Referencia Interactiva de Biología' },
        { id: 'network-map', name: 'Red Genómica Pro', description: 'Visualización Avanzada de Conceptos' },
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
      <aside className={`fixed lg:static inset-y-0 left-0 z-50 bg-white border-r border-slate-200 transform transition-all duration-500 ease-in-out ${sidebarOpen ? 'translate-x-0' : '-translate-x-full lg:translate-x-0'} flex flex-col h-screen shadow-2xl lg:shadow-none ${miniSidebar ? 'w-24' : 'w-80'}`}>

        {/* Sidebar Header */}
        <div className={`p-6 border-b border-slate-100 bg-white flex items-center ${miniSidebar ? 'justify-center' : 'justify-between'}`}>
          <div className="flex items-center gap-4">
            <div className={`w-10 h-10 bg-blue-600 rounded-2xl flex items-center justify-center shadow-lg shadow-blue-200 transition-transform hover:rotate-12 flex-shrink-0 cursor-pointer`} onClick={() => setMiniSidebar(!miniSidebar)}>
              <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
            </div>
            {!miniSidebar && (
              <div className="min-w-0 animate-in fade-in duration-300">
                <h1 className="text-sm font-black text-slate-900 uppercase tracking-tighter truncate">Genómica <span className="text-blue-600">Lab</span></h1>
                <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest truncate">Análisis Molecular</p>
              </div>
            )}
          </div>
          {!miniSidebar && (
            <button
              onClick={() => {
                if (window.innerWidth < 1024) setSidebarOpen(false)
                else setMiniSidebar(true)
              }}
              className="p-2 rounded-xl hover:bg-slate-100 transition-colors"
            >
              <svg className="w-5 h-5 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d={window.innerWidth < 1024 ? "M6 18L18 6M6 6l12 12" : "M11 19l-7-7 7-7m8 14l-7-7 7-7"} />
              </svg>
            </button>
          )}
        </div>

        {/* Workflow Progress */}
        {currentStep < 3 && !miniSidebar && (
          <div className="p-6 bg-slate-50/50 border-b border-slate-100 flex-shrink-0 animate-in fade-in">
            <div className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] mb-4 px-2">PROGRESO DEL FLUJO</div>
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
                  className={`w-full flex items-center gap-4 px-4 py-3 rounded-2xl text-sm transition-all duration-300 ${step.id === currentStep
                    ? 'bg-blue-600 text-white shadow-xl shadow-blue-200'
                    : step.id < currentStep
                      ? 'bg-blue-50 text-blue-700 hover:bg-blue-100'
                      : 'bg-slate-100/50 text-slate-400 cursor-not-allowed opacity-50'
                    }`}
                >
                  <div className={`w-7 h-7 rounded-xl flex items-center justify-center text-xs font-black ${step.id === currentStep
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
        {currentStep < 3 && !miniSidebar && (
          <div className="px-6 py-4 space-y-6">
            <div className="space-y-2">
              <h2 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-3">
                Acceso Rápido
              </h2>
              <div className="space-y-1">
                <button
                  onClick={() => {
                    setCurrentStep(3)
                    setActiveView('concept-map')
                  }}
                  className="w-full flex items-center px-4 py-3 rounded-2xl text-xs text-slate-600 hover:bg-slate-50 transition-all group"
                >
                  <span className="text-slate-400 group-hover:text-blue-600">
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d={icons['concept-map']} />
                    </svg>
                  </span>
                  <div className="ml-3 text-left">
                    <span className="block font-black uppercase tracking-tight">Mapa Conceptual</span>
                  </div>
                </button>
                <button
                  onClick={() => {
                    setCurrentStep(3)
                    setActiveView('network-map')
                  }}
                  className="w-full flex items-center px-4 py-3 rounded-2xl text-xs text-slate-600 hover:bg-slate-50 transition-all group"
                >
                  <span className="text-slate-400 group-hover:text-emerald-600">
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d={icons['network-map']} />
                    </svg>
                  </span>
                  <div className="ml-3 text-left">
                    <span className="block font-black uppercase tracking-tight">Red Genómica Pro</span>
                  </div>
                </button>
              </div>
            </div>

            <div className="p-4 bg-slate-50 rounded-2xl border border-slate-100">
              <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest text-center leading-relaxed">
                El menú completo se activará al finalizar el análisis
              </p>
            </div>
          </div>
        )}

        {currentStep === 3 && (
          <>
            <div className={`p-6 border-b border-slate-50 flex-shrink-0 transition-all ${miniSidebar ? 'px-4' : 'px-6'}`}>
              <div className="relative group">
                {!miniSidebar ? (
                  <>
                    <input
                      type="text"
                      value={sidebarSearch}
                      onChange={(e) => setSidebarSearch(e.target.value)}
                      placeholder="Filtrar..."
                      className="w-full pl-10 pr-4 py-3 text-xs font-bold border border-slate-100 bg-slate-50 rounded-2xl focus:outline-none focus:ring-4 focus:ring-blue-500/5 focus:border-blue-500/30 transition-all placeholder-slate-400"
                    />
                    <svg className="w-4 h-4 text-slate-400 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-blue-500 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                    </svg>
                  </>
                ) : (
                  <button onClick={() => setMiniSidebar(false)} className="w-full flex justify-center py-2 bg-slate-50 rounded-xl hover:bg-slate-100">
                    <svg className="w-5 h-5 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                    </svg>
                  </button>
                )}
              </div>
            </div>

            <nav className={`flex-1 overflow-y-auto ${miniSidebar ? 'px-2' : 'px-6'} py-6 sidebar-scroll min-h-0 space-y-8`}>
              {navigationSections.map((section, idx) => {
                const filteredItems = section.items.filter(item =>
                  item.name.toLowerCase().includes(sidebarSearch.toLowerCase()) ||
                  item.description.toLowerCase().includes(sidebarSearch.toLowerCase())
                )

                if (filteredItems.length === 0 && sidebarSearch) return null

                return (
                  <div key={idx} className="space-y-2">
                    {!miniSidebar && (
                      <h2 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-3">
                        {section.title}
                      </h2>
                    )}
                    {miniSidebar && idx === 0 && <div className="h-4"></div>} {/* Spacer */}

                    <div className="space-y-1">
                      {filteredItems.map((item) => (
                        <button
                          key={item.id}
                          onClick={() => {
                            setActiveView(item.id)
                            setSidebarSearch('')
                            if (window.innerWidth < 1024) setSidebarOpen(false)
                          }}
                          className={`group w-full flex items-center ${miniSidebar ? 'justify-center py-3 px-0' : 'justify-start px-4 py-3'} rounded-2xl text-xs transition-all duration-300 relative ${activeView === item.id
                            ? 'bg-blue-50 text-blue-700 shadow-sm border border-blue-100/50'
                            : 'text-slate-600 hover:bg-slate-50 hover:text-slate-900'
                            }`}
                          title={miniSidebar ? item.name : ''}
                        >
                          <span className={`flex items-center justify-center ${activeView === item.id ? 'text-blue-600' : 'text-slate-400 group-hover:text-slate-600'}`}>
                            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d={icons[item.id] || "M13 10V3L4 14h7v7l9-11h-7z"} />
                            </svg>
                          </span>

                          {!miniSidebar && (
                            <div className="ml-3 text-left">
                              <span className="block font-black uppercase tracking-tight leading-none">{item.name}</span>
                              <span className="block text-[9px] text-slate-400 font-bold mt-1 uppercase tracking-wider leading-none">{item.description}</span>
                            </div>
                          )}

                          {miniSidebar && (
                            <div className="absolute left-full top-1/2 -translate-y-1/2 ml-4 px-3 py-2 bg-slate-900 text-white text-[10px] font-bold rounded-lg opacity-0 group-hover:opacity-100 pointer-events-none transition-opacity whitespace-nowrap z-50 shadow-xl uppercase tracking-wider">
                              {item.name}
                            </div>
                          )}
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
                  <>
                    <button
                      onClick={() => setActiveView('network-map')}
                      className="hidden sm:flex items-center gap-3 px-5 py-2.5 bg-emerald-500 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-xl shadow-emerald-200 transition-all hover:bg-emerald-600 hover:-translate-y-1 active:scale-95"
                    >
                      <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
                      Mapa Pro
                    </button>
                    
                    <button
                      onClick={() => {
                        setCurrentStep(1)
                        setSelectedGenomes([])
                        setAnalysisData(null)
                        setComparisonResult(null)
                      }}
                      className="flex items-center gap-3 px-5 py-2.5 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-xl shadow-slate-200 transition-all hover:bg-blue-600 hover:-translate-y-1 active:scale-95"
                    >
                      <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M12 4v16m8-8H4" />
                      </svg>
                      Nueva Sesión
                    </button>
                  </>
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

                <div className="flex gap-4">
                  <button
                    onClick={() => {
                      setCurrentStep(3)
                      setActiveView('network-map')
                    }}
                    className="flex items-center gap-3 px-6 py-4 bg-emerald-500 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-xl shadow-emerald-200 transition-all hover:bg-emerald-600 hover:-translate-y-1 active:scale-95"
                  >
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
                    Explorar Mapa Pro
                  </button>
                </div>

                <GenomeMultiSelector
                  downloadedGenomes={downloadedGenomes}
                  selectedGenomes={selectedGenomes}
                  onSelectionChange={handleMultipleSelection}
                  onRefresh={loadDownloadedGenomes}
                />

                <div className="flex justify-end pt-6">
                  <button
                    onClick={() => {
                      if (selectedGenomes.length === 0) {
                        toast.error('Selecciona al menos un genoma')
                        return
                      }
                      setCurrentStep(2)
                    }}
                    disabled={selectedGenomes.length === 0}
                    className="flex items-center gap-4 px-8 py-3.5 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest transition-all hover:bg-blue-600 hover:-translate-y-1 shadow-xl active:scale-95 disabled:opacity-30 disabled:cursor-not-allowed"
                  >
                    <span>Siguiente Fase: {selectedGenomes.length} Genomas</span>
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M17 8l4 4m0 0l-4 4m4-4H3" />
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

                <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-12 shadow-2xl space-y-10">
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

                  <div className="flex justify-center pt-4">
                    <button
                      onClick={runAnalysis}
                      disabled={isLoading}
                      className="w-full max-w-md py-4.5 bg-slate-900 text-white font-black uppercase tracking-[0.2em] rounded-2xl transition-all hover:bg-blue-600 hover:-translate-y-1 shadow-2xl shadow-blue-200 disabled:opacity-50 flex items-center justify-center gap-6 group"
                    >
                      {isLoading ? (
                        <>
                          <div className="w-5 h-5 border-4 border-white/20 border-t-white rounded-full animate-spin"></div>
                          <span className="text-[10px] tracking-widest">Sincronizando...</span>
                        </>
                      ) : (
                        <>
                          <span className="text-[10px]">Iniciar Procesamiento Genómico</span>
                          <svg className="w-5 h-5 group-hover:translate-x-2 transition-transform" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M13 10V3L4 14h7v7l9-11h-7z" />
                          </svg>
                        </>
                      )}
                    </button>
                  </div>
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

                <div className={`bg-white rounded-[3rem] border border-slate-100 shadow-sm min-h-[800px] ${activeView === 'network-map' ? 'p-0 overflow-hidden' : 'p-10'}`}>
                  {activeView === 'dashboard' && <AnalysisDashboard analysisData={analysisData} isLoading={isLoading} status={analysisData ? 'completed' : 'idle'} />}
                  {activeView === 'codons' && <CodonVisualization codonData={analysisData?.codons} />}
                  {activeView === 'comparison' && <ComparisonResults comparisonResult={comparisonResult} selectedGenomes={selectedGenomes} />}
                  {activeView === 'genes' && <GeneStatistics geneData={analysisData?.genes} />}
                  {activeView === 'files' && <FileManager files={files} onRefresh={loadFiles} selectedGenomes={selectedGenomes} />}
                  {activeView === 'genome-map' && <GenomeViewer />}
                  {activeView === 'dogma' && <CentralDogma />}
                  {activeView === 'proteins' && <ProteinViewer />}
                  {activeView === 'codon-usage' && <CodonUsageTable />}
                  {activeView === 'gc-window' && <GCWindowViewer />}
                  {activeView === 'rna' && <RNAAnalysis />}
                  {activeView === 'blast' && <BLASTSearch genes={analysisData?.genes?.genes || []} />}
                  {activeView === 'phylo' && <PhylogeneticTree />}
                  {activeView === 'cog' && <FunctionalCategories />}
                  {activeView === 'cai' && <CAIAnalysis />}
                  {activeView === 'concept-map' && <ConceptMap />}
                  {activeView === 'network-map' && <NetworkMap />}
                  {activeView === 'chat' && <InteractiveChat hasAnalysis={!!analysisData} currentGenome={currentGenome} />}
                  {activeView === 'ai' && (
                    <AIValidation
                      validationData={aiValidation}
                      technicalValidation={analysisData?.validation}
                      isValidating={isValidatingAI}
                      onValidate={runAIValidation}
                      hasAnalysis={!!analysisData}
                      comparisonData={comparisonResult}
                      selectedGenomes={selectedGenomes}
                      currentGenome={currentGenome}
                    />
                  )}
                  {activeView === 'export' && <DataExport hasData={analysisData !== null} comparisonData={comparisonResult} currentGenome={currentGenome} selectedGenomes={selectedGenomes} />}
                </div>
              </div>
            )}
          </div>
        </main>

        <FloatingMapButton 
          activeView={activeView}
          onOpenConcept={() => {
            setCurrentStep(3)
            setActiveView('concept-map')
          }}
          onOpenNetwork={() => {
            setCurrentStep(3)
            setActiveView('network-map')
          }}
        />

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
