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

  const resultViews = [
    { id: 'dashboard', name: 'Dashboard' },
    { id: 'comparison', name: '📊 Comparación' },
    { id: 'codons', name: 'Codones' },
    { id: 'genes', name: 'Genes' },
    { id: 'filter', name: 'Filtrar Genes' },
    { id: 'files', name: 'Archivos' },
    { id: 'ai', name: 'Validación IA' },
    { id: 'export', name: 'Exportar' },
  ]

  return (
    <div className="min-h-screen bg-slate-50 flex flex-col">
      <Toaster position="top-right" />

      {/* Header */}
      <header className="bg-gradient-to-r from-teal-700 via-teal-600 to-emerald-600 text-white shadow-lg">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 bg-white/10 rounded-lg flex items-center justify-center">
                <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
              </div>
              <div>
                <h1 className="text-lg font-semibold">Análisis Genómico</h1>
                <p className="text-teal-100 text-xs">NCBI Datasets API v2</p>
              </div>
            </div>

            {/* Status Badge */}
            {currentGenome && (
              <div className="hidden sm:flex items-center gap-2 bg-white/10 px-3 py-1.5 rounded-lg text-sm">
                <span className="w-2 h-2 bg-emerald-400 rounded-full"></span>
                <span className="font-mono">{currentGenome.accession}</span>
              </div>
            )}
          </div>
        </div>
      </header>

      {/* Progress Steps */}
      <div className="bg-white border-b border-slate-200 sticky top-0 z-40">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 py-3">
          <div className="flex items-center justify-between gap-1 sm:gap-2">
            {steps.map((step, index) => (
              <button
                key={step.id}
                onClick={() => {
                  // Solo permitir ir a pasos completados o el actual
                  if (step.id <= currentStep || (step.id === 2 && downloadedGenomes.length > 0)) {
                    setCurrentStep(step.id)
                  }
                }}
                className={`flex-1 flex items-center justify-center sm:justify-start gap-2 py-2 px-2 sm:px-3 rounded-lg transition-all ${step.id === currentStep
                    ? 'bg-teal-50 text-teal-700'
                    : step.id < currentStep
                      ? 'text-teal-600 cursor-pointer hover:bg-teal-50/50'
                      : 'text-slate-400 cursor-not-allowed'
                  }`}
              >
                <div className={`w-7 h-7 sm:w-8 sm:h-8 flex-shrink-0 rounded-full flex items-center justify-center text-xs sm:text-sm font-medium ${step.id === currentStep
                    ? 'bg-teal-600 text-white'
                    : step.id < currentStep
                      ? 'bg-teal-100 text-teal-700'
                      : 'bg-slate-100 text-slate-400'
                  }`}>
                  {step.id < currentStep ? (
                    <svg className="w-3 h-3 sm:w-4 sm:h-4" fill="currentColor" viewBox="0 0 20 20">
                      <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                    </svg>
                  ) : step.id}
                </div>
                <div className="hidden md:block text-left flex-1 min-w-0">
                  <div className="text-sm font-medium truncate">{step.name}</div>
                  <div className="text-xs opacity-70 truncate">{step.description}</div>
                </div>
                {index < steps.length - 1 && (
                  <div className={`hidden lg:block flex-1 h-0.5 ml-3 ${step.id < currentStep ? 'bg-teal-300' : 'bg-slate-200'
                    }`} />
                )}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Main Content */}
      <main className="flex-1 max-w-7xl w-full mx-auto px-4 sm:px-6 py-6">

        {/* PASO 1: Seleccionar Múltiples Genomas */}
        {currentStep === 1 && (
          <div className="space-y-6">
            <div className="text-center mb-6 sm:mb-8">
              <h2 className="text-xl sm:text-2xl font-bold text-slate-800 px-4">Seleccionar Genomas para Análisis</h2>
              <p className="text-sm sm:text-base text-slate-500 mt-1 px-4">
                Selecciona uno o más genomas para analizar (se descargarán automáticamente)
              </p>
            </div>

            <GenomeMultiSelector
              downloadedGenomes={downloadedGenomes}
              selectedGenomes={selectedGenomes}
              onSelectionChange={handleMultipleSelection}
              onRefresh={loadDownloadedGenomes}
            />

            {/* Botón de continuar */}
            <div className="flex justify-end items-center pt-4 bg-white rounded-xl border border-slate-200 p-4">
              <button
                onClick={() => {
                  if (selectedGenomes.length === 0) {
                    toast.error('Selecciona al menos un genoma')
                    return
                  }
                  setCurrentStep(2)
                }}
                disabled={selectedGenomes.length === 0}
                className="px-6 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-bold hover:shadow-lg disabled:opacity-50 transition-all"
              >
                Continuar con {selectedGenomes.length} genoma{selectedGenomes.length !== 1 ? 's' : ''} →
              </button>
            </div>
          </div>
        )}

        {/* PASO 2: Descargar y Analizar */}
        {currentStep === 2 && (
          <div className="space-y-6">
            <div className="text-center mb-6 sm:mb-8">
              <h2 className="text-xl sm:text-2xl font-bold text-slate-800 px-4">Preparar Análisis</h2>
              <p className="text-sm sm:text-base text-slate-500 mt-1 px-4">
                {selectedGenomes.length === 1 
                  ? 'Analiza el genoma seleccionado' 
                  : `Análisis comparativo de ${selectedGenomes.length} genomas`}
              </p>
            </div>

            {/* Resumen de genomas seleccionados */}
            {selectedGenomes.length > 0 ? (
              <div className="bg-white rounded-xl border border-slate-200 p-6 max-w-3xl mx-auto">
                <div className="text-center mb-6">
                  <div className="inline-flex items-center gap-2 bg-teal-50 px-4 py-2 rounded-lg mb-2">
                    <span className="text-2xl">🧬</span>
                    <span className="font-bold text-teal-700">
                      {selectedGenomes.length} genoma{selectedGenomes.length !== 1 ? 's' : ''} seleccionado{selectedGenomes.length !== 1 ? 's' : ''}
                    </span>
                  </div>
                  <p className="text-slate-600 text-sm">
                    {selectedGenomes.length === 1 
                      ? 'Análisis completo de un genoma' 
                      : 'Análisis comparativo mostrará diferencias y similitudes'}
                  </p>
                </div>

                {/* Lista de genomas */}
                <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3 mb-6">
                  {selectedGenomes.map((accession, i) => {
                    const genome = downloadedGenomes.find(g => g.accession === accession)
                    return (
                      <div key={accession} className="flex items-center gap-2 p-3 bg-slate-50 rounded-lg">
                        <span className="w-6 h-6 bg-teal-600 text-white rounded-full flex items-center justify-center text-xs font-bold">
                          {i + 1}
                        </span>
                        <div className="flex-1 min-w-0">
                          <div className="font-mono text-sm text-teal-700 truncate">{accession}</div>
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
                  className="w-full py-4 bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700 text-white text-lg font-bold rounded-xl transition-all shadow-lg hover:shadow-xl disabled:opacity-50 flex items-center justify-center gap-3"
                >
                  {isLoading ? (
                    <>
                      <svg className="animate-spin h-6 w-6" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                      </svg>
                      <span>Analizando...</span>
                    </>
                  ) : (
                    <>
                      <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
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

                <p className="text-center text-xs text-slate-500 mt-3">
                  {selectedGenomes.length === 1 
                    ? 'Análisis de codones, genes, estadísticas y validación con IA' 
                    : 'Comparación de tamaño, genes, GC%, densidad génica y genes extremos'}
                </p>
              </div>
            ) : (
              <div className="bg-white rounded-xl border border-slate-200 p-12 text-center max-w-2xl mx-auto">
                <div className="w-16 h-16 mx-auto bg-amber-50 rounded-xl flex items-center justify-center mb-4">
                  <svg className="w-8 h-8 text-amber-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                  </svg>
                </div>
                <h3 className="font-semibold text-slate-700 mb-2">No hay genomas seleccionados</h3>
                <p className="text-slate-500 text-sm mb-4">Selecciona un genoma para analizar</p>
                <button
                  onClick={() => setCurrentStep(1)}
                  className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 transition-all text-sm font-medium"
                >
                  Seleccionar genoma
                </button>
              </div>
            )}

            <div className="flex justify-between pt-4 max-w-2xl mx-auto">
              <button
                onClick={() => setCurrentStep(1)}
                className="px-4 py-2 text-slate-600 hover:text-slate-800 text-sm"
              >
                ← Cambiar genoma
              </button>
            </div>
          </div>
        )}

        {/* PASO 3: Resultados */}
        {currentStep === 3 && (
          <div className="space-y-4 sm:space-y-6">
            {/* Tabs de resultados */}
            <div className="bg-white border border-slate-200 p-2 rounded-xl">
              <div className="flex flex-nowrap overflow-x-auto gap-2 pb-1 scrollbar-visible">
                {resultViews.map(view => (
                  <button
                    key={view.id}
                    onClick={() => setActiveView(view.id)}
                    className={`flex-shrink-0 px-3 sm:px-4 py-2 rounded-lg text-xs sm:text-sm font-medium transition-all whitespace-nowrap ${activeView === view.id
                        ? 'bg-teal-600 text-white'
                        : 'text-slate-600 hover:bg-slate-100'
                      }`}
                  >
                    {view.name}
                  </button>
                ))}

                {/* AI Validation Button */}
                <button
                  onClick={runAIValidation}
                  disabled={isValidatingAI || !analysisData}
                  className="flex-shrink-0 ml-auto px-3 sm:px-4 py-2 bg-gradient-to-r from-teal-500 to-emerald-500 text-white rounded-lg text-xs sm:text-sm font-medium hover:from-teal-600 hover:to-emerald-600 disabled:opacity-50 transition-all flex items-center gap-2 whitespace-nowrap"
                >
                  {isValidatingAI ? (
                    <>
                      <svg className="animate-spin h-3 w-3 sm:h-4 sm:w-4" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                      </svg>
                      <span className="hidden sm:inline">Validando...</span>
                    </>
                  ) : (
                    <>
                      <svg className="w-3 h-3 sm:w-4 sm:h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
                      </svg>
                      <span>Validar IA</span>
                    </>
                  )}
                </button>
              </div>
            </div>

            {/* Genome Info Bar */}
            {currentGenome && (
              <div className="bg-teal-50 border border-teal-100 rounded-xl p-3 sm:p-4 flex flex-col sm:flex-row sm:items-center justify-between gap-3 sm:gap-4">
                <div className="flex-1 min-w-0">
                  <div className="flex items-center gap-2 sm:gap-3 mb-1">
                    <span className="w-2 h-2 flex-shrink-0 bg-emerald-500 rounded-full"></span>
                    <span className="font-mono font-semibold text-teal-700 text-sm sm:text-base">{currentGenome.accession}</span>
                    <span className="text-slate-600 text-xs sm:text-sm truncate">{currentGenome.organism_name}</span>
                  </div>
                  {selectedGenomes.length > 1 && (
                    <p className="text-xs text-slate-500 ml-4">
                      📊 Comparación: {selectedGenomes.length} genomas • 
                      Dashboard/Codones/Genes/IA: análisis detallado de este genoma
                    </p>
                  )}
                </div>
                <button
                  onClick={() => setCurrentStep(1)}
                  className="text-xs sm:text-sm text-teal-600 hover:text-teal-700 font-medium whitespace-nowrap"
                >
                  Cambiar genoma
                </button>
              </div>
            )}

            {/* Content */}
            <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-6">
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
              {activeView === 'ai' && (
                <AIValidation
                  validationData={aiValidation}
                  isValidating={isValidatingAI}
                  onValidate={runAIValidation}
                  hasAnalysis={!!analysisData}
                  comparisonData={comparisonResult}
                  selectedGenomes={selectedGenomes}
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

            <div className="flex justify-between pt-4">
              <button
                onClick={() => {
                  setCurrentStep(2)
                  setAnalysisData(null)
                  setAiValidation(null)
                }}
                className="px-4 py-2 text-slate-600 hover:text-slate-800 text-sm"
              >
                ← Ejecutar nuevo análisis
              </button>
            </div>
          </div>
        )}
      </main>

      {/* Footer */}
      <footer className="bg-slate-800 text-slate-400 py-6 mt-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6">
          <div className="flex flex-col sm:flex-row justify-between items-center gap-4 text-sm">
            <div className="flex items-center gap-2">
              <svg className="w-5 h-5 text-teal-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
              <span>Sistema de Análisis Genómico</span>
            </div>
            <div className="flex items-center gap-4">
              <span>FastAPI + BioPython</span>
              <span>•</span>
              <span>React + Vite</span>
            </div>
          </div>
        </div>
      </footer>
    </div>
  )
}

export default App
