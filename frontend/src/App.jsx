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
    loadFiles()
    setAnalysisData(null)
    setAiValidation(null)
    setCurrentStep(3) // Avanzar al paso de análisis
    toast.success(`Genoma ${result.genome_info.accession} seleccionado`)
  }

  const runAnalysis = async () => {
    if (!currentGenome) {
      toast.error('Seleccione un genoma primero')
      return
    }

    setIsLoading(true)
    try {
      toast.loading('Ejecutando análisis genómico...', { id: 'analysis' })
      const result = await api.runCompleteAnalysis()
      setAnalysisData(result)
      setCurrentStep(4) // Avanzar a resultados
      setActiveView('dashboard')
      toast.success('Análisis completado', { id: 'analysis' })
    } catch (error) {
      console.error('Analysis error:', error)
      toast.error('Error en el análisis', { id: 'analysis' })
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
    { id: 1, name: 'Buscar Genoma', description: 'Buscar en NCBI' },
    { id: 2, name: 'Seleccionar', description: 'Elegir genoma' },
    { id: 3, name: 'Analizar', description: 'Ejecutar análisis' },
    { id: 4, name: 'Resultados', description: 'Ver datos' },
  ]

  const resultViews = [
    { id: 'dashboard', name: 'Dashboard' },
    { id: 'codons', name: 'Codones' },
    { id: 'genes', name: 'Genes' },
    { id: 'filter', name: 'Filtrar Genes' },
    { id: 'compare', name: 'Comparar' },
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

        {/* PASO 1: Buscar Genoma */}
        {currentStep === 1 && (
          <div className="space-y-6">
            <div className="text-center mb-6 sm:mb-8">
              <h2 className="text-xl sm:text-2xl font-bold text-slate-800 px-4">Buscar Genoma en NCBI</h2>
              <p className="text-sm sm:text-base text-slate-500 mt-1 px-4">Escribe el nombre del organismo o número de accesión</p>
            </div>

            <GenomeSelector
              onGenomeActivated={handleGenomeActivated}
              onGenomeDownloaded={handleGenomeDownloaded}
              mode="search"
            />

            {/* Quick Access to Downloaded */}
            {downloadedGenomes.length > 0 && (
              <div className="bg-white rounded-xl border border-slate-200 p-5 mt-6">
                <div className="flex items-center justify-between mb-4">
                  <h3 className="font-semibold text-slate-800">Genomas Descargados</h3>
                  <button
                    onClick={() => setCurrentStep(2)}
                    className="text-sm text-teal-600 hover:text-teal-700 font-medium"
                  >
                    Ver todos →
                  </button>
                </div>
                <div className="flex flex-wrap gap-2">
                  {downloadedGenomes.slice(0, 5).map(genome => (
                    <button
                      key={genome.accession}
                      onClick={() => handleGenomeActivated({ genome_info: genome })}
                      className="px-3 py-1.5 bg-slate-100 hover:bg-teal-100 text-slate-700 hover:text-teal-700 rounded-lg text-sm font-mono transition-all"
                    >
                      {genome.accession}
                    </button>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}

        {/* PASO 2: Seleccionar Genoma */}
        {currentStep === 2 && (
          <div className="space-y-6">
            <div className="text-center mb-6 sm:mb-8">
              <h2 className="text-xl sm:text-2xl font-bold text-slate-800 px-4">Seleccionar Genoma</h2>
              <p className="text-sm sm:text-base text-slate-500 mt-1 px-4">Elige el genoma que deseas analizar</p>
            </div>

            {downloadedGenomes.length === 0 ? (
              <div className="bg-white rounded-xl border border-slate-200 p-12 text-center">
                <div className="w-16 h-16 mx-auto bg-slate-100 rounded-xl flex items-center justify-center mb-4">
                  <svg className="w-8 h-8 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M20 13V6a2 2 0 00-2-2H6a2 2 0 00-2 2v7m16 0v5a2 2 0 01-2 2H6a2 2 0 01-2-2v-5m16 0h-2.586a1 1 0 00-.707.293l-2.414 2.414a1 1 0 01-.707.293h-3.172a1 1 0 01-.707-.293l-2.414-2.414A1 1 0 006.586 13H4" />
                  </svg>
                </div>
                <h3 className="font-semibold text-slate-700 mb-2">No hay genomas descargados</h3>
                <p className="text-slate-500 text-sm mb-4">Primero busca y descarga un genoma</p>
                <button
                  onClick={() => setCurrentStep(1)}
                  className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 transition-all text-sm font-medium"
                >
                  Ir a buscar
                </button>
              </div>
            ) : (
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                {downloadedGenomes.map(genome => (
                  <div
                    key={genome.accession}
                    className={`bg-white rounded-xl border p-5 transition-all cursor-pointer hover:shadow-md ${currentGenome?.accession === genome.accession
                        ? 'border-teal-500 ring-2 ring-teal-100'
                        : 'border-slate-200 hover:border-teal-300'
                      }`}
                    onClick={() => handleGenomeActivated({ genome_info: genome })}
                  >
                    <div className="flex items-start justify-between mb-3">
                      <span className="font-mono text-teal-700 font-semibold bg-teal-50 px-2 py-1 rounded">
                        {genome.accession}
                      </span>
                      {currentGenome?.accession === genome.accession && (
                        <span className="w-6 h-6 bg-teal-500 rounded-full flex items-center justify-center">
                          <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
                            <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                          </svg>
                        </span>
                      )}
                    </div>
                    <p className="text-sm text-slate-600">{genome.file_count} archivos</p>
                    <button
                      onClick={(e) => {
                        e.stopPropagation()
                        handleGenomeActivated({ genome_info: genome })
                      }}
                      className="mt-3 w-full py-2 bg-teal-600 hover:bg-teal-700 text-white rounded-lg text-sm font-medium transition-all"
                    >
                      Seleccionar y Continuar
                    </button>
                  </div>
                ))}
              </div>
            )}

            <div className="flex justify-between pt-4">
              <button
                onClick={() => setCurrentStep(1)}
                className="px-4 py-2 text-slate-600 hover:text-slate-800 text-sm"
              >
                ← Buscar más genomas
              </button>
            </div>
          </div>
        )}

        {/* PASO 3: Ejecutar Análisis */}
        {currentStep === 3 && (
          <div className="space-y-6">
            <div className="text-center mb-6 sm:mb-8">
              <h2 className="text-xl sm:text-2xl font-bold text-slate-800 px-4">Ejecutar Análisis</h2>
              <p className="text-sm sm:text-base text-slate-500 mt-1 px-4">Analiza el genoma seleccionado</p>
            </div>

            {/* Genome Info Card */}
            {currentGenome ? (
              <div className="bg-white rounded-xl border border-slate-200 p-6 max-w-2xl mx-auto">
                <div className="text-center mb-6">
                  <div className="inline-block bg-teal-50 px-4 py-2 rounded-lg mb-2">
                    <span className="font-mono text-xl text-teal-700 font-bold">{currentGenome.accession}</span>
                  </div>
                  <h3 className="text-lg font-semibold text-slate-800">{currentGenome.organism_name || 'Genoma seleccionado'}</h3>
                  {currentGenome.strain && (
                    <p className="text-slate-500 text-sm">Strain: {currentGenome.strain}</p>
                  )}
                </div>

                {/* Stats Preview */}
                {currentGenome.genome_size_mb && (
                  <div className="grid grid-cols-3 gap-4 mb-6">
                    <div className="text-center p-3 bg-slate-50 rounded-lg">
                      <div className="text-lg font-bold text-teal-700">{currentGenome.genome_size_mb} Mb</div>
                      <div className="text-xs text-slate-500">Tamaño</div>
                    </div>
                    <div className="text-center p-3 bg-slate-50 rounded-lg">
                      <div className="text-lg font-bold text-emerald-700">{currentGenome.gc_percent?.toFixed(1)}%</div>
                      <div className="text-xs text-slate-500">GC</div>
                    </div>
                    <div className="text-center p-3 bg-slate-50 rounded-lg">
                      <div className="text-lg font-bold text-slate-700">{currentGenome.gene_count?.toLocaleString() || 'N/A'}</div>
                      <div className="text-xs text-slate-500">Genes</div>
                    </div>
                  </div>
                )}

                {/* Action Button */}
                <button
                  onClick={runAnalysis}
                  disabled={isLoading}
                  className="w-full py-3 sm:py-4 bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700 disabled:from-slate-400 disabled:to-slate-500 text-white rounded-xl font-semibold text-base sm:text-lg transition-all shadow-lg flex items-center justify-center gap-2 sm:gap-3"
                >
                  {isLoading ? (
                    <>
                      <svg className="animate-spin h-4 w-4 sm:h-5 sm:w-5" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                      </svg>
                      Analizando...
                    </>
                  ) : (
                    <>
                      <svg className="w-5 h-5 sm:w-6 sm:h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z" />
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
                      </svg>
                      <span className="hidden sm:inline">Ejecutar Análisis Completo</span>
                      <span className="sm:hidden">Ejecutar Análisis</span>
                    </>
                  )}
                </button>

                <p className="text-center text-xs sm:text-sm text-slate-500 mt-3 px-2">
                  Se analizarán codones, genes y se validará contra valores de referencia
                </p>
              </div>
            ) : (
              <div className="bg-white rounded-xl border border-slate-200 p-12 text-center max-w-2xl mx-auto">
                <div className="w-16 h-16 mx-auto bg-amber-50 rounded-xl flex items-center justify-center mb-4">
                  <svg className="w-8 h-8 text-amber-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                  </svg>
                </div>
                <h3 className="font-semibold text-slate-700 mb-2">No hay genoma seleccionado</h3>
                <p className="text-slate-500 text-sm mb-4">Selecciona un genoma para analizar</p>
                <button
                  onClick={() => setCurrentStep(2)}
                  className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 transition-all text-sm font-medium"
                >
                  Seleccionar genoma
                </button>
              </div>
            )}

            <div className="flex justify-between pt-4 max-w-2xl mx-auto">
              <button
                onClick={() => setCurrentStep(2)}
                className="px-4 py-2 text-slate-600 hover:text-slate-800 text-sm"
              >
                ← Cambiar genoma
              </button>
            </div>
          </div>
        )}

        {/* PASO 4: Resultados */}
        {currentStep === 4 && (
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
                <div className="flex items-center gap-2 sm:gap-3 min-w-0 flex-1">
                  <span className="w-2 h-2 flex-shrink-0 bg-emerald-500 rounded-full"></span>
                  <span className="font-mono font-semibold text-teal-700 text-sm sm:text-base">{currentGenome.accession}</span>
                  <span className="text-slate-600 text-xs sm:text-sm truncate">{currentGenome.organism_name}</span>
                </div>
                <button
                  onClick={() => setCurrentStep(2)}
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
              {activeView === 'genes' && (
                <GeneStatistics geneData={analysisData?.genes} />
              )}
              {activeView === 'filter' && (
                <GeneFilter hasAnalysis={!!analysisData} />
              )}
              {activeView === 'compare' && (
                <GenomeComparison 
                  hasAnalysis={!!analysisData} 
                  downloadedGenomes={downloadedGenomes}
                />
              )}
              {activeView === 'files' && (
                <FileManager files={files} onRefresh={loadFiles} />
              )}
              {activeView === 'ai' && (
                <AIValidation
                  validationData={aiValidation}
                  isValidating={isValidatingAI}
                  onValidate={runAIValidation}
                  hasAnalysis={!!analysisData}
                />
              )}
              {activeView === 'export' && (
                <DataExport hasData={analysisData !== null} />
              )}
            </div>

            <div className="flex justify-between pt-4">
              <button
                onClick={() => {
                  setCurrentStep(3)
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
