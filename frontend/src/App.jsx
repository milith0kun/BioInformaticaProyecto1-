import { useState, useEffect } from 'react'
import toast from 'react-hot-toast'
import FileManager from './components/FileManager'
import AnalysisDashboard from './components/AnalysisDashboard'
import CodonVisualization from './components/CodonVisualization'
import GeneStatistics from './components/GeneStatistics'
import DataExport from './components/DataExport'
import AIValidation from './components/AIValidation'
import { api } from './services/api'

function App() {
  const [activeTab, setActiveTab] = useState('dashboard')
  const [files, setFiles] = useState([])
  const [analysisData, setAnalysisData] = useState(null)
  const [aiValidation, setAiValidation] = useState(null)
  const [isLoading, setIsLoading] = useState(false)
  const [isValidatingAI, setIsValidatingAI] = useState(false)
  const [analysisStatus, setAnalysisStatus] = useState('idle')

  // Load files on mount
  useEffect(() => {
    loadFiles()
  }, [])

  const loadFiles = async () => {
    try {
      const response = await api.detectFiles()
      setFiles(response.files)
      if (response.files.length > 0) {
        toast.success(`${response.files.length} archivos detectados`)
      }
    } catch (error) {
      console.error('Error loading files:', error)
    }
  }

  const runAnalysis = async () => {
    setIsLoading(true)
    setAnalysisStatus('running')
    
    try {
      toast.loading('Ejecutando análisis completo...', { id: 'analysis' })
      const result = await api.runCompleteAnalysis()
      setAnalysisData(result)
      setAnalysisStatus('completed')
      toast.success('Análisis completado exitosamente', { id: 'analysis' })
    } catch (error) {
      console.error('Analysis error:', error)
      setAnalysisStatus('error')
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
      toast.loading('Validando con IA (Google Gemini)...', { id: 'ai-validation' })
      const result = await api.validateWithAI()
      setAiValidation(result)
      toast.success('Validación IA completada', { id: 'ai-validation' })
    } catch (error) {
      console.error('AI Validation error:', error)
      toast.error('Error en validación IA: ' + (error.response?.data?.detail || error.message), { id: 'ai-validation' })
    } finally {
      setIsValidatingAI(false)
    }
  }

  const tabs = [
    { id: 'dashboard', name: 'Dashboard' },
    { id: 'files', name: 'Archivos' },
    { id: 'codons', name: 'Codones' },
    { id: 'genes', name: 'Genes' },
    { id: 'ai', name: 'Validación IA' },
    { id: 'export', name: 'Exportar' },
  ]

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <header className="bg-gradient-to-r from-emerald-600 to-cyan-600 text-white shadow-lg">
        <div className="max-w-7xl mx-auto px-4 py-6">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-4">
              <img src="/dna.svg" alt="DNA" className="h-12 w-12" />
              <div>
                <h1 className="text-2xl font-bold">E. coli K-12 MG1655</h1>
                <p className="text-emerald-100">Sistema de Análisis Genómico</p>
              </div>
            </div>
            <button
              onClick={runAnalysis}
              disabled={isLoading}
              className={`px-8 py-3 rounded-xl font-bold text-sm uppercase tracking-wide transition-all duration-200 ${
                isLoading
                  ? 'bg-gray-400 cursor-not-allowed text-gray-200'
                  : 'bg-white text-emerald-600 hover:bg-emerald-50 hover:shadow-xl hover:scale-105'
              }`}
            >
              {isLoading ? (
                <span className="flex items-center gap-3">
                  <svg className="animate-spin h-5 w-5" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                  </svg>
                  Analizando...
                </span>
              ) : (
                'Ejecutar Análisis'
              )}
            </button>
          </div>
        </div>
      </header>

      {/* Navigation Tabs */}
      <nav className="bg-white shadow-sm border-b border-gray-200">
        <div className="max-w-7xl mx-auto px-4">
          <div className="flex space-x-2">
            {tabs.map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`relative px-6 py-4 font-semibold text-sm uppercase tracking-wide transition-all duration-200 ${
                  activeTab === tab.id
                    ? 'text-emerald-600'
                    : 'text-gray-500 hover:text-gray-700 hover:bg-gray-50'
                }`}
              >
                {tab.name}
                {activeTab === tab.id && (
                  <div className="absolute bottom-0 left-0 right-0 h-0.5 bg-gradient-to-r from-emerald-500 to-cyan-500"></div>
                )}
              </button>
            ))}
          </div>
        </div>
      </nav>

      {/* Main Content */}
      <main className="max-w-7xl mx-auto px-4 py-8">
        {activeTab === 'dashboard' && (
          <AnalysisDashboard 
            analysisData={analysisData} 
            isLoading={isLoading}
            status={analysisStatus}
          />
        )}
        {activeTab === 'files' && (
          <FileManager 
            files={files} 
            onRefresh={loadFiles}
          />
        )}
        {activeTab === 'codons' && (
          <CodonVisualization 
            codonData={analysisData?.codons}
          />
        )}
        {activeTab === 'genes' && (
          <GeneStatistics 
            geneData={analysisData?.genes}
          />
        )}
        {activeTab === 'ai' && (
          <AIValidation
            validationData={aiValidation}
            isValidating={isValidatingAI}
            onValidate={runAIValidation}
            hasAnalysis={!!analysisData}
          />
        )}
        {activeTab === 'export' && (
          <DataExport 
            hasData={analysisData !== null}
          />
        )}
      </main>

      {/* Footer */}
      <footer className="bg-gray-800 text-gray-400 py-6 mt-12">
        <div className="max-w-7xl mx-auto px-4 text-center">
          <p>Proyecto de Bioinformática - Análisis de E. coli K-12 MG1655</p>
          <p className="text-sm mt-2">Datos: NCBI RefSeq NC_000913.3</p>
        </div>
      </footer>
    </div>
  )
}

export default App
