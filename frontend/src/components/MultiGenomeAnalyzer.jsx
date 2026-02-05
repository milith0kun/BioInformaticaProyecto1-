/**
 * MultiGenomeAnalyzer Component
 * Sistema completo para seleccionar múltiples genomas, descargarlos,
 * compararlos y mostrar resultados consolidados con genes extremos
 */
import { useState, useEffect, useMemo, useCallback } from 'react'
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Legend,
  PieChart,
  Pie,
  Cell,
  LineChart,
  Line
} from 'recharts'
import { AgGridReact } from 'ag-grid-react'
import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'
import { api } from '../services/api'
import toast from 'react-hot-toast'

// Colores para gráficos
const COLORS = ['#0d9488', '#10b981', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#14b8a6']

export default function MultiGenomeAnalyzer() {
  // Estados principales
  const [step, setStep] = useState(1) // 1: Selección, 2: Descarga, 3: Análisis, 4: Resultados
  const [selectedGenomes, setSelectedGenomes] = useState([]) // Genomas seleccionados para descargar
  const [downloadedGenomes, setDownloadedGenomes] = useState([]) // Genomas ya descargados
  const [downloadProgress, setDownloadProgress] = useState({}) // Progreso de descarga por accession
  const [comparisonResult, setComparisonResult] = useState(null) // Resultado de comparación
  const [isLoading, setIsLoading] = useState(false)
  const [activeTab, setActiveTab] = useState('summary')
  
  // Estados de búsqueda
  const [searchQuery, setSearchQuery] = useState('')
  const [searchResults, setSearchResults] = useState([])
  const [isSearching, setIsSearching] = useState(false)
  
  // Estados de filtros para genes
  const [geneFilter, setGeneFilter] = useState('all')
  const [geneSizeOrder, setGeneSizeOrder] = useState('largest')
  const [geneCount, setGeneCount] = useState(20)
  const [functionalGroups, setFunctionalGroups] = useState([])
  const [selectedGroup, setSelectedGroup] = useState('')
  
  // Cepas relacionadas predefinidas
  const relatedStrains = [
    { accession: "GCF_000005845.2", organism: "E. coli K-12 MG1655", strain: "K-12", category: "laboratory", size: "4.6 Mb" },
    { accession: "GCF_000008865.2", organism: "E. coli O157:H7 Sakai", strain: "O157:H7", category: "pathogenic", size: "5.5 Mb" },
    { accession: "GCF_000009565.2", organism: "E. coli BL21(DE3)", strain: "BL21", category: "industrial", size: "4.6 Mb" },
    { accession: "GCF_000019425.1", organism: "E. coli CFT073", strain: "CFT073", category: "pathogenic", size: "5.2 Mb" },
    { accession: "GCF_000007445.1", organism: "E. coli W3110", strain: "W3110", category: "laboratory", size: "4.6 Mb" },
    { accession: "GCF_000750555.1", organism: "E. coli Nissle 1917", strain: "Nissle", category: "probiotic", size: "5.0 Mb" },
  ]

  // Cargar genomas descargados y grupos funcionales al montar
  useEffect(() => {
    loadDownloadedGenomes()
    loadFunctionalGroups()
  }, [])

  const loadDownloadedGenomes = async () => {
    try {
      const result = await api.getDownloadedGenomes()
      setDownloadedGenomes(result.genomes || [])
    } catch (error) {
      console.error('Error:', error)
    }
  }

  const loadFunctionalGroups = async () => {
    try {
      const result = await api.getFunctionalGroups()
      setFunctionalGroups(result.groups || [])
    } catch (error) {
      console.error('Error:', error)
    }
  }

  // Búsqueda de genomas
  const searchGenomes = async () => {
    if (!searchQuery.trim() || searchQuery.length < 2) return
    
    setIsSearching(true)
    try {
      const results = await api.searchGenomes(searchQuery, 20)
      setSearchResults(results.genomes || [])
    } catch (error) {
      toast.error('Error en búsqueda: ' + error.message)
    } finally {
      setIsSearching(false)
    }
  }

  // Toggle selección de genoma
  const toggleGenomeSelection = (genome) => {
    setSelectedGenomes(prev => {
      const exists = prev.find(g => g.accession === genome.accession)
      if (exists) {
        return prev.filter(g => g.accession !== genome.accession)
      } else {
        return [...prev, genome]
      }
    })
  }

  // Seleccionar todas las cepas relacionadas
  const selectAllRelated = () => {
    setSelectedGenomes(relatedStrains.map(s => ({
      accession: s.accession,
      organism_name: s.organism,
      strain: s.strain
    })))
  }

  // Descargar todos los genomas seleccionados
  const downloadAllGenomes = async () => {
    if (selectedGenomes.length === 0) {
      toast.error('Seleccione al menos un genoma')
      return
    }

    setStep(2)
    setIsLoading(true)

    for (const genome of selectedGenomes) {
      // Verificar si ya está descargado
      const alreadyDownloaded = downloadedGenomes.some(d => d.accession === genome.accession)
      if (alreadyDownloaded) {
        setDownloadProgress(prev => ({
          ...prev,
          [genome.accession]: { status: 'completed', message: 'Ya descargado' }
        }))
        continue
      }

      setDownloadProgress(prev => ({
        ...prev,
        [genome.accession]: { status: 'downloading', message: 'Descargando...' }
      }))

      try {
        await api.downloadGenome({
          accession: genome.accession,
          include_gbff: true,
          include_gff: true,
          include_fasta: true
        })

        // Esperar a que se complete la descarga
        let completed = false
        let attempts = 0
        while (!completed && attempts < 60) {
          await new Promise(resolve => setTimeout(resolve, 2000))
          try {
            const status = await api.getGenomeDownloadStatus(genome.accession)
            if (status.status === 'completed') {
              completed = true
              setDownloadProgress(prev => ({
                ...prev,
                [genome.accession]: { status: 'completed', message: '✓ Completado' }
              }))
            } else if (status.status === 'error') {
              setDownloadProgress(prev => ({
                ...prev,
                [genome.accession]: { status: 'error', message: status.message }
              }))
              break
            }
          } catch (e) {
            attempts++
          }
          attempts++
        }
      } catch (error) {
        setDownloadProgress(prev => ({
          ...prev,
          [genome.accession]: { status: 'error', message: error.message }
        }))
      }
    }

    await loadDownloadedGenomes()
    setIsLoading(false)
    toast.success('Descargas completadas')
  }

  // Ejecutar comparación de todos los genomas
  const runFullComparison = async () => {
    setStep(3)
    setIsLoading(true)
    
    try {
      toast.loading('Analizando y comparando genomas...', { id: 'compare' })
      
      // Comparar todos los genomas descargados
      const result = await api.compareGenomes()
      setComparisonResult(result)
      
      // Cargar resumen de grupos funcionales
      try {
        const groupsSummary = await api.getGeneGroupsSummary()
        setComparisonResult(prev => ({ ...prev, groupsSummary }))
      } catch (e) {
        console.error('Error loading groups:', e)
      }
      
      setStep(4)
      toast.success(`${result.total_genomes_compared} genomas comparados`, { id: 'compare' })
    } catch (error) {
      toast.error('Error: ' + (error.response?.data?.detail || error.message), { id: 'compare' })
      setStep(2)
    } finally {
      setIsLoading(false)
    }
  }

  // Obtener genes filtrados
  const getFilteredGenes = useCallback(async () => {
    if (!comparisonResult) return

    try {
      if (geneFilter === 'size') {
        const result = await api.getGenesBySize(geneSizeOrder, geneCount)
        setComparisonResult(prev => ({ ...prev, filteredGenes: result.genes }))
      } else if (geneFilter === 'group' && selectedGroup) {
        const result = await api.getGenesByGroup(selectedGroup)
        setComparisonResult(prev => ({ ...prev, filteredGenes: result.genes }))
      }
    } catch (error) {
      console.error('Error filtering genes:', error)
    }
  }, [geneFilter, geneSizeOrder, geneCount, selectedGroup, comparisonResult])

  useEffect(() => {
    if (step === 4 && geneFilter !== 'all') {
      getFilteredGenes()
    }
  }, [geneFilter, geneSizeOrder, geneCount, selectedGroup, step])

  // Datos para gráficos
  const chartData = useMemo(() => {
    if (!comparisonResult?.genomes) return []
    return comparisonResult.genomes.map((g, i) => ({
      name: g.organism_name?.split(' ').slice(0, 3).join(' ') || g.accession,
      accession: g.accession,
      size: (g.genome_length / 1000000).toFixed(2),
      genes: g.gene_count,
      gc: g.gc_content,
      density: g.gene_density,
      avgLength: g.avg_gene_length,
      color: COLORS[i % COLORS.length]
    }))
  }, [comparisonResult])

  // Columnas para tabla de genes
  const geneColumns = [
    { field: 'locus_tag', headerName: 'Locus Tag', width: 120, sortable: true },
    { field: 'product', headerName: 'Producto', flex: 1, minWidth: 250 },
    { field: 'length', headerName: 'Longitud (bp)', width: 120, sortable: true, 
      valueFormatter: p => p.value?.toLocaleString() },
    { field: 'gc_content', headerName: 'GC%', width: 80, sortable: true,
      valueFormatter: p => p.value?.toFixed(1) },
    { field: 'strand', headerName: 'Hebra', width: 70 },
    { field: 'genome', headerName: 'Genoma', width: 150 }
  ]

  // ================== RENDERIZADO ==================

  return (
    <div className="space-y-6">
      {/* Header con progreso */}
      <div className="bg-gradient-to-r from-teal-700 to-emerald-600 rounded-2xl p-6 text-white">
        <h2 className="text-2xl font-bold mb-2">🧬 Análisis Comparativo de Genomas</h2>
        <p className="text-teal-100 mb-4">
          Seleccione múltiples genomas, descárguelos y compare sus características
        </p>
        
        {/* Stepper */}
        <div className="flex items-center justify-between mt-4">
          {[
            { num: 1, label: 'Seleccionar' },
            { num: 2, label: 'Descargar' },
            { num: 3, label: 'Analizar' },
            { num: 4, label: 'Resultados' }
          ].map((s, i) => (
            <div key={s.num} className="flex items-center">
              <div className={`w-10 h-10 rounded-full flex items-center justify-center font-bold ${
                step >= s.num ? 'bg-white text-teal-700' : 'bg-teal-600 text-teal-300'
              }`}>
                {step > s.num ? '✓' : s.num}
              </div>
              <span className={`ml-2 hidden sm:inline ${step >= s.num ? 'text-white' : 'text-teal-300'}`}>
                {s.label}
              </span>
              {i < 3 && <div className={`w-8 sm:w-16 h-1 mx-2 ${step > s.num ? 'bg-white' : 'bg-teal-600'}`} />}
            </div>
          ))}
        </div>
      </div>

      {/* PASO 1: Selección de genomas */}
      {step === 1 && (
        <div className="space-y-6">
          {/* Cepas relacionadas predefinidas */}
          <div className="bg-white rounded-xl border border-slate-200 p-5">
            <div className="flex justify-between items-center mb-4">
              <div>
                <h3 className="font-bold text-slate-800 text-lg">Cepas de E. coli Emparentadas</h3>
                <p className="text-sm text-slate-500">Seleccione múltiples cepas para comparación</p>
              </div>
              <button
                onClick={selectAllRelated}
                className="px-4 py-2 bg-teal-100 text-teal-700 rounded-lg text-sm font-medium hover:bg-teal-200"
              >
                Seleccionar Todas
              </button>
            </div>
            
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
              {relatedStrains.map(strain => {
                const isSelected = selectedGenomes.some(g => g.accession === strain.accession)
                const isDownloaded = downloadedGenomes.some(d => d.accession === strain.accession)
                
                return (
                  <button
                    key={strain.accession}
                    onClick={() => toggleGenomeSelection({
                      accession: strain.accession,
                      organism_name: strain.organism,
                      strain: strain.strain
                    })}
                    className={`p-4 rounded-xl border-2 text-left transition-all ${
                      isSelected 
                        ? 'border-teal-500 bg-teal-50' 
                        : 'border-slate-200 hover:border-teal-200 hover:bg-slate-50'
                    }`}
                  >
                    <div className="flex justify-between items-start">
                      <div>
                        <span className={`inline-block w-5 h-5 rounded border-2 mr-2 ${
                          isSelected ? 'bg-teal-500 border-teal-500' : 'border-slate-300'
                        }`}>
                          {isSelected && <span className="text-white text-xs flex items-center justify-center">✓</span>}
                        </span>
                        <span className="font-mono text-sm text-teal-700">{strain.strain}</span>
                      </div>
                      <div className="flex gap-1">
                        {isDownloaded && (
                          <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded">Descargado</span>
                        )}
                        <span className={`text-xs px-2 py-0.5 rounded ${
                          strain.category === 'pathogenic' ? 'bg-red-100 text-red-700' :
                          strain.category === 'laboratory' ? 'bg-blue-100 text-blue-700' :
                          strain.category === 'industrial' ? 'bg-amber-100 text-amber-700' :
                          'bg-green-100 text-green-700'
                        }`}>
                          {strain.category}
                        </span>
                      </div>
                    </div>
                    <p className="text-sm text-slate-700 mt-2 font-medium">{strain.organism}</p>
                    <p className="text-xs text-slate-500 mt-1">{strain.size} • {strain.accession}</p>
                  </button>
                )
              })}
            </div>
          </div>

          {/* Búsqueda personalizada */}
          <div className="bg-white rounded-xl border border-slate-200 p-5">
            <h3 className="font-bold text-slate-800 mb-4">🔍 Buscar Otros Genomas</h3>
            <div className="flex gap-3">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={(e) => e.key === 'Enter' && searchGenomes()}
                placeholder="Buscar por nombre o accession (ej: Bacillus, Salmonella, GCF_000009045.1)"
                className="flex-1 px-4 py-3 border border-slate-200 rounded-xl focus:ring-2 focus:ring-teal-500"
              />
              <button
                onClick={searchGenomes}
                disabled={isSearching}
                className="px-6 py-3 bg-teal-600 text-white rounded-xl font-medium hover:bg-teal-700 disabled:opacity-50"
              >
                {isSearching ? 'Buscando...' : 'Buscar'}
              </button>
            </div>

            {/* Resultados de búsqueda */}
            {searchResults.length > 0 && (
              <div className="mt-4 max-h-64 overflow-y-auto border border-slate-100 rounded-xl">
                {searchResults.map(result => {
                  const isSelected = selectedGenomes.some(g => g.accession === result.accession)
                  return (
                    <button
                      key={result.accession}
                      onClick={() => toggleGenomeSelection(result)}
                      className={`w-full p-3 text-left border-b border-slate-100 last:border-0 hover:bg-slate-50 ${
                        isSelected ? 'bg-teal-50' : ''
                      }`}
                    >
                      <div className="flex justify-between">
                        <span className="font-medium text-slate-800">{result.organism_name}</span>
                        <span className="font-mono text-sm text-teal-600">{result.accession}</span>
                      </div>
                      <div className="text-xs text-slate-500 mt-1">
                        {result.genome_size_mb} Mb • {result.gene_count?.toLocaleString()} genes • GC: {result.gc_percent}%
                      </div>
                    </button>
                  )
                })}
              </div>
            )}
          </div>

          {/* Resumen de selección */}
          <div className="bg-slate-50 rounded-xl p-5">
            <div className="flex justify-between items-center">
              <div>
                <h3 className="font-bold text-slate-800">
                  {selectedGenomes.length} genomas seleccionados
                </h3>
                <p className="text-sm text-slate-500">
                  {selectedGenomes.map(g => g.strain || g.accession).join(', ') || 'Ninguno seleccionado'}
                </p>
              </div>
              <button
                onClick={downloadAllGenomes}
                disabled={selectedGenomes.length === 0}
                className="px-8 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-bold hover:shadow-lg disabled:opacity-50 transition-all"
              >
                Descargar y Continuar →
              </button>
            </div>
          </div>
        </div>
      )}

      {/* PASO 2: Progreso de descarga */}
      {step === 2 && (
        <div className="bg-white rounded-xl border border-slate-200 p-6">
          <h3 className="font-bold text-slate-800 text-lg mb-4">📥 Descargando Genomas...</h3>
          
          <div className="space-y-3">
            {selectedGenomes.map(genome => {
              const progress = downloadProgress[genome.accession] || { status: 'pending', message: 'Pendiente' }
              return (
                <div key={genome.accession} className="flex items-center gap-4 p-3 bg-slate-50 rounded-lg">
                  <div className={`w-8 h-8 rounded-full flex items-center justify-center ${
                    progress.status === 'completed' ? 'bg-green-500 text-white' :
                    progress.status === 'downloading' ? 'bg-blue-500 text-white animate-pulse' :
                    progress.status === 'error' ? 'bg-red-500 text-white' :
                    'bg-slate-300 text-slate-600'
                  }`}>
                    {progress.status === 'completed' ? '✓' :
                     progress.status === 'downloading' ? '⏳' :
                     progress.status === 'error' ? '✗' : '○'}
                  </div>
                  <div className="flex-1">
                    <p className="font-medium text-slate-800">{genome.organism_name || genome.accession}</p>
                    <p className="text-sm text-slate-500">{genome.accession}</p>
                  </div>
                  <span className={`text-sm font-medium ${
                    progress.status === 'completed' ? 'text-green-600' :
                    progress.status === 'error' ? 'text-red-600' :
                    'text-blue-600'
                  }`}>
                    {progress.message}
                  </span>
                </div>
              )
            })}
          </div>

          {!isLoading && (
            <div className="mt-6 flex justify-end">
              <button
                onClick={runFullComparison}
                className="px-8 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-bold hover:shadow-lg"
              >
                Analizar y Comparar →
              </button>
            </div>
          )}
        </div>
      )}

      {/* PASO 3: Analizando */}
      {step === 3 && (
        <div className="bg-white rounded-xl border border-slate-200 p-12 text-center">
          <div className="animate-spin w-16 h-16 border-4 border-teal-500 border-t-transparent rounded-full mx-auto mb-6"></div>
          <h3 className="text-xl font-bold text-slate-800 mb-2">Analizando Genomas...</h3>
          <p className="text-slate-500">Comparando secuencias, genes y estadísticas</p>
        </div>
      )}

      {/* PASO 4: Resultados */}
      {step === 4 && comparisonResult && (
        <div className="space-y-6">
          {/* Tabs de resultados */}
          <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
            <div className="flex border-b border-slate-200">
              {[
                { id: 'summary', name: '📊 Resumen General', icon: '📊' },
                { id: 'comparison', name: '📈 Comparación', icon: '📈' },
                { id: 'genes', name: '🧬 Genes Extremos', icon: '🧬' },
                { id: 'groups', name: '🏷️ Grupos Funcionales', icon: '🏷️' }
              ].map(tab => (
                <button
                  key={tab.id}
                  onClick={() => setActiveTab(tab.id)}
                  className={`flex-1 px-4 py-3 text-sm font-medium transition-all ${
                    activeTab === tab.id
                      ? 'bg-teal-50 text-teal-700 border-b-2 border-teal-600'
                      : 'text-slate-600 hover:bg-slate-50'
                  }`}
                >
                  {tab.name}
                </button>
              ))}
            </div>

            <div className="p-6">
              {/* TAB: Resumen General */}
              {activeTab === 'summary' && (
                <div className="space-y-6">
                  {/* Stats cards */}
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div className="bg-gradient-to-br from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
                      <p className="text-teal-100 text-sm">Genomas Comparados</p>
                      <p className="text-4xl font-bold mt-1">{comparisonResult.total_genomes_compared}</p>
                    </div>
                    <div className="bg-white border border-slate-200 rounded-xl p-5">
                      <p className="text-xs text-slate-500 uppercase">Genoma Más Grande</p>
                      <p className="text-lg font-bold text-red-600 mt-1">{comparisonResult.extremes?.largest_genome}</p>
                    </div>
                    <div className="bg-white border border-slate-200 rounded-xl p-5">
                      <p className="text-xs text-slate-500 uppercase">Mayor %GC</p>
                      <p className="text-lg font-bold text-emerald-600 mt-1">{comparisonResult.extremes?.highest_gc}</p>
                    </div>
                    <div className="bg-white border border-slate-200 rounded-xl p-5">
                      <p className="text-xs text-slate-500 uppercase">Mayor Densidad</p>
                      <p className="text-lg font-bold text-blue-600 mt-1">{comparisonResult.extremes?.highest_gene_density}</p>
                    </div>
                  </div>

                  {/* Tabla comparativa */}
                  <div className="overflow-x-auto">
                    <table className="w-full">
                      <thead className="bg-slate-50">
                        <tr>
                          <th className="px-4 py-3 text-left text-xs font-bold text-slate-600 uppercase">Organismo</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Tamaño (Mb)</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Genes</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">GC%</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Densidad</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Gen Más Largo</th>
                          <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Gen Más Corto</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-100">
                        {comparisonResult.genomes?.map((g, i) => (
                          <tr key={g.accession} className="hover:bg-slate-50">
                            <td className="px-4 py-3">
                              <div>
                                <p className="font-medium text-slate-800">{g.organism_name}</p>
                                <p className="text-xs text-teal-600 font-mono">{g.accession}</p>
                              </div>
                            </td>
                            <td className="px-4 py-3 text-right text-slate-600">{(g.genome_length / 1000000).toFixed(2)}</td>
                            <td className="px-4 py-3 text-right text-slate-600">{g.gene_count?.toLocaleString()}</td>
                            <td className="px-4 py-3 text-right font-medium text-emerald-700">{g.gc_content?.toFixed(1)}%</td>
                            <td className="px-4 py-3 text-right text-slate-600">{g.gene_density?.toFixed(1)}</td>
                            <td className="px-4 py-3 text-right text-red-600 font-medium">{g.max_gene_length?.toLocaleString()} bp</td>
                            <td className="px-4 py-3 text-right text-green-600 font-medium">{g.min_gene_length?.toLocaleString()} bp</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              {/* TAB: Comparación con gráficos */}
              {activeTab === 'comparison' && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                  {/* Tamaño del genoma */}
                  <div className="bg-slate-50 rounded-xl p-4">
                    <h4 className="font-bold text-slate-800 mb-4">Tamaño del Genoma (Mb)</h4>
                    <ResponsiveContainer width="100%" height={280}>
                      <BarChart data={chartData} layout="vertical">
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis type="number" />
                        <YAxis dataKey="name" type="category" width={100} tick={{ fontSize: 10 }} />
                        <Tooltip />
                        <Bar dataKey="size" fill="#0d9488" radius={[0, 4, 4, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>

                  {/* Número de genes */}
                  <div className="bg-slate-50 rounded-xl p-4">
                    <h4 className="font-bold text-slate-800 mb-4">Número de Genes</h4>
                    <ResponsiveContainer width="100%" height={280}>
                      <BarChart data={chartData} layout="vertical">
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis type="number" />
                        <YAxis dataKey="name" type="category" width={100} tick={{ fontSize: 10 }} />
                        <Tooltip />
                        <Bar dataKey="genes" fill="#10b981" radius={[0, 4, 4, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>

                  {/* Contenido GC */}
                  <div className="bg-slate-50 rounded-xl p-4">
                    <h4 className="font-bold text-slate-800 mb-4">Contenido GC (%)</h4>
                    <ResponsiveContainer width="100%" height={280}>
                      <BarChart data={chartData}>
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="name" tick={{ fontSize: 10 }} />
                        <YAxis domain={[45, 55]} />
                        <Tooltip />
                        <Bar dataKey="gc" fill="#6366f1" radius={[4, 4, 0, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>

                  {/* Densidad génica */}
                  <div className="bg-slate-50 rounded-xl p-4">
                    <h4 className="font-bold text-slate-800 mb-4">Densidad Génica (genes/Mb)</h4>
                    <ResponsiveContainer width="100%" height={280}>
                      <BarChart data={chartData}>
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="name" tick={{ fontSize: 10 }} />
                        <YAxis />
                        <Tooltip />
                        <Bar dataKey="density" fill="#f59e0b" radius={[4, 4, 0, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                </div>
              )}

              {/* TAB: Genes Extremos */}
              {activeTab === 'genes' && (
                <div className="space-y-6">
                  {/* Filtros */}
                  <div className="flex flex-wrap gap-4 p-4 bg-slate-50 rounded-xl">
                    <div>
                      <label className="block text-sm font-medium text-slate-700 mb-1">Filtrar por</label>
                      <select
                        value={geneFilter}
                        onChange={(e) => setGeneFilter(e.target.value)}
                        className="px-3 py-2 border border-slate-200 rounded-lg"
                      >
                        <option value="all">Todos los extremos</option>
                        <option value="size">Por tamaño</option>
                        <option value="group">Por grupo funcional</option>
                      </select>
                    </div>
                    
                    {geneFilter === 'size' && (
                      <>
                        <div>
                          <label className="block text-sm font-medium text-slate-700 mb-1">Orden</label>
                          <select
                            value={geneSizeOrder}
                            onChange={(e) => setGeneSizeOrder(e.target.value)}
                            className="px-3 py-2 border border-slate-200 rounded-lg"
                          >
                            <option value="largest">Más largos</option>
                            <option value="smallest">Más cortos</option>
                          </select>
                        </div>
                        <div>
                          <label className="block text-sm font-medium text-slate-700 mb-1">Cantidad</label>
                          <select
                            value={geneCount}
                            onChange={(e) => setGeneCount(parseInt(e.target.value))}
                            className="px-3 py-2 border border-slate-200 rounded-lg"
                          >
                            <option value={10}>Top 10</option>
                            <option value={20}>Top 20</option>
                            <option value={50}>Top 50</option>
                            <option value={100}>Top 100</option>
                          </select>
                        </div>
                      </>
                    )}
                    
                    {geneFilter === 'group' && (
                      <div>
                        <label className="block text-sm font-medium text-slate-700 mb-1">Grupo</label>
                        <select
                          value={selectedGroup}
                          onChange={(e) => setSelectedGroup(e.target.value)}
                          className="px-3 py-2 border border-slate-200 rounded-lg"
                        >
                          <option value="">Seleccionar...</option>
                          {functionalGroups.map(g => (
                            <option key={g.id} value={g.id}>{g.name}</option>
                          ))}
                        </select>
                      </div>
                    )}
                  </div>

                  {/* Genes más largos y más cortos (globales) */}
                  {geneFilter === 'all' && (
                    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                      {/* Más largos */}
                      <div className="border border-slate-200 rounded-xl overflow-hidden">
                        <div className="bg-gradient-to-r from-red-50 to-orange-50 p-4 border-b">
                          <h4 className="font-bold text-slate-800">📏 Top 10 Genes Más Largos (Global)</h4>
                        </div>
                        <div className="divide-y divide-slate-100 max-h-96 overflow-y-auto">
                          {comparisonResult.longest_genes_global?.map((gene, i) => (
                            <div key={i} className="p-3 hover:bg-slate-50">
                              <div className="flex justify-between">
                                <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                                <span className="font-bold text-red-600">{gene.length?.toLocaleString()} bp</span>
                              </div>
                              <p className="text-xs text-slate-500 mt-1 line-clamp-1">{gene.product}</p>
                              <p className="text-xs text-slate-400">{gene.genome}</p>
                            </div>
                          ))}
                        </div>
                      </div>

                      {/* Más cortos */}
                      <div className="border border-slate-200 rounded-xl overflow-hidden">
                        <div className="bg-gradient-to-r from-green-50 to-teal-50 p-4 border-b">
                          <h4 className="font-bold text-slate-800">🔬 Top 10 Genes Más Cortos (Global)</h4>
                        </div>
                        <div className="divide-y divide-slate-100 max-h-96 overflow-y-auto">
                          {comparisonResult.shortest_genes_global?.map((gene, i) => (
                            <div key={i} className="p-3 hover:bg-slate-50">
                              <div className="flex justify-between">
                                <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                                <span className="font-bold text-green-600">{gene.length?.toLocaleString()} bp</span>
                              </div>
                              <p className="text-xs text-slate-500 mt-1 line-clamp-1">{gene.product}</p>
                              <p className="text-xs text-slate-400">{gene.genome}</p>
                            </div>
                          ))}
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Tabla filtrada */}
                  {(geneFilter === 'size' || geneFilter === 'group') && comparisonResult.filteredGenes && (
                    <div className="ag-theme-alpine" style={{ height: 400 }}>
                      <AgGridReact
                        rowData={comparisonResult.filteredGenes}
                        columnDefs={geneColumns}
                        defaultColDef={{ resizable: true }}
                        animateRows={true}
                      />
                    </div>
                  )}
                </div>
              )}

              {/* TAB: Grupos Funcionales */}
              {activeTab === 'groups' && (
                <div className="space-y-6">
                  <p className="text-slate-600">
                    Clasificación de genes por función basada en el producto génico anotado.
                  </p>
                  
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    {comparisonResult.groupsSummary?.groups?.map(group => (
                      <button
                        key={group.id}
                        onClick={() => {
                          setActiveTab('genes')
                          setGeneFilter('group')
                          setSelectedGroup(group.id)
                        }}
                        className="p-4 bg-white border border-slate-200 rounded-xl hover:border-teal-300 hover:bg-teal-50 transition-all text-left"
                      >
                        <div className="flex justify-between items-start mb-2">
                          <span className="text-2xl">
                            {group.id === 'metabolism' ? '⚗️' :
                             group.id === 'transport' ? '🚚' :
                             group.id === 'regulation' ? '🎛️' :
                             group.id === 'dna_rna' ? '🧬' :
                             group.id === 'protein' ? '🔬' :
                             group.id === 'cell_structure' ? '🏗️' :
                             group.id === 'stress_response' ? '🛡️' : '❓'}
                          </span>
                          <span className="text-xs bg-teal-100 text-teal-700 px-2 py-0.5 rounded-full font-bold">
                            {group.percentage?.toFixed(1)}%
                          </span>
                        </div>
                        <h4 className="font-bold text-slate-800">{group.name}</h4>
                        <p className="text-sm text-teal-600 font-medium">{group.count?.toLocaleString()} genes</p>
                        <p className="text-xs text-slate-500 mt-1">~{group.avg_length?.toFixed(0)} bp promedio</p>
                      </button>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Botón para reiniciar */}
          <div className="flex justify-center">
            <button
              onClick={() => {
                setStep(1)
                setSelectedGenomes([])
                setComparisonResult(null)
                setDownloadProgress({})
              }}
              className="px-6 py-3 border border-slate-300 text-slate-600 rounded-xl font-medium hover:bg-slate-50"
            >
              ← Iniciar Nueva Comparación
            </button>
          </div>
        </div>
      )}
    </div>
  )
}
