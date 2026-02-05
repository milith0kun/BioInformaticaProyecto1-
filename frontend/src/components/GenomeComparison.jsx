/**
 * GenomeComparison Component
 * Compara múltiples genomas y muestra resultados de análisis comparativo
 */
import { useState, useEffect, useMemo } from 'react'
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Legend,
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar
} from 'recharts'
import { api } from '../services/api'
import toast from 'react-hot-toast'

export default function GenomeComparison({ hasAnalysis, downloadedGenomes = [] }) {
  const [comparisonData, setComparisonData] = useState(null)
  const [isComparing, setIsComparing] = useState(false)
  const [relatedStrains, setRelatedStrains] = useState([])
  const [selectedStrains, setSelectedStrains] = useState([])
  const [categoryFilter, setCategoryFilter] = useState('')
  const [activeTab, setActiveTab] = useState('overview')

  // Cargar cepas relacionadas al montar
  useEffect(() => {
    loadRelatedStrains()
  }, [])

  const loadRelatedStrains = async () => {
    try {
      const result = await api.getRelatedStrains()
      setRelatedStrains(result.strains || [])
    } catch (error) {
      console.error('Error loading related strains:', error)
    }
  }

  const runComparison = async (accessions = null) => {
    setIsComparing(true)
    try {
      toast.loading('Comparando genomas...', { id: 'compare' })
      const result = await api.compareGenomes(accessions)
      setComparisonData(result)
      toast.success(`${result.total_genomes_compared} genomas comparados`, { id: 'compare' })
    } catch (error) {
      toast.error('Error al comparar: ' + (error.response?.data?.detail || error.message), { id: 'compare' })
    } finally {
      setIsComparing(false)
    }
  }

  const filteredStrains = useMemo(() => {
    if (!categoryFilter) return relatedStrains
    return relatedStrains.filter(s => s.category === categoryFilter)
  }, [relatedStrains, categoryFilter])

  // Datos para gráficos comparativos
  const chartData = useMemo(() => {
    if (!comparisonData?.genomes) return []
    return comparisonData.genomes.map(g => ({
      name: g.accession.replace('GCF_', '').replace('GCA_', '').slice(0, 8),
      fullName: g.organism_name,
      genome_size: (g.genome_length / 1000000).toFixed(2),
      genes: g.gene_count,
      gc: g.gc_content,
      density: g.gene_density,
      avg_length: g.avg_gene_length
    }))
  }, [comparisonData])

  const radarData = useMemo(() => {
    if (!comparisonData?.genomes || comparisonData.genomes.length === 0) return []
    
    const metrics = ['genome_length', 'gene_count', 'gc_content', 'gene_density', 'avg_gene_length']
    const maxValues = {}
    
    metrics.forEach(m => {
      maxValues[m] = Math.max(...comparisonData.genomes.map(g => g[m]))
    })
    
    return [
      { metric: 'Tamaño', ...Object.fromEntries(
        comparisonData.genomes.map(g => [g.accession.slice(-8), (g.genome_length / maxValues.genome_length) * 100])
      )},
      { metric: 'Genes', ...Object.fromEntries(
        comparisonData.genomes.map(g => [g.accession.slice(-8), (g.gene_count / maxValues.gene_count) * 100])
      )},
      { metric: 'GC%', ...Object.fromEntries(
        comparisonData.genomes.map(g => [g.accession.slice(-8), (g.gc_content / maxValues.gc_content) * 100])
      )},
      { metric: 'Densidad', ...Object.fromEntries(
        comparisonData.genomes.map(g => [g.accession.slice(-8), (g.gene_density / maxValues.gene_density) * 100])
      )},
      { metric: 'Long. Media', ...Object.fromEntries(
        comparisonData.genomes.map(g => [g.accession.slice(-8), (g.avg_gene_length / maxValues.avg_gene_length) * 100])
      )}
    ]
  }, [comparisonData])

  const COLORS = ['#0d9488', '#10b981', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6']

  if (!hasAnalysis && downloadedGenomes.length === 0) {
    return (
      <div className="text-center py-20">
        <div className="w-20 h-20 mx-auto bg-slate-100 rounded-2xl flex items-center justify-center mb-6">
          <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
          </svg>
        </div>
        <h2 className="text-2xl font-bold text-slate-800 mb-2">Comparación Genómica</h2>
        <p className="text-slate-500 max-w-md mx-auto">
          Descargue y analice genomas para poder compararlos entre sí. Podrá ver diferencias en tamaño, contenido GC, densidad génica y más.
        </p>
      </div>
    )
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4">
        <div>
          <h2 className="text-2xl font-bold text-slate-800">Comparación Genómica</h2>
          <p className="text-slate-500">Analice y compare múltiples genomas secuencialmente</p>
        </div>
        <button
          onClick={() => runComparison()}
          disabled={isComparing}
          className="px-6 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-medium hover:shadow-lg transition-all disabled:opacity-50 flex items-center gap-2"
        >
          {isComparing ? (
            <>
              <svg className="animate-spin h-5 w-5" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
              </svg>
              Comparando...
            </>
          ) : (
            <>
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
              </svg>
              Comparar Genomas Descargados
            </>
          )}
        </button>
      </div>

      {/* Cepas Relacionadas */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-3 mb-4">
          <div>
            <h3 className="font-semibold text-slate-800">Cepas de E. coli Emparentadas</h3>
            <p className="text-sm text-slate-500">Descargue cepas relacionadas para análisis comparativo</p>
          </div>
          <select
            value={categoryFilter}
            onChange={(e) => setCategoryFilter(e.target.value)}
            className="px-3 py-2 border border-slate-200 rounded-lg text-sm focus:ring-2 focus:ring-teal-500"
          >
            <option value="">Todas las categorías</option>
            <option value="laboratory">🧪 Laboratorio</option>
            <option value="pathogenic">⚠️ Patógenas</option>
            <option value="industrial">🏭 Industriales</option>
            <option value="probiotic">💊 Probióticas</option>
          </select>
        </div>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
          {filteredStrains.map(strain => (
            <div
              key={strain.accession}
              className="p-4 border border-slate-100 rounded-xl hover:border-teal-200 hover:bg-teal-50/30 transition-all"
            >
              <div className="flex justify-between items-start mb-2">
                <span className="font-mono text-sm text-teal-700 bg-teal-50 px-2 py-0.5 rounded">
                  {strain.strain}
                </span>
                <span className={`text-xs px-2 py-0.5 rounded-full ${
                  strain.category === 'pathogenic' ? 'bg-red-100 text-red-700' :
                  strain.category === 'laboratory' ? 'bg-blue-100 text-blue-700' :
                  strain.category === 'industrial' ? 'bg-amber-100 text-amber-700' :
                  'bg-green-100 text-green-700'
                }`}>
                  {strain.category === 'laboratory' ? '🧪 Lab' :
                   strain.category === 'pathogenic' ? '⚠️ Patógena' :
                   strain.category === 'industrial' ? '🏭 Industrial' : '💊 Probiótica'}
                </span>
              </div>
              <p className="text-xs text-slate-600 mb-2 line-clamp-2">{strain.description}</p>
              <div className="flex justify-between items-center text-xs text-slate-500">
                <span>{strain.genome_size_mb} Mb</span>
                <span className="font-mono">{strain.accession}</span>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Tabs de resultados */}
      {comparisonData && (
        <>
          <div className="flex gap-2 border-b border-slate-200">
            {[
              { id: 'overview', name: 'Resumen', icon: '📊' },
              { id: 'charts', name: 'Gráficos', icon: '📈' },
              { id: 'genes', name: 'Genes Extremos', icon: '🧬' }
            ].map(tab => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`px-4 py-3 text-sm font-medium border-b-2 transition-all ${
                  activeTab === tab.id
                    ? 'border-teal-600 text-teal-600'
                    : 'border-transparent text-slate-500 hover:text-slate-700'
                }`}
              >
                {tab.icon} {tab.name}
              </button>
            ))}
          </div>

          {/* Tab: Resumen */}
          {activeTab === 'overview' && (
            <div className="space-y-6">
              {/* Estadísticas globales */}
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <div className="bg-gradient-to-br from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
                  <p className="text-teal-100 text-sm">Genomas Comparados</p>
                  <p className="text-3xl font-bold mt-1">{comparisonData.total_genomes_compared}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                  <p className="text-xs text-slate-500 uppercase">Mayor Genoma</p>
                  <p className="text-lg font-bold text-slate-800 mt-1">{comparisonData.extremes.largest_genome}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                  <p className="text-xs text-slate-500 uppercase">Mayor GC%</p>
                  <p className="text-lg font-bold text-emerald-700 mt-1">{comparisonData.extremes.highest_gc}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                  <p className="text-xs text-slate-500 uppercase">Mayor Densidad</p>
                  <p className="text-lg font-bold text-teal-700 mt-1">{comparisonData.extremes.highest_gene_density}</p>
                </div>
              </div>

              {/* Tabla comparativa */}
              <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                <div className="p-4 border-b border-slate-100">
                  <h3 className="font-semibold text-slate-800">Comparación Detallada</h3>
                </div>
                <div className="overflow-x-auto">
                  <table className="w-full">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-4 py-3 text-left text-xs font-semibold text-slate-600 uppercase">Accession</th>
                        <th className="px-4 py-3 text-left text-xs font-semibold text-slate-600 uppercase">Organismo</th>
                        <th className="px-4 py-3 text-right text-xs font-semibold text-slate-600 uppercase">Tamaño (Mb)</th>
                        <th className="px-4 py-3 text-right text-xs font-semibold text-slate-600 uppercase">Genes</th>
                        <th className="px-4 py-3 text-right text-xs font-semibold text-slate-600 uppercase">GC%</th>
                        <th className="px-4 py-3 text-right text-xs font-semibold text-slate-600 uppercase">Densidad</th>
                        <th className="px-4 py-3 text-right text-xs font-semibold text-slate-600 uppercase">Long. Media</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {comparisonData.genomes.map((g, i) => (
                        <tr key={g.accession} className="hover:bg-slate-50">
                          <td className="px-4 py-3 font-mono text-sm text-teal-700">{g.accession}</td>
                          <td className="px-4 py-3 text-sm text-slate-700 max-w-xs truncate">{g.organism_name}</td>
                          <td className="px-4 py-3 text-sm text-right text-slate-600">{(g.genome_length / 1000000).toFixed(2)}</td>
                          <td className="px-4 py-3 text-sm text-right text-slate-600">{g.gene_count.toLocaleString()}</td>
                          <td className="px-4 py-3 text-sm text-right text-emerald-700 font-medium">{g.gc_content.toFixed(1)}%</td>
                          <td className="px-4 py-3 text-sm text-right text-slate-600">{g.gene_density.toFixed(1)}</td>
                          <td className="px-4 py-3 text-sm text-right text-slate-600">{g.avg_gene_length.toFixed(0)} bp</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>
          )}

          {/* Tab: Gráficos */}
          {activeTab === 'charts' && (
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              {/* Gráfico de barras - Tamaño */}
              <div className="bg-white rounded-xl border border-slate-200 p-5">
                <h3 className="font-semibold text-slate-800 mb-4">Tamaño del Genoma (Mb)</h3>
                <ResponsiveContainer width="100%" height={280}>
                  <BarChart data={chartData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                    <YAxis tick={{ fontSize: 11 }} />
                    <Tooltip 
                      contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                      labelFormatter={(label, payload) => payload[0]?.payload?.fullName || label}
                    />
                    <Bar dataKey="genome_size" fill="#0d9488" radius={[4, 4, 0, 0]} name="Mb" />
                  </BarChart>
                </ResponsiveContainer>
              </div>

              {/* Gráfico de barras - Genes */}
              <div className="bg-white rounded-xl border border-slate-200 p-5">
                <h3 className="font-semibold text-slate-800 mb-4">Número de Genes</h3>
                <ResponsiveContainer width="100%" height={280}>
                  <BarChart data={chartData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                    <YAxis tick={{ fontSize: 11 }} />
                    <Tooltip 
                      contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                      labelFormatter={(label, payload) => payload[0]?.payload?.fullName || label}
                    />
                    <Bar dataKey="genes" fill="#10b981" radius={[4, 4, 0, 0]} name="Genes" />
                  </BarChart>
                </ResponsiveContainer>
              </div>

              {/* Gráfico de barras - GC */}
              <div className="bg-white rounded-xl border border-slate-200 p-5">
                <h3 className="font-semibold text-slate-800 mb-4">Contenido GC (%)</h3>
                <ResponsiveContainer width="100%" height={280}>
                  <BarChart data={chartData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                    <YAxis domain={[40, 60]} tick={{ fontSize: 11 }} />
                    <Tooltip 
                      contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                      labelFormatter={(label, payload) => payload[0]?.payload?.fullName || label}
                    />
                    <Bar dataKey="gc" fill="#6366f1" radius={[4, 4, 0, 0]} name="GC%" />
                  </BarChart>
                </ResponsiveContainer>
              </div>

              {/* Gráfico de barras - Densidad */}
              <div className="bg-white rounded-xl border border-slate-200 p-5">
                <h3 className="font-semibold text-slate-800 mb-4">Densidad Génica (genes/Mb)</h3>
                <ResponsiveContainer width="100%" height={280}>
                  <BarChart data={chartData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                    <YAxis tick={{ fontSize: 11 }} />
                    <Tooltip 
                      contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                      labelFormatter={(label, payload) => payload[0]?.payload?.fullName || label}
                    />
                    <Bar dataKey="density" fill="#f59e0b" radius={[4, 4, 0, 0]} name="genes/Mb" />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* Tab: Genes Extremos */}
          {activeTab === 'genes' && (
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              {/* Genes más largos */}
              <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                <div className="p-4 border-b border-slate-100 bg-gradient-to-r from-red-50 to-orange-50">
                  <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                    <span className="text-xl">📏</span> Genes Más Largos (Global)
                  </h3>
                </div>
                <div className="divide-y divide-slate-100">
                  {comparisonData.longest_genes_global.map((gene, i) => (
                    <div key={i} className="p-3 hover:bg-slate-50">
                      <div className="flex justify-between items-start">
                        <div>
                          <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                          <p className="text-xs text-slate-500 mt-0.5 line-clamp-1">{gene.product || 'Sin descripción'}</p>
                        </div>
                        <div className="text-right">
                          <span className="font-bold text-red-600">{gene.length.toLocaleString()} bp</span>
                          <p className="text-xs text-slate-400">{gene.genome}</p>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>

              {/* Genes más cortos */}
              <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                <div className="p-4 border-b border-slate-100 bg-gradient-to-r from-green-50 to-teal-50">
                  <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                    <span className="text-xl">🔬</span> Genes Más Cortos (Global)
                  </h3>
                </div>
                <div className="divide-y divide-slate-100">
                  {comparisonData.shortest_genes_global.map((gene, i) => (
                    <div key={i} className="p-3 hover:bg-slate-50">
                      <div className="flex justify-between items-start">
                        <div>
                          <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                          <p className="text-xs text-slate-500 mt-0.5 line-clamp-1">{gene.product || 'Sin descripción'}</p>
                        </div>
                        <div className="text-right">
                          <span className="font-bold text-green-600">{gene.length.toLocaleString()} bp</span>
                          <p className="text-xs text-slate-400">{gene.genome}</p>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          )}
        </>
      )}

      {/* Info Box */}
      <div className="bg-teal-50 border border-teal-100 rounded-xl p-5">
        <div className="flex gap-3">
          <div className="w-10 h-10 bg-teal-500 rounded-xl flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
          </div>
          <div>
            <h4 className="font-semibold text-teal-900 mb-1">Comparación Secuencial</h4>
            <p className="text-sm text-teal-700">
              El sistema analiza los genomas descargados secuencialmente, comparando métricas clave como tamaño del genoma, 
              número de genes, contenido GC y densidad génica. Los genes de mayor y menor tamaño se identifican globalmente 
              para todos los genomas comparados.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
