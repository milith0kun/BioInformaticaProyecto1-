/**
 * GeneStatistics Component
 * Displays gene analysis with interactive charts and searchable table
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
  ScatterChart,
  Scatter,
  ZAxis
} from 'recharts'
import { AgGridReact } from 'ag-grid-react'
import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'
import { api } from '../services/api'

export default function GeneStatistics({ geneData }) {
  const [searchQuery, setSearchQuery] = useState('')
  const [paginatedGenes, setPaginatedGenes] = useState({ genes: [], total: 0 })
  const [currentPage, setCurrentPage] = useState(1)
  const [isLoading, setIsLoading] = useState(false)

  const loadGenes = useCallback(async (page = 1, search = '') => {
    setIsLoading(true)
    try {
      const result = await api.getGeneResults(page, 50, search)
      setPaginatedGenes(result)
      setCurrentPage(page)
    } catch (error) {
      console.error('Error loading genes:', error)
    } finally {
      setIsLoading(false)
    }
  }, [])

  useEffect(() => {
    if (geneData) loadGenes(1, '')
  }, [geneData, loadGenes])

  const handleSearch = () => loadGenes(1, searchQuery)
  const handleKeyPress = (e) => { if (e.key === 'Enter') handleSearch() }

  const columnDefs = useMemo(() => [
    { field: 'locus_tag', headerName: 'Locus Tag', width: 120, sortable: true, filter: true },
    { field: 'start', headerName: 'Inicio', width: 100, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'end', headerName: 'Fin', width: 100, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'length', headerName: 'Longitud', width: 100, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'strand', headerName: 'Hebra', width: 70, valueFormatter: p => p.value === 1 ? '+' : '-' },
    { field: 'gc_content', headerName: 'GC%', width: 70, sortable: true, valueFormatter: p => p.value?.toFixed(1) },
    { field: 'product', headerName: 'Producto', flex: 1, tooltipField: 'product' },
  ], [])

  if (!geneData) {
    return (
      <div className="text-center py-20">
        <div className="w-20 h-20 mx-auto bg-slate-100 rounded-2xl flex items-center justify-center mb-6">
          <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
          </svg>
        </div>
        <h2 className="text-xl font-bold text-slate-700 mb-2">Sin datos de genes</h2>
        <p className="text-slate-500">Ejecute el análisis completo para ver las estadísticas</p>
      </div>
    )
  }

  const sizeHistogramData = useMemo(() => {
    if (!paginatedGenes.genes.length) return []
    const bins = [
      { range: '0-300', min: 0, max: 300, count: 0 },
      { range: '300-600', min: 300, max: 600, count: 0 },
      { range: '600-900', min: 600, max: 900, count: 0 },
      { range: '900-1.2K', min: 900, max: 1200, count: 0 },
      { range: '1.2-1.5K', min: 1200, max: 1500, count: 0 },
      { range: '1.5-2K', min: 1500, max: 2000, count: 0 },
      { range: '2-3K', min: 2000, max: 3000, count: 0 },
      { range: '3K+', min: 3000, max: Infinity, count: 0 },
    ]
    paginatedGenes.genes.forEach(gene => {
      const bin = bins.find(b => gene.length >= b.min && gene.length < b.max)
      if (bin) bin.count++
    })
    return bins
  }, [paginatedGenes.genes])

  const gcScatterData = useMemo(() => {
    return paginatedGenes.genes.slice(0, 200).map(gene => ({
      x: gene.length, y: gene.gc_content, locus: gene.locus_tag
    }))
  }, [paginatedGenes.genes])

  return (
    <div className="space-y-6">
      {/* Stats Cards */}
      <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
        <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-5 text-white">
          <p className="text-teal-100 text-sm">Total de Genes</p>
          <p className="text-2xl font-bold mt-1">{geneData.total_genes.toLocaleString()}</p>
        </div>
        <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-5 text-white">
          <p className="text-emerald-100 text-sm">Total de CDS</p>
          <p className="text-2xl font-bold mt-1">{geneData.total_cds.toLocaleString()}</p>
        </div>
        <div className="bg-gradient-to-br from-slate-600 to-slate-700 rounded-xl p-5 text-white">
          <p className="text-slate-300 text-sm">Densidad Génica</p>
          <p className="text-2xl font-bold mt-1">{geneData.gene_density}</p>
          <p className="text-slate-400 text-xs">genes/Mb</p>
        </div>
        <div className="bg-gradient-to-br from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
          <p className="text-teal-100 text-sm">Contenido GC</p>
          <p className="text-2xl font-bold mt-1">{geneData.gc_content}%</p>
        </div>
      </div>

      {/* Size Statistics */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <h3 className="font-semibold text-slate-800 mb-4">Estadísticas de Tamaño de Genes</h3>
        <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
          {[
            { label: 'Media', value: geneData.size_statistics.mean.toFixed(0), unit: 'bp' },
            { label: 'Mediana', value: geneData.size_statistics.median.toFixed(0), unit: 'bp' },
            { label: 'Mínimo', value: geneData.size_statistics.min, unit: 'bp', color: 'text-teal-700' },
            { label: 'Máximo', value: geneData.size_statistics.max.toLocaleString(), unit: 'bp', color: 'text-red-600' },
            { label: 'Desv. Est.', value: geneData.size_statistics.std.toFixed(0), unit: 'bp' },
          ].map((stat, i) => (
            <div key={i} className="text-center p-3 bg-slate-50 rounded-lg">
              <p className="text-xs text-slate-500 uppercase">{stat.label}</p>
              <p className={`text-xl font-bold ${stat.color || 'text-slate-800'}`}>{stat.value}</p>
              <p className="text-xs text-slate-400">{stat.unit}</p>
            </div>
          ))}
        </div>
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Histogram */}
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">Distribución de Tamaños</h3>
          <ResponsiveContainer width="100%" height={260}>
            <BarChart data={sizeHistogramData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
              <XAxis dataKey="range" tick={{ fontSize: 10, fill: '#64748b' }} angle={-45} textAnchor="end" height={50} />
              <YAxis tick={{ fontSize: 12, fill: '#64748b' }} />
              <Tooltip
                formatter={(value) => [value, 'Genes']}
                contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
              />
              <Bar dataKey="count" fill="#14b8a6" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Scatter Plot */}
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">GC% vs Longitud</h3>
          <ResponsiveContainer width="100%" height={260}>
            <ScatterChart margin={{ bottom: 20 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
              <XAxis
                dataKey="x" type="number" name="Longitud"
                tick={{ fontSize: 12, fill: '#64748b' }}
                label={{ value: 'Longitud (bp)', position: 'bottom', offset: 0, fontSize: 11, fill: '#64748b' }}
              />
              <YAxis
                dataKey="y" type="number" name="GC%"
                tick={{ fontSize: 12, fill: '#64748b' }}
                domain={[20, 80]}
              />
              <ZAxis range={[20, 20]} />
              <Tooltip
                cursor={{ strokeDasharray: '3 3' }}
                contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
              />
              <Scatter data={gcScatterData} fill="#10b981" opacity={0.6} />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Gene Table */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4 mb-4">
          <h3 className="font-semibold text-slate-800">
            Lista de Genes ({paginatedGenes.total.toLocaleString()} total)
          </h3>
          <div className="flex items-center gap-2">
            <div className="relative">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={handleKeyPress}
                placeholder="Buscar por locus o producto..."
                className="pl-9 pr-4 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500 focus:border-teal-500 text-sm w-48 sm:w-64"
              />
              <svg className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
              </svg>
            </div>
            <button
              onClick={handleSearch}
              className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 transition-colors text-sm font-medium"
            >
              Buscar
            </button>
          </div>
        </div>

        <div className="ag-theme-alpine" style={{ height: 350 }}>
          <AgGridReact
            rowData={paginatedGenes.genes}
            columnDefs={columnDefs}
            defaultColDef={{ resizable: true }}
            animateRows={true}
            pagination={false}
            loading={isLoading}
            overlayLoadingTemplate='<span class="text-slate-500">Cargando genes...</span>'
            overlayNoRowsTemplate='<span class="text-slate-500">No se encontraron genes</span>'
          />
        </div>

        {paginatedGenes.total_pages > 1 && (
          <div className="flex justify-center items-center gap-2 mt-4">
            <button
              onClick={() => loadGenes(currentPage - 1, searchQuery)}
              disabled={currentPage === 1}
              className="px-3 py-1.5 border border-slate-200 rounded-lg disabled:opacity-50 hover:bg-slate-50 text-sm"
            >
              Anterior
            </button>
            <span className="text-sm text-slate-600">
              Página {currentPage} de {paginatedGenes.total_pages}
            </span>
            <button
              onClick={() => loadGenes(currentPage + 1, searchQuery)}
              disabled={currentPage === paginatedGenes.total_pages}
              className="px-3 py-1.5 border border-slate-200 rounded-lg disabled:opacity-50 hover:bg-slate-50 text-sm"
            >
              Siguiente
            </button>
          </div>
        )}
      </div>

      {/* Info Box */}
      <div className="bg-teal-50 border border-teal-100 rounded-xl p-4">
        <div className="flex gap-3">
          <div className="w-8 h-8 bg-teal-500 rounded-lg flex items-center justify-center flex-shrink-0">
            <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <p className="text-sm text-teal-800">
            E. coli K-12 MG1655 tiene aproximadamente 4,300 genes que cubren ~87% del genoma.
            El tamaño promedio es de ~950 pb. Los genes están distribuidos casi equitativamente entre ambas hebras.
          </p>
        </div>
      </div>
    </div>
  )
}
