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
  ZAxis,
  LineChart,
  Line
} from 'recharts'
import { MagnifyingGlassIcon, BeakerIcon } from '@heroicons/react/24/outline'
import { AgGridReact } from 'ag-grid-react'
import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'
import { api } from '../services/api'

export default function GeneStatistics({ geneData }) {
  const [searchQuery, setSearchQuery] = useState('')
  const [paginatedGenes, setPaginatedGenes] = useState({ genes: [], total: 0 })
  const [currentPage, setCurrentPage] = useState(1)
  const [isLoading, setIsLoading] = useState(false)

  // Load paginated genes
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
    if (geneData) {
      loadGenes(1, '')
    }
  }, [geneData, loadGenes])

  const handleSearch = () => {
    loadGenes(1, searchQuery)
  }

  const handleKeyPress = (e) => {
    if (e.key === 'Enter') {
      handleSearch()
    }
  }

  // AG Grid column definitions
  const columnDefs = useMemo(() => [
    { field: 'locus_tag', headerName: 'Locus Tag', width: 120, sortable: true, filter: true },
    { field: 'start', headerName: 'Inicio', width: 100, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'end', headerName: 'Fin', width: 100, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'length', headerName: 'Longitud (bp)', width: 120, sortable: true, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'strand', headerName: 'Hebra', width: 80, valueFormatter: p => p.value === 1 ? '+' : '-' },
    { field: 'gc_content', headerName: 'GC %', width: 80, sortable: true, valueFormatter: p => p.value?.toFixed(1) },
    { field: 'product', headerName: 'Producto', flex: 1, tooltipField: 'product' },
  ], [])

  if (!geneData) {
    return (
      <div className="text-center py-16">
        <BeakerIcon className="h-24 w-24 mx-auto text-gray-300 mb-6" />
        <h2 className="text-2xl font-bold text-gray-700 mb-2">
          Sin datos de genes
        </h2>
        <p className="text-gray-500">
          Ejecute el análisis completo para ver las estadísticas de genes
        </p>
      </div>
    )
  }

  // Prepare histogram data for gene sizes
  const sizeHistogramData = useMemo(() => {
    if (!paginatedGenes.genes.length) return []
    
    // Create size bins
    const bins = [
      { range: '0-300', min: 0, max: 300, count: 0 },
      { range: '300-600', min: 300, max: 600, count: 0 },
      { range: '600-900', min: 600, max: 900, count: 0 },
      { range: '900-1200', min: 900, max: 1200, count: 0 },
      { range: '1200-1500', min: 1200, max: 1500, count: 0 },
      { range: '1500-2000', min: 1500, max: 2000, count: 0 },
      { range: '2000-3000', min: 2000, max: 3000, count: 0 },
      { range: '3000+', min: 3000, max: Infinity, count: 0 },
    ]
    
    paginatedGenes.genes.forEach(gene => {
      const bin = bins.find(b => gene.length >= b.min && gene.length < b.max)
      if (bin) bin.count++
    })
    
    return bins
  }, [paginatedGenes.genes])

  // Prepare scatter data for GC content
  const gcScatterData = useMemo(() => {
    return paginatedGenes.genes.slice(0, 200).map(gene => ({
      x: gene.length,
      y: gene.gc_content,
      locus: gene.locus_tag
    }))
  }, [paginatedGenes.genes])

  return (
    <div className="space-y-6">
      {/* Stats Cards */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-6 text-white">
          <p className="text-emerald-100">Total de Genes</p>
          <p className="text-3xl font-bold mt-2">
            {geneData.total_genes.toLocaleString()}
          </p>
        </div>
        <div className="bg-gradient-to-br from-blue-500 to-blue-600 rounded-xl p-6 text-white">
          <p className="text-blue-100">Total de CDS</p>
          <p className="text-3xl font-bold mt-2">
            {geneData.total_cds.toLocaleString()}
          </p>
        </div>
        <div className="bg-gradient-to-br from-purple-500 to-purple-600 rounded-xl p-6 text-white">
          <p className="text-purple-100">Densidad Génica</p>
          <p className="text-3xl font-bold mt-2">
            {geneData.gene_density}
          </p>
          <p className="text-purple-100 text-sm">genes/Mb</p>
        </div>
        <div className="bg-gradient-to-br from-orange-500 to-orange-600 rounded-xl p-6 text-white">
          <p className="text-orange-100">Contenido GC</p>
          <p className="text-3xl font-bold mt-2">
            {geneData.gc_content}%
          </p>
        </div>
      </div>

      {/* Size Statistics */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <h3 className="font-semibold text-gray-800 mb-4">Estadísticas de Tamaño de Genes</h3>
        <div className="grid grid-cols-2 md:grid-cols-5 gap-4">
          <div className="text-center p-4 bg-gray-50 rounded-lg">
            <p className="text-sm text-gray-500">Media</p>
            <p className="text-2xl font-bold text-gray-800">
              {geneData.size_statistics.mean.toFixed(0)}
            </p>
            <p className="text-xs text-gray-400">bp</p>
          </div>
          <div className="text-center p-4 bg-gray-50 rounded-lg">
            <p className="text-sm text-gray-500">Mediana</p>
            <p className="text-2xl font-bold text-gray-800">
              {geneData.size_statistics.median.toFixed(0)}
            </p>
            <p className="text-xs text-gray-400">bp</p>
          </div>
          <div className="text-center p-4 bg-gray-50 rounded-lg">
            <p className="text-sm text-gray-500">Mínimo</p>
            <p className="text-2xl font-bold text-emerald-600">
              {geneData.size_statistics.min}
            </p>
            <p className="text-xs text-gray-400">bp</p>
          </div>
          <div className="text-center p-4 bg-gray-50 rounded-lg">
            <p className="text-sm text-gray-500">Máximo</p>
            <p className="text-2xl font-bold text-red-600">
              {geneData.size_statistics.max.toLocaleString()}
            </p>
            <p className="text-xs text-gray-400">bp</p>
          </div>
          <div className="text-center p-4 bg-gray-50 rounded-lg">
            <p className="text-sm text-gray-500">Desv. Est.</p>
            <p className="text-2xl font-bold text-gray-800">
              {geneData.size_statistics.std.toFixed(0)}
            </p>
            <p className="text-xs text-gray-400">bp</p>
          </div>
        </div>
      </div>

      {/* Charts Row */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Histogram - Gene Size Distribution */}
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">
            Distribución de Tamaños de Genes
          </h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={sizeHistogramData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
              <XAxis dataKey="range" tick={{ fontSize: 11 }} angle={-45} textAnchor="end" height={60} />
              <YAxis tick={{ fontSize: 12 }} />
              <Tooltip 
                formatter={(value) => [value, 'Genes']}
                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgba(0,0,0,0.1)' }}
              />
              <Bar dataKey="count" fill="#10b981" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Scatter Plot - GC Content vs Length */}
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">
            Contenido GC vs Longitud del Gen
          </h3>
          <ResponsiveContainer width="100%" height={300}>
            <ScatterChart margin={{ bottom: 20 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
              <XAxis 
                dataKey="x" 
                type="number" 
                name="Longitud" 
                tick={{ fontSize: 12 }}
                label={{ value: 'Longitud (bp)', position: 'bottom', offset: 0 }}
              />
              <YAxis 
                dataKey="y" 
                type="number" 
                name="GC%" 
                tick={{ fontSize: 12 }}
                domain={[20, 80]}
                label={{ value: 'GC %', angle: -90, position: 'insideLeft' }}
              />
              <ZAxis range={[20, 20]} />
              <Tooltip 
                cursor={{ strokeDasharray: '3 3' }}
                formatter={(value, name) => [
                  name === 'x' ? `${value} bp` : `${value}%`,
                  name === 'x' ? 'Longitud' : 'GC'
                ]}
                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgba(0,0,0,0.1)' }}
              />
              <Scatter data={gcScatterData} fill="#8b5cf6" opacity={0.6} />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Gene Table with Search */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4 mb-4">
          <h3 className="font-semibold text-gray-800">
            Lista de Genes ({paginatedGenes.total.toLocaleString()} total)
          </h3>
          <div className="flex items-center space-x-2">
            <div className="relative">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={handleKeyPress}
                placeholder="Buscar por locus o producto..."
                className="pl-10 pr-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-emerald-500 focus:border-emerald-500 w-64"
              />
              <MagnifyingGlassIcon className="h-5 w-5 text-gray-400 absolute left-3 top-1/2 transform -translate-y-1/2" />
            </div>
            <button
              onClick={handleSearch}
              className="px-4 py-2 bg-emerald-600 text-white rounded-lg hover:bg-emerald-700 transition-colors"
            >
              Buscar
            </button>
          </div>
        </div>

        {/* AG Grid Table */}
        <div className="ag-theme-alpine" style={{ height: 400 }}>
          <AgGridReact
            rowData={paginatedGenes.genes}
            columnDefs={columnDefs}
            defaultColDef={{
              resizable: true,
            }}
            animateRows={true}
            pagination={false}
            loading={isLoading}
            overlayLoadingTemplate='<span class="text-gray-500">Cargando genes...</span>'
            overlayNoRowsTemplate='<span class="text-gray-500">No se encontraron genes</span>'
          />
        </div>

        {/* Pagination */}
        {paginatedGenes.total_pages > 1 && (
          <div className="flex justify-center items-center space-x-2 mt-4">
            <button
              onClick={() => loadGenes(currentPage - 1, searchQuery)}
              disabled={currentPage === 1}
              className="px-3 py-1 border border-gray-300 rounded-md disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-50"
            >
              Anterior
            </button>
            <span className="text-sm text-gray-600">
              Página {currentPage} de {paginatedGenes.total_pages}
            </span>
            <button
              onClick={() => loadGenes(currentPage + 1, searchQuery)}
              disabled={currentPage === paginatedGenes.total_pages}
              className="px-3 py-1 border border-gray-300 rounded-md disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-50"
            >
              Siguiente
            </button>
          </div>
        )}
      </div>

      {/* Info Box */}
      <div className="bg-emerald-50 border border-emerald-200 rounded-xl p-4">
        <div className="flex">
          <div className="flex-shrink-0">
            <svg className="h-5 w-5 text-emerald-400" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="ml-3">
            <p className="text-sm text-emerald-700">
              E. coli K-12 MG1655 tiene aproximadamente 4,300 genes que cubren ~87% del genoma. 
              El tamaño promedio de un gen es de ~950 pb. Los genes están distribuidos casi 
              equitativamente entre las dos hebras del ADN.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
