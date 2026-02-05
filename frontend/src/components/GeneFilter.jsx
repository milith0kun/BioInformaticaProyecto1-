/**
 * GeneFilter Component
 * Filtrado avanzado de genes por tamaño, grupos funcionales y búsqueda individual
 */
import { useState, useEffect, useCallback } from 'react'
import { AgGridReact } from 'ag-grid-react'
import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'
import { api } from '../services/api'
import toast from 'react-hot-toast'

export default function GeneFilter({ hasAnalysis }) {
  const [activeFilter, setActiveFilter] = useState('size')
  const [genes, setGenes] = useState([])
  const [isLoading, setIsLoading] = useState(false)
  const [filterInfo, setFilterInfo] = useState('')
  const [functionalGroups, setFunctionalGroups] = useState([])
  const [groupsSummary, setGroupsSummary] = useState([])

  // Filtros de tamaño
  const [sizeOrder, setSizeOrder] = useState('largest')
  const [sizeCount, setSizeCount] = useState(10)

  // Filtros avanzados
  const [searchQuery, setSearchQuery] = useState('')
  const [minLength, setMinLength] = useState('')
  const [maxLength, setMaxLength] = useState('')
  const [minGC, setMinGC] = useState('')
  const [maxGC, setMaxGC] = useState('')
  const [strand, setStrand] = useState('')
  const [selectedGroup, setSelectedGroup] = useState('')

  // Paginación
  const [page, setPage] = useState(1)
  const [totalPages, setTotalPages] = useState(1)
  const [totalCount, setTotalCount] = useState(0)

  // Cargar grupos funcionales al montar
  useEffect(() => {
    if (hasAnalysis) {
      loadFunctionalGroups()
      loadGroupsSummary()
    }
  }, [hasAnalysis])

  const loadFunctionalGroups = async () => {
    try {
      const result = await api.getFunctionalGroups()
      setFunctionalGroups(result.groups || [])
    } catch (error) {
      console.error('Error loading functional groups:', error)
    }
  }

  const loadGroupsSummary = async () => {
    try {
      const result = await api.getGeneGroupsSummary()
      setGroupsSummary(result.groups || [])
    } catch (error) {
      console.error('Error loading groups summary:', error)
    }
  }

  const loadGenesBySize = async () => {
    setIsLoading(true)
    try {
      const result = await api.getGenesBySize(sizeOrder, sizeCount)
      setGenes(result.genes || [])
      setFilterInfo(result.filter_applied)
      setTotalCount(result.total_count)
    } catch (error) {
      toast.error('Error al cargar genes: ' + (error.response?.data?.detail || error.message))
    } finally {
      setIsLoading(false)
    }
  }

  const loadGenesByGroup = async (groupId) => {
    setIsLoading(true)
    try {
      const result = await api.getGenesByGroup(groupId)
      setGenes(result.genes || [])
      setFilterInfo(result.filter_applied)
      setTotalCount(result.total_count)
    } catch (error) {
      toast.error('Error al cargar genes: ' + (error.response?.data?.detail || error.message))
    } finally {
      setIsLoading(false)
    }
  }

  const searchGenesAdvanced = async () => {
    setIsLoading(true)
    try {
      const result = await api.searchGenesAdvanced({
        query: searchQuery || undefined,
        min_length: minLength ? parseInt(minLength) : undefined,
        max_length: maxLength ? parseInt(maxLength) : undefined,
        min_gc: minGC ? parseFloat(minGC) : undefined,
        max_gc: maxGC ? parseFloat(maxGC) : undefined,
        strand: strand || undefined,
        page,
        page_size: 50
      })
      setGenes(result.genes || [])
      setFilterInfo(result.filter_applied)
      setTotalCount(result.total_count)
      setTotalPages(result.total_pages)
    } catch (error) {
      toast.error('Error en búsqueda: ' + (error.response?.data?.detail || error.message))
    } finally {
      setIsLoading(false)
    }
  }

  const columnDefs = [
    { field: 'locus_tag', headerName: 'Locus Tag', width: 120, sortable: true, filter: true },
    { field: 'product', headerName: 'Producto', flex: 1, minWidth: 200, sortable: true, filter: true },
    { field: 'length', headerName: 'Longitud', width: 100, sortable: true, valueFormatter: p => `${p.value?.toLocaleString()} bp` },
    { field: 'gc_content', headerName: 'GC%', width: 80, sortable: true, valueFormatter: p => p.value?.toFixed(1) },
    { field: 'start', headerName: 'Inicio', width: 100, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'end', headerName: 'Fin', width: 100, valueFormatter: p => p.value?.toLocaleString() },
    { field: 'strand', headerName: 'Hebra', width: 70 }
  ]

  const clearFilters = () => {
    setSearchQuery('')
    setMinLength('')
    setMaxLength('')
    setMinGC('')
    setMaxGC('')
    setStrand('')
    setSelectedGroup('')
    setGenes([])
    setFilterInfo('')
    setTotalCount(0)
    setPage(1)
  }

  if (!hasAnalysis) {
    return (
      <div className="text-center py-20">
        <div className="w-20 h-20 mx-auto bg-slate-100 rounded-2xl flex items-center justify-center mb-6">
          <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
          </svg>
        </div>
        <h2 className="text-2xl font-bold text-slate-800 mb-2">Filtrado Avanzado de Genes</h2>
        <p className="text-slate-500 max-w-md mx-auto">
          Ejecute un análisis para acceder a las herramientas de filtrado por tamaño, grupos funcionales y búsqueda individual.
        </p>
      </div>
    )
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div>
        <h2 className="text-2xl font-bold text-slate-800">Filtrado Avanzado de Genes</h2>
        <p className="text-slate-500">Identifique genes por tamaño, función o características específicas</p>
      </div>

      {/* Tabs de filtro */}
      <div className="flex gap-2 bg-slate-100 p-1 rounded-xl">
        {[
          { id: 'size', name: 'Por Tamaño', icon: '📏' },
          { id: 'group', name: 'Por Grupo', icon: '🏷️' },
          { id: 'search', name: 'Búsqueda Avanzada', icon: '🔍' }
        ].map(tab => (
          <button
            key={tab.id}
            onClick={() => { setActiveFilter(tab.id); clearFilters() }}
            className={`flex-1 px-4 py-2.5 rounded-lg text-sm font-medium transition-all ${
              activeFilter === tab.id
                ? 'bg-white text-teal-700 shadow-sm'
                : 'text-slate-600 hover:text-slate-800'
            }`}
          >
            {tab.icon} {tab.name}
          </button>
        ))}
      </div>

      {/* Filtro por Tamaño */}
      {activeFilter === 'size' && (
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">Seleccionar Genes por Tamaño</h3>
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Orden</label>
              <select
                value={sizeOrder}
                onChange={(e) => setSizeOrder(e.target.value)}
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              >
                <option value="largest">📏 Más largos</option>
                <option value="smallest">🔬 Más cortos</option>
              </select>
            </div>
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Cantidad</label>
              <select
                value={sizeCount}
                onChange={(e) => setSizeCount(parseInt(e.target.value))}
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              >
                <option value={10}>Top 10</option>
                <option value={20}>Top 20</option>
                <option value={50}>Top 50</option>
                <option value={100}>Top 100</option>
              </select>
            </div>
            <div className="flex items-end">
              <button
                onClick={loadGenesBySize}
                disabled={isLoading}
                className="w-full px-4 py-2 bg-teal-600 text-white rounded-lg font-medium hover:bg-teal-700 transition-colors disabled:opacity-50"
              >
                {isLoading ? 'Cargando...' : 'Buscar'}
              </button>
            </div>
          </div>

          {/* Quick stats de tamaño */}
          <div className="mt-4 grid grid-cols-2 sm:grid-cols-4 gap-3">
            <button
              onClick={() => { setSizeOrder('largest'); setSizeCount(10); }}
              className="p-3 bg-red-50 border border-red-100 rounded-lg hover:bg-red-100 transition-all text-left"
            >
              <span className="text-2xl">📏</span>
              <p className="text-sm font-medium text-red-700 mt-1">10 Más Largos</p>
            </button>
            <button
              onClick={() => { setSizeOrder('smallest'); setSizeCount(10); }}
              className="p-3 bg-green-50 border border-green-100 rounded-lg hover:bg-green-100 transition-all text-left"
            >
              <span className="text-2xl">🔬</span>
              <p className="text-sm font-medium text-green-700 mt-1">10 Más Cortos</p>
            </button>
            <button
              onClick={() => { setSizeOrder('largest'); setSizeCount(50); }}
              className="p-3 bg-amber-50 border border-amber-100 rounded-lg hover:bg-amber-100 transition-all text-left"
            >
              <span className="text-2xl">📊</span>
              <p className="text-sm font-medium text-amber-700 mt-1">Top 50 Largos</p>
            </button>
            <button
              onClick={() => { setSizeOrder('smallest'); setSizeCount(50); }}
              className="p-3 bg-blue-50 border border-blue-100 rounded-lg hover:bg-blue-100 transition-all text-left"
            >
              <span className="text-2xl">🧪</span>
              <p className="text-sm font-medium text-blue-700 mt-1">Top 50 Cortos</p>
            </button>
          </div>
        </div>
      )}

      {/* Filtro por Grupo Funcional */}
      {activeFilter === 'group' && (
        <div className="space-y-4">
          {/* Resumen de grupos */}
          <div className="bg-white rounded-xl border border-slate-200 p-5">
            <h3 className="font-semibold text-slate-800 mb-4">Grupos Funcionales</h3>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
              {groupsSummary.map(group => (
                <button
                  key={group.id}
                  onClick={() => { setSelectedGroup(group.id); loadGenesByGroup(group.id) }}
                  className={`p-4 rounded-xl border transition-all text-left ${
                    selectedGroup === group.id
                      ? 'border-teal-500 bg-teal-50'
                      : 'border-slate-200 hover:border-teal-200 hover:bg-slate-50'
                  }`}
                >
                  <div className="flex justify-between items-start mb-2">
                    <span className="text-lg">
                      {group.id === 'metabolism' ? '⚗️' :
                       group.id === 'transport' ? '🚚' :
                       group.id === 'regulation' ? '🎛️' :
                       group.id === 'dna_rna' ? '🧬' :
                       group.id === 'protein' ? '🔬' :
                       group.id === 'cell_structure' ? '🏗️' :
                       group.id === 'stress_response' ? '🛡️' : '❓'}
                    </span>
                    <span className="text-xs bg-teal-100 text-teal-700 px-2 py-0.5 rounded-full font-medium">
                      {group.percentage.toFixed(1)}%
                    </span>
                  </div>
                  <h4 className="font-medium text-slate-800">{group.name}</h4>
                  <p className="text-sm text-slate-500 mt-0.5">{group.count.toLocaleString()} genes</p>
                  <p className="text-xs text-slate-400 mt-1">~{group.avg_length.toFixed(0)} bp promedio</p>
                </button>
              ))}
            </div>
          </div>

          {/* Descripción del grupo seleccionado */}
          {selectedGroup && (
            <div className="bg-teal-50 border border-teal-100 rounded-xl p-4">
              <p className="text-sm text-teal-800">
                <strong>Grupo seleccionado:</strong> {groupsSummary.find(g => g.id === selectedGroup)?.description}
              </p>
            </div>
          )}
        </div>
      )}

      {/* Búsqueda Avanzada */}
      {activeFilter === 'search' && (
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">Filtros Avanzados</h3>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
            {/* Búsqueda de texto */}
            <div className="lg:col-span-2">
              <label className="block text-sm font-medium text-slate-700 mb-1">Buscar (locus_tag o producto)</label>
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                placeholder="Ej: rpoB, polymerase, kinase..."
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              />
            </div>

            {/* Hebra */}
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Hebra</label>
              <select
                value={strand}
                onChange={(e) => setStrand(e.target.value)}
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              >
                <option value="">Ambas</option>
                <option value="+">+ (Forward)</option>
                <option value="-">- (Reverse)</option>
              </select>
            </div>

            {/* Rango de longitud */}
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Longitud mínima (bp)</label>
              <input
                type="number"
                value={minLength}
                onChange={(e) => setMinLength(e.target.value)}
                placeholder="Ej: 500"
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Longitud máxima (bp)</label>
              <input
                type="number"
                value={maxLength}
                onChange={(e) => setMaxLength(e.target.value)}
                placeholder="Ej: 2000"
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              />
            </div>

            {/* Rango de GC */}
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">GC mínimo (%)</label>
              <input
                type="number"
                step="0.1"
                value={minGC}
                onChange={(e) => setMinGC(e.target.value)}
                placeholder="Ej: 40"
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">GC máximo (%)</label>
              <input
                type="number"
                step="0.1"
                value={maxGC}
                onChange={(e) => setMaxGC(e.target.value)}
                placeholder="Ej: 60"
                className="w-full px-3 py-2 border border-slate-200 rounded-lg focus:ring-2 focus:ring-teal-500"
              />
            </div>
          </div>

          {/* Botones */}
          <div className="flex gap-3 mt-4">
            <button
              onClick={searchGenesAdvanced}
              disabled={isLoading}
              className="px-6 py-2 bg-teal-600 text-white rounded-lg font-medium hover:bg-teal-700 transition-colors disabled:opacity-50"
            >
              {isLoading ? 'Buscando...' : '🔍 Buscar'}
            </button>
            <button
              onClick={clearFilters}
              className="px-6 py-2 border border-slate-200 text-slate-600 rounded-lg font-medium hover:bg-slate-50 transition-colors"
            >
              Limpiar Filtros
            </button>
          </div>
        </div>
      )}

      {/* Resultados */}
      {genes.length > 0 && (
        <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
          <div className="p-4 border-b border-slate-100 flex flex-col sm:flex-row justify-between items-start sm:items-center gap-2">
            <div>
              <h3 className="font-semibold text-slate-800">Resultados</h3>
              <p className="text-sm text-slate-500">{filterInfo}</p>
            </div>
            <span className="text-sm text-teal-600 font-medium">
              {totalCount.toLocaleString()} genes encontrados
            </span>
          </div>
          
          <div className="ag-theme-alpine" style={{ height: 400 }}>
            <AgGridReact
              rowData={genes}
              columnDefs={columnDefs}
              defaultColDef={{ resizable: true }}
              animateRows={true}
              pagination={false}
            />
          </div>

          {/* Paginación para búsqueda avanzada */}
          {activeFilter === 'search' && totalPages > 1 && (
            <div className="p-4 border-t border-slate-100 flex justify-center items-center gap-2">
              <button
                onClick={() => { setPage(p => Math.max(1, p - 1)); searchGenesAdvanced() }}
                disabled={page === 1}
                className="px-3 py-1.5 border border-slate-200 rounded-lg disabled:opacity-50 hover:bg-slate-50 text-sm"
              >
                Anterior
              </button>
              <span className="text-sm text-slate-600">
                Página {page} de {totalPages}
              </span>
              <button
                onClick={() => { setPage(p => Math.min(totalPages, p + 1)); searchGenesAdvanced() }}
                disabled={page === totalPages}
                className="px-3 py-1.5 border border-slate-200 rounded-lg disabled:opacity-50 hover:bg-slate-50 text-sm"
              >
                Siguiente
              </button>
            </div>
          )}
        </div>
      )}

      {/* Info Box */}
      <div className="bg-slate-50 border border-slate-200 rounded-xl p-5">
        <div className="flex gap-3">
          <div className="w-10 h-10 bg-slate-700 rounded-xl flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
            </svg>
          </div>
          <div>
            <h4 className="font-semibold text-slate-800 mb-1">Algoritmo Optimizado</h4>
            <p className="text-sm text-slate-600">
              El filtrado utiliza algoritmos optimizados para maximizar la precisión. Los genes se clasifican 
              automáticamente en grupos funcionales basados en palabras clave del producto génico. La búsqueda 
              avanzada permite combinar múltiples criterios para identificar genes específicos.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
