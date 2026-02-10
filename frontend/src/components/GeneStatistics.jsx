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

  // Componente para renderizar la hebra con dirección clara
  const StrandCellRenderer = (props) => {
    const isForward = props.value === 1
    return (
      <span style={{
        display: 'inline-flex',
        alignItems: 'center',
        justifyContent: 'center',
        gap: '4px',
        padding: '4px 8px',
        borderRadius: '6px',
        fontSize: '10px',
        fontWeight: '600',
        background: isForward ? '#f0fdfa' : '#f5f3ff',
        color: isForward ? '#0d9488' : '#7c3aed',
        border: `1px solid ${isForward ? '#99f6e4' : '#ddd6fe'}`,
        whiteSpace: 'nowrap'
      }}>
        {isForward ? "→ 5'→3'" : "← 3'→5'"}
      </span>
    )
  }

  // Componente para renderizar codones
  const CodonCellRenderer = (props) => {
    const codon = props.value
    if (!codon) {
      return <span style={{ color: '#cbd5e1', fontStyle: 'italic', fontSize: '11px' }}>—</span>
    }

    // Colores para diferentes codones
    const isStart = ['ATG', 'GTG', 'TTG', 'CTG'].includes(codon)
    const isStop = ['TAA', 'TAG', 'TGA'].includes(codon)

    return (
      <span style={{
        fontFamily: 'monospace',
        fontSize: '11px',
        fontWeight: '600',
        padding: '2px 6px',
        borderRadius: '4px',
        background: isStart ? '#dcfce7' : isStop ? '#fee2e2' : '#f1f5f9',
        color: isStart ? '#15803d' : isStop ? '#dc2626' : '#475569'
      }}>
        {codon}
      </span>
    )
  }

  const columnDefs = useMemo(() => [
    {
      field: 'locus_tag',
      headerName: 'Locus Tag',
      width: 140,
      sortable: true,
      filter: true,
      pinned: 'left',
      cellStyle: { fontWeight: '600', color: '#0f766e', fontFamily: 'monospace', fontSize: '12px' }
    },
    {
      field: 'gene_name',
      headerName: 'Nombre Gen',
      width: 110,
      sortable: true,
      filter: true,
      valueFormatter: p => p.value || '—',
      cellStyle: params => params.value ? { color: '#0d9488', fontWeight: '500' } : { color: '#94a3b8' }
    },
    {
      field: 'start',
      headerName: 'Inicio (bp)',
      width: 110,
      sortable: true,
      valueFormatter: p => p.value?.toLocaleString(),
      cellStyle: { fontFamily: 'monospace', fontSize: '12px' }
    },
    {
      field: 'end',
      headerName: 'Fin (bp)',
      width: 110,
      sortable: true,
      valueFormatter: p => p.value?.toLocaleString(),
      cellStyle: { fontFamily: 'monospace', fontSize: '12px' }
    },
    {
      field: 'length',
      headerName: 'Longitud (bp)',
      width: 120,
      sortable: true,
      valueFormatter: p => p.value?.toLocaleString(),
      cellStyle: { fontWeight: '500', color: '#0f172a' }
    },
    {
      field: 'strand',
      headerName: 'Hebra',
      width: 95,
      sortable: true,
      cellRenderer: StrandCellRenderer,
      tooltipValueGetter: params => {
        const isForward = params.value === 1
        return isForward
          ? "Hebra Forward: 5'→3' (sentido directo)"
          : "Hebra Reverse: 3'→5' (sentido complementario)"
      }
    },
    {
      field: 'gc_content',
      headerName: 'GC%',
      width: 85,
      sortable: true,
      valueFormatter: p => p.value ? `${p.value.toFixed(1)}%` : '—',
      cellStyle: params => {
        const val = params.value
        return {
          fontWeight: '600',
          color: val > 60 ? '#dc2626' : val < 40 ? '#2563eb' : '#0f766e'
        }
      }
    },
    {
      field: 'start_codon',
      headerName: 'Inicio',
      width: 80,
      sortable: true,
      cellRenderer: CodonCellRenderer,
      tooltipValueGetter: params => {
        const codon = params.value
        if (!codon) return 'Sin codón de inicio detectado'
        const names = {
          'ATG': 'ATG (Metionina) - Inicio canónico',
          'GTG': 'GTG (Valina) - Inicio alternativo',
          'TTG': 'TTG (Leucina) - Inicio alternativo',
          'CTG': 'CTG (Leucina) - Inicio raro'
        }
        return names[codon] || `${codon} - Codón de inicio`
      }
    },
    {
      field: 'stop_codon',
      headerName: 'Parada',
      width: 80,
      sortable: true,
      cellRenderer: CodonCellRenderer,
      tooltipValueGetter: params => {
        const codon = params.value
        if (!codon) return 'Sin codón de parada detectado'
        const names = {
          'TAA': 'TAA (Ocre) - Codón de parada',
          'TAG': 'TAG (Ámbar) - Codón de parada',
          'TGA': 'TGA (Ópalo) - Codón de parada'
        }
        return names[codon] || `${codon} - Codón de parada`
      }
    },
    {
      field: 'has_introns',
      headerName: 'Intrón',
      width: 75,
      sortable: true,
      cellRenderer: params => {
        const hasIntrons = params.value
        return hasIntrons
          ? <span style={{ fontSize: '11px', fontWeight: '600', color: '#ea580c' }}>Sí</span>
          : <span style={{ fontSize: '11px', color: '#94a3b8' }}>No</span>
      },
      tooltipValueGetter: params => params.value
        ? 'Gen con intrones (splicing detectado)'
        : 'Gen sin intrones (secuencia continua)'
    },
    {
      field: 'protein_id',
      headerName: 'Protein ID',
      width: 130,
      sortable: true,
      filter: true,
      valueGetter: params => params.data.protein_id || null,
      valueFormatter: params => params.value || '—',
      cellStyle: params => ({
        fontFamily: params.value ? 'monospace' : 'inherit',
        fontSize: params.value ? '11px' : '12px',
        color: params.value ? '#7c3aed' : '#94a3b8',
        fontWeight: params.value ? '500' : 'normal',
        fontStyle: params.value ? 'normal' : 'italic'
      }),
      tooltipValueGetter: params => params.value || 'Sin ID de proteína asignado'
    },
    {
      field: 'product',
      headerName: 'Producto / Función',
      flex: 1,
      minWidth: 220,
      sortable: true,
      filter: true,
      tooltipField: 'product',
      cellStyle: { fontSize: '12px', lineHeight: '1.4' },
      valueFormatter: p => p.value || 'Proteína hipotética'
    },
  ], [])

  // Strand distribution
  const strandStats = useMemo(() => {
    if (!paginatedGenes.genes.length) return null
    const fwd = paginatedGenes.genes.filter(g => g.strand === 1).length
    const rev = paginatedGenes.genes.filter(g => g.strand === -1).length
    const total = fwd + rev
    return { fwd, rev, total, fwdPct: ((fwd / total) * 100).toFixed(1), revPct: ((rev / total) * 100).toFixed(1) }
  }, [paginatedGenes.genes])

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

        <div className="ag-theme-alpine" style={{ height: 500, width: '100%' }}>
          <AgGridReact
            rowData={paginatedGenes.genes}
            columnDefs={columnDefs}
            defaultColDef={{
              resizable: true,
              sortable: true,
              filter: false,
            }}
            animateRows={true}
            pagination={false}
            loading={isLoading}
            rowHeight={38}
            headerHeight={42}
            suppressHorizontalScroll={false}
            enableCellTextSelection={true}
            ensureDomOrder={true}
            overlayLoadingTemplate='<div style="padding:20px;"><div style="border:3px solid #14b8a6;border-top-color:transparent;border-radius:50%;width:40px;height:40px;animation:spin 0.8s linear infinite;margin:0 auto 10px;"></div><span style="color:#64748b;">Cargando genes...</span></div>'
            overlayNoRowsTemplate='<div style="padding:30px;color:#64748b;"><svg style="width:48px;height:48px;margin:0 auto 10px;display:block;opacity:0.3;" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path></svg><p style="font-weight:500;margin-bottom:5px;">No se encontraron genes</p><p style="font-size:12px;">Intenta modificar los filtros de búsqueda</p></div>'
          />
        </div>

        {/* Quick Stats under table */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 mt-4 pt-4 border-t border-slate-200">
          <div className="text-center">
            <p className="text-xs text-slate-500 mb-1">Genes en esta página</p>
            <p className="text-lg font-bold text-slate-800">{paginatedGenes.genes.length}</p>
          </div>
          <div className="text-center">
            <p className="text-xs text-slate-500 mb-1">Total de genes</p>
            <p className="text-lg font-bold text-teal-700">{paginatedGenes.total.toLocaleString()}</p>
          </div>
          <div className="text-center">
            <p className="text-xs text-slate-500 mb-1">Tamaño promedio</p>
            <p className="text-lg font-bold text-slate-800">
              {paginatedGenes.genes.length > 0
                ? Math.round(paginatedGenes.genes.reduce((sum, g) => sum + g.length, 0) / paginatedGenes.genes.length).toLocaleString()
                : '—'} bp
            </p>
          </div>
          <div className="text-center">
            <p className="text-xs text-slate-500 mb-1">GC% promedio</p>
            <p className="text-lg font-bold text-slate-800">
              {paginatedGenes.genes.length > 0
                ? (paginatedGenes.genes.reduce((sum, g) => sum + (g.gc_content || 0), 0) / paginatedGenes.genes.length).toFixed(1)
                : '—'}%
            </p>
          </div>
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

      {/* Info Box — Strand Distribution */}
      {strandStats && (
        <div className="bg-white border border-slate-200 rounded-xl p-5">
          <h3 className="font-semibold text-slate-800 mb-3">Distribución por Hebras</h3>
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <div className="bg-teal-50 border border-teal-100 rounded-lg p-4 text-center">
              <p className="text-xs text-teal-700 uppercase font-medium">→ Hebra Forward (5'→3')</p>
              <p className="text-2xl font-bold text-teal-700 mt-1">{strandStats.fwd.toLocaleString()}</p>
              <p className="text-sm text-teal-600">{strandStats.fwdPct}%</p>
            </div>
            <div className="bg-violet-50 border border-violet-100 rounded-lg p-4 text-center">
              <p className="text-xs text-violet-700 uppercase font-medium">← Hebra Reverse (3'→5')</p>
              <p className="text-2xl font-bold text-violet-700 mt-1">{strandStats.rev.toLocaleString()}</p>
              <p className="text-sm text-violet-600">{strandStats.revPct}%</p>
            </div>
            <div className="flex flex-col justify-center">
              <p className="text-xs text-slate-500 mb-2">Proporción Forward / Reverse</p>
              <div className="w-full h-4 bg-slate-100 rounded-full overflow-hidden flex">
                <div className="bg-teal-400 h-full transition-all" style={{ width: `${strandStats.fwdPct}%` }}></div>
                <div className="bg-violet-400 h-full transition-all" style={{ width: `${strandStats.revPct}%` }}></div>
              </div>
              <div className="flex justify-between mt-1 text-xs text-slate-500">
                <span>5'→3' ({strandStats.fwdPct}%)</span>
                <span>3'→5' ({strandStats.revPct}%)</span>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Categorías de Genes por Tamaño */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <h3 className="font-semibold text-slate-800 mb-4">Categorías por Tamaño de Gen</h3>
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
          {[
            {
              label: 'Muy Pequeños',
              range: '< 300 bp',
              count: paginatedGenes.genes.filter(g => g.length < 300).length,
              color: 'bg-blue-50 border-blue-200 text-blue-700'
            },
            {
              label: 'Pequeños',
              range: '300-900 bp',
              count: paginatedGenes.genes.filter(g => g.length >= 300 && g.length < 900).length,
              color: 'bg-emerald-50 border-emerald-200 text-emerald-700'
            },
            {
              label: 'Medianos',
              range: '900-1500 bp',
              count: paginatedGenes.genes.filter(g => g.length >= 900 && g.length < 1500).length,
              color: 'bg-amber-50 border-amber-200 text-amber-700'
            },
            {
              label: 'Grandes',
              range: '≥ 1500 bp',
              count: paginatedGenes.genes.filter(g => g.length >= 1500).length,
              color: 'bg-red-50 border-red-200 text-red-700'
            }
          ].map((cat, i) => (
            <div key={i} className={`${cat.color} border rounded-lg p-4 text-center`}>
              <p className="text-xs font-medium uppercase tracking-wide mb-1">{cat.label}</p>
              <p className="text-2xl font-bold mb-1">{cat.count}</p>
              <p className="text-xs opacity-75">{cat.range}</p>
            </div>
          ))}
        </div>
      </div>

      {/* Genes con Características Especiales */}
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
        <div className="bg-white rounded-lg border border-slate-200 p-4">
          <div className="flex items-start justify-between mb-2">
            <div>
              <p className="text-xs text-slate-500 uppercase font-medium">Genes Nombrados</p>
              <p className="text-2xl font-bold text-slate-800 mt-1">
                {paginatedGenes.genes.filter(g => g.gene_name).length}
              </p>
            </div>
            <div className="w-10 h-10 bg-teal-50 rounded-lg flex items-center justify-center">
              <svg className="w-5 h-5 text-teal-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 7h.01M7 3h5c.512 0 1.024.195 1.414.586l7 7a2 2 0 010 2.828l-7 7a2 2 0 01-2.828 0l-7-7A1.994 1.994 0 013 12V7a4 4 0 014-4z" />
              </svg>
            </div>
          </div>
          <p className="text-xs text-slate-600">Genes con nombre asignado</p>
        </div>

        <div className="bg-white rounded-lg border border-slate-200 p-4">
          <div className="flex items-start justify-between mb-2">
            <div>
              <p className="text-xs text-slate-500 uppercase font-medium">Con Proteína ID</p>
              <p className="text-2xl font-bold text-slate-800 mt-1">
                {paginatedGenes.genes.filter(g => g.protein_id).length}
              </p>
            </div>
            <div className="w-10 h-10 bg-purple-50 rounded-lg flex items-center justify-center">
              <svg className="w-5 h-5 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
            </div>
          </div>
          <p className="text-xs text-slate-600">Genes con ID de proteína</p>
        </div>

        <div className="bg-white rounded-lg border border-slate-200 p-4">
          <div className="flex items-start justify-between mb-2">
            <div>
              <p className="text-xs text-slate-500 uppercase font-medium">Alto contenido GC</p>
              <p className="text-2xl font-bold text-slate-800 mt-1">
                {paginatedGenes.genes.filter(g => g.gc_content > 60).length}
              </p>
            </div>
            <div className="w-10 h-10 bg-red-50 rounded-lg flex items-center justify-center">
              <svg className="w-5 h-5 text-red-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
              </svg>
            </div>
          </div>
          <p className="text-xs text-slate-600">GC% mayor a 60%</p>
        </div>

        <div className="bg-white rounded-lg border border-slate-200 p-4">
          <div className="flex items-start justify-between mb-2">
            <div>
              <p className="text-xs text-slate-500 uppercase font-medium">Bajo contenido GC</p>
              <p className="text-2xl font-bold text-slate-800 mt-1">
                {paginatedGenes.genes.filter(g => g.gc_content < 40).length}
              </p>
            </div>
            <div className="w-10 h-10 bg-blue-50 rounded-lg flex items-center justify-center">
              <svg className="w-5 h-5 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 17h8m0 0V9m0 8l-8-8-4 4-6-6" />
              </svg>
            </div>
          </div>
          <p className="text-xs text-slate-600">GC% menor a 40%</p>
        </div>
      </div>

      {/* Info Box - Resumen */}
      <div className="bg-gradient-to-r from-teal-50 to-emerald-50 border border-teal-100 rounded-xl p-5">
        <div className="flex gap-4">
          <div className="w-10 h-10 bg-teal-500 rounded-lg flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5 text-white" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="flex-1">
            <h4 className="font-semibold text-teal-900 mb-2">Resumen del Genoma</h4>
            <div className="text-sm text-teal-800 space-y-1">
              <p>
                Este genoma contiene <strong className="font-bold">{geneData.total_genes.toLocaleString()}</strong> genes,
                de los cuales <strong className="font-bold">{geneData.total_cds.toLocaleString()}</strong> son CDS (secuencias codificantes de proteínas).
              </p>
              <p>
                Contenido GC: <strong className="font-bold">{geneData.gc_content}%</strong> •
                Densidad génica: <strong className="font-bold">{geneData.gene_density} genes/Mb</strong>
              </p>
              <div className="mt-3 pt-3 border-t border-teal-200 text-xs space-y-1">
                <p>
                  <strong>Hebra Forward (→ 5'→3')</strong>: Genes codificados en la misma dirección que la secuencia de referencia.
                </p>
                <p>
                  <strong>Hebra Reverse (← 3'→5')</strong>: Genes codificados en la hebra complementaria (antisentido).
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Info Box - Codones e Intrones */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
        {/* Codones */}
        <div className="bg-white border border-slate-200 rounded-xl p-5">
          <h4 className="font-semibold text-slate-800 mb-3 flex items-center gap-2">
            <span className="w-8 h-8 bg-green-50 rounded-lg flex items-center justify-center">
              <svg className="w-4 h-4 text-green-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
              </svg>
            </span>
            Codones de Inicio y Parada
          </h4>
          <div className="space-y-3 text-sm">
            <div>
              <p className="font-medium text-slate-700 mb-2">Codones de Inicio:</p>
              <div className="space-y-1.5 ml-2">
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-green-50 text-green-700 border border-green-200">ATG</span>
                  <span className="text-xs text-slate-600">Metionina - Inicio canónico (más común)</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-green-50 text-green-700 border border-green-200">GTG</span>
                  <span className="text-xs text-slate-600">Valina - Inicio alternativo</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-green-50 text-green-700 border border-green-200">TTG</span>
                  <span className="text-xs text-slate-600">Leucina - Inicio alternativo</span>
                </div>
              </div>
            </div>
            <div>
              <p className="font-medium text-slate-700 mb-2">Codones de Parada:</p>
              <div className="space-y-1.5 ml-2">
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-red-50 text-red-700 border border-red-200">TAA</span>
                  <span className="text-xs text-slate-600">Ocre - Codón de parada</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-red-50 text-red-700 border border-red-200">TAG</span>
                  <span className="text-xs text-slate-600">Ámbar - Codón de parada</span>
                </div>
                <div className="flex items-center gap-2">
                  <span className="font-mono text-xs font-bold px-2 py-1 rounded bg-red-50 text-red-700 border border-red-200">TGA</span>
                  <span className="text-xs text-slate-600">Ópalo - Codón de parada</span>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Intrones y Exones */}
        <div className="bg-white border border-slate-200 rounded-xl p-5">
          <h4 className="font-semibold text-slate-800 mb-3 flex items-center gap-2">
            <span className="w-8 h-8 bg-purple-50 rounded-lg flex items-center justify-center">
              <svg className="w-4 h-4 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
              </svg>
            </span>
            Intrones y Exones
          </h4>
          <div className="space-y-3 text-sm text-slate-700">
            <div>
              <p className="font-medium mb-1">Exones:</p>
              <p className="text-xs text-slate-600 leading-relaxed">
                Secuencias codificantes que se mantienen en el mRNA maduro y se traducen a proteínas.
                En procariotas, la mayoría de genes son continuos (solo exones).
              </p>
            </div>
            <div>
              <p className="font-medium mb-1">Intrones:</p>
              <p className="text-xs text-slate-600 leading-relaxed">
                Secuencias no codificantes que se eliminan durante el splicing del RNA.
                Raros en procariotas, comunes en eucariotas.
              </p>
            </div>
            <div className="bg-purple-50 border border-purple-100 rounded-lg p-3 mt-3">
              <p className="text-xs text-purple-800">
                <strong>Nota:</strong> Los genes marcados con "Sí" en la columna Intrón tienen ubicaciones compuestas
                (compound locations), indicando posible splicing o secuencias discontinuas.
              </p>
            </div>
            <div className="text-xs text-slate-500 mt-2">
              <p>
                <strong>Genes con intrones:</strong> {paginatedGenes.genes.filter(g => g.has_introns).length} de {paginatedGenes.genes.length} en esta página
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
