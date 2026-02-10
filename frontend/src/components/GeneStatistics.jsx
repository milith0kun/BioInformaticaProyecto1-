/**
 * GeneStatistics Component — Clean Laboratory Edition
 * Comprehensive gene analysis with paginated list and interactive charts - 100% Spanish
 */
import { useState, useEffect, useMemo } from 'react'
import api from '../services/api'
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Cell,
  PieChart,
  Pie,
  Legend,
  ScatterChart,
  Scatter,
  ZAxis
} from 'recharts'

const COLORS = ['#2563eb', '#4f46e5', '#7c3aed', '#db2777', '#0f172a']

const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900/95 backdrop-blur-xl border border-white/10 p-4 rounded-2xl shadow-2xl animate-in fade-in zoom-in-95">
        <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">{label}</p>
        {payload.map((entry, index) => (
          <p key={index} className="text-xs font-bold text-white uppercase tracking-tight">
            {entry.name === 'count' ? 'Cantidad' : entry.name}: <span className="text-blue-200">{entry.value}</span>
          </p>
        ))}
      </div>
    )
  }
  return null
}

export default function GeneStatistics({ geneData }) {
  const [paginatedGenes, setPaginatedGenes] = useState([])
  const [page, setPage] = useState(1)
  const [totalPages, setTotalPages] = useState(1)
  const [totalGenesCount, setTotalGenesCount] = useState(0)
  const [searchTerm, setSearchTerm] = useState('')
  const [loadingTable, setLoading] = useState(false)

  const pageSize = 50

  useEffect(() => {
    if (geneData) {
      loadTableData(1, '')
    }
  }, [geneData])

  const loadTableData = async (pageNum, search) => {
    setLoading(true)
    try {
      const result = await api.getGeneResults(pageNum, pageSize, search)
      setPaginatedGenes(result?.genes || [])
      setTotalPages(result?.total_pages || 1)
      setTotalGenesCount(result?.total || 0)
    } catch (e) {
      console.warn('Genes no disponibles aún:', e.message)
      setPaginatedGenes([])
      setTotalPages(1)
      setTotalGenesCount(0)
    } finally {
      setLoading(false)
    }
  }

  const handleSearch = (e) => {
    e.preventDefault()
    setPage(1)
    loadTableData(1, searchTerm)
  }

  const handlePageChange = (newPage) => {
    if (newPage < 1 || newPage > totalPages) return
    setPage(newPage)
    loadTableData(newPage, searchTerm)
  }

  if (!geneData) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8"></div>
        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Sincronizando Dataset Génico...</p>
      </div>
    )
  }

  const lengthDistributionData = [
    { name: '0-300 pb', count: geneData.length_distribution?.['0-300'] || 0 },
    { name: '300-600 pb', count: geneData.length_distribution?.['300-600'] || 0 },
    { name: '600-900 pb', count: geneData.length_distribution?.['600-900'] || 0 },
    { name: '900-1.2K pb', count: geneData.length_distribution?.['900-1200'] || 0 },
    { name: '1.2-1.5K pb', count: geneData.length_distribution?.['1200-1500'] || 0 },
    { name: '1.5-2K pb', count: geneData.length_distribution?.['1500-2000'] || 0 },
    { name: '2-3K pb', count: geneData.length_distribution?.['2000-3000'] || 0 },
    { name: '3K+ pb', count: geneData.length_distribution?.['3000+'] || 0 },
  ]

  const scatterData = useMemo(() => {
    if (!geneData.genes) return []
    return geneData.genes.slice(0, 1000).map(g => ({
      length: g.length,
      gc: g.gc_content,
      name: g.locus_tag
    }))
  }, [geneData])

  const currentPageStats = {
    avgSize: paginatedGenes.length > 0 ? (paginatedGenes.reduce((s, g) => s + g.length, 0) / paginatedGenes.length).toFixed(0) : 0,
    avgGC: paginatedGenes.length > 0 ? (paginatedGenes.reduce((s, g) => s + g.gc_content, 0) / paginatedGenes.length).toFixed(1) : 0,
    forward: paginatedGenes.filter(g => g.strand === 1).length,
    reverse: paginatedGenes.filter(g => g.strand === -1).length,
    namedCount: paginatedGenes.filter(g => g.gene_name).length,
    proteinIdCount: paginatedGenes.filter(g => g.protein_id).length,
    highGC: paginatedGenes.filter(g => g.gc_content > 60).length,
    lowGC: paginatedGenes.filter(g => g.gc_content < 40).length,
    intronCount: paginatedGenes.filter(g => g.has_introns).length,
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Tarjetas de Métricas Principales */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-slate-600 uppercase tracking-widest mb-3 leading-none">Total de Genes</p>
          <p className="text-3xl font-black text-slate-900 tracking-tighter">{geneData.total_genes?.toLocaleString() || '0'}</p>
        </div>
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-blue-600 uppercase tracking-widest mb-3 leading-none">Total de CDS</p>
          <p className="text-3xl font-black text-blue-700 tracking-tighter">{geneData.total_cds?.toLocaleString() || '0'}</p>
        </div>
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-indigo-600 uppercase tracking-widest mb-3 leading-none">Densidad Génica</p>
          <div className="flex items-baseline gap-2">
            <p className="text-3xl font-black text-indigo-700 tracking-tighter">{geneData.gene_density || '0.00'}</p>
            <p className="text-[9px] font-bold text-slate-600 uppercase tracking-tight">genes/Mb</p>
          </div>
        </div>
        <div className="bg-slate-900 rounded-3xl p-8 text-white shadow-2xl shadow-blue-900/20 relative overflow-hidden group">
          <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 group-hover:scale-150 transition-transform duration-1000"></div>
          <p className="text-[10px] font-black text-blue-300 uppercase tracking-widest mb-3 leading-none">Contenido GC</p>
          <p className="text-3xl font-black tracking-tighter text-white">{geneData.gc_content || '0.00'}%</p>
        </div>
      </div>

      {/* Estadísticas Detalladas de Tamaño */}
      <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 p-10 opacity-5 pointer-events-none"><span className="text-6xl">📏</span></div>
        <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">Estadísticas de Tamaño de Genes</h3>
        <div className="grid grid-cols-2 md:grid-cols-5 gap-8 relative z-10">
          {[
            { label: 'Media', val: geneData.size_statistics?.mean || 0, color: 'text-slate-900' },
            { label: 'Mediana', val: geneData.size_statistics?.median || 0, color: 'text-slate-900' },
            { label: 'Mínimo', val: geneData.size_statistics?.min || 0, color: 'text-blue-700' },
            { label: 'Máximo', val: geneData.size_statistics?.max || 0, color: 'text-rose-700' },
            { label: 'Desv. Est.', val: geneData.size_statistics?.std || 0, color: 'text-slate-600' }
          ].map(stat => (
            <div key={stat.label} className="space-y-1">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest">{stat.label}</p>
              <div className="flex items-baseline gap-1">
                <p className={`text-2xl font-black tracking-tighter ${stat.color}`}>{Math.round(stat.val).toLocaleString()}</p>
                <span className="text-[9px] font-bold text-slate-400 uppercase">pb</span>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Gráficos de Análisis */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Distribución de Tamaños</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={lengthDistributionData}>
              <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
              <XAxis dataKey="name" axisLine={false} tickLine={false} tick={{ fontSize: 9, fontWeight: 900, fill: '#475569' }} dy={10} />
              <YAxis axisLine={false} tickLine={false} tick={{ fontSize: 9, fill: '#94a3b8' }} />
              <Tooltip content={<CustomTooltip />} cursor={{ fill: '#f8fafc' }} />
              <Bar dataKey="count" name="Cantidad" fill="#2563eb" radius={[8, 8, 0, 0]} barSize={32}>
                {lengthDistributionData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">GC% vs Longitud (pb)</h3>
          <ResponsiveContainer width="100%" height={300}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis type="number" dataKey="length" name="Longitud" unit=" pb" axisLine={false} tickLine={false} tick={{fontSize: 9, fill: '#475569'}} />
              <YAxis type="number" dataKey="gc" name="GC%" unit="%" axisLine={false} tickLine={false} tick={{fontSize: 9, fill: '#475569'}} domain={['auto', 'auto']} />
              <ZAxis type="category" dataKey="name" name="Locus" />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)' }} />
              <Scatter name="Genes" data={scatterData} fill="#4f46e5" fillOpacity={0.6} />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Módulo Principal: Tabla de Genes */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
        <div className="p-10 border-b border-slate-50 bg-slate-50/30 flex flex-col lg:flex-row lg:items-center justify-between gap-8">
          <div>
            <h3 className="text-xl font-black text-slate-900 uppercase tracking-tighter italic">Lista de Genes</h3>
            <p className="text-[10px] font-bold text-slate-600 uppercase tracking-widest mt-1">Sincronizado: {totalGenesCount.toLocaleString()} registros detectados</p>
          </div>
          
          <form onSubmit={handleSearch} className="relative group">
            <input 
              type="text" 
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              placeholder="Buscar por locus o producto..."
              className="pl-12 pr-24 py-3 bg-white border-2 border-slate-200 rounded-2xl text-[10px] font-black uppercase tracking-widest text-slate-900 focus:outline-none focus:border-blue-500 transition-all w-full lg:w-96 placeholder-slate-400"
            />
            <svg className="w-5 h-5 text-slate-500 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-blue-500 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeWidth={2.5} /></svg>
            <button type="submit" className="absolute right-2 top-1/2 -translate-y-1/2 px-4 py-1.5 bg-blue-600 text-white text-[9px] font-black uppercase tracking-widest rounded-xl hover:bg-blue-700 transition-all shadow-lg active:scale-95">Buscar</button>
          </form>
        </div>

        <div className="overflow-x-auto custom-scrollbar">
          <table className="w-full text-left border-collapse min-w-[1300px]">
            <thead>
              <tr className="bg-white border-b border-slate-100">
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Nombre</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Inicio</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Locus Tag</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Fin</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Extensión</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center min-w-[140px]">Hebra</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">GC%</th>
                <th className="px-4 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center">Inicio</th>
                <th className="px-4 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center">Parada</th>
                <th className="px-4 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center">Intrón</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Proteína ID</th>
                <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest min-w-[300px]">Producto</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-50">
              {loadingTable ? (
                <tr>
                  <td colSpan="12" className="px-8 py-24 text-center">
                    <div className="w-10 h-10 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mx-auto mb-6"></div>
                    <p className="text-[10px] font-black text-slate-500 uppercase tracking-[0.3em] animate-pulse">Indexando Base de Datos...</p>
                  </td>
                </tr>
              ) : (
                paginatedGenes.map((gene, i) => (
                  <tr key={i} className="hover:bg-blue-50/30 transition-colors group">
                    <td className="px-6 py-5">
                      <span className="text-[10px] font-black text-blue-700 uppercase bg-blue-50 px-2 py-1 rounded-md border border-blue-100">{gene.gene_name || '—'}</span>
                    </td>
                    <td className="px-6 py-5 text-right font-mono text-[10px] font-bold text-slate-600">{gene.start.toLocaleString()}</td>
                    <td className="px-6 py-5">
                      <span className="font-mono text-[10px] font-black text-slate-900 uppercase group-hover:text-blue-600 transition-colors">{gene.locus_tag}</span>
                    </td>
                    <td className="px-6 py-5 text-right font-mono text-[10px] font-bold text-slate-600">{gene.end.toLocaleString()}</td>
                    <td className="px-6 py-5 text-right">
                      <span className="font-mono text-xs font-black text-slate-900 tracking-tighter">{gene.length.toLocaleString()} <span className="text-[9px] text-slate-400 font-bold">pb</span></span>
                    </td>
                    <td className="px-6 py-5 text-center whitespace-nowrap">
                      <span className={`inline-block text-[9px] font-black px-3 py-1 rounded-full border shadow-sm ${gene.strand === 1 ? 'text-emerald-700 bg-emerald-50 border-emerald-200' : 'text-indigo-700 bg-indigo-50 border-indigo-200'}`}>
                        {gene.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                      </span>
                    </td>
                    <td className="px-6 py-5 text-right font-mono text-xs font-black text-blue-700">{gene.gc_content.toFixed(1)}%</td>
                    <td className="px-4 py-5 text-center">
                      <span className="font-mono text-[10px] font-black text-slate-800">{gene.start_codon || '—'}</span>
                    </td>
                    <td className="px-4 py-5 text-center">
                      <span className="font-mono text-[10px] font-black text-slate-800">{gene.stop_codon || '—'}</span>
                    </td>
                    <td className="px-4 py-5 text-center">
                      <span className={`text-[9px] font-black uppercase ${gene.has_introns ? 'text-rose-700' : 'text-slate-400'}`}>{gene.has_introns ? 'SÍ' : 'No'}</span>
                    </td>
                    <td className="px-6 py-5">
                      <span className="font-mono text-[10px] font-bold text-slate-600 group-hover:text-blue-600 transition-colors">{gene.protein_id || '—'}</span>
                    </td>
                    <td className="px-6 py-5">
                      <p className="text-[10px] font-bold text-slate-700 uppercase leading-relaxed line-clamp-1 group-hover:line-clamp-none transition-all duration-500">{gene.product || 'Proteína funcional hipotética'}</p>
                    </td>
                  </tr>
                ))
              )}
            </tbody>
          </table>
        </div>

        {/* Paginación y Resumen de Página */}
        <div className="p-10 bg-slate-50/50 border-t border-slate-100 flex flex-col lg:flex-row justify-between items-center gap-8">
          <div className="flex gap-10">
            <div className="space-y-1">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest">En esta página</p>
              <p className="text-xl font-black text-slate-900">{paginatedGenes.length}</p>
            </div>
            <div className="space-y-1 border-l border-slate-200 pl-10">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest">Tamaño promedio</p>
              <p className="text-xl font-black text-slate-900">{currentPageStats.avgSize} <span className="text-xs font-bold text-slate-500 uppercase">pb</span></p>
            </div>
            <div className="space-y-1 border-l border-slate-200 pl-10">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest">GC% promedio</p>
              <p className="text-xl font-black text-blue-700">{currentPageStats.avgGC}%</p>
            </div>
          </div>

          <div className="flex items-center gap-4">
            <button 
              onClick={() => handlePageChange(page - 1)}
              disabled={page === 1 || loadingTable}
              className="px-6 py-3 bg-white border-2 border-slate-200 rounded-2xl text-[10px] font-black uppercase tracking-widest text-slate-700 hover:bg-slate-50 disabled:opacity-30 transition-all shadow-sm active:scale-95"
            >
              Anterior
            </button>
            <div className="px-8 py-3 bg-white rounded-2xl border-2 border-slate-100 text-[10px] font-black uppercase tracking-widest shadow-inner text-slate-900">
              Página <span className="text-blue-700">{page}</span> de {totalPages}
            </div>
            <button 
              onClick={() => handlePageChange(page + 1)}
              disabled={page === totalPages || loadingTable}
              className="px-6 py-3 bg-white border-2 border-slate-100 rounded-2xl text-[10px] font-black uppercase tracking-widest text-slate-700 hover:bg-slate-50 disabled:opacity-30 transition-all shadow-sm active:scale-95"
            >
              Siguiente
            </button>
          </div>
        </div>
      </div>

      {/* Distribución por Hebras y Categorías */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm flex flex-col">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Distribución por Hebras</h3>
          <div className="flex-1 flex flex-col items-center justify-center">
            <ResponsiveContainer width="100%" height={250}>
              <PieChart>
                <Pie
                  data={[
                    { name: 'Sentido (5\'→3\')', value: currentPageStats.forward },
                    { name: 'Antisentido (3\'→5\')', value: currentPageStats.reverse }
                  ]}
                  cx="50%" cy="50%" innerRadius={60} outerRadius={90} paddingAngle={10} dataKey="value"
                >
                  <Cell fill="#2563eb" cornerRadius={10} />
                  <Cell fill="#4f46e5" cornerRadius={10} />
                </Pie>
                <Tooltip contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)' }} />
                <Legend verticalAlign="bottom" iconType="circle" />
              </PieChart>
            </ResponsiveContainer>
            <div className="w-full grid grid-cols-2 gap-4 mt-8 pt-8 border-t border-slate-50">
              <div className="text-center">
                <p className="text-[9px] font-black text-blue-600 uppercase tracking-widest mb-1">→ Hebra Forward (5\'→3\')</p>
                <p className="text-2xl font-black text-slate-900">{currentPageStats.forward}</p>
                <p className="text-[10px] font-bold text-slate-600">{((currentPageStats.forward / (paginatedGenes.length || 1)) * 100).toFixed(1)}%</p>
              </div>
              <div className="text-center border-l border-slate-100">
                <p className="text-[9px] font-black text-indigo-600 uppercase tracking-widest mb-1">← Hebra Reverse (3\'→5\')</p>
                <p className="text-2xl font-black text-slate-900">{currentPageStats.reverse}</p>
                <p className="text-[10px] font-bold text-slate-600">{((currentPageStats.reverse / (paginatedGenes.length || 1)) * 100).toFixed(1)}%</p>
              </div>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">Categorías por Tamaño de Gen</h3>
          <div className="space-y-4">
            {[
              { label: 'Muy Pequeños', desc: '< 300 bp', count: geneData.length_distribution?.['0-300'], color: 'bg-cyan-500' },
              { label: 'Pequeños', desc: '300-900 bp', count: (geneData.length_distribution?.['300-600'] || 0) + (geneData.length_distribution?.['600-900'] || 0), color: 'bg-blue-500' },
              { label: 'Medianos', desc: '900-1500 bp', count: (geneData.length_distribution?.['900-1200'] || 0) + (geneData.length_distribution?.['1200-1500'] || 0), color: 'bg-indigo-500' },
              { label: 'Grandes', desc: '≥ 1500 bp', count: (geneData.length_distribution?.['1500-2000'] || 0) + (geneData.length_distribution?.['2000-3000'] || 0) + (geneData.length_distribution?.['3000+'] || 0), color: 'bg-violet-600' }
            ].map(cat => (
              <div key={cat.label} className="p-6 bg-slate-50 rounded-2xl border border-slate-100 flex items-center justify-between group hover:bg-white hover:border-blue-200 transition-all shadow-sm">
                <div>
                  <p className="text-xs font-black text-slate-900 uppercase tracking-tight">{cat.label}</p>
                  <p className="text-[9px] font-bold text-slate-600 uppercase tracking-widest">{cat.desc}</p>
                </div>
                <div className="flex items-center gap-4">
                  <span className="text-lg font-black text-slate-900">{cat.count?.toLocaleString()}</span>
                  <div className={`w-1.5 h-1.5 rounded-full ${cat.color} group-hover:animate-pulse shadow-sm`}></div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* Estadísticas de Composición */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        {[
          { label: 'Genes Nombrados', val: currentPageStats.namedCount, desc: 'Genes con nombre asignado', icon: '🏷️' },
          { label: 'Con Proteína ID', val: currentPageStats.proteinIdCount, desc: 'Genes con ID de proteína', icon: '🔑' },
          { label: 'Alto contenido GC', val: currentPageStats.highGC, desc: 'GC% mayor a 60%', icon: '🔥' },
          { label: 'Bajo contenido GC', val: currentPageStats.lowGC, desc: 'GC% menor a 40%', icon: '❄️' }
        ].map(item => (
          <div key={item.label} className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
            <div className="flex justify-between items-start mb-4">
              <p className="text-[9px] font-black text-slate-600 uppercase tracking-widest">{item.label}</p>
              <span className="text-lg">{item.icon}</span>
            </div>
            <p className="text-3xl font-black text-slate-900 mb-2">{item.val}</p>
            <p className="text-[9px] font-bold text-slate-600 uppercase">{item.desc}</p>
          </div>
        ))}
      </div>

      {/* Resumen del Genoma y Glosario */}
      <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 space-y-12 overflow-hidden relative">
        <div className="absolute bottom-0 right-0 p-12 opacity-5 pointer-events-none"><span className="text-9xl font-black italic">GENOME</span></div>
        <div className="space-y-6 relative z-10">
          <h4 className="text-2xl font-black uppercase italic tracking-tighter text-blue-400">Resumen del Genoma</h4>
          <p className="text-sm font-medium text-slate-200 leading-relaxed max-w-4xl">
            Este genoma contiene <span className="text-white font-black">{geneData.total_genes?.toLocaleString()} genes</span>, 
            de los cuales <span className="text-white font-black">{geneData.total_cds?.toLocaleString()} son CDS</span> (secuencias codificantes de proteínas). 
            La arquitectura genómica presenta un <span className="text-blue-400 font-black">Contenido GC de {geneData.gc_content}%</span> y una 
            <span className="text-indigo-400 font-black"> Densidad génica de {geneData.gene_density} genes/Mb</span>.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-12 relative z-10">
          <div className="space-y-6">
            <h5 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.3em]">Hebra y Dirección</h5>
            <div className="space-y-4">
              <p className="text-xs text-slate-300 leading-relaxed font-medium">
                <span className="text-white font-bold block mb-1">Hebra Forward (→ 5'→3'):</span> Genes codificados en la misma dirección que la secuencia de referencia.
              </p>
              <p className="text-xs text-slate-300 leading-relaxed font-medium">
                <span className="text-white font-bold block mb-1">Hebra Reverse (← 3'→5'):</span> Genes codificados en la hebra complementaria (antisentido).
              </p>
            </div>
          </div>
          <div className="space-y-6 border-l border-white/5 pl-8">
            <h5 className="text-[10px] font-black text-indigo-400 uppercase tracking-[0.3em]">Codones de Inicio y Parada</h5>
            <div className="space-y-4">
              <p className="text-[10px] text-slate-300 font-bold uppercase"><span className="text-white">Inicio:</span> ATG (Metionina), GTG (Valina) y TTG (Leucina).</p>
              <p className="text-[10px] text-slate-300 font-bold uppercase"><span className="text-white">Parada:</span> TAA (Ocre), TAG (Ámbar) y TGA (Ópalo).</p>
            </div>
          </div>
          <div className="space-y-6 border-l border-white/5 pl-8">
            <h5 className="text-[10px] font-black text-emerald-400 uppercase tracking-[0.3em]">Intrones y Exones</h5>
            <p className="text-xs text-slate-300 leading-relaxed font-medium">
              En procariotas, la mayoría de los genes son continuos (<span className="text-white">solo exones</span>). 
              Los genes con <span className="text-emerald-400">Intrones</span> detectados indican ubicaciones compuestas o posibles eventos de splicing bacteriano.
            </p>
            <p className="text-[10px] font-black text-slate-400 uppercase">Genes con intrones: {currentPageStats.intronCount} de {paginatedGenes.length} en esta página</p>
          </div>
        </div>
      </div>
    </div>
  )
}
