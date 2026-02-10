/**
 * RNAAnalysis Component — Clean Laboratory Edition
 * Advanced analysis of tRNA and rRNA genes with charts and detailed categorization.
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
  Cell
} from 'recharts'

const COLORS = ['#2563eb', '#4f46e5', '#7c3aed', '#db2777', '#0f172a']

const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900/95 backdrop-blur-xl border border-white/10 p-4 rounded-2xl shadow-2xl animate-in fade-in zoom-in-95">
        <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">{label}</p>
        <p className="text-xs font-bold text-white uppercase tracking-tight">
          Genes: <span className="text-blue-200">{payload[0].value}</span>
        </p>
      </div>
    )
  }
  return null
}

export default function RNAAnalysis() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState('summary')

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      const result = await api.getRNAAnalysis()
      setData(result)
    } catch (e) {
      console.error("Error cargando análisis de RNA:", e)
    } finally {
      setLoading(false)
    }
  }

  const coverageChartData = useMemo(() => {
    if (!data?.amino_acid_coverage) return []
    return Object.entries(data.amino_acid_coverage)
      .sort((a, b) => b[1] - a[1])
      .map(([aa, count]) => ({ name: aa, count }))
  }, [data])

  if (loading && !data) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse">Analizando Transcripción...</p>
      </div>
    )
  }

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Header Section */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            Análisis de <span className="text-blue-600">tRNA y rRNA</span>
          </h2>
          <div className="flex flex-wrap items-center gap-4">
            <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-blue-100 shadow-sm">
              {data.organism || 'Genoma en Análisis'}
            </span>
            <p className="text-[10px] font-bold text-slate-600 uppercase tracking-widest italic">
              {data.total_trna} tRNA, {data.total_rrna} genes rRNA detectados
            </p>
          </div>
        </div>
        
        {/* Tab Navigation */}
        <div className="flex bg-slate-50 p-1.5 rounded-2xl border border-slate-100 relative z-10">
          {[
            { id: 'summary', label: 'Resumen', icon: '📊' },
            { id: 'trna', label: `tRNA (${data.total_trna})`, icon: '🧬' },
            { id: 'rrna', label: `rRNA (${data.total_rrna})`, icon: '🔬' }
          ].map(tab => (
            <button
              key={tab.id}
              onClick={() => setActiveTab(tab.id)}
              className={`flex items-center gap-3 px-6 py-2.5 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all ${
                activeTab === tab.id ? 'bg-white text-blue-600 shadow-sm border border-slate-100' : 'text-slate-400 hover:text-slate-600'
              }`}
            >
              <span>{tab.icon}</span>
              <span>{tab.label}</span>
            </button>
          ))}
        </div>
      </div>

      {activeTab === 'summary' && (
        <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
          {/* Executive Metrics */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
            <div className="bg-slate-900 rounded-3xl p-8 text-white shadow-2xl relative overflow-hidden group">
              <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 group-hover:scale-150 transition-transform duration-1000"></div>
              <p className="text-[9px] font-black text-blue-400 uppercase tracking-widest mb-3">Total tRNA</p>
              <p className="text-4xl font-black tracking-tighter">{data.total_trna}</p>
            </div>
            <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm hover:border-blue-200 transition-all">
              <p className="text-[9px] font-black text-slate-600 uppercase tracking-widest mb-3">Total rRNA</p>
              <p className="text-4xl font-black text-slate-900 tracking-tighter">{data.total_rrna}</p>
            </div>
            <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm hover:border-blue-200 transition-all">
              <p className="text-[9px] font-black text-blue-600 uppercase tracking-widest mb-3">AA Cubiertos</p>
              <div className="flex items-baseline gap-2">
                <p className="text-4xl font-black text-blue-600 tracking-tighter">{Object.keys(data.amino_acid_coverage || {}).length}</p>
                <p className="text-[10px] font-bold text-slate-400 uppercase">/ 20</p>
              </div>
            </div>
            <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm hover:border-blue-200 transition-all">
              <p className="text-[9px] font-black text-indigo-600 uppercase tracking-widest mb-3">Tipos rRNA</p>
              <p className="text-4xl font-black text-indigo-600 tracking-tighter">{Object.keys(data.rrna_types || {}).length}</p>
            </div>
          </div>

          {/* Strand Distribution */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">tRNA por Hebra</h3>
              <div className="grid grid-cols-2 gap-8">
                <div className="p-6 bg-blue-50 rounded-2xl border border-blue-100 flex flex-col items-center">
                  <p className="text-[9px] font-black text-blue-600 uppercase mb-2">→ 5'→3' Forward</p>
                  <p className="text-4xl font-black text-slate-900">{data.trna_by_strand?.forward || 0}</p>
                </div>
                <div className="p-6 bg-indigo-50 rounded-2xl border border-indigo-100 flex flex-col items-center">
                  <p className="text-[9px] font-black text-indigo-600 uppercase mb-2">← 3'→5' Reverse</p>
                  <p className="text-4xl font-black text-slate-900">{data.trna_by_strand?.reverse || 0}</p>
                </div>
              </div>
            </div>
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">rRNA por Hebra</h3>
              <div className="grid grid-cols-2 gap-8">
                <div className="p-6 bg-emerald-50 rounded-2xl border border-emerald-100 flex flex-col items-center">
                  <p className="text-[9px] font-black text-emerald-600 uppercase mb-2">→ 5'→3' Forward</p>
                  <p className="text-4xl font-black text-slate-900">{data.rrna_by_strand?.forward || 0}</p>
                </div>
                <div className="p-6 bg-rose-50 rounded-2xl border border-rose-100 flex flex-col items-center">
                  <p className="text-[9px] font-black text-rose-600 uppercase mb-2">← 3'→5' Reverse</p>
                  <p className="text-4xl font-black text-slate-900">{data.rrna_by_strand?.reverse || 0}</p>
                </div>
              </div>
            </div>
          </div>

          {/* Cobertura Chart */}
          <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
            <div className="mb-10 text-center">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-2">Cobertura de Aminoácidos por tRNA</h3>
              <p className="text-[9px] text-slate-400 font-bold uppercase tracking-widest">Número de genes tRNA identificados para cada residuo</p>
            </div>
            <ResponsiveContainer width="100%" height={350}>
              <BarChart data={coverageChartData}>
                <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
                <XAxis dataKey="name" axisLine={false} tickLine={false} tick={{ fontSize: 9, fontWeight: 900, fill: '#475569' }} />
                <YAxis axisLine={false} tickLine={false} tick={{ fontSize: 9, fill: '#94a3b8' }} />
                <Tooltip content={<CustomTooltip />} cursor={{ fill: '#f8fafc' }} />
                <Bar dataKey="count" radius={[6, 6, 0, 0]} barSize={24}>
                  {coverageChartData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>

          {/* Technical Info & Legend */}
          <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 space-y-12 overflow-hidden relative">
            <div className="absolute bottom-0 right-0 p-12 opacity-5 pointer-events-none text-9xl font-black italic">RNA</div>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-12 relative z-10">
              <div className="space-y-4">
                <h5 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.3em]">Conceptos de tRNA</h5>
                <p className="text-xs text-slate-300 leading-relaxed font-medium">
                  <span className="text-white font-bold block mb-1">RNA de transferencia (tRNA):</span> Transporta aminoácidos específicos al ribosoma durante la traducción. Cada tRNA reconoce un codón del mRNA mediante su anticodón. Los genomas bacterianos suelen optimizar el número de copias según la abundancia de codones.
                </p>
              </div>
              <div className="space-y-4 border-l border-white/5 pl-8">
                <h5 className="text-[10px] font-black text-indigo-400 uppercase tracking-[0.3em]">Conceptos de rRNA</h5>
                <p className="text-xs text-slate-300 leading-relaxed font-medium">
                  <span className="text-white font-bold block mb-1">RNA Ribosomal (rRNA):</span> Componente estructural clave del ribosoma. En procariotas, se organiza en operones que contienen las subunidades 5S, 16S y 23S. La secuencia del 16S es el estándar de oro para la filogenia bacteriana.
                </p>
              </div>
            </div>
          </div>
        </div>
      )}

      {activeTab === 'trna' && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm animate-in zoom-in-95 duration-500">
          <div className="p-8 border-b border-slate-50 bg-slate-50/30 flex justify-between items-center">
            <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Inventario Detallado de tRNA</h3>
            <span className="px-4 py-1.5 bg-blue-600 text-white text-[9px] font-black uppercase rounded-full shadow-lg">{data.total_trna} Genes</span>
          </div>
          <div className="overflow-x-auto max-h-[600px] overflow-y-auto custom-scrollbar">
            <table className="w-full text-left">
              <thead className="bg-white sticky top-0 z-10 border-b border-slate-100">
                <tr>
                  <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Identificador</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Aminoácido</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Anticodón / Nota</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Inicio</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Extensión</th>
                  <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center">Hebra</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-50">
                {data.trna_genes?.map((gene, i) => (
                  <tr key={i} className="hover:bg-blue-50/30 transition-colors group">
                    <td className="px-10 py-5 font-mono text-xs font-black text-slate-900 group-hover:text-blue-600 uppercase">{gene.locus_tag || `trna-${i+1}`}</td>
                    <td className="px-6 py-5">
                      <span className="text-[10px] font-black text-blue-700 uppercase bg-blue-50 px-2.5 py-1 rounded-md border border-blue-100">{gene.amino_acid}</span>
                    </td>
                    <td className="px-6 py-5">
                      <p className="text-[10px] text-slate-600 font-bold uppercase truncate max-w-[200px]">{gene.product || gene.anticodon}</p>
                    </td>
                    <td className="px-6 py-5 text-right font-mono text-[10px] font-bold text-slate-500">{gene.start.toLocaleString()}</td>
                    <td className="px-6 py-5 text-right font-mono text-xs font-black text-slate-900">{gene.length} pb</td>
                    <td className="px-10 py-5 text-center">
                      <span className={`text-[9px] font-black px-2.5 py-1 rounded-full border ${gene.strand === 1 ? 'text-emerald-700 bg-emerald-50 border-emerald-200' : 'text-indigo-700 bg-indigo-50 border-indigo-200'}`}>
                        {gene.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                      </span>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {activeTab === 'rrna' && (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8 animate-in zoom-in-95 duration-500">
          <div className="lg:col-span-1 space-y-6">
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Distribución de Copias</h3>
              <div className="space-y-4">
                {Object.entries(data.rrna_types || {}).map(([type, count]) => (
                  <div key={type} className="p-6 bg-slate-50 rounded-2xl border border-slate-100 flex items-center justify-between group hover:bg-white hover:border-blue-200 transition-all shadow-sm">
                    <div>
                      <p className="text-xs font-black text-slate-900 uppercase tracking-tight">{type.replace('ribosomal RNA', '').trim()}</p>
                      <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Subunidad</p>
                    </div>
                    <div className="flex items-center gap-4">
                      <span className="text-2xl font-black text-blue-600">{count}</span>
                      <div className="w-1.5 h-1.5 rounded-full bg-blue-500 group-hover:animate-ping shadow-sm"></div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>

          <div className="lg:col-span-2 bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
            <div className="p-8 border-b border-slate-50 bg-slate-50/30">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Localización Genómica de rRNA</h3>
            </div>
            <div className="overflow-x-auto max-h-[600px] overflow-y-auto custom-scrollbar">
              <table className="w-full text-left">
                <thead className="bg-white sticky top-0 z-10 border-b border-slate-100">
                  <tr>
                    <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Tipo</th>
                    <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Inicio</th>
                    <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Extensión</th>
                    <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-center">Hebra</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-50">
                  {data.rrna_genes?.map((gene, i) => (
                    <tr key={i} className="hover:bg-blue-50/30 transition-colors group">
                      <td className="px-10 py-5">
                        <p className="text-[10px] font-black text-slate-900 uppercase leading-tight">{gene.product}</p>
                        <p className="text-[9px] font-bold text-blue-600 uppercase font-mono mt-1">{gene.locus_tag}</p>
                      </td>
                      <td className="px-6 py-5 text-right font-mono text-[10px] font-bold text-slate-600">{gene.start.toLocaleString()}</td>
                      <td className="px-6 py-5 text-right font-mono text-xs font-black text-slate-900">{gene.length.toLocaleString()} pb</td>
                      <td className="px-10 py-5 text-center">
                        <span className={`text-[9px] font-black px-3 py-1 rounded-full border ${gene.strand === 1 ? 'text-emerald-700 bg-emerald-50 border-emerald-200' : 'text-indigo-700 bg-indigo-50 border-indigo-200'}`}>
                          {gene.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
