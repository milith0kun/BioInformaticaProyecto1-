/**
 * CAIAnalysis Component — Clean Laboratory Edition
 * Advanced analysis of Codon Adaptation Index (CAI) with interactive charts and technical context.
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
  ScatterChart,
  Scatter,
  ZAxis
} from 'recharts'

export default function CAIAnalysis() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState('summary')

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      // Fetch data for the top 50 expressed genes
      const result = await api.getCAIAnalysis(50)
      setData(result)
    } catch (e) {
      console.error("Error loading CAI data:", e)
    } finally {
      setLoading(false)
    }
  }

  // Memoized scatter data for performance (CAI vs Length)
  const scatterData = useMemo(() => {
    if (!data?.all_genes) return []
    // Limit to 1000 genes for visualization smoothness
    return data.all_genes.slice(0, 1000).map(g => ({
      length: g.length,
      cai: g.cai,
      name: g.gene_name || g.locus_tag
    }))
  }, [data])

  if (loading && !data) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse text-center">Calculando Índice de Adaptación...</p>
      </div>
    )
  }

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Topology Header */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            📈 CAI — <span className="text-blue-600">Codon Adaptation Index</span>
          </h2>
          <div className="flex flex-wrap items-center gap-4">
            <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-blue-100 shadow-sm">
              {data.organism || 'Escherichia coli'}
            </span>
            <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest italic">{data.total_genes?.toLocaleString()} genes analizados</p>
          </div>
        </div>
        
        <div className="flex bg-slate-50 p-1.5 rounded-2xl border border-slate-100 relative z-10 overflow-x-auto scrollbar-hide">
          {[
            { id: 'summary', label: 'Resumen', icon: '📊' },
            { id: 'top', label: 'Top Expresados', icon: '🔝' },
            { id: 'low', label: 'Baja Expresión', icon: '⬇️' },
            { id: 'scatter', label: 'CAI vs Longitud', icon: '📈' },
            { id: 'all', label: 'Todos los Genes', icon: '📋' }
          ].map(tab => (
            <button
              key={tab.id}
              onClick={() => setActiveTab(tab.id)}
              className={`flex items-center gap-3 px-6 py-2.5 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all whitespace-nowrap ${
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
        <div className="space-y-10">
          {/* Stats Grid */}
          <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-5 gap-6">
            {[
              { label: 'CAI Promedio', val: data.cai_stats?.mean, color: 'text-blue-600' },
              { label: 'Mediana', val: data.cai_stats?.median, color: 'text-indigo-600' },
              { label: 'Máximo', val: data.cai_stats?.max, color: 'text-slate-900' },
              { label: 'Mínimo', val: data.cai_stats?.min, color: 'text-rose-600' },
              { label: 'Desv. Estándar', val: data.cai_stats?.std, color: 'text-slate-400' }
            ].map((stat, i) => (
              <div key={i} className="bg-white rounded-3xl border-2 border-slate-100 p-6 shadow-sm hover:border-blue-200 transition-all">
                <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest mb-3">{stat.label}</p>
                <p className={`text-2xl font-black tracking-tighter ${stat.color}`}>{(stat.val || 0).toFixed(4)}</p>
              </div>
            ))}
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
            {/* Distribution Chart */}
            <div className="lg:col-span-2 bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm flex flex-col justify-center">
              <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-10 text-center">Distribución de CAI</h3>
              <p className="text-[9px] font-bold text-slate-300 uppercase tracking-widest mb-6 text-center">Genes agrupados por rango de adaptabilidad</p>
              <ResponsiveContainer width="100%" height={300}>
                <BarChart data={Object.entries(data.cai_distribution || {}).map(([range, count]) => ({ range, count }))}>
                  <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
                  <XAxis dataKey="range" axisLine={false} tickLine={false} tick={{ fontSize: 10, fontWeight: 900, fill: '#64748b' }} dy={10} />
                  <YAxis axisLine={false} tickLine={false} tick={{ fontSize: 10, fill: '#94a3b8' }} />
                  <Tooltip cursor={{ fill: '#f8fafc' }} contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)' }} />
                  <Bar dataKey="count" radius={[10, 10, 0, 0]} barSize={60}>
                    {Object.entries(data.cai_distribution || {}).map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={['#f43f5e', '#f59e0b', '#10b981', '#3b82f6', '#6366f1'][index]} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>

            {/* CAI Scale Context */}
            <div className="bg-slate-900 rounded-[2.5rem] p-8 text-white shadow-2xl shadow-blue-900/20 flex flex-col">
              <h4 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.2em] mb-8">Escala CAI</h4>
              <div className="space-y-6 flex-1">
                {[
                  { range: '0.8', label: 'Muy altamente expresado', desc: 'Proteínas ribosomales, chaperonas, factores de elongación', color: 'bg-indigo-500' },
                  { range: '0.6', label: 'Altamente expresado', desc: 'Proteínas de membrana, enzimas metabólicas principales', color: 'bg-blue-500' },
                  { range: '0.4', label: 'Moderadamente expresado', desc: 'Genes housekeeping, reguladores', color: 'bg-emerald-500' },
                  { range: '0.2', label: 'Baja expresión', desc: 'Genes específicos de condición, regulatorios', color: 'bg-amber-500' },
                  { range: '0.0', label: 'Muy baja expresión', desc: 'Genes horizontalmente transferidos, pseudogenes', color: 'bg-rose-500' }
                ].map(item => (
                  <div key={item.range} className="flex gap-4 group">
                    <div className={`w-10 h-10 rounded-xl ${item.color} flex items-center justify-center font-black text-xs flex-shrink-0 shadow-lg group-hover:scale-110 transition-transform`}>
                      {item.range}
                    </div>
                    <div className="space-y-1 min-w-0">
                      <p className="text-[10px] font-black uppercase tracking-tight">{item.label}</p>
                      <p className="text-[9px] font-bold text-slate-500 uppercase leading-relaxed truncate">{item.desc}</p>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {(activeTab === 'top' || activeTab === 'low' || activeTab === 'all') && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm animate-in slide-in-from-bottom-4 duration-500">
          <div className="p-8 border-b border-slate-50 bg-slate-50/50 flex items-center justify-between">
            <h3 className="text-[10px] font-black text-slate-900 uppercase tracking-[0.3em]">
              {activeTab === 'top' ? 'Genes con Alta Adaptación Genómica' : activeTab === 'low' ? 'Genes con Bajo Índice de Adaptación' : 'Inventario de Adaptación Codónica'}
            </h3>
            <span className="px-4 py-1 bg-blue-600 text-white text-[9px] font-black rounded-full shadow-lg uppercase tracking-widest">
              {activeTab === 'top' ? data.top_expressed?.length : activeTab === 'low' ? data.low_expressed?.length : data.all_genes?.length} Registros
            </span>
          </div>
          <div className="overflow-x-auto max-h-[600px] overflow-y-auto custom-scrollbar">
            <table className="w-full text-left">
              <thead className="bg-white sticky top-0 z-10 border-b border-slate-100">
                <tr>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Gen / Identificador</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Longitud</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Hebra</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Índice CAI</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Interpretación</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-50">
                {(activeTab === 'top' ? data.top_expressed : activeTab === 'low' ? data.low_expressed : data.all_genes || []).map((gene, i) => (
                  <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                    <td className="px-8 py-6">
                      <div className="space-y-1">
                        <p className="text-xs font-black text-slate-900 uppercase tracking-tight group-hover:text-blue-600 transition-colors">{gene.gene_name || gene.locus_tag}</p>
                        <p className="text-[10px] text-slate-400 font-bold uppercase truncate max-w-md">{gene.product || 'Proteína funcional hipotética'}</p>
                      </div>
                    </td>
                    <td className="px-8 py-6 text-right font-mono text-xs font-black text-slate-500">{gene.length?.toLocaleString()} bp</td>
                    <td className="px-8 py-6 text-right">
                      <span className={`text-[9px] font-black px-2 py-0.5 rounded ${gene.strand === 1 ? 'bg-blue-50 text-blue-600' : 'bg-indigo-50 text-indigo-600'}`}>
                        {gene.strand === 1 ? '5\'→3\'' : '3\'→5\''}
                      </span>
                    </td>
                    <td className="px-8 py-6 text-right">
                      <span className={`font-mono text-sm font-black tracking-tighter ${gene.cai > 0.7 ? 'text-blue-600' : gene.cai < 0.3 ? 'text-rose-600' : 'text-slate-700'}`}>
                        {(gene.cai || 0).toFixed(4)}
                      </span>
                    </td>
                    <td className="px-8 py-6">
                      <div className="flex items-center gap-4">
                        <div className="w-20 h-1 bg-slate-100 rounded-full overflow-hidden">
                          <div className={`h-full transition-all duration-1000 ${gene.cai > 0.7 ? 'bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.4)]' : gene.cai < 0.3 ? 'bg-rose-500' : 'bg-indigo-400'}`} style={{ width: `${(gene.cai || 0) * 100}%` }}></div>
                        </div>
                        <span className="text-[9px] font-black text-slate-400 uppercase tracking-widest whitespace-nowrap">
                          {gene.cai > 0.8 ? 'Ultra-Alta' : gene.cai > 0.6 ? 'Alta' : gene.cai > 0.4 ? 'Media' : 'Baja'}
                        </span>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {activeTab === 'scatter' && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm animate-in zoom-in-95 duration-500">
          <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-10 text-center">Correlación: Adaptación vs Longitud de Hebra</h3>
          <ResponsiveContainer width="100%" height={400}>
            <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis type="number" dataKey="length" name="Longitud" unit=" bp" axisLine={false} tickLine={false} tick={{fontSize: 10, fontBold: 900, fill: '#94a3b8'}} label={{ value: 'Extensión del Gen (bp)', position: 'insideBottom', offset: -10, fontSize: 9, fontWeight: 900, fill: '#cbd5e1' }} />
              <YAxis type="number" dataKey="cai" name="CAI" axisLine={false} tickLine={false} tick={{fontSize: 10, fontBold: 900, fill: '#94a3b8'}} domain={[0, 1]} label={{ value: 'Índice CAI', angle: -90, position: 'insideLeft', fontSize: 9, fontWeight: 900, fill: '#cbd5e1' }} />
              <ZAxis type="category" dataKey="name" name="Gen" />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)' }} />
              <Scatter name="Genes" data={scatterData} fill="#2563eb" fillOpacity={0.4}>
                {scatterData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={entry.cai > 0.7 ? '#2563eb' : entry.cai < 0.3 ? '#f43f5e' : '#6366f1'} />
                ))}
              </Scatter>
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      )}

      {/* Scientific Context & Glossary */}
      <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 space-y-12">
        <div className="space-y-6">
          <h4 className="text-2xl font-black uppercase italic tracking-tighter text-blue-400">Interpretación Científica del Índice CAI</h4>
          <p className="text-sm font-medium text-slate-300 leading-relaxed max-w-4xl">
            El <span className="text-white font-black">Codon Adaptation Index (CAI)</span> mide cuán adaptados están los codones de un gen al uso preferido del organismo. 
            Un valor de CAI alto indica que el gen utiliza predominantemente los codones óptimos, lo cual está fuertemente correlacionado con una <span className="text-blue-400 font-black">alta eficiencia de traducción y mayores niveles de expresión proteica</span>.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-12 pt-8 border-t border-white/5">
          <div className="space-y-4">
            <h5 className="text-[10px] font-black text-indigo-400 uppercase tracking-[0.3em]">Metodología de Referencia</h5>
            <p className="text-xs text-slate-400 leading-relaxed font-medium uppercase tracking-tight">
              Para el cálculo, se utilizan las <span className="text-white">proteínas ribosomales y factores de elongación</span> como conjunto de referencia, 
              ya que son genes constitutivamente muy expresados en todos los organismos. 
              <br/><br/>
              <span className="font-mono text-blue-400">w(codon) = freq / max_freq_sinónimo</span>
              <br/>
              <span className="font-mono text-blue-400">CAI = media geométrica de los pesos (w)</span>
            </p>
          </div>
          <div className="space-y-4">
            <h5 className="text-[10px] font-black text-rose-400 uppercase tracking-[0.3em]">Importancia en Genómica Comparativa</h5>
            <p className="text-xs text-slate-400 leading-relaxed font-medium uppercase tracking-tight">
              Los genes con un <span className="text-rose-400 font-black">CAI significativamente bajo</span> y un uso atípico de codones son candidatos principales para haber sido adquiridos mediante <span className="text-white">Transferencia Horizontal de Genes (HGT)</span> o podrían representar pseudogenes en proceso de degradación.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
