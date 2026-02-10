/**
 * ComparisonResults Component — Clean Laboratory Edition
 * Comprehensive multi-genome analysis with high-fidelity charts and metrics
 */
import { useMemo } from 'react'
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

const COLORS = ['#2563eb', '#4f46e5', '#7c3aed', '#db2777', '#0f172a', '#64748b', '#3b82f6', '#6366f1']

const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900/95 backdrop-blur-xl border border-white/10 p-4 rounded-2xl shadow-2xl animate-in fade-in zoom-in-95">
        <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">{label}</p>
        {payload.map((entry, index) => (
          <div key={index} className="flex items-center gap-3">
            <div className="w-1.5 h-1.5 rounded-full" style={{ backgroundColor: entry.color }}></div>
            <p className="text-xs font-bold text-white uppercase tracking-tight">
              {entry.name}: <span className="text-blue-200">{entry.value}</span>
            </p>
          </div>
        ))}
      </div>
    )
  }
  return null
}

const TechnicalChart = ({ title, data, dataKey, layout = "vertical", domain }) => (
  <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm transition-all hover:border-blue-200">
    <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-10 text-center">{title}</h4>
    <ResponsiveContainer width="100%" height={280}>
      <BarChart data={data} layout={layout}>
        <CartesianGrid strokeDasharray="3 3" vertical={layout === 'horizontal'} horizontal={layout === 'vertical'} stroke="#f1f5f9" />
        <XAxis 
          type={layout === 'vertical' ? 'number' : 'category'} 
          dataKey={layout === 'vertical' ? undefined : 'name'}
          axisLine={false} 
          tickLine={false} 
          tick={{ fontSize: 9, fontWeight: 700, fill: '#94a3b8' }} 
        />
        <YAxis 
          type={layout === 'vertical' ? 'category' : 'number'}
          dataKey={layout === 'vertical' ? 'name' : undefined}
          axisLine={false} 
          tickLine={false} 
          tick={{ fontSize: 9, fontWeight: 900, fill: '#1e293b' }} 
          width={layout === 'vertical' ? 100 : 40}
          domain={domain}
        />
        <Tooltip content={<CustomTooltip />} cursor={{ fill: '#f8fafc' }} />
        <Bar dataKey={dataKey} radius={layout === 'vertical' ? [0, 8, 8, 0] : [8, 8, 0, 0]} barSize={layout === 'vertical' ? 12 : 24}>
          {data.map((entry, index) => (
            <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
          ))}
        </Bar>
      </BarChart>
    </ResponsiveContainer>
  </div>
)

export default function ComparisonResults({ comparisonResult }) {
  const chartData = useMemo(() => {
    if (!comparisonResult?.genomes) return []
    return comparisonResult.genomes.map((g) => ({
      name: g.organism_name?.split(' ').slice(0, 2).join(' ') || g.accession?.substring(0, 10),
      full_name: g.organism_name,
      'Tamaño (Mb)': parseFloat(((g.genome_length || 0) / 1e6).toFixed(2)),
      'Genes': g.gene_count || 0,
      'GC%': parseFloat(g.gc_content?.toFixed(1) || 0),
      'Densidad': parseFloat(g.gene_density?.toFixed(1) || 0),
    }))
  }, [comparisonResult])

  if (!comparisonResult) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Sincronizando Dataset...</p>
      </div>
    )
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Executive Summary Header */}
      <div className="bg-slate-900 rounded-[3rem] p-10 text-white relative overflow-hidden shadow-2xl shadow-blue-900/20">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32"></div>
        <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 relative z-10">
          <div className="space-y-3">
            <h2 className="text-4xl font-black italic uppercase tracking-tighter text-blue-400">Análisis Comparativo</h2>
            <div className="flex items-center gap-4">
              <span className="px-4 py-1.5 bg-blue-500/10 text-blue-400 text-[10px] font-black uppercase tracking-widest rounded-full border border-blue-500/20">
                {comparisonResult.total_genomes_compared} Cepas Analizadas
              </span>
              <span className="text-[10px] font-bold text-slate-500 uppercase tracking-widest">
                {comparisonResult.comparison_date ? new Date(comparisonResult.comparison_date).toLocaleString() : 'Live Session'}
              </span>
            </div>
          </div>
          <div className="flex items-center gap-6">
            <div className="text-right">
              <div className="text-6xl font-black tracking-tighter text-white">{comparisonResult.total_genomes_compared}</div>
              <div className="text-[10px] font-black uppercase tracking-widest text-blue-500">Genomas</div>
            </div>
          </div>
        </div>
      </div>

      {/* Extreme Metrics Grid */}
      <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-4">
        {[
          { label: 'Más Grande', val: comparisonResult.extremes.largest_genome, icon: '🐘' },
          { label: 'Más Pequeño', val: comparisonResult.extremes.smallest_genome, icon: '🔬' },
          { label: 'Mayor GC%', val: comparisonResult.extremes.highest_gc, icon: '🔝' },
          { label: 'Menor GC%', val: comparisonResult.extremes.lowest_gc, icon: '📉' },
          { label: 'Máx Densidad', val: comparisonResult.extremes.highest_gene_density, icon: '🧬' },
          { label: 'Mín Densidad', val: comparisonResult.extremes.lowest_gene_density, icon: '🌫️' }
        ].map((item, i) => (
          <div key={i} className="bg-white border-2 border-slate-100 rounded-[2rem] p-6 shadow-sm hover:border-blue-200 transition-all group overflow-hidden relative">
            <div className="absolute -right-2 -bottom-2 text-4xl opacity-5 grayscale group-hover:grayscale-0 group-hover:scale-110 transition-all duration-500">{item.icon}</div>
            <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest mb-3 leading-none">{item.label}</p>
            <p className="text-xs font-black text-slate-900 truncate tracking-tight uppercase italic">{item.val}</p>
          </div>
        ))}
      </div>

      {/* Detailed Matrix Table */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
        <div className="p-10 border-b border-slate-50 bg-slate-50/30 flex items-center justify-between">
          <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">📊 Comparación Detallada</h3>
          <div className="flex gap-2">
            <div className="w-2 h-2 bg-blue-500 rounded-full animate-pulse"></div>
            <div className="w-2 h-2 bg-indigo-500 rounded-full animate-pulse [animation-delay:200ms]"></div>
          </div>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left border-collapse">
            <thead className="bg-white border-b border-slate-100">
              <tr>
                <th className="px-10 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest">Organismo</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Tamaño</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Genes</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">GC%</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Densidad</th>
                <th className="px-10 py-6 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Gen Mayor / Menor</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-50">
              {comparisonResult.genomes?.map((g, i) => (
                <tr key={g.accession} className="hover:bg-slate-50/80 transition-colors group">
                  <td className="px-10 py-6">
                    <div className="flex items-center gap-5">
                      <div className="w-10 h-10 rounded-2xl flex items-center justify-center text-white font-black text-xs shadow-lg group-hover:rotate-6 transition-all" style={{ backgroundColor: COLORS[i % COLORS.length] }}>
                        {g.organism_name?.charAt(0) || 'G'}
                      </div>
                      <div>
                        <p className="font-black text-slate-900 text-xs uppercase tracking-tighter leading-tight max-w-[200px] line-clamp-1">{g.organism_name}</p>
                        <p className="text-[10px] text-blue-600 font-mono tracking-widest font-black mt-1 uppercase">{g.accession}</p>
                      </div>
                    </div>
                  </td>
                  <td className="px-6 py-6 text-right font-mono text-xs font-black text-slate-600 italic">
                    {((g.genome_length || 0) / 1e6).toFixed(2)} <span className="text-[9px] text-slate-400">Mb</span>
                  </td>
                  <td className="px-6 py-6 text-right font-mono text-xs font-black text-slate-600">
                    {g.gene_count?.toLocaleString()}
                  </td>
                  <td className="px-6 py-6 text-right font-mono text-xs font-black text-blue-600">
                    {g.gc_content?.toFixed(1)}%
                  </td>
                  <td className="px-6 py-6 text-right font-mono text-xs font-black text-slate-400">
                    {g.gene_density?.toFixed(1)} <span className="text-[9px]">g/Mb</span>
                  </td>
                  <td className="px-10 py-6 text-right">
                    <div className="flex flex-col items-end gap-1">
                      <span className="text-[10px] font-black text-slate-700">{g.max_gene_length?.toLocaleString()} pb</span>
                      <span className="text-[9px] font-bold text-slate-400 uppercase tracking-tighter border-t border-slate-100 pt-1">{g.min_gene_length?.toLocaleString()} pb</span>
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Multi-Grid Analytics Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        <TechnicalChart title="Extensión Genómica (Mb)" data={chartData} dataKey="Tamaño (Mb)" layout="vertical" />
        <TechnicalChart title="Número de Genes (CDS)" data={chartData} dataKey="Genes" layout="vertical" />
        <TechnicalChart title="Contenido GC (%)" data={chartData} dataKey="GC%" layout="horizontal" domain={['dataMin - 5', 'dataMax + 5']} />
        <TechnicalChart title="Densidad Génica (g/Mb)" data={chartData} dataKey="Densidad" layout="horizontal" />
      </div>

      {/* Global Extremes: Longest and Shortest Genes */}
      <div className="grid grid-cols-1 xl:grid-cols-2 gap-10">
        {/* Longest Genes Global */}
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
          <div className="p-10 border-b border-slate-100 bg-gradient-to-br from-slate-900 to-blue-900 text-white relative">
            <div className="absolute top-0 right-0 p-8 opacity-10"><span className="text-6xl text-white">📏</span></div>
            <h4 className="text-[10px] font-black uppercase tracking-[0.3em] mb-2 text-blue-400">Genes Más Largos (Global)</h4>
            <p className="text-xs text-slate-400 font-bold uppercase tracking-widest">Top 10 Macro-estructuras detectadas</p>
          </div>
          <div className="divide-y divide-slate-50 max-h-[500px] overflow-y-auto custom-scrollbar">
            {comparisonResult.longest_genes_global?.map((gene, i) => (
              <div key={i} className="p-8 hover:bg-blue-50/30 transition-all group flex items-center justify-between gap-6">
                <div className="flex items-center gap-6 flex-1 min-w-0">
                  <div className="text-lg font-black text-slate-200 group-hover:text-blue-200 transition-colors">{(i + 1).toString().padStart(2, '0')}</div>
                  <div className="min-w-0">
                    <p className="font-mono text-sm font-black text-slate-900 group-hover:text-blue-600 transition-colors truncate uppercase italic">{gene.gene_name || gene.locus_tag || 'Unknown'}</p>
                    <p className="text-[10px] text-slate-400 font-bold uppercase tracking-tight mt-1 truncate">{gene.product || 'Proteína funcional'}</p>
                    <div className="flex items-center gap-2 mt-2">
                      <span className="text-[9px] font-black text-blue-500 uppercase px-2 py-0.5 bg-blue-50 rounded border border-blue-100">{gene.genome}</span>
                    </div>
                  </div>
                </div>
                <div className="text-right flex-shrink-0">
                  <p className="text-lg font-black text-slate-900 tracking-tighter italic">{gene.length?.toLocaleString()} <span className="text-[10px] text-slate-400 not-italic uppercase font-bold">pb</span></p>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Shortest Genes Global */}
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
          <div className="p-10 border-b border-slate-100 bg-slate-50 relative">
            <div className="absolute top-0 right-0 p-8 opacity-10"><span className="text-6xl text-slate-900">🔬</span></div>
            <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-2">Genes Más Cortos (Global)</h4>
            <p className="text-xs text-slate-500 font-bold uppercase tracking-widest">Top 10 Micro-secuencias codificantes</p>
          </div>
          <div className="divide-y divide-slate-50 max-h-[500px] overflow-y-auto custom-scrollbar">
            {comparisonResult.shortest_genes_global?.map((gene, i) => (
              <div key={i} className="p-8 hover:bg-slate-50 transition-all group flex items-center justify-between gap-6">
                <div className="flex items-center gap-6 flex-1 min-w-0">
                  <div className="text-lg font-black text-slate-200 group-hover:text-slate-400 transition-colors">{(i + 1).toString().padStart(2, '0')}</div>
                  <div className="min-w-0">
                    <p className="font-mono text-sm font-black text-slate-900 group-hover:text-slate-700 transition-colors truncate uppercase italic">{gene.gene_name || gene.locus_tag || 'Unknown'}</p>
                    <p className="text-[10px] text-slate-400 font-bold uppercase tracking-tight mt-1 truncate">{gene.product || 'Micro-proteína hipotética'}</p>
                    <div className="flex items-center gap-2 mt-2">
                      <span className="text-[9px] font-black text-slate-500 uppercase px-2 py-0.5 bg-slate-100 rounded border border-slate-200">{gene.genome}</span>
                    </div>
                  </div>
                </div>
                <div className="text-right flex-shrink-0">
                  <p className="text-lg font-black text-slate-900 tracking-tighter italic">{gene.length?.toLocaleString()} <span className="text-[10px] text-slate-400 not-italic uppercase font-bold">pb</span></p>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  )
}