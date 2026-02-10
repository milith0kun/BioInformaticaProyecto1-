/**
 * AnalysisDashboard Component
 * Main dashboard showing analysis overview and key metrics
 */

const STATUS_COLORS = {
  PASS: 'bg-emerald-50 text-emerald-700 border border-emerald-200',
  WARNING: 'bg-amber-50 text-amber-700 border border-amber-200',
  FAIL: 'bg-red-50 text-red-700 border border-red-200',
}

const STATUS_DOT = {
  PASS: 'bg-emerald-500',
  WARNING: 'bg-amber-500',
  FAIL: 'bg-red-500',
}

function StatCard({ title, value, subtitle, variant = 'blue' }) {
  const variants = {
    blue: 'border-blue-100 bg-white shadow-sm',
    indigo: 'border-indigo-100 bg-white shadow-sm',
    slate: 'border-slate-100 bg-white shadow-sm',
  }

  const textColors = {
    blue: 'text-blue-600',
    indigo: 'text-indigo-600',
    slate: 'text-slate-700',
  }

  const iconColors = {
    blue: 'bg-blue-50 text-blue-500',
    indigo: 'bg-indigo-50 text-indigo-500',
    slate: 'bg-slate-50 text-slate-500',
  }

  return (
    <div className={`rounded-3xl border-2 p-6 ${variants[variant]} transition-all hover:border-blue-200 hover:shadow-xl group`}>
      <p className="text-[10px] font-black uppercase tracking-[0.2em] text-slate-400 mb-3 group-hover:text-blue-500 transition-colors">{title}</p>
      <div className="flex items-baseline gap-2">
        <p className={`text-3xl font-black tracking-tighter ${textColors[variant]}`}>{value}</p>
        {subtitle && <p className="text-[10px] font-bold text-slate-400 uppercase tracking-tight">{subtitle}</p>}
      </div>
    </div>
  )
}

function formatNumber(num) {
  if (num >= 1000000) return (num / 1000000).toFixed(2) + 'M'
  if (num >= 1000) return (num / 1000).toFixed(1) + 'K'
  return num?.toString() || '0'
}

export default function AnalysisDashboard({ analysisData, isLoading, status }) {
  if (status === 'idle' && !analysisData) {
    return (
      <div className="text-center py-32 px-8">
        <div className="max-w-2xl mx-auto space-y-8">
          <div className="w-24 h-24 mx-auto bg-blue-50 rounded-[2rem] flex items-center justify-center border border-blue-100 shadow-inner">
            <svg className="w-10 h-10 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
            </svg>
          </div>
          <div className="space-y-4">
            <h2 className="text-4xl font-black text-slate-900 tracking-tighter uppercase italic">
              Terminal de Cómputo <span className="text-blue-600">MG1655</span>
            </h2>
            <p className="text-slate-500 font-medium max-w-md mx-auto leading-relaxed uppercase text-[10px] tracking-[0.2em]">
              Entorno listo para procesamiento de secuencias genómicas y validación molecular.
            </p>
          </div>
        </div>
      </div>
    )
  }

  if (isLoading) {
    return (
      <div className="text-center py-48 px-8">
        <div className="w-20 h-20 mx-auto relative mb-10">
          <div className="absolute inset-0 border-4 border-slate-100 rounded-full"></div>
          <div className="absolute inset-0 border-4 border-transparent border-t-blue-600 rounded-full animate-spin"></div>
        </div>
        <h2 className="text-xs font-black text-slate-900 uppercase tracking-[0.5em] animate-pulse">Sincronizando Hebras...</h2>
      </div>
    )
  }

  if (!analysisData) return null

  const { codons, genes, validation } = analysisData

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Principal Metrics Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        <StatCard 
          title="Longitud del Genoma" 
          value={`${(genes.genome_length / 1e6).toFixed(2)}M`} 
          subtitle="pares de bases" 
          variant="blue" 
        />
        <StatCard 
          title="Total de Genes" 
          value={formatNumber(genes.total_genes)} 
          subtitle={`${genes.total_cds} CDS`} 
          variant="indigo" 
        />
        <StatCard 
          title="Contenido GC" 
          value={`${genes.gc_content}%`} 
          subtitle="guanina + citosina" 
          variant="blue" 
        />
        <StatCard 
          title="Codones ATG" 
          value={formatNumber(codons.atg_count)} 
          subtitle={`${(codons.atg_count / (genes.genome_length / 1000)).toFixed(2)}/kb`} 
          variant="indigo" 
        />
      </div>

      {/* Secondary Technical Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black uppercase tracking-widest text-slate-400 mb-4">Densidad Génica</p>
          <div className="flex items-baseline gap-2">
            <p className="text-3xl font-black text-slate-900 tracking-tighter">{genes.gene_density.toFixed(2)}</p>
            <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight">genes / Mb</p>
          </div>
        </div>
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black uppercase tracking-widest text-slate-400 mb-4">Tamaño Promedio</p>
          <div className="flex items-baseline gap-2">
            <p className="text-3xl font-black text-slate-900 tracking-tighter">
              {Math.round(genes.genome_length / (genes.total_genes || 1))}
            </p>
            <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight">pb / gen</p>
          </div>
        </div>
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black uppercase tracking-widest text-slate-400 mb-4">ATG vs Genes</p>
          <div className="flex items-baseline gap-2">
            <p className="text-3xl font-black text-blue-600 tracking-tighter">
              +{codons.atg_count - genes.total_genes}
            </p>
            <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight">diferencia</p>
          </div>
        </div>
      </div>

      {/* Analytics Layout */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        <div className="lg:col-span-2 space-y-8">
          {/* Validation section */}
          <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 overflow-hidden shadow-sm">
            <div className="p-8 border-b border-slate-50 flex items-center justify-between bg-slate-50/30">
              <div className="space-y-1">
                <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">
                  {validation.validation_type === 'multi' ? 'Consenso del Grupo Analizado' : 'Validación contra Referencia'}
                </h3>
                <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest italic">
                  {validation.reference_source || 'Comparativa RefSeq NC_000913.3'}
                </p>
              </div>
              <div className={`px-4 py-1.5 rounded-full text-[9px] font-black uppercase tracking-widest ${STATUS_COLORS[validation.overall_status]}`}>
                {validation.overall_status === 'PASS' ? 'SISTEMA VALIDADO' : 'FUERA DE RANGO'}
              </div>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead className="bg-white border-b border-slate-50">
                  <tr>
                    <th className="px-8 py-4 text-[9px] font-black text-slate-400 uppercase tracking-widest">Métrica</th>
                    <th className="px-8 py-4 text-[9px] font-black text-slate-400 uppercase tracking-widest">Valor</th>
                    <th className="px-8 py-4 text-[9px] font-black text-slate-400 uppercase tracking-widest">
                      {validation.validation_type === 'multi' ? 'Media Grupo' : 'Referencia'}
                    </th>
                    <th className="px-8 py-4 text-[9px] font-black text-slate-400 uppercase tracking-widest">Desviación</th>
                    <th className="px-8 py-4 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Estado</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-50">
                  {validation.items.map((item, i) => (
                    <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                      <td className="px-8 py-5">
                        <p className="text-[10px] font-black text-slate-900 uppercase tracking-tight">{item.metric.replace(/_/g, ' ')}</p>
                      </td>
                      <td className="px-8 py-5">
                        <p className="font-mono text-xs font-black text-blue-600">
                          {(item.calculated || 0).toLocaleString()}
                        </p>
                      </td>
                      <td className="px-8 py-5">
                        <p className="font-mono text-[10px] font-bold text-slate-400">
                          {item.reference?.toLocaleString() || 'N/A'}
                        </p>
                      </td>
                      <td className="px-8 py-5">
                        <p className="font-mono text-[10px] font-black text-slate-500">
                          {item.deviation_percent || item.deviation || '0.00'}%
                        </p>
                      </td>
                      <td className="px-8 py-5 text-right">
                        <div className="flex items-center justify-end gap-2">
                          <div className={`w-1.5 h-1.5 rounded-full ${STATUS_DOT[item.status]}`}></div>
                          <span className={`text-[9px] font-black uppercase tracking-widest ${item.status === 'PASS' ? 'text-emerald-600' : item.status === 'WARNING' ? 'text-amber-600' : 'text-rose-600'}`}>
                            {item.status}
                          </span>
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>

          <div className="p-8 bg-blue-50/50 border-2 border-blue-100 rounded-[2rem] flex gap-6 items-start">
            <div className="w-12 h-12 bg-white rounded-2xl flex items-center justify-center text-xl shadow-sm flex-shrink-0">ℹ️</div>
            <div className="space-y-2">
              <h4 className="text-[10px] font-black text-blue-700 uppercase tracking-widest">Información de Referencia</h4>
              <p className="text-[11px] text-blue-900/70 font-medium leading-relaxed">
                E. coli K-12 MG1655 es la cepa de referencia de Escherichia coli, con un genoma de aproximadamente 4.6 millones de pares de bases y ~4,300 genes. Los valores de referencia corresponden a la secuencia RefSeq NC_000913.3.
              </p>
            </div>
          </div>
        </div>

        {/* Sidebar Info - Stop Codons Distribution */}
        <div className="space-y-6">
          <div className="bg-slate-900 rounded-[3rem] p-10 text-white shadow-2xl shadow-blue-900/20 relative overflow-hidden group">
            <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000"></div>
            <h4 className="text-[10px] font-black uppercase tracking-[0.2em] mb-10 text-slate-400">Distribución de Parada</h4>
            <div className="space-y-10 relative z-10">
              {Object.entries(codons?.stop_codons || {}).map(([codon, data]) => (
                <div key={codon} className="space-y-3">
                  <div className="flex items-center justify-between">
                    <span className="font-mono text-sm font-black tracking-[0.3em] text-blue-400">{codon}</span>
                    <div className="text-right">
                      <span className="text-xs font-black block">{data.count?.toLocaleString()}</span>
                      <span className="text-[9px] font-bold text-slate-500 uppercase">{(data.percentage || 0).toFixed(1)}%</span>
                    </div>
                  </div>
                  <div className="w-full h-1 bg-white/5 rounded-full overflow-hidden">
                    <div 
                      className="h-full bg-blue-500 shadow-[0_0_12px_rgba(59,130,246,0.6)] transition-all duration-1000" 
                      style={{ width: `${data.percentage || 0}%` }}
                    ></div>
                  </div>
                </div>
              ))}
            </div>
          </div>

          <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden">
             <div className="absolute top-0 right-0 p-8 opacity-5">
              <svg className="w-20 h-20" fill="currentColor" viewBox="0 0 24 24"><path d="M13 10V3L4 14h7v7l9-11h-7z"/></svg>
            </div>
            <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] mb-6">Métricas de Espacio</h4>
            <div className="space-y-6 relative z-10">
              <div className="flex justify-between items-end">
                <span className="text-[10px] font-bold text-slate-500 uppercase tracking-widest">Desviación σ</span>
                <span className="text-2xl font-black text-slate-900 tracking-tighter">{codons?.spatial_distribution?.std_dev?.toFixed(3) || '0.000'}</span>
              </div>
              <div className="flex justify-between items-end">
                <span className="text-[10px] font-bold text-slate-500 uppercase tracking-widest">Uniformidad</span>
                <span className="text-2xl font-black text-slate-900 tracking-tighter">{codons?.spatial_distribution?.uniformity_score?.toFixed(1) || '0.0'}%</span>
              </div>
              <p className="text-[10px] text-slate-400 font-medium leading-relaxed uppercase pt-4 border-t border-slate-50">
                Consistencia espacial de tripletas detectada en la hebra MG1655.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
