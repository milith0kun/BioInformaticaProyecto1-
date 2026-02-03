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

function StatCard({ title, value, subtitle, variant = 'teal' }) {
  const variants = {
    teal: 'border-teal-100 bg-gradient-to-br from-teal-50 to-white',
    emerald: 'border-emerald-100 bg-gradient-to-br from-emerald-50 to-white',
    slate: 'border-slate-100 bg-gradient-to-br from-slate-50 to-white',
  }

  const textColors = {
    teal: 'text-teal-700',
    emerald: 'text-emerald-700',
    slate: 'text-slate-700',
  }

  return (
    <div className={`rounded-xl border p-4 sm:p-5 ${variants[variant]} transition-all hover:shadow-md`}>
      <p className="text-xs uppercase tracking-wide text-slate-500 font-medium mb-1">{title}</p>
      <p className={`text-2xl sm:text-3xl font-bold ${textColors[variant]}`}>{value}</p>
      {subtitle && <p className="text-xs sm:text-sm text-slate-500 mt-1">{subtitle}</p>}
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
      <div className="text-center py-12 sm:py-20 px-4">
        <div className="max-w-xl mx-auto">
          <div className="w-16 h-16 sm:w-20 sm:h-20 mx-auto bg-teal-50 rounded-2xl flex items-center justify-center mb-4 sm:mb-6">
            <svg className="w-8 h-8 sm:w-10 sm:h-10 text-teal-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
            </svg>
          </div>
          <h2 className="text-xl sm:text-2xl font-bold text-slate-800 mb-2 sm:mb-3">
            Bienvenido al Análisis Genómico
          </h2>
          <p className="text-sm sm:text-base text-slate-600 leading-relaxed">
            Haga clic en <span className="font-semibold text-teal-600">"Ejecutar Análisis"</span> para comenzar.
            Se analizarán codones, genes y se validarán los resultados contra valores de referencia.
          </p>
        </div>
      </div>
    )
  }

  if (isLoading) {
    return (
      <div className="text-center py-12 sm:py-20 px-4">
        <div className="w-12 h-12 sm:w-16 sm:h-16 mx-auto border-4 border-teal-200 border-t-teal-600 rounded-full animate-spin mb-4 sm:mb-6"></div>
        <h2 className="text-lg sm:text-xl font-bold text-slate-800 mb-2">Analizando genoma...</h2>
        <p className="text-sm sm:text-base text-slate-500">Procesando datos genómicos</p>
      </div>
    )
  }

  if (!analysisData) {
    return (
      <div className="text-center py-12 sm:py-20 px-4">
        <div className="w-16 h-16 sm:w-20 sm:h-20 mx-auto bg-red-50 rounded-2xl flex items-center justify-center mb-4 sm:mb-6">
          <svg className="w-8 h-8 sm:w-10 sm:h-10 text-red-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
        </div>
        <h2 className="text-lg sm:text-xl font-bold text-slate-800 mb-2">Error en el análisis</h2>
        <p className="text-sm sm:text-base text-slate-500">Ocurrió un error. Por favor, inténtelo de nuevo.</p>
      </div>
    )
  }

  const { codons, genes, validation } = analysisData

  return (
    <div className="space-y-6">
      {/* Key Metrics */}
      <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
        <StatCard
          title="Longitud del Genoma"
          value={formatNumber(genes.genome_length)}
          subtitle="pares de bases"
          variant="teal"
        />
        <StatCard
          title="Total de Genes"
          value={formatNumber(genes.total_genes)}
          subtitle={`${genes.total_cds} CDS`}
          variant="emerald"
        />
        <StatCard
          title="Contenido GC"
          value={`${genes.gc_content}%`}
          subtitle="guanina + citosina"
          variant="teal"
        />
        <StatCard
          title="Codones ATG"
          value={formatNumber(codons.atg_count)}
          subtitle={`${codons.atg_density}/kb`}
          variant="emerald"
        />
      </div>

      {/* Secondary Metrics */}
      <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
        <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-5">
          <p className="text-xs uppercase tracking-wide text-slate-500 font-medium">Densidad Génica</p>
          <p className="text-2xl sm:text-3xl font-bold text-teal-700 mt-1">{genes.gene_density}</p>
          <p className="text-xs sm:text-sm text-slate-500">genes por megabase</p>
        </div>
        <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-5">
          <p className="text-xs uppercase tracking-wide text-slate-500 font-medium">Tamaño Promedio</p>
          <p className="text-2xl sm:text-3xl font-bold text-emerald-700 mt-1">{genes.size_statistics.mean.toFixed(0)}</p>
          <p className="text-xs sm:text-sm text-slate-500">pares de bases</p>
        </div>
        <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-5">
          <p className="text-xs uppercase tracking-wide text-slate-500 font-medium">ATG vs Genes</p>
          <p className="text-2xl sm:text-3xl font-bold text-slate-700 mt-1">
            {codons.gene_comparison.difference > 0 ? '+' : ''}{codons.gene_comparison.difference}
          </p>
          <p className="text-xs sm:text-sm text-slate-500">diferencia</p>
        </div>
      </div>

      {/* Stop Codons */}
      <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-6">
        <h3 className="font-semibold text-slate-800 mb-4 text-sm sm:text-base">Distribución de Codones de Parada</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 sm:gap-4">
          {Object.entries(codons.stop_codons).map(([codon, data]) => (
            <div key={codon} className="text-center p-4 bg-slate-50 rounded-xl">
              <span className="inline-block px-3 py-1 bg-white rounded-lg shadow-sm font-mono font-bold text-base sm:text-lg text-slate-800 mb-2">
                {codon}
              </span>
              <p className="text-xl sm:text-2xl font-bold text-teal-700">{formatNumber(data.count)}</p>
              <div className="mt-3 bg-slate-200 rounded-full h-2 overflow-hidden">
                <div
                  className="h-full bg-gradient-to-r from-teal-500 to-emerald-500 transition-all"
                  style={{ width: `${data.percentage}%` }}
                />
              </div>
              <p className="text-xs sm:text-sm text-slate-600 mt-1">{data.percentage.toFixed(1)}%</p>
            </div>
          ))}
        </div>
      </div>

      {/* Validation Results */}
      {validation && (
        <div className="bg-white rounded-xl border border-slate-200 p-4 sm:p-6">
          <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3 mb-4">
            <div className="flex items-center gap-2">
              <span className={`w-2 h-2 rounded-full ${STATUS_DOT[validation.overall_status]}`}></span>
              <h3 className="font-semibold text-slate-800 text-sm sm:text-base">Validación contra Referencia</h3>
            </div>
            <span className={`px-3 py-1 rounded-full text-xs sm:text-sm font-medium ${STATUS_COLORS[validation.overall_status]} w-fit`}>
              {validation.overall_status === 'PASS' && 'Validado'}
              {validation.overall_status === 'WARNING' && 'Con advertencias'}
              {validation.overall_status === 'FAIL' && 'No validado'}
            </span>
          </div>

          <div className="overflow-x-auto -mx-4 sm:mx-0 px-4 sm:px-0">
            <table className="w-full text-xs sm:text-sm min-w-[600px]">
              <thead>
                <tr className="text-left text-xs uppercase tracking-wide text-slate-500 border-b border-slate-200">
                  <th className="pb-3 font-medium whitespace-nowrap pr-2">Métrica</th>
                  <th className="pb-3 font-medium whitespace-nowrap px-2">Calculado</th>
                  <th className="pb-3 font-medium whitespace-nowrap px-2">Referencia</th>
                  <th className="pb-3 font-medium whitespace-nowrap px-2">Desviación</th>
                  <th className="pb-3 font-medium whitespace-nowrap pl-2">Estado</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-100">
                {validation.items.map((item, index) => (
                  <tr key={index} className="hover:bg-slate-50">
                    <td className="py-3 font-medium text-slate-800 pr-2">
                      <span className="line-clamp-2">
                        {item.metric.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                      </span>
                    </td>
                    <td className="py-3 text-slate-600 px-2">
                      {typeof item.calculated === 'number'
                        ? item.calculated.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.calculated}
                    </td>
                    <td className="py-3 text-slate-600 px-2">
                      {typeof item.reference === 'number'
                        ? item.reference.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.reference}
                    </td>
                    <td className="py-3 px-2">
                      <span className={`inline-block px-2 py-0.5 rounded text-xs font-medium whitespace-nowrap ${item.deviation_percent <= 5 ? 'bg-emerald-100 text-emerald-700' :
                          item.deviation_percent <= 10 ? 'bg-amber-100 text-amber-700' : 'bg-red-100 text-red-700'
                        }`}>
                        {item.deviation_percent.toFixed(2)}%
                      </span>
                    </td>
                    <td className="py-3 pl-2">
                      <div className="flex items-center gap-2">
                        <span className={`w-1.5 h-1.5 rounded-full ${STATUS_DOT[item.status]}`}></span>
                        <span className="text-xs text-slate-600 whitespace-nowrap">{item.status}</span>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Info Box */}
      <div className="bg-teal-50 border border-teal-100 rounded-xl p-4 sm:p-5">
        <div className="flex gap-3">
          <div className="w-7 h-7 sm:w-8 sm:h-8 bg-teal-500 rounded-lg flex items-center justify-center flex-shrink-0">
            <svg className="w-3.5 h-3.5 sm:w-4 sm:h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="flex-1 min-w-0">
            <h4 className="font-medium text-teal-800 mb-1 text-xs sm:text-sm">Información de Referencia</h4>
            <p className="text-xs sm:text-sm text-teal-700 leading-relaxed">
              <strong>E. coli K-12 MG1655</strong> es la cepa de referencia de Escherichia coli,
              con un genoma de aproximadamente 4.6 millones de pares de bases y ~4,300 genes.
              Los valores de referencia corresponden a la secuencia RefSeq <span className="font-mono">NC_000913.3</span>.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
