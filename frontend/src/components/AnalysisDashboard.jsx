/**
 * AnalysisDashboard Component
 * Main dashboard showing analysis overview and key metrics
 */

const STATUS_COLORS = {
  PASS: 'bg-gradient-to-r from-green-50 to-emerald-50 text-green-800 border border-green-200',
  WARNING: 'bg-gradient-to-r from-yellow-50 to-amber-50 text-yellow-800 border border-yellow-200',
  FAIL: 'bg-gradient-to-r from-red-50 to-rose-50 text-red-800 border border-red-200',
}

const STATUS_BADGE = {
  PASS: 'bg-green-500',
  WARNING: 'bg-yellow-500',
  FAIL: 'bg-red-500',
}

function StatCard({ title, value, subtitle, color = 'emerald' }) {
  const colorClasses = {
    emerald: {
      gradient: 'from-emerald-500/10 to-emerald-600/5',
      border: 'border-emerald-200',
      text: 'text-emerald-600',
      badge: 'bg-emerald-500'
    },
    blue: {
      gradient: 'from-blue-500/10 to-blue-600/5',
      border: 'border-blue-200',
      text: 'text-blue-600',
      badge: 'bg-blue-500'
    },
    purple: {
      gradient: 'from-purple-500/10 to-purple-600/5',
      border: 'border-purple-200',
      text: 'text-purple-600',
      badge: 'bg-purple-500'
    },
    orange: {
      gradient: 'from-orange-500/10 to-orange-600/5',
      border: 'border-orange-200',
      text: 'text-orange-600',
      badge: 'bg-orange-500'
    },
    cyan: {
      gradient: 'from-cyan-500/10 to-cyan-600/5',
      border: 'border-cyan-200',
      text: 'text-cyan-600',
      badge: 'bg-cyan-500'
    },
  }

  const colors = colorClasses[color]

  return (
    <div className={`group relative bg-gradient-to-br ${colors.gradient} backdrop-blur-sm rounded-2xl border ${colors.border} p-6 hover:shadow-lg hover:scale-[1.02] transition-all duration-300 overflow-hidden`}>
      <div className={`absolute top-0 right-0 w-24 h-24 ${colors.badge} opacity-5 rounded-full -mr-12 -mt-12 group-hover:scale-150 transition-transform duration-500`}></div>
      <div className="relative">
        <div className={`inline-block px-3 py-1 rounded-full text-xs font-semibold ${colors.text} bg-white/60 backdrop-blur-sm mb-3`}>
          {title}
        </div>
        <p className={`text-3xl font-bold ${colors.text} mb-1`}>{value}</p>
        {subtitle && (
          <p className="text-sm text-gray-500 font-medium">{subtitle}</p>
        )}
      </div>
    </div>
  )
}

function formatNumber(num) {
  if (num >= 1000000) {
    return (num / 1000000).toFixed(2) + 'M'
  }
  if (num >= 1000) {
    return (num / 1000).toFixed(1) + 'K'
  }
  return num?.toString() || '0'
}

export default function AnalysisDashboard({ analysisData, isLoading, status }) {
  if (status === 'idle' && !analysisData) {
    return (
      <div className="text-center py-20">
        <div className="max-w-2xl mx-auto">
          <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-gradient-to-br from-emerald-100 to-cyan-100 mb-8 animate-pulse">
            <div className="w-12 h-12 rounded-full bg-gradient-to-br from-emerald-500 to-cyan-500"></div>
          </div>
          <h2 className="text-3xl font-bold text-gray-800 mb-4">
            Bienvenido al Análisis Genómico
          </h2>
          <p className="text-gray-600 text-lg leading-relaxed">
            Haga clic en <span className="font-semibold text-emerald-600">"Ejecutar Análisis"</span> para comenzar el análisis del genoma de E. coli K-12 MG1655.
            Se analizarán codones, genes y se validarán los resultados contra valores de referencia.
          </p>
        </div>
      </div>
    )
  }

  if (isLoading) {
    return (
      <div className="text-center py-20">
        <div className="relative inline-flex items-center justify-center mb-8">
          <div className="absolute w-32 h-32 rounded-full bg-emerald-200 animate-ping opacity-20"></div>
          <div className="relative w-24 h-24 rounded-full border-4 border-emerald-200 border-t-emerald-600 animate-spin"></div>
        </div>
        <h2 className="text-3xl font-bold text-gray-800 mb-3">
          Analizando genoma...
        </h2>
        <p className="text-gray-600 text-lg">
          Procesando datos genómicos
        </p>
        <div className="mt-6 flex items-center justify-center gap-2">
          <div className="w-2 h-2 bg-emerald-500 rounded-full animate-bounce" style={{ animationDelay: '0ms' }}></div>
          <div className="w-2 h-2 bg-emerald-500 rounded-full animate-bounce" style={{ animationDelay: '150ms' }}></div>
          <div className="w-2 h-2 bg-emerald-500 rounded-full animate-bounce" style={{ animationDelay: '300ms' }}></div>
        </div>
      </div>
    )
  }

  if (!analysisData) {
    return (
      <div className="text-center py-20">
        <div className="max-w-md mx-auto">
          <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-gradient-to-br from-red-100 to-rose-100 mb-8">
            <div className="w-12 h-12 rounded-full bg-gradient-to-br from-red-500 to-rose-500"></div>
          </div>
          <h2 className="text-3xl font-bold text-gray-800 mb-4">
            Error en el análisis
          </h2>
          <p className="text-gray-600 text-lg">
            Ocurrió un error durante el análisis. Por favor, inténtelo de nuevo.
          </p>
        </div>
      </div>
    )
  }

  const { codons, genes, validation } = analysisData

  return (
    <div className="space-y-6">
      {/* Key Metrics Cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-5">
        <StatCard
          title="Longitud del Genoma"
          value={formatNumber(genes.genome_length)}
          subtitle="pares de bases"
          color="emerald"
        />
        <StatCard
          title="Total de Genes"
          value={formatNumber(genes.total_genes)}
          subtitle={`${genes.total_cds} CDS`}
          color="blue"
        />
        <StatCard
          title="Contenido GC"
          value={`${genes.gc_content}%`}
          subtitle="guanina + citosina"
          color="purple"
        />
        <StatCard
          title="Codones ATG"
          value={formatNumber(codons.atg_count)}
          subtitle={`${codons.atg_density}/kb`}
          color="cyan"
        />
      </div>

      {/* Secondary Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-5">
        <div className="group bg-white rounded-2xl shadow-sm border border-gray-100 p-6 hover:shadow-md hover:border-emerald-200 transition-all duration-300">
          <div className="flex items-center justify-between mb-3">
            <h3 className="font-semibold text-gray-700 text-sm uppercase tracking-wide">Densidad Génica</h3>
            <div className="w-2 h-2 bg-emerald-500 rounded-full group-hover:scale-150 transition-transform duration-300"></div>
          </div>
          <p className="text-4xl font-bold text-emerald-600 mb-2">
            {genes.gene_density}
          </p>
          <p className="text-gray-500 text-sm font-medium">genes por megabase</p>
        </div>

        <div className="group bg-white rounded-2xl shadow-sm border border-gray-100 p-6 hover:shadow-md hover:border-blue-200 transition-all duration-300">
          <div className="flex items-center justify-between mb-3">
            <h3 className="font-semibold text-gray-700 text-sm uppercase tracking-wide">Tamaño Promedio</h3>
            <div className="w-2 h-2 bg-blue-500 rounded-full group-hover:scale-150 transition-transform duration-300"></div>
          </div>
          <p className="text-4xl font-bold text-blue-600 mb-2">
            {genes.size_statistics.mean.toFixed(0)}
          </p>
          <p className="text-gray-500 text-sm font-medium">pares de bases</p>
        </div>

        <div className="group bg-white rounded-2xl shadow-sm border border-gray-100 p-6 hover:shadow-md hover:border-purple-200 transition-all duration-300">
          <div className="flex items-center justify-between mb-3">
            <h3 className="font-semibold text-gray-700 text-sm uppercase tracking-wide">ATG vs Genes</h3>
            <div className="w-2 h-2 bg-purple-500 rounded-full group-hover:scale-150 transition-transform duration-300"></div>
          </div>
          <p className="text-4xl font-bold text-purple-600 mb-2">
            {codons.gene_comparison.difference > 0 ? '+' : ''}{codons.gene_comparison.difference}
          </p>
          <p className="text-gray-500 text-sm font-medium">diferencia (ATG - genes anotados)</p>
        </div>
      </div>

      {/* Stop Codons */}
      <div className="bg-white rounded-2xl shadow-sm border border-gray-100 p-7 hover:shadow-md transition-all duration-300">
        <div className="flex items-center justify-between mb-6">
          <h3 className="font-bold text-gray-800 text-lg">Distribución de Codones de Parada</h3>
          <div className="flex items-center gap-2">
            <div className="w-3 h-3 bg-orange-500 rounded-full"></div>
            <span className="text-sm text-gray-500 font-medium">Stop Codons</span>
          </div>
        </div>
        <div className="grid grid-cols-3 gap-5">
          {Object.entries(codons.stop_codons).map(([codon, data]) => (
            <div key={codon} className="group text-center p-5 bg-gradient-to-br from-orange-50/50 to-amber-50/30 border border-orange-100 rounded-xl hover:shadow-md hover:scale-105 transition-all duration-300">
              <div className="inline-block px-4 py-2 bg-white rounded-lg shadow-sm mb-3 group-hover:shadow-md transition-shadow duration-300">
                <p className="text-2xl font-mono font-bold text-gray-800">{codon}</p>
              </div>
              <p className="text-4xl font-bold text-orange-600 mb-3">
                {formatNumber(data.count)}
              </p>
              <div className="relative mt-3 bg-gray-200 rounded-full h-2.5 overflow-hidden">
                <div
                  className="absolute top-0 left-0 bg-gradient-to-r from-orange-500 to-amber-500 h-full transition-all duration-700 ease-out rounded-full"
                  style={{ width: `${data.percentage}%` }}
                ></div>
              </div>
              <p className="text-sm text-gray-600 font-semibold mt-2">{data.percentage.toFixed(1)}%</p>
            </div>
          ))}
        </div>
      </div>

      {/* Validation Results */}
      {validation && (
        <div className="bg-white rounded-2xl shadow-sm border border-gray-100 p-7 hover:shadow-md transition-all duration-300">
          <div className="flex items-center justify-between mb-6">
            <div className="flex items-center gap-3">
              <div className={`w-3 h-3 rounded-full ${STATUS_BADGE[validation.overall_status]}`}></div>
              <h3 className="font-bold text-gray-800 text-lg">Validación contra Referencia</h3>
            </div>
            <span className={`px-4 py-2 rounded-full text-sm font-semibold ${
              STATUS_COLORS[validation.overall_status]
            }`}>
              {validation.overall_status === 'PASS' && 'Validado'}
              {validation.overall_status === 'WARNING' && 'Con advertencias'}
              {validation.overall_status === 'FAIL' && 'No validado'}
            </span>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full">
              <thead>
                <tr className="text-left text-xs uppercase tracking-wider text-gray-500 border-b-2 border-gray-200">
                  <th className="pb-4 font-semibold">Métrica</th>
                  <th className="pb-4 font-semibold">Calculado</th>
                  <th className="pb-4 font-semibold">Referencia</th>
                  <th className="pb-4 font-semibold">Desviación</th>
                  <th className="pb-4 font-semibold">Estado</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-100">
                {validation.items.map((item, index) => (
                  <tr key={index} className="hover:bg-gray-50 transition-colors duration-150">
                    <td className="py-4 font-semibold text-gray-800">
                      {item.metric.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                    </td>
                    <td className="py-4 text-gray-700 font-medium">
                      {typeof item.calculated === 'number'
                        ? item.calculated.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.calculated}
                    </td>
                    <td className="py-4 text-gray-700 font-medium">
                      {typeof item.reference === 'number'
                        ? item.reference.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.reference}
                    </td>
                    <td className="py-4">
                      <span className={`inline-flex items-center px-3 py-1 rounded-full text-sm font-bold ${
                        item.deviation_percent <= 5 ? 'bg-green-100 text-green-700' :
                        item.deviation_percent <= 10 ? 'bg-yellow-100 text-yellow-700' : 'bg-red-100 text-red-700'
                      }`}>
                        {item.deviation_percent.toFixed(2)}%
                      </span>
                    </td>
                    <td className="py-4">
                      <div className="flex items-center gap-2">
                        <div className={`w-2 h-2 rounded-full ${STATUS_BADGE[item.status]}`}></div>
                        <span className={`px-3 py-1 rounded-full text-xs font-semibold ${
                          STATUS_COLORS[item.status]
                        }`}>
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
      )}

      {/* Info Box */}
      <div className="relative bg-gradient-to-r from-blue-50 to-cyan-50 border-l-4 border-blue-500 rounded-2xl p-6 shadow-sm hover:shadow-md transition-all duration-300">
        <div className="flex items-start gap-4">
          <div className="flex-shrink-0 mt-1">
            <div className="w-8 h-8 rounded-full bg-blue-500 flex items-center justify-center">
              <span className="text-white font-bold text-sm">i</span>
            </div>
          </div>
          <div>
            <h4 className="text-blue-900 font-bold text-sm mb-2 uppercase tracking-wide">Información de Referencia</h4>
            <p className="text-sm text-blue-800 leading-relaxed">
              <strong className="font-bold">E. coli K-12 MG1655</strong> es la cepa de referencia de Escherichia coli,
              con un genoma de aproximadamente 4.6 millones de pares de bases y ~4,300 genes.
              Los valores de referencia corresponden a la secuencia RefSeq <span className="font-mono font-semibold">NC_000913.3</span>.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
