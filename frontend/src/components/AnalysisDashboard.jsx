/**
 * AnalysisDashboard Component
 * Main dashboard showing analysis overview and key metrics
 */
import { 
  BeakerIcon, 
  ChartBarIcon, 
  CheckCircleIcon,
  ExclamationCircleIcon,
  XCircleIcon 
} from '@heroicons/react/24/outline'

const STATUS_ICONS = {
  PASS: <CheckCircleIcon className="h-5 w-5 text-green-500" />,
  WARNING: <ExclamationCircleIcon className="h-5 w-5 text-yellow-500" />,
  FAIL: <XCircleIcon className="h-5 w-5 text-red-500" />,
}

const STATUS_COLORS = {
  PASS: 'bg-green-100 text-green-800',
  WARNING: 'bg-yellow-100 text-yellow-800',
  FAIL: 'bg-red-100 text-red-800',
}

function StatCard({ title, value, subtitle, icon, color = 'emerald' }) {
  const colorClasses = {
    emerald: 'from-emerald-500 to-emerald-600',
    blue: 'from-blue-500 to-blue-600',
    purple: 'from-purple-500 to-purple-600',
    orange: 'from-orange-500 to-orange-600',
    cyan: 'from-cyan-500 to-cyan-600',
  }

  return (
    <div className="bg-white rounded-xl shadow-sm p-6 hover:shadow-md transition-shadow">
      <div className="flex items-center justify-between">
        <div>
          <p className="text-sm font-medium text-gray-500">{title}</p>
          <p className="text-2xl font-bold text-gray-900 mt-1">{value}</p>
          {subtitle && (
            <p className="text-sm text-gray-400 mt-1">{subtitle}</p>
          )}
        </div>
        <div className={`p-3 rounded-xl bg-gradient-to-br ${colorClasses[color]}`}>
          {icon}
        </div>
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
      <div className="text-center py-16">
        <BeakerIcon className="h-24 w-24 mx-auto text-gray-300 mb-6" />
        <h2 className="text-2xl font-bold text-gray-700 mb-2">
          Bienvenido al Análisis Genómico
        </h2>
        <p className="text-gray-500 max-w-md mx-auto">
          Haga clic en "Ejecutar Análisis" para comenzar el análisis del genoma de E. coli K-12 MG1655.
          Se analizarán codones, genes y se validarán los resultados contra valores de referencia.
        </p>
      </div>
    )
  }

  if (isLoading) {
    return (
      <div className="text-center py-16">
        <div className="relative inline-block">
          <div className="h-24 w-24 rounded-full border-4 border-emerald-200 animate-spin border-t-emerald-600"></div>
          <BeakerIcon className="h-12 w-12 text-emerald-600 absolute top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2" />
        </div>
        <h2 className="text-2xl font-bold text-gray-700 mt-6 mb-2">
          Analizando genoma...
        </h2>
        <p className="text-gray-500">
          Esto puede tomar unos segundos
        </p>
      </div>
    )
  }

  if (!analysisData) {
    return (
      <div className="text-center py-16">
        <XCircleIcon className="h-24 w-24 mx-auto text-red-300 mb-6" />
        <h2 className="text-2xl font-bold text-gray-700 mb-2">
          Error en el análisis
        </h2>
        <p className="text-gray-500">
          Ocurrió un error durante el análisis. Por favor, inténtelo de nuevo.
        </p>
      </div>
    )
  }

  const { codons, genes, validation } = analysisData

  return (
    <div className="space-y-6">
      {/* Key Metrics Cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
        <StatCard
          title="Longitud del Genoma"
          value={formatNumber(genes.genome_length)}
          subtitle="pares de bases"
          icon={<ChartBarIcon className="h-6 w-6 text-white" />}
          color="emerald"
        />
        <StatCard
          title="Total de Genes"
          value={formatNumber(genes.total_genes)}
          subtitle={`${genes.total_cds} CDS`}
          icon={<BeakerIcon className="h-6 w-6 text-white" />}
          color="blue"
        />
        <StatCard
          title="Contenido GC"
          value={`${genes.gc_content}%`}
          subtitle="guanina + citosina"
          icon={<ChartBarIcon className="h-6 w-6 text-white" />}
          color="purple"
        />
        <StatCard
          title="Codones ATG"
          value={formatNumber(codons.atg_count)}
          subtitle={`${codons.atg_density}/kb`}
          icon={<BeakerIcon className="h-6 w-6 text-white" />}
          color="cyan"
        />
      </div>

      {/* Secondary Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">Densidad Génica</h3>
          <p className="text-4xl font-bold text-emerald-600">
            {genes.gene_density}
          </p>
          <p className="text-gray-500 mt-1">genes por megabase</p>
        </div>
        
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">Tamaño Promedio de Gen</h3>
          <p className="text-4xl font-bold text-blue-600">
            {genes.size_statistics.mean.toFixed(0)}
          </p>
          <p className="text-gray-500 mt-1">pares de bases</p>
        </div>

        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">Comparación ATG vs Genes</h3>
          <p className="text-4xl font-bold text-purple-600">
            {codons.gene_comparison.difference > 0 ? '+' : ''}{codons.gene_comparison.difference}
          </p>
          <p className="text-gray-500 mt-1">diferencia (ATG - genes anotados)</p>
        </div>
      </div>

      {/* Stop Codons */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <h3 className="font-semibold text-gray-800 mb-4">Distribución de Codones de Parada</h3>
        <div className="grid grid-cols-3 gap-4">
          {Object.entries(codons.stop_codons).map(([codon, data]) => (
            <div key={codon} className="text-center p-4 bg-gray-50 rounded-lg">
              <p className="text-2xl font-mono font-bold text-gray-800">{codon}</p>
              <p className="text-3xl font-bold text-orange-600 mt-2">
                {formatNumber(data.count)}
              </p>
              <div className="mt-2 bg-gray-200 rounded-full h-2 overflow-hidden">
                <div 
                  className="bg-orange-500 h-full transition-all duration-500"
                  style={{ width: `${data.percentage}%` }}
                ></div>
              </div>
              <p className="text-sm text-gray-500 mt-1">{data.percentage}%</p>
            </div>
          ))}
        </div>
      </div>

      {/* Validation Results */}
      {validation && (
        <div className="bg-white rounded-xl shadow-sm p-6">
          <div className="flex items-center justify-between mb-4">
            <h3 className="font-semibold text-gray-800">Validación contra Referencia</h3>
            <span className={`px-3 py-1 rounded-full text-sm font-medium ${
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
                <tr className="text-left text-sm text-gray-500 border-b">
                  <th className="pb-3 font-medium">Métrica</th>
                  <th className="pb-3 font-medium">Calculado</th>
                  <th className="pb-3 font-medium">Referencia</th>
                  <th className="pb-3 font-medium">Desviación</th>
                  <th className="pb-3 font-medium">Estado</th>
                </tr>
              </thead>
              <tbody>
                {validation.items.map((item, index) => (
                  <tr key={index} className="border-b last:border-b-0">
                    <td className="py-3 font-medium text-gray-800">
                      {item.metric.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                    </td>
                    <td className="py-3 text-gray-600">
                      {typeof item.calculated === 'number' 
                        ? item.calculated.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.calculated}
                    </td>
                    <td className="py-3 text-gray-600">
                      {typeof item.reference === 'number'
                        ? item.reference.toLocaleString('es-ES', { maximumFractionDigits: 2 })
                        : item.reference}
                    </td>
                    <td className="py-3">
                      <span className={`${
                        item.deviation_percent <= 5 ? 'text-green-600' :
                        item.deviation_percent <= 10 ? 'text-yellow-600' : 'text-red-600'
                      }`}>
                        {item.deviation_percent.toFixed(2)}%
                      </span>
                    </td>
                    <td className="py-3">
                      <div className="flex items-center">
                        {STATUS_ICONS[item.status]}
                        <span className={`ml-2 px-2 py-1 rounded-full text-xs font-medium ${
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
      <div className="bg-blue-50 border border-blue-200 rounded-xl p-4">
        <div className="flex">
          <div className="flex-shrink-0">
            <svg className="h-5 w-5 text-blue-400" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="ml-3">
            <p className="text-sm text-blue-700">
              <strong>E. coli K-12 MG1655</strong> es la cepa de referencia de Escherichia coli, 
              con un genoma de aproximadamente 4.6 millones de pares de bases y ~4,300 genes. 
              Los valores de referencia corresponden a la secuencia RefSeq NC_000913.3.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
