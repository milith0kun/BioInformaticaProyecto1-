/**
 * CodonVisualization Component
 * Displays codon analysis with interactive charts
 */
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  PieChart,
  Pie,
  Cell,
  Legend
} from 'recharts'

const COLORS = ['#14b8a6', '#f59e0b', '#ef4444', '#8b5cf6']

export default function CodonVisualization({ codonData }) {
  if (!codonData) {
    return (
      <div className="text-center py-20">
        <div className="w-20 h-20 mx-auto bg-slate-100 rounded-2xl flex items-center justify-center mb-6">
          <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
          </svg>
        </div>
        <h2 className="text-xl font-bold text-slate-700 mb-2">Sin datos de codones</h2>
        <p className="text-slate-500">Ejecute el análisis completo para ver la visualización</p>
      </div>
    )
  }

  const codonComparisonData = [
    { name: 'ATG', count: codonData.atg_count, fill: '#14b8a6' },
    { name: 'TAA', count: codonData.stop_codons.TAA.count, fill: '#f59e0b' },
    { name: 'TAG', count: codonData.stop_codons.TAG.count, fill: '#ef4444' },
    { name: 'TGA', count: codonData.stop_codons.TGA.count, fill: '#8b5cf6' },
  ]

  const stopCodonPieData = Object.entries(codonData.stop_codons).map(([codon, data]) => ({
    name: codon,
    value: data.count,
    percentage: data.percentage
  }))

  const totalStopCodons = Object.values(codonData.stop_codons)
    .reduce((sum, data) => sum + data.count, 0)

  return (
    <div className="space-y-6">
      {/* Stats Cards */}
      <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
        <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-5 text-white">
          <p className="text-teal-100 text-sm">Longitud del Genoma</p>
          <p className="text-2xl font-bold mt-1">
            {(codonData.genome_length / 1000000).toFixed(2)}M
          </p>
          <p className="text-teal-200 text-xs">pares de bases</p>
        </div>
        <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-5 text-white">
          <p className="text-emerald-100 text-sm">Codones ATG</p>
          <p className="text-2xl font-bold mt-1">
            {codonData.atg_count.toLocaleString()}
          </p>
          <p className="text-emerald-200 text-xs">codones de inicio</p>
        </div>
        <div className="bg-gradient-to-br from-slate-600 to-slate-700 rounded-xl p-5 text-white">
          <p className="text-slate-300 text-sm">Densidad ATG</p>
          <p className="text-2xl font-bold mt-1">{codonData.atg_density}</p>
          <p className="text-slate-400 text-xs">por kilobase</p>
        </div>
        <div className="bg-gradient-to-br from-amber-500 to-amber-600 rounded-xl p-5 text-white">
          <p className="text-amber-100 text-sm">Total Stop Codons</p>
          <p className="text-2xl font-bold mt-1">{totalStopCodons.toLocaleString()}</p>
          <p className="text-amber-200 text-xs">TAA + TAG + TGA</p>
        </div>
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Bar Chart */}
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">Comparación de Codones</h3>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={codonComparisonData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
              <XAxis dataKey="name" tick={{ fontSize: 12, fill: '#64748b' }} />
              <YAxis
                tick={{ fontSize: 12, fill: '#64748b' }}
                tickFormatter={(value) => value >= 1000 ? `${(value / 1000).toFixed(0)}K` : value}
              />
              <Tooltip
                formatter={(value) => [value.toLocaleString(), 'Conteo']}
                contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
              />
              <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                {codonComparisonData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={entry.fill} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Pie Chart */}
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <h3 className="font-semibold text-slate-800 mb-4">Distribución de Codones de Parada</h3>
          <ResponsiveContainer width="100%" height={280}>
            <PieChart>
              <Pie
                data={stopCodonPieData}
                cx="50%"
                cy="50%"
                innerRadius={50}
                outerRadius={90}
                paddingAngle={3}
                dataKey="value"
                label={({ name, percentage }) => `${name}: ${percentage}%`}
              >
                {stopCodonPieData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index + 1]} />
                ))}
              </Pie>
              <Tooltip
                formatter={(value) => [value.toLocaleString(), 'Conteo']}
                contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
              />
              <Legend />
            </PieChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Detailed Table */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <h3 className="font-semibold text-slate-800 mb-4">Análisis Detallado de Codones de Parada</h3>
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="text-left text-xs uppercase tracking-wide text-slate-500 border-b border-slate-200">
                <th className="pb-3 font-medium">Codón</th>
                <th className="pb-3 font-medium">Nombre</th>
                <th className="pb-3 font-medium">Conteo</th>
                <th className="pb-3 font-medium">Porcentaje</th>
                <th className="pb-3 font-medium">Frecuencia</th>
                <th className="pb-3 font-medium">Distribución</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(codonData.stop_codons).map(([codon, data]) => (
                <tr key={codon} className="border-b border-slate-100 last:border-0">
                  <td className="py-3">
                    <span className={`inline-flex px-2 py-0.5 rounded text-sm font-mono font-medium ${codon === 'TAA' ? 'bg-amber-100 text-amber-700' :
                        codon === 'TAG' ? 'bg-red-100 text-red-700' :
                          'bg-purple-100 text-purple-700'
                      }`}>
                      {codon}
                    </span>
                  </td>
                  <td className="py-3 text-slate-600">
                    {codon === 'TAA' && 'Ocre'}
                    {codon === 'TAG' && 'Ámbar'}
                    {codon === 'TGA' && 'Ópalo'}
                  </td>
                  <td className="py-3 font-medium text-slate-800">{data.count.toLocaleString()}</td>
                  <td className="py-3 text-slate-600">{data.percentage}%</td>
                  <td className="py-3 text-slate-600">
                    {((data.count / codonData.genome_length) * 1000).toFixed(3)}/kb
                  </td>
                  <td className="py-3 w-32">
                    <div className="bg-slate-200 rounded-full h-2 overflow-hidden">
                      <div
                        className={`h-full ${codon === 'TAA' ? 'bg-amber-500' :
                            codon === 'TAG' ? 'bg-red-500' : 'bg-purple-500'
                          }`}
                        style={{ width: `${data.percentage}%` }}
                      />
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Gene Comparison */}
      <div className="bg-gradient-to-r from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
        <h3 className="font-medium mb-4">Comparación: ATG vs Genes Anotados</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 text-center">
          <div>
            <p className="text-teal-100 text-sm">Codones ATG</p>
            <p className="text-3xl font-bold mt-1">{codonData.gene_comparison.atg_found.toLocaleString()}</p>
          </div>
          <div>
            <p className="text-teal-100 text-sm">Genes Anotados</p>
            <p className="text-3xl font-bold mt-1">{codonData.gene_comparison.annotated_genes.toLocaleString()}</p>
          </div>
          <div>
            <p className="text-teal-100 text-sm">Diferencia</p>
            <p className="text-3xl font-bold mt-1">
              {codonData.gene_comparison.difference > 0 ? '+' : ''}{codonData.gene_comparison.difference.toLocaleString()}
            </p>
          </div>
        </div>
        <p className="mt-4 text-teal-100 text-sm">
          No todos los ATG son codones de inicio. Existen ATG internos dentro de secuencias codificantes.
        </p>
      </div>

      {/* Info Box */}
      <div className="bg-amber-50 border border-amber-100 rounded-xl p-4">
        <div className="flex gap-3">
          <div className="w-8 h-8 bg-amber-500 rounded-lg flex items-center justify-center flex-shrink-0">
            <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
            </svg>
          </div>
          <p className="text-sm text-amber-800">
            <strong>TAA (Ocre)</strong> es el codón de parada más común en E. coli (~63%).
            <strong> TAG (Ámbar)</strong> es menos frecuente (~7%), y
            <strong> TGA (Ópalo)</strong> representa ~30%.
          </p>
        </div>
      </div>
    </div>
  )
}
