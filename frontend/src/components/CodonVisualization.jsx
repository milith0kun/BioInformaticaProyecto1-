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
  Legend, 
  ResponsiveContainer,
  PieChart,
  Pie,
  Cell
} from 'recharts'
import { BeakerIcon } from '@heroicons/react/24/outline'

const COLORS = ['#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4']

export default function CodonVisualization({ codonData }) {
  if (!codonData) {
    return (
      <div className="text-center py-16">
        <BeakerIcon className="h-24 w-24 mx-auto text-gray-300 mb-6" />
        <h2 className="text-2xl font-bold text-gray-700 mb-2">
          Sin datos de codones
        </h2>
        <p className="text-gray-500">
          Ejecute el análisis completo para ver la visualización de codones
        </p>
      </div>
    )
  }

  // Prepare data for bar chart (ATG vs Stop codons)
  const codonComparisonData = [
    { name: 'ATG (Inicio)', count: codonData.atg_count, fill: '#10b981' },
    { name: 'TAA', count: codonData.stop_codons.TAA.count, fill: '#f59e0b' },
    { name: 'TAG', count: codonData.stop_codons.TAG.count, fill: '#ef4444' },
    { name: 'TGA', count: codonData.stop_codons.TGA.count, fill: '#8b5cf6' },
  ]

  // Prepare data for pie chart (stop codon distribution)
  const stopCodonPieData = Object.entries(codonData.stop_codons).map(([codon, data]) => ({
    name: codon,
    value: data.count,
    percentage: data.percentage
  }))

  const totalStopCodons = Object.values(codonData.stop_codons)
    .reduce((sum, data) => sum + data.count, 0)

  return (
    <div className="space-y-6">
      {/* Header Stats */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-6 text-white">
          <p className="text-emerald-100">Longitud del Genoma</p>
          <p className="text-3xl font-bold mt-2">
            {(codonData.genome_length / 1000000).toFixed(2)}M
          </p>
          <p className="text-emerald-100 text-sm">pares de bases</p>
        </div>
        <div className="bg-gradient-to-br from-blue-500 to-blue-600 rounded-xl p-6 text-white">
          <p className="text-blue-100">Codones ATG</p>
          <p className="text-3xl font-bold mt-2">
            {codonData.atg_count.toLocaleString()}
          </p>
          <p className="text-blue-100 text-sm">codones de inicio</p>
        </div>
        <div className="bg-gradient-to-br from-purple-500 to-purple-600 rounded-xl p-6 text-white">
          <p className="text-purple-100">Densidad ATG</p>
          <p className="text-3xl font-bold mt-2">
            {codonData.atg_density}
          </p>
          <p className="text-purple-100 text-sm">por kilobase</p>
        </div>
        <div className="bg-gradient-to-br from-orange-500 to-orange-600 rounded-xl p-6 text-white">
          <p className="text-orange-100">Total Stop Codons</p>
          <p className="text-3xl font-bold mt-2">
            {totalStopCodons.toLocaleString()}
          </p>
          <p className="text-orange-100 text-sm">TAA + TAG + TGA</p>
        </div>
      </div>

      {/* Charts Row */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Bar Chart - Codon Comparison */}
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">
            Comparación de Codones
          </h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={codonComparisonData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
              <XAxis dataKey="name" tick={{ fontSize: 12 }} />
              <YAxis 
                tick={{ fontSize: 12 }}
                tickFormatter={(value) => value >= 1000 ? `${(value/1000).toFixed(0)}K` : value}
              />
              <Tooltip 
                formatter={(value) => [value.toLocaleString(), 'Conteo']}
                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgba(0,0,0,0.1)' }}
              />
              <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                {codonComparisonData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={entry.fill} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Pie Chart - Stop Codon Distribution */}
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="font-semibold text-gray-800 mb-4">
            Distribución de Codones de Parada
          </h3>
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={stopCodonPieData}
                cx="50%"
                cy="50%"
                innerRadius={60}
                outerRadius={100}
                paddingAngle={5}
                dataKey="value"
                label={({ name, percentage }) => `${name}: ${percentage}%`}
              >
                {stopCodonPieData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index + 1]} />
                ))}
              </Pie>
              <Tooltip 
                formatter={(value) => [value.toLocaleString(), 'Conteo']}
                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgba(0,0,0,0.1)' }}
              />
              <Legend />
            </PieChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Detailed Stop Codon Table */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <h3 className="font-semibold text-gray-800 mb-4">
          Análisis Detallado de Codones de Parada
        </h3>
        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="text-left text-sm text-gray-500 border-b">
                <th className="pb-3 font-medium">Codón</th>
                <th className="pb-3 font-medium">Secuencia</th>
                <th className="pb-3 font-medium">Conteo</th>
                <th className="pb-3 font-medium">Porcentaje</th>
                <th className="pb-3 font-medium">Frecuencia (por kb)</th>
                <th className="pb-3 font-medium">Distribución</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(codonData.stop_codons).map(([codon, data], index) => (
                <tr key={codon} className="border-b last:border-b-0">
                  <td className="py-4">
                    <span className={`inline-flex px-3 py-1 rounded-full text-sm font-medium ${
                      codon === 'TAA' ? 'bg-yellow-100 text-yellow-800' :
                      codon === 'TAG' ? 'bg-red-100 text-red-800' :
                      'bg-purple-100 text-purple-800'
                    }`}>
                      {codon}
                    </span>
                  </td>
                  <td className="py-4 font-mono text-gray-600">
                    {codon === 'TAA' && 'Ocre'}
                    {codon === 'TAG' && 'Ámbar'}
                    {codon === 'TGA' && 'Ópalo'}
                  </td>
                  <td className="py-4 font-semibold text-gray-800">
                    {data.count.toLocaleString()}
                  </td>
                  <td className="py-4 text-gray-600">
                    {data.percentage}%
                  </td>
                  <td className="py-4 text-gray-600">
                    {((data.count / codonData.genome_length) * 1000).toFixed(3)}
                  </td>
                  <td className="py-4 w-48">
                    <div className="bg-gray-200 rounded-full h-3 overflow-hidden">
                      <div 
                        className={`h-full transition-all duration-500 ${
                          codon === 'TAA' ? 'bg-yellow-500' :
                          codon === 'TAG' ? 'bg-red-500' :
                          'bg-purple-500'
                        }`}
                        style={{ width: `${data.percentage}%` }}
                      ></div>
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Gene Comparison */}
      <div className="bg-gradient-to-r from-cyan-500 to-blue-500 rounded-xl p-6 text-white">
        <h3 className="font-semibold mb-4">Comparación: ATG vs Genes Anotados</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className="text-center">
            <p className="text-cyan-100">Codones ATG Encontrados</p>
            <p className="text-4xl font-bold mt-2">
              {codonData.gene_comparison.atg_found.toLocaleString()}
            </p>
          </div>
          <div className="text-center">
            <p className="text-cyan-100">Genes Anotados</p>
            <p className="text-4xl font-bold mt-2">
              {codonData.gene_comparison.annotated_genes.toLocaleString()}
            </p>
          </div>
          <div className="text-center">
            <p className="text-cyan-100">Diferencia</p>
            <p className="text-4xl font-bold mt-2">
              {codonData.gene_comparison.difference > 0 ? '+' : ''}
              {codonData.gene_comparison.difference.toLocaleString()}
            </p>
          </div>
        </div>
        <p className="mt-4 text-cyan-100 text-sm">
          Nota: La diferencia se debe a que no todos los ATG son codones de inicio de genes. 
          También existen ATG internos dentro de secuencias codificantes.
        </p>
      </div>

      {/* Info Box */}
      <div className="bg-amber-50 border border-amber-200 rounded-xl p-4">
        <div className="flex">
          <div className="flex-shrink-0">
            <svg className="h-5 w-5 text-amber-400" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="ml-3">
            <h4 className="text-sm font-medium text-amber-800">Sobre los codones de parada</h4>
            <p className="text-sm text-amber-700 mt-1">
              <strong>TAA (Ocre)</strong> es el codón de parada más común en E. coli (~63%). 
              <strong> TAG (Ámbar)</strong> es el menos frecuente (~7%), mientras que 
              <strong> TGA (Ópalo)</strong> representa ~30% de los codones de terminación.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
