/**
 * ComparisonResults - Muestra resultados de comparación de múltiples genomas
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
  Legend
} from 'recharts'

const COLORS = ['#0d9488', '#10b981', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#14b8a6']

export default function ComparisonResults({ comparisonResult, selectedGenomes }) {
  // Preparar datos para gráficos
  const chartData = useMemo(() => {
    if (!comparisonResult?.genomes) return []
    return comparisonResult.genomes.map((g, i) => ({
      name: g.organism_name?.split(' ').slice(0, 2).join(' ') || g.accession.substring(0, 15),
      accession: g.accession,
      'Tamaño (Mb)': (g.genome_length / 1000000).toFixed(2),
      'Genes': g.gene_count,
      'GC%': g.gc_content?.toFixed(1),
      'Densidad': g.gene_density?.toFixed(1),
      'Gen más largo': g.max_gene_length,
      'Gen más corto': g.min_gene_length,
      color: COLORS[i % COLORS.length]
    }))
  }, [comparisonResult])

  if (!comparisonResult) {
    return (
      <div className="text-center py-12 text-slate-500">
        <p>No hay datos de comparación disponibles</p>
        <p className="text-sm mt-1">Selecciona múltiples genomas y ejecuta el análisis</p>
      </div>
    )
  }

  return (
    <div className="space-y-6">
      {/* Resumen */}
      <div className="bg-gradient-to-r from-teal-600 to-emerald-600 rounded-xl p-6 text-white">
        <div className="flex items-center justify-between mb-4">
          <div>
            <h2 className="text-2xl font-bold">Análisis Comparativo</h2>
            <p className="text-teal-100">
              {comparisonResult.total_genomes_compared} genomas analizados • {comparisonResult.comparison_date}
            </p>
          </div>
          <div className="text-right">
            <div className="text-4xl font-bold">{comparisonResult.total_genomes_compared}</div>
            <div className="text-teal-100 text-sm">genomas</div>
          </div>
        </div>
      </div>

      {/* Métricas extremas */}
      {comparisonResult.extremes && (
        <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-3">
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Más Grande</p>
            <p className="text-sm font-bold text-red-600 mt-1 truncate">{comparisonResult.extremes.largest_genome}</p>
          </div>
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Más Pequeño</p>
            <p className="text-sm font-bold text-green-600 mt-1 truncate">{comparisonResult.extremes.smallest_genome}</p>
          </div>
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Mayor GC%</p>
            <p className="text-sm font-bold text-purple-600 mt-1 truncate">{comparisonResult.extremes.highest_gc}</p>
          </div>
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Menor GC%</p>
            <p className="text-sm font-bold text-orange-600 mt-1 truncate">{comparisonResult.extremes.lowest_gc}</p>
          </div>
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Mayor Densidad</p>
            <p className="text-sm font-bold text-blue-600 mt-1 truncate">{comparisonResult.extremes.highest_gene_density}</p>
          </div>
          <div className="bg-white border border-slate-200 rounded-xl p-4">
            <p className="text-xs text-slate-500 uppercase font-medium">Menor Densidad</p>
            <p className="text-sm font-bold text-slate-600 mt-1 truncate">{comparisonResult.extremes.lowest_gene_density}</p>
          </div>
        </div>
      )}

      {/* Tabla comparativa */}
      <div className="bg-white border border-slate-200 rounded-xl overflow-hidden">
        <div className="p-4 border-b border-slate-200 bg-slate-50">
          <h3 className="font-bold text-slate-800">📊 Comparación Detallada</h3>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full">
            <thead className="bg-slate-50 border-b border-slate-200">
              <tr>
                <th className="px-4 py-3 text-left text-xs font-bold text-slate-600 uppercase">Organismo</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Tamaño</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Genes</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">GC%</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Densidad</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Gen Mayor</th>
                <th className="px-4 py-3 text-right text-xs font-bold text-slate-600 uppercase">Gen Menor</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-100">
              {comparisonResult.genomes?.map((g, i) => (
                <tr key={g.accession} className="hover:bg-slate-50">
                  <td className="px-4 py-3">
                    <div className="flex items-center gap-2">
                      <span 
                        className="w-3 h-3 rounded-full flex-shrink-0" 
                        style={{ backgroundColor: COLORS[i % COLORS.length] }}
                      />
                      <div>
                        <p className="font-medium text-slate-800 text-sm">{g.organism_name || g.accession}</p>
                        <p className="text-xs text-teal-600 font-mono">{g.accession}</p>
                      </div>
                    </div>
                  </td>
                  <td className="px-4 py-3 text-right text-sm text-slate-600">
                    {(g.genome_length / 1000000).toFixed(2)} Mb
                  </td>
                  <td className="px-4 py-3 text-right text-sm text-slate-600">
                    {g.gene_count?.toLocaleString()}
                  </td>
                  <td className="px-4 py-3 text-right text-sm font-medium text-emerald-700">
                    {g.gc_content?.toFixed(1)}%
                  </td>
                  <td className="px-4 py-3 text-right text-sm text-slate-600">
                    {g.gene_density?.toFixed(1)} g/Mb
                  </td>
                  <td className="px-4 py-3 text-right text-sm text-red-600 font-medium">
                    {g.max_gene_length?.toLocaleString()} bp
                  </td>
                  <td className="px-4 py-3 text-right text-sm text-green-600 font-medium">
                    {g.min_gene_length?.toLocaleString()} bp
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Gráficos */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Tamaño del genoma */}
        <div className="bg-white border border-slate-200 rounded-xl p-4">
          <h4 className="font-bold text-slate-800 mb-4">Tamaño del Genoma (Mb)</h4>
          <ResponsiveContainer width="100%" height={250}>
            <BarChart data={chartData} layout="vertical">
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" />
              <YAxis dataKey="name" type="category" width={80} tick={{ fontSize: 10 }} />
              <Tooltip />
              <Bar dataKey="Tamaño (Mb)" fill="#0d9488" radius={[0, 4, 4, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Número de genes */}
        <div className="bg-white border border-slate-200 rounded-xl p-4">
          <h4 className="font-bold text-slate-800 mb-4">Número de Genes</h4>
          <ResponsiveContainer width="100%" height={250}>
            <BarChart data={chartData} layout="vertical">
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" />
              <YAxis dataKey="name" type="category" width={80} tick={{ fontSize: 10 }} />
              <Tooltip />
              <Bar dataKey="Genes" fill="#10b981" radius={[0, 4, 4, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Contenido GC */}
        <div className="bg-white border border-slate-200 rounded-xl p-4">
          <h4 className="font-bold text-slate-800 mb-4">Contenido GC (%)</h4>
          <ResponsiveContainer width="100%" height={250}>
            <BarChart data={chartData}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="name" tick={{ fontSize: 10 }} />
              <YAxis domain={['dataMin - 2', 'dataMax + 2']} />
              <Tooltip />
              <Bar dataKey="GC%" fill="#6366f1" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Densidad génica */}
        <div className="bg-white border border-slate-200 rounded-xl p-4">
          <h4 className="font-bold text-slate-800 mb-4">Densidad Génica (genes/Mb)</h4>
          <ResponsiveContainer width="100%" height={250}>
            <BarChart data={chartData}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="name" tick={{ fontSize: 10 }} />
              <YAxis />
              <Tooltip />
              <Bar dataKey="Densidad" fill="#f59e0b" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Genes extremos globales */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Genes más largos */}
        {comparisonResult.longest_genes_global && comparisonResult.longest_genes_global.length > 0 && (
          <div className="bg-white border border-slate-200 rounded-xl overflow-hidden">
            <div className="p-4 border-b border-slate-200 bg-gradient-to-r from-red-50 to-orange-50">
              <h4 className="font-bold text-slate-800">📏 Genes Más Largos (Global)</h4>
              <p className="text-xs text-slate-500">Top 10 genes más largos de todos los genomas</p>
            </div>
            <div className="divide-y divide-slate-100 max-h-80 overflow-y-auto">
              {comparisonResult.longest_genes_global.map((gene, i) => (
                <div key={i} className="p-3 hover:bg-slate-50">
                  <div className="flex justify-between items-start">
                    <div className="flex-1 min-w-0">
                      <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                      <p className="text-xs text-slate-500 mt-0.5 truncate">{gene.product}</p>
                      <p className="text-xs text-slate-400">{gene.genome}</p>
                    </div>
                    <span className="font-bold text-red-600 ml-2">{gene.length?.toLocaleString()} bp</span>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Genes más cortos */}
        {comparisonResult.shortest_genes_global && comparisonResult.shortest_genes_global.length > 0 && (
          <div className="bg-white border border-slate-200 rounded-xl overflow-hidden">
            <div className="p-4 border-b border-slate-200 bg-gradient-to-r from-green-50 to-teal-50">
              <h4 className="font-bold text-slate-800">🔬 Genes Más Cortos (Global)</h4>
              <p className="text-xs text-slate-500">Top 10 genes más cortos de todos los genomas</p>
            </div>
            <div className="divide-y divide-slate-100 max-h-80 overflow-y-auto">
              {comparisonResult.shortest_genes_global.map((gene, i) => (
                <div key={i} className="p-3 hover:bg-slate-50">
                  <div className="flex justify-between items-start">
                    <div className="flex-1 min-w-0">
                      <span className="font-mono text-sm text-teal-700">{gene.locus_tag}</span>
                      <p className="text-xs text-slate-500 mt-0.5 truncate">{gene.product}</p>
                      <p className="text-xs text-slate-400">{gene.genome}</p>
                    </div>
                    <span className="font-bold text-green-600 ml-2">{gene.length?.toLocaleString()} bp</span>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
