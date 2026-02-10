/**
 * CodonVisualization Component — Clean Laboratory Edition
 * Visualizes codon distribution, bias, and technical metrics with interactive charts.
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
  Legend,
  AreaChart,
  Area
} from 'recharts'

const COLORS = ['#2563eb', '#4f46e5', '#7c3aed', '#db2777']

function formatLargeNumber(num) {
  if (!num) return '0'
  if (num >= 1000000) return (num / 1000000).toFixed(2) + 'M'
  if (num >= 1000) return (num / 1000).toLocaleString('es-ES', { minimumFractionDigits: 3 })
  return num.toLocaleString('es-ES')
}

const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900/95 backdrop-blur-xl border border-white/10 p-4 rounded-2xl shadow-2xl animate-in fade-in zoom-in-95">
        <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">{label}</p>
        <p className="text-xs font-bold text-white uppercase tracking-tight">
          Conteo: <span className="text-blue-200">{payload[0].value.toLocaleString()}</span>
        </p>
      </div>
    )
  }
  return null
}

export default function CodonVisualization({ codonData }) {
  if (!codonData) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-600 uppercase tracking-widest">Sincronizando Dataset de Codones...</p>
      </div>
    )
  }

  const codonComparisonData = [
    { name: 'ATG (Inicio)', count: codonData.atg_count || 0, fill: '#2563eb' },
    { name: 'TAA (Ocre)', count: codonData.stop_codons?.TAA?.count || 0, fill: '#4f46e5' },
    { name: 'TAG (Ámbar)', count: codonData.stop_codons?.TAG?.count || 0, fill: '#7c3aed' },
    { name: 'TGA (Ópalo)', count: codonData.stop_codons?.TGA?.count || 0, fill: '#db2777' },
  ]

  const totalStopCodons = Object.values(codonData.stop_codons || {})
    .reduce((sum, d) => sum + (d.count || 0), 0)

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Principal Metrics Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-slate-600 uppercase tracking-widest mb-3 leading-none">Extensión de Hebra</p>
          <p className="text-3xl font-black text-slate-900 tracking-tighter">{formatLargeNumber(codonData.genome_length)}</p>
          <p className="text-[9px] font-bold text-slate-500 uppercase tracking-tight mt-1">pares de bases</p>
        </div>
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-blue-600 uppercase tracking-widest mb-3 leading-none">Iniciadores ATG</p>
          <p className="text-3xl font-black text-blue-700 tracking-tighter">{formatLargeNumber(codonData.atg_count)}</p>
          <p className="text-[9px] font-bold text-slate-500 uppercase tracking-tight mt-1">sitios de iniciación</p>
        </div>
        <div className="bg-white rounded-3xl border-2 border-slate-100 p-8 shadow-sm group hover:border-blue-200 transition-all">
          <p className="text-[10px] font-black text-indigo-600 uppercase tracking-widest mb-3 leading-none">Densidad ATG</p>
          <p className="text-3xl font-black text-indigo-700 tracking-tighter">{codonData.atg_density || 0}</p>
          <p className="text-[9px] font-bold text-slate-500 uppercase tracking-tight mt-1">codones por kb</p>
        </div>
        <div className="bg-slate-900 rounded-3xl p-8 text-white shadow-2xl shadow-blue-900/20 relative overflow-hidden group">
          <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000"></div>
          <p className="text-[10px] font-black text-blue-300 uppercase tracking-widest mb-3 leading-none">Codones de Parada</p>
          <p className="text-3xl font-black tracking-tighter text-white">{formatLargeNumber(totalStopCodons)}</p>
          <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight mt-1">Total TAA/TAG/TGA</p>
        </div>
      </div>

      {/* Analytics Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Frecuencia de Codones Críticos</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={codonComparisonData}>
              <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
              <XAxis dataKey="name" axisLine={false} tickLine={false} tick={{ fontSize: 9, fontWeight: 900, fill: '#475569' }} dy={10} />
              <YAxis axisLine={false} tickLine={false} tick={{ fontSize: 9, fill: '#94a3b8' }} />
              <Tooltip content={<CustomTooltip />} cursor={{ fill: '#f8fafc' }} />
              <Bar dataKey="count" radius={[10, 10, 0, 0]} barSize={40}>
                {codonComparisonData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={entry.fill} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm flex flex-col items-center">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Distribución de Terminación</h3>
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={Object.entries(codonData.stop_codons || {}).map(([name, d]) => ({ name, value: d.count || 0 }))}
                cx="50%" cy="50%" innerRadius={60} outerRadius={90} paddingAngle={8} dataKey="value"
              >
                {Object.keys(codonData.stop_codons || {}).map((_, i) => (
                  <Cell key={`cell-${i}`} fill={COLORS[i % COLORS.length]} cornerRadius={10} />
                ))}
              </Pie>
              <Tooltip contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)', fontSize: 11 }} />
              <Legend 
                verticalAlign="bottom" 
                iconType="circle"
                formatter={(value) => {
                  const d = codonData.stop_codons[value]
                  return <span className="text-[10px] font-black text-slate-700 uppercase ml-2">{value}: {d?.percentage?.toFixed(1)}%</span>
                }}
              />
            </PieChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Detailed Analysis Table */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
        <div className="p-8 border-b border-slate-50 bg-slate-50/30">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Métricas de Terminación de Lectura</h3>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left">
            <thead>
              <tr className="bg-white border-b border-slate-100">
                <th className="px-10 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest">Codón</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest">Identidad</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest text-right">Conteo</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest text-right">Porcentaje</th>
                <th className="px-6 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest text-right">Densidad</th>
                <th className="px-10 py-6 text-[10px] font-black text-slate-500 uppercase tracking-widest text-center">Dispersión Espacial</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-50">
              {Object.entries(codonData.stop_codons || {}).map(([codon, d]) => (
                <tr key={codon} className="hover:bg-slate-50/50 transition-colors group">
                  <td className="px-10 py-6 font-mono text-sm font-black text-blue-600 tracking-[0.2em]">{codon}</td>
                  <td className="px-6 py-6">
                    <span className="text-[10px] font-black text-slate-800 uppercase tracking-widest px-3 py-1 bg-slate-100 rounded-lg">
                      {codon === 'TAA' ? 'Ocre' : codon === 'TAG' ? 'Ámbar' : 'Ópalo'}
                    </span>
                  </td>
                  <td className="px-6 py-6 text-right text-xs font-black text-slate-900">{d.count?.toLocaleString()}</td>
                  <td className="px-6 py-6 text-right">
                    <div className="flex flex-col items-end gap-1">
                      <span className="text-xs font-black text-blue-700">{(d.percentage || 0).toFixed(2)}%</span>
                      <div className="w-16 h-1 bg-slate-100 rounded-full overflow-hidden">
                        <div className="h-full bg-blue-600" style={{ width: `${d.percentage}%` }}></div>
                      </div>
                    </div>
                  </td>
                  <td className="px-6 py-6 text-right font-mono text-xs font-bold text-slate-600">{(d.density_per_kb || 0).toFixed(3)}/kb</td>
                  <td className="px-10 py-6">
                    <div className="w-28 h-8 mx-auto grayscale group-hover:grayscale-0 transition-all opacity-40 group-hover:opacity-100">
                      <ResponsiveContainer width="100%" height={100}>
                        <AreaChart data={(d.spatial_distribution?.count_per_window || [0,5,2,8,3,6]).map(v => ({v}))}>
                          <Area type="monotone" dataKey="v" stroke={COLORS[Object.keys(codonData.stop_codons).indexOf(codon)]} fill={COLORS[Object.keys(codonData.stop_codons).indexOf(codon)]} fillOpacity={0.2} strokeWidth={1.5} dot={false} />
                        </AreaChart>
                      </ResponsiveContainer>
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Comparison: ATG vs Genes */}
      {codonData.gene_comparison && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm">
          <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">Análisis Comparativo: Iniciadores vs Anotación</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div className="bg-slate-50 p-8 rounded-[2rem] border border-slate-100 hover:border-blue-200 transition-all shadow-sm">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-3 leading-none">Total ATG Detectados</p>
              <p className="text-3xl font-black text-slate-900 tracking-tighter">{codonData.gene_comparison.atg_found?.toLocaleString()}</p>
            </div>
            <div className="bg-slate-50 p-8 rounded-[2rem] border border-slate-100 hover:border-blue-200 transition-all shadow-sm">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-3 leading-none">Genes Anotados (Ref)</p>
              <p className="text-3xl font-black text-blue-700 tracking-tighter">{codonData.gene_comparison.annotated_genes?.toLocaleString()}</p>
            </div>
            <div className="bg-slate-50 p-8 rounded-[2rem] border border-slate-100 hover:border-blue-200 transition-all shadow-sm">
              <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-3 leading-none">Diferencia Neta</p>
              <p className={`text-3xl font-black tracking-tighter ${codonData.gene_comparison.difference > 0 ? 'text-rose-700' : 'text-emerald-700'}`}>
                {codonData.gene_comparison.difference > 0 ? '+' : ''}{codonData.gene_comparison.difference?.toLocaleString()}
              </p>
            </div>
          </div>
          <div className="mt-8 p-8 bg-blue-50/50 border-2 border-blue-100 rounded-[2rem] flex gap-6 items-start">
            <div className="text-2xl text-blue-600 font-black">!</div>
            <div className="space-y-2">
              <h4 className="text-[10px] font-black text-blue-700 uppercase tracking-widest">Interpretación del Sesgo</h4>
              <p className="text-[11px] font-medium text-blue-900/70 leading-relaxed uppercase tracking-tight">
                La disparidad entre ATG detectados y genes anotados se debe a que el triplete ATG también codifica para Metionina interna en las proteínas. No todas las instancias de ATG actúan como codones de inicio funcionales.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Scientific Context Footer */}
      <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 flex flex-col md:flex-row gap-10 items-center overflow-hidden relative">
        <div className="absolute top-0 left-0 w-full h-full opacity-5 pointer-events-none text-9xl font-black italic">STOP</div>
        <div className="w-16 h-16 bg-blue-600 rounded-3xl flex items-center justify-center text-3xl shadow-lg relative z-10 font-black italic">!</div>
        <p className="text-sm font-medium leading-relaxed italic text-slate-300 relative z-10 max-w-4xl">
          "En el genoma de E. coli MG1655, el codón <span className="text-white font-black">TAA (Ocre)</span> es el terminador predominante (~63%), reflejando una fuerte presión selectiva para una terminación eficiente. <span className="text-white font-black">TAG (Ámbar)</span> y <span className="text-white font-black">TGA (Ópalo)</span> se utilizan con menor frecuencia, a menudo en genes con niveles de expresión moderados o regulados."
        </p>
      </div>
    </div>
  )
}