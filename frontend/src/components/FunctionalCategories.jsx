/**
 * FunctionalCategories Component — Clean Laboratory Edition
 * COG Functional classification with interactive charts and detailed table.
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
  PieChart,
  Pie,
  Legend
} from 'recharts'

const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900/95 backdrop-blur-xl border border-white/10 p-4 rounded-2xl shadow-2xl animate-in fade-in zoom-in-95">
        <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">{label}</p>
        <p className="text-xs font-bold text-white uppercase tracking-tight">
          Genes: <span className="text-blue-200">{payload[0].value}</span>
        </p>
      </div>
    )
  }
  return null
}

export default function FunctionalCategories() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState('visual')

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      const result = await api.getFunctionalCategories()
      setData(result)
    } catch (e) {
      console.error("Error cargando categorías COG:", e)
    } finally {
      setLoading(false)
    }
  }

  const chartData = useMemo(() => {
    if (!data?.categories) return []
    return data.categories.map(c => ({
      code: c.code,
      name: c.name,
      count: c.count,
      color: c.color
    }))
  }, [data])

  const pieData = useMemo(() => {
    if (!data) return []
    return [
      { name: 'Clasificados', value: data.categorized, fill: '#2563eb' },
      { name: 'Sin Clasificar', value: data.uncategorized, fill: '#cbd5e1' }
    ]
  }, [data])

  if (loading && !data) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse text-center">Catalogando Funciones...</p>
      </div>
    )
  }

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Header Section */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            Categorías <span className="text-blue-600">Funcionales (COG)</span>
          </h2>
          <div className="flex flex-wrap items-center gap-4">
            <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-blue-100 shadow-sm">
              {data.organism || 'Genoma en Análisis'}
            </span>
            <p className="text-[10px] font-bold text-slate-600 uppercase tracking-widest italic">
              {data.categorized} de {data.total_cds} CDS clasificados ({((data.categorized/data.total_cds)*100).toFixed(1)}%)
            </p>
          </div>
        </div>
        
        {/* Tab Navigation */}
        <div className="flex bg-slate-50 p-1.5 rounded-2xl border border-slate-100 relative z-10">
          {[
            { id: 'visual', label: 'Gráficos', icon: '📊' },
            { id: 'table', label: 'Tabla Detallada', icon: '📋' }
          ].map(tab => (
            <button
              key={tab.id}
              onClick={() => setActiveTab(tab.id)}
              className={`flex items-center gap-3 px-6 py-2.5 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all ${
                activeTab === tab.id ? 'bg-white text-blue-600 shadow-sm border border-slate-100' : 'text-slate-400 hover:text-slate-600'
              }`}
            >
              <span>{tab.icon}</span>
              <span>{tab.label}</span>
            </button>
          ))}
        </div>
      </div>

      {/* Summary Header Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        {[
          { label: 'Total CDS', val: data.total_cds, icon: 'bg-slate-900', text: 'text-white' },
          { label: 'Clasificados', val: data.categorized, icon: 'bg-blue-600', text: 'text-blue-600' },
          { label: 'Sin Clasificar', val: data.uncategorized, icon: 'bg-slate-400', text: 'text-slate-400' },
          { label: 'Categorías', val: data.categories?.length, icon: 'bg-indigo-600', text: 'text-indigo-600' }
        ].map((item, i) => (
          <div key={i} className={`${item.icon === 'bg-slate-900' ? 'bg-slate-900 text-white shadow-blue-900/20 shadow-2xl' : 'bg-white border-2 border-slate-100 shadow-sm'} rounded-[2rem] p-8 group hover:border-blue-200 transition-all relative overflow-hidden`}>
            {item.icon === 'bg-slate-900' && <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000"></div>}
            <p className={`text-[9px] font-black uppercase tracking-widest mb-3 ${item.icon === 'bg-slate-900' ? 'text-blue-400' : 'text-slate-400'}`}>{item.label}</p>
            <p className="text-3xl font-black tracking-tighter">{item.val?.toLocaleString()}</p>
          </div>
        ))}
      </div>

      {activeTab === 'visual' && (
        <div className="space-y-8 animate-in slide-in-from-bottom-4 duration-700">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
            {/* Pie Chart: Coverage */}
            <div className="lg:col-span-1 bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm flex flex-col items-center">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Cobertura de Clasificación</h3>
              <ResponsiveContainer width="100%" height={250}>
                <PieChart>
                  <Pie data={pieData} innerRadius={60} outerRadius={80} paddingAngle={5} dataKey="value">
                    {pieData.map((entry, index) => <Cell key={`cell-${index}`} fill={entry.fill} cornerRadius={10} />)}
                  </Pie>
                  <Tooltip contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)' }} />
                  <Legend verticalAlign="bottom" iconType="circle" />
                </PieChart>
              </ResponsiveContainer>
              <div className="mt-8 pt-8 border-t border-slate-50 w-full text-center">
                <p className="text-4xl font-black text-blue-600 tracking-tighter">{((data.categorized/data.total_cds)*100).toFixed(1)}%</p>
                <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest mt-2">Éxito de Catalogación</p>
              </div>
            </div>

            {/* Bar Chart: Distribution */}
            <div className="lg:col-span-2 bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10 text-center">Genes por Categoría Principal</h3>
              <ResponsiveContainer width="100%" height={350}>
                <BarChart data={chartData.slice(0, 12)}>
                  <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
                  <XAxis dataKey="code" axisLine={false} tickLine={false} tick={{ fontSize: 10, fontWeight: 900, fill: '#1e293b' }} />
                  <YAxis axisLine={false} tickLine={false} tick={{ fontSize: 9, fill: '#94a3b8' }} />
                  <Tooltip content={<CustomTooltip />} cursor={{ fill: '#f8fafc' }} />
                  <Bar dataKey="count" radius={[6, 6, 0, 0]} barSize={35}>
                    {chartData.map((entry, index) => <Cell key={`cell-${index}`} fill={entry.color} />)}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Detailed Distribution Cards */}
          <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm">
            <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em] mb-10">Distribución de Roles Biológicos</h3>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {chartData.map((cat) => (
                <div key={cat.code} className="p-6 bg-slate-50 rounded-[2rem] border border-slate-100 group hover:border-blue-200 hover:bg-white transition-all shadow-sm">
                  <div className="flex items-start justify-between mb-4">
                    <div className="flex items-center gap-4">
                      <div className="w-10 h-10 rounded-2xl flex items-center justify-center text-xs font-black text-white shadow-lg shadow-blue-900/10" style={{ backgroundColor: cat.color }}>
                        {cat.code}
                      </div>
                      <div className="space-y-1">
                        <h5 className="text-[10px] font-black text-slate-900 uppercase tracking-tight leading-tight max-w-[180px]">{cat.name}</h5>
                        <p className="text-[9px] font-bold text-slate-400">Categoría COG</p>
                      </div>
                    </div>
                    <span className="text-sm font-black text-slate-900">{cat.count}</span>
                  </div>
                  <div className="w-full h-1 bg-white rounded-full overflow-hidden">
                    <div className="h-full transition-all duration-1000 shadow-sm" style={{ backgroundColor: cat.color, width: `${(cat.count / data.total_cds) * 100 * 5}%` }}></div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      )}

      {activeTab === 'table' && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm animate-in zoom-in-95 duration-500">
          <div className="p-8 border-b border-slate-50 bg-slate-50/30 flex justify-between items-center">
            <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Diccionario Funcional Detallado</h3>
            <div className="flex gap-4">
              <span className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Dataset: {data.organism}</span>
            </div>
          </div>
          <div className="overflow-x-auto max-h-[700px] overflow-y-auto custom-scrollbar">
            <table className="w-full text-left">
              <thead className="bg-white sticky top-0 z-10 border-b border-slate-100">
                <tr>
                  <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Código</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest">Categoría Funcional</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Cantidad</th>
                  <th className="px-6 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Porcentaje</th>
                  <th className="px-10 py-6 text-[9px] font-black text-slate-500 uppercase tracking-widest text-right">Distribución</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-50">
                {data.categories?.map((cat, i) => (
                  <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                    <td className="px-10 py-5">
                      <div className="w-8 h-8 rounded-xl flex items-center justify-center text-white font-black text-[10px] shadow-sm" style={{ backgroundColor: cat.color }}>
                        {cat.code}
                      </div>
                    </td>
                    <td className="px-6 py-5">
                      <p className="text-[10px] font-black text-slate-900 uppercase tracking-tight">{cat.name}</p>
                    </td>
                    <td className="px-6 py-5 text-right font-mono text-xs font-black text-slate-700">{cat.count.toLocaleString()}</td>
                    <td className="px-6 py-5 text-right font-mono text-xs font-black text-blue-600">
                      {((cat.count / data.total_cds) * 100).toFixed(1)}%
                    </td>
                    <td className="px-10 py-5">
                      <div className="flex items-center justify-end gap-4">
                        <div className="w-24 h-1 bg-slate-100 rounded-full overflow-hidden">
                          <div className="h-full transition-all duration-1000 shadow-sm" style={{ backgroundColor: cat.color, width: `${(cat.count / data.total_cds) * 100 * 5}%` }}></div>
                        </div>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Glossary & Technical Definition */}
      <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 space-y-8 overflow-hidden relative">
        <div className="absolute bottom-0 right-0 p-12 opacity-5 pointer-events-none text-9xl font-black italic">COG</div>
        <div className="space-y-4 relative z-10">
          <h4 className="text-xl font-black uppercase tracking-widest text-blue-400">Sistema de Clasificación COG</h4>
          <p className="text-sm text-slate-300 leading-relaxed font-medium max-w-4xl">
            El sistema <span className="text-white font-black">COG (Clusters of Orthologous Groups)</span> es un modelo desarrollado por el NCBI para la anotación funcional de proteínas basada en relaciones evolutivas. Agrupa secuencias ortólogas que comparten funciones biológicas fundamentales a través de distintos linajes genómicos.
          </p>
        </div>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-10 relative z-10 pt-6 border-t border-white/5">
          <p className="text-xs text-slate-400 leading-relaxed font-medium">
            <span className="text-blue-400 font-bold block mb-2">Metodología de Análisis:</span> La clasificación se realiza mediante el procesamiento de palabras clave en las descripciones de productos génicos, mapeándolas contra descriptores estandarizados de roles metabólicos, regulatorios y estructurales.
          </p>
          <p className="text-xs text-slate-400 leading-relaxed font-medium">
            <span className="text-indigo-400 font-bold block mb-2">Categorías Críticas:</span> El motor resalta funciones esenciales como <span className="text-white">Traducción (J)</span>, <span className="text-white">Transcripción (K)</span> y <span className="text-white">Metabolismo Energético (C)</span>, fundamentales para la viabilidad celular del MG1655.
          </p>
        </div>
      </div>
    </div>
  )
}
