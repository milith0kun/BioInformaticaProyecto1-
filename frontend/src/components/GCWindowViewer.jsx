/**
 * GCWindowViewer Component — Clean Laboratory Edition
 * Full GC sliding window analysis with Content, Skew, and combined views
 */
import { useState, useEffect, useMemo } from 'react'
import api from '../services/api'
import {
  AreaChart,
  Area,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  LineChart,
  Line,
  ReferenceLine
} from 'recharts'

export default function GCWindowViewer() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [windowSize, setWindowSize] = useState(5000)
  const [step, setStep] = useState(1000)
  const [viewMode, setViewMode] = useState('both') // 'gc', 'skew', 'both'

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      const result = await api.getGCSlidingWindow(windowSize, step)
      setData(result)
    } catch (e) {
      console.error(e)
    } finally {
      setLoading(false)
    }
  }

  // Compute statistics from data
  const stats = useMemo(() => {
    if (!data?.gc_windows?.length) return null

    const backendStats = data.gc_stats || {}
    return {
      average: (data.average_gc || 0).toFixed(2),
      min: (backendStats.min_gc || 0).toFixed(2),
      max: (backendStats.max_gc || 0).toFixed(2),
      stdDev: (backendStats.std_gc || 0).toFixed(2),
      genomeSize: data.genome_length || 0,
      windowCount: data.total_windows || data.gc_windows.length,
      organism: data.organism || 'Genoma activo',
      accession: data.accession || ''
    }
  }, [data])

  // Format position for X axis
  const formatPosition = (val) => {
    if (val == null) return ''
    if (val >= 1e6) return `${(val / 1e6).toFixed(2)}M`
    if (val >= 1e3) return `${(val / 1e3).toFixed(0)}K`
    return val.toString()
  }

  if (loading && !data) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse">Escaneando Secuencia...</p>
      </div>
    )
  }

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Header + Controls */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="flex flex-col gap-8 relative z-10">
          {/* Title */}
          <div className="space-y-2">
            <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
              📊 Análisis GC — <span className="text-blue-600">Ventana Deslizante</span>
            </h2>
            <p className="text-[10px] font-bold text-slate-400 uppercase tracking-[0.2em]">
              {stats?.organism} {stats?.accession ? `(${stats.accession})` : ''} — {stats?.windowCount?.toLocaleString() || 0} ventanas
            </p>
          </div>

          {/* Controls Row */}
          <div className="flex flex-wrap items-end gap-6">
            <div className="space-y-2">
              <label className="text-[9px] font-black text-slate-400 uppercase tracking-widest px-2">Ventana:</label>
              <select
                value={windowSize}
                onChange={(e) => setWindowSize(Number(e.target.value))}
                className="block w-28 px-4 py-2.5 text-xs font-bold bg-slate-50 border-2 border-slate-100 rounded-xl focus:outline-none focus:border-blue-500 transition-all"
              >
                <option value={1000}>1 kb</option>
                <option value={2000}>2 kb</option>
                <option value={5000}>5 kb</option>
                <option value={10000}>10 kb</option>
                <option value={20000}>20 kb</option>
              </select>
            </div>
            <div className="space-y-2">
              <label className="text-[9px] font-black text-slate-400 uppercase tracking-widest px-2">Paso:</label>
              <select
                value={step}
                onChange={(e) => setStep(Number(e.target.value))}
                className="block w-28 px-4 py-2.5 text-xs font-bold bg-slate-50 border-2 border-slate-100 rounded-xl focus:outline-none focus:border-blue-500 transition-all"
              >
                <option value={500}>0.5 kb</option>
                <option value={1000}>1 kb</option>
                <option value={2000}>2 kb</option>
                <option value={5000}>5 kb</option>
              </select>
            </div>
            <button
              onClick={loadData}
              disabled={loading}
              className="px-8 py-2.5 bg-blue-600 text-white rounded-xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-700 transition-all shadow-lg shadow-blue-200 disabled:opacity-50"
            >
              {loading ? 'Calculando...' : 'Recalcular'}
            </button>
          </div>

          {/* View Mode Tabs */}
          <div className="flex bg-slate-50 p-1.5 rounded-2xl border border-slate-100 w-fit">
            {[
              { id: 'gc', label: '📈 GC Content' },
              { id: 'skew', label: '📉 GC Skew' },
              { id: 'both', label: '📊 Ambos' }
            ].map(tab => (
              <button
                key={tab.id}
                onClick={() => setViewMode(tab.id)}
                className={`px-6 py-2.5 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all ${viewMode === tab.id
                  ? 'bg-white text-blue-600 shadow-sm border border-slate-100'
                  : 'text-slate-400 hover:text-slate-600'
                  }`}
              >
                {tab.label}
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Statistics Cards */}
      {stats && (
        <div className="grid grid-cols-2 md:grid-cols-5 gap-4">
          {[
            { label: 'Promedio GC', value: `${stats.average}%`, color: 'text-blue-600' },
            { label: 'GC Mínimo', value: `${stats.min}%`, color: 'text-emerald-600' },
            { label: 'GC Máximo', value: `${stats.max}%`, color: 'text-rose-600' },
            { label: 'Desv. Estándar', value: `${stats.stdDev}%`, color: 'text-amber-600' },
            { label: 'Genoma', value: `${(stats.genomeSize / 1e6).toFixed(2)} Mb`, color: 'text-indigo-600' }
          ].map((stat, i) => (
            <div key={i} className="bg-white rounded-[2rem] border-2 border-slate-100 p-6 shadow-sm text-center">
              <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest mb-2">{stat.label}</p>
              <p className={`text-2xl font-black ${stat.color}`}>{stat.value}</p>
            </div>
          ))}
        </div>
      )}

      {/* GC Content Chart */}
      {(viewMode === 'gc' || viewMode === 'both') && (
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <div className="mb-6 space-y-2">
            <h3 className="text-sm font-black text-slate-900 uppercase tracking-tight">%GC a lo largo del genoma</h3>
            <p className="text-[10px] font-medium text-slate-400 leading-relaxed">
              Contenido GC calculado con ventana deslizante de {windowSize >= 1000 ? `${windowSize / 1000} kb` : `${windowSize} pb`} y paso de {step >= 1000 ? `${(step / 1000).toFixed(1)} kb` : `${step} pb`}. La línea punteada indica el GC promedio ({stats?.average || 0}%).
            </p>
          </div>
          <ResponsiveContainer width="100%" height={350}>
            <AreaChart data={data.gc_windows || []}>
              <defs>
                <linearGradient id="colorGc" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#2563eb" stopOpacity={0.15} />
                  <stop offset="95%" stopColor="#2563eb" stopOpacity={0} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
              <XAxis
                dataKey="position"
                tickFormatter={formatPosition}
                tick={{ fontSize: 9, fill: '#94a3b8', fontWeight: 700 }}
                axisLine={{ stroke: '#e2e8f0' }}
                tickLine={false}
                interval={Math.max(1, Math.floor((data.gc_windows?.length || 1) / 20))}
                label={{ value: 'Posición genómica', position: 'insideBottom', offset: -5, style: { fontSize: 10, fill: '#94a3b8', fontWeight: 700, textTransform: 'uppercase' } }}
              />
              <YAxis
                domain={['dataMin - 5', 'dataMax + 5']}
                axisLine={false}
                tickLine={false}
                tick={{ fontSize: 10, fill: '#94a3b8', fontWeight: 700 }}
                label={{ value: 'GC%', angle: -90, position: 'insideLeft', style: { fontSize: 10, fill: '#94a3b8', fontWeight: 700 } }}
              />
              <Tooltip
                contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)', fontSize: 11 }}
                labelFormatter={(val) => `Posición: ${(val || 0).toLocaleString()} pb`}
                formatter={(value) => [`${value?.toFixed(2)}%`, 'GC']}
              />
              <Area type="monotone" dataKey="gc" stroke="#2563eb" strokeWidth={2} fillOpacity={1} fill="url(#colorGc)" dot={false} />
              <ReferenceLine
                y={parseFloat(stats?.average || 0)}
                stroke="#94a3b8"
                strokeDasharray="5 5"
                strokeWidth={2}
                label={{ value: `Promedio: ${stats?.average || 0}%`, position: 'right', style: { fontSize: 10, fill: '#64748b', fontWeight: 700 } }}
              />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      )}

      {/* GC Skew Chart */}
      {(viewMode === 'skew' || viewMode === 'both') && (
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
          <div className="mb-6 space-y-2">
            <h3 className="text-sm font-black text-slate-900 uppercase tracking-tight">GC Skew (G−C)/(G+C)</h3>
            <p className="text-[10px] font-medium text-slate-400 leading-relaxed">
              GC Skew indica la asimetría en distribución de G y C entre las dos hebras del DNA. Valores positivos = exceso G (hebra líder), valores negativos = exceso C (hebra rezagada). El punto de cruce indica el origen/terminación de replicación.
            </p>
          </div>
          <ResponsiveContainer width="100%" height={300}>
            <AreaChart data={data.gc_skew || []}>
              <defs>
                <linearGradient id="colorSkewPos" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#4f46e5" stopOpacity={0.15} />
                  <stop offset="95%" stopColor="#4f46e5" stopOpacity={0} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
              <XAxis
                dataKey="position"
                tickFormatter={formatPosition}
                tick={{ fontSize: 9, fill: '#94a3b8', fontWeight: 700 }}
                axisLine={{ stroke: '#e2e8f0' }}
                tickLine={false}
                interval={Math.max(1, Math.floor((data.gc_skew?.length || 1) / 20))}
              />
              <YAxis
                axisLine={false}
                tickLine={false}
                tick={{ fontSize: 10, fill: '#94a3b8', fontWeight: 700 }}
                label={{ value: 'GC Skew', angle: -90, position: 'insideLeft', style: { fontSize: 10, fill: '#94a3b8', fontWeight: 700 } }}
              />
              <Tooltip
                contentStyle={{ borderRadius: '1.5rem', border: 'none', boxShadow: '0 20px 25px -5px rgb(0 0 0 / 0.1)', fontSize: 11 }}
                labelFormatter={(val) => `Posición: ${(val || 0).toLocaleString()} pb`}
                formatter={(value) => [value?.toFixed(4), 'Skew']}
              />
              <ReferenceLine y={0} stroke="#cbd5e1" strokeWidth={2} />
              <Area type="monotone" dataKey="skew" stroke="#4f46e5" strokeWidth={2} fillOpacity={0.1} fill="url(#colorSkewPos)" dot={false} />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      )}

      {/* Educational Footer */}
      <div className="bg-slate-50 rounded-[2.5rem] border border-slate-100 p-10">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
          <div className="space-y-3">
            <h4 className="text-[10px] font-black text-slate-600 uppercase tracking-widest">Contenido GC (%GC):</h4>
            <p className="text-xs text-slate-500 font-medium leading-relaxed">
              Proporción de bases Guanina + Citosina. Varía a lo largo del genoma según islas CpG, genes, y regiones codificantes vs. no codificantes.
            </p>
          </div>
          <div className="space-y-3">
            <h4 className="text-[10px] font-black text-slate-600 uppercase tracking-widest">GC Skew:</h4>
            <p className="text-xs text-slate-500 font-medium leading-relaxed">
              (G−C)/(G+C) calculado por ventana. Los cambios de signo del GC skew acumulativo indican el origen (oriC) y término de replicación (ter) en genomas circulares.
            </p>
          </div>
        </div>
        <div className="mt-6 pt-6 border-t border-slate-200 flex items-center justify-between">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">
            Ventana actual: {windowSize >= 1000 ? `${windowSize / 1000} kb` : `${windowSize} pb`}, Paso: {step >= 1000 ? `${(step / 1000).toFixed(1)} kb` : `${step} pb`}, Puntos: {stats?.windowCount?.toLocaleString() || 0}
          </p>
          <span className="text-[9px] font-black text-blue-500 uppercase tracking-widest animate-pulse">LABORATORY SYSTEM ONLINE</span>
        </div>
      </div>
    </div>
  )
}