/**
 * GCWindowViewer Component
 * Interactive GC content sliding window visualization with
 * Recharts line graphs for GC% and GC Skew across the genome
 */
import { useState, useEffect } from 'react'
import {
    LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip,
    ResponsiveContainer, ReferenceLine, Area, AreaChart
} from 'recharts'
import api from '../services/api'

export default function GCWindowViewer() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [windowSize, setWindowSize] = useState(5000)
    const [step, setStep] = useState(1000)
    const [viewMode, setViewMode] = useState('gc') // 'gc' | 'skew' | 'both'
    const [zoomRange, setZoomRange] = useState(null)

    useEffect(() => { loadData() }, [])

    const loadData = async (ws, st) => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getGCSlidingWindow(ws || windowSize, st || step)
            setData(result)
        } catch (e) {
            setError(e.response?.data?.detail || 'Error cargando datos GC')
        } finally {
            setLoading(false)
        }
    }

    const handleRecalculate = () => {
        loadData(windowSize, step)
    }

    // Sample data for display (avoid rendering too many points)
    const sampledGC = data?.gc_windows
        ? (data.gc_windows.length > 2000
            ? data.gc_windows.filter((_, i) => i % Math.ceil(data.gc_windows.length / 2000) === 0)
            : data.gc_windows
        ).map(w => ({
            pos: w.position,
            posLabel: `${(w.position / 1e6).toFixed(2)}M`,
            gc: w.gc
        }))
        : []

    const sampledSkew = data?.gc_skew
        ? (data.gc_skew.length > 2000
            ? data.gc_skew.filter((_, i) => i % Math.ceil(data.gc_skew.length / 2000) === 0)
            : data.gc_skew
        ).map(w => ({
            pos: w.position,
            posLabel: `${(w.position / 1e6).toFixed(2)}M`,
            skew: w.skew,
            skewColor: w.skew >= 0 ? '#22c55e' : '#ef4444'
        }))
        : []

    if (error) {
        return (
            <div className="bg-amber-50 border border-amber-200 rounded-xl p-6 text-center">
                <p className="text-amber-700 font-medium">⚠️ {error}</p>
                <p className="text-amber-600 text-sm mt-2">Active un genoma primero para analizar el contenido GC.</p>
            </div>
        )
    }

    return (
        <div className="space-y-6">
            <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
                <div>
                    <h2 className="text-xl font-bold text-slate-800">📊 Análisis GC — Ventana Deslizante</h2>
                    <p className="text-sm text-slate-500">
                        {data ? `${data.organism} (${data.accession}) — ${data.total_windows.toLocaleString()} ventanas` : 'Cargando...'}
                    </p>
                </div>

                <div className="flex items-center gap-3 flex-wrap">
                    <div className="flex items-center gap-1.5">
                        <label className="text-xs text-slate-500">Ventana:</label>
                        <select value={windowSize} onChange={(e) => setWindowSize(parseInt(e.target.value))}
                            className="px-2 py-1.5 border border-slate-200 rounded-lg text-xs">
                            <option value={1000}>1 kb</option>
                            <option value={2000}>2 kb</option>
                            <option value={5000}>5 kb</option>
                            <option value={10000}>10 kb</option>
                            <option value={20000}>20 kb</option>
                        </select>
                    </div>
                    <div className="flex items-center gap-1.5">
                        <label className="text-xs text-slate-500">Paso:</label>
                        <select value={step} onChange={(e) => setStep(parseInt(e.target.value))}
                            className="px-2 py-1.5 border border-slate-200 rounded-lg text-xs">
                            <option value={500}>500 bp</option>
                            <option value={1000}>1 kb</option>
                            <option value={2000}>2 kb</option>
                            <option value={5000}>5 kb</option>
                        </select>
                    </div>
                    <button onClick={handleRecalculate}
                        className="px-4 py-1.5 bg-teal-600 text-white rounded-lg text-xs font-medium hover:bg-teal-700">
                        Recalcular
                    </button>
                </div>
            </div>

            {/* View Mode Tabs */}
            <div className="flex gap-2">
                {[
                    { id: 'gc', label: 'GC Content', icon: '📈' },
                    { id: 'skew', label: 'GC Skew', icon: '📉' },
                    { id: 'both', label: 'Ambos', icon: '📊' },
                ].map(v => (
                    <button key={v.id} onClick={() => setViewMode(v.id)}
                        className={`px-3 py-2 rounded-lg text-xs font-medium transition-all ${viewMode === v.id
                            ? 'bg-teal-600 text-white'
                            : 'bg-white border border-slate-200 text-slate-600 hover:border-teal-300'
                            }`}>
                        {v.icon} {v.label}
                    </button>
                ))}
            </div>

            {loading ? (
                <div className="text-center py-16">
                    <div className="w-12 h-12 border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-4"></div>
                    <p className="text-slate-500">Calculando GC sliding window...</p>
                    <p className="text-xs text-slate-400 mt-1">Ventana: {windowSize.toLocaleString()} bp, Paso: {step.toLocaleString()} bp</p>
                </div>
            ) : data ? (
                <>
                    {/* Stats Cards */}
                    <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
                        <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-4 text-white">
                            <p className="text-teal-100 text-xs">Promedio GC</p>
                            <p className="text-2xl font-bold">{data.average_gc}%</p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-4">
                            <p className="text-xs text-slate-500">GC Mínimo</p>
                            <p className="text-xl font-bold text-red-600">{data.gc_stats.min_gc}%</p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-4">
                            <p className="text-xs text-slate-500">GC Máximo</p>
                            <p className="text-xl font-bold text-green-600">{data.gc_stats.max_gc}%</p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-4">
                            <p className="text-xs text-slate-500">Desv. Estándar</p>
                            <p className="text-xl font-bold text-slate-800">{data.gc_stats.std_gc}%</p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-4">
                            <p className="text-xs text-slate-500">Genoma</p>
                            <p className="text-xl font-bold text-slate-800">{(data.genome_length / 1e6).toFixed(2)} Mb</p>
                        </div>
                    </div>

                    {/* GC Content Chart */}
                    {(viewMode === 'gc' || viewMode === 'both') && (
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <h3 className="font-semibold text-slate-800 mb-1">%GC a lo largo del genoma</h3>
                            <p className="text-xs text-slate-500 mb-4">
                                Contenido GC calculado con ventana deslizante de {(windowSize / 1000).toFixed(0)} kb y paso de {(step / 1000).toFixed(1)} kb.
                                La línea punteada indica el GC promedio ({data.average_gc}%).
                            </p>
                            <ResponsiveContainer width="100%" height={300}>
                                <AreaChart data={sampledGC}>
                                    <defs>
                                        <linearGradient id="gcGrad" x1="0" y1="0" x2="0" y2="1">
                                            <stop offset="5%" stopColor="#14b8a6" stopOpacity={0.3} />
                                            <stop offset="95%" stopColor="#14b8a6" stopOpacity={0.05} />
                                        </linearGradient>
                                    </defs>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                    <XAxis
                                        dataKey="posLabel"
                                        tickCount={10}
                                        tick={{ fontSize: 10, fill: '#64748b' }}
                                        label={{ value: 'Posición genómica', position: 'bottom', offset: -5, fontSize: 11, fill: '#94a3b8' }}
                                    />
                                    <YAxis
                                        domain={['auto', 'auto']}
                                        tick={{ fontSize: 10, fill: '#64748b' }}
                                        label={{ value: 'GC%', angle: -90, position: 'insideLeft', fontSize: 11, fill: '#94a3b8' }}
                                    />
                                    <Tooltip
                                        contentStyle={{ borderRadius: '10px', border: '1px solid #e2e8f0', fontSize: '12px' }}
                                        formatter={(value) => [`${value}%`, 'GC']}
                                    />
                                    <ReferenceLine y={data.average_gc} stroke="#94a3b8" strokeDasharray="5 5"
                                        label={{ value: `Promedio: ${data.average_gc}%`, position: 'right', fontSize: 10, fill: '#94a3b8' }} />
                                    <Area type="monotone" dataKey="gc" stroke="#14b8a6" fill="url(#gcGrad)" strokeWidth={1.5} dot={false} />
                                </AreaChart>
                            </ResponsiveContainer>
                        </div>
                    )}

                    {/* GC Skew Chart */}
                    {(viewMode === 'skew' || viewMode === 'both') && (
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <h3 className="font-semibold text-slate-800 mb-1">GC Skew (G−C)/(G+C)</h3>
                            <p className="text-xs text-slate-500 mb-4">
                                GC Skew indica la asimetría en distribución de G y C entre las dos hebras del DNA.
                                Valores <span className="text-green-600 font-medium">positivos</span> = exceso G (hebra líder),
                                valores <span className="text-red-600 font-medium">negativos</span> = exceso C (hebra rezagada).
                                El punto de cruce indica el origen/terminación de replicación.
                            </p>
                            <ResponsiveContainer width="100%" height={250}>
                                <AreaChart data={sampledSkew}>
                                    <defs>
                                        <linearGradient id="skewGradUp" x1="0" y1="0" x2="0" y2="1">
                                            <stop offset="5%" stopColor="#22c55e" stopOpacity={0.3} />
                                            <stop offset="95%" stopColor="#22c55e" stopOpacity={0} />
                                        </linearGradient>
                                        <linearGradient id="skewGradDown" x1="0" y1="1" x2="0" y2="0">
                                            <stop offset="5%" stopColor="#ef4444" stopOpacity={0.3} />
                                            <stop offset="95%" stopColor="#ef4444" stopOpacity={0} />
                                        </linearGradient>
                                    </defs>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                    <XAxis
                                        dataKey="posLabel"
                                        tickCount={10}
                                        tick={{ fontSize: 10, fill: '#64748b' }}
                                    />
                                    <YAxis
                                        tick={{ fontSize: 10, fill: '#64748b' }}
                                        domain={['auto', 'auto']}
                                        label={{ value: 'GC Skew', angle: -90, position: 'insideLeft', fontSize: 11, fill: '#94a3b8' }}
                                    />
                                    <Tooltip
                                        contentStyle={{ borderRadius: '10px', border: '1px solid #e2e8f0', fontSize: '12px' }}
                                        formatter={(value) => [value.toFixed(4), 'Skew']}
                                    />
                                    <ReferenceLine y={0} stroke="#94a3b8" strokeDasharray="3 3" />
                                    <Area type="monotone" dataKey="skew" stroke="#6366f1" strokeWidth={1.5} dot={false}
                                        fill="url(#skewGradUp)" />
                                </AreaChart>
                            </ResponsiveContainer>
                        </div>
                    )}

                    {/* Info */}
                    <div className="bg-teal-50 border border-teal-100 rounded-xl p-4">
                        <div className="flex gap-3">
                            <div className="w-8 h-8 bg-teal-500 rounded-lg flex items-center justify-center flex-shrink-0">
                                <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
                                    <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
                                </svg>
                            </div>
                            <div className="text-sm text-teal-800 space-y-1">
                                <p><strong>Contenido GC (%GC)</strong>: Proporción de bases Guanina + Citosina. Varía a lo largo del genoma según islas CpG, genes, y regiones codificantes vs. no codificantes.</p>
                                <p><strong>GC Skew</strong>: (G−C)/(G+C) calculado por ventana. Los cambios de signo del GC skew acumulativo indican el <strong>origen (oriC)</strong> y <strong>término de replicación (ter)</strong> en genomas circulares.</p>
                                <p className="text-xs text-teal-700">Ventana actual: <strong>{(windowSize / 1000).toFixed(0)} kb</strong>, Paso: <strong>{(step / 1000).toFixed(1)} kb</strong>, Puntos: <strong>{data.total_windows.toLocaleString()}</strong></p>
                            </div>
                        </div>
                    </div>
                </>
            ) : null}
        </div>
    )
}
