/**
 * FunctionalCategories Component
 * Classify genes into COG-like functional categories
 * with horizontal bar charts and expandable gene lists
 */
import { useState, useEffect, useMemo } from 'react'
import {
    BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
    ResponsiveContainer, PieChart, Pie, Cell
} from 'recharts'
import api from '../services/api'

export default function FunctionalCategories() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [expandedCat, setExpandedCat] = useState(null)
    const [viewMode, setViewMode] = useState('chart') // 'chart' | 'table'

    useEffect(() => { loadData() }, [])

    const loadData = async () => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getFunctionalCategories()
            setData(result)
        } catch (e) {
            setError(e.response?.data?.detail || 'Error clasificando genes')
        } finally {
            setLoading(false)
        }
    }

    const chartData = useMemo(() => {
        if (!data?.categories) return []
        return data.categories.map(c => ({
            name: `${c.code} - ${c.name}`,
            code: c.code,
            count: c.count,
            fill: c.color,
            shortName: c.name.length > 25 ? c.name.slice(0, 22) + '...' : c.name,
        }))
    }, [data])

    const pieData = useMemo(() => {
        if (!data?.categories) return []
        const top = data.categories.slice(0, 8)
        const others = data.categories.slice(8)
        const result = top.map(c => ({ name: c.name, value: c.count, color: c.color }))
        if (others.length > 0) {
            result.push({
                name: 'Otros',
                value: others.reduce((s, c) => s + c.count, 0) + (data.uncategorized || 0),
                color: '#94a3b8'
            })
        }
        return result
    }, [data])

    if (error) {
        return (
            <div className="bg-amber-50 border border-amber-200 rounded-xl p-6 text-center">
                <p className="text-amber-700 font-medium">⚠️ {error}</p>
                <p className="text-amber-600 text-sm mt-2">Active un genoma primero.</p>
            </div>
        )
    }

    if (loading) {
        return (
            <div className="text-center py-16">
                <div className="w-12 h-12 border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-4"></div>
                <p className="text-slate-500">Clasificando genes en categorías funcionales...</p>
            </div>
        )
    }

    if (!data) return null

    const categorizedPct = data.total_cds > 0 ? ((data.categorized / data.total_cds) * 100).toFixed(1) : 0

    return (
        <div className="space-y-6">
            <div>
                <h2 className="text-xl font-bold text-slate-800">📋 Categorías Funcionales (COG)</h2>
                <p className="text-sm text-slate-500">
                    {data.organism} — {data.categorized.toLocaleString()} de {data.total_cds.toLocaleString()} CDS clasificados ({categorizedPct}%)
                </p>
            </div>

            {/* Stats */}
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-4 text-white">
                    <p className="text-teal-100 text-xs">Total CDS</p>
                    <p className="text-2xl font-bold">{data.total_cds.toLocaleString()}</p>
                </div>
                <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-4 text-white">
                    <p className="text-emerald-100 text-xs">Clasificados</p>
                    <p className="text-2xl font-bold">{data.categorized.toLocaleString()}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Sin Clasificar</p>
                    <p className="text-2xl font-bold text-slate-600">{data.uncategorized.toLocaleString()}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Categorías</p>
                    <p className="text-2xl font-bold text-slate-800">{data.categories.length}</p>
                </div>
            </div>

            {/* Coverage bar */}
            <div className="bg-white rounded-xl border border-slate-200 p-4">
                <div className="flex justify-between items-center mb-2">
                    <span className="text-xs text-slate-500">Cobertura de clasificación</span>
                    <span className="text-xs font-bold text-slate-700">{categorizedPct}%</span>
                </div>
                <div className="w-full h-4 bg-slate-100 rounded-full overflow-hidden flex">
                    {data.top_categories.map((cat, i) => (
                        <div key={i} className="h-full transition-all"
                            style={{
                                width: `${(cat.count / data.total_cds * 100)}%`,
                                backgroundColor: cat.color,
                            }}
                            title={`${cat.code}: ${cat.name} (${cat.count})`}
                        ></div>
                    ))}
                </div>
                <div className="flex flex-wrap gap-2 mt-2">
                    {data.top_categories.slice(0, 6).map((cat, i) => (
                        <span key={i} className="flex items-center gap-1 text-[10px] text-slate-600">
                            <span className="w-2 h-2 rounded" style={{ backgroundColor: cat.color }}></span>
                            {cat.code}
                        </span>
                    ))}
                </div>
            </div>

            {/* View toggle */}
            <div className="flex gap-2">
                {[
                    { id: 'chart', label: '📊 Gráficos' },
                    { id: 'table', label: '📋 Tabla Detallada' },
                ].map(v => (
                    <button key={v.id} onClick={() => setViewMode(v.id)}
                        className={`px-3 py-2 rounded-lg text-xs font-medium transition-all ${viewMode === v.id
                            ? 'bg-teal-600 text-white'
                            : 'bg-white border border-slate-200 text-slate-600 hover:border-teal-300'
                            }`}>
                        {v.label}
                    </button>
                ))}
            </div>

            {viewMode === 'chart' && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                    {/* Horizontal Bar Chart */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-4 text-sm">Genes por Categoría</h3>
                        <ResponsiveContainer width="100%" height={chartData.length * 30 + 40}>
                            <BarChart data={chartData} layout="vertical" margin={{ left: 10, right: 20 }}>
                                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" horizontal={false} />
                                <XAxis type="number" tick={{ fontSize: 10, fill: '#64748b' }} />
                                <YAxis type="category" dataKey="code" tick={{ fontSize: 11, fill: '#334155', fontWeight: 600 }} width={25} />
                                <Tooltip
                                    formatter={(value, _, props) => [value, props.payload.name]}
                                    contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0', fontSize: '12px' }}
                                />
                                <Bar dataKey="count" radius={[0, 4, 4, 0]}>
                                    {chartData.map((entry, i) => (
                                        <Cell key={i} fill={entry.fill} />
                                    ))}
                                </Bar>
                            </BarChart>
                        </ResponsiveContainer>
                    </div>

                    {/* Pie Chart */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-4 text-sm">Distribución</h3>
                        <ResponsiveContainer width="100%" height={300}>
                            <PieChart>
                                <Pie data={pieData} dataKey="value" nameKey="name"
                                    cx="50%" cy="50%" outerRadius={110} innerRadius={60} paddingAngle={2}>
                                    {pieData.map((entry, i) => (
                                        <Cell key={i} fill={entry.color} />
                                    ))}
                                </Pie>
                                <Tooltip formatter={(value) => [value, 'Genes']} />
                            </PieChart>
                        </ResponsiveContainer>
                        <div className="grid grid-cols-2 gap-1 mt-2">
                            {pieData.map((p, i) => (
                                <span key={i} className="flex items-center gap-1.5 text-[10px] text-slate-600">
                                    <span className="w-2 h-2 rounded flex-shrink-0" style={{ backgroundColor: p.color }}></span>
                                    <span className="truncate">{p.name}</span>
                                    <span className="ml-auto font-bold">{p.value}</span>
                                </span>
                            ))}
                        </div>
                    </div>
                </div>
            )}

            {/* Detailed Table */}
            {viewMode === 'table' && (
                <div className="space-y-2">
                    {data.categories.map((cat) => (
                        <div key={cat.code} className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                            <button onClick={() => setExpandedCat(expandedCat === cat.code ? null : cat.code)}
                                className="w-full flex items-center gap-3 p-4 text-left hover:bg-slate-50 transition-colors">
                                <span className="w-8 h-8 rounded-lg flex items-center justify-center text-white font-bold text-xs"
                                    style={{ backgroundColor: cat.color }}>
                                    {cat.code}
                                </span>
                                <div className="flex-1 min-w-0">
                                    <p className="text-sm font-medium text-slate-800">{cat.name}</p>
                                    <p className="text-xs text-slate-500">{cat.count} genes ({(cat.count / data.total_cds * 100).toFixed(1)}%)</p>
                                </div>
                                {/* Mini bar */}
                                <div className="w-32 h-2 bg-slate-100 rounded-full overflow-hidden hidden sm:block">
                                    <div className="h-full rounded-full" style={{
                                        width: `${Math.min(100, (cat.count / (data.categories[0]?.count || 1)) * 100)}%`,
                                        backgroundColor: cat.color
                                    }}></div>
                                </div>
                                <span className="font-bold text-slate-800 text-sm w-12 text-right">{cat.count}</span>
                                <span className={`text-slate-400 transition-transform ${expandedCat === cat.code ? 'rotate-180' : ''}`}>▼</span>
                            </button>

                            {expandedCat === cat.code && cat.sample_genes && (
                                <div className="border-t border-slate-100 p-4 bg-slate-50/50">
                                    <p className="text-xs text-slate-500 mb-2">
                                        Mostrando {Math.min(20, cat.sample_genes.length)} de {cat.count} genes
                                    </p>
                                    <div className="grid grid-cols-1 sm:grid-cols-2 gap-1">
                                        {cat.sample_genes.map((g, i) => (
                                            <div key={i} className="flex items-baseline gap-2 text-xs py-0.5">
                                                <code className="text-teal-700 font-mono bg-teal-50 px-1 rounded">{g.locus_tag}</code>
                                                {g.gene_name && <span className="font-bold text-slate-700">{g.gene_name}</span>}
                                                <span className="text-slate-500 truncate">{g.product}</span>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                            )}
                        </div>
                    ))}
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
                        <p><strong>COG (Clusters of Orthologous Groups)</strong>: Sistema de clasificación funcional desarrollado por NCBI. Agrupa proteínas en categorías según su función (metabolismo, transporte, regulación, etc.).</p>
                        <p className="text-xs text-teal-700">
                            La clasificación se realiza por análisis de keywords en las descripciones de productos génicos.
                            Las categorías principales son: <strong>J</strong> (Traducción), <strong>K</strong> (Transcripción),
                            <strong> L</strong> (Replicación), <strong>C</strong> (Energía), <strong>E</strong> (Aminoácidos),
                            <strong> G</strong> (Carbohidratos), <strong>M</strong> (Membrana), <strong>T</strong> (Señales).
                        </p>
                    </div>
                </div>
            </div>
        </div>
    )
}
