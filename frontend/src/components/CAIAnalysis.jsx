/**
 * CAIAnalysis Component
 * Codon Adaptation Index analysis per gene
 * Identifies highly vs lowly expressed genes
 */
import { useState, useEffect, useMemo } from 'react'
import {
    BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
    ResponsiveContainer, ScatterChart, Scatter, Cell,
    AreaChart, Area
} from 'recharts'
import api from '../services/api'

export default function CAIAnalysis() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [viewTab, setViewTab] = useState('overview')
    const [searchTerm, setSearchTerm] = useState('')
    const [sortBy, setSortBy] = useState('cai_desc')

    useEffect(() => { loadData() }, [])

    const loadData = async () => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getCAIAnalysis(50)
            setData(result)
        } catch (e) {
            setError(e.response?.data?.detail || 'Error calculando CAI')
        } finally {
            setLoading(false)
        }
    }

    // Distribution chart data
    const distData = useMemo(() => {
        if (!data?.cai_distribution) return []
        return Object.entries(data.cai_distribution).map(([range, count]) => ({
            range,
            count,
            fill: range === '0.8-1.0' ? '#14b8a6' :
                range === '0.6-0.8' ? '#22c55e' :
                    range === '0.4-0.6' ? '#f59e0b' :
                        range === '0.2-0.4' ? '#f97316' : '#ef4444'
        }))
    }, [data])

    // Scatter data (CAI vs Length)
    const scatterData = useMemo(() => {
        if (!data?.all_genes) return []
        return data.all_genes.slice(0, 500).map(g => ({
            x: g.length,
            y: g.cai,
            name: g.locus_tag,
            gene: g.gene_name,
            product: g.product,
        }))
    }, [data])

    // Filtered/sorted gene list for table
    const filteredGenes = useMemo(() => {
        if (!data?.all_genes) return []
        let genes = [...data.all_genes]

        if (searchTerm) {
            const q = searchTerm.toLowerCase()
            genes = genes.filter(g =>
                g.locus_tag.toLowerCase().includes(q) ||
                g.gene_name?.toLowerCase().includes(q) ||
                g.product?.toLowerCase().includes(q)
            )
        }

        switch (sortBy) {
            case 'cai_desc': genes.sort((a, b) => b.cai - a.cai); break
            case 'cai_asc': genes.sort((a, b) => a.cai - b.cai); break
            case 'length_desc': genes.sort((a, b) => b.length - a.length); break
            case 'name': genes.sort((a, b) => a.locus_tag.localeCompare(b.locus_tag)); break
        }

        return genes.slice(0, 200)
    }, [data, searchTerm, sortBy])

    // CAI color
    const caiColor = (cai) => {
        if (cai >= 0.8) return '#14b8a6'
        if (cai >= 0.6) return '#22c55e'
        if (cai >= 0.4) return '#f59e0b'
        if (cai >= 0.2) return '#f97316'
        return '#ef4444'
    }

    const caiLabel = (cai) => {
        if (cai >= 0.8) return 'Muy alto'
        if (cai >= 0.6) return 'Alto'
        if (cai >= 0.4) return 'Medio'
        if (cai >= 0.2) return 'Bajo'
        return 'Muy bajo'
    }

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
                <p className="text-slate-500">Calculando Codon Adaptation Index...</p>
                <p className="text-xs text-slate-400 mt-1">Analizando uso de codones de cada gen</p>
            </div>
        )
    }

    if (!data) return null

    return (
        <div className="space-y-6">
            <div>
                <h2 className="text-xl font-bold text-slate-800">📈 CAI — Codon Adaptation Index</h2>
                <p className="text-sm text-slate-500">
                    {data.organism} — {data.total_genes.toLocaleString()} genes analizados
                </p>
            </div>

            {/* Tabs */}
            <div className="flex gap-2 flex-wrap">
                {[
                    { id: 'overview', label: '📊 Resumen' },
                    { id: 'top', label: '🔝 Top Expresados' },
                    { id: 'low', label: '⬇️ Baja Expresión' },
                    { id: 'scatter', label: '📈 CAI vs Longitud' },
                    { id: 'table', label: '📋 Todos los Genes' },
                ].map(t => (
                    <button key={t.id} onClick={() => setViewTab(t.id)}
                        className={`px-3 py-2 rounded-lg text-xs font-medium transition-all ${viewTab === t.id
                            ? 'bg-teal-600 text-white'
                            : 'bg-white border border-slate-200 text-slate-600 hover:border-teal-300'
                            }`}>
                        {t.label}
                    </button>
                ))}
            </div>

            {/* Stats Cards */}
            <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
                <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-4 text-white">
                    <p className="text-teal-100 text-xs">CAI Promedio</p>
                    <p className="text-2xl font-bold">{data.cai_stats.mean}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Mediana</p>
                    <p className="text-xl font-bold text-slate-800">{data.cai_stats.median}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Máximo</p>
                    <p className="text-xl font-bold text-green-600">{data.cai_stats.max}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Mínimo</p>
                    <p className="text-xl font-bold text-red-600">{data.cai_stats.min}</p>
                </div>
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <p className="text-xs text-slate-500">Desv. Estándar</p>
                    <p className="text-xl font-bold text-slate-700">{data.cai_stats.std}</p>
                </div>
            </div>

            {/* OVERVIEW */}
            {viewTab === 'overview' && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                    {/* Distribution */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-1 text-sm">Distribución de CAI</h3>
                        <p className="text-xs text-slate-500 mb-4">Genes agrupados por rango de CAI</p>
                        <ResponsiveContainer width="100%" height={250}>
                            <BarChart data={distData}>
                                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                <XAxis dataKey="range" tick={{ fontSize: 11, fill: '#334155' }} />
                                <YAxis tick={{ fontSize: 10, fill: '#64748b' }} />
                                <Tooltip
                                    formatter={(value) => [value, 'Genes']}
                                    contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                                />
                                <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                                    {distData.map((entry, i) => (
                                        <Cell key={i} fill={entry.fill} />
                                    ))}
                                </Bar>
                            </BarChart>
                        </ResponsiveContainer>
                    </div>

                    {/* CAI Scale Legend */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-4 text-sm">Escala CAI</h3>
                        <div className="space-y-3">
                            {[
                                { range: '0.8 — 1.0', label: 'Muy altamente expresado', color: '#14b8a6', desc: 'Proteínas ribosomales, chaperonas, factores de elongación' },
                                { range: '0.6 — 0.8', label: 'Altamente expresado', color: '#22c55e', desc: 'Proteínas de membrana, enzimas metabólicas principales' },
                                { range: '0.4 — 0.6', label: 'Moderadamente expresado', color: '#f59e0b', desc: 'Genes housekeeping, reguladores' },
                                { range: '0.2 — 0.4', label: 'Baja expresión', color: '#f97316', desc: 'Genes específicos de condición, regulatorios' },
                                { range: '0.0 — 0.2', label: 'Muy baja expresión', color: '#ef4444', desc: 'Genes horizontalmente transferidos, pseudogenes' },
                            ].map((item, i) => (
                                <div key={i} className="flex items-start gap-3">
                                    <div className="w-10 h-6 rounded flex-shrink-0 flex items-center justify-center text-white text-[10px] font-bold"
                                        style={{ backgroundColor: item.color }}>
                                        {item.range.split(' — ')[0]}
                                    </div>
                                    <div>
                                        <p className="text-sm font-medium text-slate-800">{item.label}</p>
                                        <p className="text-xs text-slate-500">{item.desc}</p>
                                    </div>
                                </div>
                            ))}
                        </div>
                    </div>
                </div>
            )}

            {/* TOP EXPRESSED */}
            {viewTab === 'top' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-4 text-sm">
                        🔝 Top {data.top_expressed.length} Genes con CAI más Alto (Altamente Expresados)
                    </h3>
                    <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                            <thead className="bg-slate-50">
                                <tr>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">#</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Locus</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Gen</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">CAI</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Hebra</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Producto</th>
                                </tr>
                            </thead>
                            <tbody>
                                {data.top_expressed.map((g, i) => (
                                    <tr key={i} className={`border-b border-slate-100 ${i % 2 ? 'bg-slate-50/50' : ''}`}>
                                        <td className="px-3 py-2 text-xs text-slate-400">{i + 1}</td>
                                        <td className="px-3 py-2 font-mono text-xs text-teal-700">{g.locus_tag}</td>
                                        <td className="px-3 py-2 font-bold text-xs text-slate-800">{g.gene_name || '—'}</td>
                                        <td className="px-3 py-2">
                                            <span className="font-mono text-xs font-bold px-1.5 py-0.5 rounded"
                                                style={{ color: caiColor(g.cai), backgroundColor: caiColor(g.cai) + '20' }}>
                                                {g.cai}
                                            </span>
                                        </td>
                                        <td className="px-3 py-2 text-xs text-slate-600">{g.length.toLocaleString()} bp</td>
                                        <td className="px-3 py-2">
                                            <span className={`text-[10px] font-medium px-1.5 py-0.5 rounded ${g.strand === 1 ? 'bg-teal-50 text-teal-600' : 'bg-violet-50 text-violet-600'}`}>
                                                {g.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                                            </span>
                                        </td>
                                        <td className="px-3 py-2 text-xs text-slate-600 max-w-[200px] truncate" title={g.product}>{g.product}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            {/* LOW EXPRESSED */}
            {viewTab === 'low' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-4 text-sm">
                        ⬇️ {data.low_expressed.length} Genes con CAI más Bajo (Baja Expresión / Transferencia Horizontal)
                    </h3>
                    <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                            <thead className="bg-slate-50">
                                <tr>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">#</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Locus</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Gen</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">CAI</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Producto</th>
                                </tr>
                            </thead>
                            <tbody>
                                {data.low_expressed.map((g, i) => (
                                    <tr key={i} className={`border-b border-slate-100 ${i % 2 ? 'bg-slate-50/50' : ''}`}>
                                        <td className="px-3 py-2 text-xs text-slate-400">{data.total_genes - data.low_expressed.length + i + 1}</td>
                                        <td className="px-3 py-2 font-mono text-xs text-red-700">{g.locus_tag}</td>
                                        <td className="px-3 py-2 font-bold text-xs text-slate-800">{g.gene_name || '—'}</td>
                                        <td className="px-3 py-2">
                                            <span className="font-mono text-xs font-bold px-1.5 py-0.5 rounded"
                                                style={{ color: caiColor(g.cai), backgroundColor: caiColor(g.cai) + '20' }}>
                                                {g.cai}
                                            </span>
                                        </td>
                                        <td className="px-3 py-2 text-xs text-slate-600">{g.length.toLocaleString()} bp</td>
                                        <td className="px-3 py-2 text-xs text-slate-600 max-w-[250px] truncate" title={g.product}>{g.product}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            {/* SCATTER: CAI vs Length */}
            {viewTab === 'scatter' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-1 text-sm">CAI vs Longitud del Gen</h3>
                    <p className="text-xs text-slate-500 mb-4">
                        Genes cortos con CAI alto suelen ser ribosomales. Genes largos con CAI bajo pueden ser de transferencia horizontal.
                    </p>
                    <ResponsiveContainer width="100%" height={400}>
                        <ScatterChart margin={{ top: 10, right: 20, bottom: 30, left: 10 }}>
                            <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                            <XAxis type="number" dataKey="x" name="Longitud"
                                tick={{ fontSize: 10, fill: '#64748b' }}
                                label={{ value: 'Longitud (bp)', position: 'bottom', offset: 10, fontSize: 11 }}
                            />
                            <YAxis type="number" dataKey="y" name="CAI" domain={[0, 1]}
                                tick={{ fontSize: 10, fill: '#64748b' }}
                                label={{ value: 'CAI', angle: -90, position: 'insideLeft', fontSize: 11 }}
                            />
                            <Tooltip
                                content={({ active, payload }) => {
                                    if (!active || !payload?.[0]) return null
                                    const d = payload[0].payload
                                    return (
                                        <div className="bg-white border border-slate-200 rounded-lg p-2 shadow-lg text-xs">
                                            <p className="font-bold text-slate-800">{d.name} {d.gene ? `(${d.gene})` : ''}</p>
                                            <p className="text-slate-600">{d.product?.slice(0, 60)}</p>
                                            <p className="mt-1">CAI: <strong style={{ color: caiColor(d.y) }}>{d.y}</strong> | Longitud: <strong>{d.x.toLocaleString()} bp</strong></p>
                                        </div>
                                    )
                                }}
                            />
                            <Scatter data={scatterData}>
                                {scatterData.map((entry, i) => (
                                    <Cell key={i} fill={caiColor(entry.y)} fillOpacity={0.6} r={3} />
                                ))}
                            </Scatter>
                        </ScatterChart>
                    </ResponsiveContainer>
                </div>
            )}

            {/* FULL TABLE */}
            {viewTab === 'table' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <div className="flex items-center gap-3 mb-4">
                        <input
                            type="text"
                            value={searchTerm}
                            onChange={(e) => setSearchTerm(e.target.value)}
                            placeholder="Buscar gen, locus, producto..."
                            className="flex-1 px-3 py-2 border border-slate-200 rounded-lg text-sm outline-none focus:ring-2 focus:ring-teal-300"
                        />
                        <select value={sortBy} onChange={(e) => setSortBy(e.target.value)}
                            className="px-3 py-2 border border-slate-200 rounded-lg text-xs">
                            <option value="cai_desc">CAI ↓</option>
                            <option value="cai_asc">CAI ↑</option>
                            <option value="length_desc">Longitud ↓</option>
                            <option value="name">Locus Tag</option>
                        </select>
                    </div>
                    <p className="text-xs text-slate-500 mb-2">
                        Mostrando {filteredGenes.length} de {data.total_genes.toLocaleString()} genes
                    </p>
                    <div className="overflow-x-auto max-h-[500px] overflow-y-auto">
                        <table className="w-full text-sm">
                            <thead className="bg-slate-50 sticky top-0">
                                <tr>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Locus</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Gen</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">CAI</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Nivel</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Hebra</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Producto</th>
                                </tr>
                            </thead>
                            <tbody>
                                {filteredGenes.map((g, i) => (
                                    <tr key={i} className={`border-b border-slate-100 hover:bg-teal-50/30 ${i % 2 ? 'bg-slate-50/50' : ''}`}>
                                        <td className="px-3 py-1.5 font-mono text-xs">{g.locus_tag}</td>
                                        <td className="px-3 py-1.5 font-bold text-xs">{g.gene_name || '—'}</td>
                                        <td className="px-3 py-1.5">
                                            <span className="font-mono text-xs font-bold px-1.5 py-0.5 rounded"
                                                style={{ color: caiColor(g.cai), backgroundColor: caiColor(g.cai) + '15' }}>
                                                {g.cai}
                                            </span>
                                        </td>
                                        <td className="px-3 py-1.5 text-[10px]" style={{ color: caiColor(g.cai) }}>
                                            {caiLabel(g.cai)}
                                        </td>
                                        <td className="px-3 py-1.5 text-xs text-slate-600">{g.length.toLocaleString()}</td>
                                        <td className="px-3 py-1.5">
                                            <span className={`text-[10px] px-1 py-0.5 rounded ${g.strand === 1 ? 'bg-teal-50 text-teal-600' : 'bg-violet-50 text-violet-600'}`}>
                                                {g.strand === 1 ? '→' : '←'}
                                            </span>
                                        </td>
                                        <td className="px-3 py-1.5 text-xs text-slate-500 max-w-[200px] truncate" title={g.product}>{g.product}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
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
                        <p><strong>CAI (Codon Adaptation Index)</strong>: Mide cuán adaptados están los codones de un gen al uso preferido del organismo. Un CAI alto indica que el gen usa los codones óptimos → mayor expresión.</p>
                        <p className="text-xs text-teal-700">
                            <strong>Referencia</strong>: Se usan las proteínas ribosomales y factores de elongación como set de referencia (genes altamente expresados en todos los organismos).
                            <strong> w(codon)</strong> = freq / max_freq_sinónimo. <strong>CAI</strong> = media geométrica de los w de todos los codones del gen.
                        </p>
                        <p className="text-xs text-teal-700">
                            Genes con CAI bajo y uso atípico de codones son candidatos a <strong>transferencia horizontal</strong> (HGT).
                        </p>
                    </div>
                </div>
            </div>
        </div>
    )
}
