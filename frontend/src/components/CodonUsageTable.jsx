/**
 * CodonUsageTable Component
 * Complete 64-codon usage analysis with RSCU, amino acid grouping, and visualizations
 */
import { useState, useEffect, useMemo } from 'react'
import {
    BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
    ResponsiveContainer, ScatterChart, Scatter, RadarChart,
    Radar, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Cell
} from 'recharts'
import api from '../services/api'

const AA_COLORS = {
    Ala: '#14b8a6', Arg: '#f59e0b', Asn: '#3b82f6', Asp: '#ef4444',
    Cys: '#f97316', Gln: '#06b6d4', Glu: '#dc2626', Gly: '#22c55e',
    His: '#8b5cf6', Ile: '#ec4899', Leu: '#0ea5e9', Lys: '#d946ef',
    Met: '#10b981', Phe: '#6366f1', Pro: '#eab308', Ser: '#2dd4bf',
    Thr: '#a855f7', Trp: '#84cc16', Tyr: '#f43f5e', Val: '#0891b2',
}

export default function CodonUsageTable() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [viewMode, setViewMode] = useState('table') // table | aminoacid | charts
    const [sortBy, setSortBy] = useState('codon') // codon | count | rscu | amino_acid
    const [sortDir, setSortDir] = useState('asc')
    const [filterAA, setFilterAA] = useState('all')

    useEffect(() => {
        loadData()
    }, [])

    const loadData = async () => {
        setLoading(true)
        try {
            const result = await api.getCompleteCodonUsage()
            setData(result)
        } catch (e) {
            console.error('Error loading codon usage:', e)
        } finally {
            setLoading(false)
        }
    }

    const sortedTable = useMemo(() => {
        if (!data?.codon_table) return []
        let filtered = [...data.codon_table]

        if (filterAA !== 'all') {
            filtered = filtered.filter(c => c.amino_acid === filterAA)
        }

        filtered.sort((a, b) => {
            let valA, valB
            switch (sortBy) {
                case 'count': valA = a.count; valB = b.count; break
                case 'rscu': valA = a.rscu; valB = b.rscu; break
                case 'frequency': valA = a.frequency; valB = b.frequency; break
                case 'amino_acid': valA = a.amino_acid; valB = b.amino_acid; break
                default: valA = a.codon; valB = b.codon
            }
            if (typeof valA === 'string') {
                return sortDir === 'asc' ? valA.localeCompare(valB) : valB.localeCompare(valA)
            }
            return sortDir === 'asc' ? valA - valB : valB - valA
        })

        return filtered
    }, [data, sortBy, sortDir, filterAA])

    const aminoAcidChartData = useMemo(() => {
        if (!data?.amino_acid_usage) return []
        return Object.entries(data.amino_acid_usage)
            .map(([aa, info]) => ({
                name: `${info.letter} (${aa})`,
                fraction: info.fraction,
                total: info.total_count,
                codons: info.codons.length,
                fill: AA_COLORS[aa] || '#94a3b8'
            }))
            .sort((a, b) => b.fraction - a.fraction)
    }, [data])

    const rscuChartData = useMemo(() => {
        if (!data?.codon_table) return []
        return data.codon_table
            .filter(c => c.amino_acid !== 'Stop')
            .map(c => ({
                codon: c.codon,
                rscu: c.rscu,
                count: c.count,
                aa: c.amino_acid,
                fill: AA_COLORS[c.amino_acid] || '#94a3b8'
            }))
    }, [data])

    const toggleSort = (field) => {
        if (sortBy === field) {
            setSortDir(d => d === 'asc' ? 'desc' : 'asc')
        } else {
            setSortBy(field)
            setSortDir('desc')
        }
    }

    const getRSCUColor = (rscu) => {
        if (rscu >= 1.5) return 'bg-teal-100 text-teal-800'
        if (rscu >= 1.0) return 'bg-emerald-50 text-emerald-700'
        if (rscu >= 0.5) return 'bg-amber-50 text-amber-700'
        return 'bg-red-50 text-red-700'
    }

    if (loading) {
        return (
            <div className="text-center py-20">
                <div className="w-16 h-16 mx-auto border-4 border-teal-200 border-t-teal-600 rounded-full animate-spin mb-6"></div>
                <h2 className="text-xl font-bold text-slate-800 mb-2">Calculando uso de codones...</h2>
                <p className="text-slate-500">Analizando todas las secuencias CDS del genoma</p>
            </div>
        )
    }

    if (!data) {
        return (
            <div className="text-center py-20">
                <p className="text-slate-500">No hay datos de uso de codones. Execute el análisis primero.</p>
            </div>
        )
    }

    return (
        <div className="space-y-6">
            {/* Summary Stats */}
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-4">
                <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-4 text-white">
                    <p className="text-teal-100 text-xs">Total Codones</p>
                    <p className="text-2xl font-bold">{data.total_codons?.toLocaleString()}</p>
                </div>
                <div className="bg-gradient-to-br from-emerald-500 to-emerald-600 rounded-xl p-4 text-white">
                    <p className="text-emerald-100 text-xs">GC3 Content</p>
                    <p className="text-2xl font-bold">{data.gc3_content}%</p>
                    <p className="text-emerald-200 text-[10px]">3ra posición del codón</p>
                </div>
                <div className="bg-gradient-to-br from-violet-500 to-violet-600 rounded-xl p-4 text-white">
                    <p className="text-violet-100 text-xs">Nc (Codones Efectivos)</p>
                    <p className="text-2xl font-bold">{data.effective_number_of_codons}</p>
                    <p className="text-violet-200 text-[10px]">Wright (1990)</p>
                </div>
                <div className="bg-gradient-to-br from-slate-600 to-slate-700 rounded-xl p-4 text-white">
                    <p className="text-slate-300 text-xs">Aminoácidos Codificados</p>
                    <p className="text-2xl font-bold">20</p>
                    <p className="text-slate-400 text-[10px]">+ 3 stop codons</p>
                </div>
            </div>

            {/* View Mode Tabs */}
            <div className="flex gap-2 flex-wrap">
                {[
                    { id: 'table', label: 'Tabla de Codones', icon: '📋' },
                    { id: 'aminoacid', label: 'Por Aminoácido', icon: '🧪' },
                    { id: 'charts', label: 'Visualizaciones', icon: '📊' },
                ].map(tab => (
                    <button
                        key={tab.id}
                        onClick={() => setViewMode(tab.id)}
                        className={`px-4 py-2 rounded-lg text-sm font-medium transition-all ${viewMode === tab.id
                                ? 'bg-teal-600 text-white shadow-md'
                                : 'bg-white text-slate-600 border border-slate-200 hover:border-teal-300'
                            }`}
                    >
                        {tab.icon} {tab.label}
                    </button>
                ))}
            </div>

            {/* ==================== CODON TABLE ==================== */}
            {viewMode === 'table' && (
                <div className="bg-white rounded-xl border border-slate-200">
                    {/* Filter bar */}
                    <div className="px-4 py-3 border-b border-slate-200 flex items-center gap-3 flex-wrap">
                        <span className="text-xs text-slate-500">Filtrar por aminoácido:</span>
                        <select
                            value={filterAA}
                            onChange={(e) => setFilterAA(e.target.value)}
                            className="px-3 py-1.5 border border-slate-200 rounded-lg text-sm"
                        >
                            <option value="all">Todos (64 codones)</option>
                            {Object.entries(data.amino_acid_usage || {}).map(([aa, info]) => (
                                <option key={aa} value={aa}>{info.letter} — {info.name} ({aa})</option>
                            ))}
                            <option value="Stop">Stop codons</option>
                        </select>
                    </div>

                    <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                            <thead>
                                <tr className="text-xs uppercase tracking-wide text-slate-500 border-b border-slate-200">
                                    {[
                                        { key: 'codon', label: 'Codón' },
                                        { key: 'amino_acid', label: 'Aminoácido' },
                                        { key: 'count', label: 'Conteo' },
                                        { key: 'frequency', label: 'Freq/1000' },
                                        { key: 'rscu', label: 'RSCU' },
                                    ].map(col => (
                                        <th
                                            key={col.key}
                                            onClick={() => toggleSort(col.key)}
                                            className="px-4 py-3 font-medium cursor-pointer hover:bg-slate-50 whitespace-nowrap select-none"
                                        >
                                            {col.label}
                                            {sortBy === col.key && (
                                                <span className="ml-1">{sortDir === 'asc' ? '↑' : '↓'}</span>
                                            )}
                                        </th>
                                    ))}
                                    <th className="px-4 py-3 font-medium whitespace-nowrap">RSCU Barra</th>
                                </tr>
                            </thead>
                            <tbody>
                                {sortedTable.map((row, i) => (
                                    <tr key={i} className="border-b border-slate-50 hover:bg-slate-50">
                                        <td className="px-4 py-2.5">
                                            <span className="font-mono font-bold text-slate-800 bg-slate-100 px-2 py-0.5 rounded">
                                                {row.codon}
                                            </span>
                                        </td>
                                        <td className="px-4 py-2.5">
                                            <span className="flex items-center gap-2">
                                                <span
                                                    className="w-3 h-3 rounded-full flex-shrink-0"
                                                    style={{ backgroundColor: AA_COLORS[row.amino_acid] || '#94a3b8' }}
                                                ></span>
                                                <span className="text-slate-700">{row.amino_acid_name}</span>
                                                <span className="text-slate-400 text-xs font-mono">({row.letter})</span>
                                            </span>
                                        </td>
                                        <td className="px-4 py-2.5 font-mono text-slate-700">{row.count.toLocaleString()}</td>
                                        <td className="px-4 py-2.5 font-mono text-slate-600">{row.frequency.toFixed(2)}</td>
                                        <td className="px-4 py-2.5">
                                            <span className={`inline-block px-2 py-0.5 rounded text-xs font-mono font-semibold ${getRSCUColor(row.rscu)}`}>
                                                {row.rscu.toFixed(3)}
                                            </span>
                                        </td>
                                        <td className="px-4 py-2.5 w-32">
                                            <div className="bg-slate-200 rounded-full h-2 overflow-hidden">
                                                <div
                                                    className="h-full rounded-full transition-all"
                                                    style={{
                                                        width: `${Math.min(row.rscu / 3 * 100, 100)}%`,
                                                        backgroundColor: AA_COLORS[row.amino_acid] || '#94a3b8'
                                                    }}
                                                />
                                            </div>
                                        </td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>

                    {/* RSCU Legend */}
                    <div className="px-4 py-3 border-t border-slate-200 flex flex-wrap gap-3 text-xs text-slate-500">
                        <span className="font-medium">RSCU:</span>
                        <span className="flex items-center gap-1"><span className="w-3 h-2 bg-teal-100 rounded"></span>≥1.5 (preferido)</span>
                        <span className="flex items-center gap-1"><span className="w-3 h-2 bg-emerald-50 rounded border"></span>1.0–1.5 (normal)</span>
                        <span className="flex items-center gap-1"><span className="w-3 h-2 bg-amber-50 rounded border"></span>0.5–1.0 (menos usado)</span>
                        <span className="flex items-center gap-1"><span className="w-3 h-2 bg-red-50 rounded border"></span>&lt;0.5 (poco usado)</span>
                    </div>
                </div>
            )}

            {/* ==================== BY AMINO ACID ==================== */}
            {viewMode === 'aminoacid' && data.amino_acid_usage && (
                <div className="space-y-4">
                    {Object.entries(data.amino_acid_usage).sort((a, b) => b[1].total_count - a[1].total_count).map(([aa, info]) => (
                        <div key={aa} className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                            <div className="flex items-center justify-between px-4 py-3 bg-gradient-to-r from-slate-50 to-white border-b border-slate-100">
                                <div className="flex items-center gap-3">
                                    <span
                                        className="w-8 h-8 rounded-lg flex items-center justify-center text-white font-bold text-sm"
                                        style={{ backgroundColor: AA_COLORS[aa] || '#94a3b8' }}
                                    >
                                        {info.letter}
                                    </span>
                                    <div>
                                        <h4 className="font-semibold text-slate-800">{info.name} ({aa})</h4>
                                        <p className="text-xs text-slate-500">{info.fraction}% del total • {info.total_count.toLocaleString()} codones</p>
                                    </div>
                                </div>
                            </div>

                            <div className="p-4">
                                <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-4 gap-3">
                                    {info.codons.sort((a, b) => b.rscu - a.rscu).map((c, i) => {
                                        const isPreferred = c.rscu >= 1.3
                                        return (
                                            <div
                                                key={i}
                                                className={`rounded-lg p-3 border ${isPreferred ? 'bg-teal-50 border-teal-200' : 'bg-slate-50 border-slate-200'}`}
                                            >
                                                <div className="flex items-center justify-between mb-2">
                                                    <span className="font-mono font-bold text-lg text-slate-800">{c.codon}</span>
                                                    {isPreferred && <span className="text-[10px] text-teal-600 font-semibold">★ PREFERIDO</span>}
                                                </div>
                                                <div className="space-y-1 text-xs">
                                                    <div className="flex justify-between">
                                                        <span className="text-slate-500">Conteo</span>
                                                        <span className="font-mono font-medium">{c.count.toLocaleString()}</span>
                                                    </div>
                                                    <div className="flex justify-between">
                                                        <span className="text-slate-500">RSCU</span>
                                                        <span className={`font-mono font-semibold ${isPreferred ? 'text-teal-700' : 'text-slate-700'}`}>
                                                            {c.rscu.toFixed(3)}
                                                        </span>
                                                    </div>
                                                    <div className="flex justify-between">
                                                        <span className="text-slate-500">Freq/1000</span>
                                                        <span className="font-mono">{c.frequency_per_thousand}</span>
                                                    </div>
                                                </div>
                                                {/* RSCU bar */}
                                                <div className="mt-2 bg-slate-200 rounded-full h-1.5 overflow-hidden">
                                                    <div
                                                        className="h-full rounded-full"
                                                        style={{
                                                            width: `${Math.min(c.rscu / 3 * 100, 100)}%`,
                                                            backgroundColor: AA_COLORS[aa] || '#94a3b8'
                                                        }}
                                                    />
                                                </div>
                                            </div>
                                        )
                                    })}
                                </div>
                            </div>
                        </div>
                    ))}
                </div>
            )}

            {/* ==================== CHARTS ==================== */}
            {viewMode === 'charts' && (
                <div className="space-y-6">
                    {/* Amino acid composition */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-4">Composición de Aminoácidos (%)</h3>
                        <ResponsiveContainer width="100%" height={350}>
                            <BarChart data={aminoAcidChartData} margin={{ bottom: 50 }}>
                                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                <XAxis
                                    dataKey="name"
                                    tick={{ fontSize: 10, fill: '#64748b' }}
                                    angle={-45}
                                    textAnchor="end"
                                    interval={0}
                                />
                                <YAxis tick={{ fontSize: 11, fill: '#64748b' }} label={{ value: '%', position: 'insideLeft' }} />
                                <Tooltip
                                    formatter={(value, name) => [`${value}%`, 'Fracción']}
                                    contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                                />
                                <Bar dataKey="fraction" radius={[3, 3, 0, 0]}>
                                    {aminoAcidChartData.map((entry, i) => (
                                        <Cell key={i} fill={entry.fill} />
                                    ))}
                                </Bar>
                            </BarChart>
                        </ResponsiveContainer>
                    </div>

                    {/* RSCU distribution */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-4">Distribución RSCU por Codón</h3>
                        <ResponsiveContainer width="100%" height={300}>
                            <BarChart data={rscuChartData.slice(0, 61)} margin={{ bottom: 50 }}>
                                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                <XAxis
                                    dataKey="codon"
                                    tick={{ fontSize: 8, fill: '#64748b' }}
                                    angle={-90}
                                    textAnchor="end"
                                    interval={0}
                                />
                                <YAxis tick={{ fontSize: 11, fill: '#64748b' }} />
                                <Tooltip
                                    formatter={(value, name) => [value.toFixed(3), 'RSCU']}
                                    contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                                    labelFormatter={(label) => `Codón: ${label}`}
                                />
                                {/* Reference line at RSCU=1 */}
                                <Bar dataKey="rscu" radius={[2, 2, 0, 0]}>
                                    {rscuChartData.map((entry, i) => (
                                        <Cell key={i} fill={entry.fill} opacity={0.8} />
                                    ))}
                                </Bar>
                            </BarChart>
                        </ResponsiveContainer>
                        <p className="text-xs text-slate-500 mt-2 text-center">
                            RSCU = 1.0 indica uso equitativo entre sinónimos. Valores &gt; 1.0 indican preferencia.
                        </p>
                    </div>

                    {/* Bio Info Box */}
                    <div className="bg-gradient-to-r from-teal-700 to-emerald-700 rounded-xl p-5 text-white">
                        <h3 className="font-semibold mb-3 flex items-center gap-2">
                            <span className="text-2xl">📊</span>
                            Interpretación del Uso de Codones
                        </h3>
                        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-sm text-teal-100">
                            <div>
                                <h4 className="font-medium text-white mb-1">Nc (Codones Efectivos) = {data.effective_number_of_codons}</h4>
                                <p className="text-xs">
                                    {data.effective_number_of_codons < 35
                                        ? 'Sesgo de codones alto. Indica genes altamente expresados con preferencia fuerte.'
                                        : data.effective_number_of_codons < 50
                                            ? 'Sesgo moderado. Típico de genomas bacterianos con expresión balanceada.'
                                            : 'Sesgo bajo. Uso relativamente uniforme de codones sinónimos.'}
                                </p>
                            </div>
                            <div>
                                <h4 className="font-medium text-white mb-1">GC3 = {data.gc3_content}%</h4>
                                <p className="text-xs">
                                    {data.gc3_content > 60
                                        ? 'Alto contenido GC en 3ra posición. Correlacionado con alta expresión génica.'
                                        : data.gc3_content > 40
                                            ? 'GC3 moderado, consistente con el contenido GC global del genoma.'
                                            : 'Bajo GC3. Posible presión mutacional por composición de bases.'}
                                </p>
                            </div>
                        </div>
                    </div>
                </div>
            )}
        </div>
    )
}
