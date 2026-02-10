/**
 * RNAAnalysis Component
 * Displays tRNA and rRNA analysis from GenBank file
 * Including amino acid coverage, strand distribution, and gene tables
 */
import { useState, useEffect, useMemo } from 'react'
import {
    BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
    ResponsiveContainer, PieChart, Pie, Cell
} from 'recharts'
import api from '../services/api'

const COLORS = ['#14b8a6', '#10b981', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6',
    '#ec4899', '#06b6d4', '#f97316', '#22c55e', '#a855f7', '#3b82f6']

export default function RNAAnalysis() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [viewTab, setViewTab] = useState('overview')

    useEffect(() => { loadData() }, [])

    const loadData = async () => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getRNAAnalysis()
            setData(result)
        } catch (e) {
            setError(e.response?.data?.detail || 'Error cargando datos RNA')
        } finally {
            setLoading(false)
        }
    }

    // AA coverage chart data
    const aaCoverageData = useMemo(() => {
        if (!data?.amino_acid_coverage) return []
        return Object.entries(data.amino_acid_coverage)
            .sort((a, b) => b[1] - a[1])
            .map(([aa, count]) => ({ aa, count }))
    }, [data])

    // rRNA types chart data
    const rrnaPieData = useMemo(() => {
        if (!data?.rrna_types) return []
        return Object.entries(data.rrna_types).map(([type, count]) => ({ name: type, value: count }))
    }, [data])

    if (error) {
        return (
            <div className="bg-amber-50 border border-amber-200 rounded-xl p-6 text-center">
                <p className="text-amber-700 font-medium">⚠️ {error}</p>
                <p className="text-amber-600 text-sm mt-2">Active un genoma primero para analizar tRNA/rRNA.</p>
            </div>
        )
    }

    if (loading) {
        return (
            <div className="text-center py-16">
                <div className="w-12 h-12 border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-4"></div>
                <p className="text-slate-500">Analizando genes tRNA y rRNA...</p>
            </div>
        )
    }

    if (!data) return null

    return (
        <div className="space-y-6">
            <div>
                <h2 className="text-xl font-bold text-slate-800">🧪 Análisis de tRNA y rRNA</h2>
                <p className="text-sm text-slate-500">{data.organism} — {data.total_trna} tRNA, {data.total_rrna} rRNA genes</p>
            </div>

            {/* Tabs */}
            <div className="flex gap-2">
                {[
                    { id: 'overview', label: 'Resumen', icon: '📊' },
                    { id: 'trna', label: `tRNA (${data.total_trna})`, icon: '🧬' },
                    { id: 'rrna', label: `rRNA (${data.total_rrna})`, icon: '🔬' },
                ].map(t => (
                    <button key={t.id} onClick={() => setViewTab(t.id)}
                        className={`px-4 py-2 rounded-lg text-sm font-medium transition-all ${viewTab === t.id
                            ? 'bg-teal-600 text-white'
                            : 'bg-white border border-slate-200 text-slate-600 hover:border-teal-300'
                            }`}>
                        {t.icon} {t.label}
                    </button>
                ))}
            </div>

            {/* OVERVIEW */}
            {viewTab === 'overview' && (
                <>
                    {/* Stats */}
                    <div className="grid grid-cols-2 sm:grid-cols-4 gap-4">
                        <div className="bg-gradient-to-br from-teal-500 to-teal-600 rounded-xl p-5 text-white">
                            <p className="text-teal-100 text-sm">Total tRNA</p>
                            <p className="text-3xl font-bold mt-1">{data.total_trna}</p>
                        </div>
                        <div className="bg-gradient-to-br from-violet-500 to-violet-600 rounded-xl p-5 text-white">
                            <p className="text-violet-100 text-sm">Total rRNA</p>
                            <p className="text-3xl font-bold mt-1">{data.total_rrna}</p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <p className="text-xs text-slate-500">Aminoácidos cubiertos</p>
                            <p className="text-3xl font-bold text-slate-800 mt-1">
                                {Object.keys(data.amino_acid_coverage).length}
                                <span className="text-sm font-normal text-slate-400 ml-1">/ 20</span>
                            </p>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <p className="text-xs text-slate-500">Tipos rRNA</p>
                            <p className="text-3xl font-bold text-slate-800 mt-1">{Object.keys(data.rrna_types).length}</p>
                        </div>
                    </div>

                    {/* Strand Distribution */}
                    <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <h3 className="font-semibold text-slate-800 mb-3 text-sm">tRNA por Hebra</h3>
                            <div className="flex items-center gap-4">
                                <div className="flex-1">
                                    <div className="flex items-center gap-2 mb-1.5">
                                        <span className="w-3 h-3 bg-teal-400 rounded"></span>
                                        <span className="text-sm text-slate-600">→ 5'→3' Forward</span>
                                        <span className="ml-auto font-bold text-slate-800">{data.trna_by_strand.forward}</span>
                                    </div>
                                    <div className="flex items-center gap-2">
                                        <span className="w-3 h-3 bg-violet-400 rounded"></span>
                                        <span className="text-sm text-slate-600">← 3'→5' Reverse</span>
                                        <span className="ml-auto font-bold text-slate-800">{data.trna_by_strand.reverse}</span>
                                    </div>
                                    <div className="w-full h-3 bg-slate-100 rounded-full overflow-hidden flex mt-3">
                                        <div className="bg-teal-400 h-full" style={{ width: `${(data.trna_by_strand.forward / data.total_trna * 100)}%` }}></div>
                                        <div className="bg-violet-400 h-full" style={{ width: `${(data.trna_by_strand.reverse / data.total_trna * 100)}%` }}></div>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <h3 className="font-semibold text-slate-800 mb-3 text-sm">rRNA por Hebra</h3>
                            <div className="flex items-center gap-4">
                                <div className="flex-1">
                                    <div className="flex items-center gap-2 mb-1.5">
                                        <span className="w-3 h-3 bg-teal-400 rounded"></span>
                                        <span className="text-sm text-slate-600">→ 5'→3' Forward</span>
                                        <span className="ml-auto font-bold text-slate-800">{data.rrna_by_strand.forward}</span>
                                    </div>
                                    <div className="flex items-center gap-2">
                                        <span className="w-3 h-3 bg-violet-400 rounded"></span>
                                        <span className="text-sm text-slate-600">← 3'→5' Reverse</span>
                                        <span className="ml-auto font-bold text-slate-800">{data.rrna_by_strand.reverse}</span>
                                    </div>
                                    {data.total_rrna > 0 && (
                                        <div className="w-full h-3 bg-slate-100 rounded-full overflow-hidden flex mt-3">
                                            <div className="bg-teal-400 h-full" style={{ width: `${(data.rrna_by_strand.forward / data.total_rrna * 100)}%` }}></div>
                                            <div className="bg-violet-400 h-full" style={{ width: `${(data.rrna_by_strand.reverse / data.total_rrna * 100)}%` }}></div>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </div>
                    </div>

                    {/* AA Coverage */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="font-semibold text-slate-800 mb-1">Cobertura de Aminoácidos por tRNA</h3>
                        <p className="text-xs text-slate-500 mb-4">Número de genes tRNA para cada aminoácido. Los aminoácidos más comunes tienen más copias de tRNA.</p>
                        <ResponsiveContainer width="100%" height={280}>
                            <BarChart data={aaCoverageData}>
                                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                <XAxis dataKey="aa" tick={{ fontSize: 11, fill: '#334155', fontWeight: 600 }} />
                                <YAxis tick={{ fontSize: 10, fill: '#64748b' }} />
                                <Tooltip
                                    formatter={(value) => [value, 'Genes tRNA']}
                                    contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0' }}
                                />
                                <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                                    {aaCoverageData.map((_, i) => (
                                        <Cell key={i} fill={COLORS[i % COLORS.length]} />
                                    ))}
                                </Bar>
                            </BarChart>
                        </ResponsiveContainer>
                    </div>

                    {/* rRNA Types Pie */}
                    {rrnaPieData.length > 0 && (
                        <div className="bg-white rounded-xl border border-slate-200 p-5">
                            <h3 className="font-semibold text-slate-800 mb-4">Tipos de rRNA</h3>
                            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                                <ResponsiveContainer width="100%" height={250}>
                                    <PieChart>
                                        <Pie data={rrnaPieData} dataKey="value" nameKey="name" cx="50%" cy="50%"
                                            outerRadius={90} innerRadius={50} paddingAngle={3}>
                                            {rrnaPieData.map((_, i) => (
                                                <Cell key={i} fill={COLORS[i % COLORS.length]} />
                                            ))}
                                        </Pie>
                                        <Tooltip />
                                    </PieChart>
                                </ResponsiveContainer>
                                <div className="space-y-2 flex flex-col justify-center">
                                    {rrnaPieData.map((r, i) => (
                                        <div key={i} className="flex items-center gap-3 text-sm">
                                            <span className="w-3 h-3 rounded flex-shrink-0" style={{ backgroundColor: COLORS[i % COLORS.length] }}></span>
                                            <span className="text-slate-700">{r.name}</span>
                                            <span className="ml-auto font-bold text-slate-800">{r.value}</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>
                    )}
                </>
            )}

            {/* tRNA TABLE */}
            {viewTab === 'trna' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-4">
                        Genes tRNA ({data.total_trna})
                    </h3>
                    <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                            <thead className="bg-slate-50">
                                <tr>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Locus</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Producto</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Aminoácido</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Inicio</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Fin</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Hebra</th>
                                </tr>
                            </thead>
                            <tbody>
                                {data.trna_genes.map((trna, i) => (
                                    <tr key={i} className={`border-b border-slate-100 ${i % 2 === 0 ? '' : 'bg-slate-50/50'}`}>
                                        <td className="px-3 py-2 font-mono text-xs">{trna.locus_tag}</td>
                                        <td className="px-3 py-2 text-xs">{trna.product}</td>
                                        <td className="px-3 py-2 font-medium text-teal-700">{trna.amino_acid}</td>
                                        <td className="px-3 py-2 font-mono text-xs">{trna.start.toLocaleString()}</td>
                                        <td className="px-3 py-2 font-mono text-xs">{trna.end.toLocaleString()}</td>
                                        <td className="px-3 py-2 text-xs">{trna.length} bp</td>
                                        <td className="px-3 py-2">
                                            <span className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-lg text-[10px] font-medium ${trna.strand === 1
                                                ? 'bg-teal-50 text-teal-600'
                                                : 'bg-violet-50 text-violet-600'
                                                }`}>
                                                {trna.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                                            </span>
                                        </td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            {/* rRNA TABLE */}
            {viewTab === 'rrna' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-4">
                        Genes rRNA ({data.total_rrna})
                    </h3>
                    <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                            <thead className="bg-slate-50">
                                <tr>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Locus</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Producto</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Inicio</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Fin</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Hebra</th>
                                </tr>
                            </thead>
                            <tbody>
                                {data.rrna_genes.map((rrna, i) => (
                                    <tr key={i} className={`border-b border-slate-100 ${i % 2 === 0 ? '' : 'bg-slate-50/50'}`}>
                                        <td className="px-3 py-2 font-mono text-xs">{rrna.locus_tag}</td>
                                        <td className="px-3 py-2 text-xs font-medium text-violet-700">{rrna.product}</td>
                                        <td className="px-3 py-2 font-mono text-xs">{rrna.start.toLocaleString()}</td>
                                        <td className="px-3 py-2 font-mono text-xs">{rrna.end.toLocaleString()}</td>
                                        <td className="px-3 py-2 text-xs">{rrna.length.toLocaleString()} bp</td>
                                        <td className="px-3 py-2">
                                            <span className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-lg text-[10px] font-medium ${rrna.strand === 1
                                                ? 'bg-teal-50 text-teal-600'
                                                : 'bg-violet-50 text-violet-600'
                                                }`}>
                                                {rrna.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                                            </span>
                                        </td>
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
                        <p><strong>tRNA</strong>: RNA de transferencia — transporta aminoácidos al ribosoma durante la traducción. Cada tRNA reconoce un codón específico del mRNA mediante apareamiento de bases con su anticodón.</p>
                        <p><strong>rRNA</strong>: RNA ribosomal — componente estructural y catalítico de los ribosomas. En procariotas: 5S, 16S y 23S rRNA. El 16S rRNA es la base de la taxonomía bacteriana.</p>
                        <p className="text-xs text-teal-700">Los genomas bacterianos típicos tienen 1-15 operones rRNA y 60-90 genes tRNA.</p>
                    </div>
                </div>
            </div>
        </div>
    )
}
