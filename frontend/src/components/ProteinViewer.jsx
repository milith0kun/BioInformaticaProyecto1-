/**
 * ProteinViewer Component — Enhanced
 * Rich protein browser with NCBI links, amino acid composition,
 * hydrophobicity profile, and detailed gene/strand information
 */
import { useState, useEffect, useMemo } from 'react'
import api from '../services/api'
import ProteinStructureViewer from './ProteinStructureViewer'

const AA_PROPERTIES = {
    A: { name: 'Alanina', group: 'Hidrofóbico', color: '#4ade80' },
    I: { name: 'Isoleucina', group: 'Hidrofóbico', color: '#4ade80' },
    L: { name: 'Leucina', group: 'Hidrofóbico', color: '#4ade80' },
    M: { name: 'Metionina', group: 'Hidrofóbico', color: '#4ade80' },
    F: { name: 'Fenilalanina', group: 'Hidrofóbico', color: '#4ade80' },
    W: { name: 'Triptófano', group: 'Hidrofóbico', color: '#4ade80' },
    V: { name: 'Valina', group: 'Hidrofóbico', color: '#4ade80' },
    P: { name: 'Prolina', group: 'Especial', color: '#fbbf24' },
    G: { name: 'Glicina', group: 'Especial', color: '#fbbf24' },
    C: { name: 'Cisteína', group: 'Especial', color: '#fbbf24' },
    D: { name: 'Ácido aspártico', group: 'Cargado (−)', color: '#f87171' },
    E: { name: 'Ácido glutámico', group: 'Cargado (−)', color: '#f87171' },
    K: { name: 'Lisina', group: 'Cargado (+)', color: '#fb923c' },
    R: { name: 'Arginina', group: 'Cargado (+)', color: '#fb923c' },
    H: { name: 'Histidina', group: 'Cargado (+)', color: '#fb923c' },
    N: { name: 'Asparagina', group: 'Polar', color: '#60a5fa' },
    Q: { name: 'Glutamina', group: 'Polar', color: '#60a5fa' },
    S: { name: 'Serina', group: 'Polar', color: '#60a5fa' },
    T: { name: 'Treonina', group: 'Polar', color: '#60a5fa' },
    Y: { name: 'Tirosina', group: 'Polar', color: '#60a5fa' },
    '*': { name: 'Stop', group: 'Stop', color: '#ef4444' },
}

export default function ProteinViewer() {
    const [proteins, setProteins] = useState([])
    const [loading, setLoading] = useState(false)
    const [page, setPage] = useState(1)
    const [totalPages, setTotalPages] = useState(0)
    const [total, setTotal] = useState(0)
    const [search, setSearch] = useState('')
    const [searchInput, setSearchInput] = useState('')
    const [selectedProtein, setSelectedProtein] = useState(null)
    const [detailLoading, setDetailLoading] = useState(false)
    const pageSize = 30

    useEffect(() => { loadProteins() }, [page, search])

    const loadProteins = async () => {
        setLoading(true)
        try {
            const data = await api.getProteins(page, pageSize, search)
            setProteins(data.proteins || [])
            setTotal(data.total)
            setTotalPages(data.total_pages)
        } catch (e) {
            console.error('Error loading proteins:', e)
        } finally {
            setLoading(false)
        }
    }

    const handleSearch = (e) => {
        e.preventDefault()
        setPage(1)
        setSearch(searchInput)
    }

    const selectProtein = async (protein) => {
        setDetailLoading(true)
        setSelectedProtein(protein)
        try {
            const detail = await api.getProteinDetail(protein.protein_id || protein.locus_tag)
            setSelectedProtein(detail)
        } catch (e) {
            console.error('Error loading protein detail:', e)
        } finally {
            setDetailLoading(false)
        }
    }

    const getSizeCategory = (length) => {
        if (length < 100) return { label: 'Pequeña', color: 'bg-blue-100 text-blue-700', desc: '<100 aa' }
        if (length < 300) return { label: 'Media', color: 'bg-emerald-100 text-emerald-700', desc: '100-300 aa' }
        if (length < 600) return { label: 'Grande', color: 'bg-amber-100 text-amber-700', desc: '300-600 aa' }
        return { label: 'Muy grande', color: 'bg-red-100 text-red-700', desc: '>600 aa' }
    }

    // Amino acid composition analysis
    const aaComposition = useMemo(() => {
        if (!selectedProtein?.full_sequence) return null
        const seq = selectedProtein.full_sequence
        const counts = {}
        for (const aa of seq) {
            counts[aa] = (counts[aa] || 0) + 1
        }
        // Group counts
        const groups = { 'Hidrofóbico': 0, 'Polar': 0, 'Cargado (+)': 0, 'Cargado (−)': 0, 'Especial': 0 }
        Object.entries(counts).forEach(([aa, count]) => {
            const prop = AA_PROPERTIES[aa]
            if (prop && groups[prop.group] !== undefined) {
                groups[prop.group] += count
            }
        })
        const totalAA = seq.length
        return {
            counts,
            total: totalAA,
            groups: Object.entries(groups).map(([group, count]) => ({
                group,
                count,
                pct: ((count / totalAA) * 100).toFixed(1)
            })),
            topAA: Object.entries(counts)
                .sort((a, b) => b[1] - a[1])
                .slice(0, 6)
                .map(([aa, count]) => ({
                    aa,
                    count,
                    pct: ((count / totalAA) * 100).toFixed(1),
                    ...(AA_PROPERTIES[aa] || { name: aa, group: 'Otro', color: '#94a3b8' })
                }))
        }
    }, [selectedProtein])

    return (
        <div className="space-y-6">
            {/* Header + Search */}
            <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
                <div>
                    <h2 className="text-xl font-bold text-slate-800">🔬 Proteínas del Genoma</h2>
                    <p className="text-sm text-slate-500">{total.toLocaleString()} proteínas detectadas en secuencias CDS</p>
                </div>
                <form onSubmit={handleSearch} className="flex items-center gap-2 w-full sm:w-auto">
                    <div className="relative flex-1 sm:w-72">
                        <input
                            type="text"
                            value={searchInput}
                            onChange={(e) => setSearchInput(e.target.value)}
                            placeholder="Buscar: gen, proteína, producto, locus..."
                            className="w-full pl-9 pr-4 py-2 border border-slate-200 rounded-xl text-sm focus:outline-none focus:ring-2 focus:ring-teal-400"
                        />
                        <svg className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                        </svg>
                    </div>
                    <button type="submit" className="px-4 py-2 bg-teal-600 text-white rounded-xl text-sm font-medium hover:bg-teal-700">
                        Buscar
                    </button>
                </form>
            </div>

            {/* Protein Grid */}
            {loading ? (
                <div className="text-center py-16">
                    <div className="w-12 h-12 border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-4"></div>
                    <p className="text-slate-500">Extrayendo proteínas del GenBank...</p>
                </div>
            ) : (
                <>
                    <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
                        {proteins.map((protein, i) => {
                            const size = getSizeCategory(protein.length)
                            return (
                                <button
                                    key={i}
                                    onClick={() => selectProtein(protein)}
                                    className={`text-left bg-white rounded-xl border p-4 transition-all hover:shadow-md group ${selectedProtein?.protein_id === protein.protein_id
                                        ? 'border-teal-400 ring-2 ring-teal-100'
                                        : 'border-slate-200 hover:border-teal-300'
                                        }`}
                                >
                                    <div className="flex items-start justify-between gap-2 mb-2">
                                        <div className="min-w-0">
                                            <h4 className="font-mono text-sm font-semibold text-slate-800 truncate">
                                                {protein.protein_id || protein.locus_tag}
                                            </h4>
                                            {protein.gene_name && (
                                                <span className="text-xs text-teal-600 font-medium">{protein.gene_name}</span>
                                            )}
                                        </div>
                                        <span className={`px-2 py-0.5 rounded text-[10px] font-medium flex-shrink-0 ${size.color}`}>
                                            {size.label}
                                        </span>
                                    </div>
                                    <p className="text-xs text-slate-500 line-clamp-2 mb-2">{protein.product || 'Proteína hipotética'}</p>
                                    <div className="flex items-center gap-3 text-[10px] text-slate-400">
                                        <span className="font-medium">{protein.length} aa</span>
                                        <span>{(protein.molecular_weight_approx / 1000).toFixed(1)} kDa</span>
                                        <span className={`flex items-center gap-0.5 font-medium px-1.5 py-0.5 rounded ${protein.strand === 1
                                            ? 'bg-teal-50 text-teal-600'
                                            : 'bg-violet-50 text-violet-600'
                                            }`}>
                                            {protein.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                                        </span>
                                    </div>
                                    <div className="text-[10px] text-slate-400 mt-1 font-mono">
                                        {protein.start?.toLocaleString()} — {protein.end?.toLocaleString()} bp
                                    </div>
                                </button>
                            )
                        })}
                    </div>

                    {/* Pagination */}
                    <div className="flex items-center justify-between">
                        <p className="text-xs text-slate-500">
                            Mostrando {((page - 1) * pageSize) + 1}–{Math.min(page * pageSize, total)} de {total.toLocaleString()}
                        </p>
                        <div className="flex items-center gap-2">
                            <button
                                onClick={() => setPage(p => Math.max(1, p - 1))}
                                disabled={page === 1}
                                className="px-3 py-1.5 bg-white border border-slate-200 rounded-lg text-sm disabled:opacity-40 hover:bg-slate-50"
                            >
                                ← Anterior
                            </button>
                            <span className="text-sm text-slate-600 px-3">{page} / {totalPages}</span>
                            <button
                                onClick={() => setPage(p => Math.min(totalPages, p + 1))}
                                disabled={page >= totalPages}
                                className="px-3 py-1.5 bg-white border border-slate-200 rounded-lg text-sm disabled:opacity-40 hover:bg-slate-50"
                            >
                                Siguiente →
                            </button>
                        </div>
                    </div>
                </>
            )}

            {/* ===== Protein Detail Panel (ahora con visualización de estructuras) ===== */}
            {selectedProtein && (
                <ProteinStructureViewer
                    protein={selectedProtein}
                    onClose={() => setSelectedProtein(null)}
                />
            )}
        </div>
    )
}
