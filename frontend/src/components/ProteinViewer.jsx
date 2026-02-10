/**
 * ProteinViewer Component — Redesigned
 * Modern protein browser with infinite scrolling and enhanced navigation
 */
import { useState, useEffect, useMemo, useCallback } from 'react'
import api from '../services/api'
import ProteinStructureViewer from './ProteinStructureViewer'

const AA_PROPERTIES = {
    A: { name: 'Alanina', group: 'Hidrofóbico', color: '#10b981' },
    I: { name: 'Isoleucina', group: 'Hidrofóbico', color: '#10b981' },
    L: { name: 'Leucina', group: 'Hidrofóbico', color: '#10b981' },
    M: { name: 'Metionina', group: 'Hidrofóbico', color: '#10b981' },
    F: { name: 'Fenilalanina', group: 'Hidrofóbico', color: '#10b981' },
    W: { name: 'Triptófano', group: 'Hidrofóbico', color: '#10b981' },
    V: { name: 'Valina', group: 'Hidrofóbico', color: '#10b981' },
    P: { name: 'Prolina', group: 'Especial', color: '#f59e0b' },
    G: { name: 'Glicina', group: 'Especial', color: '#f59e0b' },
    C: { name: 'Cisteína', group: 'Especial', color: '#f59e0b' },
    D: { name: 'Ácido aspártico', group: 'Cargado (−)', color: '#ef4444' },
    E: { name: 'Ácido glutámico', group: 'Cargado (−)', color: '#ef4444' },
    K: { name: 'Lisina', group: 'Cargado (+)', color: '#3b82f6' },
    R: { name: 'Arginina', group: 'Cargado (+)', color: '#3b82f6' },
    H: { name: 'Histidina', group: 'Cargado (+)', color: '#3b82f6' },
    N: { name: 'Asparagina', group: 'Polar', color: '#8b5cf6' },
    Q: { name: 'Glutamina', group: 'Polar', color: '#8b5cf6' },
    S: { name: 'Serina', group: 'Polar', color: '#8b5cf6' },
    T: { name: 'Treonina', group: 'Polar', color: '#8b5cf6' },
    Y: { name: 'Tirosina', group: 'Polar', color: '#8b5cf6' },
    '*': { name: 'Stop', group: 'Stop', color: '#ef4444' },
}

export default function ProteinViewer() {
    const [proteins, setProteins] = useState([])
    const [loading, setLoading] = useState(false)
    const [loadingMore, setLoadingMore] = useState(false)
    const [total, setTotal] = useState(0)
    const [page, setPage] = useState(1)
    const [totalPages, setTotalPages] = useState(1)
    const [search, setSearch] = useState('')
    const [searchInput, setSearchInput] = useState('')
    const [selectedProteinIndex, setSelectedProteinIndex] = useState(null)
    const [detailLoading, setDetailLoading] = useState(false)

    const pageSize = 50

    // Initial load
    useEffect(() => {
        loadProteins(1, true)
    }, [search])

    const loadProteins = async (pageNum, isInitial = false) => {
        if (isInitial) setLoading(true)
        else setLoadingMore(true)
        
        try {
            const data = await api.getProteins(pageNum, pageSize, search)
            if (isInitial) {
                setProteins(data.proteins || [])
            } else {
                setProteins(prev => [...prev, ...(data.proteins || [])])
            }
            setTotal(data.total_count || 0)
            setTotalPages(data.total_pages || 1)
            setPage(pageNum)
        } catch (e) {
            console.error('Error loading proteins:', e)
        } finally {
            setLoading(false)
            setLoadingMore(false)
        }
    }

    const handleLoadMore = () => {
        if (page < totalPages && !loadingMore) {
            loadProteins(page + 1)
        }
    }

    const handleSearch = (e) => {
        e.preventDefault()
        setSearch(searchInput)
        setSelectedProteinIndex(null)
    }

    const selectProtein = async (index) => {
        setDetailLoading(true)
        setSelectedProteinIndex(index)
        try {
            const protein = proteins[index]
            const detail = await api.getProteinDetail(protein.protein_id || protein.locus_tag)
            // Update the protein in the list with full details
            const updatedProteins = [...proteins]
            updatedProteins[index] = detail
            setProteins(updatedProteins)
        } catch (e) {
            console.error('Error loading protein detail:', e)
        } finally {
            setDetailLoading(false)
        }
    }

    const navigateProtein = useCallback((direction) => {
        if (selectedProteinIndex === null) return

        const newIndex = direction === 'next'
            ? Math.min(selectedProteinIndex + 1, proteins.length - 1)
            : Math.max(selectedProteinIndex - 1, 0)

        if (newIndex !== selectedProteinIndex) {
            selectProtein(newIndex)
        }
    }, [selectedProteinIndex, proteins.length])

    const getSizeCategory = (length) => {
        if (length < 100) return { label: 'Pequeña', color: 'from-cyan-500 to-blue-500', ring: 'ring-cyan-400/30', text: 'text-cyan-400' }
        if (length < 300) return { label: 'Media', color: 'from-emerald-500 to-teal-500', ring: 'ring-emerald-400/30', text: 'text-emerald-400' }
        if (length < 600) return { label: 'Grande', color: 'from-amber-500 to-orange-500', ring: 'ring-amber-400/30', text: 'text-amber-400' }
        return { label: 'Muy Grande', color: 'from-rose-500 to-pink-500', ring: 'ring-rose-400/30', text: 'text-rose-400' }
    }

    const selectedProtein = selectedProteinIndex !== null ? proteins[selectedProteinIndex] : null

    return (
        <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900">
            <div className="max-w-7xl mx-auto p-6 space-y-8">
                {/* Header Section */}
                <div className="flex flex-col lg:flex-row items-start lg:items-center justify-between gap-6 bg-slate-800/40 p-8 rounded-3xl border border-slate-700/50 backdrop-blur-md">
                    <div className="space-y-3">
                        <h2 className="text-3xl lg:text-5xl font-black tracking-tight bg-gradient-to-r from-cyan-400 via-blue-500 to-indigo-500 bg-clip-text text-transparent">
                            PROTEOME BROWSER
                        </h2>
                        <div className="flex items-center gap-3">
                            <div className="px-3 py-1 bg-cyan-500/10 rounded-full border border-cyan-500/20">
                                <p className="text-cyan-400 text-xs font-bold uppercase tracking-widest">
                                    {total.toLocaleString()} Proteínas detectadas
                                </p>
                            </div>
                            {search && (
                                <div className="px-3 py-1 bg-amber-500/10 rounded-full border border-amber-500/20">
                                    <p className="text-amber-400 text-xs font-bold uppercase tracking-widest">
                                        Filtrado: {search}
                                    </p>
                                </div>
                            )}
                        </div>
                    </div>

                    {/* Search Bar */}
                    <form onSubmit={handleSearch} className="w-full lg:w-auto">
                        <div className="relative group">
                            <input
                                type="text"
                                value={searchInput}
                                onChange={(e) => setSearchInput(e.target.value)}
                                placeholder="Buscar gen, ID o función..."
                                className="w-full lg:w-[400px] pl-14 pr-4 py-4 bg-slate-900/80 border-2 border-slate-700/50 rounded-2xl text-slate-200 placeholder-slate-500 text-sm focus:outline-none focus:border-cyan-500/50 focus:ring-4 focus:ring-cyan-500/10 transition-all shadow-inner"
                            />
                            <svg className="w-6 h-6 text-slate-500 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-cyan-400 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                            </svg>
                            <button type="submit" className="absolute right-2.5 top-1/2 -translate-y-1/2 px-6 py-2 bg-gradient-to-r from-cyan-600 to-blue-600 text-white rounded-xl text-xs font-black uppercase tracking-widest hover:brightness-110 shadow-lg transition-all active:scale-95">
                                BUSCAR
                            </button>
                        </div>
                    </form>
                </div>

                {/* Loading State */}
                {loading ? (
                    <div className="flex flex-col items-center justify-center py-40">
                        <div className="relative w-32 h-32 mb-8">
                            <div className="absolute inset-0 border-[6px] border-slate-800 rounded-full"></div>
                            <div className="absolute inset-0 border-[6px] border-transparent border-t-cyan-500 rounded-full animate-spin"></div>
                            <div className="absolute inset-4 border-[6px] border-transparent border-t-indigo-500 rounded-full animate-spin" style={{ animationDirection: 'reverse', animationDuration: '1.5s' }}></div>
                        </div>
                        <p className="text-slate-400 text-sm font-black uppercase tracking-[0.2em] animate-pulse">Sincronizando Proteoma...</p>
                    </div>
                ) : (
                    <div className="space-y-12">
                        {/* Protein Grid */}
                        <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 2xl:grid-cols-4 gap-6">
                            {proteins.map((protein, index) => {
                                const size = getSizeCategory(protein.length)
                                const isSelected = selectedProteinIndex === index

                                return (
                                    <button
                                        key={`${protein.protein_id}-${index}`}
                                        onClick={() => selectProtein(index)}
                                        className={`group relative text-left bg-slate-800/30 backdrop-blur-sm rounded-3xl border-2 transition-all duration-500 overflow-hidden hover:-translate-y-2 ${
                                            isSelected
                                                ? `border-cyan-500/50 bg-slate-800/80 shadow-2xl shadow-cyan-500/20 ring-4 ${size.ring}`
                                                : 'border-slate-700/30 hover:border-slate-600/50 hover:bg-slate-800/50'
                                        }`}
                                    >
                                        <div className="p-6 space-y-4 relative z-10">
                                            {/* Header */}
                                            <div className="flex items-start justify-between gap-3">
                                                <div className="min-w-0 flex-1">
                                                    <h4 className="font-mono text-lg font-black text-slate-100 truncate group-hover:text-cyan-400 transition-colors">
                                                        {protein.protein_id || protein.locus_tag}
                                                    </h4>
                                                    <div className="flex flex-wrap gap-2 mt-2">
                                                        {protein.gene_name && (
                                                            <span className="px-2 py-0.5 bg-cyan-500/10 text-cyan-400 text-[10px] font-black uppercase tracking-tighter rounded-md border border-cyan-500/20">
                                                                GEN: {protein.gene_name}
                                                            </span>
                                                        )}
                                                        <span className={`px-2 py-0.5 bg-gradient-to-r ${size.color} text-white text-[10px] font-black uppercase tracking-tighter rounded-md shadow-md`}>
                                                            {size.label}
                                                        </span>
                                                    </div>
                                                </div>
                                            </div>

                                            {/* Product Description */}
                                            <p className="text-xs text-slate-400 font-medium line-clamp-2 leading-relaxed h-8">
                                                {protein.product || 'Proteína hipotética'}
                                            </p>

                                            {/* Visual Bar */}
                                            <div className="w-full h-1.5 bg-slate-900/50 rounded-full overflow-hidden">
                                                <div 
                                                    className={`h-full bg-gradient-to-r ${size.color} transition-all duration-1000 group-hover:w-full`}
                                                    style={{ width: `${Math.min((protein.length / 1000) * 100, 100)}%` }}
                                                />
                                            </div>

                                            {/* Stats Grid */}
                                            <div className="grid grid-cols-2 gap-4 pt-4 border-t border-slate-700/30">
                                                <div className="space-y-1">
                                                    <p className="text-[10px] text-slate-500 font-bold uppercase tracking-widest">Longitud</p>
                                                    <p className="font-mono text-sm text-slate-200">
                                                        <span className="font-black text-cyan-400">{protein.length}</span> <span className="text-[10px] text-slate-500">aa</span>
                                                    </p>
                                                </div>
                                                <div className="space-y-1">
                                                    <p className="text-[10px] text-slate-500 font-bold uppercase tracking-widest">Peso Mol.</p>
                                                    <p className="font-mono text-sm text-slate-200">
                                                        <span className="font-black text-indigo-400">{(protein.molecular_weight_approx / 1000).toFixed(1)}</span> <span className="text-[10px] text-slate-500">kDa</span>
                                                    </p>
                                                </div>
                                            </div>
                                        </div>

                                        {/* Background Decor */}
                                        <div className={`absolute -right-4 -bottom-4 w-24 h-24 bg-gradient-to-br ${size.color} blur-3xl opacity-0 group-hover:opacity-20 transition-opacity duration-700`}></div>
                                    </button>
                                )
                            })}
                        </div>

                        {/* Load More Button */}
                        {page < totalPages && (
                            <div className="flex justify-center py-10">
                                <button
                                    onClick={handleLoadMore}
                                    disabled={loadingMore}
                                    className="group relative px-12 py-5 bg-slate-800 border-2 border-slate-700 rounded-2xl overflow-hidden transition-all hover:border-cyan-500/50 active:scale-95 disabled:opacity-50"
                                >
                                    <div className={`absolute inset-0 bg-gradient-to-r from-cyan-600/20 to-indigo-600/20 translate-y-full group-hover:translate-y-0 transition-transform duration-500`}></div>
                                    <div className="relative flex items-center gap-4">
                                        {loadingMore ? (
                                            <>
                                                <div className="w-5 h-5 border-2 border-slate-400 border-t-cyan-400 rounded-full animate-spin"></div>
                                                <span className="text-sm font-black text-slate-400 uppercase tracking-widest">Procesando...</span>
                                            </>
                                        ) : (
                                            <>
                                                <svg className="w-5 h-5 text-cyan-400 group-hover:rotate-180 transition-transform duration-700" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M19 9l-7 7-7-7" />
                                                </svg>
                                                <span className="text-sm font-black text-slate-200 uppercase tracking-[0.2em]">Cargar más proteínas</span>
                                            </>
                                        )}
                                    </div>
                                </button>
                            </div>
                        )}

                        {/* Empty State */}
                        {proteins.length === 0 && !loading && (
                            <div className="flex flex-col items-center justify-center py-40 bg-slate-800/20 rounded-3xl border-2 border-dashed border-slate-700/50">
                                <div className="text-8xl mb-6 grayscale opacity-30">🧬</div>
                                <p className="text-slate-200 text-xl font-black uppercase tracking-widest mb-2">Sin Resultados</p>
                                <p className="text-slate-500 text-sm font-medium">No se encontraron secuencias con los criterios especificados</p>
                                <button 
                                    onClick={() => { setSearchInput(''); setSearch(''); }}
                                    className="mt-8 text-cyan-400 text-xs font-black uppercase tracking-widest hover:underline"
                                >
                                    Limpiar Filtros
                                </button>
                            </div>
                        )}
                    </div>
                )}
            </div>

            {/* Protein Detail Modal */}
            {selectedProtein && (
                <ProteinStructureViewer
                    protein={selectedProtein}
                    onClose={() => setSelectedProteinIndex(null)}
                    onNavigate={navigateProtein}
                    currentIndex={selectedProteinIndex}
                    totalProteins={proteins.length}
                    loading={detailLoading}
                />
            )}
        </div>
    )
}
