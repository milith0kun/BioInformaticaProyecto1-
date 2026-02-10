/**
 * ProteinViewer Component — Zen Laboratory Edition
 * High-performance protein browser with traditional pagination and refined aesthetics
 */
import { useState, useEffect, useCallback } from 'react'
import api from '../services/api'
import ProteinStructureViewer from './ProteinStructureViewer'

export default function ProteinViewer() {
    const [proteins, setProteins] = useState([])
    const [loading, setLoading] = useState(false)
    const [total, setTotal] = useState(0)
    const [page, setPage] = useState(1)
    const [totalPages, setTotalPages] = useState(1)
    const [search, setSearch] = useState('')
    const [searchInput, setSearchInput] = useState('')
    const [selectedProteinIndex, setSelectedProteinIndex] = useState(null)
    const [detailLoading, setDetailLoading] = useState(false)

    const pageSize = 24 // Optimized grid layout

    // Initial load and page/search changes
    useEffect(() => {
        loadProteins(page)
    }, [page, search])

    const loadProteins = async (pageNum) => {
        setLoading(true)
        try {
            const data = await api.getProteins(pageNum, pageSize, search)
            setProteins(data.proteins || [])
            setTotal(data.total || 0)
            setTotalPages(data.total_pages || 1)
        } catch (e) {
            console.error('Error loading proteins:', e)
        } finally {
            setLoading(false)
        }
    }

    const handleSearch = (e) => {
        e.preventDefault()
        setSearch(searchInput)
        setPage(1)
        setSelectedProteinIndex(null)
    }

    const selectProtein = async (index) => {
        setDetailLoading(true)
        setSelectedProteinIndex(index)
        try {
            const protein = proteins[index]
            const detail = await api.getProteinDetail(protein.protein_id || protein.locus_tag)
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
        if (length < 100) return { label: 'Micro', color: 'from-cyan-500 to-blue-500', ring: 'ring-cyan-500/20' }
        if (length < 300) return { label: 'Pequeña', color: 'from-blue-500 to-indigo-500', ring: 'ring-blue-500/20' }
        if (length < 600) return { label: 'Media', color: 'from-indigo-500 to-violet-500', ring: 'ring-indigo-500/20' }
        return { label: 'Grande', color: 'from-violet-500 to-fuchsia-600', ring: 'ring-violet-500/20' }
    }

    const selectedProtein = selectedProteinIndex !== null ? proteins[selectedProteinIndex] : null

    return (
        <div className="min-h-screen bg-[#f1f5f9]">
            <div className="max-w-7xl mx-auto p-8 lg:p-12 space-y-12">
                
                {/* Dashboard Header */}
                <div className="flex flex-col lg:flex-row items-start lg:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border border-slate-200 backdrop-blur-2xl shadow-sm relative overflow-hidden">
                    <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
                    
                    <div className="space-y-4 relative z-10">
                        <h2 className="text-4xl lg:text-6xl font-black tracking-tighter text-slate-900 uppercase italic">
                            PROTEOMA <span className="text-blue-600">LAB</span>
                        </h2>
                        <div className="flex items-center gap-4">
                            <div className="px-4 py-1.5 bg-blue-50 rounded-full border border-blue-100">
                                <p className="text-blue-600 text-[10px] font-black uppercase tracking-[0.2em]">
                                    {total.toLocaleString()} SECUENCIAS
                                </p>
                            </div>
                            <div className="w-1.5 h-1.5 rounded-full bg-emerald-500 shadow-[0_0_8px_rgba(16,185,129,0.5)]"></div>
                            <span className="text-[10px] font-bold text-slate-400 uppercase tracking-widest">Sincronizado</span>
                        </div>
                    </div>

                    {/* Search Engine */}
                    <form onSubmit={handleSearch} className="w-full lg:w-auto relative z-10">
                        <div className="relative group">
                            <input
                                type="text"
                                value={searchInput}
                                onChange={(e) => setSearchInput(e.target.value)}
                                placeholder="Filtrar por Gen o ID..."
                                className="w-full lg:w-[450px] pl-16 pr-6 py-5 bg-slate-50 border-2 border-slate-100 rounded-3xl text-slate-900 placeholder-slate-400 text-sm focus:outline-none focus:border-blue-500/50 focus:ring-8 focus:ring-blue-500/5 transition-all shadow-inner font-medium"
                            />
                            <svg className="w-6 h-6 text-slate-400 absolute left-6 top-1/2 -translate-y-1/2 group-focus-within:text-blue-600 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                            </svg>
                            <button type="submit" className="absolute right-3 top-1/2 -translate-y-1/2 px-8 py-2.5 bg-blue-600 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-700 shadow-xl transition-all active:scale-95">
                                BUSCAR
                            </button>
                        </div>
                    </form>
                </div>

                {/* Content Grid */}
                {loading ? (
                    <div className="flex flex-col items-center justify-center py-48">
                        <div className="w-24 h-24 relative mb-8">
                            <div className="absolute inset-0 border-4 border-slate-200 rounded-full"></div>
                            <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin"></div>
                        </div>
                        <p className="text-slate-400 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse text-center">Indexando Base de Datos...</p>
                    </div>
                ) : (
                    <div className="space-y-16">
                        <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 2xl:grid-cols-4 gap-8">
                            {proteins.map((protein, index) => {
                                const size = getSizeCategory(protein.length)
                                const isSelected = selectedProteinIndex === index

                                return (
                                    <button
                                        key={`${protein.protein_id}-${index}`}
                                        onClick={() => selectProtein(index)}
                                        className={`group relative text-left bg-white rounded-[2.5rem] border-2 transition-all duration-700 overflow-hidden hover:-translate-y-3 ${
                                            isSelected
                                                ? `border-blue-500 bg-blue-50/30 shadow-2xl shadow-blue-500/10 ring-8 ring-blue-500/5`
                                                : 'border-white shadow-md hover:shadow-xl hover:border-slate-100'
                                        }`}
                                    >
                                        <div className="p-8 space-y-5 relative z-10">
                                            {/* Header */}
                                            <div className="flex items-start justify-between">
                                                <div className="min-w-0 flex-1">
                                                    <h4 className="font-mono text-xl font-black text-slate-900 truncate group-hover:text-blue-600 transition-colors tracking-tighter uppercase">
                                                        {protein.protein_id || protein.locus_tag}
                                                    </h4>
                                                    <div className="flex flex-wrap gap-2 mt-3">
                                                        {protein.gene_name && (
                                                            <span className="px-2.5 py-1 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-lg border border-blue-100">
                                                                {protein.gene_name}
                                                            </span>
                                                        )}
                                                        <span className={`px-2.5 py-1 bg-gradient-to-r ${size.color} text-white text-[9px] font-black uppercase tracking-widest rounded-lg shadow-lg`}>
                                                            {size.label}
                                                        </span>
                                                    </div>
                                                </div>
                                            </div>

                                            {/* Description */}
                                            <p className="text-[11px] text-slate-500 font-medium line-clamp-2 leading-relaxed min-h-[2.5rem]">
                                                {protein.product || 'Proteína funcional hipotética'}
                                            </p>

                                            {/* Geometric Progress */}
                                            <div className="w-full h-1 bg-slate-100 rounded-full overflow-hidden">
                                                <div 
                                                    className={`h-full bg-gradient-to-r ${size.color} transition-all duration-1000 group-hover:w-full shadow-[0_0_8px_rgba(59,130,246,0.2)]`}
                                                    style={{ width: `${Math.min((protein.length / 1000) * 100, 100)}%` }}
                                                />
                                            </div>

                                            {/* Analytics Footer */}
                                            <div className="flex items-center justify-between pt-6 border-t border-slate-50">
                                                <div className="space-y-1">
                                                    <p className="text-[9px] text-slate-400 font-black uppercase tracking-widest">Residuos</p>
                                                    <p className="font-mono text-sm text-slate-900 font-black">{protein.length}</p>
                                                </div>
                                                <div className="space-y-1 text-right">
                                                    <p className="text-[9px] text-slate-400 font-black uppercase tracking-widest">Masa</p>
                                                    <p className="font-mono text-sm text-slate-900 font-black">{(protein.molecular_weight_approx / 1000).toFixed(1)} <span className="text-[10px] text-slate-400">kDa</span></p>
                                                </div>
                                            </div>
                                        </div>
                                    </button>
                                )
                            })}
                        </div>

                        {/* Traditional Pagination */}
                        <div className="flex flex-col md:flex-row items-center justify-center gap-8 py-12 border-t border-slate-200">
                            <div className="flex items-center gap-2">
                                <button
                                    onClick={() => setPage(p => Math.max(1, p - 1))}
                                    disabled={page === 1}
                                    className="p-4 bg-white hover:bg-slate-50 rounded-2xl border border-slate-200 disabled:opacity-20 transition-all group shadow-sm"
                                >
                                    <svg className="w-5 h-5 text-slate-400 group-hover:text-blue-600 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M15 19l-7-7 7-7" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5}/></svg>
                                </button>
                                
                                <div className="flex gap-2">
                                    {[...Array(Math.min(5, totalPages))].map((_, i) => {
                                        let pageNum = 1
                                        if (totalPages <= 5) pageNum = i + 1
                                        else if (page <= 3) pageNum = i + 1
                                        else if (page >= totalPages - 2) pageNum = totalPages - 4 + i
                                        else pageNum = page - 2 + i

                                        if (pageNum <= 0 || pageNum > totalPages) return null

                                        return (
                                            <button
                                                key={pageNum}
                                                onClick={() => setPage(pageNum)}
                                                className={`w-12 h-12 rounded-2xl font-black text-xs transition-all border ${
                                                    page === pageNum
                                                        ? 'bg-blue-600 border-blue-500 text-white shadow-xl shadow-blue-900/20 scale-110'
                                                        : 'bg-white border-slate-200 text-slate-500 hover:text-blue-600 hover:bg-slate-50 shadow-sm'
                                                }`}
                                            >
                                                {pageNum}
                                            </button>
                                        )
                                    })}
                                </div>

                                <button
                                    onClick={() => setPage(p => Math.min(totalPages, p + 1))}
                                    disabled={page === totalPages}
                                    className="p-4 bg-white hover:bg-slate-50 rounded-2xl border border-slate-200 disabled:opacity-20 transition-all group shadow-sm"
                                >
                                    <svg className="w-5 h-5 text-slate-400 group-hover:text-blue-600 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M9 5l7 7-7 7" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5}/></svg>
                                </button>
                            </div>
                            
                            <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">
                                PÁGINA {page} DE {totalPages}
                            </p>
                        </div>
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