/**
 * ProteinViewer Component — Clean Laboratory Edition
 * High-performance protein browser with size filtering, advanced search and 3D visualization.
 */
import { useState, useEffect, useCallback, useMemo } from 'react'
import api from '../services/api'
import ProteinStructureViewer from './ProteinStructureViewer'

const SIZE_FILTERS = [
  { id: 'all', label: 'Todas', min: 0, max: 10000, color: 'text-slate-500' },
  { id: 'micro', label: 'Micro', min: 0, max: 100, color: 'text-cyan-500' },
  { id: 'small', label: 'Pequeñas', min: 101, max: 300, color: 'text-blue-500' },
  { id: 'medium', label: 'Medias', min: 301, max: 600, color: 'text-indigo-500' },
  { id: 'large', label: 'Grandes', min: 601, max: 10000, color: 'text-violet-600' },
]

export default function ProteinViewer() {
    const [allProteins, setAllProteins] = useState([])
    const [loading, setLoading] = useState(false)
    const [page, setPage] = useState(1)
    const [search, setSearch] = useState('')
    const [searchInput, setSearchInput] = useState('')
    const [sizeFilter, setSizeFilter] = useState('all')
    const [selectedProteinIndex, setSelectedProteinIndex] = useState(null)
    const [detailLoading, setDetailLoading] = useState(false)

    const pageSize = 24

    // Load ALL available proteins once (limited to 5000 by backend)
    useEffect(() => {
        loadData()
    }, [])

    const loadData = async () => {
        setLoading(true)
        try {
            const data = await api.getProteins(1, 5000, '')
            setAllProteins(data.proteins || [])
        } catch (e) {
            console.error('Error loading proteins:', e)
        } finally {
            setLoading(false)
        }
    }

    // Advanced filtering logic
    const filteredProteins = useMemo(() => {
        const currentRange = SIZE_FILTERS.find(f => f.id === sizeFilter)
        const term = search.toLowerCase().trim()

        return allProteins.filter(p => {
            const matchesSearch = term === '' || 
                (p.product?.toLowerCase().includes(term)) ||
                (p.gene_name?.toLowerCase().includes(term)) ||
                (p.protein_id?.toLowerCase().includes(term)) ||
                (p.locus_tag?.toLowerCase().includes(term))
            
            const matchesSize = p.length >= currentRange.min && p.length <= currentRange.max
            
            return matchesSearch && matchesSize
        })
    }, [allProteins, search, sizeFilter])

    const totalPages = Math.max(1, Math.ceil(filteredProteins.length / pageSize))
    const paginatedProteins = useMemo(() => {
        const start = (page - 1) * pageSize
        return filteredProteins.slice(start, start + pageSize)
    }, [filteredProteins, page])

    const handleSearch = (e) => {
        e.preventDefault()
        setSearch(searchInput)
        setPage(1)
        setSelectedProteinIndex(null)
    }

    const selectProtein = async (globalIndex) => {
        const protein = paginatedProteins[globalIndex]
        if (!protein) return
        
        setDetailLoading(true)
        setSelectedProteinIndex(globalIndex)
        
        try {
            const detail = await api.getProteinDetail(protein.protein_id || protein.locus_tag)
            setAllProteins(prev => prev.map(p => 
                (p.protein_id === detail.protein_id || p.locus_tag === detail.locus_tag) ? detail : p
            ))
        } catch (e) {
            console.error('Error loading protein detail:', e)
        } finally {
            setDetailLoading(false)
        }
    }

    const navigateProtein = useCallback((direction) => {
        if (selectedProteinIndex === null) return
        const newIndex = direction === 'next'
            ? Math.min(selectedProteinIndex + 1, paginatedProteins.length - 1)
            : Math.max(selectedProteinIndex - 1, 0)
        
        if (newIndex !== selectedProteinIndex) {
            selectProtein(newIndex)
        }
    }, [selectedProteinIndex, paginatedProteins])

    // FIX: Using correct variable name 'selectedProteinIndex'
    const currentSelectedProtein = selectedProteinIndex !== null ? paginatedProteins[selectedProteinIndex] : null

    return (
        <div className="space-y-10 animate-in fade-in duration-1000">
            
            {/* 1. Dashboard Header & Search */}
            <div className="flex flex-col lg:flex-row items-stretch lg:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
                <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
                
                <div className="space-y-4 relative z-10">
                    <h2 className="text-4xl font-black tracking-tighter text-slate-900 uppercase italic leading-none">
                        Proteoma <span className="text-blue-600">Lab</span>
                    </h2>
                    <div className="flex items-center gap-4">
                        <div className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[10px] font-black uppercase tracking-widest rounded-full border border-blue-100 shadow-sm">
                            {allProteins.length.toLocaleString()} Secuencias Detectadas
                        </div>
                        {search && (
                            <span className="text-[10px] font-bold text-slate-600 uppercase tracking-widest">
                                Filtrando: <span className="text-blue-600">{filteredProteins.length}</span> resultados
                            </span>
                        )}
                    </div>
                </div>

                <form onSubmit={handleSearch} className="relative z-10 w-full lg:w-[500px]">
                    <div className="relative group">
                        <input
                            type="text"
                            value={searchInput}
                            onChange={(e) => setSearchInput(e.target.value)}
                            placeholder="Buscar por Producto, Gen o ID..."
                            className="w-full pl-12 pr-32 py-4 bg-slate-50 border-2 border-slate-100 rounded-2xl text-slate-900 text-[10px] font-black uppercase tracking-widest focus:border-blue-500/50 transition-all placeholder-slate-400 shadow-inner"
                        />
                        <svg className="w-5 h-5 text-slate-400 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-blue-600 transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                        </svg>
                        <button type="submit" className="absolute right-2 top-1/2 -translate-y-1/2 px-6 py-2 bg-slate-900 text-white rounded-xl text-[9px] font-black uppercase tracking-widest hover:bg-blue-600 transition-all active:scale-95">
                            Buscar
                        </button>
                    </div>
                </form>
            </div>

            {/* 2. Size Filters Tabs (Centered) */}
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-2 flex flex-wrap justify-center gap-3 shadow-sm">
                {SIZE_FILTERS.map(f => (
                    <button
                        key={f.id}
                        onClick={() => { setSizeFilter(f.id); setPage(1); }}
                        className={`flex items-center gap-3 px-8 py-3 rounded-2xl text-[10px] font-black uppercase tracking-widest transition-all ${
                            sizeFilter === f.id 
                                ? 'bg-blue-600 text-white shadow-xl shadow-blue-900/20' 
                                : 'text-slate-500 hover:bg-slate-50 hover:text-slate-900'
                        }`}
                    >
                        <span>{f.label}</span>
                        <span className={`ml-2 px-2 py-0.5 rounded-lg text-[8px] ${sizeFilter === f.id ? 'bg-white/20 text-white' : 'bg-slate-100 text-slate-400'}`}>
                            {allProteins.filter(p => p.length >= f.min && p.length <= f.max).length}
                        </span>
                    </button>
                ))}
            </div>

            {/* 3. Content Grid */}
            {loading ? (
                <div className="flex flex-col items-center justify-center py-48">
                    <div className="w-16 h-16 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
                    <p className="text-slate-600 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse text-center">Indexando Proteoma...</p>
                </div>
            ) : filteredProteins.length === 0 ? (
                <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-20 text-center space-y-6">
                    <h3 className="text-xl font-black text-slate-900 uppercase tracking-tighter">Sin Resultados</h3>
                    <p className="text-slate-600 text-xs font-bold uppercase max-w-sm mx-auto leading-relaxed">No se encontraron proteínas que coincidan con los criterios de búsqueda y tamaño seleccionados.</p>
                    <button onClick={() => { setSearch(''); setSearchInput(''); setSizeFilter('all'); }} className="px-8 py-3 bg-blue-50 text-blue-600 rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-100 transition-all shadow-sm">Limpiar Filtros</button>
                </div>
            ) : (
                <div className="space-y-12">
                    <div className="grid grid-cols-1 sm:grid-cols-2 xl:grid-cols-3 2xl:grid-cols-4 gap-6">
                        {paginatedProteins.map((protein, index) => (
                            <button
                                key={`${protein.protein_id}-${index}`}
                                onClick={() => selectProtein(index)}
                                className="group relative text-left bg-white rounded-[2.5rem] border-2 border-slate-100 hover:border-blue-200 transition-all duration-500 overflow-hidden hover:-translate-y-2 hover:shadow-2xl hover:shadow-blue-900/5 active:scale-[0.98]"
                            >
                                <div className="p-8 space-y-5 relative z-10">
                                    <div className="flex items-start justify-between gap-4">
                                        <div className="min-w-0 flex-1">
                                            <h4 className="font-mono text-lg font-black text-slate-900 truncate group-hover:text-blue-600 transition-colors tracking-tighter uppercase leading-none">
                                                {protein.gene_name || protein.locus_tag}
                                            </h4>
                                            <p className="text-[10px] font-black text-blue-600 uppercase tracking-widest mt-2">{protein.protein_id}</p>
                                        </div>
                                    </div>

                                    <p className="text-[11px] text-slate-600 font-bold uppercase line-clamp-2 leading-relaxed min-h-[2.5rem]">
                                        {protein.product || 'Proteína funcional hipotética'}
                                    </p>

                                    <div className="w-full h-1 bg-slate-50 rounded-full overflow-hidden">
                                        <div 
                                            className={`h-full bg-blue-600 transition-all duration-1000 shadow-[0_0_8px_rgba(37,99,235,0.3)]`}
                                            style={{ width: `${Math.min((protein.length / 1000) * 100, 100)}%` }}
                                        />
                                    </div>

                                    <div className="flex items-center justify-between pt-6 border-t border-slate-50">
                                        <div className="space-y-1">
                                            <p className="text-[8px] text-slate-500 font-black uppercase tracking-widest">Residuos</p>
                                            <p className="font-mono text-xs text-slate-900 font-black">{protein.length} aa</p>
                                        </div>
                                        <div className="space-y-1 text-right">
                                            <p className="text-[8px] text-slate-500 font-black uppercase tracking-widest">Masa</p>
                                            <p className="font-mono text-xs text-slate-900 font-black">{((protein.molecular_weight_approx || 0) / 1000).toFixed(1)} <span className="text-[10px] text-slate-400 uppercase font-bold">kDa</span></p>
                                        </div>
                                    </div>
                                </div>
                            </button>
                        ))}
                    </div>

                    <div className="flex flex-col md:flex-row items-center justify-between gap-8 py-10 px-10 bg-white rounded-[2.5rem] border-2 border-slate-100 shadow-sm">
                        <div className="flex items-center gap-4">
                            <button
                                onClick={() => setPage(p => Math.max(1, p - 1))}
                                disabled={page === 1}
                                className="p-4 bg-slate-50 hover:bg-slate-100 rounded-2xl border border-slate-100 disabled:opacity-20 transition-all group active:scale-95 shadow-sm"
                            >
                                <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M15 19l-7-7 7-7" strokeLinecap="round" strokeLinejoin="round" strokeWidth={3}/></svg>
                            </button>
                            
                            <div className="px-8 py-3 bg-white rounded-2xl border-2 border-slate-100 text-[10px] font-black uppercase tracking-widest shadow-inner text-slate-900">
                                Página <span className="text-blue-700">{page}</span> de {totalPages}
                            </div>

                            <button
                                onClick={() => setPage(p => Math.min(totalPages, p + 1))}
                                disabled={page === totalPages}
                                className="p-4 bg-slate-50 hover:bg-slate-100 rounded-2xl border border-slate-100 disabled:opacity-20 transition-all group active:scale-95 shadow-sm"
                            >
                                <svg className="w-5 h-5 text-slate-600" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M9 5l7 7-7 7" strokeLinecap="round" strokeLinejoin="round" strokeWidth={3}/></svg>
                            </button>
                        </div>
                        
                        <div className="flex gap-10">
                            <div className="text-center">
                                <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-1">Mostrando</p>
                                <p className="text-xl font-black text-slate-900">{paginatedProteins.length}</p>
                            </div>
                            <div className="text-center border-l border-slate-100 pl-10">
                                <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-1">Total Filtrado</p>
                                <p className="text-xl font-black text-blue-600">{filteredProteins.length}</p>
                            </div>
                        </div>
                    </div>
                </div>
            )}

            {currentSelectedProtein && (
                <ProteinStructureViewer
                    protein={currentSelectedProtein}
                    onClose={() => setSelectedProteinIndex(null)}
                    onNavigate={navigateProtein}
                    currentIndex={selectedProteinIndex}
                    totalProteins={paginatedProteins.length}
                    loading={detailLoading}
                />
            )}
        </div>
    )
}
