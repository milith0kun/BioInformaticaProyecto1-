import { useState, useEffect, useRef } from 'react'
import toast from 'react-hot-toast'
import { api } from '../services/api'

// Debounce hook
function useDebounce(value, delay) {
    const [debouncedValue, setDebouncedValue] = useState(value)

    useEffect(() => {
        const handler = setTimeout(() => setDebouncedValue(value), delay)
        return () => clearTimeout(handler)
    }, [value, delay])

    return debouncedValue
}

export default function GenomeSelector({ onGenomeActivated, onGenomeDownloaded, mode = 'search' }) {
    const [searchQuery, setSearchQuery] = useState('')
    const [searchResults, setSearchResults] = useState([])
    const [selectedGenome, setSelectedGenome] = useState(null)
    const [downloadStatus, setDownloadStatus] = useState(null)
    const [isSearching, setIsSearching] = useState(false)
    const [isLoadingInfo, setIsLoadingInfo] = useState(false)
    const [isDownloading, setIsDownloading] = useState(false)
    const [showSuggestions, setShowSuggestions] = useState(false)

    const searchInputRef = useRef(null)
    const suggestionsRef = useRef(null)
    const debouncedQuery = useDebounce(searchQuery, 300)

    useEffect(() => {
        if (debouncedQuery.length >= 2) {
            performSearch(debouncedQuery)
        } else {
            setSearchResults([])
        }
    }, [debouncedQuery])

    useEffect(() => {
        let interval
        if (isDownloading && selectedGenome) {
            interval = setInterval(async () => {
                try {
                    const status = await api.getGenomeDownloadStatus(selectedGenome.accession)
                    setDownloadStatus(status)

                    if (status.status === 'completed') {
                        setIsDownloading(false)
                        toast.success(`Genoma descargado`, { id: 'genome-download' })
                        if (onGenomeDownloaded) onGenomeDownloaded()
                    } else if (status.status === 'error') {
                        setIsDownloading(false)
                        toast.error(status.message, { id: 'genome-download' })
                    }
                } catch (error) {
                    console.error('Status error:', error)
                }
            }, 2000)
        }
        return () => clearInterval(interval)
    }, [isDownloading, selectedGenome])

    useEffect(() => {
        const handleClickOutside = (event) => {
            if (suggestionsRef.current && !suggestionsRef.current.contains(event.target) &&
                searchInputRef.current && !searchInputRef.current.contains(event.target)) {
                setShowSuggestions(false)
            }
        }
        document.addEventListener('mousedown', handleClickOutside)
        return () => document.removeEventListener('mousedown', handleClickOutside)
    }, [])

    const isAccessionNumber = (input) => /^GC[AF]_\d+(\.\d+)?$/i.test(input.trim())

    const extractAccessionFromUrl = (input) => {
        const match = input.match(/GC[AF]_\d+(\.\d+)?/i)
        return match ? match[0] : null
    }

    const performSearch = async (query) => {
        if (!query.trim()) return

        const accession = extractAccessionFromUrl(query) || (isAccessionNumber(query) ? query.trim() : null)

        if (accession) {
            setIsLoadingInfo(true)
            try {
                const info = await api.getGenomeInfo(accession)
                setSearchResults([{
                    accession: info.accession,
                    organism_name: info.organism_name,
                    strain: info.strain,
                    assembly_level: info.assembly_level,
                    genome_size_mb: info.genome_size_mb,
                    gc_percent: info.gc_percent,
                    gene_count: info.gene_count,
                    is_reference: info.is_reference
                }])
                setShowSuggestions(true)
            } catch (error) {
                console.log('Genome not found:', accession)
                setSearchResults([])
            } finally {
                setIsLoadingInfo(false)
            }
        } else {
            setIsSearching(true)
            try {
                const result = await api.searchGenomes(query, 15)
                setSearchResults(result.results || [])
                setShowSuggestions(result.results?.length > 0)
            } catch (error) {
                // 404 means no results found, not an error
                console.log('No results for:', query)
                setSearchResults([])
                setShowSuggestions(false)
            } finally {
                setIsSearching(false)
            }
        }
    }

    const handleSelectGenome = async (genome) => {
        setSelectedGenome(genome)
        setShowSuggestions(false)
        setSearchQuery('')

        // Obtener info completa
        setIsLoadingInfo(true)
        try {
            const info = await api.getGenomeInfo(genome.accession)
            setSelectedGenome(info)
        } catch (error) {
            console.error('Error getting genome info:', error)
        } finally {
            setIsLoadingInfo(false)
        }
    }

    const handleDownload = async () => {
        if (!selectedGenome) return

        setIsDownloading(true)
        setDownloadStatus({ status: 'starting', message: 'Iniciando...' })
        toast.loading('Descargando desde NCBI...', { id: 'genome-download' })

        try {
            await api.downloadGenome({
                accession: selectedGenome.accession,
                include_gbff: true,
                include_gff: true,
                include_fasta: true,
                include_protein: false,
                include_cds: false,
                include_rna: false
            })
        } catch (error) {
            toast.error('Error al descargar', { id: 'genome-download' })
            setIsDownloading(false)
        }
    }

    const handleActivateAfterDownload = async () => {
        if (!selectedGenome) return

        try {
            toast.loading('Activando...', { id: 'activate' })
            const result = await api.activateGenome(selectedGenome.accession)
            toast.success(`Genoma listo`, { id: 'activate' })
            if (onGenomeActivated) onGenomeActivated(result)
        } catch (error) {
            toast.error('Error al activar', { id: 'activate' })
        }
    }

    const generateApiUrl = (accession) =>
        `https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/download`

    return (
        <div className="space-y-6">
            {/* Search Bar */}
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
                <div className="relative">
                    <div className="absolute left-4 top-1/2 -translate-y-1/2 text-slate-400">
                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                        </svg>
                    </div>
                    <input
                        ref={searchInputRef}
                        type="text"
                        value={searchQuery}
                        onChange={(e) => setSearchQuery(e.target.value)}
                        onFocus={() => searchResults.length > 0 && setShowSuggestions(true)}
                        placeholder="Escribe: Escherichia coli, Bacillus, GCF_000005845.2..."
                        className="w-full pl-12 pr-12 py-4 text-lg border-2 border-slate-200 rounded-xl focus:border-teal-500 focus:ring-4 focus:ring-teal-50 transition-all"
                    />

                    {(isSearching || isLoadingInfo) && (
                        <div className="absolute right-4 top-1/2 -translate-y-1/2">
                            <svg className="animate-spin h-5 w-5 text-teal-500" viewBox="0 0 24 24">
                                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                            </svg>
                        </div>
                    )}

                    {/* Suggestions Dropdown */}
                    {showSuggestions && searchResults.length > 0 && (
                        <div
                            ref={suggestionsRef}
                            className="absolute z-50 w-full mt-2 bg-white rounded-xl shadow-xl border border-slate-200 max-h-80 overflow-y-auto"
                        >
                            <div className="px-4 py-3 bg-slate-50 border-b border-slate-100 sticky top-0">
                                <span className="text-sm text-slate-600 font-medium">
                                    {searchResults.length} resultados
                                </span>
                            </div>

                            {searchResults.map((genome, index) => (
                                <div
                                    key={genome.accession}
                                    onClick={() => handleSelectGenome(genome)}
                                    className={`p-4 cursor-pointer transition-all hover:bg-teal-50 border-b border-slate-50 last:border-0 ${index === 0 ? 'bg-teal-50/30' : ''
                                        }`}
                                >
                                    <div className="flex justify-between items-start gap-4">
                                        <div className="flex-1 min-w-0">
                                            <div className="flex items-center gap-2 flex-wrap">
                                                <span className="font-mono text-sm text-teal-700 bg-teal-100 px-2 py-0.5 rounded font-medium">
                                                    {genome.accession}
                                                </span>
                                                {genome.is_reference && (
                                                    <span className="text-xs bg-emerald-100 text-emerald-700 px-2 py-0.5 rounded-full">
                                                        Reference
                                                    </span>
                                                )}
                                            </div>
                                            <h4 className="font-medium text-slate-800 mt-1 truncate">{genome.organism_name}</h4>
                                            {genome.strain && (
                                                <p className="text-sm text-slate-500 truncate">Strain: {genome.strain}</p>
                                            )}
                                        </div>
                                        <div className="text-right text-sm flex-shrink-0">
                                            <div className="font-semibold text-teal-600">{genome.genome_size_mb} Mb</div>
                                            <div className="text-xs text-slate-400">{genome.assembly_level}</div>
                                        </div>
                                    </div>
                                </div>
                            ))}
                        </div>
                    )}
                </div>

                {/* Hint Text */}
                {searchQuery.length > 0 && searchQuery.length < 2 && (
                    <p className="mt-2 text-sm text-slate-400">Escribe al menos 2 caracteres...</p>
                )}

                {/* Examples */}
                <div className="mt-4 flex flex-wrap gap-2">
                    <span className="text-sm text-slate-500">Prueba:</span>
                    {['Escherichia coli', 'Bacillus subtilis', 'GCF_000005845.2'].map(ex => (
                        <button
                            key={ex}
                            onClick={() => setSearchQuery(ex)}
                            className="px-3 py-1 bg-slate-100 hover:bg-teal-100 text-slate-600 hover:text-teal-700 rounded-lg text-sm transition-all"
                        >
                            {ex}
                        </button>
                    ))}
                </div>
            </div>

            {/* Selected Genome Card */}
            {selectedGenome && (
                <div className="bg-white rounded-xl shadow-sm border border-teal-200 p-6">
                    <div className="flex items-start justify-between mb-4">
                        <div>
                            <div className="flex items-center gap-3 mb-2">
                                <span className="font-mono text-xl text-teal-700 bg-teal-50 px-3 py-1 rounded-lg font-bold">
                                    {selectedGenome.accession}
                                </span>
                                {selectedGenome.is_reference && (
                                    <span className="bg-emerald-50 text-emerald-700 px-2 py-1 rounded-full text-xs font-medium">
                                        Reference
                                    </span>
                                )}
                            </div>
                            <h3 className="text-xl font-bold text-slate-800">{selectedGenome.organism_name}</h3>
                            {selectedGenome.strain && (
                                <p className="text-slate-600">Strain: {selectedGenome.strain}</p>
                            )}
                        </div>
                        <button
                            onClick={() => setSelectedGenome(null)}
                            className="p-2 text-slate-400 hover:text-slate-600 hover:bg-slate-100 rounded-lg transition-all"
                        >
                            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                            </svg>
                        </button>
                    </div>

                    {/* Stats */}
                    <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 mb-4">
                        {[
                            { label: 'Tamaño', value: `${selectedGenome.genome_size_mb} Mb`, color: 'teal' },
                            { label: 'GC %', value: `${selectedGenome.gc_percent?.toFixed(1)}%`, color: 'emerald' },
                            { label: 'Genes', value: selectedGenome.gene_count?.toLocaleString() || 'N/A', color: 'slate' },
                            { label: 'Nivel', value: selectedGenome.assembly_level, color: 'slate' },
                        ].map((stat, i) => (
                            <div key={i} className="bg-slate-50 rounded-lg p-3 text-center">
                                <div className="text-xs text-slate-500 uppercase">{stat.label}</div>
                                <div className={`text-lg font-bold ${stat.color === 'teal' ? 'text-teal-700' : stat.color === 'emerald' ? 'text-emerald-700' : 'text-slate-700'}`}>
                                    {stat.value}
                                </div>
                            </div>
                        ))}
                    </div>

                    {/* Links */}
                    <div className="flex flex-wrap gap-4 text-sm mb-4">
                        <a
                            href={`https://www.ncbi.nlm.nih.gov/datasets/genome/${selectedGenome.accession}/`}
                            target="_blank"
                            rel="noopener noreferrer"
                            className="text-teal-600 hover:text-teal-800 flex items-center gap-1"
                        >
                            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
                            </svg>
                            Ver en NCBI
                        </a>
                        <button
                            onClick={() => {
                                navigator.clipboard.writeText(generateApiUrl(selectedGenome.accession))
                                toast.success('URL copiada')
                            }}
                            className="text-slate-500 hover:text-slate-700 flex items-center gap-1"
                        >
                            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 16H6a2 2 0 01-2-2V6a2 2 0 012-2h8a2 2 0 012 2v2m-6 12h8a2 2 0 002-2v-8a2 2 0 00-2-2h-8a2 2 0 00-2 2v8a2 2 0 002 2z" />
                            </svg>
                            Copiar API URL
                        </button>
                    </div>

                    {/* Download Progress */}
                    {isDownloading && downloadStatus && (
                        <div className="bg-teal-50 rounded-lg p-4 mb-4">
                            <div className="flex items-center gap-3">
                                <svg className="w-5 h-5 text-teal-600 animate-bounce" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M9 19l3 3m0 0l3-3m-3 3V10" />
                                </svg>
                                <div>
                                    <div className="font-medium text-teal-800">Descargando...</div>
                                    <div className="text-sm text-teal-600">{downloadStatus.message}</div>
                                </div>
                            </div>
                        </div>
                    )}

                    {/* Actions */}
                    <div className="flex flex-col sm:flex-row gap-3">
                        {downloadStatus?.status === 'completed' ? (
                            <button
                                onClick={handleActivateAfterDownload}
                                className="flex-1 py-3 bg-gradient-to-r from-emerald-500 to-teal-500 hover:from-emerald-600 hover:to-teal-600 text-white rounded-xl font-semibold transition-all shadow-lg flex items-center justify-center gap-2"
                            >
                                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                                </svg>
                                Usar este genoma
                            </button>
                        ) : (
                            <button
                                onClick={handleDownload}
                                disabled={isDownloading}
                                className="flex-1 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700 disabled:from-slate-400 disabled:to-slate-500 text-white rounded-xl font-semibold transition-all shadow-lg flex items-center justify-center gap-2"
                            >
                                {isDownloading ? (
                                    <>
                                        <svg className="animate-spin h-5 w-5" viewBox="0 0 24 24">
                                            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                                            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                                        </svg>
                                        Descargando...
                                    </>
                                ) : (
                                    <>
                                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                                        </svg>
                                        Descargar de NCBI
                                    </>
                                )}
                            </button>
                        )}
                    </div>
                </div>
            )}

            {/* Help Box */}
            {!selectedGenome && (
                <div className="bg-teal-50 border border-teal-100 rounded-xl p-5">
                    <h4 className="font-semibold text-teal-800 mb-3">¿Cómo buscar?</h4>
                    <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 text-sm">
                        <div className="flex items-start gap-2">
                            <span className="w-6 h-6 bg-teal-500 text-white rounded text-xs flex items-center justify-center font-medium">1</span>
                            <div>
                                <div className="font-medium text-teal-800">Por nombre</div>
                                <div className="text-teal-600">Ej: "Escherichia coli"</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-2">
                            <span className="w-6 h-6 bg-teal-500 text-white rounded text-xs flex items-center justify-center font-medium">2</span>
                            <div>
                                <div className="font-medium text-teal-800">Por accession</div>
                                <div className="text-teal-600">Ej: "GCF_000005845.2"</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-2">
                            <span className="w-6 h-6 bg-teal-500 text-white rounded text-xs flex items-center justify-center font-medium">3</span>
                            <div>
                                <div className="font-medium text-teal-800">Por URL</div>
                                <div className="text-teal-600">Pega una URL de NCBI</div>
                            </div>
                        </div>
                    </div>
                </div>
            )}
        </div>
    )
}
