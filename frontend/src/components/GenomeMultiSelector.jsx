/**
 * GenomeMultiSelector - Componente para selección múltiple de genomas
 * Unificado con la paleta "Clean Laboratory"
 */
import { useState, useEffect } from 'react'
import { api } from '../services/api'
import toast from 'react-hot-toast'

export default function GenomeMultiSelector({ 
  downloadedGenomes, 
  selectedGenomes, 
  onSelectionChange,
  onRefresh 
}) {
  const [searchQuery, setSearchQuery] = useState('')
  const [searchResults, setSearchResults] = useState([])
  const [isSearching, setIsSearching] = useState(false)
  const [isDownloadingBatch, setIsDownloadingBatch] = useState(false)
  const [downloadProgress, setDownloadProgress] = useState({})
  const [searchLimit, setSearchLimit] = useState(50)
  const [showSuggestions, setShowSuggestions] = useState(true)

  const [relatedStrains, setRelatedStrains] = useState([])

  useEffect(() => {
    loadInitialStrains()
  }, [])

  const loadInitialStrains = async () => {
    try {
      const popular = await api.getPopularGenomes()
      if (popular && popular.length > 0) {
        setRelatedStrains(popular.slice(0, 12))
      } else {
        // Fallback to defaults if API fails or returns empty
        const defaultStrains = [
          { accession: "GCF_000005845.2", organism_name: "E. coli K-12 MG1655", category: "laboratory" },
          { accession: "GCF_000008865.2", organism_name: "E. coli O157:H7 Sakai", category: "pathogenic" },
          { accession: "GCF_000009565.2", organism_name: "E. coli BL21(DE3)", category: "industrial" },
          { accession: "GCF_000019425.1", organism_name: "E. coli CFT073", category: "pathogenic" },
          { accession: "GCF_000007445.1", organism_name: "E. coli W3110", category: "laboratory" },
          { accession: "GCF_000750555.1", organism_name: "E. coli Nissle 1917", category: "probiotic" },
        ]
        setRelatedStrains(defaultStrains)
      }
    } catch (e) {
      console.error("Error loading popular genomes:", e)
    }
  }

  const shuffleStrains = () => {
    const shuffled = [...relatedStrains].sort(() => Math.random() - 0.5)
    setRelatedStrains(shuffled)
  }

  const searchSuggestions = [
    { term: "Escherichia coli", desc: "E. coli" },
    { term: "Salmonella enterica", desc: "Salmonella" },
    { term: "Bacillus subtilis", desc: "Bacillus" },
    { term: "Staphylococcus aureus", desc: "Staph" },
    { term: "Klebsiella pneumoniae", desc: "Klebsiella" },
    { term: "Pseudomonas aeruginosa", desc: "Pseudomonas" },
  ]

  const searchGenomes = async (customQuery = null) => {
    const query = customQuery || searchQuery
    if (!query.trim() || query.length < 3) return
    setIsSearching(true)
    setShowSuggestions(false)
    if (customQuery) setSearchQuery(customQuery)
    try {
      const results = await api.searchGenomes(query, searchLimit)
      setSearchResults(results.results || [])
    } catch (error) {
      setSearchResults([])
    } finally { setIsSearching(false) }
  }

  const toggleSelection = (accession) => {
    const isCurrentlySelected = selectedGenomes.includes(accession)
    const newSelection = isCurrentlySelected
      ? selectedGenomes.filter(a => a !== accession)
      : [...selectedGenomes, accession]
    onSelectionChange(newSelection)
  }

  const isDownloaded = (accession) => downloadedGenomes.some(g => g.accession === accession)
  const isSelected = (accession) => selectedGenomes.includes(accession)

  return (
    <div className="space-y-10 animate-in fade-in duration-700">
      {/* Search Header */}
      <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-48 h-48 bg-blue-500/5 blur-3xl -mr-24 -mt-24"></div>
        <div className="relative z-10 space-y-6">
          <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">Repositorio NCBI</h4>
          <div className="flex gap-4">
            <div className="flex-1 relative group">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={(e) => e.key === 'Enter' && searchGenomes()}
                placeholder="Accession o Nombre de Cepa..."
                className="w-full pl-12 pr-6 py-4 bg-slate-50 border-2 border-slate-100 rounded-2xl text-slate-900 text-sm focus:outline-none focus:border-blue-500/50 transition-all font-medium"
              />
              <svg className="w-5 h-5 text-slate-400 absolute left-4 top-1/2 -translate-y-1/2 group-focus-within:text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeWidth={2.5} /></svg>
            </div>
            <button onClick={() => searchGenomes()} className="px-10 bg-blue-600 text-white rounded-2xl font-black uppercase tracking-widest text-[10px] hover:bg-blue-500 shadow-xl transition-all">Buscar</button>
          </div>
          
          <div className="flex flex-wrap gap-2">
            {searchSuggestions.map(s => (
              <button key={s.term} onClick={() => searchGenomes(s.term)} className="px-4 py-1.5 bg-slate-100 hover:bg-blue-50 text-[10px] font-bold text-slate-500 hover:text-blue-600 rounded-full transition-all border border-transparent hover:border-blue-100 uppercase">
                {s.desc}
              </button>
            ))}
          </div>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-10">
        {/* Available Selection */}
        <div className="space-y-6">
          <div className="flex items-center justify-between px-4">
            <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em]">
              {searchResults.length > 0 ? `Resultados NCBI (${searchResults.length})` : 'Genomas Sugeridos'}
            </h4>
            <div className="flex items-center gap-4">
              {searchResults.length > 0 && (
                <button 
                  onClick={() => setSearchResults([])} 
                  className="text-[9px] font-black text-rose-500 uppercase tracking-widest hover:underline"
                >
                  Limpiar Búsqueda
                </button>
              )}
              <button onClick={shuffleStrains} className="text-[9px] font-black text-blue-600 uppercase tracking-widest hover:underline">Refrescar</button>
            </div>
          </div>
          <div className="grid grid-cols-1 gap-3 max-h-[600px] overflow-y-auto pr-2 custom-scrollbar pb-10">
            {isSearching ? (
              <div className="py-20 text-center animate-pulse">
                <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.5em]">Consultando NCBI...</p>
              </div>
            ) : (searchResults.length > 0 ? searchResults : relatedStrains).map(s => {
              const selected = isSelected(s.accession)
              return (
                <button
                  key={s.accession}
                  onClick={() => toggleSelection(s.accession)}
                  className={`w-full flex items-center justify-between p-5 rounded-[2rem] border-2 transition-all duration-500 group ${
                    selected ? 'bg-blue-600 border-blue-500 text-white shadow-xl shadow-blue-200' : 'bg-white border-slate-100 hover:border-blue-200 text-slate-900 shadow-sm'
                  }`}
                >
                  <div className="text-left flex-1 min-w-0 pr-4">
                    <p className={`font-mono text-[9px] font-black mb-1 ${selected ? 'text-blue-100' : 'text-blue-600'}`}>{s.accession}</p>
                    <p className="text-[11px] font-black uppercase truncate tracking-tight">{s.organism_name || s.organism}</p>
                    {s.strain && <p className={`text-[9px] font-bold ${selected ? 'text-blue-200' : 'text-slate-400'}`}>CEPA: {s.strain}</p>}
                  </div>
                  <div className={`w-8 h-8 rounded-xl flex items-center justify-center transition-all flex-shrink-0 ${selected ? 'bg-white/20' : 'bg-slate-50 text-slate-300 group-hover:bg-blue-50 group-hover:text-blue-500'}`}>
                    {selected ? <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M5 13l4 4L19 7" strokeWidth={3} strokeLinecap="round" strokeLinejoin="round"/></svg> : <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M12 4v16m8-8H4" strokeWidth={3} strokeLinecap="round" strokeLinejoin="round"/></svg>}
                  </div>
                </button>
              )
            })}
          </div>
        </div>

        {/* Workspace Cohort */}
        <div className="space-y-6">
          <div className="flex items-center justify-between px-4">
            <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em]">Entorno de Trabajo</h4>
            <span className="text-[9px] font-black text-white bg-slate-900 px-3 py-1 rounded-full uppercase">{selectedGenomes.length} Genomas</span>
          </div>
          <div className="bg-slate-50/50 rounded-[2.5rem] border-2 border-dashed border-slate-200 p-8 min-h-[400px]">
            {selectedGenomes.length === 0 ? (
              <div className="flex flex-col items-center justify-center h-full text-center space-y-4 py-20 opacity-40">
                <p className="text-[10px] font-black uppercase tracking-[0.2em]">Esperando Selección...</p>
              </div>
            ) : (
              <div className="space-y-3">
                {selectedGenomes.map(acc => (
                  <div key={acc} className="flex items-center justify-between p-4 bg-white rounded-2xl border border-slate-100 shadow-sm">
                    <div className="flex items-center gap-4">
                      <div className="w-1.5 h-1.5 bg-blue-500 rounded-full shadow-[0_0_8px_rgba(59,130,246,0.5)]"></div>
                      <span className="font-mono text-[10px] font-black text-slate-700 tracking-widest">{acc}</span>
                    </div>
                    <button onClick={() => toggleSelection(acc)} className="p-2 text-slate-300 hover:text-rose-500 transition-all"><svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M6 18L18 6M6 6l12 12" strokeWidth={2.5}/></svg></button>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}