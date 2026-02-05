/**
 * GenomeMultiSelector - Componente para selección múltiple de genomas
 * Permite seleccionar 1 o más genomas, descargarlos en grupo y analizarlos
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

  // Lista extendida de cepas de E. coli
  const allEcoliStrains = [
    { accession: "GCF_000005845.2", organism: "E. coli K-12 MG1655", strain: "K-12", category: "laboratory" },
    { accession: "GCF_000008865.2", organism: "E. coli O157:H7 Sakai", strain: "O157:H7", category: "pathogenic" },
    { accession: "GCF_000009565.2", organism: "E. coli BL21(DE3)", strain: "BL21", category: "industrial" },
    { accession: "GCF_000019425.1", organism: "E. coli CFT073", strain: "CFT073", category: "pathogenic" },
    { accession: "GCF_000007445.1", organism: "E. coli W3110", strain: "W3110", category: "laboratory" },
    { accession: "GCF_000750555.1", organism: "E. coli Nissle 1917", strain: "Nissle", category: "probiotic" },
    { accession: "GCF_000013265.1", organism: "E. coli O127:H6 E2348/69", strain: "E2348/69", category: "pathogenic" },
    { accession: "GCF_000017985.1", organism: "E. coli SE11", strain: "SE11", category: "laboratory" },
    { accession: "GCF_000026245.1", organism: "E. coli ETEC H10407", strain: "H10407", category: "pathogenic" },
    { accession: "GCF_000010385.1", organism: "E. coli O103:H2 12009", strain: "O103:H2", category: "pathogenic" },
    { accession: "GCF_000025165.1", organism: "E. coli DH1", strain: "DH1", category: "laboratory" },
    { accession: "GCF_000167715.1", organism: "E. coli C", strain: "C-ATCC 8739", category: "laboratory" },
  ]

  // Seleccionar aleatoriamente 6 cepas para mostrar
  const [relatedStrains, setRelatedStrains] = useState(() => {
    const shuffled = [...allEcoliStrains].sort(() => Math.random() - 0.5)
    return shuffled.slice(0, 6)
  })

  // Función para cambiar las cepas mostradas
  const shuffleStrains = () => {
    const shuffled = [...allEcoliStrains].sort(() => Math.random() - 0.5)
    setRelatedStrains(shuffled.slice(0, 6))
  }

  // Sugerencias de búsqueda populares
  const searchSuggestions = [
    { term: "Escherichia coli", icon: "🦠", desc: "E. coli" },
    { term: "Salmonella enterica", icon: "🦠", desc: "Salmonella" },
    { term: "Bacillus subtilis", icon: "🦠", desc: "Bacillus" },
    { term: "Staphylococcus aureus", icon: "🦠", desc: "Staphylococcus" },
    { term: "Pseudomonas aeruginosa", icon: "🦠", desc: "Pseudomonas" },
    { term: "Mycobacterium tuberculosis", icon: "🦠", desc: "Tuberculosis" },
    { term: "Streptococcus pneumoniae", icon: "🦠", desc: "Streptococcus" },
    { term: "Klebsiella pneumoniae", icon: "🦠", desc: "Klebsiella" },
  ]

  // Buscar genomas con límite configurable
  const searchGenomes = async (customQuery = null, customLimit = null) => {
    const query = customQuery || searchQuery
    const limit = customLimit || searchLimit
    
    if (!query.trim() || query.length < 3) {
      if (query.length === 0) {
        setSearchResults([])
        setShowSuggestions(true)
      }
      return
    }
    
    setIsSearching(true)
    setShowSuggestions(false)
    
    // Si es una sugerencia, actualizar el input
    if (customQuery) {
      setSearchQuery(customQuery)
    }
    
    try {
      const results = await api.searchGenomes(query, limit)
      setSearchResults(results.results || [])
      if (customQuery && results.results && results.results.length > 0) {
        toast.success(`${results.results.length} genomas encontrados para "${query}"`)
      }
    } catch (error) {
      console.error('Search error:', error)
      setSearchResults([])
      if (customQuery) {
        toast.error('No se encontraron genomas. Intenta con otro término.')
      }
    } finally {
      setIsSearching(false)
    }
  }

  // Búsqueda automática mientras escribe (debounce)
  useEffect(() => {
    if (searchQuery.length === 0) {
      setSearchResults([])
      setShowSuggestions(true)
      return
    }

    if (searchQuery.length < 3) {
      return
    }

    const timeoutId = setTimeout(async () => {
      setIsSearching(true)
      setShowSuggestions(false)
      
      try {
        const results = await api.searchGenomes(searchQuery, searchLimit)
        setSearchResults(results.results || [])
      } catch (error) {
        console.error('Search error:', error)
        setSearchResults([])
      } finally {
        setIsSearching(false)
      }
    }, 600) // Espera 600ms después de que el usuario deja de escribir

    return () => clearTimeout(timeoutId)
  }, [searchQuery, searchLimit])

  // Limpiar búsqueda
  const clearSearch = () => {
    setSearchResults([])
    setSearchQuery('')
    setShowSuggestions(true)
  }

  // Extraer número base del accession (sin prefijo GCA/GCF)
  const getBaseAccession = (accession) => {
    // GCF_000005845.2 → 000005845.2, GCA_000005845.2 → 000005845.2
    return accession.replace(/^GC[AF]_/, '')
  }

  // Verificar si un genoma ya está seleccionado (considerando GCA/GCF equivalentes)
  const isEquivalentSelected = (accession) => {
    const baseAcc = getBaseAccession(accession)
    return selectedGenomes.some(selected => getBaseAccession(selected) === baseAcc)
  }

  // Obtener el accession equivalente ya seleccionado
  const getEquivalentSelected = (accession) => {
    const baseAcc = getBaseAccession(accession)
    return selectedGenomes.find(selected => getBaseAccession(selected) === baseAcc)
  }

  // Toggle selección de genoma (funciona para descargados y no descargados)
  const toggleSelection = (accession) => {
    const isCurrentlySelected = selectedGenomes.includes(accession)
    
    // Check for equivalent selection (GCA vs GCF)
    if (!isCurrentlySelected) {
      const equivalent = getEquivalentSelected(accession)
      if (equivalent && equivalent !== accession) {
        toast.error(`Ya tienes seleccionado ${equivalent} (mismo genoma con diferente prefijo)`, { duration: 3000 })
        return
      }
    }
    
    const newSelection = isCurrentlySelected
      ? selectedGenomes.filter(a => a !== accession)
      : [...selectedGenomes, accession]
    
    console.log('Toggle selection:', {
      accession,
      wasSelected: isCurrentlySelected,
      newSelection
    })
    
    onSelectionChange(newSelection)
    
    // Feedback visual
    if (!isCurrentlySelected) {
      toast.success(`${accession} seleccionado`, { duration: 1000 })
    }
  }

  // Seleccionar todas las cepas (de las 12 disponibles)
  const selectAllStrains = () => {
    const allAccessions = allEcoliStrains.map(s => s.accession)
    onSelectionChange(allAccessions)
  }

  // Descargar todos los genomas seleccionados que no estén descargados
  const downloadSelectedGenomes = async () => {
    const toDownload = selectedGenomes.filter(acc => !isDownloaded(acc))
    
    if (toDownload.length === 0) {
      toast.success('Todos los genomas seleccionados ya están descargados')
      return
    }

    setIsDownloadingBatch(true)
    toast.loading(`Descargando ${toDownload.length} genomas...`, { id: 'batch' })

    for (const accession of toDownload) {
      setDownloadProgress(prev => ({ ...prev, [accession]: 'downloading' }))
      
      try {
        await api.downloadGenome({
          accession,
          include_gbff: true,
          include_gff: true,
          include_fasta: true
        })

        // Esperar a que termine
        let completed = false
        let attempts = 0
        while (!completed && attempts < 60) {
          await new Promise(resolve => setTimeout(resolve, 2000))
          try {
            const status = await api.getGenomeDownloadStatus(accession)
            if (status.status === 'completed') {
              completed = true
              setDownloadProgress(prev => ({ ...prev, [accession]: 'completed' }))
            } else if (status.status === 'error') {
              setDownloadProgress(prev => ({ ...prev, [accession]: 'error' }))
              break
            }
          } catch (e) {
            attempts++
          }
          attempts++
        }
      } catch (error) {
        setDownloadProgress(prev => ({ ...prev, [accession]: 'error' }))
      }
    }

    await onRefresh()
    setIsDownloadingBatch(false)
    toast.success(`${toDownload.length} genomas descargados`, { id: 'batch' })
  }

  const isDownloaded = (accession) => downloadedGenomes.some(g => g.accession === accession)
  const isSelected = (accession) => selectedGenomes.includes(accession)
  const isEquivalentAlreadySelected = (accession) => {
    // Check if equivalent (GCA/GCF) is already selected
    if (selectedGenomes.includes(accession)) return false // It's the same, not equivalent
    return isEquivalentSelected(accession)
  }
  const getDownloadStatus = (accession) => downloadProgress[accession] || null

  // Contar cuántos están descargados y cuántos faltan
  const downloadedCount = selectedGenomes.filter(acc => isDownloaded(acc)).length
  const toDownloadCount = selectedGenomes.length - downloadedCount

  // Función para eliminar genoma
  const handleDeleteGenome = async (accession) => {
    try {
      toast.loading(`Eliminando ${accession}...`, { id: 'delete' })
      await api.deleteGenome(accession)
      toast.success(`${accession} eliminado`, { id: 'delete' })
      
      // Actualizar lista de descargados
      if (onRefresh) {
        await onRefresh()
      }
      
      // Quitar de selección si estaba seleccionado
      if (selectedGenomes.includes(accession)) {
        onSelectionChange(selectedGenomes.filter(acc => acc !== accession))
      }
    } catch (error) {
      toast.error(`Error al eliminar: ${error.message}`, { id: 'delete' })
    }
  }

  return (
    <div className="space-y-6">
      {/* Resumen y acciones */}
      <div className="bg-gradient-to-r from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
        <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
          <div className="flex-1">
            <h3 className="font-bold text-xl mb-1">
              {selectedGenomes.length === 0 
                ? '🧬 Selecciona genomas' 
                : `✓ ${selectedGenomes.length} genoma${selectedGenomes.length > 1 ? 's' : ''} seleccionado${selectedGenomes.length > 1 ? 's' : ''}`}
            </h3>
            <p className="text-teal-100 text-sm">
              {selectedGenomes.length === 0 
                ? 'Marca las casillas de los genomas que quieres analizar'
                : selectedGenomes.length === 1 
                  ? `${downloadedCount}/1 descargado • Análisis individual`
                  : `${downloadedCount}/${selectedGenomes.length} descargados • Análisis comparativo`}
            </p>
          </div>
          
          {selectedGenomes.length > 0 && (
            <div className="flex gap-2">
              {toDownloadCount > 0 && (
                <button
                  onClick={downloadSelectedGenomes}
                  disabled={isDownloadingBatch}
                  className="px-4 py-2 bg-white text-teal-700 rounded-lg font-medium hover:bg-teal-50 disabled:opacity-50 flex items-center gap-2"
                >
                  {isDownloadingBatch ? (
                    <>
                      <svg className="animate-spin h-4 w-4" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                      </svg>
                      Descargando...
                    </>
                  ) : (
                    <>
                      <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M9 19l3 3m0 0l3-3m-3 3V10" />
                      </svg>
                      Descargar {toDownloadCount}
                    </>
                  )}
                </button>
              )}
              <button
                onClick={() => onSelectionChange([])}
                className="px-4 py-2 bg-white/20 hover:bg-white/30 rounded-lg font-medium text-sm"
              >
                Limpiar
              </button>
            </div>
          )}
        </div>
      </div>

      {/* Genomas descargados - Gestión */}
      {downloadedGenomes.length > 0 && (
        <div className="bg-white rounded-xl border border-slate-200 p-5">
          <div className="flex items-center justify-between mb-4">
            <h3 className="font-semibold text-slate-800">
              📦 Genomas Descargados ({downloadedGenomes.length})
            </h3>
            <div className="flex items-center gap-4">
              <button
                onClick={() => {
                  // Seleccionar todos los descargados que no tengan equivalente seleccionado
                  downloadedGenomes.forEach(genome => {
                    if (!selectedGenomes.includes(genome.accession) && !isEquivalentAlreadySelected(genome.accession)) {
                      setSelectedGenomes(prev => [...prev, genome.accession])
                    }
                  })
                }}
                className="text-xs text-teal-600 hover:text-teal-700 font-medium"
              >
                Seleccionar todos
              </button>
              <span className="text-xs text-slate-500">
                Haz clic para seleccionar
              </span>
            </div>
          </div>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3 max-h-96 overflow-y-auto">
            {downloadedGenomes.map(genome => {
              const selected = isSelected(genome.accession)
              const hasEquivalent = isEquivalentAlreadySelected(genome.accession)
              const equivalentAccession = hasEquivalent ? getEquivalentSelected(genome.accession) : null
              
              return (
                <div 
                  key={genome.accession} 
                  onClick={() => toggleSelection(genome.accession)}
                  className={`flex items-center gap-2 p-3 rounded-lg border transition-all ${
                    hasEquivalent
                      ? 'bg-orange-50/50 border-orange-300 cursor-not-allowed opacity-60'
                      : selected
                        ? 'bg-teal-50 border-teal-500 shadow-md cursor-pointer'
                        : 'bg-slate-50 border-slate-200 hover:border-teal-300 cursor-pointer'
                  }`}
                >
                  {/* Checkbox */}
                  <div className={`w-5 h-5 rounded border-2 flex items-center justify-center flex-shrink-0 transition-all ${
                    hasEquivalent
                      ? 'bg-orange-200 border-orange-400'
                      : selected 
                        ? 'bg-teal-500 border-teal-500' 
                        : 'border-slate-300'
                  }`}>
                    {selected && !hasEquivalent && (
                      <svg className="w-3 h-3 text-white" fill="currentColor" viewBox="0 0 20 20">
                        <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                      </svg>
                    )}
                    {hasEquivalent && (
                      <svg className="w-3 h-3 text-orange-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01" />
                      </svg>
                    )}
                  </div>
                  
                  <div className="flex-1 min-w-0">
                    <div className={`font-mono text-sm font-semibold truncate ${
                      hasEquivalent ? 'text-orange-600' : selected ? 'text-teal-700' : 'text-slate-700'
                    }`}>
                      {genome.accession}
                    </div>
                    {hasEquivalent ? (
                      <div className="text-xs text-orange-600 truncate">
                        ≈ {equivalentAccession}
                      </div>
                    ) : genome.organism_name ? (
                      <div className="text-xs text-slate-500 truncate">
                        {genome.organism_name}
                      </div>
                    ) : null}
                  </div>
                  <button
                    onClick={(e) => {
                      e.stopPropagation() // Evitar que se dispare toggleSelection
                      handleDeleteGenome(genome.accession)
                    }}
                    className="flex-shrink-0 p-2 text-red-600 hover:bg-red-50 rounded-lg transition-colors"
                    title="Eliminar genoma"
                  >
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
                    </svg>
                  </button>
                </div>
              )
            })}
          </div>
        </div>
      )}

      {/* Cepas de E. coli predefinidas */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <div className="flex items-center justify-between mb-4">
          <div>
            <h3 className="font-semibold text-slate-800">Cepas de E. coli para Comparar</h3>
            <p className="text-sm text-slate-500">Selecciona 1 o más cepas • Descarga las que necesites</p>
          </div>
          <div className="flex gap-2">
            <button
              onClick={shuffleStrains}
              className="text-sm text-slate-600 hover:text-slate-800 font-medium flex items-center gap-1"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
              Otras cepas
            </button>
            <button
              onClick={selectAllStrains}
              className="text-sm text-teal-600 hover:text-teal-700 font-medium"
            >
              Seleccionar todas (12)
            </button>
          </div>
        </div>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
          {relatedStrains.map(strain => {
            const downloaded = isDownloaded(strain.accession)
            const selected = isSelected(strain.accession)
            const downloadStatus = getDownloadStatus(strain.accession)
            const hasEquivalent = isEquivalentAlreadySelected(strain.accession)
            const equivalentAccession = hasEquivalent ? getEquivalentSelected(strain.accession) : null
            
            return (
              <div
                key={strain.accession}
                onClick={() => toggleSelection(strain.accession)}
                className={`p-3 rounded-xl border-2 text-left transition-all ${
                  hasEquivalent
                    ? 'border-orange-300 bg-orange-50/50 cursor-not-allowed opacity-60'
                    : selected 
                      ? 'border-teal-500 bg-teal-50 shadow-md cursor-pointer' 
                      : 'border-slate-200 hover:border-teal-300 hover:shadow-sm cursor-pointer'
                }`}
              >
                <div className="flex items-start gap-2 mb-2">
                  <div className={`w-5 h-5 rounded border-2 flex items-center justify-center flex-shrink-0 mt-0.5 transition-all ${
                    hasEquivalent
                      ? 'bg-orange-200 border-orange-400'
                      : selected 
                        ? 'bg-teal-500 border-teal-500' 
                        : 'border-slate-300'
                  }`}>
                    {selected && !hasEquivalent && (
                      <svg className="w-3 h-3 text-white" fill="currentColor" viewBox="0 0 20 20">
                        <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                      </svg>
                    )}
                    {hasEquivalent && (
                      <svg className="w-3 h-3 text-orange-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                      </svg>
                    )}
                  </div>
                  
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center justify-between">
                      <span className={`font-mono text-sm font-medium ${hasEquivalent ? 'text-orange-600' : 'text-teal-700'}`}>{strain.strain}</span>
                      {hasEquivalent ? (
                        <span className="text-xs px-1.5 py-0.5 rounded bg-orange-100 text-orange-700">
                          ≈ {equivalentAccession?.slice(-8)}
                        </span>
                      ) : (
                        <span className={`text-xs px-1.5 py-0.5 rounded ${
                          strain.category === 'pathogenic' ? 'bg-red-100 text-red-700' :
                          strain.category === 'laboratory' ? 'bg-blue-100 text-blue-700' :
                          strain.category === 'industrial' ? 'bg-amber-100 text-amber-700' :
                          'bg-green-100 text-green-700'
                        }`}>
                          {strain.category}
                        </span>
                      )}
                    </div>
                    <p className="text-xs text-slate-600 mt-1 truncate">{strain.organism}</p>
                  </div>
                </div>
                
                <div className="text-xs font-medium">
                  {hasEquivalent ? (
                    <span className="text-orange-600 flex items-center gap-1">
                      <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01" />
                      </svg>
                      Equivalente ya seleccionado
                    </span>
                  ) : downloaded ? (
                    <span className="text-green-600 flex items-center gap-1">
                      <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
                        <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                      </svg>
                      Descargado y listo
                    </span>
                  ) : downloadStatus === 'downloading' ? (
                    <span className="text-blue-600 flex items-center gap-1">
                      <svg className="animate-spin h-3 w-3" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                      </svg>
                      Descargando...
                    </span>
                  ) : downloadStatus === 'completed' ? (
                    <span className="text-green-600">✓ Descarga completa</span>
                  ) : downloadStatus === 'error' ? (
                    <span className="text-red-600">✗ Error en descarga</span>
                  ) : selected ? (
                    <span className="text-slate-500">Se descargará</span>
                  ) : (
                    <span className="text-slate-400">No descargado</span>
                  )}
                </div>
              </div>
            )
          })}
        </div>
      </div>

      {/* Búsqueda personalizada */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <div className="flex items-center justify-between mb-4">
          <div>
            <h3 className="font-semibold text-slate-800">🔍 Buscar Cualquier Genoma en NCBI</h3>
            <p className="text-sm text-slate-500 mt-1">
              Busca por nombre completo, género, especie o términos comunes
            </p>
          </div>
          {searchResults.length > 0 && (
            <button
              onClick={clearSearch}
              className="text-sm text-slate-500 hover:text-slate-700 flex items-center gap-1"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
              Limpiar
            </button>
          )}
        </div>
        
        <div className="space-y-3 mb-4">
          <div className="flex gap-3">
            <div className="flex-1 relative">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={(e) => e.key === 'Enter' && searchGenomes()}
              placeholder="Ej: coli, Salmonella, Bacillus, GCF_000005845..."
                className="w-full px-4 py-3 border border-slate-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-teal-500 text-sm"
              />
              {isSearching && (
                <div className="absolute right-3 top-1/2 -translate-y-1/2">
                  <svg className="animate-spin h-5 w-5 text-teal-600" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                </div>
              )}
            </div>
            <button
              onClick={() => searchGenomes()}
              disabled={isSearching || !searchQuery.trim()}
              className="px-6 py-3 bg-teal-600 text-white rounded-lg font-medium hover:bg-teal-700 disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
            >
              {isSearching ? (
                <>
                  <svg className="animate-spin h-4 w-4" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                  Buscando...
                </>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                  </svg>
                  Buscar
                </>
              )}
            </button>
          </div>

          {/* Indicador de búsqueda dinámica */}
          {searchQuery.length > 0 && searchQuery.length < 3 && (
            <div className="text-xs text-slate-500 mb-2">
              💡 Escribe al menos 3 caracteres para buscar automáticamente
            </div>
          )}

          {/* Límite de resultados */}
          <div className="flex items-center gap-3 text-sm">
            <span className="text-slate-600">Resultados máximos:</span>
            <div className="flex gap-2">
              {[20, 50, 100].map(limit => (
                <button
                  key={limit}
                  onClick={() => setSearchLimit(limit)}
                  className={`px-3 py-1 rounded-lg font-medium transition-all ${
                    searchLimit === limit
                      ? 'bg-teal-600 text-white'
                      : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                  }`}
                >
                  {limit}
                </button>
              ))}
            </div>
          </div>
        </div>

        {/* Sugerencias de búsqueda */}
        {showSuggestions && searchResults.length === 0 && !isSearching && (
          <div className="mb-4">
            <p className="text-xs text-slate-500 mb-3 font-medium">HAZ CLIC PARA BUSCAR:</p>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
              {searchSuggestions.map(suggestion => (
                <button
                  key={suggestion.term}
                  onClick={() => searchGenomes(suggestion.term)}
                  className="p-3 bg-slate-50 hover:bg-teal-50 border border-slate-200 hover:border-teal-300 rounded-lg text-left transition-all group"
                >
                  <div className="flex items-center gap-2 mb-1">
                    <span className="text-lg">{suggestion.icon}</span>
                    <span className="text-xs font-medium text-slate-700 group-hover:text-teal-700">
                      {suggestion.desc}
                    </span>
                  </div>
                  <p className="text-xs text-slate-500 italic truncate">{suggestion.term}</p>
                </button>
              ))}
            </div>
            <p className="text-xs text-blue-600 mt-3 bg-blue-50 p-2 rounded border border-blue-200">
              💡 <strong>Puedes buscar por:</strong> Nombres completos (Escherichia coli), 
              términos cortos (coli, salmonella), o números de accesión (GCF_...)
            </p>
          </div>
        )}

        {/* Resultados de búsqueda */}
        {searchResults.length > 0 && (
          <div className="space-y-2">
            <div className="flex items-center justify-between text-sm text-slate-600 mb-3 pb-2 border-b border-slate-200">
              <span className="font-medium">
                📊 {searchResults.length} genomas encontrados
                {searchResults.length >= searchLimit && ` (límite: ${searchLimit})`}
              </span>
              <span className="text-xs text-slate-500">Haz clic para seleccionar</span>
            </div>
            
            <div className="max-h-96 overflow-y-auto space-y-2 pr-2">
              {searchResults.map(result => {
                const downloaded = isDownloaded(result.accession)
                const selected = isSelected(result.accession)
                const hasEquivalent = isEquivalentAlreadySelected(result.accession)
                const equivalentAccession = hasEquivalent ? getEquivalentSelected(result.accession) : null
                
                return (
                  <div
                    key={result.accession}
                    onClick={() => toggleSelection(result.accession)}
                    className={`p-4 rounded-lg border-2 transition-all ${
                      hasEquivalent
                        ? 'border-orange-300 bg-orange-50/50 cursor-not-allowed opacity-60'
                        : selected
                          ? 'border-teal-500 bg-teal-50 shadow-md cursor-pointer'
                          : 'border-slate-200 hover:border-teal-300 hover:bg-slate-50 cursor-pointer'
                    }`}
                  >
                    <div className="flex items-start gap-3">
                      <div className={`w-5 h-5 rounded border-2 flex items-center justify-center flex-shrink-0 mt-0.5 transition-all ${
                        hasEquivalent
                          ? 'bg-orange-200 border-orange-400'
                          : selected 
                            ? 'bg-teal-500 border-teal-500' 
                            : 'border-slate-300'
                      }`}>
                        {selected && !hasEquivalent && (
                          <svg className="w-3 h-3 text-white" fill="currentColor" viewBox="0 0 20 20">
                            <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                          </svg>
                        )}
                        {hasEquivalent && (
                          <svg className="w-3 h-3 text-orange-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                          </svg>
                        )}
                      </div>
                      
                      <div className="flex-1 min-w-0">
                        <div className="flex items-start justify-between gap-2 mb-1">
                          <div className="flex items-center gap-2 flex-wrap">
                            <span className={`font-mono text-sm font-semibold ${hasEquivalent ? 'text-orange-600' : 'text-teal-700'}`}>{result.accession}</span>
                            {hasEquivalent && (
                              <span className="text-xs bg-orange-100 text-orange-700 px-2 py-0.5 rounded-full font-medium">
                                ⚠️ Equivalente a {equivalentAccession}
                              </span>
                            )}
                            {result.is_reference && !hasEquivalent && (
                              <span className="text-xs bg-amber-100 text-amber-700 px-2 py-0.5 rounded-full font-medium">
                                ⭐ Referencia
                              </span>
                            )}
                            {downloaded && !hasEquivalent && (
                              <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded-full font-medium">
                                ✓ Descargado
                              </span>
                            )}
                            {selected && !downloaded && !hasEquivalent && (
                              <span className="text-xs bg-blue-100 text-blue-700 px-2 py-0.5 rounded-full font-medium">
                                Se descargará
                              </span>
                            )}
                          </div>
                        </div>
                        
                        <p className="text-sm text-slate-800 font-medium mb-1">
                          {result.organism_name}
                          {result.strain && <span className="text-slate-500"> • {result.strain}</span>}
                        </p>
                        
                        <div className="flex items-center gap-3 text-xs text-slate-500 mt-2 flex-wrap">
                          {result.genome_size_mb && (
                            <span className="flex items-center gap-1">
                              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 7v10c0 2 1 3 3 3h10c2 0 3-1 3-3V7c0-2-1-3-3-3H7C5 4 4 5 4 7z" />
                              </svg>
                              {result.genome_size_mb} Mb
                            </span>
                          )}
                          {result.assembly_level && (
                            <span className="flex items-center gap-1">
                              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
                              </svg>
                              {result.assembly_level}
                            </span>
                          )}
                          {result.assembly_name && (
                            <span className="flex items-center gap-1 truncate">
                              <svg className="w-3 h-3 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 7h.01M7 3h5c.512 0 1.024.195 1.414.586l7 7a2 2 0 010 2.828l-7 7a2 2 0 01-2.828 0l-7-7A1.994 1.994 0 013 12V7a4 4 0 014-4z" />
                              </svg>
                              {result.assembly_name}
                            </span>
                          )}
                          {result.submission_date && (
                            <span className="flex items-center gap-1">
                              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 7V3m8 4V3m-9 8h10M5 21h14a2 2 0 002-2V7a2 2 0 00-2-2H5a2 2 0 00-2 2v12a2 2 0 002 2z" />
                              </svg>
                              {new Date(result.submission_date).getFullYear()}
                            </span>
                          )}
                        </div>
                      </div>
                    </div>
                  </div>
                )
              })}
            </div>
            
            {searchResults.length >= searchLimit && (
              <div className="mt-3 p-3 bg-blue-50 border border-blue-200 rounded-lg">
                <p className="text-sm text-blue-700">
                  💡 <strong>Mostrando {searchLimit} resultados.</strong> Aumenta el límite o refina tu búsqueda para ver más genomas.
                </p>
              </div>
            )}
          </div>
        )}
        
        {!isSearching && searchQuery && searchResults.length === 0 && !showSuggestions && (
          <div className="text-center py-8 text-slate-500">
            <svg className="w-12 h-12 mx-auto mb-3 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
            </svg>
            <p className="text-sm font-medium">No se encontraron genomas para "{searchQuery}"</p>
            <div className="mt-4 bg-blue-50 border border-blue-200 rounded-lg p-4 text-left max-w-md mx-auto">
              <p className="text-xs font-semibold text-blue-800 mb-2">✓ Intenta buscar por:</p>
              <ul className="text-xs text-blue-700 space-y-1">
                <li>• <strong>Términos comunes:</strong> coli, salmonella, bacillus</li>
                <li>• <strong>Nombres completos:</strong> Escherichia coli, Bacillus subtilis</li>
                <li>• <strong>Número de accesión:</strong> GCF_000005845</li>
              </ul>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
