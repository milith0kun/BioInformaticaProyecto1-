/**
 * Protein3DViewer Component
 * Visualizador 3D integrado de estructuras de proteínas usando 3Dmol.js
 * Soporta descarga automática desde PDB y AlphaFold
 */
import { useEffect, useRef, useState } from 'react'

// 3Dmol.js se carga globalmente desde CDN (ver index.html)
const get3Dmol = () => {
  if (typeof window !== 'undefined' && window.$3Dmol) {
    return window.$3Dmol
  }
  return null
}

const RCSB_PDB_API = 'https://www.rcsb.org/structure'
const ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/entry'
const PDB_FILE_API = 'https://files.rcsb.org/download'
const ALPHAFOLD_FILE_API = 'https://alphafold.ebi.ac.uk/files'

export default function Protein3DViewer({ proteinId, proteinSequence, productName }) {
  const viewerRef = useRef(null)
  const containerRef = useRef(null)
  const [viewer, setViewer] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [structureLoaded, setStructureLoaded] = useState(false)
  const [structureSource, setStructureSource] = useState(null) // 'pdb', 'alphafold', or 'predicted'
  const [searchQuery, setSearchQuery] = useState('')
  const [viewStyle, setViewStyle] = useState('cartoon')
  const [colorScheme, setColorScheme] = useState('spectrum')

  // Inicializar el visualizador 3Dmol
  useEffect(() => {
    if (containerRef.current && !viewer) {
      const $3Dmol = get3Dmol()

      if (!$3Dmol) {
        console.warn('3Dmol.js no está disponible aún, reintentando...')
        // Reintentar después de un corto delay
        const timer = setTimeout(() => {
          const retry = get3Dmol()
          if (retry && containerRef.current) {
            try {
              const config = { backgroundColor: '#1e293b' }
              const glviewer = retry.createViewer(containerRef.current, config)
              setViewer(glviewer)
              viewerRef.current = glviewer
            } catch (err) {
              console.error('Error inicializando visualizador 3D:', err)
              setError('El visualizador 3D no está disponible. La librería 3Dmol.js no se cargó correctamente.')
            }
          } else {
            setError('El visualizador 3D no está disponible. Por favor, recargue la página.')
          }
        }, 500)
        return () => clearTimeout(timer)
      }

      try {
        const config = { backgroundColor: '#1e293b' }
        const glviewer = $3Dmol.createViewer(containerRef.current, config)
        setViewer(glviewer)
        viewerRef.current = glviewer
      } catch (err) {
        console.error('Error inicializando visualizador 3D:', err)
        setError('No se pudo inicializar el visualizador 3D')
      }
    }

    // Cleanup cuando se desmonta el componente
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear()
        } catch (err) {
          console.error('Error limpiando visualizador:', err)
        }
      }
    }
  }, [])

  // Buscar y cargar estructura desde PDB
  const searchPDB = async (query, silent = false) => {
    setLoading(true)
    setError(null)
    setStructureLoaded(false)

    try {
      // Buscar en RCSB PDB
      const searchUrl = `https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify({
        query: {
          type: 'terminal',
          service: 'text',
          parameters: {
            attribute: 'struct.title',
            operator: 'contains_words',
            value: query
          }
        },
        return_type: 'entry',
        request_options: {
          results_content_type: ['experimental'],
          return_all_hits: false,
          pager: { start: 0, rows: 1 }
        }
      }))}`

      const searchResponse = await fetch(searchUrl)
      if (!searchResponse.ok) throw new Error('No se pudo buscar en PDB')

      const searchData = await searchResponse.json()

      if (!searchData.result_set || searchData.result_set.length === 0) {
        throw new Error('No se encontraron estructuras en PDB')
      }

      const pdbId = searchData.result_set[0].identifier
      await loadPDBStructure(pdbId)

    } catch (err) {
      if (!silent) {
        console.log('ℹ️ No se encontró en PDB, intentando AlphaFold...')
      }
      // Intentar con AlphaFold
      await tryAlphaFold(query, silent)
    }
  }

  // Cargar estructura desde PDB por ID
  const loadPDBStructure = async (pdbId) => {
    if (!viewer) return

    try {
      setLoading(true)
      setError(null)

      const pdbUrl = `${PDB_FILE_API}/${pdbId}.pdb`
      const response = await fetch(pdbUrl)

      if (!response.ok) throw new Error('No se pudo descargar la estructura PDB')

      const pdbData = await response.text()

      viewer.clear()
      viewer.addModel(pdbData, 'pdb')
      applyStyle()
      viewer.zoomTo()
      viewer.render()

      setStructureLoaded(true)
      setStructureSource('pdb')
      setLoading(false)
    } catch (err) {
      console.error('Error cargando PDB:', err)
      setError('No se pudo cargar la estructura desde PDB')
      setLoading(false)
    }
  }

  // Intentar con AlphaFold
  const tryAlphaFold = async (uniprotId, silent = false) => {
    if (!viewer) return

    try {
      setLoading(true)
      setError(null)

      // Intentar descargar desde AlphaFold usando UniProt ID
      const alphafoldUrl = `${ALPHAFOLD_FILE_API}/AF-${uniprotId}-F1-model_v4.pdb`
      const response = await fetch(alphafoldUrl, { mode: 'cors' })

      if (!response.ok) throw new Error('No disponible en AlphaFold')

      const pdbData = await response.text()

      viewer.clear()
      viewer.addModel(pdbData, 'pdb')
      applyStyle()
      viewer.zoomTo()
      viewer.render()

      setStructureLoaded(true)
      setStructureSource('alphafold')
      setLoading(false)
    } catch (err) {
      if (!silent) {
        console.log('ℹ️ No se encontró en AlphaFold, generando estructura predicha...')
      }
      // Si falla todo, generar estructura predicha simple
      generatePredictedStructure()
    }
  }

  // Generar estructura predicha simple basada en secuencia
  const generatePredictedStructure = () => {
    if (!viewer || !proteinSequence) {
      setError('No hay datos suficientes para generar una estructura')
      setLoading(false)
      return
    }

    try {
      // Crear un modelo simple en espiral alfa-hélice
      let pdbContent = 'HEADER    PREDICTED STRUCTURE\n'
      pdbContent += 'TITLE     ESTRUCTURA PREDICHA (APROXIMADA)\n'

      const sequence = proteinSequence.substring(0, Math.min(200, proteinSequence.length))

      // Generar coordenadas en espiral (alfa hélice)
      const radius = 2.3  // Radio de la hélice
      const rise = 1.5    // Altura por residuo
      const angle = 100   // Ángulo por residuo (grados)

      sequence.split('').forEach((aa, i) => {
        const theta = (angle * i * Math.PI) / 180
        const x = radius * Math.cos(theta)
        const y = radius * Math.sin(theta)
        const z = i * rise

        // Formato PDB: ATOM, número, tipo, residuo, cadena, número residuo, x, y, z
        pdbContent += `ATOM  ${String(i + 1).padStart(5, ' ')}  CA  ${aa.padEnd(3, ' ')} A${String(i + 1).padStart(4, ' ')}    `
        pdbContent += `${x.toFixed(3).padStart(8, ' ')}${y.toFixed(3).padStart(8, ' ')}${z.toFixed(3).padStart(8, ' ')}\n`
      })

      pdbContent += 'END\n'

      viewer.clear()
      viewer.addModel(pdbContent, 'pdb')
      applyStyle()
      viewer.zoomTo()
      viewer.render()

      setStructureLoaded(true)
      setStructureSource('predicted')
      setLoading(false)
    } catch (err) {
      console.error('Error generando estructura:', err)
      setError('No se pudo generar estructura predicha')
      setLoading(false)
    }
  }

  // Aplicar estilo de visualización
  const applyStyle = () => {
    if (!viewer) return

    const models = viewer.getModel()
    if (!models) return

    // Limpiar estilos previos
    viewer.setStyle({}, {})

    // Aplicar nuevo estilo
    const styleConfig = {}

    switch (viewStyle) {
      case 'cartoon':
        styleConfig.cartoon = { color: colorScheme }
        break
      case 'sphere':
        styleConfig.sphere = { colorscheme: colorScheme, scale: 0.3 }
        break
      case 'stick':
        styleConfig.stick = { colorscheme: colorScheme }
        break
      case 'line':
        styleConfig.line = { colorscheme: colorScheme }
        break
      case 'cross':
        styleConfig.cross = { colorscheme: colorScheme }
        break
      default:
        styleConfig.cartoon = { color: colorScheme }
    }

    viewer.setStyle({}, styleConfig)
    viewer.render()
  }

  // Aplicar estilo cuando cambian las opciones
  useEffect(() => {
    if (structureLoaded && viewer) {
      try {
        applyStyle()
      } catch (err) {
        console.error('Error aplicando estilo:', err)
      }
    }
  }, [viewStyle, colorScheme, structureLoaded, viewer])

  // Búsqueda automática cuando se proporciona protein ID
  useEffect(() => {
    if (proteinId && viewer && proteinSequence) {
      setSearchQuery(proteinId)

      // Detectar si es un ID de RefSeq/NCBI (WP_, YP_, NP_, etc.)
      const isRefSeqId = /^(WP_|YP_|NP_|XP_|AP_)/.test(proteinId)

      if (isRefSeqId) {
        // IDs de RefSeq raramente tienen estructura en PDB/AlphaFold
        // Generar predicción directamente
        console.log(`🧬 ID de RefSeq detectado (${proteinId}), generando estructura predicha...`)
        generatePredictedStructure()
      } else {
        // Intentar buscar en PDB/AlphaFold solo para IDs que podrían existir
        searchPDB(proteinId)
      }
    } else if (!proteinId && proteinSequence && viewer) {
      // Si no hay protein_id pero sí secuencia, generar predicción automáticamente
      generatePredictedStructure()
    }
  }, [proteinId, proteinSequence, viewer])

  const handleSearch = (e) => {
    e.preventDefault()
    if (searchQuery.trim()) {
      searchPDB(searchQuery.trim())
    }
  }

  const handleLoadPredicted = () => {
    generatePredictedStructure()
  }

  return (
    <div className="space-y-4">
      {/* Controles */}
      <div className="bg-white border border-slate-200 rounded-lg p-4">
        <div className="flex flex-col gap-3">
          {/* Búsqueda */}
          <form onSubmit={handleSearch} className="flex gap-2">
            <input
              type="text"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              placeholder="PDB ID, UniProt ID o nombre de proteína..."
              className="flex-1 px-3 py-2 border border-slate-200 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-teal-400"
            />
            <button
              type="submit"
              disabled={loading || !searchQuery.trim()}
              className="px-4 py-2 bg-teal-600 text-white rounded-lg hover:bg-teal-700 disabled:opacity-50 text-sm font-medium transition-colors"
            >
              {loading ? 'Buscando...' : 'Buscar'}
            </button>
            {proteinSequence && (
              <button
                type="button"
                onClick={handleLoadPredicted}
                disabled={loading}
                className="px-4 py-2 bg-purple-600 text-white rounded-lg hover:bg-purple-700 disabled:opacity-50 text-sm font-medium transition-colors"
              >
                Generar Predicción
              </button>
            )}
          </form>

          {/* Opciones de visualización */}
          {structureLoaded && (
            <div className="flex flex-wrap gap-3">
              <div className="flex items-center gap-2">
                <label className="text-xs font-medium text-slate-600">Estilo:</label>
                <select
                  value={viewStyle}
                  onChange={(e) => setViewStyle(e.target.value)}
                  className="px-2 py-1 border border-slate-200 rounded text-xs focus:outline-none focus:ring-2 focus:ring-teal-400"
                >
                  <option value="cartoon">Cartoon (Cinta)</option>
                  <option value="sphere">Esferas</option>
                  <option value="stick">Palos</option>
                  <option value="line">Líneas</option>
                  <option value="cross">Cruz</option>
                </select>
              </div>

              <div className="flex items-center gap-2">
                <label className="text-xs font-medium text-slate-600">Color:</label>
                <select
                  value={colorScheme}
                  onChange={(e) => setColorScheme(e.target.value)}
                  className="px-2 py-1 border border-slate-200 rounded text-xs focus:outline-none focus:ring-2 focus:ring-teal-400"
                >
                  <option value="spectrum">Espectro</option>
                  <option value="chain">Por cadena</option>
                  <option value="ss">Estructura 2ª</option>
                  <option value="residue">Por residuo</option>
                  <option value="hydrophobicity">Hidrofobicidad</option>
                  <option value="charge">Carga</option>
                </select>
              </div>

              {structureSource && (
                <div className="ml-auto flex items-center gap-2">
                  <span className="text-xs text-slate-500">Fuente:</span>
                  <span className={`text-xs font-medium px-2 py-1 rounded ${
                    structureSource === 'pdb'
                      ? 'bg-blue-100 text-blue-700'
                      : structureSource === 'alphafold'
                      ? 'bg-green-100 text-green-700'
                      : 'bg-purple-100 text-purple-700'
                  }`}>
                    {structureSource === 'pdb' && 'PDB Experimental'}
                    {structureSource === 'alphafold' && 'AlphaFold (IA)'}
                    {structureSource === 'predicted' && 'Predicción Simple'}
                  </span>
                </div>
              )}
            </div>
          )}
        </div>
      </div>

      {/* Visualizador 3D */}
      <div className="bg-slate-900 rounded-lg border border-slate-700 overflow-hidden">
        <div
          ref={containerRef}
          style={{
            width: '100%',
            height: '500px',
            position: 'relative'
          }}
        />

        {/* Overlay de estado */}
        {!structureLoaded && !loading && (
          <div className="absolute inset-0 flex items-center justify-center bg-slate-900/95">
            <div className="text-center p-8">
              <svg className="w-16 h-16 text-slate-600 mx-auto mb-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
              <h3 className="text-lg font-semibold text-slate-300 mb-2">Sin estructura cargada</h3>
              <p className="text-sm text-slate-500">
                Busca una estructura por PDB ID, UniProt ID o genera una predicción
              </p>
            </div>
          </div>
        )}

        {loading && (
          <div className="absolute inset-0 flex items-center justify-center bg-slate-900/95">
            <div className="text-center">
              <div className="w-16 h-16 border-4 border-teal-500 border-t-transparent rounded-full animate-spin mx-auto mb-4"></div>
              <p className="text-slate-300 text-sm">Cargando estructura 3D...</p>
            </div>
          </div>
        )}

        {error && (
          <div className="absolute top-4 left-4 right-4 bg-red-50 border border-red-200 rounded-lg p-3">
            <div className="flex items-start gap-2">
              <svg className="w-5 h-5 text-red-600 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" clipRule="evenodd" />
              </svg>
              <div className="flex-1">
                <p className="text-sm font-medium text-red-800">{error}</p>
                <p className="text-xs text-red-600 mt-1">
                  Intenta con otro ID o genera una predicción simple
                </p>
              </div>
              <button
                onClick={() => setError(null)}
                className="text-red-600 hover:text-red-800"
              >
                <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
                  <path fillRule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clipRule="evenodd" />
                </svg>
              </button>
            </div>
          </div>
        )}
      </div>

      {/* Información sobre predicción automática para RefSeq */}
      {structureSource === 'predicted' && proteinId && /^(WP_|YP_|NP_|XP_|AP_)/.test(proteinId) && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
          <div className="flex items-start gap-3">
            <svg className="w-5 h-5 text-amber-600 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
            <div className="flex-1">
              <h4 className="text-sm font-semibold text-amber-900 mb-1">Estructura Predicha - ID RefSeq</h4>
              <p className="text-xs text-amber-800">
                Esta proteína tiene un ID de RefSeq ({proteinId}), que típicamente no tiene estructura experimental en PDB ni AlphaFold.
                Se ha generado una <strong>predicción simple</strong> basada en hélices α de la secuencia primaria.
                Para estructuras más precisas, busca directamente por PDB ID o UniProt ID.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Información */}
      <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
        <h4 className="text-sm font-semibold text-blue-900 mb-2">💡 Cómo usar el visualizador 3D</h4>
        <div className="text-xs text-blue-800 space-y-1">
          <p>• <strong>Buscar:</strong> PDB ID (ej: 1MBO), UniProt ID (ej: P69905) o nombre de proteína</p>
          <p>• <strong>Rotar:</strong> Click izquierdo + arrastrar</p>
          <p>• <strong>Zoom:</strong> Rueda del mouse o pinch en móvil</p>
          <p>• <strong>Mover:</strong> Click derecho + arrastrar (o Ctrl + click izquierdo)</p>
          <p>• <strong>Predicción:</strong> Genera estructura aproximada (hélice α) de la secuencia</p>
          <p className="text-amber-700 mt-2">
            ⚠️ <strong>Nota:</strong> IDs de RefSeq (WP_, YP_, NP_) generan predicciones automáticamente
          </p>
        </div>
      </div>
    </div>
  )
}
