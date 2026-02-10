/**
 * MolstarProteinViewer Component — Enhanced
 * 3D protein structure visualizer with Mol*
 * Supports tertiary and quaternary structures
 */
import { useState, useEffect, useMemo, useRef } from 'react'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/build/viewer/molstar.css';

// Molstar Viewer Component with Enhanced Configuration
const Molstar = ({ pdbUrl }) => {
  const parentRef = useRef(null);
  const pluginRef = useRef(null);

  useEffect(() => {
    async function init() {
      if (!parentRef.current) return;

      // Enhanced spec to ensure quaternary structures are displayed
      const spec = DefaultPluginUISpec();
      spec.layout = {
          initial: {
              isExpanded: false,
              showControls: true,
              controlsDisplay: 'reactive'
          }
      };

      // Clear container to avoid "already passed to createRoot" warning
      parentRef.current.innerHTML = '';

      // Create plugin instance
      const plugin = await createPluginUI(parentRef.current, spec);
      pluginRef.current = plugin;

      if (pdbUrl) {
        try {
          // Load structure
          const data = await plugin.builders.data.download(
            { url: pdbUrl, isBinary: false },
            { state: { isGhost: true } }
          );

          // Parse as PDB
          const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');

          // IMPORTANT: Use 'all-models' preset to show ALL chains/assemblies
          // This ensures quaternary structures are properly displayed
          await plugin.builders.structure.hierarchy.applyPreset(
            trajectory,
            'all-models', // Changed from 'default' to show quaternary structures
            {
              representationPreset: 'auto', // Automatically choose best representation
              showUnitcell: false,
              showWater: false,
              showHet: true,
              ignoreHydrogens: true
            }
          );

        } catch (err) {
          console.error('Error loading structure in Molstar:', err);
        }
      }
    }

    init();

    return () => {
      if (pluginRef.current) {
        pluginRef.current.dispose();
        pluginRef.current = null;
      }
    };
  }, [pdbUrl]);

  return <div ref={parentRef} style={{ width: '100%', height: '100%', minHeight: '500px' }} />;
};

const PDB_FILE_API = 'https://files.rcsb.org/download'
const ALPHAFOLD_FILE_API = 'https://alphafold.ebi.ac.uk/files'

export default function MolstarProteinViewer({ proteinId, proteinSequence, productName }) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [structureUrl, setStructureUrl] = useState(null)
  const [structureSource, setStructureSource] = useState(null)
  const [searchQuery, setSearchQuery] = useState('')
  const [manualPdbId, setManualPdbId] = useState('')

  const hasAttemptedLoad = useRef(false)

  // Cleanup blob URLs
  useEffect(() => {
    return () => {
      if (structureUrl && structureUrl.startsWith('blob:')) {
        URL.revokeObjectURL(structureUrl)
      }
    }
  }, [structureUrl])

  // Auto-detect ID types
  const isRefSeqId = useMemo(() => {
    if (!proteinId) return false
    return /^(WP_|YP_|NP_|XP_|AP_)/.test(proteinId)
  }, [proteinId])

  const isPdbId = useMemo(() => {
    if (!proteinId) return false
    return /^[A-Z0-9]{4}$/i.test(proteinId)
  }, [proteinId])

  const isAlphaFoldId = useMemo(() => {
    if (!proteinId) return false
    return /^AF-/.test(proteinId)
  }, [proteinId])

  // Predict structure using ESMFold
  const predictStructure = async () => {
    try {
      setLoading(true)
      setError(null)

      if (!proteinSequence || proteinSequence.length < 10) {
        throw new Error('Secuencia de proteína no disponible o demasiado corta')
      }

      if (proteinSequence.length > 400) {
        throw new Error(`Secuencia demasiado larga (${proteinSequence.length} aa). ESMFold soporta máximo 400 aminoácidos.`)
      }

      console.log(`🧬 Prediciendo estructura para ${proteinId} (${proteinSequence.length} aa)...`)

      const response = await fetch(`/api/ncbi/protein/predict-structure/${proteinId}`, {
        method: 'POST'
      })

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}))
        throw new Error(errorData.detail || `Error ${response.status} prediciendo estructura`)
      }

      const pdbContent = await response.text()
      const blob = new Blob([pdbContent], { type: 'chemical/x-pdb' })
      const pdbUrl = URL.createObjectURL(blob)

      setStructureUrl(pdbUrl)
      setStructureSource('esmfold')
      setError(null)
      console.log('✅ Estructura predicha con ESMFold')
      return true
    } catch (err) {
      setError(err.message)
      console.error('Error prediciendo estructura:', err)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Try AlphaFold DB
  const tryAlphaFold = async (uniprotId) => {
    try {
      setLoading(true)
      setError(null)

      const searchUrl = `https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`
      const response = await fetch(searchUrl)

      if (!response.ok) throw new Error('No encontrado en AlphaFold')

      const data = await response.json()

      if (data && data[0]) {
        const entry = data[0]
        const pdbUrl = entry.pdbUrl || `${ALPHAFOLD_FILE_API}/AF-${uniprotId}-F1-model_v4.pdb`

        setStructureUrl(pdbUrl)
        setStructureSource('alphafold')
        setError(null)
        console.log('✅ Estructura de AlphaFold cargada:', pdbUrl)
        return true
      }

      throw new Error('No se encontró predicción en AlphaFold')
    } catch (err) {
      console.warn('AlphaFold no disponible:', err.message)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Search AlphaFold DB by RefSeq
  const searchAlphaFoldByRefSeq = async () => {
    try {
      setLoading(true)
      setError(null)

      const searchUrl = `https://www.alphafold.ebi.ac.uk/api/prediction/${proteinId}`
      const response = await fetch(searchUrl)

      if (response.ok) {
        const data = await response.json()
        if (data && data[0]) {
          const entry = data[0]
          const pdbUrl = entry.pdbUrl || `${ALPHAFOLD_FILE_API}/AF-${proteinId}-F1-model_v4.pdb`

          setStructureUrl(pdbUrl)
          setStructureSource('alphafold')
          setError(null)
          console.log('✅ Estructura encontrada en AlphaFold DB')
          return true
        }
      }

      throw new Error('No encontrado en AlphaFold DB')
    } catch (err) {
      console.warn('AlphaFold DB búsqueda fallida:', err.message)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Load from PDB
  const loadFromPDB = async (pdbId) => {
    try {
      setLoading(true)
      setError(null)

      const cleanPdbId = pdbId.trim().toUpperCase()
      if (!/^[A-Z0-9]{4}$/i.test(cleanPdbId)) {
        throw new Error(`ID de PDB inválido: "${pdbId}". Debe tener exactamente 4 caracteres (ej: 1CRN, 4HHB)`)
      }

      const pdbUrl = `${PDB_FILE_API}/${cleanPdbId}.pdb`

      const response = await fetch(pdbUrl, { method: 'HEAD' })
      if (!response.ok) {
        throw new Error(`Estructura PDB ${cleanPdbId} no encontrada (${response.status})`)
      }

      setStructureUrl(pdbUrl)
      setStructureSource('pdb')
      setError(null)
      console.log('✅ Estructura PDB cargada:', pdbUrl)
      return true
    } catch (err) {
      setError(err.message)
      console.error('Error cargando PDB:', err)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Search PDB by name
  const searchPDB = async (query) => {
    try {
      setLoading(true)
      setError(null)

      const sanitizedQuery = query.trim().slice(0, 100)
      if (sanitizedQuery.length < 2) {
        throw new Error('La búsqueda debe tener al menos 2 caracteres')
      }

      const searchPayload = {
        query: {
          type: 'terminal',
          service: 'text',
          parameters: {
            attribute: 'struct.title',
            operator: 'contains_phrase',
            value: sanitizedQuery
          }
        },
        return_type: 'entry',
        request_options: {
          results_content_type: ['experimental'],
          sort: [{ sort_by: 'score', direction: 'desc' }],
          pager: { start: 0, rows: 5 }
        }
      }

      const searchUrl = `https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(searchPayload))}`
      const response = await fetch(searchUrl)

      if (!response.ok) {
        throw new Error(`Error en búsqueda PDB: ${response.status} ${response.statusText}`)
      }

      const data = await response.json()

      if (data.result_set && data.result_set.length > 0) {
        const firstResult = data.result_set[0].identifier
        console.log('🔍 Encontrado en PDB:', firstResult)
        await loadFromPDB(firstResult)
        return true
      }

      throw new Error('No se encontraron resultados en PDB para esta búsqueda')
    } catch (err) {
      setError(err.message)
      console.error('Error buscando en PDB:', err)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Auto-load on mount
  useEffect(() => {
    if (!proteinId) return
    if (hasAttemptedLoad.current) return
    hasAttemptedLoad.current = true

    const autoLoad = async () => {
      if (isPdbId) {
        await loadFromPDB(proteinId)
        return
      }

      if (isAlphaFoldId) {
        const afId = proteinId.replace('AF-', '').split('-')[0]
        await tryAlphaFold(afId)
        return
      }

      if (isRefSeqId && proteinSequence) {
        const seqLength = proteinSequence.length

        // Try AlphaFold DB first (supports large proteins)
        console.log(`🧬 ID de RefSeq detectado (${proteinId}, ${seqLength} aa). Buscando en AlphaFold DB...`)
        const foundInAlphaFold = await searchAlphaFoldByRefSeq()
        if (foundInAlphaFold) return

        // If not found and <= 400 aa, predict with ESMFold
        if (seqLength <= 400) {
          console.log(`🧬 Prediciendo estructura con ESMFold (${seqLength} aa)...`)
          await predictStructure()
          return
        }

        // Large proteins: show message
        console.log(`⚠️ Proteína grande (${seqLength} aa). ESMFold tiene límite de 400 aa.`)
        setError(`Proteína grande (${seqLength} aa). Usa la búsqueda manual en PDB o busca por nombre en AlphaFold DB.`)
        return
      }

      if (productName) {
        console.log(`🔍 Buscando estructura para: ${productName}`)
        await searchPDB(productName)
      }
    }

    autoLoad()

    return () => {
      hasAttemptedLoad.current = false
    }
  }, [proteinId, productName, proteinSequence, isPdbId, isAlphaFoldId, isRefSeqId])

  // Handlers
  const handleSearch = async (e) => {
    e.preventDefault()
    if (!searchQuery.trim()) return
    await searchPDB(searchQuery.trim())
  }

  const handleManualPdbLoad = async (e) => {
    e.preventDefault()
    if (!manualPdbId.trim()) return
    await loadFromPDB(manualPdbId.trim())
  }

  return (
    <div className="space-y-6">
      {/* Dynamic Prediction/Search HUD */}
      {!structureUrl && !loading && (
        <div className="bg-slate-800/20 backdrop-blur-md rounded-[2rem] border border-white/5 p-10 text-center animate-in fade-in zoom-in-95 duration-700">
          <div className="w-24 h-24 bg-gradient-to-br from-cyan-500/20 to-blue-500/20 rounded-full flex items-center justify-center mx-auto mb-8 border border-cyan-500/30">
            <span className="text-4xl animate-pulse">🧬</span>
          </div>
          
          <h3 className="text-xl font-black text-white uppercase tracking-tighter mb-4">Resolución de Estructura Requerida</h3>
          <p className="text-slate-400 text-sm max-w-md mx-auto mb-10 leading-relaxed">
            No hay datos 3D en caché para esta proteína. Elige un método de resolución para inicializar el visor molecular.
          </p>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 max-w-2xl mx-auto">
            {/* ESMFold Option */}
            <button
              onClick={predictStructure}
              disabled={!proteinSequence || proteinSequence.length > 400}
              className="group relative p-6 bg-slate-900/50 rounded-3xl border border-white/5 text-left transition-all hover:border-cyan-500/50 hover:bg-slate-800 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              <div className="flex items-center gap-4 mb-3">
                <div className="p-3 bg-cyan-500/10 rounded-xl text-cyan-400 font-black text-xs">IA</div>
                <h4 className="text-xs font-black text-slate-200 uppercase tracking-widest">Predicción ESMFold</h4>
              </div>
              <p className="text-[10px] text-slate-500 font-medium">
                Genera un modelo 3D de novo usando ESMFold (máx. 400 aa).
              </p>
              {proteinSequence?.length > 400 && (
                <p className="text-[9px] text-amber-500 mt-2 font-bold uppercase tracking-tighter">Secuencia muy larga (límite 400aa)</p>
              )}
            </button>

            {/* PDB Option */}
            <div className="group relative p-6 bg-slate-900/50 rounded-3xl border border-white/5 text-left transition-all hover:border-indigo-500/50 hover:bg-slate-800">
              <div className="flex items-center gap-4 mb-3">
                <div className="p-3 bg-indigo-500/10 rounded-xl text-indigo-400 font-black text-xs">EXP</div>
                <h4 className="text-xs font-black text-slate-200 uppercase tracking-widest">Importación PDB</h4>
              </div>
              <form onSubmit={handleManualPdbLoad} className="flex gap-2">
                <input
                  type="text"
                  value={manualPdbId}
                  onChange={(e) => setManualPdbId(e.target.value)}
                  placeholder="Ingresa ID de PDB..."
                  className="flex-1 bg-slate-950 border border-white/5 rounded-xl px-3 py-2 text-[10px] font-mono text-cyan-400 focus:outline-none focus:border-indigo-500"
                  maxLength={4}
                />
                <button type="submit" className="p-2 bg-indigo-600 rounded-xl text-white">
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M19 13l-7 7-7-7m14-8l-7 7-7-7" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}/></svg>
                </button>
              </form>
            </div>
          </div>
        </div>
      )}

      {/* Loading Overlay */}
      {loading && (
        <div className="h-[600px] bg-slate-950/50 backdrop-blur-sm rounded-[2.5rem] flex flex-col items-center justify-center border border-white/5">
          <div className="relative w-20 h-20 mb-6">
            <div className="absolute inset-0 border-4 border-slate-800 rounded-full"></div>
            <div className="absolute inset-0 border-4 border-transparent border-t-cyan-500 rounded-full animate-spin"></div>
          </div>
          <p className="text-slate-500 text-[10px] font-black uppercase tracking-[0.3em] animate-pulse">Procesando Malla Molecular...</p>
        </div>
      )}

      {/* Error HUD */}
      {error && !structureUrl && !loading && (
        <div className="p-6 bg-amber-500/10 rounded-2xl border border-amber-500/20 text-center">
          <p className="text-amber-300 text-xs font-bold uppercase tracking-widest mb-4">Error de Resolución: {error}</p>
          <button onClick={() => {setError(null); setSearchQuery('');}} className="text-[10px] text-cyan-400 font-black uppercase hover:underline">Reintentar</button>
        </div>
      )}

      {/* Main Viewer Canvas */}
      {structureUrl && (
        <div className="relative group rounded-[2rem] overflow-hidden border border-white/10 shadow-3xl bg-black">
          {/* Viewer HUD Overlay */}
          <div className="absolute top-6 left-6 z-10 flex gap-3 pointer-events-none">
            <div className="px-4 py-2 bg-slate-900/80 backdrop-blur-xl border border-white/5 rounded-full flex items-center gap-3">
              <div className="w-2 h-2 bg-emerald-500 rounded-full shadow-[0_0_8px_rgba(16,185,129,0.5)]"></div>
              <span className="text-[9px] font-black text-white uppercase tracking-widest">
                {structureSource === 'pdb' ? 'Modelo Experimental' : structureSource === 'alphafold' ? 'Predicción AlphaFold' : 'Predicción ESMFold'}
              </span>
            </div>
          </div>

          {/* Action Tools */}
          <div className="absolute top-6 right-6 z-10 flex gap-2">
            <button 
              onClick={() => setStructureUrl(null)} 
              className="p-3 bg-slate-900/80 backdrop-blur-xl border border-white/5 rounded-2xl text-slate-400 hover:text-white transition-all pointer-events-auto shadow-xl"
              title="Cambiar Método de Resolución"
            >
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}/></svg>
            </button>
          </div>

          <div style={{ width: '100%', height: '600px' }}>
            <Molstar pdbUrl={structureUrl} />
          </div>
        </div>
      )}
    </div>
  )
}
