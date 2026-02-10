/**
 * MolstarProteinViewer Component
 * Visualizador 3D avanzado de estructuras de proteínas usando Mol*
 * Soporta PDB, AlphaFold y estructuras predichas
 */
import { useState, useEffect, useMemo, useRef } from 'react'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/build/viewer/molstar.css'; 

// Componente Wrapper para Molstar
const Molstar = ({ pdbUrl }) => {
  const parentRef = useRef(null);
  const pluginRef = useRef(null);

  useEffect(() => {
    async function init() {
      if (!parentRef.current) return;
      
      const spec = DefaultPluginUISpec();
      spec.layout = {
          initial: {
              isExpanded: false,
              showControls: true
          }
      };

      const plugin = await createPluginUI(parentRef.current, spec);
      pluginRef.current = plugin;

      if (pdbUrl) {
        const data = await plugin.builders.data.download({ url: pdbUrl }, { state: { isGhost: true } });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
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

  return <div ref={parentRef} style={{ width: '100%', height: '100%' }} />;
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

  // Detectar tipo de ID automáticamente
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

  // Intentar cargar estructura desde AlphaFold
  const tryAlphaFold = async (uniprotId) => {
    try {
      setLoading(true)
      setError(null)

      // Buscar en AlphaFold
      const searchUrl = `https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`
      const response = await fetch(searchUrl)

      if (!response.ok) {
        throw new Error('No encontrado en AlphaFold')
      }

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

  // Cargar estructura desde PDB
  const loadFromPDB = async (pdbId) => {
    try {
      setLoading(true)
      setError(null)

      const pdbUrl = `${PDB_FILE_API}/${pdbId}.pdb`

      // Verificar que existe
      const response = await fetch(pdbUrl, { method: 'HEAD' })
      if (!response.ok) {
        throw new Error(`Estructura PDB ${pdbId} no encontrada`)
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

  // Búsqueda en PDB por nombre de proteína
  const searchPDB = async (query) => {
    try {
      setLoading(true)
      setError(null)

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
          sort: [{ sort_by: 'score', direction: 'desc' }],
          pager: { start: 0, rows: 5 }
        }
      }))}`

      const response = await fetch(searchUrl)
      const data = await response.json()

      if (data.result_set && data.result_set.length > 0) {
        const firstResult = data.result_set[0].identifier
        console.log('🔍 Encontrado en PDB:', firstResult)
        await loadFromPDB(firstResult)
        return true
      }

      throw new Error('No se encontraron resultados en PDB')
    } catch (err) {
      setError(err.message)
      console.error('Error buscando en PDB:', err)
      return false
    } finally {
      setLoading(false)
    }
  }

  // Efecto para cargar automáticamente al montar
  useEffect(() => {
    if (!proteinId) return

    const autoLoad = async () => {
      // Si es un ID de PDB válido, cargar directamente
      if (isPdbId) {
        await loadFromPDB(proteinId)
        return
      }

      // Si es un ID de AlphaFold, cargar directamente
      if (isAlphaFoldId) {
        const afId = proteinId.replace('AF-', '').split('-')[0]
        await tryAlphaFold(afId)
        return
      }

      // Si es RefSeq, intentar buscar por nombre de producto
      if (isRefSeqId && productName) {
        console.log(`🧬 ID de RefSeq detectado (${proteinId}), buscando por nombre: ${productName}`)
        await searchPDB(productName)
        return
      }

      // Por defecto, intentar buscar por ID como query
      if (productName) {
        await searchPDB(productName)
      }
    }

    autoLoad()
  }, [proteinId, productName, isPdbId, isAlphaFoldId, isRefSeqId])

  // Manejadores de eventos
  const handleSearch = async (e) => {
    e.preventDefault()
    if (!searchQuery.trim()) return
    await searchPDB(searchQuery.trim())
  }

  const handleManualPdbLoad = async (e) => {
    e.preventDefault()
    if (!manualPdbId.trim()) return
    await loadFromPDB(manualPdbId.trim().toUpperCase())
  }

  return (
    <div className="space-y-4">
      {/* Información de la proteína */}
      <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
        <div className="grid grid-cols-2 gap-3 text-sm">
          <div>
            <span className="text-slate-500 font-medium">ID Proteína:</span>
            <p className="text-slate-800 font-mono text-xs mt-1">{proteinId || 'No disponible'}</p>
          </div>
          <div>
            <span className="text-slate-500 font-medium">Producto:</span>
            <p className="text-slate-800 text-xs mt-1">{productName || 'No disponible'}</p>
          </div>
          <div>
            <span className="text-slate-500 font-medium">Longitud:</span>
            <p className="text-slate-800 font-mono text-xs mt-1">{proteinSequence?.length || 0} aa</p>
          </div>
          {structureSource && (
            <div>
              <span className="text-slate-500 font-medium">Fuente:</span>
              <p className="text-slate-800 text-xs mt-1 capitalize">
                {structureSource === 'pdb' && '🔬 PDB (Experimental)'}
                {structureSource === 'alphafold' && '🤖 AlphaFold (Predicción)'}
              </p>
            </div>
          )}
        </div>
      </div>

      {/* Controles de búsqueda */}
      <div className="space-y-3">
        {/* Búsqueda por nombre */}
        <form onSubmit={handleSearch} className="flex gap-2">
          <input
            type="text"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            placeholder="Buscar en PDB por nombre (ej: 'insulin')"
            className="flex-1 px-3 py-2 border border-slate-300 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
          <button
            type="submit"
            disabled={loading || !searchQuery.trim()}
            className="px-4 py-2 bg-blue-600 text-white rounded-lg text-sm font-medium hover:bg-blue-700 disabled:bg-slate-300 disabled:cursor-not-allowed transition-colors"
          >
            🔍 Buscar
          </button>
        </form>

        {/* Cargar PDB ID directo */}
        <form onSubmit={handleManualPdbLoad} className="flex gap-2">
          <input
            type="text"
            value={manualPdbId}
            onChange={(e) => setManualPdbId(e.target.value)}
            placeholder="Cargar PDB ID directo (ej: 1CRN, 4HHB)"
            className="flex-1 px-3 py-2 border border-slate-300 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-teal-500"
            maxLength={4}
          />
          <button
            type="submit"
            disabled={loading || !manualPdbId.trim()}
            className="px-4 py-2 bg-teal-600 text-white rounded-lg text-sm font-medium hover:bg-teal-700 disabled:bg-slate-300 disabled:cursor-not-allowed transition-colors"
          >
            📥 Cargar
          </button>
        </form>
      </div>

      {/* Estado de carga y errores */}
      {loading && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 text-center">
          <div className="w-8 h-8 border-3 border-blue-200 border-t-blue-600 rounded-full animate-spin mx-auto mb-2"></div>
          <p className="text-blue-700 text-sm font-medium">Cargando estructura 3D...</p>
        </div>
      )}

      {error && !structureUrl && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
          <p className="text-amber-800 text-sm">
            <span className="font-semibold">⚠️ Aviso:</span> {error}
          </p>
          <p className="text-amber-700 text-xs mt-2">
            Intenta buscar por nombre de proteína o cargar un PDB ID directamente.
          </p>
        </div>
      )}

      {/* Visualizador Mol* */}
      {structureUrl ? (
        <div className="border-2 border-slate-300 rounded-xl overflow-hidden bg-slate-900">
          <div className="bg-slate-800 px-4 py-2 border-b border-slate-700">
            <p className="text-slate-200 text-sm font-medium">
              🧬 Visualizador Mol* - {structureSource === 'pdb' ? 'Estructura Experimental (PDB)' : 'Predicción AlphaFold'}
            </p>
          </div>
          <div style={{ width: '100%', height: '600px' }}>
            <Molstar pdbUrl={structureUrl} />
          </div>
          <div className="bg-slate-800 px-4 py-2 border-t border-slate-700">
            <p className="text-slate-400 text-xs">
              URL: <a href={structureUrl} target="_blank" rel="noopener noreferrer" className="text-blue-400 hover:underline font-mono">{structureUrl}</a>
            </p>
          </div>
        </div>
      ) : (
        !loading && (
          <div className="bg-slate-100 border-2 border-dashed border-slate-300 rounded-xl p-12 text-center">
            <div className="text-6xl mb-4">🧬</div>
            <p className="text-slate-700 font-medium mb-2">Visualizador Mol* listo</p>
            <p className="text-slate-500 text-sm mb-4">
              Busca una proteína por nombre o carga un PDB ID directamente para visualizar su estructura 3D.
            </p>
            <div className="text-xs text-slate-400 space-y-1">
              <p>• Estructuras experimentales desde RCSB PDB</p>
              <p>• Predicciones desde AlphaFold Database</p>
              <p>• Controles interactivos de rotación, zoom y estilo</p>
            </div>
          </div>
        )
      )}

      {/* Información adicional */}
      {structureUrl && (
        <div className="bg-teal-50 border border-teal-200 rounded-lg p-4">
          <p className="text-teal-900 font-semibold text-sm mb-2">💡 Controles de Mol*:</p>
          <ul className="text-teal-800 text-xs space-y-1 list-disc list-inside">
            <li><strong>Click izquierdo + arrastrar:</strong> Rotar estructura</li>
            <li><strong>Rueda del ratón:</strong> Zoom in/out</li>
            <li><strong>Click derecho + arrastrar:</strong> Mover (pan)</li>
            <li><strong>Menú superior:</strong> Cambiar estilos de visualización (cartoon, surface, ball-stick)</li>
          </ul>
        </div>
      )}
    </div>
  )
}
