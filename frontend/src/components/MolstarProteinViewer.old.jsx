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

  // Ref para prevenir cargas duplicadas (React StrictMode)
  const hasAttemptedLoad = useRef(false)

  // Cleanup blob URLs
  useEffect(() => {
    return () => {
      if (structureUrl && structureUrl.startsWith('blob:')) {
        URL.revokeObjectURL(structureUrl)
      }
    }
  }, [structureUrl])

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

  // Predecir estructura usando ESMFold (local backend)
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

      // Get PDB content
      const pdbContent = await response.text()

      // Create blob URL for the PDB data
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

  // Intentar cargar estructura desde AlphaFold DB
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

  // Buscar proteína en AlphaFold DB por UniProt ID (extraído de RefSeq)
  const searchAlphaFoldByRefSeq = async () => {
    try {
      setLoading(true)
      setError(null)

      console.log(`🔍 Buscando en AlphaFold DB para RefSeq: ${proteinId}`)

      // Intentar buscar directamente con el RefSeq ID
      // AlphaFold DB indexa algunas proteínas por RefSeq
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

  // Cargar estructura desde PDB
  const loadFromPDB = async (pdbId) => {
    try {
      setLoading(true)
      setError(null)

      // Validar formato de PDB ID (debe ser 4 caracteres alfanuméricos)
      const cleanPdbId = pdbId.trim().toUpperCase()
      if (!/^[A-Z0-9]{4}$/i.test(cleanPdbId)) {
        throw new Error(`ID de PDB inválido: "${pdbId}". Debe tener exactamente 4 caracteres (ej: 1CRN, 4HHB)`)
      }

      const pdbUrl = `${PDB_FILE_API}/${cleanPdbId}.pdb`

      // Verificar que existe
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

  // Búsqueda en PDB por nombre de proteína
  const searchPDB = async (query) => {
    try {
      setLoading(true)
      setError(null)

      // Sanitizar query: eliminar caracteres especiales y limitar longitud
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

  // Efecto para cargar automáticamente al montar
  useEffect(() => {
    if (!proteinId) return

    // Prevenir cargas duplicadas (React StrictMode en desarrollo monta componentes dos veces)
    if (hasAttemptedLoad.current) return
    hasAttemptedLoad.current = true

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

      // Si es RefSeq, intentar estrategias según tamaño de proteína
      if (isRefSeqId && proteinSequence) {
        const seqLength = proteinSequence.length

        // Primero, intentar AlphaFold DB (puede tener proteínas grandes ya predichas)
        console.log(`🧬 ID de RefSeq detectado (${proteinId}, ${seqLength} aa). Buscando en AlphaFold DB...`)
        const foundInAlphaFold = await searchAlphaFoldByRefSeq()
        if (foundInAlphaFold) {
          return
        }

        // Si no está en AlphaFold DB y es <400 aa, predecir con ESMFold
        if (seqLength <= 400) {
          console.log(`🧬 Prediciendo estructura con ESMFold (${seqLength} aa)...`)
          await predictStructure()
          return
        }

        // Proteínas grandes: mostrar mensaje informativo
        console.log(`⚠️ Proteína grande (${seqLength} aa). ESMFold tiene límite de 400 aa.`)
        setError(`Proteína grande (${seqLength} aa). Usa la búsqueda manual en PDB o busca por nombre en AlphaFold DB.`)
        return
      }

      // Por defecto, intentar buscar por nombre de producto si está disponible
      if (productName) {
        console.log(`🔍 Buscando estructura para: ${productName}`)
        await searchPDB(productName)
      }
    }

    autoLoad()

    // Cleanup: resetear flag cuando cambie el proteinId
    return () => {
      hasAttemptedLoad.current = false
    }
  }, [proteinId, productName, proteinSequence, isPdbId, isAlphaFoldId, isRefSeqId])

  // Manejadores de eventos
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
                {structureSource === 'esmfold' && '🧬 ESMFold (Predicción Local)'}
              </p>
            </div>
          )}
        </div>
      </div>

      {/* Controles de búsqueda */}
      <div className="space-y-3">
        {/* Predicción de estructura con ESMFold */}
        {proteinSequence && (
          <div className={`rounded-lg p-4 border-2 ${
            proteinSequence.length > 400
              ? 'bg-gradient-to-r from-amber-50 to-orange-50 border-amber-300'
              : 'bg-gradient-to-r from-purple-50 to-blue-50 border-purple-300'
          }`}>
            <div className="flex items-start justify-between gap-4">
              <div className="flex-1">
                <h3 className={`font-bold text-sm mb-1 ${
                  proteinSequence.length > 400 ? 'text-amber-900' : 'text-purple-900'
                }`}>
                  🧬 Predicción de Estructura 3D con ESMFold
                </h3>
                <p className={`text-xs mb-2 ${
                  proteinSequence.length > 400 ? 'text-amber-700' : 'text-purple-700'
                }`}>
                  {proteinSequence.length > 0 && `Longitud: ${proteinSequence.length} aminoácidos`}
                </p>
                {proteinSequence.length > 400 ? (
                  <div className="bg-amber-100 border border-amber-300 rounded px-2 py-1.5 mt-2">
                    <p className="text-amber-900 text-xs font-medium">
                      ⚠️ Proteína muy grande. ESMFold solo acepta hasta 400 aa.
                    </p>
                    <p className="text-amber-800 text-xs mt-1">
                      Usa las búsquedas manuales abajo o busca en AlphaFold DB (soporta proteínas grandes).
                    </p>
                  </div>
                ) : proteinSequence.length < 10 ? (
                  <p className="text-red-700 text-xs font-medium">
                    ⚠️ Proteína demasiado corta (mínimo 10 aa requeridos)
                  </p>
                ) : (
                  <p className="text-green-700 text-xs font-medium">
                    ✓ Tamaño compatible con ESMFold
                  </p>
                )}
              </div>
              <button
                onClick={predictStructure}
                disabled={loading || !proteinSequence || proteinSequence.length > 400 || proteinSequence.length < 10}
                className="px-5 py-2.5 bg-gradient-to-r from-purple-600 to-blue-600 text-white rounded-lg text-sm font-bold hover:from-purple-700 hover:to-blue-700 disabled:from-slate-300 disabled:to-slate-400 disabled:cursor-not-allowed transition-all shadow-md hover:shadow-lg whitespace-nowrap"
              >
                {loading ? '⏳ Prediciendo...' : '🚀 Predecir'}
              </button>
            </div>
          </div>
        )}

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
          <div className="mt-3 space-y-2 text-xs text-amber-700">
            <p className="font-medium">Alternativas para proteínas grandes:</p>
            <ul className="list-disc list-inside space-y-1 ml-2">
              <li>Busca el nombre de la proteína en el campo de búsqueda (puede estar en PDB)</li>
              <li>Busca en <a href={`https://alphafold.ebi.ac.uk/search/text/${productName || proteinId}`} target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline font-medium">AlphaFold DB</a> y copia el PDB ID</li>
              <li>Busca en <a href={`https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22value%22%3A%22${productName || proteinId}%22%7D%7D%7D`} target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline font-medium">RCSB PDB</a></li>
            </ul>
          </div>
        </div>
      )}

      {/* Visualizador Mol* */}
      {structureUrl ? (
        <div className="border-2 border-slate-300 rounded-xl overflow-hidden bg-slate-900">
          <div className="bg-slate-800 px-4 py-2 border-b border-slate-700">
            <p className="text-slate-200 text-sm font-medium">
              🧬 Visualizador Mol* - {
                structureSource === 'pdb' ? 'Estructura Experimental (PDB)' :
                structureSource === 'alphafold' ? 'Predicción AlphaFold' :
                'Predicción ESMFold'
              }
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
              <p>• Predicción local con ESMFold AI (recomendado)</p>
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
