/**
 * MolstarProteinViewer Component — Clean Laboratory Edition
 * 3D protein structure visualizer with Mol*
 * Supports tertiary and quaternary structures with refined light aesthetic
 */
import { useState, useEffect, useMemo, useRef } from 'react'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import 'molstar/build/viewer/molstar.css';

// Internal Molstar Component with strict lifecycle management
const Molstar = ({ pdbUrl }) => {
  const parentRef = useRef(null);
  const pluginRef = useRef(null);
  const isInitializing = useRef(false);

  useEffect(() => {
    let isMounted = true;
    
    async function init() {
      // Prevent concurrent initialization
      if (!parentRef.current || isInitializing.current) return;
      isInitializing.current = true;

      // Ensure full cleanup of previous instances
      if (pluginRef.current) {
        try {
          pluginRef.current.dispose();
        } catch (e) {}
        pluginRef.current = null;
      }

      // Absolute clear of container to prevent React 18 root conflicts
      if (parentRef.current) {
        parentRef.current.innerHTML = '';
      }

      const spec = DefaultPluginUISpec();
      spec.layout = {
          initial: {
              isExpanded: false,
              showControls: true,
              controlsDisplay: 'reactive'
          }
      };

      try {
        const plugin = await createPluginUI(parentRef.current, spec);
        
        if (!isMounted) {
          plugin.dispose();
          isInitializing.current = false;
          return;
        }

        pluginRef.current = plugin;

        if (pdbUrl && isMounted) {
          const data = await plugin.builders.data.download(
            { url: pdbUrl, isBinary: false },
            { state: { isGhost: true } }
          );

          const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');

          await plugin.builders.structure.hierarchy.applyPreset(
            trajectory,
            'all-models',
            {
              representationPreset: 'auto',
              showUnitcell: false,
              showWater: false,
              showHet: true,
              ignoreHydrogens: true
            }
          );
        }
      } catch (err) {
        if (isMounted) console.error('Error loading structure in Molstar:', err);
      } finally {
        isInitializing.current = false;
      }
    }

    // Delay init slightly to ensure DOM is ready and previous cleanups finished
    const timer = setTimeout(init, 50);

    return () => {
      isMounted = false;
      clearTimeout(timer);
      if (pluginRef.current) {
        const p = pluginRef.current;
        pluginRef.current = null;
        // Schedule disposal to happen after React's immediate lifecycle
        setTimeout(() => p.dispose(), 0);
      }
    };
  }, [pdbUrl]);

  return <div ref={parentRef} style={{ width: '100%', height: '100%', minHeight: '500px' }} />;
};

const PDB_FILE_API = 'https://files.rcsb.org/download'

export default function MolstarProteinViewer({ proteinId, proteinSequence, productName }) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [structureUrl, setStructureUrl] = useState(null)
  const [structureSource, setStructureSource] = useState(null)
  const [manualPdbId, setManualPdbId] = useState('')

  useEffect(() => {
    return () => {
      if (structureUrl && structureUrl.startsWith('blob:')) {
        URL.revokeObjectURL(structureUrl)
      }
    }
  }, [structureUrl])

  const predictStructure = async () => {
    try {
      setLoading(true)
      setError(null)

      if (!proteinSequence || proteinSequence.length < 10) {
        throw new Error('Secuencia demasiado corta')
      }

      const response = await fetch(`/api/ncbi/protein/predict-structure/${proteinId}`, {
        method: 'POST'
      })

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}))
        throw new Error(errorData.detail || `Error ${response.status}`)
      }

      const pdbContent = await response.text()
      const blob = new Blob([pdbContent], { type: 'chemical/x-pdb' })
      const pdbUrl = URL.createObjectURL(blob)

      setStructureUrl(pdbUrl)
      setStructureSource('esmfold')
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  const loadFromPDB = async (pdbId) => {
    try {
      setLoading(true)
      setError(null)
      const cleanPdbId = pdbId.trim().toUpperCase()
      if (!/^[A-Z0-9]{4}$/i.test(cleanPdbId)) {
        throw new Error(`ID PDB inválido: ${pdbId}`)
      }
      const pdbUrl = `${PDB_FILE_API}/${cleanPdbId}.pdb`
      const response = await fetch(pdbUrl, { method: 'HEAD' })
      if (!response.ok) throw new Error(`PDB ${cleanPdbId} no encontrado`)
      setStructureUrl(pdbUrl)
      setStructureSource('pdb')
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  const handleManualPdbLoad = (e) => {
    e.preventDefault()
    if (manualPdbId.trim()) loadFromPDB(manualPdbId.trim())
  }

  return (
    <div className="space-y-6">
      {/* HUD de Resolución */}
      {!structureUrl && !loading && (
        <div className="bg-white rounded-[2.5rem] border border-slate-200 p-12 text-center shadow-sm animate-in fade-in zoom-in-95 duration-700">
          <div className="w-28 h-28 bg-blue-50 rounded-full flex items-center justify-center mx-auto mb-10 border border-blue-100 relative">
            <span className="text-5xl animate-pulse">🧬</span>
          </div>
          
          <h3 className="text-2xl font-black text-slate-900 uppercase tracking-tighter mb-4">Resolución Estructural</h3>
          <p className="text-slate-500 text-sm max-w-lg mx-auto mb-12 leading-relaxed font-medium">
            Para visualizar la arquitectura 3D, el sistema debe resolver la geometría molecular mediante IA o datos experimentales.
          </p>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-3xl mx-auto">
            <button
              onClick={predictStructure}
              disabled={!proteinSequence || proteinSequence.length > 1000}
              className={`group p-8 rounded-[2rem] border transition-all duration-500 text-left ${
                proteinSequence?.length > 1000 
                  ? 'bg-slate-50 border-slate-100 opacity-50 cursor-not-allowed' 
                  : 'bg-white border-slate-200 hover:border-blue-400 hover:shadow-xl'
              }`}
            >
              <div className="flex items-center gap-4 mb-4">
                <div className="p-3 bg-blue-50 rounded-2xl font-black text-xs text-blue-600">IA</div>
                <h4 className="text-xs font-black text-slate-700 uppercase tracking-widest">ESMFold</h4>
              </div>
              <p className="text-xs text-slate-400 font-medium leading-relaxed">
                Predicción de novo mediante redes neuronales (Hasta 1000 aa).
              </p>
              {proteinSequence?.length > 1000 && (
                <p className="mt-4 text-[10px] font-black text-rose-500 uppercase">Proteína demasiado grande ({proteinSequence.length} aa)</p>
              )}
            </button>

            <div className="p-8 bg-white rounded-[2rem] border border-slate-200 text-left transition-all duration-500 hover:border-indigo-400 hover:shadow-xl">
              <div className="flex items-center gap-4 mb-4">
                <div className="p-3 bg-indigo-50 rounded-2xl font-black text-xs text-indigo-600">DB</div>
                <h4 className="text-xs font-black text-slate-700 uppercase tracking-widest">PDB Import</h4>
              </div>
              <p className="text-xs text-slate-400 font-medium leading-relaxed mb-6">Importación de estructuras experimentales mediante ID oficial.</p>
              <form onSubmit={handleManualPdbLoad} className="flex gap-2">
                <input
                  type="text"
                  value={manualPdbId}
                  onChange={(e) => setManualPdbId(e.target.value)}
                  placeholder="ID PDB..."
                  className="flex-1 bg-slate-50 border border-slate-200 rounded-xl px-4 py-3 text-xs font-mono text-blue-600 focus:outline-none focus:border-blue-400 placeholder-slate-300"
                  maxLength={10}
                />
                <button type="submit" className="px-4 bg-blue-600 text-white rounded-xl shadow-lg active:scale-95 transition-all">
                  <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5}/></svg>
                </button>
              </form>
            </div>
          </div>
        </div>
      )}

      {loading && (
        <div className="h-[650px] bg-slate-50 rounded-[3rem] flex flex-col items-center justify-center border border-slate-200">
          <div className="relative w-20 h-20 mb-8">
            <div className="absolute inset-0 border-4 border-slate-200 rounded-full"></div>
            <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin"></div>
          </div>
          <p className="text-slate-400 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse">Procesando Malla Molecular</p>
        </div>
      )}

      {error && !structureUrl && !loading && (
        <div className="p-8 bg-rose-50 rounded-[2rem] border border-rose-100 text-center max-w-2xl mx-auto">
          <h4 className="text-sm font-black text-rose-600 uppercase tracking-widest mb-2">Error de Resolución</h4>
          <p className="text-xs text-rose-400 font-medium mb-6 leading-relaxed">{error}</p>
          <div className="flex justify-center gap-4">
            <button onClick={() => setError(null)} className="px-6 py-2 bg-white border border-rose-200 rounded-xl text-[10px] text-rose-600 font-black uppercase tracking-widest hover:bg-rose-100 transition-all">Limpiar</button>
            <button onClick={() => window.open(`https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22value%22%3A%22${productName || proteinId}%22%7D%7D%2C%22return_type%22%3A%22entry%22%7D`, '_blank')} className="px-6 py-2 bg-blue-600 text-white rounded-xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-700 transition-all">Buscar en RCSB</button>
          </div>
        </div>
      )}

      {structureUrl && (
        <div className="relative group rounded-[3rem] overflow-hidden border border-slate-200 shadow-2xl bg-slate-900">
          <div className="absolute top-8 left-8 z-10 flex gap-3 pointer-events-none">
            <div className="px-5 py-2.5 bg-white/90 backdrop-blur-xl border border-slate-200 rounded-full flex items-center gap-4 shadow-xl">
              <div className="w-2.5 h-2.5 bg-emerald-500 rounded-full shadow-[0_0_12px_rgba(16,185,129,0.4)] animate-pulse"></div>
              <span className="text-[10px] font-black text-slate-700 uppercase tracking-[0.15em]">
                {structureSource === 'pdb' ? 'Modelo Experimental' : 'Predicción Molecular'}
              </span>
              <div className="w-px h-3 bg-slate-200"></div>
              <span className="text-[9px] font-bold text-slate-400">{proteinSequence?.length} aa</span>
            </div>
          </div>

          <div className="absolute top-8 right-8 z-10 flex gap-3">
            <button 
              onClick={() => setStructureUrl(null)} 
              className="w-12 h-12 bg-white/90 backdrop-blur-xl border border-slate-200 rounded-2xl text-slate-400 hover:text-blue-600 transition-all flex items-center justify-center shadow-xl group/btn"
              title="Reiniciar Resolución"
            >
              <svg className="w-5 h-5 group-hover/btn:rotate-180 transition-transform duration-700" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5}/></svg>
            </button>
          </div>

          <div style={{ width: '100%', height: '700px' }}>
            <Molstar pdbUrl={structureUrl} />
          </div>

          <div className="absolute bottom-8 left-1/2 -translate-x-1/2 px-8 py-3 bg-white/80 backdrop-blur-xl rounded-full border border-slate-200 opacity-0 group-hover:opacity-100 transition-all duration-700 flex gap-12 shadow-2xl pointer-events-none">
            <div className="flex items-center gap-3">
              <span className="text-[10px] font-black text-blue-600 uppercase tracking-widest">Rotar</span>
              <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">LMB</span>
            </div>
            <div className="flex items-center gap-3">
              <span className="text-[10px] font-black text-blue-600 uppercase tracking-widest">Mover</span>
              <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">RMB</span>
            </div>
            <div className="flex items-center gap-3">
              <span className="text-[10px] font-black text-blue-600 uppercase tracking-widest">Zoom</span>
              <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Rueda</span>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
