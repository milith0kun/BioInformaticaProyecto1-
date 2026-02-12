/**
 * ProteinStructureViewer — Production Molecular Analysis Suite
 * Uses Mol* for ALL 3D visualization with different representation presets.
 * Backend provides real PDB-parsed secondary structure (HELIX/SHEET records)
 * and chain analysis for quaternary structure.
 */
import { useState, useEffect, useRef } from 'react'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18'
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec'
import 'molstar/build/viewer/molstar.css'
import { api } from '../services/api'
import ProteinPipeline from './ProteinPipeline'

// ─── Amino acid properties ───
const AA_PROPERTIES = {
    A: { name: 'Alanina', group: 'Hidrofóbico', color: '#10b981' },
    I: { name: 'Isoleucina', group: 'Hidrofóbico', color: '#10b981' },
    L: { name: 'Leucina', group: 'Hidrofóbico', color: '#10b981' },
    M: { name: 'Metionina', group: 'Hidrofóbico', color: '#10b981' },
    F: { name: 'Fenilalanina', group: 'Hidrofóbico', color: '#10b981' },
    W: { name: 'Triptófano', group: 'Hidrofóbico', color: '#10b981' },
    V: { name: 'Valina', group: 'Hidrofóbico', color: '#10b981' },
    P: { name: 'Prolina', group: 'Especial', color: '#f59e0b' },
    G: { name: 'Glicina', group: 'Especial', color: '#f59e0b' },
    C: { name: 'Cisteína', group: 'Especial', color: '#f59e0b' },
    D: { name: 'Ác. aspártico', group: 'Cargado (−)', color: '#ef4444' },
    E: { name: 'Ác. glutámico', group: 'Cargado (−)', color: '#ef4444' },
    K: { name: 'Lisina', group: 'Cargado (+)', color: '#3b82f6' },
    R: { name: 'Arginina', group: 'Cargado (+)', color: '#3b82f6' },
    H: { name: 'Histidina', group: 'Cargado (+)', color: '#3b82f6' },
    N: { name: 'Asparagina', group: 'Polar', color: '#8b5cf6' },
    Q: { name: 'Glutamina', group: 'Polar', color: '#8b5cf6' },
    S: { name: 'Serina', group: 'Polar', color: '#8b5cf6' },
    T: { name: 'Treonina', group: 'Polar', color: '#8b5cf6' },
    Y: { name: 'Tirosina', group: 'Polar', color: '#8b5cf6' },
    '*': { name: 'Stop', group: 'Stop', color: '#ef4444' },
}

// ─── Mol* Viewer with Preset Support ───
function MolstarViewer({ pdbUrl, preset, height = 600 }) {
    const parentRef = useRef(null)
    const pluginRef = useRef(null)
    const isInitializing = useRef(false)

    useEffect(() => {
        let isMounted = true

        async function init() {
            if (!parentRef.current || isInitializing.current || !pdbUrl) return
            isInitializing.current = true

            if (pluginRef.current) {
                try { pluginRef.current.dispose() } catch (e) { }
                pluginRef.current = null
            }
            if (parentRef.current) parentRef.current.innerHTML = ''

            const spec = DefaultPluginUISpec()
            spec.layout = {
                initial: {
                    isExpanded: false,
                    showControls: true,
                    controlsDisplay: 'reactive'
                }
            }

            try {
                const plugin = await createPluginUI(parentRef.current, spec)
                if (!isMounted) { plugin.dispose(); isInitializing.current = false; return }
                pluginRef.current = plugin

                if (pdbUrl && isMounted) {
                    const isBinary = pdbUrl.toLowerCase().endsWith('.bcif')
                    const data = await plugin.builders.data.download(
                        { url: pdbUrl, isBinary },
                        { state: { isGhost: true } }
                    )

                    const format = isBinary ? 'mmcif' : 'pdb'
                    const trajectory = await plugin.builders.structure.parseTrajectory(data, format)

                    if (preset === 'cartoon' || preset === 'surface') {
                        // Manual representation for secondary/quaternary views
                        const model = await plugin.builders.structure.createModel(trajectory)
                        const structure = await plugin.builders.structure.createStructure(model)

                        // Create polymer component
                        const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer')

                        if (polymer) {
                            if (preset === 'cartoon') {
                                // SECONDARY STRUCTURE: Cartoon colored by SS type
                                // (helices=red, sheets=blue, coils=gray — standard bioinformatics view)
                                await plugin.builders.structure.representation.addRepresentation(
                                    polymer,
                                    {
                                        type: 'cartoon',
                                        color: 'secondary-structure-type',
                                        size: 'uniform'
                                    }
                                )
                            } else if (preset === 'surface') {
                                // QUATERNARY STRUCTURE: Cartoon + transparent surface, colored by chain
                                await plugin.builders.structure.representation.addRepresentation(
                                    polymer,
                                    {
                                        type: 'cartoon',
                                        color: 'chain-id',
                                        size: 'uniform'
                                    }
                                )
                                await plugin.builders.structure.representation.addRepresentation(
                                    polymer,
                                    {
                                        type: 'gaussian-surface',
                                        color: 'chain-id',
                                        typeParams: { alpha: 0.25 }
                                    }
                                )
                            }
                        }

                        // Also show ligands/het atoms if any
                        const ligand = await plugin.builders.structure.tryCreateComponentStatic(structure, 'ligand')
                        if (ligand) {
                            await plugin.builders.structure.representation.addRepresentation(
                                ligand,
                                { type: 'ball-and-stick', color: 'element-symbol' }
                            )
                        }
                    } else {
                        // AUTO / DEFAULT: Standard Mol* preset
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
                        )
                    }
                }
            } catch (err) {
                if (isMounted) console.error('Error loading structure in Molstar:', err)
            } finally {
                isInitializing.current = false
            }
        }

        const timer = setTimeout(init, 80)
        return () => {
            isMounted = false
            clearTimeout(timer)
            if (pluginRef.current) {
                const p = pluginRef.current
                pluginRef.current = null
                setTimeout(() => { try { p.dispose() } catch (e) { } }, 0)
            }
        }
    }, [pdbUrl, preset])

    return <div ref={parentRef} style={{ width: '100%', height, minHeight: 400 }} />
}

// ─── Structure Info Badge ───
function SourceBadge({ source, seqLength }) {
    const labels = {
        pdb: 'Experimental (RCSB PDB)',
        alphafold: 'AlphaFold DB',
        esmfold: 'ESMFold (Predicción IA)',
        unknown: 'Estructura Resuelta'
    }
    return (
        <div className="flex items-center gap-3 px-5 py-2.5 bg-white/90 backdrop-blur-xl border border-slate-200 rounded-full shadow-xl">
            <div className="w-2.5 h-2.5 bg-emerald-500 rounded-full shadow-[0_0_12px_rgba(16,185,129,0.4)] animate-pulse" />
            <span className="text-[10px] font-black text-slate-700 uppercase tracking-[0.15em]">
                {labels[source] || labels.unknown}
            </span>
            <div className="w-px h-3 bg-slate-200" />
            <span className="text-[9px] font-bold text-slate-400">{seqLength} aa</span>
        </div>
    )
}

// ─── Main Component ───
export default function ProteinStructureViewer({
    protein,
    onClose,
    onNavigate,
    currentIndex,
    totalProteins,
    loading: detailLoading
}) {
    const [activeTab, setActiveTab] = useState('overview')
    const [structureUrl, setStructureUrl] = useState(null)
    const [structureSource, setStructureSource] = useState(null)
    const [structureAnalysis, setStructureAnalysis] = useState(null)
    const [pipelineData, setPipelineData] = useState(null)
    const [resolving, setResolving] = useState(false)
    const [error, setError] = useState(null)

    // Keyboard navigation
    useEffect(() => {
        const handleKeyDown = (e) => {
            if (e.key === 'Escape') onClose()
            if (e.key === 'ArrowLeft') onNavigate?.('prev')
            if (e.key === 'ArrowRight') onNavigate?.('next')
        }
        window.addEventListener('keydown', handleKeyDown)
        return () => window.removeEventListener('keydown', handleKeyDown)
    }, [onClose, onNavigate])

    // Reset when protein changes
    useEffect(() => {
        setStructureUrl(null)
        setStructureSource(null)
        setStructureAnalysis(null)
        setPipelineData(null)
        setError(null)
        setActiveTab('overview')
    }, [protein?.protein_id])

    // Cleanup blob URLs
    useEffect(() => {
        return () => {
            if (structureUrl?.startsWith('blob:')) URL.revokeObjectURL(structureUrl)
        }
    }, [structureUrl])

    // ── Resolve structure ──
    const resolveStructure = async () => {
        if (!protein?.protein_id) return
        setResolving(true)
        setError(null)

        try {
            const response = await fetch(`/api/ncbi/protein/predict-structure/${protein.protein_id}`, {
                method: 'POST'
            })
            if (!response.ok) {
                const errData = await response.json().catch(() => ({}))
                throw new Error(errData.detail || `Error ${response.status}`)
            }
            const pdbContent = await response.text()
            const blob = new Blob([pdbContent], { type: 'chemical/x-pdb' })
            const url = URL.createObjectURL(blob)
            const source = response.headers.get('X-Structure-Source') || 'esmfold'
            setStructureUrl(url)
            setStructureSource(source)

            // Fetch full structural analysis and pipeline data in background
            try {
                const [analysis, pipeline] = await Promise.allSettled([
                    api.getProteinStructureAnalysis(protein.protein_id),
                    api.getProteinStructurePipeline(protein.protein_id)
                ])
                if (analysis.status === 'fulfilled') setStructureAnalysis(analysis.value)
                if (pipeline.status === 'fulfilled') setPipelineData(pipeline.value)
            } catch (analysisErr) {
                console.warn('Structure analysis partial:', analysisErr)
            }
        } catch (err) {
            setError(err.message)
        } finally {
            setResolving(false)
        }
    }

    // ── Amino acid composition ──
    const composition = () => {
        if (!protein?.full_sequence) return null
        const seq = protein.full_sequence
        const counts = {}
        for (const aa of seq) counts[aa] = (counts[aa] || 0) + 1

        const groups = { 'Hidrofóbico': 0, 'Polar': 0, 'Cargado (+)': 0, 'Cargado (−)': 0, 'Especial': 0 }
        Object.entries(counts).forEach(([aa, count]) => {
            const prop = AA_PROPERTIES[aa]
            if (prop && groups[prop.group] !== undefined) groups[prop.group] += count
        })

        return {
            total: seq.length,
            groups: Object.entries(groups).map(([group, count]) => ({
                group, count, pct: ((count / seq.length) * 100).toFixed(1)
            })),
            topAA: Object.entries(counts)
                .sort((a, b) => b[1] - a[1])
                .slice(0, 5)
                .map(([aa, count]) => ({
                    aa, count, pct: ((count / seq.length) * 100).toFixed(1),
                    ...(AA_PROPERTIES[aa] || { name: aa, group: 'Otro', color: '#94a3b8' })
                }))
        }
    }
    const aaComp = composition()

    if (!protein) return null

    const ss = structureAnalysis?.secondary_structure
    const quat = structureAnalysis?.quaternary
    const quality = structureAnalysis?.quality

    // Tabs definition
    const tabs = [
        { id: 'overview', label: 'Estructura 3D' },
        ...(structureUrl ? [
            { id: 'pipeline', label: '🧬 Pipeline Completo' },
            { id: 'secondary', label: `Secundaria${ss ? ` (${ss.stats.helix_count}H/${ss.stats.sheet_count}E)` : ''}` },
            { id: 'quaternary', label: `Cuaternaria${quat ? ` (${quat.oligomeric_state})` : ''}` },
        ] : []),
        { id: 'primary', label: 'Secuencia' },
        { id: 'composition', label: 'Residuos' },
    ]

    return (
        <div
            className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/40 backdrop-blur-sm p-4 overflow-hidden"
            onClick={(e) => e.target === e.currentTarget && onClose()}
        >
            <div className="relative w-full h-full max-w-7xl max-h-screen bg-[#f8fafc] rounded-3xl shadow-2xl overflow-hidden border border-slate-200 flex flex-col" onClick={e => e.stopPropagation()}>

                {/* Header */}
                <div className="relative bg-white/80 backdrop-blur-2xl border-b border-slate-200 px-10 py-6 flex-shrink-0">
                    <div className="flex items-center justify-between gap-8">
                        <button onClick={onClose} className="flex items-center gap-3 px-5 py-2.5 bg-slate-100 hover:bg-slate-200 text-slate-600 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 group shadow-sm">
                            <svg className="w-5 h-5 group-hover:rotate-90 transition-transform duration-300" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" /></svg>
                            <span className="text-xs font-black uppercase tracking-widest">Cerrar</span>
                        </button>

                        <div className="flex-1 min-w-0 text-center space-y-1">
                            <h3 className="font-mono text-xl font-black text-blue-600 tracking-tighter truncate uppercase">
                                {protein.protein_id || protein.locus_tag}
                            </h3>
                            <p className="text-xs font-bold text-slate-400 uppercase tracking-[0.2em] truncate">{protein.product || 'Proteína hipotética'}</p>
                        </div>

                        {onNavigate && (
                            <div className="flex items-center gap-3">
                                <button onClick={() => onNavigate('prev')} disabled={currentIndex === 0} className="p-3 bg-slate-100 hover:bg-slate-200 text-slate-500 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 disabled:opacity-20 shadow-sm">
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M15 19l-7-7 7-7" /></svg>
                                </button>
                                <div className="px-5 py-2.5 bg-white rounded-2xl border border-slate-200 text-xs font-black font-mono text-blue-600 shadow-inner">
                                    {currentIndex + 1} <span className="text-slate-300 px-1">/</span> {totalProteins}
                                </div>
                                <button onClick={() => onNavigate('next')} disabled={currentIndex === totalProteins - 1} className="p-3 bg-slate-100 hover:bg-slate-200 text-slate-500 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 disabled:opacity-20 shadow-sm">
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M9 5l7 7-7 7" /></svg>
                                </button>
                            </div>
                        )}
                    </div>

                    <div className="flex items-center justify-center gap-8 mt-6">
                        <div className="flex items-center gap-2 px-4 py-1.5 bg-blue-50 rounded-full border border-blue-100">
                            <div className="w-1.5 h-1.5 rounded-full bg-blue-500" />
                            <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest"><span className="text-blue-600">{protein.length}</span> Residuos</span>
                        </div>
                        <div className="flex items-center gap-2 px-4 py-1.5 bg-indigo-50 rounded-full border border-indigo-100">
                            <div className="w-1.5 h-1.5 rounded-full bg-indigo-500" />
                            <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest"><span className="text-indigo-600">{((protein.molecular_weight_approx || 0) / 1000).toFixed(1)}</span> kDa</span>
                        </div>
                        {protein.gene_name && (
                            <div className="flex items-center gap-2 px-4 py-1.5 bg-emerald-50 rounded-full border border-emerald-100">
                                <div className="w-1.5 h-1.5 rounded-full bg-emerald-500" />
                                <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest">Gen: <span className="text-emerald-600">{protein.gene_name}</span></span>
                            </div>
                        )}
                    </div>
                </div>

                {/* Tabs */}
                <div className="bg-slate-50/50 backdrop-blur-2xl border-b border-slate-200 px-10 flex-shrink-0">
                    <div className="flex gap-10 overflow-x-auto scrollbar-hide">
                        {tabs.map(tab => (
                            <button
                                key={tab.id}
                                onClick={() => setActiveTab(tab.id)}
                                className={`group relative py-6 text-[10px] font-black uppercase tracking-[0.25em] transition-all whitespace-nowrap ${activeTab === tab.id ? 'text-blue-600' : 'text-slate-500 hover:text-slate-900'
                                    }`}
                            >
                                <span>{tab.label}</span>
                                {activeTab === tab.id && (
                                    <div className="absolute bottom-0 left-0 right-0 h-0.5 bg-blue-600 shadow-[0_0_12px_rgba(59,130,246,0.3)]" />
                                )}
                            </button>
                        ))}
                    </div>
                </div>

                {/* Content */}
                <div className="flex-1 overflow-y-auto p-8 lg:p-12 space-y-12 custom-scrollbar min-h-0">
                    {detailLoading ? (
                        <div className="flex flex-col items-center justify-center py-48">
                            <div className="relative w-24 h-24 mb-8">
                                <div className="absolute inset-0 border-4 border-slate-200 rounded-full" />
                                <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin" />
                            </div>
                            <p className="text-slate-400 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse text-center">Cargando datos moleculares</p>
                        </div>
                    ) : (
                        <>
                            {/* ══════ TAB: 3D Structure ══════ */}
                            {activeTab === 'overview' && (
                                <div className="space-y-10 animate-in fade-in slide-in-from-bottom-6 duration-1000">
                                    {!structureUrl && !resolving && (
                                        <div className="bg-white rounded-[2.5rem] border border-slate-200 p-12 text-center shadow-sm">
                                            <div className="w-28 h-28 bg-blue-50 rounded-full flex items-center justify-center mx-auto mb-10 border border-blue-100">
                                                <span className="text-5xl animate-pulse">🧬</span>
                                            </div>
                                            <h3 className="text-2xl font-black text-slate-900 uppercase tracking-tighter mb-4">Resolución Estructural</h3>
                                            <p className="text-slate-500 text-sm max-w-lg mx-auto mb-12 leading-relaxed font-medium">
                                                El sistema buscará en <strong>AlphaFold DB</strong>, <strong>RCSB PDB</strong> (experimental) y <strong>ESMFold</strong> (predicción IA) automáticamente.
                                            </p>
                                            <button
                                                onClick={resolveStructure}
                                                disabled={!protein?.full_sequence}
                                                className="px-10 py-4 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-600 transition-all active:scale-95 shadow-2xl disabled:opacity-30"
                                            >
                                                ▶ Resolver Estructura 3D
                                            </button>
                                            {protein?.full_sequence?.length > 800 && (
                                                <p className="mt-4 text-[10px] font-black text-amber-500 uppercase">Proteína grande ({protein.full_sequence.length} aa): Puede demorar</p>
                                            )}
                                        </div>
                                    )}

                                    {resolving && (
                                        <div className="h-[500px] bg-slate-50 rounded-[3rem] flex flex-col items-center justify-center border border-slate-200">
                                            <div className="relative w-20 h-20 mb-8">
                                                <div className="absolute inset-0 border-4 border-slate-200 rounded-full" />
                                                <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin" />
                                            </div>
                                            <p className="text-slate-400 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse">Buscando en AlphaFold / PDB / ESMFold...</p>
                                        </div>
                                    )}

                                    {error && (
                                        <div className="p-8 bg-rose-50 rounded-[2rem] border border-rose-100 text-center max-w-2xl mx-auto">
                                            <h4 className="text-sm font-black text-rose-600 uppercase tracking-widest mb-2">Error de Resolución</h4>
                                            <p className="text-xs text-rose-400 font-medium mb-6 leading-relaxed">{error}</p>
                                            <div className="flex justify-center gap-4">
                                                <button onClick={() => setError(null)} className="px-6 py-2 bg-white border border-rose-200 rounded-xl text-[10px] text-rose-600 font-black uppercase">Reintentar</button>
                                                <button onClick={() => window.open(`https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22value%22%3A%22${protein.product || protein.protein_id}%22%7D%7D%2C%22return_type%22%3A%22entry%22%7D`, '_blank')} className="px-6 py-2 bg-blue-600 text-white rounded-xl text-[10px] font-black uppercase">Buscar en RCSB</button>
                                            </div>
                                        </div>
                                    )}

                                    {structureUrl && (
                                        <>
                                            <div className="relative group rounded-[3rem] overflow-hidden border border-slate-200 shadow-2xl bg-slate-900">
                                                <div className="absolute top-8 left-8 z-10 pointer-events-none">
                                                    <SourceBadge source={structureSource} seqLength={protein?.full_sequence?.length || protein.length} />
                                                </div>
                                                <div className="absolute top-8 right-8 z-10">
                                                    <button onClick={() => { setStructureUrl(null); setStructureAnalysis(null) }} className="w-12 h-12 bg-white/90 border border-slate-200 rounded-2xl text-slate-400 hover:text-blue-600 transition-all flex items-center justify-center shadow-xl" title="Reiniciar">
                                                        <svg className="w-5 h-5 hover:rotate-180 transition-transform duration-700" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} /></svg>
                                                    </button>
                                                </div>
                                                <MolstarViewer pdbUrl={structureUrl} preset="auto" height={700} />
                                            </div>

                                            {/* Quick stats based on real structural analysis */}
                                            {structureAnalysis && (
                                                <div className="grid grid-cols-1 md:grid-cols-4 gap-6">
                                                    <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                                        <div className="flex items-center gap-4 mb-4">
                                                            <div className="w-10 h-10 bg-rose-50 rounded-2xl flex items-center justify-center border border-rose-100 text-rose-600 font-black text-xs">2°</div>
                                                            <h5 className="text-[10px] font-black text-slate-700 uppercase tracking-[0.2em]">Secundaria</h5>
                                                        </div>
                                                        <p className="text-[11px] text-slate-500 leading-relaxed font-medium">
                                                            {ss?.stats.helix_count || 0} hélices α ({ss?.stats.helix_pct || 0}%) · {ss?.stats.sheet_count || 0} láminas β ({ss?.stats.sheet_pct || 0}%)
                                                        </p>
                                                    </div>
                                                    <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                                        <div className="flex items-center gap-4 mb-4">
                                                            <div className="w-10 h-10 bg-blue-50 rounded-2xl flex items-center justify-center border border-blue-100 text-blue-600 font-black text-xs">3°</div>
                                                            <h5 className="text-[10px] font-black text-slate-700 uppercase tracking-[0.2em]">Terciaria</h5>
                                                        </div>
                                                        <p className="text-[11px] text-slate-500 leading-relaxed font-medium">
                                                            {structureAnalysis.total_atoms?.toLocaleString()} átomos · Fuente: {structureSource}
                                                        </p>
                                                    </div>
                                                    <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                                        <div className="flex items-center gap-4 mb-4">
                                                            <div className="w-10 h-10 bg-indigo-50 rounded-2xl flex items-center justify-center border border-indigo-100 text-indigo-600 font-black text-xs">4°</div>
                                                            <h5 className="text-[10px] font-black text-slate-700 uppercase tracking-[0.2em]">Cuaternaria</h5>
                                                        </div>
                                                        <p className="text-[11px] text-slate-500 leading-relaxed font-medium">
                                                            {quat?.oligomeric_state} · {quat?.num_chains} cadena{quat?.num_chains > 1 ? 's' : ''}
                                                        </p>
                                                    </div>
                                                    <div className="bg-gradient-to-br from-blue-600 to-indigo-700 rounded-[2rem] p-8 text-white shadow-xl relative overflow-hidden group">
                                                        <div className="absolute top-0 right-0 w-32 h-32 bg-white/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000" />
                                                        <p className="text-[10px] font-black uppercase tracking-widest mb-3 opacity-70 relative z-10">
                                                            {quality?.is_plddt ? 'pLDDT' : 'Resolución'}
                                                        </p>
                                                        <p className="text-xs font-bold leading-relaxed relative z-10">
                                                            {quality?.is_plddt
                                                                ? 'Los colores B-factor representan confianza pLDDT de la predicción IA.'
                                                                : structureAnalysis.resolution
                                                                    ? `${structureAnalysis.resolution} Å de resolución experimental.`
                                                                    : 'Use los controles de Mol* para explorar la geometría molecular.'}
                                                        </p>
                                                    </div>
                                                </div>
                                            )}
                                        </>
                                    )}
                                </div>
                            )}

                            {/* ══════ TAB: Structure Pipeline (Primary→Quaternary) ══════ */}
                            {activeTab === 'pipeline' && pipelineData && (
                                <div className="animate-in fade-in slide-in-from-bottom-6 duration-700">
                                    <ProteinPipeline
                                        pipelineData={pipelineData}
                                        pdbUrl={structureUrl}
                                        structureSource={structureSource}
                                    />
                                </div>
                            )}

                            {/* ══════ TAB: Secondary Structure (Mol* cartoon + real data) ══════ */}
                            {activeTab === 'secondary' && structureUrl && (
                                <div className="space-y-10 animate-in fade-in slide-in-from-right-8 duration-700">
                                    {/* Mol* in cartoon mode for secondary structure */}
                                    <div className="rounded-[3rem] overflow-hidden border border-slate-200 shadow-2xl bg-slate-900 relative">
                                        <div className="absolute top-6 left-6 z-10 px-4 py-2 bg-white/90 rounded-full border border-slate-200 shadow-lg pointer-events-none">
                                            <span className="text-[10px] font-black text-slate-700 uppercase tracking-widest">Representación Cartoon — Estructura Secundaria</span>
                                        </div>
                                        <MolstarViewer pdbUrl={structureUrl} preset="cartoon" height={500} />
                                    </div>

                                    {/* Real HELIX/SHEET data from PDB */}
                                    {ss && (
                                        <>
                                            <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                                                {[
                                                    { label: 'Hélices α', count: ss.stats.helix_count, residues: ss.stats.helix_residues, pct: ss.stats.helix_pct, color: '#e11d48', bg: '#fff1f2', border: '#fecdd3' },
                                                    { label: 'Láminas β', count: ss.stats.sheet_count, residues: ss.stats.sheet_residues, pct: ss.stats.sheet_pct, color: '#2563eb', bg: '#eff6ff', border: '#bfdbfe' },
                                                    { label: 'Ovillos / Coil', count: '-', residues: ss.stats.coil_residues, pct: ss.stats.coil_pct, color: '#64748b', bg: '#f8fafc', border: '#e2e8f0' },
                                                ].map(s => (
                                                    <div key={s.label} className="rounded-3xl p-8 border text-center shadow-sm" style={{ backgroundColor: s.bg, borderColor: s.border }}>
                                                        <p className="text-[10px] font-black uppercase tracking-[0.2em] mb-3" style={{ color: s.color }}>{s.label}</p>
                                                        <p className="text-5xl font-black mb-2" style={{ color: s.color }}>{s.pct}%</p>
                                                        <p className="text-[10px] text-slate-400 font-bold">{s.residues} residuos{s.count !== '-' ? ` · ${s.count} segmentos` : ''}</p>
                                                    </div>
                                                ))}
                                            </div>

                                            {/* Helices table */}
                                            {ss.helices.length > 0 && (
                                                <div className="bg-white rounded-2xl border border-slate-200 overflow-hidden shadow-sm">
                                                    <div className="px-6 py-4 border-b border-slate-100">
                                                        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Hélices α del PDB ({ss.helices.length})</p>
                                                    </div>
                                                    <div className="max-h-60 overflow-y-auto">
                                                        <table className="w-full text-xs">
                                                            <thead className="bg-slate-50 sticky top-0">
                                                                <tr>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">ID</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Cadena</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Inicio</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Fin</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Largo</th>
                                                                </tr>
                                                            </thead>
                                                            <tbody>
                                                                {ss.helices.map((h, i) => (
                                                                    <tr key={i} className="border-t border-slate-50 hover:bg-rose-50/50">
                                                                        <td className="px-4 py-2 font-mono font-bold text-rose-600">{h.id}</td>
                                                                        <td className="px-4 py-2 font-bold text-slate-600">{h.init_chain || 'A'}</td>
                                                                        <td className="px-4 py-2 font-mono">{h.init_res_name} {h.init_seq_num}</td>
                                                                        <td className="px-4 py-2 font-mono">{h.end_res_name} {h.end_seq_num}</td>
                                                                        <td className="px-4 py-2 font-mono font-bold text-rose-600">{h.length || (h.end_seq_num - h.init_seq_num + 1)}</td>
                                                                    </tr>
                                                                ))}
                                                            </tbody>
                                                        </table>
                                                    </div>
                                                </div>
                                            )}

                                            {/* Sheets table */}
                                            {ss.sheets.length > 0 && (
                                                <div className="bg-white rounded-2xl border border-slate-200 overflow-hidden shadow-sm">
                                                    <div className="px-6 py-4 border-b border-slate-100">
                                                        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Láminas β del PDB ({ss.sheets.length} hebras)</p>
                                                    </div>
                                                    <div className="max-h-60 overflow-y-auto">
                                                        <table className="w-full text-xs">
                                                            <thead className="bg-slate-50 sticky top-0">
                                                                <tr>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Lámina</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Hebra</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Cadena</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Inicio</th>
                                                                    <th className="px-4 py-2 text-left text-[9px] font-black text-slate-400 uppercase">Fin</th>
                                                                </tr>
                                                            </thead>
                                                            <tbody>
                                                                {ss.sheets.map((s, i) => (
                                                                    <tr key={i} className="border-t border-slate-50 hover:bg-blue-50/50">
                                                                        <td className="px-4 py-2 font-mono font-bold text-blue-600">{s.sheet_id}</td>
                                                                        <td className="px-4 py-2 font-bold text-slate-600">{s.strand}</td>
                                                                        <td className="px-4 py-2 font-bold text-slate-600">{s.init_chain || 'A'}</td>
                                                                        <td className="px-4 py-2 font-mono">{s.init_res_name} {s.init_seq_num}</td>
                                                                        <td className="px-4 py-2 font-mono">{s.end_res_name} {s.end_seq_num}</td>
                                                                    </tr>
                                                                ))}
                                                            </tbody>
                                                        </table>
                                                    </div>
                                                </div>
                                            )}
                                        </>
                                    )}

                                    {!ss && (
                                        <div className="bg-slate-50 rounded-3xl p-12 text-center border border-slate-100">
                                            <div className="w-8 h-8 border-4 border-slate-200 border-t-blue-500 rounded-full animate-spin mx-auto mb-4" />
                                            <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Analizando registros HELIX/SHEET del PDB...</p>
                                        </div>
                                    )}
                                </div>
                            )}

                            {/* ══════ TAB: Quaternary Structure ══════ */}
                            {activeTab === 'quaternary' && structureUrl && (
                                <div className="space-y-10 animate-in fade-in slide-in-from-bottom-8 duration-700">
                                    {/* Mol* with surface/chain coloring */}
                                    <div className="rounded-[3rem] overflow-hidden border border-slate-200 shadow-2xl bg-slate-900 relative">
                                        <div className="absolute top-6 left-6 z-10 px-4 py-2 bg-white/90 rounded-full border border-slate-200 shadow-lg pointer-events-none">
                                            <span className="text-[10px] font-black text-slate-700 uppercase tracking-widest">
                                                Vista por Cadenas — {quat?.oligomeric_state || 'Analizando...'}
                                            </span>
                                        </div>
                                        <MolstarViewer pdbUrl={structureUrl} preset="surface" height={500} />
                                    </div>

                                    {quat && (
                                        <>
                                            {/* Oligomeric state card */}
                                            <div className="bg-gradient-to-br from-slate-900 via-slate-800 to-blue-900 rounded-3xl p-10 text-white relative overflow-hidden shadow-2xl">
                                                <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[80px] rounded-full" />
                                                <div className="relative z-10 flex items-start justify-between">
                                                    <div className="space-y-3">
                                                        <p className="text-[10px] font-black uppercase tracking-[0.3em] text-blue-300/70">Estado Oligomérico</p>
                                                        <h3 className="text-4xl font-black uppercase tracking-tighter">{quat.oligomeric_state}</h3>
                                                        <p className="text-sm text-slate-300 max-w-md leading-relaxed">
                                                            {quat.is_complex
                                                                ? `Complejo multiproteico con ${quat.num_chains} cadenas distintas detectadas en el archivo PDB.`
                                                                : 'Estructura monomérica — la proteína funciona como unidad individual en esta resolución.'}
                                                        </p>
                                                    </div>
                                                    <div className="text-right space-y-1">
                                                        <p className="text-6xl font-black text-blue-400/50">{quat.num_chains}</p>
                                                        <p className="text-[9px] font-black uppercase tracking-widest text-slate-400">cadena{quat.num_chains > 1 ? 's' : ''}</p>
                                                    </div>
                                                </div>
                                            </div>

                                            {/* Chain details */}
                                            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
                                                {quat.chains.map((chain, i) => {
                                                    const chainColors = ['#3b82f6', '#e11d48', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899']
                                                    const color = chainColors[i % chainColors.length]
                                                    return (
                                                        <div key={chain.chain_id} className="bg-white rounded-2xl border border-slate-200 p-6 shadow-sm hover:shadow-lg transition-all">
                                                            <div className="flex items-center gap-4 mb-4">
                                                                <div className="w-10 h-10 rounded-2xl flex items-center justify-center font-black text-white text-sm" style={{ backgroundColor: color }}>
                                                                    {chain.chain_id}
                                                                </div>
                                                                <div>
                                                                    <p className="text-sm font-black text-slate-900">Cadena {chain.chain_id}</p>
                                                                    <p className="text-[10px] text-slate-400 font-bold">{chain.residue_count} residuos · {chain.atom_count.toLocaleString()} átomos</p>
                                                                </div>
                                                            </div>
                                                            <div className="space-y-2">
                                                                <div className="flex justify-between text-[10px] font-bold">
                                                                    <span className="text-slate-400 uppercase">Rango</span>
                                                                    <span className="font-mono text-slate-600">{chain.min_residue}–{chain.max_residue}</span>
                                                                </div>
                                                                <div className="flex justify-between text-[10px] font-bold">
                                                                    <span className="text-slate-400 uppercase">{quality?.is_plddt ? 'pLDDT promedio' : 'B-factor promedio'}</span>
                                                                    <span className="font-mono" style={{ color: chain.avg_bfactor > 70 ? '#16a34a' : chain.avg_bfactor > 50 ? '#ca8a04' : '#dc2626' }}>
                                                                        {chain.avg_bfactor}
                                                                    </span>
                                                                </div>
                                                            </div>
                                                        </div>
                                                    )
                                                })}
                                            </div>
                                        </>
                                    )}
                                </div>
                            )}

                            {/* ══════ TAB: Primary Sequence ══════ */}
                            {activeTab === 'primary' && (
                                <div className="space-y-8 animate-in fade-in zoom-in-95 duration-500">
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <div className="flex items-center justify-between mb-10">
                                            <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter">Mapa de Secuencia Polipeptídica</h4>
                                        </div>
                                        <div className="bg-slate-50 rounded-3xl p-8 border border-slate-100 shadow-inner">
                                            <div className="font-mono text-sm leading-[2.2] tracking-[0.2em] break-all select-all text-slate-700">
                                                {protein.full_sequence?.split('').map((aa, i) => {
                                                    const prop = AA_PROPERTIES[aa]
                                                    // Color by SS assignment if available
                                                    const ssType = ss?.assignment?.[String(i + 1)]
                                                    const ssBorder = ssType === 'H' ? '2px solid #e11d48' : ssType === 'E' ? '2px solid #2563eb' : 'none'
                                                    return (
                                                        <span
                                                            key={i}
                                                            style={{ color: prop?.color || '#64748b', borderBottom: ssBorder }}
                                                            title={`${aa}${i + 1}: ${prop?.name || aa}${ssType ? ` (${ssType === 'H' ? 'Hélice' : 'Lámina'})` : ''}`}
                                                            className="inline-block hover:scale-150 hover:bg-blue-100 hover:shadow-xl transition-all cursor-crosshair px-0.5 rounded"
                                                        >{aa}</span>
                                                    )
                                                })}
                                            </div>
                                        </div>
                                        <div className="flex flex-wrap gap-8 mt-10 p-6 bg-slate-50/50 rounded-2xl border border-slate-100">
                                            {[
                                                { color: '#10b981', label: 'Hidrofóbico', desc: 'A, I, L, M, F, W, V' },
                                                { color: '#8b5cf6', label: 'Polar', desc: 'N, Q, S, T, Y' },
                                                { color: '#3b82f6', label: 'Básico (+)', desc: 'K, R, H' },
                                                { color: '#ef4444', label: 'Ácido (−)', desc: 'D, E' },
                                                { color: '#f59e0b', label: 'Especial', desc: 'P, G, C' },
                                            ].map(item => (
                                                <div key={item.label} className="space-y-1">
                                                    <div className="flex items-center gap-2">
                                                        <div className="w-2 h-2 rounded-full shadow-sm" style={{ backgroundColor: item.color }} />
                                                        <span className="text-[10px] font-black uppercase tracking-widest text-slate-600">{item.label}</span>
                                                    </div>
                                                    <p className="text-[9px] text-slate-400 font-mono pl-4">{item.desc}</p>
                                                </div>
                                            ))}
                                            {ss && (
                                                <>
                                                    <div className="w-px bg-slate-200" />
                                                    <div className="space-y-1">
                                                        <div className="flex items-center gap-2">
                                                            <div className="w-6 h-0.5 bg-rose-500 rounded" />
                                                            <span className="text-[10px] font-black uppercase tracking-widest text-slate-600">Hélice α</span>
                                                        </div>
                                                    </div>
                                                    <div className="space-y-1">
                                                        <div className="flex items-center gap-2">
                                                            <div className="w-6 h-0.5 bg-blue-600 rounded" />
                                                            <span className="text-[10px] font-black uppercase tracking-widest text-slate-600">Lámina β</span>
                                                        </div>
                                                    </div>
                                                </>
                                            )}
                                        </div>
                                    </div>
                                </div>
                            )}

                            {/* ══════ TAB: Composition ══════ */}
                            {activeTab === 'composition' && aaComp && (
                                <div className="space-y-10 animate-in fade-in slide-in-from-bottom-8 duration-700">
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter mb-10">Perfil Fisicoquímico</h4>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-5 gap-6">
                                            {aaComp.groups.map(g => (
                                                <div key={g.group} className="bg-slate-50/50 rounded-3xl p-6 border border-slate-100 hover:border-blue-200 transition-all hover:-translate-y-1 shadow-sm">
                                                    <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest mb-4">{g.group}</p>
                                                    <div className="flex items-baseline gap-1 mb-4">
                                                        <p className="text-3xl font-black text-slate-900">{g.pct}</p>
                                                        <p className="text-sm font-black text-slate-400">%</p>
                                                    </div>
                                                    <div className="w-full h-1 bg-slate-200 rounded-full overflow-hidden">
                                                        <div className="h-full bg-gradient-to-r from-blue-500 to-indigo-500 transition-all duration-1000" style={{ width: `${g.pct}%` }} />
                                                    </div>
                                                    <p className="text-[10px] font-bold text-slate-400 mt-4">{g.count} Residuos</p>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter mb-10">Aminoácidos Dominantes</h4>
                                        <div className="grid grid-cols-2 lg:grid-cols-5 gap-6">
                                            {aaComp.topAA.map((aa) => (
                                                <div key={aa.aa} className="bg-slate-50/50 rounded-3xl p-8 border border-slate-100 text-center relative overflow-hidden group hover:border-blue-200 transition-all shadow-sm">
                                                    <div className="text-5xl font-mono font-black mb-4 opacity-60 group-hover:opacity-100 group-hover:scale-110 transition-all duration-500" style={{ color: aa.color }}>{aa.aa}</div>
                                                    <p className="text-2xl font-black text-slate-900 mb-1">{aa.pct}%</p>
                                                    <p className="text-[10px] font-black uppercase tracking-widest text-slate-400">{aa.name}</p>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                </div>
                            )}
                        </>
                    )}
                </div>

                {/* Footer */}
                <div className="bg-white/80 backdrop-blur-md border-t border-slate-200 px-6 py-4 flex-shrink-0">
                    <div className="flex items-center justify-center gap-8 text-[10px] text-slate-400 font-bold uppercase tracking-widest">
                        <span className="flex items-center gap-3">
                            <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">ESC</kbd> Cerrar
                        </span>
                        <span className="flex items-center gap-3">
                            <div className="flex gap-1">
                                <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">←</kbd>
                                <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">→</kbd>
                            </div>
                            Navegar
                        </span>
                    </div>
                </div>
            </div>
        </div>
    )
}