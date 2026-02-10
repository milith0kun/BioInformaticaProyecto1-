/**
 * ProteinStructureViewer Component — Fully Redesigned
 * Fullscreen modal with navigation, responsive design, and modern UI
 */
import { useState, useEffect } from 'react'
import MolstarProteinViewer from './MolstarProteinViewer'

// Amino acid properties
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
    D: { name: 'Ácido aspártico', group: 'Cargado (−)', color: '#ef4444' },
    E: { name: 'Ácido glutámico', group: 'Cargado (−)', color: '#ef4444' },
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

// Secondary structure prediction
const predictSecondaryStructure = (sequence) => {
    const helixPropensity = { A: 1.45, E: 1.53, L: 1.34, M: 1.20, K: 1.07, F: 1.12, Q: 1.17, W: 1.14, I: 1.00, V: 0.90, D: 0.98, H: 1.24, R: 0.79, T: 0.82, S: 0.79, C: 0.77, Y: 0.61, N: 0.73, G: 0.53, P: 0.59 }
    const sheetPropensity = { V: 1.87, I: 1.60, F: 1.43, Y: 1.39, W: 1.19, L: 1.22, T: 1.20, C: 1.30, M: 1.07, Q: 0.98, R: 0.90, H: 0.80, A: 0.97, N: 0.65, G: 0.81, K: 0.74, S: 0.72, E: 0.26, D: 0.45, P: 0.62 }

    const structure = []
    const windowSize = 7

    for (let i = 0; i < sequence.length; i++) {
        let helixScore = 0
        let sheetScore = 0
        let count = 0

        for (let j = Math.max(0, i - Math.floor(windowSize / 2)); j <= Math.min(sequence.length - 1, i + Math.floor(windowSize / 2)); j++) {
            const aa = sequence[j]
            if (helixPropensity[aa]) {
                helixScore += helixPropensity[aa]
                sheetScore += sheetPropensity[aa]
                count++
            }
        }

        if (count > 0) {
            helixScore /= count
            sheetScore /= count
        }

        if (helixScore > 1.1 && helixScore > sheetScore) {
            structure.push({ type: 'helix', confidence: Math.min(helixScore, 2) - 1 })
        } else if (sheetScore > 1.1 && sheetScore > helixScore) {
            structure.push({ type: 'sheet', confidence: Math.min(sheetScore, 2) - 1 })
        } else {
            structure.push({ type: 'loop', confidence: 0 })
        }
    }

    return structure
}

export default function ProteinStructureViewer({
    protein,
    onClose,
    onNavigate,
    currentIndex,
    totalProteins,
    loading
}) {
    const [activeTab, setActiveTab] = useState('tertiary')
    const [secondaryStructure, setSecondaryStructure] = useState(null)

    useEffect(() => {
        if (protein?.full_sequence) {
            const structure = predictSecondaryStructure(protein.full_sequence)
            setSecondaryStructure(structure)
        }
    }, [protein])

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

    // Calculate amino acid composition
    const composition = () => {
        if (!protein?.full_sequence) return null
        const seq = protein.full_sequence
        const counts = {}
        for (const aa of seq) {
            counts[aa] = (counts[aa] || 0) + 1
        }

        const groups = { 'Hidrofóbico': 0, 'Polar': 0, 'Cargado (+)': 0, 'Cargado (−)': 0, 'Especial': 0 }
        Object.entries(counts).forEach(([aa, count]) => {
            const prop = AA_PROPERTIES[aa]
            if (prop && groups[prop.group] !== undefined) {
                groups[prop.group] += count
            }
        })

        const totalAA = seq.length
        return {
            counts,
            total: totalAA,
            groups: Object.entries(groups).map(([group, count]) => ({
                group,
                count,
                pct: ((count / totalAA) * 100).toFixed(1)
            })),
            topAA: Object.entries(counts)
                .sort((a, b) => b[1] - a[1])
                .slice(0, 5)
                .map(([aa, count]) => ({
                    aa,
                    count,
                    pct: ((count / totalAA) * 100).toFixed(1),
                    ...(AA_PROPERTIES[aa] || { name: aa, group: 'Otro', color: '#94a3b8' })
                }))
        }
    }

    const aaComp = composition()

    // Calculate secondary structure regions
    const secondaryRegions = () => {
        if (!secondaryStructure || !protein?.full_sequence) return null

        const regions = []
        let current = { type: secondaryStructure[0].type, start: 0, end: 0 }

        for (let i = 1; i < secondaryStructure.length; i++) {
            if (secondaryStructure[i].type === current.type) {
                current.end = i
            } else {
                regions.push({ ...current })
                current = { type: secondaryStructure[i].type, start: i, end: i }
            }
        }
        regions.push(current)

        const totalLength = protein.full_sequence.length
        const stats = {
            helix: ((regions.filter(r => r.type === 'helix').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / totalLength) * 100).toFixed(1),
            sheet: ((regions.filter(r => r.type === 'sheet').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / totalLength) * 100).toFixed(1),
            loop: ((regions.filter(r => r.type === 'loop').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / totalLength) * 100).toFixed(1)
        }

        return { regions, stats }
    }

    const secRegions = secondaryRegions()

    if (!protein) return null

    return (
        <div 
            className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 backdrop-blur-sm p-4 overflow-hidden"
            onClick={(e) => e.target === e.currentTarget && onClose()}
        >
            <div className="relative w-full h-full max-w-7xl max-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 rounded-3xl shadow-2xl overflow-hidden border border-slate-700/50 flex flex-col" onClick={e => e.stopPropagation()}>

                {/* Header with Navigation */}
                <div className="relative bg-gradient-to-r from-slate-800/90 to-slate-900/90 backdrop-blur-md border-b border-slate-700/50 px-6 py-4 flex-shrink-0">
                    <div className="flex items-center justify-between gap-4">
                        {/* Close Button */}
                        <button
                            onClick={onClose}
                            className="flex items-center gap-2 px-4 py-2 bg-slate-700/50 hover:bg-slate-700 text-slate-300 hover:text-white rounded-xl transition-all group"
                        >
                            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                            </svg>
                            <span className="text-sm font-medium">Cerrar</span>
                        </button>

                        {/* Protein Info */}
                        <div className="flex-1 min-w-0 text-center">
                            <h3 className="font-mono text-lg font-bold text-cyan-400 truncate">
                                {protein.protein_id || protein.locus_tag}
                            </h3>
                            <p className="text-sm text-slate-400 truncate">{protein.product || 'Proteína hipotética'}</p>
                        </div>

                        {/* Navigation Controls */}
                        {onNavigate && (
                            <div className="flex items-center gap-2">
                                <button
                                    onClick={() => onNavigate('prev')}
                                    disabled={currentIndex === 0}
                                    className="p-2 bg-slate-700/50 hover:bg-slate-700 text-slate-300 hover:text-white rounded-xl transition-all disabled:opacity-30 disabled:cursor-not-allowed"
                                    title="Proteína anterior (←)"
                                >
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
                                    </svg>
                                </button>

                                <div className="px-3 py-1.5 bg-slate-700/30 rounded-lg text-sm font-mono text-slate-300">
                                    {currentIndex + 1} / {totalProteins}
                                </div>

                                <button
                                    onClick={() => onNavigate('next')}
                                    disabled={currentIndex === totalProteins - 1}
                                    className="p-2 bg-slate-700/50 hover:bg-slate-700 text-slate-300 hover:text-white rounded-xl transition-all disabled:opacity-30 disabled:cursor-not-allowed"
                                    title="Siguiente proteína (→)"
                                >
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                                    </svg>
                                </button>
                            </div>
                        )}
                    </div>

                    {/* Breadcrumb / Stats Bar */}
                    <div className="flex items-center gap-4 mt-3 text-xs text-slate-400 flex-wrap">
                        <div className="flex items-center gap-1.5">
                            <svg className="w-4 h-4 text-cyan-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                            </svg>
                            <span className="font-mono text-cyan-400">{protein.length}</span> aminoácidos
                        </div>
                        <div className="flex items-center gap-1.5">
                            <svg className="w-4 h-4 text-purple-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 6l3 1m0 0l-3 9a5.002 5.002 0 006.001 0M6 7l3 9M6 7l6-2m6 2l3-1m-3 1l-3 9a5.002 5.002 0 006.001 0M18 7l3 9m-3-9l-6-2m0-2v2m0 16V5m0 16H9m3 0h3" />
                            </svg>
                            <span className="font-mono text-purple-400">{(protein.molecular_weight_approx / 1000).toFixed(1)}</span> kDa
                        </div>
                        {protein.gene_name && (
                            <div className="px-2 py-0.5 bg-emerald-500/20 text-emerald-400 rounded-md font-medium">
                                Gen: {protein.gene_name}
                            </div>
                        )}
                    </div>
                </div>

                {/* Tabs */}
                <div className="bg-slate-900/80 backdrop-blur-xl border-b border-white/5 px-8 flex-shrink-0 sticky top-0 z-10">
                    <div className="flex gap-8 overflow-x-auto scrollbar-hide">
                        {[
                            { id: 'tertiary', label: 'ESTRUCTURA 3D Y 4D', icon: '💎' },
                            { id: 'primary', label: 'SECUENCIA PRIMARIA', icon: '🧬' },
                            { id: 'secondary', label: 'PREDICCIÓN SECUNDARIA', icon: '🌀' },
                            { id: 'composition', label: 'ANÁLISIS DE RESIDUOS', icon: '📊' },
                        ].map(tab => (
                            <button
                                key={tab.id}
                                onClick={() => setActiveTab(tab.id)}
                                className={`group relative py-6 text-[10px] font-black uppercase tracking-[0.2em] transition-all whitespace-nowrap ${
                                    activeTab === tab.id
                                        ? 'text-cyan-400'
                                        : 'text-slate-500 hover:text-slate-300'
                                }`}
                            >
                                <div className="flex items-center gap-3">
                                    <span className={`text-base transition-transform duration-500 group-hover:scale-125 ${activeTab === tab.id ? 'scale-110' : ''}`}>
                                        {tab.icon}
                                    </span>
                                    <span>{tab.label}</span>
                                </div>
                                {activeTab === tab.id && (
                                    <div className="absolute bottom-0 left-0 right-0 h-1 bg-gradient-to-r from-cyan-500 to-indigo-500 rounded-full shadow-[0_0_12px_rgba(34,211,238,0.5)]"></div>
                                )}
                            </button>
                        ))}
                    </div>
                </div>

                {/* Content Area - Scrollable */}
                <div className="flex-1 overflow-y-auto p-8 space-y-10 custom-scrollbar">
                    {loading && (
                        <div className="flex flex-col items-center justify-center py-40">
                            <div className="relative w-20 h-20 mb-6">
                                <div className="absolute inset-0 border-4 border-slate-800 rounded-full"></div>
                                <div className="absolute inset-0 border-4 border-transparent border-t-cyan-500 rounded-full animate-spin"></div>
                            </div>
                            <p className="text-slate-500 text-[10px] font-black uppercase tracking-[0.3em] animate-pulse">Calculando Geometría Molecular...</p>
                        </div>
                    )}

                    {!loading && activeTab === 'tertiary' && (
                        <div className="space-y-8 animate-in fade-in slide-in-from-bottom-4 duration-700">
                            {/* 3D Viewer Main Card - Full Width */}
                            <div className="w-full bg-slate-950 rounded-[2.5rem] border border-white/5 overflow-hidden shadow-2xl relative group min-h-[600px]">
                                <MolstarProteinViewer
                                    proteinId={protein.protein_id}
                                    proteinSequence={protein.full_sequence}
                                    productName={protein.product}
                                />
                                
                                {/* Minimal HUD overlay for interaction hints */}
                                <div className="absolute bottom-8 left-1/2 -translate-x-1/2 flex gap-10 px-8 py-3 bg-slate-900/40 backdrop-blur-xl rounded-full border border-white/5 opacity-0 group-hover:opacity-100 transition-all duration-500 pointer-events-none">
                                    <div className="flex items-center gap-2">
                                        <span className="text-[10px] font-black text-cyan-400 uppercase tracking-widest">LMB</span>
                                        <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Rotar</span>
                                    </div>
                                    <div className="flex items-center gap-2">
                                        <span className="text-[10px] font-black text-cyan-400 uppercase tracking-widest">RMB</span>
                                        <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Mover</span>
                                    </div>
                                    <div className="flex items-center gap-2">
                                        <span className="text-[10px] font-black text-cyan-400 uppercase tracking-widest">Rueda</span>
                                        <span className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Zoom</span>
                                    </div>
                                </div>
                            </div>

                            {/* Structural Context - Horizontal Layout Below Viewer */}
                            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
                                <div className="bg-slate-900/50 rounded-3xl border border-white/5 p-6">
                                    <div className="flex items-center gap-3 mb-3">
                                        <div className="w-8 h-8 bg-indigo-500/20 rounded-xl flex items-center justify-center border border-indigo-500/30 text-indigo-400 font-bold text-xs">3°</div>
                                        <h5 className="text-[10px] font-black text-slate-300 uppercase tracking-widest">Plegamiento Terciario</h5>
                                    </div>
                                    <p className="text-[10px] text-slate-500 leading-relaxed font-medium">
                                        Conformación 3D completa mapeando núcleos hidrofóbicos y puentes iónicos.
                                    </p>
                                </div>

                                <div className="bg-slate-900/50 rounded-3xl border border-white/5 p-6">
                                    <div className="flex items-center gap-3 mb-3">
                                        <div className="w-8 h-8 bg-cyan-500/20 rounded-xl flex items-center justify-center border border-cyan-500/30 text-cyan-400 font-bold text-xs">4°</div>
                                        <h5 className="text-[10px] font-black text-slate-300 uppercase tracking-widest">Complejo Cuaternario</h5>
                                    </div>
                                    <p className="text-[10px] text-slate-500 leading-relaxed font-medium">
                                        Detección de ensamblaje de subunidades. Las cadenas se colorean automáticamente.
                                    </p>
                                </div>

                                <div className="bg-slate-900/50 rounded-3xl border border-white/5 p-6">
                                    <h5 className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em] mb-4">Estado de Resolución</h5>
                                    <div className="space-y-3">
                                        <div className="flex justify-between text-[9px] font-bold">
                                            <span className="text-slate-400 uppercase">Precisión Geometría</span>
                                            <span className="text-emerald-400">ALTA</span>
                                        </div>
                                        <div className="w-full h-1 bg-slate-950 rounded-full overflow-hidden">
                                            <div className="h-full bg-emerald-500 w-[92%] shadow-[0_0_8px_rgba(16,185,129,0.4)]"></div>
                                        </div>
                                    </div>
                                </div>

                                <div className="bg-gradient-to-br from-indigo-600 to-blue-700 rounded-3xl p-6 text-white shadow-xl shadow-indigo-900/20">
                                    <p className="text-[10px] font-black uppercase tracking-widest mb-2 opacity-80">Consejo Pro</p>
                                    <p className="text-xs font-bold leading-tight">Usa el menú superior en Mol* y cambia a 'Spacefill' para ver superficies.</p>
                                </div>
                            </div>
                        </div>
                    )}

                    {!loading && activeTab === 'primary' && (
                        <div className="space-y-8 animate-in fade-in zoom-in-95 duration-500">
                            <div className="bg-slate-900/50 rounded-[2.5rem] border border-white/5 p-10 shadow-2xl">
                                <div className="flex items-center justify-between mb-10">
                                    <h4 className="text-lg font-black text-slate-100 uppercase tracking-tighter">Mapa de Secuencia Polipeptídica</h4>
                                    <button className="px-4 py-2 bg-slate-800 hover:bg-slate-700 text-[10px] font-black uppercase tracking-widest rounded-xl transition-all border border-white/5">
                                        Descargar FASTA
                                    </button>
                                </div>
                                
                                <div className="bg-slate-950/80 rounded-3xl p-8 border border-white/5 inner-shadow">
                                    <div className="font-mono text-sm leading-[2.2] tracking-[0.2em] break-all select-all">
                                        {protein.full_sequence?.split('').map((aa, i) => {
                                            const prop = AA_PROPERTIES[aa]
                                            return (
                                                <span
                                                    key={i}
                                                    style={{ color: prop?.color || '#94a3b8' }}
                                                    title={`${aa}${i + 1}: ${prop?.name || aa}`}
                                                    className="inline-block hover:scale-150 hover:bg-white/10 hover:shadow-xl transition-all cursor-crosshair px-0.5 rounded"
                                                >
                                                    {aa}
                                                </span>
                                            )
                                        })}
                                    </div>
                                </div>

                                {/* Legend - Designer Style */}
                                <div className="flex flex-wrap gap-8 mt-10 p-6 bg-slate-900/30 rounded-2xl border border-white/5">
                                    {[
                                        { color: '#10b981', label: 'Hidrofóbico', desc: 'A, I, L, M, F, W, V' },
                                        { color: '#8b5cf6', label: 'Polar', desc: 'N, Q, S, T, Y' },
                                        { color: '#3b82f6', label: 'Básico (+)', desc: 'K, R, H' },
                                        { color: '#ef4444', label: 'Ácido (−)', desc: 'D, E' },
                                        { color: '#f59e0b', label: 'Especial', desc: 'P, G, C' },
                                    ].map(item => (
                                        <div key={item.label} className="space-y-1">
                                            <div className="flex items-center gap-2">
                                                <div className="w-2 h-2 rounded-full shadow-lg" style={{ backgroundColor: item.color }}></div>
                                                <span className="text-[10px] font-black uppercase tracking-widest text-slate-300">{item.label}</span>
                                            </div>
                                            <p className="text-[9px] text-slate-500 font-mono pl-4">{item.desc}</p>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>
                    )}

                    {!loading && activeTab === 'secondary' && secRegions && (
                        <div className="space-y-8 animate-in fade-in slide-in-from-right-8 duration-700">
                            <div className="bg-slate-900/50 rounded-[2.5rem] border border-white/5 p-10 shadow-2xl">
                                <h4 className="text-lg font-black text-slate-100 uppercase tracking-tighter mb-10">Propensión de Estructura Secundaria</h4>
                                
                                {/* Visual representation - Pro Design */}
                                <div className="relative h-32 bg-slate-950 rounded-2xl border border-white/5 overflow-hidden shadow-inner flex">
                                    {secRegions.regions.map((region, i) => {
                                        const width = ((region.end - region.start + 1) / protein.full_sequence.length) * 100
                                        const colors = {
                                            helix: 'from-rose-500 to-pink-600 shadow-[0_0_15px_rgba(244,63,94,0.3)]',
                                            sheet: 'from-cyan-500 to-blue-600 shadow-[0_0_15px_rgba(6,182,212,0.3)]',
                                            loop: 'from-slate-700 to-slate-800'
                                        }
                                        return (
                                            <div
                                                key={i}
                                                className={`h-full bg-gradient-to-b ${colors[region.type]} transition-all hover:brightness-125 border-x border-black/20 flex items-center justify-center group/region relative`}
                                                style={{ width: `${width}%` }}
                                            >
                                                {width > 2 && (
                                                    <span className="text-[8px] font-black text-white/40 uppercase rotate-90 scale-75 lg:rotate-0 lg:scale-100">
                                                        {region.type === 'helix' ? 'H' : region.type === 'sheet' ? 'L' : 'O'}
                                                    </span>
                                                )}
                                                {/* Tooltip on hover */}
                                                <div className="absolute -top-12 left-1/2 -translate-x-1/2 px-3 py-1 bg-slate-800 rounded-lg text-[10px] font-bold text-white opacity-0 group-hover/region:opacity-100 pointer-events-none whitespace-nowrap z-20 border border-white/10 shadow-xl transition-opacity">
                                                    {region.type === 'helix' ? 'Hélice Alfa' : region.type === 'sheet' ? 'Lámina Beta' : 'Ovillo Aleatorio'}
                                                    <span className="ml-2 text-slate-500">[{region.start + 1}-{region.end + 1}]</span>
                                                </div>
                                            </div>
                                        )
                                    })}
                                </div>

                                {/* Modern Stats Grid */}
                                <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mt-10">
                                    {[
                                        { label: 'Hélices Alfa', val: secRegions.stats.helix, color: 'from-rose-500/20 to-pink-500/20', border: 'border-rose-500/30', text: 'text-rose-400' },
                                        { label: 'Láminas Beta', val: secRegions.stats.sheet, color: 'from-cyan-500/20 to-blue-500/20', border: 'border-cyan-500/30', text: 'text-cyan-400' },
                                        { label: 'Ovillos Aleatorios', val: secRegions.stats.loop, color: 'from-slate-600/20 to-slate-500/20', border: 'border-slate-500/30', text: 'text-slate-400' }
                                    ].map(stat => (
                                        <div key={stat.label} className={`bg-gradient-to-br ${stat.color} rounded-3xl p-6 border ${stat.border} flex flex-col items-center text-center`}>
                                            <p className="text-[10px] font-black uppercase tracking-[0.2em] text-slate-400 mb-2">{stat.label}</p>
                                            <p className={`text-4xl font-black ${stat.text}`}>{stat.val}%</p>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>
                    )}

                    {!loading && activeTab === 'composition' && aaComp && (
                        <div className="space-y-10 animate-in fade-in slide-in-from-bottom-8 duration-700">
                            {/* Detailed Residue Distribution */}
                            <div className="bg-slate-900/50 rounded-[2.5rem] border border-white/5 p-10 shadow-2xl">
                                <h4 className="text-lg font-black text-slate-100 uppercase tracking-tighter mb-10">Perfil Fisicoquímico</h4>
                                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-5 gap-6">
                                    {aaComp.groups.map(g => (
                                        <div key={g.group} className="bg-slate-950/50 rounded-3xl p-6 border border-white/5 hover:border-cyan-500/30 transition-all hover:-translate-y-1">
                                            <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest mb-4">{g.group}</p>
                                            <div className="flex items-baseline gap-1 mb-4">
                                                <p className="text-3xl font-black text-white">{g.pct}</p>
                                                <p className="text-sm font-black text-slate-600">%</p>
                                            </div>
                                            <div className="w-full h-1 bg-slate-900 rounded-full overflow-hidden">
                                                <div
                                                    className="h-full bg-gradient-to-r from-cyan-500 to-indigo-500 transition-all duration-1000 shadow-[0_0_8px_rgba(34,211,238,0.3)]"
                                                    style={{ width: `${g.pct}%` }}
                                                />
                                            </div>
                                            <p className="text-[10px] font-bold text-slate-500 mt-4">{g.count} Residuos</p>
                                        </div>
                                    ))}
                                </div>
                            </div>

                            {/* Most Frequent Residues */}
                            <div className="bg-slate-900/50 rounded-[2.5rem] border border-white/5 p-10 shadow-2xl">
                                <h4 className="text-lg font-black text-slate-100 uppercase tracking-tighter mb-10">Aminoácidos Dominantes</h4>
                                <div className="grid grid-cols-2 lg:grid-cols-5 gap-6">
                                    {aaComp.topAA.map((aa, i) => (
                                        <div key={aa.aa} className="bg-slate-950/50 rounded-3xl p-8 border border-white/5 text-center relative overflow-hidden group">
                                            <div className={`absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-transparent via-cyan-500 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-500`}></div>
                                            <div className="text-5xl font-mono font-black mb-4 opacity-80 group-hover:scale-110 transition-transform duration-500" style={{ color: aa.color }}>
                                                {aa.aa}
                                            </div>
                                            <p className="text-2xl font-black text-white mb-1">{aa.pct}%</p>
                                            <p className="text-[10px] font-black uppercase tracking-widest text-slate-500">{aa.name}</p>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>
                    )}
                </div>

                {/* Footer - Keyboard shortcuts hint */}
                <div className="bg-slate-800/50 backdrop-blur-sm border-t border-slate-700/50 px-6 py-3 flex-shrink-0">
                    <div className="flex items-center justify-center gap-6 text-xs text-slate-500">
                        <span className="flex items-center gap-2">
                            <kbd className="px-2 py-1 bg-slate-700/50 rounded border border-slate-600 text-slate-300">ESC</kbd>
                            Cerrar
                        </span>
                        <span className="flex items-center gap-2">
                            <kbd className="px-2 py-1 bg-slate-700/50 rounded border border-slate-600 text-slate-300">←</kbd>
                            <kbd className="px-2 py-1 bg-slate-700/50 rounded border border-slate-600 text-slate-300">→</kbd>
                            Navegar
                        </span>
                    </div>
                </div>
            </div>
        </div>
    )
}
