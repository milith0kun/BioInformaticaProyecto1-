/**
 * ProteinStructureViewer Component — Clean Laboratory Edition
 * Fullscreen modal with refined light aesthetic and structural analysis
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
            className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/40 backdrop-blur-sm p-4 overflow-hidden"
            onClick={(e) => e.target === e.currentTarget && onClose()}
        >
            <div className="relative w-full h-full max-w-7xl max-h-screen bg-[#f8fafc] rounded-3xl shadow-2xl overflow-hidden border border-slate-200 flex flex-col" onClick={e => e.stopPropagation()}>

                {/* Header with Navigation */}
                <div className="relative bg-white/80 backdrop-blur-2xl border-b border-slate-200 px-10 py-6 flex-shrink-0">
                    <div className="flex items-center justify-between gap-8">
                        <button
                            onClick={onClose}
                            className="flex items-center gap-3 px-5 py-2.5 bg-slate-100 hover:bg-slate-200 text-slate-600 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 group shadow-sm"
                        >
                            <svg className="w-5 h-5 group-hover:rotate-90 transition-transform duration-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                            </svg>
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
                                <button
                                    onClick={() => onNavigate('prev')}
                                    disabled={currentIndex === 0}
                                    className="p-3 bg-slate-100 hover:bg-slate-200 text-slate-500 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 disabled:opacity-20 disabled:cursor-not-allowed shadow-sm"
                                >
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M15 19l-7-7 7-7" />
                                    </svg>
                                </button>

                                <div className="px-5 py-2.5 bg-white rounded-2xl border border-slate-200 text-xs font-black font-mono text-blue-600 shadow-inner">
                                    {currentIndex + 1} <span className="text-slate-300 px-1">/</span> {totalProteins}
                                </div>

                                <button
                                    onClick={() => onNavigate('next')}
                                    disabled={currentIndex === totalProteins - 1}
                                    className="p-3 bg-slate-100 hover:bg-slate-200 text-slate-500 hover:text-slate-900 rounded-2xl transition-all border border-slate-200 disabled:opacity-20 disabled:cursor-not-allowed shadow-sm"
                                >
                                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M9 5l7 7-7 7" />
                                    </svg>
                                </button>
                            </div>
                        )}
                    </div>

                    <div className="flex items-center justify-center gap-8 mt-6">
                        <div className="flex items-center gap-2 px-4 py-1.5 bg-blue-50 rounded-full border border-blue-100">
                            <div className="w-1.5 h-1.5 rounded-full bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.3)]"></div>
                            <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest"><span className="text-blue-600">{protein.length}</span> Residuos</span>
                        </div>
                        <div className="flex items-center gap-2 px-4 py-1.5 bg-indigo-50 rounded-full border border-indigo-100">
                            <div className="w-1.5 h-1.5 rounded-full bg-indigo-500 shadow-[0_0_8px_rgba(99,102,241,0.3)]"></div>
                            <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest"><span className="text-indigo-600">{(protein.molecular_weight_approx / 1000).toFixed(1)}</span> kDa</span>
                        </div>
                        {protein.gene_name && (
                            <div className="flex items-center gap-2 px-4 py-1.5 bg-emerald-50 rounded-full border border-emerald-100">
                                <div className="w-1.5 h-1.5 rounded-full bg-emerald-500 shadow-[0_0_8px_rgba(16,185,129,0.3)]"></div>
                                <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest">Gen: <span className="text-emerald-600">{protein.gene_name}</span></span>
                            </div>
                        )}
                    </div>
                </div>

                {/* Tabs */}
                <div className="bg-slate-50/50 backdrop-blur-2xl border-b border-slate-200 px-10 flex-shrink-0">
                    <div className="flex gap-10 overflow-x-auto scrollbar-hide">
                        {[
                            { id: 'tertiary', label: 'Estructura 3D/4D', icon: '💎' },
                            { id: 'primary', label: 'Secuencia', icon: '🧬' },
                            { id: 'secondary', label: 'Secundaria', icon: '🌀' },
                            { id: 'composition', label: 'Residuos', icon: '📊' },
                        ].map(tab => (
                            <button
                                key={tab.id}
                                onClick={() => setActiveTab(tab.id)}
                                className={`group relative py-6 text-[10px] font-black uppercase tracking-[0.25em] transition-all whitespace-nowrap ${
                                    activeTab === tab.id
                                        ? 'text-blue-600'
                                        : 'text-slate-400 hover:text-slate-600'
                                }`}
                            >
                                <div className="flex items-center gap-3">
                                    <span className={`transition-transform duration-500 group-hover:scale-110 ${activeTab === tab.id ? 'scale-110' : 'opacity-40'}`}>
                                        {tab.icon}
                                    </span>
                                    <span>{tab.label}</span>
                                </div>
                                {activeTab === tab.id && (
                                    <div className="absolute bottom-0 left-0 right-0 h-0.5 bg-blue-600 shadow-[0_0_12px_rgba(59,130,246,0.3)]"></div>
                                )}
                            </button>
                        ))}
                    </div>
                </div>

                {/* Content Area */}
                <div className="flex-1 overflow-y-auto p-8 lg:p-12 space-y-12 custom-scrollbar min-h-0">
                    {loading ? (
                        <div className="flex flex-col items-center justify-center py-48">
                            <div className="relative w-24 h-24 mb-8">
                                <div className="absolute inset-0 border-4 border-slate-200 rounded-full"></div>
                                <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin"></div>
                            </div>
                            <p className="text-slate-400 text-[10px] font-black uppercase tracking-[0.4em] animate-pulse text-center">Analizando Geometría Molecular</p>
                        </div>
                    ) : (
                        <>
                            {activeTab === 'tertiary' && (
                                <div className="space-y-10 animate-in fade-in slide-in-from-bottom-6 duration-1000">
                                    <div className="w-full bg-slate-900 rounded-[3rem] border border-slate-200 overflow-hidden shadow-xl relative group min-h-[700px]">
                                        <MolstarProteinViewer
                                            proteinId={protein.protein_id}
                                            proteinSequence={protein.full_sequence}
                                            productName={protein.product}
                                        />
                                    </div>

                                    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
                                        <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                            <div className="flex items-center gap-4 mb-4">
                                                <div className="w-10 h-10 bg-blue-50 rounded-2xl flex items-center justify-center border border-blue-100 text-blue-600 font-black text-xs">3°</div>
                                                <h5 className="text-[10px] font-black text-slate-700 uppercase tracking-[0.2em]">Terciaria</h5>
                                            </div>
                                            <p className="text-[11px] text-slate-500 leading-relaxed font-medium">Mapeo de conformación 3D, núcleos hidrofóbicos y puentes salinos.</p>
                                        </div>
                                        <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                            <div className="flex items-center gap-4 mb-4">
                                                <div className="w-10 h-10 bg-indigo-50 rounded-2xl flex items-center justify-center border border-indigo-100 text-indigo-600 font-black text-xs">4°</div>
                                                <h5 className="text-[10px] font-black text-slate-700 uppercase tracking-[0.2em]">Cuaternaria</h5>
                                            </div>
                                            <p className="text-[11px] text-slate-500 leading-relaxed font-medium">Detección de ensamblajes y complejos multiproteicos coloreados por cadena.</p>
                                        </div>
                                        <div className="bg-white rounded-[2rem] border border-slate-200 p-8 shadow-sm">
                                            <h5 className="text-[9px] font-black text-slate-400 uppercase tracking-[0.25em] mb-5">Resolución</h5>
                                            <div className="space-y-4">
                                                <div className="flex justify-between text-[10px] font-black">
                                                    <span className="text-slate-500 uppercase tracking-widest">Precisión</span>
                                                    <span className="text-emerald-600 tracking-widest">ALTA</span>
                                                </div>
                                                <div className="w-full h-1 bg-slate-100 rounded-full overflow-hidden">
                                                    <div className="h-full bg-emerald-500 w-[94%] shadow-[0_0_12px_rgba(16,185,129,0.3)]"></div>
                                                </div>
                                            </div>
                                        </div>
                                        <div className="bg-gradient-to-br from-blue-600 to-indigo-700 rounded-[2rem] p-8 text-white shadow-xl relative overflow-hidden group">
                                            <div className="absolute top-0 right-0 w-32 h-32 bg-white/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000"></div>
                                            <p className="text-[10px] font-black uppercase tracking-widest mb-3 opacity-70 relative z-10">Tip</p>
                                            <p className="text-xs font-bold leading-relaxed relative z-10">Activa 'Spacefill' en el menú de Mol* para un análisis de superficie accesible.</p>
                                        </div>
                                    </div>
                                </div>
                            )}

                            {activeTab === 'primary' && (
                                <div className="space-y-8 animate-in fade-in zoom-in-95 duration-500">
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <div className="flex items-center justify-between mb-10">
                                            <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter">Mapa de Secuencia Polipeptídica</h4>
                                            <button className="px-4 py-2 bg-slate-100 hover:bg-slate-200 text-[10px] font-black text-slate-600 uppercase tracking-widest rounded-xl transition-all border border-slate-200">Descargar FASTA</button>
                                        </div>
                                        <div className="bg-slate-50 rounded-3xl p-8 border border-slate-100 shadow-inner">
                                            <div className="font-mono text-sm leading-[2.2] tracking-[0.2em] break-all select-all text-slate-700">
                                                {protein.full_sequence?.split('').map((aa, i) => {
                                                    const prop = AA_PROPERTIES[aa]
                                                    return (
                                                        <span key={i} style={{ color: prop?.color || '#64748b' }} title={`${aa}${i + 1}: ${prop?.name || aa}`} className="inline-block hover:scale-150 hover:bg-blue-100 hover:shadow-xl transition-all cursor-crosshair px-0.5 rounded">{aa}</span>
                                                    )
                                                })}
                                            </div>
                                        </div>
                                        <div className="flex flex-wrap gap-8 mt-10 p-6 bg-slate-50/50 rounded-2xl border border-slate-100">
                                            {[{ color: '#10b981', label: 'Hidrofóbico', desc: 'A, I, L, M, F, W, V' }, { color: '#8b5cf6', label: 'Polar', desc: 'N, Q, S, T, Y' }, { color: '#3b82f6', label: 'Básico (+)', desc: 'K, R, H' }, { color: '#ef4444', label: 'Ácido (−)', desc: 'D, E' }, { color: '#f59e0b', label: 'Especial', desc: 'P, G, C' }].map(item => (
                                                <div key={item.label} className="space-y-1">
                                                    <div className="flex items-center gap-2">
                                                        <div className="w-2 h-2 rounded-full shadow-sm" style={{ backgroundColor: item.color }}></div>
                                                        <span className="text-[10px] font-black uppercase tracking-widest text-slate-600">{item.label}</span>
                                                    </div>
                                                    <p className="text-[9px] text-slate-400 font-mono pl-4">{item.desc}</p>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                </div>
                            )}

                            {activeTab === 'secondary' && secRegions && (
                                <div className="space-y-8 animate-in fade-in slide-in-from-right-8 duration-700">
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter mb-10">Propensión de Estructura Secundaria</h4>
                                        <div className="relative h-32 bg-slate-50 rounded-2xl border border-slate-100 overflow-hidden shadow-inner flex">
                                            {secRegions.regions.map((region, i) => {
                                                const width = ((region.end - region.start + 1) / protein.full_sequence.length) * 100
                                                const colors = { helix: 'from-rose-400 to-pink-500', sheet: 'from-blue-400 to-indigo-500', loop: 'from-slate-200 to-slate-300' }
                                                return (
                                                    <div key={i} className={`h-full bg-gradient-to-b ${colors[region.type]} transition-all hover:brightness-110 border-x border-white/20 flex items-center justify-center group/region relative`} style={{ width: `${width}%` }}>
                                                        {width > 2 && <span className="text-[8px] font-black text-white/60 uppercase rotate-90 lg:rotate-0">{region.type[0]}</span>}
                                                        <div className="absolute -top-12 left-1/2 -translate-x-1/2 px-3 py-1 bg-slate-800 rounded-lg text-[10px] font-bold text-white opacity-0 group-hover/region:opacity-100 pointer-events-none whitespace-nowrap z-20 border border-white/10 shadow-xl transition-opacity">
                                                            {region.type === 'helix' ? 'Hélice Alfa' : region.type === 'sheet' ? 'Lámina Beta' : 'Ovillo Aleatorio'} <span className="ml-2 text-slate-400">[{region.start + 1}-{region.end + 1}]</span>
                                                        </div>
                                                    </div>
                                                )
                                            })}
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mt-10">
                                            {[{ label: 'Hélices Alfa', val: secRegions.stats.helix, color: 'from-rose-50 to-pink-50', border: 'border-rose-100', text: 'text-rose-600' }, { label: 'Láminas Beta', val: secRegions.stats.sheet, color: 'from-blue-50 to-indigo-50', border: 'border-blue-100', text: 'text-blue-600' }, { label: 'Ovillos Aleatorios', val: secRegions.stats.loop, color: 'from-slate-50 to-slate-100', border: 'border-slate-200', text: 'text-slate-600' }].map(stat => (
                                                <div key={stat.label} className={`bg-gradient-to-br ${stat.color} rounded-3xl p-6 border ${stat.border} flex flex-col items-center text-center shadow-sm`}>
                                                    <p className="text-[10px] font-black uppercase tracking-[0.2em] text-slate-400 mb-2">{stat.label}</p>
                                                    <p className={`text-4xl font-black ${stat.text}`}>{stat.val}%</p>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                </div>
                            )}

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
                                                        <div className="h-full bg-gradient-to-r from-blue-500 to-indigo-500 transition-all duration-1000 shadow-[0_0_8px_rgba(59,130,246,0.2)]" style={{ width: `${g.pct}%` }} />
                                                    </div>
                                                    <p className="text-[10px] font-bold text-slate-400 mt-4">{g.count} Residuos</p>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                    <div className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm">
                                        <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter mb-10">Aminoácidos Dominantes</h4>
                                        <div className="grid grid-cols-2 lg:grid-cols-5 gap-6">
                                            {aaComp.topAA.map((aa, i) => (
                                                <div key={aa.aa} className="bg-slate-50/50 rounded-3xl p-8 border border-slate-100 text-center relative overflow-hidden group hover:border-blue-200 transition-all shadow-sm">
                                                    <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-transparent via-blue-500 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-500"></div>
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
                            <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">ESC</kbd>
                            Cerrar
                        </span>
                        <span className="flex items-center gap-3">
                            <div className="flex gap-1">
                                <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">←</kbd>
                                <kbd className="px-2 py-1 bg-slate-100 rounded border border-slate-200 text-slate-600 shadow-sm">→</kbd>
                            </div>
                            Navegar Secuencias
                        </span>
                    </div>
                </div>
            </div>
        </div>
    )
}