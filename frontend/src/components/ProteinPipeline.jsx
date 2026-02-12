/**
 * ProteinPipeline — Complete Protein Structure Hierarchy Visualization
 * Shows the biological process: Primary → Secondary → Tertiary → Quaternary
 * Like a textbook diagram but interactive and with real data.
 */
import { useState, useEffect, useRef } from 'react'
import { createPluginUI } from 'molstar/lib/mol-plugin-ui/react18'
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec'
import { MolScriptBuilder as Q } from 'molstar/lib/mol-script/language/builder'
import 'molstar/build/viewer/molstar.css'

// ─── Color Maps ───
const GROUP_COLORS = {
    hydrophobic: '#10b981',
    polar: '#8b5cf6',
    positive: '#3b82f6',
    negative: '#ef4444',
    special: '#f59e0b',
    unknown: '#94a3b8'
}

const GROUP_LABELS = {
    hydrophobic: 'Hidrofóbico',
    polar: 'Polar',
    positive: 'Cargado (+)',
    negative: 'Cargado (−)',
    special: 'Especial'
}

const SS_COLORS = { H: '#e11d48', E: '#2563eb', C: '#94a3b8' }
const SS_LABELS = { H: 'α-Hélice', E: 'β-Lámina', C: 'Ovillo' }

// ─── MolstarMini: Compact Mol* viewer for pipeline ───
function MolstarMini({ pdbUrl, preset, height = 400, highlightResidues = [] }) {
    const parentRef = useRef(null)
    const pluginRef = useRef(null)
    const initRef = useRef(false)

    useEffect(() => {
        let alive = true
        async function init() {
            if (!parentRef.current || initRef.current || !pdbUrl) return
            initRef.current = true
            if (pluginRef.current) { try { pluginRef.current.dispose() } catch (e) { } }
            if (parentRef.current) parentRef.current.innerHTML = ''

            const spec = DefaultPluginUISpec()
            spec.layout = { initial: { isExpanded: false, showControls: false, controlsDisplay: 'reactive' } }

            try {
                const plugin = await createPluginUI(parentRef.current, spec)
                if (!alive) { plugin.dispose(); initRef.current = false; return }
                pluginRef.current = plugin

                const data = await plugin.builders.data.download({ url: pdbUrl, isBinary: false }, { state: { isGhost: true } })
                const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb')

                const model = await plugin.builders.structure.createModel(trajectory)
                const structure = await plugin.builders.structure.createStructure(model)
                const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer')

                if (polymer) {
                    if (preset === 'cartoon-ss') {
                        await plugin.builders.structure.representation.addRepresentation(polymer, {
                            type: 'cartoon', color: 'secondary-structure-type', size: 'uniform'
                        })
                    } else if (preset === 'cartoon-chain') {
                        await plugin.builders.structure.representation.addRepresentation(polymer, {
                            type: 'cartoon', color: 'chain-id', size: 'uniform'
                        })
                        await plugin.builders.structure.representation.addRepresentation(polymer, {
                            type: 'gaussian-surface', color: 'chain-id', typeParams: { alpha: 0.2 }
                        })
                    } else {
                        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'all-models', {
                            representationPreset: 'auto', showUnitcell: false, showWater: false
                        })
                    }

                    // Highlight Active Sites / Functional Residues
                    if (highlightResidues && highlightResidues.length > 0) {
                        // Create a selection query for the specified residues
                        const residueQuery = Q.struct.generator.atomGroups({
                            'residue-test': Q.pred.inList(Q.struct.atom.resno(), highlightResidues),
                            'group-by': Q.struct.atom.resno()
                        })

                        // Create a component for the highlighted residues
                        const highlightedComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
                            structure,
                            residueQuery,
                            'highlighted-residues',
                            { label: 'Residuos Destacados' }
                        )

                        if (highlightedComponent) {
                            // Add a ball-and-stick representation to the highlighted component
                            await plugin.builders.structure.representation.addRepresentation(highlightedComponent, {
                                type: 'ball-and-stick',
                                color: 'uniform',
                                colorParams: { value: 0x8b5cf6 }, // Purple color
                                size: 'uniform',
                                sizeParams: { value: 0.2 }
                            })
                        }
                    }
                }
            } catch (e) { console.error('MolstarMini error:', e) }
            finally { initRef.current = false }
        }
        const t = setTimeout(init, 100)
        return () => { alive = false; clearTimeout(t); if (pluginRef.current) { const p = pluginRef.current; pluginRef.current = null; setTimeout(() => { try { p.dispose() } catch (e) { } }, 0) } }
    }, [pdbUrl, preset, highlightResidues])

    return <div ref={parentRef} style={{ width: '100%', height, minHeight: 300 }} />
}

// ─── Pipeline Arrow (SVG transition between levels) ───
function PipelineArrow({ label, sublabel, color = '#3b82f6', delay = 0 }) {
    return (
        <div className="flex flex-col items-center py-6 animate-in fade-in duration-700" style={{ animationDelay: `${delay}ms` }}>
            <svg width="60" height="80" viewBox="0 0 60 80" fill="none">
                <defs>
                    <linearGradient id={`arrow-${label.replace(/\s/g, '')}`} x1="30" y1="0" x2="30" y2="80" gradientUnits="userSpaceOnUse">
                        <stop stopColor={color} stopOpacity="0.15" />
                        <stop offset="1" stopColor={color} stopOpacity="0.8" />
                    </linearGradient>
                </defs>
                <rect x="26" y="0" width="8" height="50" rx="4" fill={`url(#arrow-${label.replace(/\s/g, '')})`}>
                    <animate attributeName="height" values="30;50;30" dur="2s" repeatCount="indefinite" />
                </rect>
                <polygon points="10,50 50,50 30,78" fill={color} opacity="0.7">
                    <animate attributeName="opacity" values="0.5;0.9;0.5" dur="2s" repeatCount="indefinite" />
                </polygon>
            </svg>
            <div className="text-center mt-2">
                <p className="text-[10px] font-black uppercase tracking-[0.2em]" style={{ color }}>{label}</p>
                {sublabel && <p className="text-[9px] text-slate-400 font-medium mt-0.5">{sublabel}</p>}
            </div>
        </div>
    )
}

// ─── Level Header ───
function LevelHeader({ number, title, subtitle, color, icon }) {
    return (
        <div className="flex items-center gap-6 mb-8">
            <div className="w-16 h-16 rounded-3xl flex items-center justify-center text-white font-black text-xl shadow-lg relative overflow-hidden" style={{ backgroundColor: color }}>
                <div className="absolute inset-0 bg-white/10" />
                <span className="relative z-10">{icon}</span>
            </div>
            <div>
                <p className="text-[9px] font-black uppercase tracking-[0.3em] mb-1" style={{ color }}>NIVEL {number}</p>
                <h3 className="text-2xl font-black text-slate-900 uppercase tracking-tighter">{title}</h3>
                <p className="text-xs text-slate-500 font-medium mt-1 max-w-lg">{subtitle}</p>
            </div>
        </div>
    )
}

// ─── Primary Structure: Animated Bead Chain ───
function PrimaryBeadChain({ residues, maxDisplay = 80 }) {
    const [hoveredIdx, setHoveredIdx] = useState(null)
    const display = residues.slice(0, maxDisplay)
    const truncated = residues.length > maxDisplay

    return (
        <div className="space-y-6">
            <div className="bg-slate-50/80 rounded-3xl p-8 border border-slate-100 shadow-inner relative overflow-hidden">
                <div className="absolute top-4 right-6 flex items-center gap-2">
                    <span className="text-[8px] font-black text-slate-300 uppercase tracking-widest">NH₂</span>
                    <div className="w-16 h-px bg-gradient-to-r from-slate-300 to-transparent" />
                </div>
                <div className="flex flex-wrap gap-1 items-center">
                    <span className="text-[10px] font-black text-emerald-600 mr-2 py-1">NH₂—</span>
                    {display.map((r, i) => {
                        const color = GROUP_COLORS[r.group] || GROUP_COLORS.unknown
                        const isHovered = hoveredIdx === i
                        return (
                            <div
                                key={i}
                                className="relative group"
                                onMouseEnter={() => setHoveredIdx(i)}
                                onMouseLeave={() => setHoveredIdx(null)}
                            >
                                {/* Peptide bond line */}
                                {i > 0 && <div className="absolute -left-1 top-1/2 w-2 h-0.5 bg-slate-300 -translate-y-1/2 z-0" />}
                                {/* Amino acid bead */}
                                <div
                                    className="relative z-10 w-8 h-8 rounded-full flex items-center justify-center text-white font-mono text-[11px] font-black shadow-md cursor-crosshair transition-all duration-200"
                                    style={{
                                        backgroundColor: color,
                                        transform: isHovered ? 'scale(1.6)' : 'scale(1)',
                                        zIndex: isHovered ? 50 : 10,
                                        boxShadow: isHovered ? `0 0 20px ${color}60` : undefined,
                                        animationDelay: `${i * 15}ms`
                                    }}
                                    title={`${r.aa}${r.pos}: ${r.name} (${GROUP_LABELS[r.group] || r.group})`}
                                >
                                    {r.aa}
                                </div>
                                {/* Tooltip */}
                                {isHovered && (
                                    <div className="absolute -top-16 left-1/2 -translate-x-1/2 bg-slate-900 text-white px-3 py-2 rounded-xl shadow-2xl z-[60] whitespace-nowrap pointer-events-none">
                                        <p className="text-[10px] font-black">{r.aa}{r.pos} — {r.name}</p>
                                        <p className="text-[8px] text-slate-300">{GROUP_LABELS[r.group]} · H: {r.hydropathy}</p>
                                    </div>
                                )}
                            </div>
                        )
                    })}
                    {truncated && <span className="text-xs font-black text-slate-400 ml-2">...+{residues.length - maxDisplay} más</span>}
                    <span className="text-[10px] font-black text-rose-600 ml-2 py-1">—COOH</span>
                </div>
            </div>
            {/* Legend */}
            <div className="flex flex-wrap gap-4">
                {Object.entries(GROUP_COLORS).filter(([k]) => k !== 'unknown').map(([group, color]) => (
                    <div key={group} className="flex items-center gap-2">
                        <div className="w-3 h-3 rounded-full shadow-sm" style={{ backgroundColor: color }} />
                        <span className="text-[9px] font-black text-slate-500 uppercase tracking-widest">{GROUP_LABELS[group]}</span>
                    </div>
                ))}
            </div>
        </div>
    )
}

// ─── Secondary Structure: Ribbon Diagram + Stats ───
function SecondaryRibbon({ assignment, helices, sheets, stats, hbonds, sequence }) {
    const canvasRef = useRef(null)
    const seqLen = sequence?.length || 0

    // Draw ribbon diagram on canvas
    useEffect(() => {
        if (!canvasRef.current || !assignment || seqLen === 0) return
        const canvas = canvasRef.current
        const ctx = canvas.getContext('2d')
        const W = canvas.width = canvas.parentElement.clientWidth
        const H = canvas.height = 120

        ctx.clearRect(0, 0, W, H)

        const padding = 40
        const drawW = W - padding * 2
        const scale = drawW / seqLen
        const midY = H / 2

        // Draw backbone line
        ctx.strokeStyle = '#e2e8f0'
        ctx.lineWidth = 2
        ctx.beginPath()
        ctx.moveTo(padding, midY)
        ctx.lineTo(W - padding, midY)
        ctx.stroke()

        // Draw SS elements
        for (let i = 1; i <= seqLen; i++) {
            const ss = assignment[String(i)]
            const x = padding + (i - 1) * scale
            const w = Math.max(scale, 1.5)

            if (ss === 'H') {
                // Helix: rounded rectangle above center
                ctx.fillStyle = '#e11d4860'
                ctx.strokeStyle = '#e11d48'
                ctx.lineWidth = 1.5
                const hh = 30
                ctx.beginPath()
                ctx.roundRect(x, midY - hh, w + 0.5, hh, 2)
                ctx.fill()
                ctx.stroke()
            } else if (ss === 'E') {
                // Sheet: angular shape below center
                ctx.fillStyle = '#2563eb50'
                ctx.strokeStyle = '#2563eb'
                ctx.lineWidth = 1.5
                const sh = 24
                ctx.beginPath()
                ctx.roundRect(x, midY, w + 0.5, sh, 2)
                ctx.fill()
                ctx.stroke()
            } else {
                // Coil: thin line at center
                ctx.fillStyle = '#94a3b8'
                ctx.fillRect(x, midY - 1, w + 0.5, 2)
            }
        }

        // Draw H-bond indicators (subset)
        if (hbonds && hbonds.length > 0) {
            ctx.strokeStyle = '#f59e0b40'
            ctx.lineWidth = 0.5
            ctx.setLineDash([2, 3])
            const shown = hbonds.slice(0, 60)
            for (const hb of shown) {
                const x1 = padding + (hb.donor - 1) * scale + scale / 2
                const x2 = padding + (hb.acceptor - 1) * scale + scale / 2
                const arcH = hb.ss_type === 'H' ? -20 : 16
                ctx.beginPath()
                ctx.moveTo(x1, midY)
                ctx.quadraticCurveTo((x1 + x2) / 2, midY + arcH, x2, midY)
                ctx.stroke()
            }
            ctx.setLineDash([])
        }

        // Labels
        ctx.fillStyle = '#64748b'
        ctx.font = 'bold 9px system-ui'
        ctx.textAlign = 'left'
        ctx.fillText('N', padding - 15, midY + 4)
        ctx.textAlign = 'right'
        ctx.fillText('C', W - padding + 15, midY + 4)

        // Residue numbers
        ctx.fillStyle = '#94a3b8'
        ctx.font = '8px monospace'
        ctx.textAlign = 'center'
        const step = Math.max(Math.floor(seqLen / 10), 1)
        for (let i = 1; i <= seqLen; i += step) {
            const x = padding + (i - 1) * scale
            ctx.fillText(String(i), x, H - 5)
            ctx.fillStyle = '#e2e8f020'
            ctx.fillRect(x, 10, 0.5, H - 20)
            ctx.fillStyle = '#94a3b8'
        }
    }, [assignment, seqLen, hbonds])

    return (
        <div className="space-y-8">
            {/* Canvas ribbon diagram */}
            <div className="bg-white rounded-3xl border border-slate-200 p-6 shadow-sm overflow-hidden">
                <div className="flex items-center justify-between mb-4">
                    <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Mapa de Estructura Secundaria</p>
                    <div className="flex items-center gap-4">
                        <div className="flex items-center gap-1.5"><div className="w-4 h-3 rounded-sm bg-rose-500/40 border border-rose-500" /><span className="text-[8px] font-bold text-slate-500">α-Hélice</span></div>
                        <div className="flex items-center gap-1.5"><div className="w-4 h-3 rounded-sm bg-blue-500/40 border border-blue-500" /><span className="text-[8px] font-bold text-slate-500">β-Lámina</span></div>
                        <div className="flex items-center gap-1.5"><div className="w-4 h-1 bg-slate-400 rounded" /><span className="text-[8px] font-bold text-slate-500">Ovillo</span></div>
                        <div className="flex items-center gap-1.5"><div className="w-4 h-px bg-amber-400" style={{ borderTop: '1px dashed #f59e0b' }} /><span className="text-[8px] font-bold text-slate-500">Enlace H</span></div>
                    </div>
                </div>
                <canvas ref={canvasRef} className="w-full" style={{ height: 120 }} />
            </div>

            {/* Stats cards */}
            {stats && (
                <div className="grid grid-cols-3 gap-4">
                    {[
                        { label: 'α-Hélices', count: stats.helix_count, pct: stats.helix_pct, residues: stats.helix_residues, color: '#e11d48', bg: '#fff1f2', border: '#fecdd3' },
                        { label: 'β-Láminas', count: stats.sheet_count, pct: stats.sheet_pct, residues: stats.sheet_residues, color: '#2563eb', bg: '#eff6ff', border: '#bfdbfe' },
                        { label: 'Ovillos', count: '—', pct: stats.coil_pct, residues: stats.coil_residues, color: '#64748b', bg: '#f8fafc', border: '#e2e8f0' },
                    ].map(s => (
                        <div key={s.label} className="rounded-2xl p-6 border text-center" style={{ backgroundColor: s.bg, borderColor: s.border }}>
                            <p className="text-4xl font-black mb-1" style={{ color: s.color }}>{s.pct}%</p>
                            <p className="text-[10px] font-black uppercase tracking-widest mb-1" style={{ color: s.color }}>{s.label}</p>
                            <p className="text-[9px] text-slate-400 font-bold">{s.residues} res{s.count !== '—' ? ` · ${s.count} seg` : ''}</p>
                        </div>
                    ))}
                </div>
            )}
        </div>
    )
}


// ══════════════════════════════════════════════════════════════
// MAIN PIPELINE COMPONENT
// ══════════════════════════════════════════════════════════════
export default function ProteinPipeline({ pipelineData, pdbUrl, structureSource }) {
    if (!pipelineData) return null

    const { primary, secondary, tertiary, quaternary, quality, gene_info } = pipelineData
    const [showActiveSites, setShowActiveSites] = useState(false)

    return (
        <div className="space-y-4 animate-in fade-in duration-1000">

            {/* ════════ LEVEL 1: PRIMARY STRUCTURE ════════ */}
            <section className="bg-white rounded-[2.5rem] border border-slate-200 p-10 shadow-sm animate-in slide-in-from-bottom-8 duration-700" style={{ animationDelay: '0ms' }}>
                <LevelHeader
                    number="1"
                    title="Estructura Primaria"
                    subtitle="La secuencia lineal de aminoácidos unidos por enlaces peptídicos. Define la identidad y todas las propiedades de la proteína."
                    color="#10b981"
                    icon="🔗"
                />

                <PrimaryBeadChain residues={primary.residue_properties} maxDisplay={100} />

                <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mt-8">
                    <div className="bg-emerald-50/50 rounded-2xl p-5 border border-emerald-100 text-center">
                        <p className="text-2xl font-black text-emerald-700">{primary.length}</p>
                        <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Aminoácidos</p>
                    </div>
                    <div className="bg-blue-50/50 rounded-2xl p-5 border border-blue-100 text-center">
                        <p className="text-2xl font-black text-blue-700">{primary.molecular_weight_kda}</p>
                        <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">kDa (peso mol.)</p>
                    </div>
                    <div className="bg-indigo-50/50 rounded-2xl p-5 border border-indigo-100 text-center">
                        <p className="text-2xl font-black text-indigo-700">{primary.net_charge > 0 ? '+' : ''}{primary.net_charge}</p>
                        <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Carga neta</p>
                    </div>
                    <div className="bg-slate-50/50 rounded-2xl p-5 border border-slate-100 text-center">
                        <p className="text-2xl font-black text-slate-700">{gene_info.gene_name || gene_info.locus_tag || '—'}</p>
                        <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Gen codificante</p>
                    </div>
                </div>

                {/* Functional Domains Visualization */}
                {pipelineData.functional?.conserved_domains?.length > 0 && (
                    <div className="mt-8 border-t border-slate-100 pt-8 animate-in fade-in duration-700 delay-300">
                        <div className="flex items-center gap-3 mb-4">
                            <div className="w-8 h-8 rounded-full bg-indigo-100 flex items-center justify-center text-indigo-600">
                                <span className="text-xs">🧩</span>
                            </div>
                            <div>
                                <h4 className="text-sm font-black text-slate-700">Dominios Conservados</h4>
                                <p className="text-[10px] text-slate-400 font-medium">Regiones funcionales identificadas en la secuencia</p>
                            </div>
                        </div>
                        <div className="relative h-12 bg-slate-100 rounded-xl overflow-hidden w-full">
                            {/* Base line */}
                            <div className="absolute top-1/2 left-0 right-0 h-0.5 bg-slate-300 -translate-y-1/2" />

                            {/* Domains */}
                            {pipelineData.functional.conserved_domains.map((domain, i) => {
                                const startPct = (domain.start / primary.length) * 100
                                const widthPct = ((domain.end - domain.start) / primary.length) * 100
                                const colors = ['#f472b6', '#34d399', '#60a5fa', '#a78bfa']
                                const color = colors[i % colors.length]

                                return (
                                    <div
                                        key={i}
                                        className="absolute top-2 h-8 rounded-lg flex items-center justify-center text-[9px] font-black text-white shadow-sm hover:scale-105 transition-transform cursor-help z-10"
                                        style={{
                                            left: `${startPct}%`,
                                            width: `${widthPct}%`,
                                            backgroundColor: color,
                                            minWidth: '20px'
                                        }}
                                        title={`${domain.name}: ${domain.start}-${domain.end}\n${domain.description}`}
                                    >
                                        <span className="truncate px-2">{domain.name}</span>
                                    </div>
                                )
                            })}
                        </div>
                    </div>
                )}
            </section>

            {/* Arrow: Primary → Secondary */}
            <PipelineArrow
                label="Enlaces de Hidrógeno"
                sublabel="Interacciones C=O···H-N del esqueleto peptídico"
                color="#e11d48"
                delay={200}
            />

            {/* ════════ LEVEL 2: SECONDARY STRUCTURE ════════ */}
            <section className="bg-white rounded-[2.5rem] border border-slate-200 shadow-sm overflow-hidden animate-in slide-in-from-bottom-8 duration-700" style={{ animationDelay: '300ms' }}>
                <div className="p-10">
                    <LevelHeader
                        number="2"
                        title="Estructura Secundaria"
                        subtitle="Patrones locales de plegamiento estabilizados por enlaces de hidrógeno entre el esqueleto peptídico: α-hélices (i→i+4) y β-láminas."
                        color="#e11d48"
                        icon="🌀"
                    />
                </div>

                {secondary ? (
                    <div className="grid grid-cols-1 lg:grid-cols-2 gap-0">
                        {/* Left: Data visualization */}
                        <div className="p-10 pt-0 space-y-6">
                            <SecondaryRibbon
                                assignment={secondary.assignment}
                                helices={secondary.helices}
                                sheets={secondary.sheets}
                                stats={secondary.stats}
                                hbonds={secondary.hydrogen_bonds}
                                sequence={primary.sequence}
                            />
                            {secondary.total_hbonds > 0 && (
                                <div className="bg-amber-50/50 rounded-2xl p-5 border border-amber-100 flex items-center gap-4">
                                    <div className="w-10 h-10 bg-amber-100 rounded-xl flex items-center justify-center text-amber-600 font-black text-sm">H</div>
                                    <div>
                                        <p className="text-sm font-black text-amber-800">{secondary.total_hbonds} enlaces de hidrógeno</p>
                                        <p className="text-[10px] text-amber-600 font-medium">Responsables de estabilizar hélices α y láminas β</p>
                                    </div>
                                </div>
                            )}
                        </div>
                        {/* Right: Mol* cartoon colored by SS type */}
                        {pdbUrl && (
                            <div className="border-l border-slate-100 bg-slate-900 relative">
                                <div className="absolute top-4 left-4 z-10 px-3 py-1.5 bg-white/90 rounded-full border border-slate-200 shadow-lg pointer-events-none">
                                    <span className="text-[8px] font-black text-slate-600 uppercase tracking-widest">Mol* — Cartoon / SS-Type</span>
                                </div>
                                <MolstarMini pdbUrl={pdbUrl} preset="cartoon-ss" height={400} />
                            </div>
                        )}
                    </div>
                ) : (
                    <div className="p-10 pt-0">
                        <div className="bg-slate-50 rounded-2xl p-8 text-center border border-slate-100">
                            <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">
                                Resuelve la estructura 3D primero para ver datos reales de estructura secundaria
                            </p>
                        </div>
                    </div>
                )}
            </section>

            {/* Arrow: Secondary → Tertiary */}
            <PipelineArrow
                label="Interacciones No Covalentes"
                sublabel="Puentes disulfuro, interacciones hidrofóbicas, puentes salinos, Van der Waals"
                color="#6366f1"
                delay={400}
            />

            {/* ════════ LEVEL 3: TERTIARY STRUCTURE ════════ */}
            <section className="bg-white rounded-[2.5rem] border border-slate-200 shadow-sm overflow-hidden animate-in slide-in-from-bottom-8 duration-700" style={{ animationDelay: '500ms' }}>
                <div className="p-10">
                    <LevelHeader
                        number="3"
                        title="Estructura Terciaria"
                        subtitle="Plegamiento 3D completo de la cadena polipeptídica. Estabilizada por puentes disulfuro (S-S), interacciones hidrofóbicas, puentes salinos y fuerzas de Van der Waals."
                        color="#6366f1"
                        icon="🧬"
                    />
                </div>

                {pdbUrl ? (
                    <div className="grid grid-cols-1 lg:grid-cols-3 gap-0">
                        {/* Mol* 3D viewer */}
                        <div className="lg:col-span-2 bg-slate-900 border-t border-slate-200 relative">
                            <div className="absolute top-4 left-4 z-10 px-3 py-1.5 bg-white/90 rounded-full border border-slate-200 shadow-lg pointer-events-none">
                                <span className="text-[8px] font-black text-slate-600 uppercase tracking-widest">Mol* — Estructura 3D Completa</span>
                            </div>
                            <MolstarMini
                                pdbUrl={pdbUrl}
                                preset="auto"
                                height={500}
                                highlightResidues={showActiveSites ? (pipelineData.functional?.active_sites?.map(s => s.residue) || []) : []}
                            />
                        </div>
                        {/* Info panel */}
                        <div className="border-l border-slate-100 p-8 space-y-6 bg-slate-50/30">
                            {tertiary && (
                                <>
                                    {/* Functional Sites Toggle */}
                                    {pipelineData.functional?.active_sites?.length > 0 && (
                                        <div className="bg-indigo-50 border border-indigo-100 rounded-2xl p-4">
                                            <div className="flex items-center justify-between mb-3">
                                                <div className="flex items-center gap-2">
                                                    <span className="text-lg">⚡</span>
                                                    <h4 className="text-xs font-black text-indigo-900 uppercase tracking-wider">Sitios Activos</h4>
                                                </div>
                                                <button
                                                    onClick={() => setShowActiveSites(!showActiveSites)}
                                                    className={`relative inline-flex h-5 w-9 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2 ${showActiveSites ? 'bg-indigo-600' : 'bg-slate-300'}`}
                                                >
                                                    <span
                                                        className={`${showActiveSites ? 'translate-x-5' : 'translate-x-1'} inline-block h-3 w-3 transform rounded-full bg-white transition-transform`}
                                                    />
                                                </button>
                                            </div>

                                            {showActiveSites && (
                                                <div className="space-y-2 animate-in fade-in slide-in-from-top-2 duration-300">
                                                    {pipelineData.functional.active_sites.map((site, i) => (
                                                        <div key={i} className="flex items-center gap-2 bg-white/60 rounded-lg p-2 text-[10px] border border-indigo-100/50">
                                                            <span className="px-1.5 py-0.5 bg-indigo-100 text-indigo-700 rounded font-mono font-bold">{site.aa}{site.residue}</span>
                                                            <span className="text-slate-600 truncate flex-1">{site.description}</span>
                                                        </div>
                                                    ))}
                                                </div>
                                            )}
                                            {!showActiveSites && (
                                                <p className="text-[9px] text-indigo-600/70 font-medium">
                                                    {pipelineData.functional.active_sites.length} sitios identificados. Activa para visualizar en 3D.
                                                </p>
                                            )}
                                        </div>
                                    )}

                                    <div className="space-y-4">
                                        <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Fuente</p>
                                        <div className="flex items-center gap-3">
                                            <div className="w-2.5 h-2.5 bg-emerald-500 rounded-full shadow-[0_0_8px_rgba(16,185,129,0.4)] animate-pulse" />
                                            <span className="text-xs font-black text-slate-700 uppercase">{tertiary.source}</span>
                                        </div>
                                        <div className="text-[10px] text-slate-500 font-medium space-y-1">
                                            <p>{tertiary.total_atoms?.toLocaleString()} átomos totales</p>
                                            <p>{tertiary.total_residues} residuos resueltos</p>
                                            {tertiary.resolution && <p>Resolución: {tertiary.resolution} Å</p>}
                                        </div>
                                    </div>

                                    {/* Disulfide bonds */}
                                    {tertiary.disulfide_bonds?.length > 0 && (
                                        <div className="space-y-3">
                                            <p className="text-[9px] font-black text-amber-600 uppercase tracking-widest">Puentes Disulfuro (S-S)</p>
                                            <div className="space-y-2">
                                                {tertiary.disulfide_bonds.slice(0, 6).map((ss, i) => (
                                                    <div key={i} className="flex items-center gap-2 text-[10px] font-mono">
                                                        <span className="px-2 py-1 bg-amber-50 border border-amber-200 rounded-lg font-black text-amber-700">Cys{ss.residue1}</span>
                                                        <span className="text-amber-400">⟷</span>
                                                        <span className="px-2 py-1 bg-amber-50 border border-amber-200 rounded-lg font-black text-amber-700">Cys{ss.residue2}</span>
                                                        {ss.predicted && <span className="text-[8px] text-slate-300 italic">pred</span>}
                                                    </div>
                                                ))}
                                            </div>
                                        </div>
                                    )}

                                    {/* Domains */}
                                    {tertiary.domains?.length > 0 && (
                                        <div className="space-y-3">
                                            <p className="text-[9px] font-black text-indigo-500 uppercase tracking-widest">Dominios Estructurales</p>
                                            <div className="space-y-2">
                                                {tertiary.domains.map((d, i) => (
                                                    <div key={i} className="bg-indigo-50/50 border border-indigo-100 rounded-xl p-3">
                                                        <p className="text-[10px] font-black text-indigo-700">{d.name}</p>
                                                        <p className="text-[9px] text-slate-400 font-mono">{d.start}–{d.end} ({d.end - d.start + 1} aa)</p>
                                                    </div>
                                                ))}
                                            </div>
                                        </div>
                                    )}

                                    {/* Quality / pLDDT */}
                                    {quality && (
                                        <div className="bg-gradient-to-br from-indigo-600 to-blue-700 rounded-2xl p-5 text-white">
                                            <p className="text-[9px] font-black uppercase tracking-widest opacity-60 mb-2">
                                                {quality.is_plddt ? 'Confianza pLDDT' : 'B-Factor promedio'}
                                            </p>
                                            <p className="text-3xl font-black">{quality.avg_confidence}</p>
                                            <p className="text-[9px] opacity-60 mt-1">
                                                {quality.is_plddt
                                                    ? quality.avg_confidence > 70 ? 'Alta confianza' : 'Confianza moderada'
                                                    : 'Flexibilidad molecular'}
                                            </p>
                                        </div>
                                    )}
                                </>
                            )}
                        </div>
                    </div>
                ) : (
                    <div className="p-10 pt-0">
                        <div className="bg-slate-50 rounded-2xl p-8 text-center border border-slate-100">
                            <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">
                                Resuelve la estructura 3D para visualizar el plegamiento terciario
                            </p>
                        </div>
                    </div>
                )}
            </section>

            {/* Arrow: Tertiary → Quaternary */}
            <PipelineArrow
                label="Ensamblaje de Subunidades"
                sublabel="Interacciones proteína-proteína entre cadenas polipeptídicas"
                color="#0891b2"
                delay={600}
            />

            {/* ════════ LEVEL 4: QUATERNARY STRUCTURE ════════ */}
            <section className="bg-white rounded-[2.5rem] border border-slate-200 shadow-sm overflow-hidden animate-in slide-in-from-bottom-8 duration-700" style={{ animationDelay: '700ms' }}>
                <div className="p-10">
                    <LevelHeader
                        number="4"
                        title="Estructura Cuaternaria"
                        subtitle="Ensamblaje de múltiples cadenas polipeptídicas (subunidades) formando un complejo proteico funcional. Ej: hemoglobina (α₂β₂)."
                        color="#0891b2"
                        icon="🏗️"
                    />
                </div>

                {quaternary ? (
                    <div className="grid grid-cols-1 lg:grid-cols-2 gap-0">
                        {/* Left: Chain info */}
                        <div className="p-10 pt-0 space-y-6">
                            {/* Oligomeric state banner */}
                            <div className="bg-gradient-to-r from-slate-900 via-slate-800 to-cyan-900 rounded-2xl p-8 text-white relative overflow-hidden">
                                <div className="absolute top-0 right-0 w-40 h-40 bg-cyan-500/10 blur-[60px] rounded-full" />
                                <div className="relative z-10 flex items-center justify-between">
                                    <div>
                                        <p className="text-[9px] font-black uppercase tracking-[0.3em] text-cyan-300/70 mb-1">Estado Oligomérico</p>
                                        <h4 className="text-3xl font-black uppercase tracking-tighter">{quaternary.oligomeric_state}</h4>
                                        <p className="text-[10px] text-slate-300 mt-2">
                                            {quaternary.is_complex
                                                ? `Complejo con ${quaternary.num_chains} subunidades identificadas en la estructura PDB.`
                                                : 'Esta proteína funciona como una unidad monomérica individual.'}
                                        </p>
                                    </div>
                                    <div className="text-right">
                                        <p className="text-5xl font-black text-cyan-400/40">{quaternary.num_chains}</p>
                                        <p className="text-[8px] font-black uppercase tracking-widest text-slate-400">cadena{quaternary.num_chains > 1 ? 's' : ''}</p>
                                    </div>
                                </div>
                            </div>

                            {/* Chain details */}
                            <div className="space-y-3">
                                {quaternary.chains?.map((chain, i) => {
                                    const chainColors = ['#3b82f6', '#e11d48', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899']
                                    const c = chainColors[i % chainColors.length]
                                    return (
                                        <div key={chain.chain_id} className="flex items-center gap-4 bg-slate-50/50 rounded-xl p-4 border border-slate-100 hover:border-cyan-200 transition-all">
                                            <div className="w-10 h-10 rounded-xl flex items-center justify-center text-white font-black text-sm" style={{ backgroundColor: c }}>
                                                {chain.chain_id}
                                            </div>
                                            <div className="flex-1">
                                                <p className="text-xs font-black text-slate-800">Cadena {chain.chain_id}</p>
                                                <p className="text-[9px] text-slate-400 font-medium">{chain.residue_count} res · {chain.atom_count?.toLocaleString()} átomos · Rango {chain.min_residue}–{chain.max_residue}</p>
                                            </div>
                                            <div className="text-right">
                                                <p className="text-sm font-black" style={{ color: chain.avg_bfactor > 70 ? '#16a34a' : chain.avg_bfactor > 50 ? '#ca8a04' : '#dc2626' }}>
                                                    {chain.avg_bfactor}
                                                </p>
                                                <p className="text-[8px] text-slate-400 font-bold">{quality?.is_plddt ? 'pLDDT' : 'B-factor'}</p>
                                            </div>
                                        </div>
                                    )
                                })}
                            </div>

                            {/* Note for monomers */}
                            {!quaternary.is_complex && (
                                <div className="bg-cyan-50/50 rounded-2xl p-5 border border-cyan-100">
                                    <p className="text-[10px] text-cyan-700 font-medium leading-relaxed">
                                        <strong>Nota:</strong> La estructura cuaternaria aplica a proteínas con múltiples subunidades
                                        (ej: hemoglobina α₂β₂, colágeno triple hélice, cápsides virales). Esta proteína es monomérica
                                        en la estructura resuelta.
                                    </p>
                                </div>
                            )}
                        </div>

                        {/* Right: Mol* chain-colored view */}
                        {pdbUrl && (
                            <div className="border-l border-slate-100 bg-slate-900 relative">
                                <div className="absolute top-4 left-4 z-10 px-3 py-1.5 bg-white/90 rounded-full border border-slate-200 shadow-lg pointer-events-none">
                                    <span className="text-[8px] font-black text-slate-600 uppercase tracking-widest">Mol* — Coloreado por Cadena</span>
                                </div>
                                <MolstarMini pdbUrl={pdbUrl} preset="cartoon-chain" height={400} />
                            </div>
                        )}
                    </div>
                ) : (
                    <div className="p-10 pt-0">
                        <div className="bg-slate-50 rounded-2xl p-8 text-center border border-slate-100">
                            <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">
                                Resuelve la estructura 3D para analizar la organización cuaternaria
                            </p>
                        </div>
                    </div>
                )}
            </section>

            {/* Summary footer */}
            <div className="bg-gradient-to-r from-slate-900 via-slate-800 to-blue-900 rounded-3xl p-10 text-white relative overflow-hidden shadow-2xl">
                <div className="absolute inset-0 bg-[radial-gradient(circle_at_30%_50%,rgba(59,130,246,0.1),transparent_70%)]" />
                <div className="relative z-10 text-center space-y-4">
                    <p className="text-[9px] font-black uppercase tracking-[0.4em] text-blue-300/60">Resumen de Organización Estructural</p>
                    <div className="flex items-center justify-center gap-6 flex-wrap">
                        <div className="flex items-center gap-3">
                            <div className="w-3 h-3 bg-emerald-500 rounded-full" />
                            <span className="text-xs font-bold">{primary.length} aa</span>
                        </div>
                        {secondary && (
                            <div className="flex items-center gap-3">
                                <div className="w-3 h-3 bg-rose-500 rounded-full" />
                                <span className="text-xs font-bold">{secondary.stats?.helix_pct}% hélice · {secondary.stats?.sheet_pct}% lámina</span>
                            </div>
                        )}
                        {tertiary && (
                            <div className="flex items-center gap-3">
                                <div className="w-3 h-3 bg-indigo-500 rounded-full" />
                                <span className="text-xs font-bold">{tertiary.total_atoms?.toLocaleString()} átomos</span>
                            </div>
                        )}
                        {quaternary && (
                            <div className="flex items-center gap-3">
                                <div className="w-3 h-3 bg-cyan-500 rounded-full" />
                                <span className="text-xs font-bold">{quaternary.oligomeric_state}</span>
                            </div>
                        )}
                    </div>
                    <p className="text-[10px] text-slate-300 max-w-2xl mx-auto leading-relaxed">
                        {gene_info.product || 'Proteína hipotética'} — codificada por el gen <strong>{gene_info.gene_name || gene_info.locus_tag}</strong>
                        {gene_info.start ? ` (posición ${gene_info.start}–${gene_info.end})` : ''}
                    </p>
                </div>
            </div>
        </div>
    )
}
