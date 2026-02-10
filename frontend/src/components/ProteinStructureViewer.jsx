/**
 * ProteinStructureViewer Component
 * Visualización completa de estructuras de proteínas:
 * - Estructura Primaria: Secuencia de aminoácidos
 * - Estructura Secundaria: Hélices alfa, láminas beta, loops
 * - Estructura Terciaria: Estructura 3D completa
 * - Estructura Cuaternaria: Complejos multiproteicos
 */
import { useState, useEffect, useRef } from 'react'
import MolstarProteinViewer from './MolstarProteinViewer'

// Propiedades de aminoácidos
const AA_PROPERTIES = {
    A: { name: 'Alanina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    I: { name: 'Isoleucina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    L: { name: 'Leucina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    M: { name: 'Metionina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    F: { name: 'Fenilalanina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    W: { name: 'Triptófano', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    V: { name: 'Valina', group: 'Hidrofóbico', color: '#4ade80', type: 'nonpolar' },
    P: { name: 'Prolina', group: 'Especial', color: '#fbbf24', type: 'special' },
    G: { name: 'Glicina', group: 'Especial', color: '#fbbf24', type: 'special' },
    C: { name: 'Cisteína', group: 'Especial', color: '#fbbf24', type: 'sulfur' },
    D: { name: 'Ácido aspártico', group: 'Cargado (−)', color: '#f87171', type: 'acidic' },
    E: { name: 'Ácido glutámico', group: 'Cargado (−)', color: '#f87171', type: 'acidic' },
    K: { name: 'Lisina', group: 'Cargado (+)', color: '#fb923c', type: 'basic' },
    R: { name: 'Arginina', group: 'Cargado (+)', color: '#fb923c', type: 'basic' },
    H: { name: 'Histidina', group: 'Cargado (+)', color: '#fb923c', type: 'basic' },
    N: { name: 'Asparagina', group: 'Polar', color: '#60a5fa', type: 'polar' },
    Q: { name: 'Glutamina', group: 'Polar', color: '#60a5fa', type: 'polar' },
    S: { name: 'Serina', group: 'Polar', color: '#60a5fa', type: 'polar' },
    T: { name: 'Treonina', group: 'Polar', color: '#60a5fa', type: 'polar' },
    Y: { name: 'Tirosina', group: 'Polar', color: '#60a5fa', type: 'polar' },
    '*': { name: 'Stop', group: 'Stop', color: '#ef4444', type: 'stop' },
}

// Predicción simple de estructura secundaria basada en propensión de aminoácidos
const predictSecondaryStructure = (sequence) => {
    // Propensión de aminoácidos para formar hélices alfa (basado en datos experimentales)
    const helixPropensity = { A: 1.45, E: 1.53, L: 1.34, M: 1.20, K: 1.07, F: 1.12, Q: 1.17, W: 1.14, I: 1.00, V: 0.90, D: 0.98, H: 1.24, R: 0.79, T: 0.82, S: 0.79, C: 0.77, Y: 0.61, N: 0.73, G: 0.53, P: 0.59 }

    // Propensión para formar láminas beta
    const sheetPropensity = { V: 1.87, I: 1.60, F: 1.43, Y: 1.39, W: 1.19, L: 1.22, T: 1.20, C: 1.30, M: 1.07, Q: 0.98, R: 0.90, H: 0.80, A: 0.97, N: 0.65, G: 0.81, K: 0.74, S: 0.72, E: 0.26, D: 0.45, P: 0.62 }

    const structure = []
    const windowSize = 7

    for (let i = 0; i < sequence.length; i++) {
        let helixScore = 0
        let sheetScore = 0
        let count = 0

        // Ventana deslizante para suavizar predicción
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

        // Determinar estructura
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

export default function ProteinStructureViewer({ protein, onClose }) {
    const [activeTab, setActiveTab] = useState('primary')
    const [secondaryStructure, setSecondaryStructure] = useState(null)
    const [showLegend, setShowLegend] = useState(true)
    const viewerRef = useRef(null)
    const viewer3DRef = useRef(null)

    useEffect(() => {
        if (protein?.full_sequence) {
            // Predecir estructura secundaria
            const structure = predictSecondaryStructure(protein.full_sequence)
            setSecondaryStructure(structure)
        }
    }, [protein])

    // Análisis de composición
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
                .slice(0, 6)
                .map(([aa, count]) => ({
                    aa,
                    count,
                    pct: ((count / totalAA) * 100).toFixed(1),
                    ...(AA_PROPERTIES[aa] || { name: aa, group: 'Otro', color: '#94a3b8' })
                }))
        }
    }

    const aaComp = composition()

    // Renderizar estructura secundaria
    const renderSecondaryStructure = () => {
        if (!secondaryStructure || !protein?.full_sequence) return null

        // Agrupar estructuras consecutivas del mismo tipo
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

        return (
            <div className="space-y-4">
                <div className="bg-slate-50 rounded-lg p-4">
                    <h4 className="text-sm font-semibold text-slate-700 mb-2">Predicción de Estructura Secundaria</h4>
                    <p className="text-xs text-slate-600 mb-4">
                        Basado en propensión de aminoácidos. Las regiones en rojo son hélices alfa, azul son láminas beta, y gris son loops/giros.
                    </p>

                    {/* Visualización lineal mejorada */}
                    <div className="relative h-16 bg-white rounded-lg border-2 border-slate-200 overflow-hidden shadow-sm">
                        {regions.map((region, i) => {
                            const width = ((region.end - region.start + 1) / protein.full_sequence.length) * 100
                            const left = (region.start / protein.full_sequence.length) * 100
                            const colors = {
                                helix: 'bg-gradient-to-r from-red-400 to-red-500',
                                sheet: 'bg-gradient-to-r from-blue-400 to-blue-500',
                                loop: 'bg-gradient-to-r from-slate-300 to-slate-400'
                            }
                            const borderColors = {
                                helix: 'border-red-600',
                                sheet: 'border-blue-600',
                                loop: 'border-slate-500'
                            }
                            return (
                                <div
                                    key={i}
                                    className={`absolute h-full ${colors[region.type] || 'bg-slate-300'} border-r ${borderColors[region.type]} transition-all hover:brightness-110 cursor-pointer`}
                                    style={{ left: `${left}%`, width: `${width}%` }}
                                    title={`${region.type === 'helix' ? 'Hélice α' : region.type === 'sheet' ? 'Lámina β' : 'Loop'} (${region.start + 1}-${region.end + 1}) - ${region.end - region.start + 1} aa`}
                                />
                            )
                        })}
                    </div>

                    {/* Gráfico de distribución por posición */}
                    <div className="mt-4">
                        <h5 className="text-xs font-semibold text-slate-700 mb-2">Distribución a lo Largo de la Secuencia</h5>
                        <div className="bg-white rounded-lg border border-slate-200 p-3">
                            <svg width="100%" height="80" className="overflow-visible">
                                {secondaryStructure.map((ss, i) => {
                                    const x = (i / protein.full_sequence.length) * 100
                                    const heights = { helix: 60, sheet: 40, loop: 20 }
                                    const colors = { helix: '#f87171', sheet: '#60a5fa', loop: '#cbd5e1' }
                                    return (
                                        <rect
                                            key={i}
                                            x={`${x}%`}
                                            y={80 - heights[ss.type]}
                                            width={`${100 / protein.full_sequence.length}%`}
                                            height={heights[ss.type]}
                                            fill={colors[ss.type]}
                                            opacity="0.7"
                                        >
                                            <title>{`Posición ${i + 1}: ${ss.type === 'helix' ? 'Hélice α' : ss.type === 'sheet' ? 'Lámina β' : 'Loop'}`}</title>
                                        </rect>
                                    )
                                })}
                                {/* Eje horizontal */}
                                <line x1="0" y1="80" x2="100%" y2="80" stroke="#94a3b8" strokeWidth="1" />
                                {/* Etiquetas */}
                                <text x="2%" y="75" fontSize="10" fill="#64748b">1</text>
                                <text x="48%" y="75" fontSize="10" fill="#64748b" textAnchor="middle">{Math.floor(protein.full_sequence.length / 2)}</text>
                                <text x="96%" y="75" fontSize="10" fill="#64748b" textAnchor="end">{protein.full_sequence.length}</text>
                            </svg>
                        </div>
                    </div>

                    {/* Estadísticas */}
                    <div className="grid grid-cols-3 gap-3 mt-4">
                        <div className="bg-white rounded-lg p-3 border border-slate-200">
                            <div className="flex items-center gap-2 mb-1">
                                <div className="w-3 h-3 bg-red-400 rounded"></div>
                                <span className="text-xs font-medium text-slate-700">Hélices α</span>
                            </div>
                            <p className="text-lg font-bold text-slate-800">
                                {((regions.filter(r => r.type === 'helix').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / protein.full_sequence.length) * 100).toFixed(1)}%
                            </p>
                        </div>
                        <div className="bg-white rounded-lg p-3 border border-slate-200">
                            <div className="flex items-center gap-2 mb-1">
                                <div className="w-3 h-3 bg-blue-400 rounded"></div>
                                <span className="text-xs font-medium text-slate-700">Láminas β</span>
                            </div>
                            <p className="text-lg font-bold text-slate-800">
                                {((regions.filter(r => r.type === 'sheet').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / protein.full_sequence.length) * 100).toFixed(1)}%
                            </p>
                        </div>
                        <div className="bg-white rounded-lg p-3 border border-slate-200">
                            <div className="flex items-center gap-2 mb-1">
                                <div className="w-3 h-3 bg-slate-300 rounded"></div>
                                <span className="text-xs font-medium text-slate-700">Loops/Giros</span>
                            </div>
                            <p className="text-lg font-bold text-slate-800">
                                {((regions.filter(r => r.type === 'loop').reduce((sum, r) => sum + (r.end - r.start + 1), 0) / protein.full_sequence.length) * 100).toFixed(1)}%
                            </p>
                        </div>
                    </div>

                    {/* Lista de regiones */}
                    <div className="mt-4">
                        <h5 className="text-xs font-semibold text-slate-700 mb-2">Regiones Detectadas</h5>
                        <div className="max-h-48 overflow-y-auto space-y-1">
                            {regions.filter(r => r.end - r.start >= 3).map((region, i) => (
                                <div key={i} className="flex items-center gap-2 text-xs bg-white p-2 rounded border border-slate-100">
                                    <div className={`w-2 h-2 rounded-full ${region.type === 'helix' ? 'bg-red-400' : region.type === 'sheet' ? 'bg-blue-400' : 'bg-slate-300'}`}></div>
                                    <span className="font-medium text-slate-700">
                                        {region.type === 'helix' ? 'Hélice α' : region.type === 'sheet' ? 'Lámina β' : 'Loop'}
                                    </span>
                                    <span className="text-slate-500">
                                        Residuos {region.start + 1}-{region.end + 1}
                                    </span>
                                    <span className="text-slate-400">
                                        ({region.end - region.start + 1} aa)
                                    </span>
                                </div>
                            ))}
                        </div>
                    </div>
                </div>

                {/* Leyenda de estructuras */}
                <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                    <h5 className="text-sm font-semibold text-blue-900 mb-2">Información sobre Estructuras Secundarias</h5>
                    <div className="space-y-2 text-xs text-blue-800">
                        <div className="flex items-start gap-2">
                            <div className="w-2 h-2 bg-red-400 rounded-full mt-1 flex-shrink-0"></div>
                            <div>
                                <span className="font-medium">Hélices α:</span> Estructura estabilizada por puentes de hidrógeno entre C=O y N-H. Común en proteínas estructurales.
                            </div>
                        </div>
                        <div className="flex items-start gap-2">
                            <div className="w-2 h-2 bg-blue-400 rounded-full mt-1 flex-shrink-0"></div>
                            <div>
                                <span className="font-medium">Láminas β:</span> Cadenas extendidas con puentes de hidrógeno laterales. Común en enzimas y proteínas de unión.
                            </div>
                        </div>
                        <div className="flex items-start gap-2">
                            <div className="w-2 h-2 bg-slate-300 rounded-full mt-1 flex-shrink-0"></div>
                            <div>
                                <span className="font-medium">Loops/Giros:</span> Regiones flexibles que conectan hélices y láminas. Importantes para función y movimiento.
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        )
    }

    if (!protein) return null

    return (
        <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
            {/* Header */}
            <div className="bg-gradient-to-r from-slate-800 to-teal-800 px-6 py-4 text-white">
                <div className="flex items-start justify-between">
                    <div>
                        <h3 className="font-bold text-lg">{protein.protein_id || protein.locus_tag}</h3>
                        <p className="text-teal-200 text-sm mt-1">{protein.product || 'Proteína hipotética'}</p>
                        {protein.gene_name && (
                            <p className="text-teal-300 text-xs mt-1">Gen: {protein.gene_name}</p>
                        )}
                    </div>
                    <button onClick={onClose} className="text-slate-300 hover:text-white transition-colors">
                        <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                        </svg>
                    </button>
                </div>
            </div>

            {/* Tabs */}
            <div className="border-b border-slate-200 bg-slate-50 px-6">
                <div className="flex gap-1 overflow-x-auto">
                    {[
                        { id: 'primary', label: 'Estructura Primaria', icon: '🔗' },
                        { id: 'secondary', label: 'Estructura Secundaria', icon: '🌀' },
                        { id: 'tertiary', label: 'Estructura Terciaria (3D)', icon: '🧬' },
                        { id: 'composition', label: 'Composición', icon: '📊' },
                        { id: 'properties', label: 'Propiedades', icon: '⚗️' },
                    ].map(tab => (
                        <button
                            key={tab.id}
                            onClick={() => setActiveTab(tab.id)}
                            className={`px-4 py-3 text-sm font-medium transition-colors border-b-2 ${
                                activeTab === tab.id
                                    ? 'border-teal-600 text-teal-700'
                                    : 'border-transparent text-slate-600 hover:text-slate-800'
                            }`}
                        >
                            <span className="mr-2">{tab.icon}</span>
                            {tab.label}
                        </button>
                    ))}
                </div>
            </div>

            {/* Content */}
            <div className="p-6 max-h-[600px] overflow-y-auto">
                {activeTab === 'primary' && (
                    <div className="space-y-4">
                        <div>
                            <h4 className="text-sm font-semibold text-slate-700 mb-2">Secuencia de Aminoácidos</h4>
                            <p className="text-xs text-slate-600 mb-3">
                                Secuencia primaria de {protein.full_sequence?.length || 0} aminoácidos.
                                Coloreados por grupo químico.
                            </p>
                            <div className="bg-slate-900 rounded-lg p-4 overflow-x-auto">
                                <p className="font-mono text-sm break-all leading-loose tracking-wide">
                                    {protein.full_sequence?.split('').map((aa, i) => {
                                        const prop = AA_PROPERTIES[aa]
                                        return (
                                            <span
                                                key={i}
                                                style={{ color: prop?.color || '#94a3b8' }}
                                                title={`${aa}${i + 1} — ${prop?.name || aa} (${prop?.group || '?'})`}
                                                className="hover:bg-slate-800 cursor-help"
                                            >
                                                {aa}
                                            </span>
                                        )
                                    })}
                                </p>
                            </div>

                            {/* Leyenda */}
                            {showLegend && (
                                <div className="flex flex-wrap gap-3 mt-3 text-xs">
                                    <span className="flex items-center gap-1.5">
                                        <span className="w-3 h-3 rounded-full" style={{ backgroundColor: '#4ade80' }}></span>
                                        <span className="text-slate-600">Hidrofóbico (A,I,L,M,F,W,V)</span>
                                    </span>
                                    <span className="flex items-center gap-1.5">
                                        <span className="w-3 h-3 rounded-full" style={{ backgroundColor: '#60a5fa' }}></span>
                                        <span className="text-slate-600">Polar (N,Q,S,T,Y)</span>
                                    </span>
                                    <span className="flex items-center gap-1.5">
                                        <span className="w-3 h-3 rounded-full" style={{ backgroundColor: '#fb923c' }}></span>
                                        <span className="text-slate-600">Cargado + (K,R,H)</span>
                                    </span>
                                    <span className="flex items-center gap-1.5">
                                        <span className="w-3 h-3 rounded-full" style={{ backgroundColor: '#f87171' }}></span>
                                        <span className="text-slate-600">Cargado − (D,E)</span>
                                    </span>
                                    <span className="flex items-center gap-1.5">
                                        <span className="w-3 h-3 rounded-full" style={{ backgroundColor: '#fbbf24' }}></span>
                                        <span className="text-slate-600">Especial (G,C,P)</span>
                                    </span>
                                </div>
                            )}
                        </div>

                        {/* Información de enlaces peptídicos */}
                        <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                            <h5 className="text-sm font-semibold text-blue-900 mb-3">Enlaces Peptídicos</h5>
                            <p className="text-xs text-blue-800 mb-3">
                                La estructura primaria está formada por <span className="font-bold">enlaces peptídicos</span> entre
                                el grupo carboxilo (-COOH) de un aminoácido y el grupo amino (-NH₂) del siguiente,
                                liberando una molécula de agua (H₂O).
                            </p>

                            <div className="bg-white rounded-lg p-3 border border-blue-200">
                                <div className="mb-2">
                                    <span className="text-xs font-semibold text-blue-900">Total de enlaces peptídicos:</span>
                                    <span className="ml-2 text-sm font-bold text-blue-700">
                                        {(protein.full_sequence?.length || 1) - 1} enlaces
                                    </span>
                                </div>

                                <div className="mt-3 space-y-2">
                                    <p className="text-xs text-blue-700 font-medium">Reacción de condensación:</p>
                                    <div className="text-xs font-mono bg-blue-50 rounded p-2 overflow-x-auto">
                                        <div className="whitespace-nowrap text-blue-800">
                                            H₂N-CHR₁-COOH + H₂N-CHR₂-COOH → H₂N-CHR₁-CO-NH-CHR₂-COOH + H₂O
                                        </div>
                                    </div>

                                    <p className="text-xs text-blue-700 font-medium mt-3">Cadena polipeptídica resultante:</p>
                                    <div className="text-xs font-mono bg-blue-50 rounded p-2 overflow-x-auto">
                                        <div className="whitespace-nowrap text-blue-800">
                                            H₂N—CHR₁—<span className="font-bold text-red-600">CO—NH</span>—CHR₂—<span className="font-bold text-red-600">CO—NH</span>—CHR₃—<span className="font-bold text-red-600">CO—NH</span>—CHRₙ—COOH
                                        </div>
                                    </div>

                                    <div className="mt-3 flex items-start gap-2 bg-amber-50 border border-amber-200 rounded p-2">
                                        <div className="text-amber-600 font-bold text-xs mt-0.5">⚡</div>
                                        <p className="text-xs text-amber-800">
                                            <span className="font-medium">Enlaces destacados en rojo:</span> El enlace peptídico
                                            <span className="font-mono font-bold"> CO—NH </span> es planar y rígido debido a la
                                            resonancia del orbital π. Esta rigidez limita la rotación y determina la conformación de la proteína.
                                        </p>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                )}

                {activeTab === 'secondary' && renderSecondaryStructure()}

                {activeTab === 'tertiary' && (
                    <div className="space-y-4">
                        <div>
                            <h4 className="text-sm font-semibold text-slate-700 mb-2">Visualización 3D de la Estructura Terciaria</h4>
                            <p className="text-xs text-slate-600 mb-4">
                                Visualizador Mol* de última generación con integración de PDB y AlphaFold.
                                Busca la proteína por nombre, carga un PDB ID directamente, o visualiza estructuras predichas.
                            </p>
                            <MolstarProteinViewer
                                proteinId={protein.protein_id}
                                proteinSequence={protein.full_sequence}
                                productName={protein.product}
                            />
                        </div>

                        {/* Información educativa */}
                        <div className="bg-purple-50 border border-purple-200 rounded-lg p-4">
                            <h5 className="text-sm font-semibold text-purple-900 mb-2">Sobre la Estructura Terciaria</h5>
                            <div className="space-y-2 text-xs text-purple-800">
                                <p>
                                    <span className="font-medium">Estructura Terciaria:</span> Es el plegamiento tridimensional
                                    completo de una cadena polipeptídica. Esta estructura determina la función biológica de la proteína.
                                </p>
                                <p>
                                    <span className="font-medium">Fuerzas estabilizadoras:</span> Puentes disulfuro (S-S),
                                    interacciones hidrofóbicas, puentes de hidrógeno, fuerzas de van der Waals, e interacciones iónicas.
                                </p>
                                <p>
                                    <span className="font-medium">Fuentes de datos:</span> PDB contiene estructuras experimentales
                                    determinadas por cristalografía de rayos X, RMN o criomicroscopía electrónica. AlphaFold proporciona
                                    predicciones de alta precisión basadas en aprendizaje profundo.
                                </p>
                            </div>
                        </div>
                    </div>
                )}

                {activeTab === 'composition' && aaComp && (
                    <div className="space-y-4">
                        <div>
                            <h4 className="text-sm font-semibold text-slate-700 mb-3">Composición por Grupos</h4>
                            <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
                                {aaComp.groups.map(g => (
                                    <div key={g.group} className="bg-slate-50 rounded-lg p-3 border border-slate-200">
                                        <p className="text-xs text-slate-600 mb-1">{g.group}</p>
                                        <p className="text-2xl font-bold text-slate-800">{g.pct}%</p>
                                        <div className="w-full h-2 bg-slate-200 rounded-full mt-2">
                                            <div className="h-full rounded-full transition-all"
                                                style={{
                                                    width: `${g.pct}%`,
                                                    backgroundColor: g.group === 'Hidrofóbico' ? '#4ade80'
                                                        : g.group === 'Polar' ? '#60a5fa'
                                                            : g.group.includes('+') ? '#fb923c'
                                                                : g.group.includes('−') ? '#f87171'
                                                                    : '#fbbf24'
                                                }}
                                            />
                                        </div>
                                        <p className="text-xs text-slate-500 mt-1">{g.count} residuos</p>
                                    </div>
                                ))}
                            </div>
                        </div>

                        <div>
                            <h4 className="text-sm font-semibold text-slate-700 mb-3">Top Aminoácidos</h4>
                            <div className="grid grid-cols-2 sm:grid-cols-3 gap-2">
                                {aaComp.topAA.map(aa => (
                                    <div key={aa.aa} className="bg-white border border-slate-200 rounded-lg p-3">
                                        <div className="flex items-center justify-between mb-1">
                                            <span className="font-mono text-xl font-bold" style={{ color: aa.color }}>{aa.aa}</span>
                                            <span className="text-lg font-bold text-slate-800">{aa.pct}%</span>
                                        </div>
                                        <p className="text-xs text-slate-600">{aa.name}</p>
                                        <p className="text-xs text-slate-500">{aa.count} residuos</p>
                                    </div>
                                ))}
                            </div>
                        </div>
                    </div>
                )}

                {activeTab === 'properties' && (
                    <div className="space-y-4">
                        <div className="grid grid-cols-2 sm:grid-cols-3 gap-4">
                            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
                                <span className="text-xs text-slate-600">Longitud</span>
                                <p className="text-2xl font-bold text-slate-800">{protein.length} aa</p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
                                <span className="text-xs text-slate-600">Peso Molecular</span>
                                <p className="text-2xl font-bold text-slate-800">
                                    {(protein.molecular_weight_approx / 1000).toFixed(1)} kDa
                                </p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
                                <span className="text-xs text-slate-600">Posición Genómica</span>
                                <p className="text-sm font-mono font-medium text-slate-800">
                                    {protein.start?.toLocaleString()} - {protein.end?.toLocaleString()}
                                </p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200">
                                <span className="text-xs text-slate-600">Longitud DNA</span>
                                <p className="text-2xl font-bold text-slate-800">
                                    {protein.start && protein.end
                                        ? `${Math.abs(protein.end - protein.start).toLocaleString()} bp`
                                        : '—'}
                                </p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-4 border border-slate-200 col-span-2">
                                <span className="text-xs text-slate-600">Hebra / Dirección</span>
                                <p className="font-medium mt-1">
                                    <span className={`inline-flex items-center gap-2 px-3 py-1.5 rounded-lg ${protein.strand === 1
                                        ? 'bg-teal-100 text-teal-700'
                                        : 'bg-violet-100 text-violet-700'
                                        }`}>
                                        {protein.strand === 1 ? '→ 5\'→3\' Forward' : '← 3\'→5\' Reverse'}
                                    </span>
                                </p>
                            </div>
                        </div>

                        {/* Enlaces externos */}
                        <div>
                            <h4 className="text-sm font-semibold text-slate-700 mb-3">Enlaces Externos</h4>
                            <div className="flex flex-wrap gap-2">
                                {protein.protein_id && (
                                    <a href={`https://www.ncbi.nlm.nih.gov/protein/${protein.protein_id}`}
                                        target="_blank" rel="noopener noreferrer"
                                        className="inline-flex items-center gap-2 px-4 py-2 bg-teal-50 border border-teal-200 rounded-lg text-sm text-teal-700 hover:bg-teal-100 font-medium transition-colors">
                                        NCBI Protein
                                    </a>
                                )}
                                {protein.protein_id && (
                                    <a href={`https://www.uniprot.org/uniprotkb?query=${protein.protein_id}`}
                                        target="_blank" rel="noopener noreferrer"
                                        className="inline-flex items-center gap-2 px-4 py-2 bg-purple-50 border border-purple-200 rounded-lg text-sm text-purple-700 hover:bg-purple-100 font-medium transition-colors">
                                        UniProt
                                    </a>
                                )}
                                {protein.protein_id && (
                                    <a href={`https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22value%22%3A%22${protein.protein_id}%22%7D%7D%7D`}
                                        target="_blank" rel="noopener noreferrer"
                                        className="inline-flex items-center gap-2 px-4 py-2 bg-blue-50 border border-blue-200 rounded-lg text-sm text-blue-700 hover:bg-blue-100 font-medium transition-colors">
                                        PDB (Estructura 3D)
                                    </a>
                                )}
                                {protein.protein_id && (
                                    <a href={`https://alphafold.ebi.ac.uk/search/text/${protein.protein_id}`}
                                        target="_blank" rel="noopener noreferrer"
                                        className="inline-flex items-center gap-2 px-4 py-2 bg-green-50 border border-green-200 rounded-lg text-sm text-green-700 hover:bg-green-100 font-medium transition-colors">
                                        AlphaFold
                                    </a>
                                )}
                            </div>
                        </div>

                        {/* Nota sobre estructura terciaria */}
                        <div className="bg-teal-50 border border-teal-200 rounded-lg p-4">
                            <h5 className="text-sm font-semibold text-teal-900 mb-2">Visualización 3D Integrada</h5>
                            <p className="text-xs text-teal-800 mb-3">
                                Hemos integrado un visualizador 3D interactivo para esta proteína. Puedes explorar su estructura
                                terciaria con controles interactivos y diferentes estilos de visualización.
                            </p>
                            <button
                                onClick={() => setActiveTab('tertiary')}
                                className="inline-flex items-center gap-2 px-4 py-2 bg-teal-600 text-white rounded-lg text-sm hover:bg-teal-700 font-medium transition-colors"
                            >
                                Ver Estructura 3D
                            </button>
                        </div>
                    </div>
                )}
            </div>
        </div>
    )
}
