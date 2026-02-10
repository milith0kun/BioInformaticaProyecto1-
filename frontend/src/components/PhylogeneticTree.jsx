/**
 * PhylogeneticTree Component
 * Visualize evolutionary relationships between compared genomes
 * Uses Canvas to render a dendrogram from UPGMA clustering
 */
import { useState, useEffect, useRef, useCallback } from 'react'
import api from '../services/api'

export default function PhylogeneticTree() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [showMatrix, setShowMatrix] = useState(false)
    const canvasRef = useRef(null)

    useEffect(() => { loadTree() }, [])

    const loadTree = async () => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getPhylogeneticTree()
            setData(result)
        } catch (e) {
            const detail = e.response?.data?.detail || 'Error generando árbol filogenético'
            setError(detail)
        } finally {
            setLoading(false)
        }
    }

    // Draw dendrogram on canvas
    const drawTree = useCallback(() => {
        if (!data?.tree || !canvasRef.current) return

        const canvas = canvasRef.current
        const ctx = canvas.getContext('2d')
        const dpr = window.devicePixelRatio || 1

        // Dynamic width calculation
        const containerWidth = canvas.parentElement.clientWidth || 800
        const displayW = Math.max(800, containerWidth - 40)

        // Handle single genome case
        if (data.single_genome) {
            const displayH = 200
            canvas.width = displayW * dpr
            canvas.height = displayH * dpr
            canvas.style.width = displayW + 'px'
            canvas.style.height = displayH + 'px'
            ctx.scale(dpr, dpr)
            ctx.fillStyle = '#ffffff'
            ctx.fillRect(0, 0, displayW, displayH)
            ctx.fillStyle = '#64748b'
            ctx.font = '700 14px system-ui, sans-serif'
            ctx.textAlign = 'center'
            ctx.fillText('Genoma Único — Comparación no disponible', displayW / 2, displayH / 2)
            return
        }

        const baseSpacing = (data.labels || []).length <= 5 ? 100 : (data.labels || []).length <= 10 ? 80 : 60
        const displayH = Math.max(400, (data.labels || []).length * baseSpacing)

        canvas.width = displayW * dpr
        canvas.height = displayH * dpr
        canvas.style.width = displayW + 'px'
        canvas.style.height = displayH + 'px'
        ctx.scale(dpr, dpr)

        // Clean Laboratory Background
        ctx.fillStyle = '#ffffff'
        ctx.fillRect(0, 0, displayW, displayH)

        const marginLeft = 80
        const marginRight = 350 // Ample space for long labels and Mb info
        const marginTop = 80
        const marginBottom = 80
        const treeW = displayW - marginLeft - marginRight
        const treeH = displayH - marginTop - marginBottom

        // Get leaf order from tree traversal to prevent crossing branches
        const orderedLeaves = []
        const getLeafOrder = (node) => {
            if (!node) return
            if (!node.children || node.children.length === 0) {
                orderedLeaves.push(node.name)
            } else {
                // Draw higher nodes first for cleaner look
                const sortedChildren = [...node.children].sort((a, b) => (b.height || 0) - (a.height || 0))
                sortedChildren.forEach(getLeafOrder)
            }
        }
        getLeafOrder(data.tree)

        const n = orderedLeaves.length
        const leafSpacing = treeH / (n - 1 || 1)

        // Map labels to y positions using tree order
        const leafY = {}
        orderedLeaves.forEach((label, i) => {
            leafY[label] = marginTop + i * leafSpacing
        })

        // Scale x position: root at marginLeft, leaves at marginLeft + treeW
        const getMaxHeight = (node) => {
            if (!node || !node.children || node.children.length === 0) return 0
            return Math.max(node.height || 0, ...node.children.map(getMaxHeight))
        }
        const maxH = getMaxHeight(data.tree) || 1

        // Root is at height maxH, leaves at height 0
        const xScale = (h) => marginLeft + (1 - h / maxH) * treeW

        const COLORS = ['#3b82f6', '#6366f1', '#8b5cf6', '#ec4899', '#f43f5e']

        // Recursive draw
        const drawNode = (node) => {
            if (!node) return { y: 0, x: 0 }
            if (!node.children || node.children.length === 0) {
                return { y: leafY[node.name] || 0, x: marginLeft + treeW }
            }

            const childPositions = node.children.map(drawNode)
            const nodeX = xScale(node.height || 0)
            const minY = Math.min(...childPositions.map(p => p.y))
            const maxY = Math.max(...childPositions.map(p => p.y))
            const nodeY = (minY + maxY) / 2

            // Draw vertical connector
            ctx.beginPath()
            ctx.strokeStyle = '#cbd5e1'
            ctx.lineWidth = 2.5
            ctx.lineCap = 'round'
            ctx.lineJoin = 'round'
            ctx.moveTo(nodeX, minY)
            ctx.lineTo(nodeX, maxY)
            ctx.stroke()

            // Draw horizontal branches
            childPositions.forEach((pos, i) => {
                ctx.beginPath()
                ctx.strokeStyle = '#cbd5e1'
                ctx.moveTo(nodeX, pos.y)
                ctx.lineTo(pos.x, pos.y)
                ctx.stroke()

                // Branch length label - Only show if significant and enough space
                const child = node.children[i]
                if (child && child.branch_length !== undefined && child.branch_length > 0.0001) {
                    const midX = (nodeX + pos.x) / 2
                    ctx.fillStyle = '#64748b'
                    ctx.font = 'italic 9px font-mono'
                    ctx.textAlign = 'center'
                    ctx.fillText((child.branch_length || 0).toFixed(4), midX, pos.y - 6)
                }
            })

            return { y: nodeY, x: nodeX }
        }

        drawNode(data.tree)

        // Draw leaf info and dots using tree order
        ctx.textAlign = 'left'
        ctx.textBaseline = 'middle'
        orderedLeaves.forEach((label, i) => {
            const y = leafY[label] || 0
            const x = marginLeft + treeW

            // Connection dot
            ctx.beginPath()
            ctx.fillStyle = COLORS[i % COLORS.length]
            ctx.arc(x, y, 6, 0, Math.PI * 2)
            ctx.fill()
            ctx.strokeStyle = '#ffffff'
            ctx.lineWidth = 2
            ctx.stroke()

            // Label - Use 900 for extra bold
            ctx.fillStyle = '#1e293b'
            ctx.font = '900 13px system-ui, sans-serif'
            ctx.fillText(label || 'Unknown', x + 20, y - 10)

            // Subtitle info
            const genome = data.genomes?.find(g => g.name === label)
            if (genome) {
                ctx.fillStyle = '#64748b'
                ctx.font = '500 11px system-ui, sans-serif'
                ctx.fillText(
                    `${genome.accession} • GC: ${genome.gc}% • ${(genome.length / 1e6 || 0).toFixed(2)} Mb`,
                    x + 20, y + 10
                )
            }
        })

        // Header Decoration
        ctx.fillStyle = '#3b82f6'
        ctx.fillRect(marginLeft, 30, 50, 5)
        ctx.fillStyle = '#0f172a'
        ctx.font = '900 18px system-ui, sans-serif'
        ctx.textAlign = 'left'
        ctx.fillText('ANÁLISIS FILOGENÉTICO UPGMA', marginLeft, 60)

        // Scale bar
        const scaleLen = 100
        const scaleVal = (maxH * (100 / treeW)).toFixed(3)
        const scaleY = displayH - 50
        ctx.beginPath()
        ctx.strokeStyle = '#64748b'
        ctx.lineWidth = 1.5
        ctx.moveTo(marginLeft, scaleY)
        ctx.lineTo(marginLeft + scaleLen, scaleY)
        ctx.moveTo(marginLeft, scaleY - 6)
        ctx.lineTo(marginLeft, scaleY + 6)
        ctx.moveTo(marginLeft + scaleLen, scaleY - 6)
        ctx.lineTo(marginLeft + scaleLen, scaleY + 6)
        ctx.stroke()
        ctx.fillStyle = '#64748b'
        ctx.font = 'bold 11px font-mono'
        ctx.textAlign = 'center'
        ctx.fillText(`${scaleVal} dist.`, marginLeft + scaleLen / 2, scaleY + 25)

    }, [data])

    useEffect(() => { drawTree() }, [drawTree])

    if (error) {
        const hasNoGenomes = error.includes("Encontrados: 0")
        const hasSingleGenome = error.includes("Encontrados: 1")

        return (
            <div className="space-y-4">
                <h2 className="text-xl font-bold text-slate-800">🌳 Árbol Filogenético</h2>
                <div className={`${hasNoGenomes ? 'bg-red-50 border-red-200' : 'bg-blue-50 border-blue-200'} border rounded-xl p-6 text-center`}>
                    <div className="text-6xl mb-4">
                        {hasNoGenomes ? '📂' : hasSingleGenome ? '🧬' : '⚠️'}
                    </div>
                    <p className={`${hasNoGenomes ? 'text-red-700' : 'text-blue-700'} font-medium text-lg mb-2`}>
                        {hasNoGenomes && 'No hay genomas descargados'}
                        {hasSingleGenome && 'Solo hay 1 genoma descargado'}
                        {!hasNoGenomes && !hasSingleGenome && error}
                    </p>
                    <p className={`${hasNoGenomes ? 'text-red-600' : 'text-blue-600'} text-sm mb-4`}>
                        {hasNoGenomes && 'Descargue al menos 1 genoma para ver su información, o 2+ genomas para construir un árbol filogenético.'}
                        {hasSingleGenome && 'Descargue al menos 1 genoma adicional para comparar y construir el árbol filogenético.'}
                        {!hasNoGenomes && !hasSingleGenome && 'Cargue más genomas desde la pestaña de archivos para construir el árbol.'}
                    </p>
                    <button
                        onClick={() => window.location.href = '#files'}
                        className={`px-6 py-2 ${hasNoGenomes ? 'bg-red-600 hover:bg-red-700' : 'bg-blue-600 hover:bg-blue-700'} text-white rounded-lg font-medium transition-colors`}
                    >
                        Ir a Descargar Genomas
                    </button>
                </div>
            </div>
        )
    }

    if (loading) {
        return (
            <div className="flex flex-col items-center justify-center py-48">
                <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
                <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse text-center">Calculando Distancias...</p>
                <p className="text-[9px] font-bold text-slate-400 mt-4 uppercase tracking-widest">Comparando arquitectura genómica</p>
            </div>
        )
    }

    if (!data) return null

    const genomesCount = data.labels?.length || 0
    const isSingleGenome = data.single_genome || genomesCount === 1

    return (
        <div className="space-y-10 animate-in fade-in duration-1000">
            <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
                <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
                <div className="space-y-4 relative z-10">
                    <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
                        Árbol <span className="text-blue-600">Filogenético</span>
                    </h2>
                    <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest">
                        {isSingleGenome ? (
                            <>
                                <span className="text-blue-600 font-black">1 GENOMA</span> CARGADO — {data.method}
                            </>
                        ) : (
                            <>
                                <span className="text-blue-600 font-black">{genomesCount} GENOMAS</span> COMPARADOS — {data.method}
                            </>
                        )}
                    </p>
                </div>
                <div className="flex gap-3 relative z-10">
                    {!isSingleGenome && (
                        <button onClick={() => setShowMatrix(!showMatrix)}
                            className={`px-6 py-2.5 rounded-2xl text-[9px] font-black uppercase tracking-widest transition-all shadow-xl active:scale-95 flex items-center gap-2 ${showMatrix
                                ? 'bg-blue-600 text-white'
                                : 'bg-white border-2 border-slate-100 text-slate-600 hover:border-blue-200'
                                }`}>
                            <span>📊</span> Matriz Distancias
                        </button>
                    )}
                    <button onClick={loadTree}
                        className="px-6 py-2.5 bg-slate-900 text-white rounded-2xl text-[9px] font-black uppercase tracking-widest transition-all hover:bg-blue-600 shadow-xl active:scale-95 flex items-center gap-2">
                        <span>🔄</span> Recalcular
                    </button>
                </div>
            </div>

            {/* Info banner for single genome */}
            {isSingleGenome && (
                <div className="bg-blue-50 border-2 border-blue-100 rounded-[2rem] p-8 flex items-center gap-6">
                    <div className="w-12 h-12 bg-white rounded-2xl flex items-center justify-center text-2xl shadow-sm">💡</div>
                    <p className="text-[10px] font-bold text-blue-700 uppercase tracking-widest leading-relaxed">
                        Se requieren al menos <span className="font-black">2 genomas</span> para construir un árbol comparativo.
                        Descarga genomas adicionales para habilitar el motor de clustering.
                    </p>
                </div>
            )}

            {/* Genome Summary Cards */}
            <div className={`grid gap-6 ${genomesCount === 1 ? 'grid-cols-1 max-w-md mx-auto' :
                genomesCount === 2 ? 'grid-cols-1 sm:grid-cols-2' :
                    'grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4'
                }`}>
                {data.genomes?.map((g, i) => (
                    <div
                        key={i}
                        className="bg-white rounded-[2rem] border-2 border-slate-100 p-6 flex items-center gap-5 hover:border-blue-200 hover:shadow-xl transition-all group"
                    >
                        <div
                            className="w-12 h-12 rounded-2xl flex items-center justify-center text-white font-black text-sm shadow-lg group-hover:rotate-12 transition-transform shrink-0"
                            style={{ backgroundColor: ['#2563eb', '#4f46e5', '#7c3aed', '#db2777', '#e11d48', '#0891b2', '#059669', '#ea580c'][i % 8] }}
                        >
                            {g.name.charAt(0)}
                        </div>
                        <div className="flex-1 min-w-0">
                            <p className="font-black text-slate-900 truncate text-xs uppercase tracking-tight group-hover:text-blue-600 transition-colors" title={g.name}>
                                {g.name}
                            </p>
                            <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest mt-1">
                                {g.accession} | GC: {g.gc}% | {(g.length / 1e6).toFixed(2)} Mb | {g.gene_count} genes
                            </p>
                        </div>
                    </div>
                ))}
            </div>

            {/* Dendrogram Canvas */}
            <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 overflow-x-auto shadow-sm relative">
                <div className="absolute top-8 left-10 flex items-center gap-3">
                    <div className="w-2 h-2 bg-blue-500 rounded-full animate-pulse"></div>
                    <span className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Dendrograma Visual</span>
                </div>
                <canvas ref={canvasRef} className="mx-auto" />
            </div>

            {/* Distance Matrix */}
            {showMatrix && data.distance_matrix && !isSingleGenome && (
                <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm animate-in fade-in slide-in-from-bottom-4 duration-500">
                    <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-8">Matriz de Distancia Molecular</h3>
                    <div className="overflow-x-auto custom-scrollbar">
                        <table className="w-full text-[10px] border-separate border-spacing-1">
                            <thead>
                                <tr>
                                    <th className="p-3"></th>
                                    {data.labels.map((l, i) => (
                                        <th key={i} className="p-3 text-slate-400 font-black text-left uppercase tracking-tighter truncate max-w-[120px]" title={l}>
                                            {l.split(' ').slice(0, 2).join(' ')}
                                        </th>
                                    ))}
                                </tr>
                            </thead>
                            <tbody>
                                {data.distance_matrix.map((row, i) => (
                                    <tr key={i}>
                                        <td className="p-3 font-black text-slate-700 uppercase tracking-tighter truncate max-w-[120px]" title={data.labels[i]}>
                                            {data.labels[i].split(' ').slice(0, 2).join(' ')}
                                        </td>
                                        {row.map((val, j) => (
                                            <td key={j} className="p-3 text-center font-mono rounded-xl transition-all hover:scale-110"
                                                style={{
                                                    backgroundColor: i === j ? '#f8fafc' :
                                                        `rgba(37, 99, 235, ${Math.max(0.05, 1 - (val || 0) * 2)})`,
                                                    color: (val || 0) < 0.3 ? '#ffffff' : '#475569',
                                                    fontWeight: (val || 0) < 0.3 ? '900' : '500'
                                                }}>
                                                {i === j ? '—' : (val || 0).toFixed(4)}
                                            </td>
                                        ))}
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            {/* Technical Methodology Info */}
            <div className="bg-slate-900 rounded-[2.5rem] p-10 text-white shadow-2xl relative overflow-hidden group shadow-blue-900/20">
                <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32 rounded-full"></div>
                <div className="flex flex-col md:flex-row gap-10 relative z-10">
                    <div className="w-16 h-16 bg-blue-600 rounded-3xl flex items-center justify-center text-3xl shadow-lg flex-shrink-0">📊</div>
                    <div className="space-y-4 flex-1">
                        <h4 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.3em]">Metodología UPGMA</h4>
                        <p className="text-sm font-medium text-slate-300 leading-relaxed max-w-4xl">
                            Dendrograma UPGMA: Clustering jerárquico no ponderado por pares (Unweighted Pair Group Method with Arithmetic mean). Las distancias se calculan combinando:
                        </p>
                        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 pt-4">
                            {[
                                { label: 'Contenido GC', pct: '30%', desc: 'Diferencia en % GC entre genomas' },
                                { label: 'Tamaño Genoma', pct: '10%', desc: 'Ratio de tamaños' },
                                { label: 'Contenido Génico', pct: '30%', desc: 'Índice Jaccard de genes compartidos' },
                                { label: 'Productos Génicos', pct: '30%', desc: 'Similitud de las funciones génicas' }
                            ].map(item => (
                                <div key={item.label} className="p-4 bg-white/5 rounded-2xl border border-white/10">
                                    <div className="flex justify-between items-baseline mb-2">
                                        <span className="text-[9px] font-black text-blue-400 uppercase">{item.label}</span>
                                        <span className="text-xs font-black">{item.pct}</span>
                                    </div>
                                    <p className="text-[9px] text-slate-500 font-bold uppercase">{item.desc}</p>
                                </div>
                            ))}
                        </div>
                        <p className="text-[10px] text-slate-500 font-bold italic pt-4">
                            Nota: Para análisis filogenético riguroso se recomienda usar alineamiento de 16S rRNA o genes housekeeping.
                        </p>
                    </div>
                </div>
                <div className="mt-8 pt-6 border-t border-white/5 flex items-center justify-end">
                    <span className="text-[9px] font-black text-blue-500 uppercase tracking-widest animate-pulse">LABORATORY SYSTEM ONLINE</span>
                </div>
            </div>
        </div>
    )
}