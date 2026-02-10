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

        // Handle single genome case
        if (data.single_genome) {
            const canvas = canvasRef.current
            const ctx = canvas.getContext('2d')
            const dpr = window.devicePixelRatio || 1

            const displayW = canvas.parentElement.clientWidth - 40
            const displayH = 200

            canvas.width = displayW * dpr
            canvas.height = displayH * dpr
            canvas.style.width = displayW + 'px'
            canvas.style.height = displayH + 'px'
            ctx.scale(dpr, dpr)

            // Clear
            ctx.fillStyle = '#f8fafc'
            ctx.fillRect(0, 0, displayW, displayH)

            // Draw single genome indicator
            ctx.fillStyle = '#334155'
            ctx.font = 'bold 14px system-ui, sans-serif'
            ctx.textAlign = 'center'
            ctx.fillText('Genoma Único - Sin comparación disponible', displayW / 2, displayH / 2 - 20)

            ctx.fillStyle = '#64748b'
            ctx.font = '12px system-ui, sans-serif'
            ctx.fillText(`${data.labels[0]}`, displayW / 2, displayH / 2 + 10)

            return
        }

        const canvas = canvasRef.current
        const ctx = canvas.getContext('2d')
        const dpr = window.devicePixelRatio || 1

        const displayW = canvas.parentElement.clientWidth - 40
        // Adjust height based on number of genomes (more space for more genomes)
        const baseSpacing = data.labels.length <= 5 ? 80 : data.labels.length <= 10 ? 60 : 45
        const displayH = Math.max(400, data.labels.length * baseSpacing)

        canvas.width = displayW * dpr
        canvas.height = displayH * dpr
        canvas.style.width = displayW + 'px'
        canvas.style.height = displayH + 'px'
        ctx.scale(dpr, dpr)

        // Clear
        ctx.fillStyle = '#f8fafc'
        ctx.fillRect(0, 0, displayW, displayH)

        const marginLeft = 40
        const marginRight = 220 // space for labels
        const marginTop = 30
        const marginBottom = 30
        const treeW = displayW - marginLeft - marginRight
        const treeH = displayH - marginTop - marginBottom

        const leaves = data.labels
        const n = leaves.length
        const leafSpacing = treeH / (n + 1)

        // Assign y positions to leaves
        const leafY = {}
        leaves.forEach((label, i) => {
            leafY[label] = marginTop + (i + 1) * leafSpacing
        })

        // Get max height for scaling
        const getMaxHeight = (node) => {
            if (!node.children || node.children.length === 0) return 0
            return Math.max(node.height || 0, ...node.children.map(getMaxHeight))
        }
        const maxH = getMaxHeight(data.tree) || 1

        // Scale x position based on height
        const xScale = (h) => marginLeft + (h / maxH) * treeW

        // Colors for branches
        const BRANCH_COLORS = ['#14b8a6', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6',
            '#ec4899', '#06b6d4', '#10b981', '#f97316', '#3b82f6']

        let colorIdx = 0

        // Recursive draw
        const drawNode = (node, depth = 0) => {
            if (!node.children || node.children.length === 0) {
                // Leaf node
                const y = leafY[node.name]
                if (y === undefined) return { y, x: marginLeft }
                return { y, x: marginLeft }
            }

            const childPositions = node.children.map(child => drawNode(child, depth + 1))

            // Node position
            const nodeX = xScale(node.height || 0)
            const minY = Math.min(...childPositions.map(p => p.y))
            const maxY = Math.max(...childPositions.map(p => p.y))
            const nodeY = (minY + maxY) / 2

            const branchColor = BRANCH_COLORS[colorIdx % BRANCH_COLORS.length]
            colorIdx++

            // Draw vertical connector
            ctx.beginPath()
            ctx.strokeStyle = branchColor
            ctx.lineWidth = 2
            ctx.moveTo(nodeX, minY)
            ctx.lineTo(nodeX, maxY)
            ctx.stroke()

            // Draw horizontal branches to children
            childPositions.forEach((pos, i) => {
                ctx.beginPath()
                ctx.strokeStyle = branchColor
                ctx.lineWidth = 2
                ctx.moveTo(nodeX, pos.y)
                ctx.lineTo(pos.x, pos.y)
                ctx.stroke()

                // Distance label on branch
                const child = node.children[i]
                if (child.branch_length !== undefined) {
                    const midX = (nodeX + pos.x) / 2
                    ctx.fillStyle = '#94a3b8'
                    ctx.font = '9px monospace'
                    ctx.textAlign = 'center'
                    ctx.fillText(child.branch_length.toFixed(3), midX, pos.y - 5)
                }
            })

            return { y: nodeY, x: nodeX }
        }

        drawNode(data.tree)

        // Draw leaf labels
        ctx.textAlign = 'left'
        ctx.textBaseline = 'middle'
        leaves.forEach((label, i) => {
            const y = leafY[label]

            // Draw dot
            ctx.beginPath()
            ctx.fillStyle = BRANCH_COLORS[i % BRANCH_COLORS.length]
            ctx.arc(marginLeft - 5, y, 4, 0, Math.PI * 2)
            ctx.fill()

            // Draw label
            ctx.fillStyle = '#1e293b'
            ctx.font = 'bold 12px system-ui, sans-serif'
            ctx.fillText(label, marginLeft + 5, y)

            // Draw genome info
            const genome = data.genomes?.find(g => g.name === label)
            if (genome) {
                ctx.fillStyle = '#64748b'
                ctx.font = '10px system-ui, sans-serif'
                ctx.fillText(
                    `${genome.accession} | GC: ${genome.gc}% | ${(genome.length / 1e6).toFixed(2)} Mb | ${genome.gene_count} genes`,
                    marginLeft + 5, y + 16
                )
            }
        })

        // Title
        ctx.fillStyle = '#334155'
        ctx.font = 'bold 14px system-ui, sans-serif'
        ctx.textAlign = 'left'
        ctx.fillText('Dendrograma UPGMA', marginLeft, 18)

        // Scale bar
        const scaleLen = treeW * 0.2
        const scaleVal = (maxH * 0.2).toFixed(3)
        const scaleY = displayH - 15
        ctx.beginPath()
        ctx.strokeStyle = '#64748b'
        ctx.lineWidth = 1.5
        ctx.moveTo(marginLeft, scaleY)
        ctx.lineTo(marginLeft + scaleLen, scaleY)
        ctx.stroke()
        ctx.fillStyle = '#64748b'
        ctx.font = '10px monospace'
        ctx.textAlign = 'center'
        ctx.fillText(scaleVal, marginLeft + scaleLen / 2, scaleY - 5)

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
            <div className="text-center py-16">
                <div className="w-12 h-12 border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-4"></div>
                <p className="text-slate-500">Calculando distancias filogenéticas...</p>
                <p className="text-xs text-slate-400 mt-1">Comparando genomas por GC, tamaño, y contenido génico</p>
            </div>
        )
    }

    if (!data) return null

    const genomesCount = data.labels?.length || 0
    const isSingleGenome = data.single_genome || genomesCount === 1

    return (
        <div className="space-y-6">
            <div className="flex items-start justify-between">
                <div>
                    <h2 className="text-xl font-bold text-slate-800">🌳 Árbol Filogenético</h2>
                    <p className="text-sm text-slate-500">
                        {isSingleGenome ? (
                            <>
                                <span className="font-medium text-blue-600">1 genoma</span> cargado — {data.method}
                            </>
                        ) : (
                            <>
                                <span className="font-medium text-teal-600">{genomesCount} genomas</span> comparados — {data.method}
                            </>
                        )}
                    </p>
                </div>
                <div className="flex gap-2">
                    {!isSingleGenome && (
                        <button onClick={() => setShowMatrix(!showMatrix)}
                            className={`px-3 py-1.5 rounded-lg text-xs font-medium transition-all ${showMatrix
                                ? 'bg-teal-600 text-white'
                                : 'bg-white border border-slate-200 text-slate-600 hover:border-teal-300'
                                }`}>
                            📊 Matriz Distancias
                        </button>
                    )}
                    <button onClick={loadTree}
                        className="px-3 py-1.5 bg-white border border-slate-200 rounded-lg text-xs font-medium text-slate-600 hover:border-teal-300">
                        🔄 Recalcular
                    </button>
                </div>
            </div>

            {/* Info banner for single genome */}
            {isSingleGenome && (
                <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                    <p className="text-blue-800 text-sm">
                        <span className="font-semibold">💡 Información:</span> Se necesitan al menos 2 genomas para construir un árbol filogenético comparativo.
                        Descargue genomas adicionales desde la pestaña de archivos para habilitar la comparación.
                    </p>
                </div>
            )}

            {/* Genome Summary Cards */}
            <div className={`grid gap-2 ${
                genomesCount === 1 ? 'grid-cols-1 max-w-md mx-auto' :
                genomesCount === 2 ? 'grid-cols-1 sm:grid-cols-2' :
                genomesCount <= 6 ? 'grid-cols-1 sm:grid-cols-2 lg:grid-cols-3' :
                genomesCount <= 12 ? 'grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4' :
                'grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 xl:grid-cols-6'
            }`}>
                {data.genomes?.map((g, i) => (
                    <div
                        key={i}
                        className={`bg-white rounded-lg border border-slate-200 flex items-center gap-3 hover:shadow-md transition-shadow ${
                            genomesCount > 12 ? 'p-2' : 'p-4'
                        }`}
                    >
                        <div
                            className={`rounded-lg flex items-center justify-center text-white font-bold flex-shrink-0 ${
                                genomesCount > 12 ? 'w-8 h-8 text-xs' : 'w-10 h-10 text-sm'
                            }`}
                            style={{ backgroundColor: ['#14b8a6', '#6366f1', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#06b6d4', '#10b981', '#f97316', '#3b82f6'][i % 10] }}
                        >
                            {g.name.charAt(0)}
                        </div>
                        <div className="flex-1 min-w-0">
                            <p
                                className={`font-semibold text-slate-800 truncate ${
                                    genomesCount > 12 ? 'text-xs' : 'text-sm'
                                }`}
                                title={g.name}
                            >
                                {g.name}
                            </p>
                            <p className={`text-slate-500 ${genomesCount > 12 ? 'text-[9px]' : 'text-[10px]'}`}>
                                {genomesCount > 12 ? (
                                    // Compact view for many genomes
                                    <>GC: {g.gc}% | {(g.length / 1e6).toFixed(1)} Mb</>
                                ) : (
                                    // Full view for fewer genomes
                                    <>{g.accession} | GC: {g.gc}% | {(g.length / 1e6).toFixed(2)} Mb | {g.gene_count} genes</>
                                )}
                            </p>
                        </div>
                    </div>
                ))}
            </div>

            {/* Summary stats for many genomes */}
            {genomesCount > 12 && (
                <div className="bg-slate-50 border border-slate-200 rounded-lg p-4">
                    <p className="text-sm text-slate-600">
                        <span className="font-semibold">📊 Resumen:</span> {genomesCount} genomas |
                        GC promedio: {(data.genomes.reduce((sum, g) => sum + g.gc, 0) / genomesCount).toFixed(2)}% |
                        Tamaño promedio: {(data.genomes.reduce((sum, g) => sum + g.length, 0) / genomesCount / 1e6).toFixed(2)} Mb |
                        Genes totales: {data.genomes.reduce((sum, g) => sum + g.gene_count, 0).toLocaleString()}
                    </p>
                </div>
            )}

            {/* Dendrogram Canvas */}
            <div className="bg-white rounded-xl border border-slate-200 p-5 overflow-x-auto">
                <canvas ref={canvasRef} />
            </div>

            {/* Distance Matrix */}
            {showMatrix && data.distance_matrix && !isSingleGenome && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <div className="flex items-center justify-between mb-4">
                        <h3 className="font-semibold text-slate-800">Matriz de Distancias Genéticas</h3>
                        {genomesCount > 8 && (
                            <span className="text-xs text-slate-500 bg-slate-100 px-2 py-1 rounded">
                                {genomesCount}×{genomesCount} matriz
                            </span>
                        )}
                    </div>
                    <div className="overflow-x-auto max-h-96 overflow-y-auto">
                        <table className="text-xs">
                            <thead>
                                <tr>
                                    <th className="px-2 py-1"></th>
                                    {data.labels.map((l, i) => (
                                        <th key={i} className="px-2 py-1 text-slate-600 font-medium text-left max-w-[120px] truncate" title={l}>
                                            {l.split(' ').slice(0, 2).join(' ')}
                                        </th>
                                    ))}
                                </tr>
                            </thead>
                            <tbody>
                                {data.distance_matrix.map((row, i) => (
                                    <tr key={i}>
                                        <td className="px-2 py-1 font-medium text-slate-600 max-w-[120px] truncate" title={data.labels[i]}>
                                            {data.labels[i].split(' ').slice(0, 2).join(' ')}
                                        </td>
                                        {row.map((val, j) => (
                                            <td key={j} className="px-2 py-1 text-center font-mono"
                                                style={{
                                                    backgroundColor: i === j ? '#f1f5f9' :
                                                        `rgba(20, 184, 166, ${Math.max(0.05, 1 - val * 3)})`,
                                                    color: val < 0.2 ? '#134e4a' : '#334155'
                                                }}>
                                                {i === j ? '—' : val.toFixed(4)}
                                            </td>
                                        ))}
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            {/* Info */}
            <div className="bg-teal-50 border border-teal-100 rounded-xl p-4">
                <div className="flex gap-3">
                    <div className="w-8 h-8 bg-teal-500 rounded-lg flex items-center justify-center flex-shrink-0">
                        <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
                            <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
                        </svg>
                    </div>
                    <div className="text-sm text-teal-800 space-y-1">
                        <p><strong>Dendrograma UPGMA</strong>: Clustering jerárquico no ponderado por pares (Unweighted Pair Group Method with Arithmetic mean). Las distancias se calculan combinando:</p>
                        <ul className="text-xs text-teal-700 list-disc list-inside space-y-0.5">
                            <li><strong>Contenido GC</strong> (30%): Diferencia en % GC entre genomas</li>
                            <li><strong>Tamaño del genoma</strong> (10%): Ratio de tamaños</li>
                            <li><strong>Contenido génico</strong> (30%): Índice Jaccard de genes compartidos</li>
                            <li><strong>Productos génicos</strong> (30%): Similitud de las funciones génicas</li>
                        </ul>
                        <p className="text-xs text-teal-700">Nota: Para análisis filogenético riguroso se recomienda usar alineamiento de 16S rRNA o genes housekeeping.</p>
                    </div>
                </div>
            </div>
        </div>
    )
}
