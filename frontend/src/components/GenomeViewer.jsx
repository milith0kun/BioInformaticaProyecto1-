/**
 * GenomeViewer Component — Enhanced
 * Interactive circular genome map with zoom, pan, gene selection
 * Sequence viewer with colored bases and gene annotation overlay
 * Position search with intuitive guidance
 */
import { useState, useEffect, useRef, useCallback } from 'react'
import api from '../services/api'

const BASE_COLORS = {
    A: { bg: '#22c55e', label: 'Adenina' },
    T: { bg: '#ef4444', label: 'Timina' },
    G: { bg: '#f59e0b', label: 'Guanina' },
    C: { bg: '#3b82f6', label: 'Citosina' },
}

export default function GenomeViewer() {
    const [mapData, setMapData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [selectedGene, setSelectedGene] = useState(null)
    const [hoveredGene, setHoveredGene] = useState(null)
    const [sequenceData, setSequenceData] = useState(null)
    const [seqStart, setSeqStart] = useState(0)
    const [seqLoading, setSeqLoading] = useState(false)
    const [searchPos, setSearchPos] = useState('')
    const [posResult, setPosResult] = useState(null)
    const [viewMode, setViewMode] = useState('map')
    const [zoom, setZoom] = useState(1)
    const [geneSearch, setGeneSearch] = useState('')
    const [seqWindowSize, setSeqWindowSize] = useState(600)
    const canvasRef = useRef(null)

    useEffect(() => { loadMapData() }, [])

    const loadMapData = async () => {
        setLoading(true)
        try {
            const data = await api.getGenomeMapData()
            setMapData(data)
        } catch (e) {
            console.error('Error loading genome map:', e)
        } finally {
            setLoading(false)
        }
    }

    const loadSequence = async (start = 0) => {
        setSeqLoading(true)
        try {
            const end = start + seqWindowSize
            const data = await api.getSequenceSegment(start, end)
            setSequenceData(data)
            setSeqStart(start)
        } catch (e) {
            console.error('Error loading sequence:', e)
        } finally {
            setSeqLoading(false)
        }
    }

    const searchPosition = async () => {
        const pos = parseInt(searchPos)
        if (isNaN(pos) || pos < 0) return
        try {
            const result = await api.getGeneAtPosition(pos)
            setPosResult(result)
        } catch (e) {
            console.error('Error searching position:', e)
        }
    }

    const selectGene = async (gene) => {
        setSelectedGene(gene)
        try {
            const detail = await api.getGeneDetail(gene.locus_tag)
            setSelectedGene({ ...gene, ...detail })
        } catch (e) { /* Keep basic info */ }
    }

    // Filtered gene list
    const filteredGenes = mapData?.genes?.filter(g => {
        if (!geneSearch) return true
        const q = geneSearch.toLowerCase()
        return (g.gene_name || '').toLowerCase().includes(q) ||
            (g.locus_tag || '').toLowerCase().includes(q) ||
            (g.product || '').toLowerCase().includes(q)
    }) || []

    // ===== DRAW CIRCULAR MAP =====
    const drawGenomeMap = useCallback(() => {
        if (!canvasRef.current || !mapData) return

        const canvas = canvasRef.current
        const ctx = canvas.getContext('2d')
        const dpr = window.devicePixelRatio || 1

        const displayWidth = canvas.parentElement.clientWidth
        const baseHeight = 500
        const displayHeight = Math.min(displayWidth, baseHeight) * zoom

        canvas.width = displayWidth * dpr
        canvas.height = displayHeight * dpr
        canvas.style.width = displayWidth + 'px'
        canvas.style.height = displayHeight + 'px'
        ctx.scale(dpr, dpr)

        const centerX = displayWidth / 2
        const centerY = displayHeight / 2
        const baseRadius = Math.min(centerX, centerY) - 70
        const radius = baseRadius

        // Clear
        ctx.clearRect(0, 0, displayWidth, displayHeight)
        ctx.fillStyle = '#0f172a'
        ctx.fillRect(0, 0, displayWidth, displayHeight)

        // Subtle grid circles
        for (let r = radius - 40; r <= radius + 40; r += 20) {
            ctx.beginPath()
            ctx.arc(centerX, centerY, r, 0, Math.PI * 2)
            ctx.strokeStyle = 'rgba(255,255,255,0.03)'
            ctx.lineWidth = 1
            ctx.stroke()
        }

        // GC Content ring (outer) — enhanced with smooth gradient
        if (mapData.gc_windows) {
            const gcRadius = radius + 30
            const avgGC = mapData.average_gc || 50

            mapData.gc_windows.forEach((w, i) => {
                const nextW = mapData.gc_windows[i + 1] || mapData.gc_windows[0]
                const startAngle = (w.angle - 90) * Math.PI / 180
                const endAngle = (nextW.angle - 90) * Math.PI / 180

                const deviation = w.gc - avgGC
                const barHeight = Math.min(Math.abs(deviation) * 2, 18)

                if (deviation >= 0) {
                    ctx.strokeStyle = `rgba(34, 197, 94, ${Math.min(Math.abs(deviation) / 8, 0.9)})`
                } else {
                    ctx.strokeStyle = `rgba(239, 68, 68, ${Math.min(Math.abs(deviation) / 8, 0.9)})`
                }

                ctx.beginPath()
                ctx.lineWidth = 3
                const r = gcRadius + (deviation >= 0 ? barHeight : -barHeight)
                ctx.arc(centerX, centerY, r, startAngle, endAngle)
                ctx.stroke()
            })

            // GC average line
            ctx.beginPath()
            ctx.arc(centerX, centerY, gcRadius, 0, Math.PI * 2)
            ctx.strokeStyle = 'rgba(255,255,255,0.15)'
            ctx.lineWidth = 0.5
            ctx.setLineDash([4, 4])
            ctx.stroke()
            ctx.setLineDash([])
        }

        // Main chromosome ring — glowing
        ctx.beginPath()
        ctx.arc(centerX, centerY, radius, 0, Math.PI * 2)
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.4)'
        ctx.lineWidth = 2
        ctx.stroke()

        // Inner ring for reverse strand
        ctx.beginPath()
        ctx.arc(centerX, centerY, radius - 22, 0, Math.PI * 2)
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.2)'
        ctx.lineWidth = 1
        ctx.stroke()

        // Draw genes — forward on outer, reverse on inner
        if (mapData.genes) {
            mapData.genes.forEach(gene => {
                const startAngle = (gene.angle_start - 90) * Math.PI / 180
                const endAngle = (gene.angle_end - 90) * Math.PI / 180

                if (Math.abs(endAngle - startAngle) < 0.001) return

                const isHovered = hoveredGene?.locus_tag === gene.locus_tag
                const isSelected = selectedGene?.locus_tag === gene.locus_tag

                const isForward = gene.strand === 1
                const color = isForward ? '#14b8a6' : '#a78bfa'
                const r = isForward ? radius + 2 : radius - 22

                ctx.beginPath()
                ctx.arc(centerX, centerY, r, startAngle, endAngle)

                if (isSelected) {
                    ctx.strokeStyle = '#fbbf24'
                    ctx.lineWidth = 7
                    ctx.shadowColor = '#fbbf24'
                    ctx.shadowBlur = 8
                } else if (isHovered) {
                    ctx.strokeStyle = '#34d399'
                    ctx.lineWidth = 5
                    ctx.shadowColor = '#34d399'
                    ctx.shadowBlur = 6
                } else {
                    ctx.strokeStyle = color
                    ctx.lineWidth = isForward ? 5 : 4
                    ctx.shadowColor = 'transparent'
                    ctx.shadowBlur = 0
                }

                ctx.globalAlpha = isHovered || isSelected ? 1 : 0.75
                ctx.stroke()
                ctx.globalAlpha = 1
                ctx.shadowBlur = 0
            })
        }

        // Scale markers  
        if (mapData.genome_length) {
            const tickInterval = mapData.genome_length > 3000000 ? 500000 : 100000
            for (let pos = 0; pos < mapData.genome_length; pos += tickInterval) {
                const angle = (pos / mapData.genome_length * 360 - 90) * Math.PI / 180
                const x1 = centerX + (radius + 50) * Math.cos(angle)
                const y1 = centerY + (radius + 50) * Math.sin(angle)
                const x2 = centerX + (radius + 56) * Math.cos(angle)
                const y2 = centerY + (radius + 56) * Math.sin(angle)

                ctx.beginPath()
                ctx.moveTo(x1, y1)
                ctx.lineTo(x2, y2)
                ctx.strokeStyle = 'rgba(148, 163, 184, 0.5)'
                ctx.lineWidth = 1
                ctx.stroke()

                const labelX = centerX + (radius + 64) * Math.cos(angle)
                const labelY = centerY + (radius + 64) * Math.sin(angle)
                ctx.fillStyle = 'rgba(148, 163, 184, 0.7)'
                ctx.font = `${zoom > 1.3 ? 10 : 8}px system-ui`
                ctx.textAlign = 'center'
                ctx.textBaseline = 'middle'
                ctx.fillText(`${(pos / 1000000).toFixed(1)}M`, labelX, labelY)
            }
        }

        // Center info — styled
        const organism = mapData.organism || 'Genoma'
        const orgParts = organism.split(' ')

        ctx.textAlign = 'center'
        ctx.textBaseline = 'middle'

        // Organism name
        ctx.fillStyle = '#f8fafc'
        ctx.font = 'bold 15px system-ui'
        if (orgParts.length >= 2) {
            ctx.font = 'italic bold 15px system-ui'
            ctx.fillText(orgParts[0], centerX, centerY - 30)
            ctx.font = 'italic 13px system-ui'
            ctx.fillText(orgParts.slice(1).join(' ').substring(0, 20), centerX, centerY - 12)
        } else {
            ctx.fillText(organism.substring(0, 20), centerX, centerY - 20)
        }

        // Stats
        ctx.font = '12px system-ui'
        ctx.fillStyle = '#94a3b8'
        ctx.fillText(`${(mapData.genome_length / 1e6).toFixed(2)} Mb`, centerX, centerY + 8)
        ctx.fillText(`${mapData.total_genes?.toLocaleString()} genes`, centerX, centerY + 24)

        // GC% with color
        const gc = mapData.average_gc || 0
        ctx.fillStyle = gc > 50 ? '#34d399' : '#fbbf24'
        ctx.font = 'bold 13px system-ui'
        ctx.fillText(`GC: ${gc}%`, centerX, centerY + 42)

    }, [mapData, hoveredGene, selectedGene, zoom])

    useEffect(() => { drawGenomeMap() }, [drawGenomeMap])

    useEffect(() => {
        const handleResize = () => drawGenomeMap()
        window.addEventListener('resize', handleResize)
        return () => window.removeEventListener('resize', handleResize)
    }, [drawGenomeMap])

    // Canvas mouse handler for hover
    const handleCanvasMouseMove = useCallback((e) => {
        if (!canvasRef.current || !mapData?.genes) return
        const canvas = canvasRef.current
        const rect = canvas.getBoundingClientRect()
        const mouseX = e.clientX - rect.left
        const mouseY = e.clientY - rect.top

        const centerX = rect.width / 2
        const centerY = rect.height / 2
        const baseRadius = Math.min(centerX, centerY) - 70

        const dx = mouseX - centerX
        const dy = mouseY - centerY
        const dist = Math.sqrt(dx * dx + dy * dy)
        let angle = Math.atan2(dy, dx) * 180 / Math.PI + 90
        if (angle < 0) angle += 360

        // Check if mouse is near a gene arc
        if (dist > baseRadius - 35 && dist < baseRadius + 20) {
            const found = mapData.genes.find(g => {
                let s = g.angle_start, en = g.angle_end
                if (s > en) return angle >= s || angle <= en
                return angle >= s && angle <= en
            })
            if (found) {
                if (hoveredGene?.locus_tag !== found.locus_tag) setHoveredGene(found)
                return
            }
        }
        if (hoveredGene) setHoveredGene(null)
    }, [mapData, hoveredGene])

    const handleCanvasClick = useCallback((e) => {
        if (hoveredGene) selectGene(hoveredGene)
    }, [hoveredGene])

    if (loading) {
        return (
            <div className="text-center py-20">
                <div className="w-16 h-16 mx-auto border-4 border-teal-200 border-t-teal-600 rounded-full animate-spin mb-6"></div>
                <h2 className="text-xl font-bold text-slate-800 mb-2">Cargando mapa genómico...</h2>
                <p className="text-slate-500">Procesando genes, contenido GC y anotaciones</p>
            </div>
        )
    }

    return (
        <div className="space-y-6">
            {/* View Mode Tabs */}
            <div className="flex gap-2 flex-wrap">
                {[
                    { id: 'map', label: 'Mapa Circular', icon: '🧬' },
                    { id: 'sequence', label: 'Visor de Secuencia', icon: '🔤' },
                    { id: 'position', label: 'Buscar Posición', icon: '📍' },
                ].map(tab => (
                    <button
                        key={tab.id}
                        onClick={() => {
                            setViewMode(tab.id)
                            if (tab.id === 'sequence' && !sequenceData) loadSequence(0)
                        }}
                        className={`px-4 py-2 rounded-lg text-sm font-medium transition-all ${viewMode === tab.id
                            ? 'bg-teal-600 text-white shadow-md'
                            : 'bg-white text-slate-600 border border-slate-200 hover:border-teal-300'
                            }`}
                    >
                        {tab.icon} {tab.label}
                    </button>
                ))}
            </div>

            {/* ==================== CIRCULAR MAP ==================== */}
            {viewMode === 'map' && mapData && (
                <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                    {/* Canvas Map */}
                    <div className="lg:col-span-2 bg-slate-900 rounded-xl border border-slate-700 p-4">
                        <div className="flex items-center justify-between mb-3">
                            <h3 className="font-semibold text-white text-sm">
                                Mapa Genómico — <em>{mapData.organism}</em>
                            </h3>
                            {/* Zoom Controls */}
                            <div className="flex items-center gap-2">
                                <button
                                    onClick={() => setZoom(z => Math.max(0.6, z - 0.15))}
                                    className="w-7 h-7 bg-slate-700 hover:bg-slate-600 text-white rounded flex items-center justify-center text-sm"
                                >−</button>
                                <span className="text-xs text-slate-400 w-12 text-center font-mono">{Math.round(zoom * 100)}%</span>
                                <button
                                    onClick={() => setZoom(z => Math.min(2, z + 0.15))}
                                    className="w-7 h-7 bg-slate-700 hover:bg-slate-600 text-white rounded flex items-center justify-center text-sm"
                                >+</button>
                                <button
                                    onClick={() => setZoom(1)}
                                    className="px-2 py-1 bg-slate-700 hover:bg-slate-600 text-slate-300 rounded text-xs"
                                >Reset</button>
                            </div>
                        </div>
                        <div className="overflow-auto" style={{ maxHeight: zoom > 1 ? '700px' : 'none' }}>
                            <canvas
                                ref={canvasRef}
                                className="w-full cursor-crosshair"
                                onMouseMove={handleCanvasMouseMove}
                                onClick={handleCanvasClick}
                            />
                        </div>

                        {/* Hover tooltip */}
                        {hoveredGene && (
                            <div className="mt-2 flex items-center gap-3 bg-slate-800 rounded-lg px-3 py-2 text-xs text-slate-200">
                                <span className={`w-2 h-2 rounded-full ${hoveredGene.strand === 1 ? 'bg-teal-400' : 'bg-violet-400'}`}></span>
                                <span className="font-mono font-bold">{hoveredGene.gene_name || hoveredGene.locus_tag}</span>
                                <span className="text-slate-400">{hoveredGene.product || 'Hipotético'}</span>
                                <span className="ml-auto text-slate-500">{hoveredGene.start?.toLocaleString()}–{hoveredGene.end?.toLocaleString()}</span>
                                <span className="text-slate-500">
                                    {hoveredGene.strand === 1 ? '5\'→3\' Forward' : '3\'→5\' Reverse'}
                                </span>
                            </div>
                        )}

                        {/* Legend */}
                        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 mt-4 text-xs text-slate-300">
                            <div className="flex items-center gap-2 bg-slate-800 rounded-lg px-3 py-2">
                                <span className="w-4 h-1 bg-teal-400 rounded"></span>
                                <span>Hebra + (5'→3' Forward)</span>
                            </div>
                            <div className="flex items-center gap-2 bg-slate-800 rounded-lg px-3 py-2">
                                <span className="w-4 h-1 bg-violet-400 rounded"></span>
                                <span>Hebra − (3'→5' Reverse)</span>
                            </div>
                            <div className="flex items-center gap-2 bg-slate-800 rounded-lg px-3 py-2">
                                <span className="w-4 h-1 bg-green-400 rounded"></span>
                                <span>GC alto (&gt; promedio)</span>
                            </div>
                            <div className="flex items-center gap-2 bg-slate-800 rounded-lg px-3 py-2">
                                <span className="w-4 h-1 bg-red-400 rounded"></span>
                                <span>GC bajo (&lt; promedio)</span>
                            </div>
                        </div>
                    </div>

                    {/* Gene List Sidebar */}
                    <div className="bg-white rounded-xl border border-slate-200 p-4 flex flex-col max-h-[650px]">
                        <h3 className="font-semibold text-slate-800 mb-2 text-sm">
                            Genes ({mapData.genes?.length?.toLocaleString() || 0})
                        </h3>
                        <input
                            type="text"
                            value={geneSearch}
                            onChange={(e) => setGeneSearch(e.target.value)}
                            placeholder="Filtrar genes..."
                            className="w-full px-3 py-2 mb-2 border border-slate-200 rounded-lg text-xs focus:outline-none focus:border-teal-400"
                        />
                        <div className="overflow-y-auto flex-1 space-y-1">
                            {filteredGenes.slice(0, 300).map((gene, i) => (
                                <button
                                    key={i}
                                    onClick={() => selectGene(gene)}
                                    onMouseEnter={() => setHoveredGene(gene)}
                                    onMouseLeave={() => setHoveredGene(null)}
                                    className={`w-full text-left px-3 py-2 rounded-lg text-xs transition-all ${selectedGene?.locus_tag === gene.locus_tag
                                        ? 'bg-teal-50 border border-teal-300'
                                        : 'hover:bg-slate-50'
                                        }`}
                                >
                                    <div className="flex items-center justify-between">
                                        <span className="font-mono font-medium text-slate-800">
                                            {gene.gene_name || gene.locus_tag}
                                        </span>
                                        <span className={`flex items-center gap-1 text-[10px] font-medium px-1.5 py-0.5 rounded ${gene.strand === 1
                                            ? 'bg-teal-50 text-teal-600'
                                            : 'bg-violet-50 text-violet-600'
                                            }`}>
                                            {gene.strand === 1 ? '→ 5\'→3\'' : '← 3\'→5\''}
                                        </span>
                                    </div>
                                    <p className="text-slate-500 truncate mt-0.5">{gene.product || 'Hipotético'}</p>
                                    <p className="text-slate-400 mt-0.5 font-mono">
                                        {gene.start?.toLocaleString()} — {gene.end?.toLocaleString()}
                                    </p>
                                </button>
                            ))}
                            {filteredGenes.length > 300 && (
                                <p className="text-center text-xs text-slate-400 py-2">
                                    +{(filteredGenes.length - 300).toLocaleString()} genes más. Use el filtro para refinar.
                                </p>
                            )}
                        </div>
                    </div>
                </div>
            )}

            {/* Selected Gene Detail Panel */}
            {selectedGene && (
                <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                    <div className="bg-gradient-to-r from-slate-800 to-teal-800 px-5 py-4 text-white">
                        <div className="flex items-start justify-between">
                            <div>
                                <h3 className="font-bold text-lg">
                                    {selectedGene.gene_name || selectedGene.locus_tag}
                                    <span className="text-sm font-normal text-teal-300 ml-2">
                                        ({selectedGene.locus_tag})
                                    </span>
                                </h3>
                                <p className="text-teal-200 text-sm">{selectedGene.product}</p>
                            </div>
                            <button onClick={() => setSelectedGene(null)} className="text-slate-300 hover:text-white">
                                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                                </svg>
                            </button>
                        </div>
                    </div>

                    <div className="p-5 space-y-4">
                        <div className="grid grid-cols-2 sm:grid-cols-5 gap-3 text-sm">
                            <div className="bg-slate-50 rounded-lg p-3">
                                <span className="text-xs text-slate-500">Posición</span>
                                <p className="font-mono font-medium">{selectedGene.start?.toLocaleString()} — {selectedGene.end?.toLocaleString()}</p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-3">
                                <span className="text-xs text-slate-500">Longitud</span>
                                <p className="font-mono font-bold text-lg">{selectedGene.length?.toLocaleString()} bp</p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-3">
                                <span className="text-xs text-slate-500">GC Content</span>
                                <p className="font-mono font-bold text-lg">{selectedGene.gc_content}%</p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-3">
                                <span className="text-xs text-slate-500">Hebra</span>
                                <p className="font-medium">
                                    <span className={`inline-flex items-center gap-1.5 px-2 py-1 rounded-lg text-sm ${selectedGene.strand === 1
                                        ? 'bg-teal-100 text-teal-700'
                                        : 'bg-violet-100 text-violet-700'
                                        }`}>
                                        {selectedGene.strand === 1 ? '→ 5\'→3\' Forward' : '← 3\'→5\' Reverse'}
                                    </span>
                                </p>
                            </div>
                            <div className="bg-slate-50 rounded-lg p-3">
                                <span className="text-xs text-slate-500">Protein ID</span>
                                <p className="font-mono text-sm">
                                    {selectedGene.protein_id ? (
                                        <a href={`https://www.ncbi.nlm.nih.gov/protein/${selectedGene.protein_id}`}
                                            target="_blank" rel="noopener noreferrer"
                                            className="text-teal-600 hover:underline">
                                            {selectedGene.protein_id}
                                        </a>
                                    ) : '—'}
                                </p>
                            </div>
                        </div>

                        {/* DNA Sequence */}
                        {selectedGene.sequence && (
                            <div>
                                <h4 className="text-xs uppercase tracking-wide text-slate-500 font-medium mb-2 flex items-center gap-2">
                                    DNA Sequence
                                    <span className="text-[10px] text-slate-400 font-normal normal-case">
                                        ({selectedGene.strand === 1 ? '5\'→3\' coding strand' : '3\'→5\' reverse complement shown as 5\'→3\''})
                                    </span>
                                </h4>
                                <div className="bg-slate-900 rounded-lg p-3 overflow-x-auto">
                                    <p className="font-mono text-xs break-all leading-relaxed">
                                        {selectedGene.sequence.substring(0, 240).split('').map((base, j) => (
                                            <span key={j} style={{ color: BASE_COLORS[base]?.bg || '#94a3b8' }}>{base}</span>
                                        ))}
                                        {selectedGene.sequence.length > 240 && <span className="text-slate-500">...({selectedGene.sequence.length - 240} más)</span>}
                                    </p>
                                </div>
                            </div>
                        )}

                        {/* Translation */}
                        {selectedGene.translation && (
                            <div>
                                <h4 className="text-xs uppercase tracking-wide text-slate-500 font-medium mb-2">
                                    Proteína ({selectedGene.translation.length} aa, ~{(selectedGene.translation.length * 110 / 1000).toFixed(1)} kDa)
                                </h4>
                                <div className="bg-slate-800 rounded-lg p-3 overflow-x-auto">
                                    <p className="font-mono text-xs break-all leading-relaxed">
                                        {selectedGene.translation.substring(0, 150).split('').map((aa, j) => {
                                            const hydrophobic = 'AILMFWVP'.includes(aa)
                                            const charged = 'DEKRH'.includes(aa)
                                            const polar = 'NQSTY'.includes(aa)
                                            const color = charged ? '#ff6b6b' : hydrophobic ? '#4ade80' : polar ? '#60a5fa' : '#fbbf24'
                                            return <span key={j} style={{ color }}>{aa}</span>
                                        })}
                                        {selectedGene.translation.length > 150 && <span className="text-slate-500">...</span>}
                                    </p>
                                </div>
                            </div>
                        )}

                        {/* NCBI Links */}
                        <div className="flex flex-wrap gap-2">
                            {selectedGene.protein_id && (
                                <a href={`https://www.ncbi.nlm.nih.gov/protein/${selectedGene.protein_id}`}
                                    target="_blank" rel="noopener noreferrer"
                                    className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-teal-50 border border-teal-200 rounded-lg text-xs text-teal-700 hover:bg-teal-100">
                                    🔗 NCBI Protein
                                </a>
                            )}
                            {selectedGene.gene_name && (
                                <a href={`https://www.ncbi.nlm.nih.gov/gene/?term=${selectedGene.gene_name}`}
                                    target="_blank" rel="noopener noreferrer"
                                    className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-blue-50 border border-blue-200 rounded-lg text-xs text-blue-700 hover:bg-blue-100">
                                    🔗 NCBI Gene
                                </a>
                            )}
                            {selectedGene.db_xref?.map((ref, i) => (
                                <span key={i} className="px-2 py-1 bg-slate-100 text-slate-600 rounded text-xs font-mono">{ref}</span>
                            ))}
                        </div>
                    </div>
                </div>
            )}

            {/* ==================== SEQUENCE VIEWER ==================== */}
            {viewMode === 'sequence' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <div className="flex flex-col sm:flex-row items-start sm:items-center gap-4 mb-4">
                        <h3 className="font-semibold text-slate-800">🔤 Visor de Secuencia Genómica</h3>
                        <div className="flex items-center gap-2 flex-wrap">
                            <div className="flex items-center gap-1">
                                <label className="text-xs text-slate-500">Posición:</label>
                                <input
                                    type="number"
                                    value={seqStart}
                                    onChange={(e) => setSeqStart(parseInt(e.target.value) || 0)}
                                    className="w-28 px-2.5 py-1.5 border border-slate-200 rounded-lg text-sm font-mono focus:outline-none focus:border-teal-400"
                                    min={0}
                                />
                            </div>
                            <div className="flex items-center gap-1">
                                <label className="text-xs text-slate-500">Bases:</label>
                                <select
                                    value={seqWindowSize}
                                    onChange={(e) => setSeqWindowSize(parseInt(e.target.value))}
                                    className="px-2 py-1.5 border border-slate-200 rounded-lg text-sm"
                                >
                                    <option value={300}>300 bp</option>
                                    <option value={600}>600 bp</option>
                                    <option value={1200}>1,200 bp</option>
                                    <option value={3000}>3,000 bp</option>
                                </select>
                            </div>
                            <button
                                onClick={() => loadSequence(seqStart)}
                                className="px-3 py-1.5 bg-teal-600 text-white rounded-lg text-sm font-medium hover:bg-teal-700"
                            >
                                Ir
                            </button>
                        </div>
                        {sequenceData && (
                            <div className="flex gap-1.5 ml-auto">
                                <button
                                    onClick={() => loadSequence(Math.max(0, seqStart - seqWindowSize))}
                                    disabled={seqStart === 0}
                                    className="px-3 py-1.5 bg-slate-100 text-slate-600 rounded-lg text-sm hover:bg-slate-200 disabled:opacity-40"
                                >
                                    ← Anterior
                                </button>
                                <button
                                    onClick={() => loadSequence(seqStart + seqWindowSize)}
                                    className="px-3 py-1.5 bg-slate-100 text-slate-600 rounded-lg text-sm hover:bg-slate-200"
                                >
                                    Siguiente →
                                </button>
                            </div>
                        )}
                    </div>

                    {/* Help box */}
                    {!sequenceData && !seqLoading && (
                        <div className="bg-teal-50 border border-teal-100 rounded-xl p-4 mb-4">
                            <h4 className="font-medium text-teal-800 text-sm mb-2">💡 ¿Cómo usar el visor?</h4>
                            <ul className="text-xs text-teal-700 space-y-1">
                                <li>• Ingrese una <strong>posición genómica</strong> (ej: 0, 500000, 1000000) y haga clic en "Ir"</li>
                                <li>• Seleccione cuántas bases desea visualizar (300 a 3,000)</li>
                                <li>• Las bases están coloreadas: <span className="font-bold" style={{ color: '#22c55e' }}>A</span> <span className="font-bold" style={{ color: '#ef4444' }}>T</span> <span className="font-bold" style={{ color: '#f59e0b' }}>G</span> <span className="font-bold" style={{ color: '#3b82f6' }}>C</span></li>
                                <li>• Los genes en la región se muestran como etiquetas con su hebra</li>
                                <li>• Use "← Anterior" y "Siguiente →" para navegar</li>
                            </ul>
                        </div>
                    )}

                    {seqLoading ? (
                        <div className="text-center py-10">
                            <div className="w-8 h-8 border-2 border-teal-200 border-t-teal-600 rounded-full animate-spin mx-auto mb-3"></div>
                            <p className="text-sm text-slate-500">Cargando secuencia...</p>
                        </div>
                    ) : sequenceData ? (
                        <>
                            {/* Genes in view */}
                            {sequenceData.genes && sequenceData.genes.length > 0 && (
                                <div className="mb-4">
                                    <p className="text-xs text-slate-500 mb-2">Genes en esta región:</p>
                                    <div className="flex flex-wrap gap-2">
                                        {sequenceData.genes.map((g, i) => (
                                            <span key={i} className={`px-2.5 py-1.5 rounded-lg text-xs font-medium border ${g.strand === 1
                                                ? 'bg-teal-50 text-teal-700 border-teal-200'
                                                : 'bg-violet-50 text-violet-700 border-violet-200'
                                                }`}>
                                                {g.strand === 1 ? '→' : '←'} {g.gene_name || g.locus_tag}
                                                <span className="text-slate-400 ml-1 font-mono text-[10px]">
                                                    {g.start.toLocaleString()}–{g.end.toLocaleString()}
                                                </span>
                                                <span className="text-slate-400 ml-1">
                                                    ({g.strand === 1 ? '5\'→3\'' : '3\'→5\''})
                                                </span>
                                            </span>
                                        ))}
                                    </div>
                                </div>
                            )}

                            {/* Sequence display */}
                            <div className="bg-slate-900 rounded-xl p-4 overflow-x-auto font-mono text-xs">
                                {sequenceData.formatted?.map((line, i) => (
                                    <div key={i} className="flex gap-3 leading-relaxed">
                                        <span className="text-slate-500 w-16 text-right select-none flex-shrink-0">
                                            {line.position.toLocaleString()}
                                        </span>
                                        <span className="tracking-wider">
                                            {line.sequence.split('').map((base, j) => (
                                                <span key={j} style={{ color: BASE_COLORS[base]?.bg || '#94a3b8' }}>{base}</span>
                                            ))}
                                        </span>
                                        <span className="text-slate-600 text-[10px] select-none flex-shrink-0">
                                            GC:{line.gc_count}
                                        </span>
                                    </div>
                                ))}
                            </div>

                            <div className="mt-3 flex items-center justify-between text-xs text-slate-500 flex-wrap gap-2">
                                <span>Región: {sequenceData.start.toLocaleString()} — {sequenceData.end.toLocaleString()} ({(sequenceData.end - sequenceData.start).toLocaleString()} bp)</span>
                                <span>GC: {sequenceData.gc_content}% | Genoma total: {(sequenceData.genome_length / 1e6).toFixed(2)} Mb</span>
                            </div>

                            {/* Base color legend */}
                            <div className="mt-3 flex items-center gap-4 text-xs text-slate-500">
                                {Object.entries(BASE_COLORS).map(([base, { bg, label }]) => (
                                    <span key={base} className="flex items-center gap-1">
                                        <span className="w-2.5 h-2.5 rounded" style={{ backgroundColor: bg }}></span>
                                        <span style={{ color: bg }} className="font-mono font-bold">{base}</span>
                                        <span>{label}</span>
                                    </span>
                                ))}
                            </div>
                        </>
                    ) : null}
                </div>
            )}

            {/* ==================== POSITION SEARCH ==================== */}
            {viewMode === 'position' && (
                <div className="bg-white rounded-xl border border-slate-200 p-5">
                    <h3 className="font-semibold text-slate-800 mb-2">📍 Buscar Gen por Posición Genómica</h3>
                    <p className="text-xs text-slate-500 mb-4">
                        Ingrese una posición en pares de bases (bp) para encontrar qué gen se ubica en esa coordenada del genoma.
                        {mapData && (
                            <span className="ml-1">Rango válido: <strong className="font-mono">0 — {mapData.genome_length?.toLocaleString()}</strong></span>
                        )}
                    </p>

                    <div className="flex items-center gap-3 mb-6">
                        <div className="relative flex-1">
                            <input
                                type="number"
                                value={searchPos}
                                onChange={(e) => setSearchPos(e.target.value)}
                                onKeyDown={(e) => e.key === 'Enter' && searchPosition()}
                                placeholder="Ej: 500000, 1500000, 3000000..."
                                className="w-full px-4 py-3 pl-9 border-2 border-slate-200 rounded-xl text-sm font-mono focus:border-teal-400 focus:outline-none"
                                min={0}
                                max={mapData?.genome_length}
                            />
                            <svg className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
                            </svg>
                        </div>
                        <button
                            onClick={searchPosition}
                            className="px-6 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-medium hover:from-teal-700 hover:to-emerald-700 transition-all"
                        >
                            Buscar
                        </button>
                    </div>

                    {/* Quick position buttons */}
                    {mapData && (
                        <div className="flex flex-wrap gap-2 mb-4">
                            <span className="text-xs text-slate-500 py-1">Posiciones rápidas:</span>
                            {[0, 100000, 500000, 1000000, 2000000, 3000000, 4000000].filter(p => p < (mapData.genome_length || 0)).map(p => (
                                <button
                                    key={p}
                                    onClick={() => { setSearchPos(String(p)); }}
                                    className="px-2 py-1 bg-slate-100 hover:bg-teal-50 border border-slate-200 hover:border-teal-300 rounded text-xs font-mono transition-all"
                                >
                                    {p === 0 ? 'Inicio' : `${(p / 1e6).toFixed(1)}M`}
                                </button>
                            ))}
                        </div>
                    )}

                    {posResult && (
                        <div className={`rounded-xl p-5 border ${posResult.found ? 'bg-emerald-50 border-emerald-200' : 'bg-amber-50 border-amber-200'}`}>
                            {posResult.found ? (
                                <>
                                    <div className="flex items-center gap-2 mb-3">
                                        <span className="text-emerald-600 text-lg">✓</span>
                                        <h4 className="font-bold text-slate-800">Gen encontrado en posición {posResult.position?.toLocaleString()}</h4>
                                    </div>
                                    <div className="grid grid-cols-2 sm:grid-cols-3 gap-4 text-sm">
                                        <div>
                                            <span className="text-xs text-slate-500">Locus Tag</span>
                                            <p className="font-mono font-medium">{posResult.locus_tag}</p>
                                        </div>
                                        <div>
                                            <span className="text-xs text-slate-500">Gen</span>
                                            <p className="font-medium">{posResult.gene_name || '—'}</p>
                                        </div>
                                        <div>
                                            <span className="text-xs text-slate-500">Producto</span>
                                            <p className="text-sm">{posResult.product || 'Hipotético'}</p>
                                        </div>
                                        <div>
                                            <span className="text-xs text-slate-500">Hebra</span>
                                            <p className={`font-medium ${posResult.strand === 1 ? 'text-teal-600' : 'text-violet-600'}`}>
                                                {posResult.strand === 1 ? '→ 5\'→3\' Forward' : '← 3\'→5\' Reverse'}
                                            </p>
                                        </div>
                                        <div>
                                            <span className="text-xs text-slate-500">Rango</span>
                                            <p className="font-mono text-xs">{posResult.start?.toLocaleString()} — {posResult.end?.toLocaleString()}</p>
                                        </div>
                                        <div>
                                            <span className="text-xs text-slate-500">GC%</span>
                                            <p className="font-mono">{posResult.gc_content}%</p>
                                        </div>
                                    </div>
                                    {posResult.protein_id && (
                                        <div className="mt-3">
                                            <a href={`https://www.ncbi.nlm.nih.gov/protein/${posResult.protein_id}`}
                                                target="_blank" rel="noopener noreferrer"
                                                className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-teal-100 border border-teal-200 rounded-lg text-xs text-teal-700 hover:bg-teal-200">
                                                🔗 NCBI Protein: {posResult.protein_id}
                                            </a>
                                        </div>
                                    )}
                                </>
                            ) : (
                                <>
                                    <h4 className="font-bold text-slate-800 mb-2">
                                        {posResult.message || 'Región intergénica — no hay gen en esta posición'}
                                    </h4>
                                    {posResult.nearest_genes && (
                                        <div className="mt-3">
                                            <p className="text-xs text-slate-600 mb-2">Genes más cercanos:</p>
                                            <div className="space-y-2">
                                                {posResult.nearest_genes.map((g, i) => (
                                                    <div key={i} className="flex items-center gap-3 text-sm bg-white p-2 rounded-lg border border-slate-100">
                                                        <span className={`px-2 py-0.5 rounded text-xs font-medium ${g.direction === 'upstream'
                                                            ? 'bg-blue-100 text-blue-700'
                                                            : 'bg-orange-100 text-orange-700'
                                                            }`}>
                                                            {g.direction === 'upstream' ? '↑' : '↓'} {g.distance?.toLocaleString()} bp
                                                        </span>
                                                        <span className="font-mono text-xs font-medium">{g.gene_name || g.locus_tag}</span>
                                                        <span className="text-slate-400 text-xs">
                                                            {g.start?.toLocaleString()}–{g.end?.toLocaleString()}
                                                        </span>
                                                    </div>
                                                ))}
                                            </div>
                                        </div>
                                    )}
                                </>
                            )}
                        </div>
                    )}
                </div>
            )}
        </div>
    )
}
