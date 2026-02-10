/**
 * GenomeViewer Component — Clean Laboratory Edition
 * Comprehensive genome visualization with Circular Map, Sequence Inspection, and Position Lookup.
 */
import { useState, useEffect, useRef, useMemo } from 'react'
import api from '../services/api'
import toast from 'react-hot-toast'

export default function GenomeViewer() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [activeTab, setActiveTab] = useState('circular')
  const [hoveredItem, setHoveredItem] = useState(null)
  const [searchTerm, setSearchTerm] = useState('')
  const [zoom, setZoom] = useState(100)
  const [visibleLayers, setVisibility] = useState({
    forward: true,
    reverse: true,
    gc: true,
    skew: true
  })

  // Sequence Viewer States
  const [sequenceData, setSequenceData] = useState(null)
  const [seqLoading, setSeqLoading] = useState(false)
  const [seqConfig, setSeqConfig] = useState({ start: 0, length: 600 })

  // Position Search State
  const [searchPos, setSearchPos] = useState('')
  const [searchResult, setSearchResult] = useState(null)
  const [searchLoading, setSearchLoading] = useState(false)

  const containerRef = useRef(null)

  useEffect(() => {
    loadMapData()
  }, [])

  useEffect(() => {
    if (activeTab === 'sequence') {
      loadSequence(seqConfig.start, seqConfig.start + seqConfig.length)
    }
  }, [activeTab, seqConfig])

  const loadMapData = async () => {
    setLoading(true)
    try {
      const result = await api.getGenomeMapData()
      setData(result)
    } catch (e) {
      setError('Error cargando mapa genómico')
    } finally {
      setLoading(false)
    }
  }

  const loadSequence = async (start, end) => {
    setSeqLoading(true)
    try {
      const result = await api.getSequenceSegment(start, end)
      setSequenceData(result)
    } catch (e) {
      console.error("Error cargando secuencia:", e)
    } finally {
      setSeqLoading(false)
    }
  }

  const handleJump = (e) => {
    e.preventDefault()
    const pos = parseInt(e.target.elements.pos.value)
    const len = parseInt(e.target.elements.len.value)
    if (!isNaN(pos) && !isNaN(len)) {
      setSeqConfig({ start: Math.max(0, pos), length: Math.min(len, 5000) })
    }
  }

  const navigateSeq = (direction) => {
    const step = Math.floor(seqConfig.length / 2)
    const newStart = direction === 'next'
      ? seqConfig.start + step
      : Math.max(0, seqConfig.start - step)
    setSeqConfig(prev => ({ ...prev, start: newStart }))
  }

  const performPositionSearch = async (pos) => {
    const targetPos = parseInt(pos)
    if (isNaN(targetPos)) return

    setSearchLoading(true)
    try {
      const result = await api.getGeneAtPosition(targetPos)
      setSearchResult(result)
      if (!result) toast.error('Posición fuera de rango')
    } catch (e) {
      toast.error('Error en la búsqueda')
    } finally {
      setSearchLoading(false)
    }
  }

  // Mathematics for Circular Map
  const polarToCartesian = (centerX, centerY, radius, angleInDegrees) => {
    const angleInRadians = (angleInDegrees - 90) * Math.PI / 180.0
    return {
      x: centerX + (radius * Math.cos(angleInRadians)),
      y: centerY + (radius * Math.sin(angleInRadians))
    }
  }

  const describeArc = (x, y, radius, startAngle, endAngle) => {
    const start = polarToCartesian(x, y, radius, endAngle)
    const end = polarToCartesian(x, y, radius, startAngle)
    const largeArcFlag = endAngle - startAngle <= 180 ? "0" : "1"
    return ["M", start.x, start.y, "A", radius, radius, 0, largeArcFlag, 0, end.x, end.y].join(" ")
  }

  const filteredGenes = useMemo(() => {
    if (!data?.genes) return []
    if (!searchTerm) return data.genes
    const s = searchTerm.toLowerCase()
    return data.genes.filter(g => (g.gene_name?.toLowerCase().includes(s)) || (g.locus_tag?.toLowerCase().includes(s)))
  }, [data, searchTerm])

  const layers = useMemo(() => {
    if (!data) return null
    const cx = 50, cy = 50
    return {
      forward: (data.genes || []).filter(g => g.strand === 1).slice(0, 800).map((g, i) => {
        const isHovered = hoveredItem?.locus_tag === g.locus_tag
        return (
          <path
            key={`f-${i}`}
            d={describeArc(cx, cy, 42, g.angle_start, g.angle_end)}
            fill="none"
            stroke={isHovered ? "#2563eb" : "#3b82f6"}
            strokeWidth={isHovered ? "6" : "3"}
            className="cursor-pointer transition-all opacity-80 hover:opacity-100"
            onMouseEnter={() => setHoveredItem({ ...g, type: 'Sentido (+)' })}
            onMouseLeave={() => setHoveredItem(null)}
          />
        )
      }),
      reverse: (data.genes || []).filter(g => g.strand === -1).slice(0, 800).map((g, i) => {
        const isHovered = hoveredItem?.locus_tag === g.locus_tag
        return (
          <path
            key={`r-${i}`}
            d={describeArc(cx, cy, 38, g.angle_start, g.angle_end)}
            fill="none"
            stroke={isHovered ? "#4f46e5" : "#6366f1"}
            strokeWidth={isHovered ? "6" : "3"}
            className="cursor-pointer transition-all opacity-80 hover:opacity-100"
            onMouseEnter={() => setHoveredItem({ ...g, type: 'Antisentido (-)' })}
            onMouseLeave={() => setHoveredItem(null)}
          />
        )
      }),
      gc: (data.gc_windows || []).map((w, i, arr) => {
        const nextAngle = arr[i + 1]?.angle || 360
        const isHigh = w.gc > (data.average_gc || 0)
        return (
          <path
            key={`gc-${i}`}
            d={describeArc(cx, cy, 47, w.angle, nextAngle)}
            fill="none"
            stroke={isHigh ? '#94a3b8' : '#e2e8f0'}
            strokeWidth="4"
            className="cursor-pointer hover:stroke-blue-300 transition-all"
            onMouseEnter={() => setHoveredItem({ 
              type: isHigh ? 'GC alto (> promedio)' : 'GC bajo (< promedio)',
              gc_content: (w.gc || 0).toFixed(2),
              start: w.position || 0,
              product: `Ventana de contenido GC detectada en la posición ${(w.position || 0).toLocaleString()} pb.`
            })}
            onMouseLeave={() => setHoveredItem(null)}
          />
        )
      }),
      skew: (data.gc_skew || []).map((s, i, arr) => {
        const nextAngle = arr[i + 1]?.angle || 360
        const color = s.skew > 0 ? '#f43f5e' : '#fb7185'
        return (
          <path
            key={`sk-${i}`}
            d={describeArc(cx, cy, 33, s.angle, nextAngle)}
            fill="none"
            stroke={color}
            strokeWidth="3"
            strokeOpacity={Math.min(Math.abs(s.skew) * 5 + 0.2, 1)}
            className="cursor-pointer hover:stroke-rose-400 transition-all"
            onMouseEnter={() => setHoveredItem({ 
              type: 'GC Skew Picos',
              gc_content: (s.skew || 0).toFixed(4),
              start: s.position || 0,
              product: `Variación de asimetría GC (Skew) detectada en la posición ${(s.position || 0).toLocaleString()} pb.`
            })}
            onMouseLeave={() => setHoveredItem(null)}
          />
        )
      })
    }
  }, [data, hoveredItem])

  if (loading) return (
    <div className="flex flex-col items-center justify-center py-48">
      <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
      <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse">Cartografía Molecular...</p>
    </div>
  )

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-700">
      {/* Header */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            Explorador <span className="text-blue-600">Genómico</span>
          </h2>
          <div className="flex items-center gap-4">
            <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-blue-100">{data.accession}</span>
            <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest italic">{data.organism}</p>
          </div>
        </div>

        <div className="flex bg-slate-50 p-1.5 rounded-2xl border border-slate-100 relative z-10">
          {[
            { id: 'circular', label: 'Mapa Circular', icon: '🧬' },
            { id: 'sequence', label: 'Visor de Secuencia', icon: '🔤' },
            { id: 'search', label: 'Buscar Posición', icon: '📍' }
          ].map(tab => (
            <button key={tab.id} onClick={() => setActiveTab(tab.id)} className={`flex items-center gap-3 px-6 py-2.5 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all ${activeTab === tab.id ? 'bg-white text-blue-600 shadow-sm border border-slate-100' : 'text-slate-400 hover:text-slate-600'}`}>
              <span>{tab.icon}</span><span>{tab.label}</span>
            </button>
          ))}
        </div>
      </div>

      {activeTab === 'circular' && (
        <div className="grid grid-cols-1 xl:grid-cols-4 gap-8">
          <div className="xl:col-span-3 space-y-8">
            <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden flex flex-col min-h-[750px]">
              {/* Controls */}
              <div className="absolute top-8 left-8 flex items-center gap-6 z-20 bg-white/80 backdrop-blur-md px-6 py-3 rounded-full border border-slate-100 shadow-xl">
                <div className="flex items-center gap-3 pr-6 border-r border-slate-100">
                  <div className="w-2 h-2 bg-blue-500 rounded-full animate-pulse shadow-[0_0_8px_rgba(59,130,246,0.5)]"></div>
                  <span className="text-[9px] font-black text-slate-900 uppercase tracking-widest">Atlas Interactivo</span>
                </div>
                <div className="flex items-center gap-2">
                  <button onClick={() => setZoom(z => Math.max(50, z - 10))} className="w-8 h-8 rounded-xl bg-slate-50 text-slate-600 font-black shadow-sm">−</button>
                  <span className="text-[10px] font-black text-blue-600 w-12 text-center">{zoom}%</span>
                  <button onClick={() => setZoom(z => Math.min(300, z + 10))} className="w-8 h-8 rounded-xl bg-slate-50 text-slate-600 font-black shadow-sm">+</button>
                </div>
              </div>

              <div className="flex-1 flex items-center justify-center relative">
                <div style={{ transform: `scale(${zoom / 100})`, transition: 'transform 0.3s ease-out' }} className="w-full max-w-[600px] aspect-square">
                  <svg viewBox="0 0 100 100" className="w-full h-full drop-shadow-2xl">
                    {visibleLayers.gc && layers.gc}
                    {visibleLayers.forward && layers.forward}
                    {visibleLayers.reverse && layers.reverse}
                    {visibleLayers.skew && layers.skew}
                    <circle cx="50" cy="50" r="28" fill="white" className="shadow-inner" />
                    <text x="50" y="48" textAnchor="middle" className="text-[4px] font-black fill-slate-900 uppercase italic">MG1655</text>
                    <text x="50" y="54" textAnchor="middle" className="text-[3px] font-bold fill-blue-500 uppercase tracking-[0.3em]">ATLAS</text>
                  </svg>
                </div>
              </div>

              {/* Legend Toggles */}
              <div className="bg-slate-50/50 rounded-3xl p-6 flex flex-wrap gap-8 items-center justify-center border border-slate-100">
                {[
                  { id: 'forward', label: 'Hebra +', color: 'bg-blue-500' },
                  { id: 'reverse', label: 'Hebra −', color: 'bg-indigo-500' },
                  { id: 'gc', label: 'Contenido GC', color: 'bg-slate-400' },
                  { id: 'skew', label: 'Skew Picos', color: 'bg-rose-500' }
                ].map(layer => (
                  <button key={layer.id} onClick={() => setVisibility(v => ({ ...v, [layer.id]: !v[layer.id] }))} className={`flex items-center gap-3 transition-all ${visibleLayers[layer.id] ? 'opacity-100 scale-105' : 'opacity-30 grayscale scale-95'}`}>
                    <div className={`w-3 h-3 rounded-full ${layer.color} shadow-sm`}></div>
                    <span className="text-[10px] font-black text-slate-700 uppercase tracking-widest">{layer.label}</span>
                  </button>
                ))}
              </div>
            </div>

            {/* Bottom Info Panel */}
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm min-h-[120px] transition-all duration-500">
              {hoveredItem ? (
                <div className="animate-in fade-in slide-in-from-bottom-2 duration-300">
                  <div className="flex items-center justify-between mb-4">
                    <div className="flex items-center gap-4">
                      <span className="px-3 py-1 bg-blue-600 text-white text-[9px] font-black uppercase tracking-widest rounded-lg">{hoveredItem.type}</span>
                      <h4 className="font-mono text-xl font-black text-slate-900 uppercase tracking-tighter">{hoveredItem.gene_name || hoveredItem.locus_tag || 'Región Genómica'}</h4>
                    </div>
                    <div className="text-right">
                      <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">
                        {hoveredItem.type === 'GC Skew Picos' ? 'Valor Skew' : 'Contenido GC'}
                      </p>
                      <p className={`text-lg font-black ${hoveredItem.type === 'GC Skew Picos' ? 'text-rose-600' : 'text-blue-600'}`}>
                        {hoveredItem.gc_content}{hoveredItem.type.includes('GC') && !hoveredItem.type.includes('Skew') ? '%' : ''}
                      </p>
                    </div>
                  </div>
                  <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
                    <div className="md:col-span-2">
                      <p className="text-xs font-medium text-slate-500 leading-relaxed uppercase tracking-tight line-clamp-2">{hoveredItem.product || 'Sin descripción funcional disponible.'}</p>
                    </div>
                    <div className="flex justify-end items-center gap-6">
                      <div>
                        <p className="text-[8px] font-black text-slate-400 uppercase tracking-widest">Coordenadas</p>
                        <p className="text-xs font-bold text-slate-700 font-mono">{(hoveredItem.start || 0).toLocaleString()} — {(hoveredItem.end || 0).toLocaleString()}</p>
                      </div>
                      {hoveredItem.length && (
                        <div className="border-l border-slate-100 pl-6">
                          <p className="text-[8px] font-black text-slate-400 uppercase tracking-widest">Longitud</p>
                          <p className="text-xs font-bold text-slate-700">{(hoveredItem.length || 0).toLocaleString()} bp</p>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              ) : (
                <div className="h-full flex items-center justify-center opacity-30 italic">
                  <p className="text-xs font-bold text-slate-400 uppercase tracking-[0.2em]">Pase el cursor sobre un elemento del mapa para ver los detalles técnicos</p>
                </div>
              )}
            </div>
          </div>

          {/* Sidebar */}
          <div className="space-y-6">
            <div className="bg-slate-900 rounded-[2.5rem] p-8 text-white shadow-2xl relative overflow-hidden">
              <div className="relative z-10 space-y-6">
                <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em]">Índice de Genes</h4>
                <div className="relative">
                  <input type="text" value={searchTerm} onChange={(e) => setSearchTerm(e.target.value)} placeholder="Filtrar..." className="w-full bg-white/5 border border-white/10 rounded-2xl px-10 py-3 text-xs font-bold text-white placeholder-slate-600 focus:outline-none focus:border-blue-500/50" />
                  <svg className="w-4 h-4 text-slate-600 absolute left-4 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeWidth={2.5} /></svg>
                </div>
              </div>
            </div>
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 shadow-sm overflow-hidden flex flex-col max-h-[600px]">
              <div className="flex-1 overflow-y-auto p-4 space-y-3 custom-scrollbar">
                {filteredGenes.slice(0, 100).map((gene, i) => (
                  <div key={i} className={`p-5 rounded-2xl border transition-all group cursor-pointer ${hoveredItem?.locus_tag === gene.locus_tag ? 'bg-blue-600 border-blue-500 text-white shadow-lg' : 'bg-slate-50 border-slate-100 hover:bg-blue-50'}`}
                    onMouseEnter={() => setHoveredItem({ ...gene, type: gene.strand === 1 ? 'Sentido (+)' : 'Antisentido (-)' })}
                    onMouseLeave={() => setHoveredItem(null)}
                    onClick={() => {
                      setSeqConfig({ start: gene.start || 0, length: 600 })
                      setActiveTab('sequence')
                    }}
                  >
                    <div className="flex justify-between items-start mb-2">
                      <p className={`font-mono text-xs font-black uppercase tracking-tighter ${hoveredItem?.locus_tag === gene.locus_tag ? 'text-white' : 'text-slate-900 group-hover:text-blue-600'}`}>{gene.gene_name || gene.locus_tag}</p>
                      <span className={`text-[9px] font-black px-2 py-0.5 rounded ${gene.strand === 1 ? 'bg-blue-100 text-blue-600' : 'bg-indigo-100 text-indigo-600'}`}>{gene.strand === 1 ? '→' : '←'}</span>
                    </div>
                    <p className={`text-[10px] font-bold uppercase truncate mb-3 ${hoveredItem?.locus_tag === gene.locus_tag ? 'text-blue-100' : 'text-slate-400'}`}>{gene.product || 'Hipotético'}</p>
                    <div className={`flex justify-between items-center border-t pt-3 ${hoveredItem?.locus_tag === gene.locus_tag ? 'border-white/10' : 'border-slate-200/50'}`}>
                      <span className={`text-[9px] font-black ${hoveredItem?.locus_tag === gene.locus_tag ? 'text-white/60' : 'text-slate-500'}`}>{(gene.start || 0).toLocaleString()} — {(gene.end || 0).toLocaleString()}</span>
                      <span className={`text-[9px] font-bold uppercase italic ${hoveredItem?.locus_tag === gene.locus_tag ? 'text-white/40' : 'text-slate-300'}`}>pb</span>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {activeTab === 'sequence' && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm animate-in slide-in-from-right-4 duration-500">
          {/* Seq Header */}
          <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-8 mb-10 pb-8 border-b border-slate-50">
            <div className="flex items-center gap-6">
              <div className="w-14 h-14 bg-blue-600 rounded-2xl flex items-center justify-center text-white text-2xl shadow-xl shadow-blue-200">🔤</div>
              <div>
                <h3 className="text-xl font-black text-slate-900 uppercase tracking-tighter italic">Visor de Secuencia Genómica</h3>
                <p className="text-[10px] font-bold text-slate-400 uppercase tracking-widest mt-1">Interfaz de Inspección Nucleotídica</p>
              </div>
            </div>

            <form onSubmit={handleJump} className="flex flex-wrap items-center gap-4 bg-slate-50 p-4 rounded-3xl border border-slate-100">
              <div className="flex items-center gap-3 px-4 py-2 bg-white rounded-2xl border border-slate-100 shadow-sm">
                <span className="text-[9px] font-black text-slate-400 uppercase">Posición:</span>
                <input name="pos" type="number" defaultValue={seqConfig.start} className="w-24 text-xs font-black focus:outline-none" />
              </div>
              <div className="flex items-center gap-3 px-4 py-2 bg-white rounded-2xl border border-slate-100 shadow-sm">
                <span className="text-[9px] font-black text-slate-400 uppercase">Bases:</span>
                <input name="len" type="number" defaultValue={seqConfig.length} className="w-16 text-xs font-black focus:outline-none" />
              </div>
              <button type="submit" className="px-8 py-3 bg-blue-600 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-700 transition-all shadow-lg shadow-blue-200">Ir</button>
            </form>
          </div>

          <div className="grid grid-cols-1 xl:grid-cols-4 gap-10">
            {/* Main Sequence View */}
            <div className="xl:col-span-3 space-y-8">
              <div className="flex justify-between items-center mb-4">
                <button onClick={() => navigateSeq('prev')} className="px-6 py-2.5 bg-slate-100 hover:bg-slate-200 rounded-xl text-[9px] font-black uppercase tracking-widest text-slate-600 transition-all">← Anterior</button>
                <button onClick={() => navigateSeq('next')} className="px-6 py-2.5 bg-slate-100 hover:bg-slate-200 rounded-xl text-[9px] font-black uppercase tracking-widest text-slate-600 transition-all">Siguiente →</button>
              </div>

              <div className="bg-slate-900 rounded-[2.5rem] p-10 shadow-2xl relative overflow-hidden min-h-[400px]">
                {seqLoading ? (
                  <div className="absolute inset-0 flex flex-col items-center justify-center bg-slate-900/50 backdrop-blur-sm z-10">
                    <div className="w-10 h-10 border-4 border-slate-700 border-t-blue-500 rounded-full animate-spin mb-4"></div>
                    <p className="text-[9px] font-black text-slate-500 uppercase tracking-widest">Leyendo Hebra...</p>
                  </div>
                ) : null}

                <div className="font-mono text-sm leading-[2.2] relative z-0">
                  {sequenceData?.formatted?.map((line, i) => (
                    <div key={i} className="flex gap-8 group hover:bg-white/5 p-1 rounded transition-colors">
                      <span className="text-blue-500/40 text-[10px] font-black w-20 text-right select-none">{line.position.toLocaleString()}</span>
                      <div className="flex-1 break-all tracking-[0.3em] font-black text-slate-300">
                        {line.sequence.split('').map((b, j) => (
                          <span key={j} className={b === 'A' ? 'text-emerald-400' : b === 'T' ? 'text-rose-400' : b === 'G' ? 'text-blue-400' : b === 'C' ? 'text-amber-400' : ''}>{b}</span>
                        ))}
                      </div>
                    </div>
                  ))}
                </div>
              </div>

              {/* Stats Footer */}
              <div className="flex flex-col md:flex-row justify-between items-center p-8 bg-slate-50 rounded-[2rem] border border-slate-100 gap-6">
                <div className="flex gap-10">
                  <div className="space-y-1">
                    <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Región Visualizada</p>
                    <p className="text-sm font-black text-slate-900 uppercase italic">{(seqConfig.start || 0).toLocaleString()} — {((seqConfig.start || 0) + (sequenceData?.length || 0)).toLocaleString()} <span className="text-[9px] text-slate-400 not-italic">({sequenceData?.length || 0} bp)</span></p>
                  </div>
                  <div className="space-y-1 border-l border-slate-200 pl-10">
                    <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Contenido GC Local</p>
                    <p className="text-sm font-black text-blue-600">{sequenceData?.gc_content || 0}%</p>
                  </div>
                </div>
                <div className="text-right">
                  <p className="text-[9px] font-black text-slate-400 uppercase tracking-widest">Genoma Total</p>
                  <p className="text-sm font-black text-slate-500">{((data?.genome_length || 0) / 1e6).toFixed(2)} Mb</p>
                </div>
              </div>

              {/* Atomic Legend */}
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                {[
                  { label: 'Adenina', b: 'A', color: 'bg-emerald-400' },
                  { label: 'Timina', b: 'T', color: 'bg-rose-400' },
                  { label: 'Guanina', b: 'G', color: 'bg-blue-400' },
                  { label: 'Citosina', b: 'C', color: 'bg-amber-400' }
                ].map(l => (
                  <div key={l.b} className="flex items-center gap-4 p-4 bg-white border-2 border-slate-100 rounded-2xl shadow-sm">
                    <div className={`w-8 h-8 rounded-xl ${l.color} flex items-center justify-center text-white font-black shadow-lg`}>{l.b}</div>
                    <span className="text-[10px] font-black text-slate-700 uppercase tracking-widest">{l.label}</span>
                  </div>
                ))}
              </div>
            </div>

            {/* Region Genes Sidebar */}
            <div className="space-y-6">
              <div className="bg-slate-900 rounded-[2.5rem] p-8 text-white relative overflow-hidden">
                <h4 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.2em] relative z-10">Genes en esta región</h4>
                <p className="text-[9px] text-slate-500 font-bold uppercase mt-2 relative z-10">({sequenceData?.genes?.length || 0} detectados)</p>
              </div>
              <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 shadow-sm overflow-hidden flex flex-col max-h-[600px] overflow-y-auto custom-scrollbar">
                <div className="p-4 space-y-3">
                  {sequenceData?.genes?.map((g, i) => (
                    <div key={i} className="p-5 bg-slate-50 rounded-2xl border border-slate-100 group transition-all">
                      <div className="flex justify-between items-start mb-2">
                        <p className="font-mono text-xs font-black text-slate-900 uppercase">→ {g.gene_name || g.locus_tag}</p>
                        <span className="text-[8px] font-black px-2 py-0.5 bg-blue-100 text-blue-600 rounded">({g.strand === 1 ? "5'→3'" : "3'→5'"})</span>
                      </div>
                      <p className="text-[9px] font-black text-slate-500 uppercase tracking-tighter">{(g.start || 0).toLocaleString()} — {(g.end || 0).toLocaleString()}</p>
                    </div>
                  ))}
                  {!sequenceData?.genes?.length && <p className="text-[10px] text-slate-300 font-black text-center py-10 uppercase tracking-widest">Región Intergénica</p>}
                </div>
              </div>
            </div>
          </div>
        </div>
      )}

      {activeTab === 'search' && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-12 shadow-sm min-h-[600px] flex flex-col space-y-12">
          <div className="max-w-2xl mx-auto w-full space-y-10">
            <div className="text-center space-y-4">
              <div className="w-20 h-20 bg-blue-50 rounded-[2rem] flex items-center justify-center mx-auto shadow-inner">
                <span className="text-4xl">📍</span>
              </div>
              <h3 className="text-3xl font-black text-slate-900 uppercase tracking-tighter">📍 Buscar Gen por Posición Genómica</h3>
              <p className="text-slate-500 font-medium leading-relaxed text-sm">Ingrese una posición en pares de bases (bp) para encontrar qué gen se ubica en esa coordenada del genoma.</p>
              <div className="p-4 bg-slate-50 rounded-2xl border border-slate-100">
                <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Rango válido: <span className="text-blue-600 font-black">0 — {data.genome_length?.toLocaleString()}</span></p>
              </div>
            </div>

            <div className="space-y-6">
              <div className="flex gap-4 p-2 bg-slate-50 rounded-3xl border border-slate-100 shadow-inner">
                <input
                  type="number"
                  value={searchPos}
                  onChange={(e) => setSearchPos(e.target.value)}
                  onKeyDown={(e) => e.key === 'Enter' && performPositionSearch(searchPos)}
                  placeholder="Ej: 500000, 1500000, 3000000..."
                  className="flex-1 px-8 py-4 bg-transparent border-none text-sm font-black text-slate-900 placeholder-slate-400 focus:ring-0"
                />
                <button
                  onClick={() => performPositionSearch(searchPos)}
                  disabled={searchLoading}
                  className="px-10 py-4 bg-blue-600 text-white rounded-2xl text-xs font-black uppercase tracking-widest hover:bg-blue-700 transition-all shadow-xl active:scale-95 disabled:opacity-50"
                >
                  {searchLoading ? 'Buscando...' : 'Buscar'}
                </button>
              </div>

              <div className="flex flex-wrap items-center justify-center gap-3">
                <span className="text-[10px] font-black text-slate-400 uppercase tracking-widest mr-2">Posiciones rápidas:</span>
                {[
                  { label: 'Inicio', val: 1 },
                  { label: '0.1M', val: 100000 },
                  { label: '0.5M', val: 500000 },
                  { label: '1.0M', val: 1000000 },
                  { label: '2.0M', val: 2000000 },
                  { label: '3.0M', val: 3000000 },
                  { label: '4.0M', val: 4000000 }
                ].map(pos => (
                  <button
                    key={pos.label}
                    onClick={() => {
                      setSearchPos(pos.val)
                      performPositionSearch(pos.val)
                    }}
                    className="px-4 py-2 bg-white border border-slate-200 rounded-xl text-[10px] font-black text-slate-600 hover:border-blue-500 hover:text-blue-600 transition-all shadow-sm"
                  >
                    {pos.label}
                  </button>
                ))}
              </div>
            </div>

            {searchResult && (
              <div className="animate-in fade-in zoom-in-95 duration-500 bg-white rounded-[2.5rem] border-2 border-blue-100 p-10 shadow-xl shadow-blue-500/5">
                {searchResult.type === 'intergenic' ? (
                  <div className="text-center space-y-4">
                    <span className="px-4 py-1.5 bg-slate-100 text-slate-500 text-[10px] font-black uppercase rounded-full">Región Intergénica</span>
                    <p className="text-lg font-black text-slate-900">Coordenada {(parseInt(searchPos) || 0).toLocaleString()} bp</p>
                    <p className="text-sm text-slate-500 font-medium leading-relaxed">{searchResult.message}</p>
                  </div>
                ) : (
                  <div className="space-y-8">
                    <div className="flex items-start justify-between">
                      <div className="space-y-2">
                        <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[10px] font-black uppercase rounded-full border border-blue-100">Gen Detectado</span>
                        <h4 className="text-4xl font-black text-slate-900 tracking-tighter uppercase italic">{searchResult.gene_name || searchResult.locus_tag}</h4>
                      </div>
                      <div className="text-right">
                        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Hebra</p>
                        <p className={`text-xl font-black ${searchResult.strand === 1 ? 'text-emerald-600' : 'text-indigo-600'}`}>
                          {searchResult.strand === 1 ? 'Forward (+)' : 'Reverse (-)'}
                        </p>
                      </div>
                    </div>

                    <div className="grid grid-cols-2 gap-8 py-8 border-y border-slate-50">
                      <div>
                        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-2">Ubicación</p>
                        <p className="text-sm font-bold text-slate-700 font-mono">{(searchResult.start ?? 0).toLocaleString()} — {(searchResult.end ?? 0).toLocaleString()} bp</p>
                      </div>
                      <div>
                        <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-2">Locus Tag</p>
                        <p className="text-sm font-bold text-slate-700 font-mono">{searchResult.locus_tag}</p>
                      </div>
                    </div>

                    <div>
                      <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-3">Descripción Funcional</p>
                      <p className="text-sm font-medium text-slate-600 leading-relaxed uppercase tracking-tight">{searchResult.product}</p>
                    </div>

                    <button
                      onClick={() => {
                        setSeqConfig({ start: searchResult.start ?? 0, length: 600 })
                        setActiveTab('sequence')
                      }}
                      className="w-full py-4 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-600 transition-all shadow-lg"
                    >
                      Analizar Secuencia del Gen
                    </button>
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  )
}
