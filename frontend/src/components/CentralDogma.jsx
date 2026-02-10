/**
 * CentralDogma Component — Clean Laboratory Edition
 * Advanced molecular flow visualization: DNA -> mRNA -> Protein
 */
import { useState, useEffect } from 'react'
import api from '../services/api'

const AA_COLORS = {
  // Hidrofóbicos
  'A': 'bg-emerald-500', 'L': 'bg-emerald-500', 'I': 'bg-emerald-500', 'V': 'bg-emerald-500', 
  'M': 'bg-emerald-500', 'F': 'bg-emerald-500', 'W': 'bg-emerald-500', 'P': 'bg-emerald-500',
  // Polares / Cargados +
  'S': 'bg-blue-500', 'T': 'bg-blue-500', 'N': 'bg-blue-500', 'Q': 'bg-blue-500', 
  'Y': 'bg-blue-500', 'C': 'bg-blue-500', 'G': 'bg-blue-500', 'K': 'bg-indigo-600', 
  'R': 'bg-indigo-600', 'H': 'bg-indigo-600',
  // Cargados -
  'D': 'bg-rose-500', 'E': 'bg-rose-500',
  // Stop
  '*': 'bg-slate-900'
}

const BASE_COLORS = {
  'A': 'text-emerald-500',
  'T': 'text-rose-500',
  'G': 'text-blue-500',
  'C': 'text-amber-500',
  'U': 'text-indigo-500'
}

export default function CentralDogma() {
  const [locusTag, setLocusTag] = useState('')
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [geneList, setGeneList] = useState([])
  const [showFullTable, setShowFullTable] = useState(false)

  useEffect(() => {
    loadInitialGenes()
  }, [])

  const loadInitialGenes = async () => {
    try {
      const result = await api.getGeneResults(1, 24)
      setGeneList(result?.genes || [])
    } catch (e) {
      console.warn('Genes no disponibles aún:', e.message)
      setGeneList([])
    }
  }

  const handleSearch = async (tag) => {
    const target = tag || locusTag
    if (!target.trim()) return
    setLoading(true)
    setError(null)
    setData(null)
    setShowFullTable(false)
    try {
      const result = await api.getCentralDogma(target.trim())
      setData(result)
    } catch (err) {
      setError('Locus Tag no identificado en el genoma activo.')
    } finally { setLoading(false) }
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      
      {/* 1. Buscador y Selector Superior (Aprovechando Ancho) */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="flex flex-col lg:flex-row lg:items-center gap-10 relative z-10">
          <div className="space-y-2 flex-shrink-0">
            <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
              Dogma <span className="text-blue-600">Central</span>
            </h2>
            <p className="text-[10px] font-bold text-slate-500 uppercase tracking-[0.2em]">DNA → mRNA → Proteína</p>
          </div>

          <div className="flex-1 space-y-4">
            <form onSubmit={(e) => { e.preventDefault(); handleSearch(); }} className="flex gap-3">
              <input
                type="text"
                value={locusTag}
                onChange={(e) => setLocusTag(e.target.value)}
                placeholder="Ingrese locus_tag (ej: b0001)..."
                className="flex-1 max-w-md px-6 py-3 bg-slate-50 border-2 border-slate-100 rounded-2xl text-slate-900 text-xs font-bold focus:border-blue-500 focus:outline-none transition-all placeholder-slate-400"
              />
              <button type="submit" className="px-8 py-3 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest hover:bg-blue-600 transition-all shadow-lg active:scale-95">Buscar</button>
            </form>
            
            <div className="flex items-center gap-4">
              <span className="text-[9px] font-black text-slate-400 uppercase tracking-widest flex-shrink-0">Selección Rápida:</span>
              <div className="flex gap-2 overflow-x-auto pb-2 scrollbar-hide">
                {geneList.slice(0, 12).map(g => (
                  <button
                    key={g.locus_tag}
                    onClick={() => handleSearch(g.locus_tag)}
                    className={`px-4 py-1.5 rounded-xl border-2 text-[10px] font-black uppercase whitespace-nowrap transition-all ${
                      data?.locus_tag === g.locus_tag ? 'bg-blue-600 border-blue-600 text-white' : 'bg-white border-slate-100 text-slate-500 hover:border-blue-200 hover:text-blue-600'
                    }`}
                  >
                    {g.gene_name || g.locus_tag}
                  </button>
                ))}
              </div>
            </div>
          </div>
        </div>
      </div>

      {!data && !loading && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-20 text-center flex flex-col items-center justify-center min-h-[500px]">
          <div className="w-24 h-24 bg-blue-50 rounded-full flex items-center justify-center text-5xl mb-8 animate-bounce">🧬</div>
          <h3 className="text-2xl font-black text-slate-900 uppercase tracking-tighter mb-4">Simulador de Transmisión Molecular</h3>
          <p className="text-slate-500 max-w-md mx-auto text-sm leading-relaxed font-medium">Seleccione un gen para visualizar el proceso completo de flujo de información genética con alta precisión científica.</p>
        </div>
      )}

      {loading && (
        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-20 flex flex-col items-center justify-center min-h-[500px]">
          <div className="w-16 h-16 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
          <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse">Secuenciando Dogma...</p>
        </div>
      )}

      {error && (
        <div className="bg-rose-50 rounded-[2.5rem] border-2 border-rose-100 p-10 text-center text-rose-600 font-black uppercase text-xs tracking-widest">
          {error}
        </div>
      )}

      {data && !loading && (
        <div className="space-y-10 animate-in zoom-in-95 duration-700">
          
          {/* 2. GRID Superior: Metadata y Diagrama de Flujo */}
          <div className="grid grid-cols-1 xl:grid-cols-3 gap-8">
            {/* Metadata del Gen */}
            <div className="xl:col-span-1 bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm relative overflow-hidden flex flex-col justify-between">
              <div className="absolute top-0 right-0 p-8 opacity-5 pointer-events-none text-8xl font-black text-slate-200">⬡</div>
              <div className="space-y-6 relative z-10">
                <div>
                  <h2 className="text-4xl font-black text-slate-900 tracking-tighter uppercase italic leading-none">{data.gene_name || data.locus_tag}</h2>
                  <p className="text-[10px] font-black text-blue-600 uppercase tracking-widest mt-2">({data.locus_tag})</p>
                  <p className="text-[11px] font-bold text-slate-500 uppercase tracking-tight mt-4 leading-relaxed">{data.product}</p>
                </div>
                
                <div className="grid grid-cols-2 gap-y-6 gap-x-10 pt-6 border-t border-slate-50">
                  {[
                    { label: 'Posición', val: `${data.start}-${data.end}`, icon: '📍' },
                    { label: 'Hebra', val: data.strand_label, icon: '🔗' },
                    { label: 'DNA', val: `${data.length_bp} bp`, icon: '🧬' },
                    { label: 'Proteína', val: `${data.length_aa} aa`, icon: '🧪' },
                    { label: 'GC%', val: `${data.gc_content}%`, icon: '📊' },
                    { label: 'Masa', val: `${data.molecular_weight_kda} kDa`, icon: '⚖️' }
                  ].map(stat => (
                    <div key={stat.label}>
                      <p className="text-[8px] font-black text-slate-400 uppercase tracking-widest mb-1">{stat.label}</p>
                      <p className="text-xs font-black text-slate-800 uppercase flex items-center gap-2">
                        {stat.val}
                      </p>
                    </div>
                  ))}
                </div>
              </div>

              <div className="mt-10 pt-8 border-t border-slate-50 flex flex-wrap gap-4 relative z-10">
                {data.ncbi_protein_url && <a href={data.ncbi_protein_url} target="_blank" className="px-4 py-2 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-xl hover:bg-blue-600 hover:text-white transition-all border border-blue-100">🔗 NCBI Protein</a>}
                {data.ncbi_gene_url && <a href={data.ncbi_gene_url} target="_blank" className="px-4 py-2 bg-slate-50 text-slate-600 text-[9px] font-black uppercase tracking-widest rounded-xl hover:bg-slate-900 hover:text-white transition-all border border-slate-100">🔗 NCBI Gene</a>}
              </div>
            </div>

            {/* Visual Flow Diagram */}
            <div className="xl:col-span-2 bg-slate-900 rounded-[3rem] p-10 text-white relative overflow-hidden shadow-2xl shadow-blue-900/20">
              <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32"></div>
              <h3 className="text-[10px] font-black text-blue-400 uppercase tracking-[0.3em] mb-12 text-center">Protocolo de Flujo Molecular</h3>
              <div className="flex flex-col md:flex-row items-center justify-between gap-10 relative z-10 h-full max-h-[300px]">
                <div className="text-center space-y-4 flex-1">
                  <p className="text-[10px] font-black uppercase text-blue-400">DNA</p>
                  <div className="bg-white/5 border border-white/10 rounded-3xl p-6 hover:bg-white/10 transition-all cursor-default">
                    <p className="text-sm font-black">Doble Hebra</p>
                    <p className="text-[9px] font-bold opacity-40 uppercase tracking-widest mt-1">5'→3' / 3'→5'</p>
                  </div>
                </div>
                <div className="text-blue-500 animate-pulse text-2xl">→</div>
                <div className="text-center space-y-4 flex-1">
                  <p className="text-[10px] font-black uppercase text-indigo-400">Transcripción</p>
                  <div className="bg-indigo-500/20 border border-indigo-500/30 rounded-3xl p-6 hover:bg-indigo-500/30 transition-all cursor-default">
                    <p className="text-sm font-black">RNA Polimerasa</p>
                    <p className="text-[9px] font-bold opacity-40 uppercase tracking-widest mt-1">mRNA 5'→3'</p>
                  </div>
                </div>
                <div className="text-indigo-500 animate-pulse text-2xl">→</div>
                <div className="text-center space-y-4 flex-1">
                  <p className="text-[10px] font-black uppercase text-emerald-400">Traducción</p>
                  <div className="bg-emerald-500/20 border border-emerald-500/30 rounded-3xl p-6 hover:bg-emerald-500/30 transition-all cursor-default">
                    <p className="text-sm font-black">Ribosoma</p>
                    <p className="text-[9px] font-bold opacity-40 uppercase tracking-widest mt-1">Proteína N→C</p>
                  </div>
                </div>
              </div>
            </div>
          </div>

          {/* 3. Alineamiento de Secuencias (Aprovechando Ancho) */}
          <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-12 shadow-sm space-y-16 overflow-hidden">
            <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Alineamiento de Transmisión Genética</h3>
            
            <div className="space-y-12">
              {/* DNA Sense */}
              <div className="space-y-4">
                <p className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em]">DNA Hebra Codificante (Sense) <span className="text-blue-600 ml-4">5' → 3'</span></p>
                <div className="bg-slate-50 rounded-2xl p-8 font-mono text-sm leading-relaxed break-all font-black text-slate-800 tracking-[0.3em] shadow-inner border border-slate-100">
                  <span className="text-blue-600 mr-4 text-xs">5'</span>
                  {data.dna_sense_5to3.split('').map((b, i) => <span key={i} className={BASE_COLORS[b]}>{b}</span>)}
                  <span className="text-blue-600 ml-4 text-xs">3'</span>
                </div>
              </div>

              {/* Interaction lines */}
              <div className="flex justify-center font-mono text-slate-200 text-xs tracking-[0.3em] h-4 select-none overflow-hidden">
                {"|".repeat(Math.min(data.dna_sense_5to3.length, 200))}
              </div>

              {/* DNA Template */}
              <div className="space-y-4">
                <p className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em]">DNA Hebra Molde (Template) <span className="text-indigo-600 ml-4">3' → 5'</span></p>
                <div className="bg-slate-50 rounded-2xl p-8 font-mono text-sm leading-relaxed break-all font-black text-slate-800 tracking-[0.3em] shadow-inner border border-slate-100">
                  <span className="text-indigo-600 mr-4 text-xs">3'</span>
                  {data.dna_template_3to5.split('').map((b, i) => <span key={i} className={BASE_COLORS[b]}>{b}</span>)}
                  <span className="text-indigo-600 ml-4 text-xs">5'</span>
                </div>
              </div>

              <div className="flex flex-col items-center gap-3 py-4">
                <div className="w-px h-12 bg-gradient-to-b from-blue-500 to-indigo-500 animate-pulse"></div>
                <p className="text-[9px] font-black text-blue-600 uppercase tracking-widest bg-blue-50 px-4 py-1.5 rounded-full">Transcripción Molecular</p>
              </div>

              {/* mRNA */}
              <div className="space-y-4">
                <p className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em]">mRNA <span className="text-indigo-400 ml-4">5' → 3'</span></p>
                <div className="bg-indigo-50/30 rounded-2xl p-8 font-mono text-sm leading-relaxed break-all font-black text-slate-800 tracking-[0.3em] shadow-inner border border-indigo-100">
                  <span className="text-indigo-600 mr-4 text-xs">5'</span>
                  {data.mrna_5to3.split('').map((b, i) => <span key={i} className={b === 'U' ? 'text-indigo-500' : BASE_COLORS[b]}>{b}</span>)}
                  <span className="text-indigo-600 ml-4 text-xs">3'</span>
                </div>
              </div>

              <div className="flex flex-col items-center gap-3 py-4">
                <div className="w-px h-12 bg-gradient-to-b from-indigo-500 to-emerald-500 animate-pulse"></div>
                <p className="text-[9px] font-black text-emerald-600 uppercase tracking-widest bg-emerald-50 px-4 py-1.5 rounded-full">Traducción Ribosómica</p>
              </div>

              {/* Proteína (NH2...COOH) */}
              <div className="space-y-4">
                <p className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em]">Cadena Polipeptídica <span className="text-emerald-500 ml-4">NH₂ → COOH</span></p>
                <div className="bg-emerald-50/30 rounded-2xl p-8 font-mono text-sm leading-relaxed break-all font-black text-slate-800 tracking-[0.5em] shadow-inner border border-emerald-100 relative">
                  <div className="absolute top-4 right-8 px-3 py-1 bg-emerald-500 text-white text-[8px] font-black uppercase rounded-lg shadow-lg">Proteína</div>
                  <span className="text-emerald-600 mr-4 text-xs">NH₂</span>
                  {data.full_protein.split('').map((aa, i) => (
                    <span key={i} className={`px-0.5 rounded ${AA_COLORS[aa] || 'bg-slate-400'} text-white mx-0.5 shadow-sm`}>{aa}</span>
                  ))}
                  <span className="text-emerald-600 ml-4 text-xs">COOH</span>
                </div>
              </div>
            </div>
          </div>

          {/* 4. Tabla de Codones */}
          <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
            <div className="p-10 border-b border-slate-50 bg-slate-50/30 flex justify-between items-center">
              <h3 className="text-[10px] font-black text-slate-600 uppercase tracking-[0.3em]">Alineamiento Codón por Codón</h3>
              <button onClick={() => setShowFullTable(!showFullTable)} className="px-8 py-3 bg-slate-900 text-white rounded-2xl text-[9px] font-black uppercase tracking-widest hover:bg-blue-600 transition-all shadow-xl active:scale-95">
                {showFullTable ? 'Ver menos' : `Ver todos (${data.total_codons})`}
              </button>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead className="bg-white border-b border-slate-100">
                  <tr className="text-[9px] font-black text-slate-500 uppercase tracking-widest">
                    <th className="px-10 py-6">#</th>
                    <th className="px-6 py-6 text-blue-600">DNA 5'→3'</th>
                    <th className="px-6 py-6 text-indigo-600">DNA 3'→5'</th>
                    <th className="px-6 py-6 text-indigo-400">mRNA</th>
                    <th className="px-6 py-6 text-emerald-600 text-center">AA</th>
                    <th className="px-10 py-6">Identidad Química</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-50">
                  {(showFullTable ? data.codons : data.codons.slice(0, 24)).map((c, i) => (
                    <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                      <td className="px-10 py-5 font-mono text-[10px] text-slate-400 font-black">{c.position}</td>
                      <td className="px-6 py-5 font-mono text-sm font-black text-blue-600 tracking-widest">{c.dna_codon}</td>
                      <td className="px-6 py-5 font-mono text-sm font-black text-indigo-600 tracking-widest">{c.template_codon}</td>
                      <td className="px-6 py-5 font-mono text-sm font-black text-indigo-400 tracking-widest">{c.rna_codon}</td>
                      <td className="px-6 py-5 text-center">
                        <span className={`w-10 h-10 rounded-xl flex items-center justify-center text-white font-mono text-sm font-black shadow-lg mx-auto transition-transform group-hover:scale-110 ${AA_COLORS[c.amino_acid] || 'bg-slate-400'}`}>
                          {c.amino_acid}
                        </span>
                      </td>
                      <td className="px-10 py-5">
                        <span className={`text-[10px] font-black uppercase ${c.is_start ? 'text-emerald-600' : c.is_stop ? 'text-rose-600' : 'text-slate-700'}`}>
                          {c.is_start && '🟢 '}
                          {c.is_stop && '🛑 '}
                          {c.amino_acid_name}
                          {c.is_start && ' (INICIO)'}
                          {c.is_stop && ' (PARADA)'}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>

          {/* 5. Leyendas Bioquímicas */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h4 className="text-[10px] font-black text-slate-600 uppercase tracking-widest mb-8">Código de Color: Bases</h4>
              <div className="grid grid-cols-2 md:grid-cols-3 gap-4">
                {[
                  { b: 'A', label: 'Adenina', c: 'bg-emerald-400' },
                  { b: 'T', label: 'Timina', c: 'bg-rose-400' },
                  { b: 'G', label: 'Guanina', c: 'bg-blue-400' },
                  { b: 'C', label: 'Citosina', c: 'bg-amber-400' },
                  { b: 'U', label: 'Uracilo', c: 'bg-indigo-400' }
                ].map(item => (
                  <div key={item.b} className="flex items-center gap-4 p-3 bg-slate-50 rounded-2xl border border-slate-100">
                    <div className={`w-8 h-8 rounded-xl ${item.c} flex items-center justify-center text-white font-black text-xs shadow-md`}>{item.b}</div>
                    <span className="text-[10px] font-black text-slate-700 uppercase">{item.label}</span>
                  </div>
                ))}
              </div>
            </div>
            <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm">
              <h4 className="text-[10px] font-black text-slate-600 uppercase tracking-widest mb-8">Propiedades de Aminoácidos</h4>
              <div className="space-y-6">
                {[
                  { label: 'Hidrofóbicos', c: 'bg-emerald-500', desc: 'A, L, I, V, M, F, W, P', text: 'Residuos no polares que tienden a ocultarse del agua.' },
                  { label: 'Polares / Básicos', c: 'bg-blue-500', desc: 'S, T, N, Q, Y, C, G, K, R, H', text: 'Residuos con afinidad por el agua o carga positiva.' },
                  { label: 'Cargados Negativos', c: 'bg-rose-500', desc: 'D, E', text: 'Residuos con carga ácida bajo condiciones fisiológicas.' }
                ].map(item => (
                  <div key={item.label} className="flex gap-6 group">
                    <div className={`w-3 h-3 rounded-full ${item.c} mt-1 flex-shrink-0 group-hover:scale-125 transition-transform shadow-sm`}></div>
                    <div className="space-y-1">
                      <p className="text-[10px] font-black text-slate-900 uppercase tracking-tight">{item.label}</p>
                      <p className="text-[9px] font-bold text-blue-600 font-mono">{item.desc}</p>
                      <p className="text-[9px] text-slate-500 font-medium leading-relaxed">{item.text}</p>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}