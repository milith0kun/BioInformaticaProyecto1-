/**
 * BLASTSearch Component — Clean Laboratory Edition
 * Advanced NCBI BLAST integration with sequence inspection and educational context
 */
import { useState, useEffect } from 'react'
import { api } from '../services/api'
import toast from 'react-hot-toast'

export default function BLASTSearch({ genes = [] }) {
  const [sequence, setSequence] = useState('')
  const [program, setProgram] = useState('blastn')
  const [database, setDatabase] = useState('nt')
  const [maxHits, setMaxHits] = useState(10)
  const [loading, setLoading] = useState(false)
  const [results, setResults] = useState(null)
  const [rid, setRid] = useState(null)
  const [polling, setPolling] = useState(false)

  // Character count
  const charCount = sequence.trim().length

  const handleClear = () => {
    setSequence('')
    setResults(null)
    setRid(null)
  }

  const handleLoadGene = (e) => {
    const locusTag = e.target.value
    if (!locusTag) return
    
    // In a real scenario, we might need to fetch the actual sequence from the backend
    // Since we don't have all sequences in the 'genes' prop, we'll show a placeholder
    // or if the backend supports it, fetch by locus_tag.
    // For now, we'll try to find if we have it or just toast.
    toast.loading('Cargando secuencia del gen...', { id: 'load-gene', duration: 1000 })
    
    // Simulando la carga de secuencia (en una implementación real llamaríamos a api.getGeneDetail o similar)
    api.getCentralDogma(locusTag).then(data => {
      if (data && data.dna_sequence) {
        setSequence(data.dna_sequence)
        setProgram('blastn')
        toast.success('Secuencia cargada correctamente', { id: 'load-gene' })
      }
    }).catch(() => {
      toast.error('No se pudo cargar la secuencia completa', { id: 'load-gene' })
    })
  }

  const runBlast = async (e) => {
    e.preventDefault()
    if (charCount < 10) {
      toast.error('Secuencia demasiado corta (mínimo 10 caracteres)')
      return
    }

    setLoading(true)
    setResults(null)
    setRid(null)

    try {
      toast.loading('Enviando consulta a NCBI...', { id: 'blast' })
      const response = await api.submitBLAST(sequence.trim(), program, database, maxHits)
      
      if (response.rid) {
        setRid(response.rid)
        toast.loading(`Consulta aceptada. RID: ${response.rid}. Esperando resultados...`, { id: 'blast' })
        startPolling(response.rid)
      } else {
        throw new Error('No se recibió ID de rastreo (RID)')
      }
    } catch (err) {
      console.error(err)
      toast.error('Error al iniciar BLAST: ' + (err.response?.data?.detail || err.message), { id: 'blast' })
      setLoading(false)
    }
  }

  const startPolling = (requestid) => {
    setPolling(true)
    let attempts = 0
    const interval = setInterval(async () => {
      attempts++
      try {
        const data = await api.getBLASTResults(requestid)
        
        if (data.status === 'READY') {
          clearInterval(interval)
          setResults(data)
          setLoading(false)
          setPolling(false)
          toast.success('Búsqueda BLAST completada con éxito', { id: 'blast' })
        } else if (data.status === 'ERROR') {
          clearInterval(interval)
          setLoading(false)
          setPolling(false)
          toast.error('NCBI informó un error en el procesamiento', { id: 'blast' })
        } else if (attempts > 30) { // 30 * 5s = 2.5 minutes
          clearInterval(interval)
          setLoading(false)
          setPolling(false)
          toast.error('Tiempo de espera agotado. Intente revisar el RID más tarde.', { id: 'blast' })
        }
      } catch (e) {
        console.error('Polling error:', e)
      }
    }, 5000)
  }

  const programs = [
    { id: 'blastn', name: 'blastn', desc: 'Nucleotide → Nucleotide' },
    { id: 'blastp', name: 'blastp', desc: 'Protein → Protein' },
    { id: 'blastx', name: 'blastx', desc: 'Translated Nucleotide → Protein' },
    { id: 'tblastn', name: 'tblastn', desc: 'Protein → Translated Nucleotide' }
  ]

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Header Module */}
      <div className="bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="relative z-10 flex flex-col md:flex-row md:items-center justify-between gap-8">
          <div className="space-y-2">
            <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
              🔍 BLAST — <span className="text-blue-600">Búsqueda de Secuencias</span>
            </h2>
            <p className="text-[10px] font-bold text-slate-400 uppercase tracking-[0.2em]">Compare secuencias contra las bases de datos de NCBI sin salir de la aplicación</p>
          </div>
        </div>
      </div>

      <div className="grid grid-cols-1 xl:grid-cols-3 gap-8">
        {/* Input Section */}
        <div className="xl:col-span-2 space-y-8">
          <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm space-y-8">
            <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
              <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">Secuencia de Consulta</h3>
              <div className="flex items-center gap-4 bg-slate-50 px-4 py-2 rounded-2xl border border-slate-100 shadow-inner">
                <span className="text-[9px] font-black text-slate-400 uppercase">Cargar gen:</span>
                <select 
                  onChange={handleLoadGene}
                  className="bg-transparent text-[10px] font-black text-blue-600 uppercase focus:outline-none max-w-[150px] truncate"
                >
                  <option value="">Seleccionar gen...</option>
                  {genes.map(g => (
                    <option key={g.locus_tag} value={g.locus_tag}>{g.gene_name || g.locus_tag}</option>
                  ))}
                </select>
              </div>
            </div>

            <div className="relative">
              <textarea
                value={sequence}
                onChange={(e) => setSequence(e.target.value)}
                placeholder="Pegue su secuencia aquí (DNA o proteína, formato FASTA o texto plano)..."
                className="w-full h-64 p-8 bg-slate-50 border-2 border-slate-100 rounded-[2.5rem] font-mono text-sm font-black tracking-widest text-slate-700 focus:outline-none focus:border-blue-500/50 transition-all shadow-inner placeholder-slate-300 uppercase"
              />
              <div className="absolute bottom-6 right-8 flex items-center gap-6">
                <span className="text-[10px] font-black text-slate-400 uppercase tracking-widest">{charCount} caracteres</span>
                <button 
                  onClick={handleClear}
                  className="text-[9px] font-black text-rose-500 uppercase tracking-widest hover:text-rose-600 transition-colors"
                >
                  Limpiar
                </button>
              </div>
            </div>

            <div className="p-6 bg-blue-50/50 rounded-2xl border border-blue-100">
              <p className="text-[9px] font-black text-blue-400 uppercase tracking-widest mb-2">Ejemplo:</p>
              <p className="font-mono text-[10px] text-blue-700/60 break-all leading-relaxed">
                ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA
              </p>
            </div>
          </div>

          {/* Results Module */}
          {loading && (
            <div className="bg-slate-900 rounded-[3rem] p-20 text-center space-y-8 shadow-2xl shadow-blue-900/40 relative overflow-hidden">
              <div className="absolute inset-0 bg-[radial-gradient(circle_at_center,_var(--tw-gradient-stops))] from-blue-500/10 via-transparent to-transparent"></div>
              <div className="w-20 h-20 border-4 border-white/5 border-t-blue-500 rounded-full animate-spin mx-auto relative z-10"></div>
              <div className="space-y-2 relative z-10">
                <p className="text-[10px] font-black text-blue-400 uppercase tracking-[0.4em] animate-pulse">Procesando Alineamiento Global</p>
                <p className="text-[9px] font-bold text-slate-500 uppercase tracking-widest">NCBI está analizando su secuencia... RID: {rid || 'Asignando'}</p>
              </div>
            </div>
          )}

          {results && (
            <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm animate-in slide-in-from-bottom-8 duration-700">
              <div className="p-8 border-b border-slate-50 bg-slate-50/50 flex items-center justify-between">
                <div>
                  <h3 className="text-[10px] font-black text-slate-900 uppercase tracking-[0.3em]">Resultados de Similitud Detectados</h3>
                  <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest mt-1">Fuente: NCBI Blast Database</p>
                </div>
                <div className="px-5 py-2 bg-blue-600 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-lg">
                  {results.hits?.length || 0} Aciertos
                </div>
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-left">
                  <thead>
                    <tr className="bg-white border-b border-slate-100">
                      <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Identificador</th>
                      <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Score</th>
                      <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">E-Value</th>
                      <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Identidad</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-50">
                    {results.hits?.map((hit, i) => (
                      <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                        <td className="px-8 py-6">
                          <div className="space-y-1">
                            <p className="font-mono text-xs font-black text-blue-600 group-hover:text-blue-700 transition-colors uppercase">{hit.accession}</p>
                            <p className="text-[10px] text-slate-500 font-bold uppercase truncate max-w-md">{hit.title}</p>
                          </div>
                        </td>
                        <td className="px-8 py-6 text-right font-mono text-xs font-black text-slate-700">{hit.score}</td>
                        <td className="px-8 py-6 text-right font-mono text-xs font-black text-rose-600">{hit.evalue}</td>
                        <td className="px-8 py-6 text-right">
                          <div className="flex items-center justify-end gap-4">
                            <span className="text-xs font-black text-slate-900">{(hit.identity_pct || 0).toFixed(1)}%</span>
                            <div className="w-20 h-1 bg-slate-100 rounded-full overflow-hidden">
                              <div className="h-full bg-blue-500 transition-all duration-1000 shadow-[0_0_8px_rgba(59,130,246,0.4)]" style={{ width: `${hit.identity_pct}%` }}></div>
                            </div>
                          </div>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>

        {/* Configuration Sidebar */}
        <div className="space-y-8">
          <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-8 shadow-sm space-y-10">
            <div className="space-y-6">
              <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-2">Programa BLAST</h4>
              <div className="grid grid-cols-1 gap-3">
                {programs.map(p => (
                  <button
                    key={p.id}
                    onClick={() => setProgram(p.id)}
                    className={`p-5 rounded-2xl border-2 text-left transition-all ${
                      program === p.id 
                        ? 'border-blue-500 bg-blue-50 shadow-xl shadow-blue-500/5' 
                        : 'border-slate-50 bg-slate-50 hover:border-slate-200'
                    }`}
                  >
                    <p className={`text-xs font-black uppercase tracking-widest ${program === p.id ? 'text-blue-600' : 'text-slate-700'}`}>{p.name}</p>
                    <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight mt-1">{p.desc}</p>
                  </button>
                ))}
              </div>
            </div>

            <div className="space-y-6">
              <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-2">Base de Datos</h4>
              <div className="space-y-2">
                {[
                  { id: 'nt', label: 'nt (GenBank core nucleotides)' },
                  { id: 'refseq_genomic', label: 'refseq_genomic' },
                  { id: 'est', label: 'est (Expressed Sequence Tags)' }
                ].map(db => (
                  <button
                    key={db.id}
                    onClick={() => setDatabase(db.id)}
                    className={`w-full p-4 rounded-xl border-2 text-left text-[10px] font-black transition-all ${
                      database === db.id 
                        ? 'border-blue-200 bg-blue-50/50 text-blue-600' 
                        : 'border-slate-50 text-slate-400 hover:border-slate-100 hover:text-slate-600'
                    }`}
                  >
                    {db.label}
                  </button>
                ))}
              </div>
            </div>

            <div className="space-y-6">
              <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em] px-2">Opciones</h4>
              <div className="p-6 bg-slate-50 rounded-[2rem] border border-slate-100 flex items-center justify-between">
                <span className="text-[10px] font-black text-slate-500 uppercase tracking-widest">Max Hits</span>
                <input 
                  type="number" 
                  value={maxHits} 
                  onChange={(e) => setMaxHits(parseInt(e.target.value))}
                  className="w-16 bg-white border-2 border-slate-100 rounded-xl px-3 py-1.5 text-xs font-black text-blue-600 focus:outline-none focus:border-blue-500 transition-all"
                />
              </div>
            </div>

            <button 
              onClick={runBlast}
              disabled={loading}
              className="w-full py-6 bg-slate-900 text-white rounded-[2rem] font-black uppercase tracking-[0.25em] text-[10px] shadow-2xl hover:bg-blue-600 transition-all active:scale-95 disabled:opacity-50 group flex items-center justify-center gap-4"
            >
              <span>🔍 Ejecutar BLAST</span>
              <div className="w-1 h-1 rounded-full bg-white animate-pulse"></div>
            </button>
          </div>

          <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-10 shadow-sm space-y-6">
            <h4 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.2em]">Contexto Científico</h4>
            <p className="text-[11px] text-slate-500 font-medium leading-relaxed uppercase tracking-tight">
              BLAST (Basic Local Alignment Search Tool) compara secuencias de nucleótidos o proteínas contra las bases de datos de NCBI para encontrar regiones de similitud.
            </p>
            <div className="space-y-4 pt-4 border-t border-slate-50">
              <div className="space-y-1">
                <p className="text-[9px] font-black text-slate-900 uppercase">E-value:</p>
                <p className="text-[9px] text-slate-400 font-bold uppercase leading-relaxed">Probabilidad de encontrar un hit igual o mejor por azar. Valores menores indican mayor significancia.</p>
              </div>
              <div className="space-y-1">
                <p className="text-[9px] font-black text-slate-900 uppercase">Score:</p>
                <p className="text-[9px] text-slate-400 font-bold uppercase leading-relaxed">Puntuación del alineamiento (bits). Mayor score = mejor alineamiento.</p>
              </div>
              <div className="space-y-1">
                <p className="text-[9px] font-black text-slate-900 uppercase">Identidad:</p>
                <p className="text-[9px] text-slate-400 font-bold uppercase leading-relaxed">Porcentaje de bases/residuos idénticos entre query y subject.</p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}