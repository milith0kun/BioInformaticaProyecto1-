/**
 * BLASTSearch Component
 * Submit BLAST searches against NCBI databases directly from the app.
 * Supports blastn/blastp, shows results with alignment stats.
 */
import { useState, useEffect, useRef } from 'react'
import api from '../services/api'

export default function BLASTSearch() {
    const [sequence, setSequence] = useState('')
    const [program, setProgram] = useState('blastn')
    const [database, setDatabase] = useState('nt')
    const [maxHits, setMaxHits] = useState(10)
    const [submitting, setSubmitting] = useState(false)
    const [polling, setPolling] = useState(false)
    const [rid, setRid] = useState(null)
    const [results, setResults] = useState(null)
    const [error, setError] = useState(null)
    const [proteins, setProteins] = useState([])
    const pollRef = useRef(null)

    // Load available proteins for quick sequence selection
    useEffect(() => {
        loadProteins()
        return () => { if (pollRef.current) clearInterval(pollRef.current) }
    }, [])

    const loadProteins = async () => {
        try {
            const data = await api.getProteins()
            setProteins(data?.proteins?.slice(0, 30) || [])
        } catch (e) { /* ok */ }
    }

    const handleSubmit = async () => {
        const cleanSeq = sequence.replace(/[^A-Za-z\n]/g, '').replace(/\n/g, '')
        if (cleanSeq.length < 10) {
            setError('La secuencia debe tener al menos 10 bases/aminoácidos')
            return
        }

        setSubmitting(true)
        setError(null)
        setResults(null)
        setRid(null)

        try {
            const resp = await api.submitBLAST(cleanSeq, program, database, maxHits)
            if (resp.rid) {
                setRid(resp.rid)
                setSubmitting(false)
                startPolling(resp.rid)
            } else {
                setError('No se recibió RID de NCBI')
                setSubmitting(false)
            }
        } catch (e) {
            setError(e.response?.data?.detail || 'Error enviando BLAST')
            setSubmitting(false)
        }
    }

    const startPolling = (blastRid) => {
        setPolling(true)
        let attempts = 0
        pollRef.current = setInterval(async () => {
            attempts++
            try {
                const resp = await api.getBLASTResults(blastRid)
                if (resp.status === 'COMPLETE') {
                    setResults(resp)
                    setPolling(false)
                    clearInterval(pollRef.current)
                } else if (attempts > 60) { // 5 min timeout
                    setPolling(false)
                    clearInterval(pollRef.current)
                    setResults({
                        status: 'TIMEOUT',
                        blast_url: `https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${blastRid}&FORMAT_TYPE=HTML`,
                        rid: blastRid
                    })
                }
            } catch (e) {
                // continue polling
            }
        }, 5000)
    }

    const loadGeneSequence = async (locusTag) => {
        try {
            const detail = await api.getCentralDogma(locusTag)
            if (detail?.dna_sequence) {
                setSequence(detail.dna_sequence)
                setProgram('blastn')
                setDatabase('nt')
            } else if (detail?.protein_sequence) {
                setSequence(detail.protein_sequence)
                setProgram('blastp')
                setDatabase('nr')
            }
        } catch (e) {
            setError('Error cargando secuencia del gen')
        }
    }

    const PROGRAMS = [
        { id: 'blastn', label: 'blastn', desc: 'Nucleotide → Nucleotide' },
        { id: 'blastp', label: 'blastp', desc: 'Protein → Protein' },
        { id: 'blastx', label: 'blastx', desc: 'Translated Nucleotide → Protein' },
        { id: 'tblastn', label: 'tblastn', desc: 'Protein → Translated Nucleotide' },
    ]

    const DATABASES = {
        blastn: [
            { id: 'nt', label: 'nt (GenBank core nucleotides)' },
            { id: 'refseq_genomic', label: 'refseq_genomic' },
            { id: 'est', label: 'est (Expressed Sequence Tags)' },
        ],
        blastp: [
            { id: 'nr', label: 'nr (Non-redundant proteins)' },
            { id: 'swissprot', label: 'Swiss-Prot' },
            { id: 'refseq_protein', label: 'RefSeq Proteins' },
        ],
        blastx: [
            { id: 'nr', label: 'nr (Non-redundant proteins)' },
            { id: 'swissprot', label: 'Swiss-Prot' },
        ],
        tblastn: [
            { id: 'nt', label: 'nt (GenBank core nucleotides)' },
            { id: 'refseq_genomic', label: 'refseq_genomic' },
        ],
    }

    return (
        <div className="space-y-6">
            <div>
                <h2 className="text-xl font-bold text-slate-800">🔍 BLAST — Búsqueda de Secuencias</h2>
                <p className="text-sm text-slate-500">
                    Compare secuencias contra las bases de datos de NCBI sin salir de la aplicación
                </p>
            </div>

            {/* Input Area */}
            <div className="bg-white rounded-xl border border-slate-200 p-5">
                <div className="flex items-center justify-between mb-3">
                    <h3 className="font-semibold text-slate-800 text-sm">Secuencia de Consulta</h3>
                    {proteins.length > 0 && (
                        <div className="flex items-center gap-2">
                            <span className="text-xs text-slate-500">Cargar gen:</span>
                            <select
                                onChange={(e) => e.target.value && loadGeneSequence(e.target.value)}
                                className="px-2 py-1 border border-slate-200 rounded-lg text-xs max-w-[200px]"
                                defaultValue=""
                            >
                                <option value="">Seleccionar gen...</option>
                                {proteins.map(p => (
                                    <option key={p.locus_tag} value={p.locus_tag}>
                                        {p.locus_tag} — {(p.product || '').slice(0, 40)}
                                    </option>
                                ))}
                            </select>
                        </div>
                    )}
                </div>

                <textarea
                    value={sequence}
                    onChange={(e) => setSequence(e.target.value)}
                    placeholder="Pegue su secuencia aquí (DNA o proteína, formato FASTA o texto plano)...&#10;&#10;Ejemplo:&#10;ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
                    className="w-full h-36 p-3 font-mono text-xs bg-slate-50 border border-slate-200 rounded-lg resize-y focus:ring-2 focus:ring-teal-500 focus:border-teal-500 outline-none"
                />
                <div className="flex items-center justify-between mt-2 text-xs text-slate-500">
                    <span>{sequence.replace(/[^A-Za-z]/g, '').length} caracteres</span>
                    <button onClick={() => setSequence('')} className="text-red-500 hover:text-red-600">Limpiar</button>
                </div>
            </div>

            {/* Config */}
            <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <label className="text-xs font-medium text-slate-600 block mb-2">Programa BLAST</label>
                    <div className="space-y-1.5">
                        {PROGRAMS.map(p => (
                            <label key={p.id} className={`flex items-start gap-2 p-2 rounded-lg cursor-pointer transition-all ${program === p.id ? 'bg-teal-50 border border-teal-200' : 'hover:bg-slate-50'}`}>
                                <input type="radio" name="program" value={p.id} checked={program === p.id}
                                    onChange={() => {
                                        setProgram(p.id)
                                        setDatabase(DATABASES[p.id]?.[0]?.id || 'nt')
                                    }}
                                    className="mt-0.5 accent-teal-600" />
                                <div>
                                    <span className="text-sm font-medium text-slate-800">{p.label}</span>
                                    <p className="text-[10px] text-slate-500">{p.desc}</p>
                                </div>
                            </label>
                        ))}
                    </div>
                </div>

                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <label className="text-xs font-medium text-slate-600 block mb-2">Base de Datos</label>
                    <div className="space-y-1.5">
                        {(DATABASES[program] || []).map(db => (
                            <label key={db.id} className={`flex items-center gap-2 p-2 rounded-lg cursor-pointer transition-all ${database === db.id ? 'bg-teal-50 border border-teal-200' : 'hover:bg-slate-50'}`}>
                                <input type="radio" name="database" value={db.id} checked={database === db.id}
                                    onChange={() => setDatabase(db.id)} className="accent-teal-600" />
                                <span className="text-sm text-slate-700">{db.label}</span>
                            </label>
                        ))}
                    </div>
                </div>

                <div className="bg-white rounded-xl border border-slate-200 p-4">
                    <label className="text-xs font-medium text-slate-600 block mb-2">Opciones</label>
                    <div className="space-y-3">
                        <div>
                            <label className="text-xs text-slate-500">Max Hits</label>
                            <select value={maxHits} onChange={(e) => setMaxHits(parseInt(e.target.value))}
                                className="w-full px-2 py-1.5 border border-slate-200 rounded-lg text-sm mt-1">
                                <option value={5}>5</option>
                                <option value={10}>10</option>
                                <option value={20}>20</option>
                                <option value={50}>50</option>
                            </select>
                        </div>
                    </div>

                    <button
                        onClick={handleSubmit}
                        disabled={submitting || polling || sequence.length < 10}
                        className="w-full mt-4 px-4 py-3 bg-teal-600 text-white rounded-lg font-medium text-sm
                            hover:bg-teal-700 disabled:opacity-50 disabled:cursor-not-allowed transition-all"
                    >
                        {submitting ? '⏳ Enviando...' : polling ? '🔄 Esperando resultados...' : '🔍 Ejecutar BLAST'}
                    </button>
                </div>
            </div>

            {error && (
                <div className="bg-red-50 border border-red-200 rounded-xl p-4 text-red-700 text-sm">
                    ❌ {error}
                </div>
            )}

            {/* Polling Status */}
            {polling && (
                <div className="bg-amber-50 border border-amber-200 rounded-xl p-6 text-center">
                    <div className="w-12 h-12 border-3 border-amber-200 border-t-amber-600 rounded-full animate-spin mx-auto mb-4"></div>
                    <p className="font-semibold text-amber-800">Búsqueda BLAST en progreso...</p>
                    <p className="text-sm text-amber-700 mt-1">RID: <code className="bg-amber-100 px-1 rounded">{rid}</code></p>
                    <p className="text-xs text-amber-600 mt-2">
                        Las búsquedas BLAST pueden tomar entre 30 segundos y varios minutos.
                        Los resultados se actualizarán automáticamente.
                    </p>
                    <a href={`https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${rid}&FORMAT_TYPE=HTML`}
                        target="_blank" rel="noopener noreferrer"
                        className="inline-block mt-3 text-xs text-amber-700 underline hover:text-amber-900">
                        Ver en NCBI BLAST →
                    </a>
                </div>
            )}

            {/* Results */}
            {results && (
                <div className="space-y-4">
                    {results.status === 'TIMEOUT' ? (
                        <div className="bg-amber-50 border border-amber-200 rounded-xl p-5 text-center">
                            <p className="font-semibold text-amber-800">⏰ Tiempo de espera excedido</p>
                            <p className="text-sm text-amber-700 mt-2">
                                La búsqueda BLAST continúa en NCBI. Consulte los resultados directamente:
                            </p>
                            <a href={results.blast_url} target="_blank" rel="noopener noreferrer"
                                className="inline-block mt-3 px-4 py-2 bg-amber-600 text-white rounded-lg text-sm hover:bg-amber-700">
                                🔗 Ver resultados en NCBI BLAST
                            </a>
                        </div>
                    ) : (
                        <>
                            <div className="flex items-center justify-between">
                                <h3 className="font-semibold text-slate-800">
                                    Resultados ({results.total_hits} hits)
                                </h3>
                                <a href={results.blast_url} target="_blank" rel="noopener noreferrer"
                                    className="text-xs text-teal-600 hover:text-teal-800 underline">
                                    Ver en NCBI BLAST →
                                </a>
                            </div>

                            {results.hits && results.hits.length > 0 ? (
                                <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
                                    <div className="overflow-x-auto">
                                        <table className="w-full text-sm">
                                            <thead className="bg-slate-50">
                                                <tr>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">#</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Accession</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500 min-w-[200px]">Descripción</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Score</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">E-value</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Identidad</th>
                                                    <th className="text-left px-3 py-2 text-xs font-medium text-slate-500">Longitud</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {results.hits.map((hit, i) => (
                                                    <tr key={i} className={`border-b border-slate-100 ${i % 2 ? 'bg-slate-50/50' : ''} hover:bg-teal-50/50 transition-colors`}>
                                                        <td className="px-3 py-2 text-xs text-slate-400">{i + 1}</td>
                                                        <td className="px-3 py-2">
                                                            <a href={hit.url} target="_blank" rel="noopener noreferrer"
                                                                className="text-teal-600 hover:text-teal-800 font-mono text-xs underline">
                                                                {hit.accession}
                                                            </a>
                                                        </td>
                                                        <td className="px-3 py-2 text-xs text-slate-700">
                                                            {hit.title?.slice(0, 120)}{hit.title?.length > 120 ? '...' : ''}
                                                        </td>
                                                        <td className="px-3 py-2 font-mono text-xs font-bold text-slate-800">
                                                            {typeof hit.score === 'number' ? hit.score.toFixed(1) : hit.score}
                                                        </td>
                                                        <td className="px-3 py-2 font-mono text-xs">
                                                            <span className={`px-1.5 py-0.5 rounded ${hit.evalue < 1e-50 ? 'bg-green-100 text-green-700' :
                                                                hit.evalue < 1e-10 ? 'bg-teal-100 text-teal-700' :
                                                                    'bg-amber-100 text-amber-700'
                                                                }`}>
                                                                {typeof hit.evalue === 'number' ? hit.evalue.toExponential(1) : hit.evalue}
                                                            </span>
                                                        </td>
                                                        <td className="px-3 py-2 text-xs">
                                                            <span className={`font-bold ${hit.identity_pct >= 99 ? 'text-green-600' :
                                                                hit.identity_pct >= 90 ? 'text-teal-600' :
                                                                    hit.identity_pct >= 70 ? 'text-amber-600' : 'text-red-600'
                                                                }`}>
                                                                {hit.identity_pct}%
                                                            </span>
                                                        </td>
                                                        <td className="px-3 py-2 font-mono text-xs text-slate-600">
                                                            {hit.align_len?.toLocaleString()} bp
                                                        </td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            ) : (
                                <div className="bg-slate-50 rounded-xl p-8 text-center">
                                    <p className="text-slate-500">No se encontraron hits significativos.</p>
                                    <a href={results.blast_url} target="_blank" rel="noopener noreferrer"
                                        className="text-sm text-teal-600 underline mt-2 inline-block">
                                        Ver resultados completos en NCBI →
                                    </a>
                                </div>
                            )}
                        </>
                    )}
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
                        <p><strong>BLAST</strong> (Basic Local Alignment Search Tool) compara secuencias de nucleótidos o proteínas contra las bases de datos de NCBI para encontrar regiones de similitud.</p>
                        <p className="text-xs text-teal-700">
                            <strong>E-value</strong>: Probabilidad de encontrar un hit igual o mejor por azar. Valores menores indican mayor significancia.
                            <strong> Score</strong>: Puntuación del alineamiento (bits). Mayor score = mejor alineamiento.
                            <strong> Identidad</strong>: Porcentaje de bases/residuos idénticos entre query y subject.
                        </p>
                    </div>
                </div>
            </div>
        </div>
    )
}
