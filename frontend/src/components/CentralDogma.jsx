/**
 * CentralDogma Component
 * Shows the complete flow: DNA (5'→3' / 3'→5') → mRNA (transcription) → Protein (translation)
 * Interactive visualization of codon-by-codon translation with NCBI references
 */
import { useState, useEffect } from 'react'
import { api } from '../services/api'

// Amino acid color palette by biochemical property
const AA_COLORS = {
    // Hydrophobic
    A: '#4a9eff', L: '#4a9eff', I: '#4a9eff', V: '#4a9eff', M: '#4a9eff', F: '#4a9eff', W: '#4a9eff', P: '#4a9eff',
    // Polar uncharged
    S: '#50c878', T: '#50c878', N: '#50c878', Q: '#50c878', Y: '#50c878', C: '#50c878', G: '#50c878',
    // Positively charged
    K: '#ff6b6b', R: '#ff6b6b', H: '#ff6b6b',
    // Negatively charged
    D: '#ffd93d', E: '#ffd93d',
    // Special
    '*': '#999',
}

const BASE_COLORS = {
    A: '#22c55e', T: '#ef4444', G: '#f59e0b', C: '#3b82f6',
    U: '#a855f7',
}

export default function CentralDogma() {
    const [data, setData] = useState(null)
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [searchTag, setSearchTag] = useState('')
    const [proteins, setProteins] = useState([])
    const [selectedGene, setSelectedGene] = useState(null)
    const [showFullProtein, setShowFullProtein] = useState(false)
    const [visibleCodons, setVisibleCodons] = useState(20)

    // Load first proteins for selection
    useEffect(() => {
        loadProteins()
    }, [])

    const loadProteins = async () => {
        try {
            const result = await api.getProteins(1, 20)
            setProteins(result.proteins || [])
        } catch (e) {
            console.error('Error loading proteins:', e)
        }
    }

    const loadDogmaData = async (locusTag) => {
        setLoading(true)
        setError(null)
        try {
            const result = await api.getCentralDogma(locusTag)
            setData(result)
            setSelectedGene(locusTag)
        } catch (e) {
            setError(e.response?.data?.detail || 'Error cargando datos del dogma central')
        } finally {
            setLoading(false)
        }
    }

    const handleSearch = (e) => {
        e.preventDefault()
        if (searchTag.trim()) {
            loadDogmaData(searchTag.trim())
        }
    }

    const renderColoredDNA = (seq, isRNA = false) => {
        if (!seq) return null
        return seq.split('').map((base, i) => (
            <span key={i} style={{ color: BASE_COLORS[base] || '#666', fontWeight: 500 }}>
                {base}
            </span>
        ))
    }

    const renderColoredProtein = (seq) => {
        if (!seq) return null
        return seq.split('').map((aa, i) => (
            <span
                key={i}
                style={{
                    color: '#fff',
                    backgroundColor: AA_COLORS[aa] || '#888',
                    borderRadius: '2px',
                    padding: '0 2px',
                    margin: '0 0.5px',
                    fontSize: '0.8rem',
                }}
            >
                {aa}
            </span>
        ))
    }

    return (
        <div className="space-y-6">
            {/* Header */}
            <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4">
                <div>
                    <h2 className="text-xl font-bold text-slate-800 flex items-center gap-2">
                        🧬 Dogma Central de la Biología Molecular
                    </h2>
                    <p className="text-sm text-slate-500 mt-1">
                        DNA → mRNA (Transcripción) → Proteína (Traducción)
                    </p>
                </div>
            </div>

            {/* Gene Selector */}
            <div className="bg-white rounded-xl border border-slate-200 p-5">
                <h3 className="text-sm font-semibold text-slate-700 mb-3">Seleccionar Gen</h3>

                {/* Search by locus tag */}
                <form onSubmit={handleSearch} className="flex gap-2 mb-4">
                    <input
                        type="text"
                        value={searchTag}
                        onChange={(e) => setSearchTag(e.target.value)}
                        placeholder="Ingrese locus_tag (ej: b0001)"
                        className="flex-1 px-3 py-2 border border-slate-200 rounded-lg text-sm focus:outline-none focus:border-teal-400"
                    />
                    <button
                        type="submit"
                        disabled={loading}
                        className="px-4 py-2 bg-teal-600 text-white rounded-lg text-sm font-medium hover:bg-teal-700 disabled:opacity-50"
                    >
                        Buscar
                    </button>
                </form>

                {/* Quick select from proteins */}
                {proteins.length > 0 && (
                    <div>
                        <p className="text-xs text-slate-500 mb-2">O seleccione un gen:</p>
                        <div className="flex flex-wrap gap-2 max-h-32 overflow-y-auto">
                            {proteins.slice(0, 30).map((p) => (
                                <button
                                    key={p.locus_tag}
                                    onClick={() => {
                                        setSearchTag(p.locus_tag)
                                        loadDogmaData(p.locus_tag)
                                    }}
                                    className={`px-2 py-1 rounded text-xs border transition-all ${selectedGene === p.locus_tag
                                            ? 'bg-teal-600 text-white border-teal-600'
                                            : 'bg-slate-50 text-slate-600 border-slate-200 hover:border-teal-300'
                                        }`}
                                >
                                    <span className="font-mono font-medium">{p.gene_name || p.locus_tag}</span>
                                    {p.product && (
                                        <span className="text-slate-400 ml-1">({p.product.substring(0, 20)})</span>
                                    )}
                                </button>
                            ))}
                        </div>
                    </div>
                )}
            </div>

            {/* Loading */}
            {loading && (
                <div className="text-center py-12">
                    <div className="w-12 h-12 mx-auto border-3 border-teal-200 border-t-teal-600 rounded-full animate-spin mb-4"></div>
                    <p className="text-slate-500">Extrayendo datos del dogma central...</p>
                </div>
            )}

            {/* Error */}
            {error && (
                <div className="bg-red-50 border border-red-200 rounded-xl p-4 text-red-700 text-sm">
                    {error}
                </div>
            )}

            {/* Results */}
            {data && !loading && (
                <div className="space-y-5">
                    {/* Gene Info Card */}
                    <div className="bg-gradient-to-r from-slate-800 to-slate-900 rounded-xl p-5 text-white">
                        <div className="flex flex-col sm:flex-row justify-between gap-4">
                            <div>
                                <h3 className="text-lg font-bold flex items-center gap-2">
                                    <span className="text-teal-400">⬡</span>
                                    {data.gene_name || data.locus_tag}
                                    <span className="text-sm font-normal text-slate-400">({data.locus_tag})</span>
                                </h3>
                                <p className="text-slate-300 text-sm mt-1">{data.product}</p>
                            </div>
                            <div className="flex flex-wrap gap-3 text-xs">
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">Posición</div>
                                    <div className="font-mono font-bold">{data.start?.toLocaleString()}-{data.end?.toLocaleString()}</div>
                                </div>
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">Hebra</div>
                                    <div className="font-bold">{data.strand_label}</div>
                                </div>
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">DNA</div>
                                    <div className="font-mono font-bold">{data.length_bp?.toLocaleString()} bp</div>
                                </div>
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">Proteína</div>
                                    <div className="font-mono font-bold">{data.length_aa} aa</div>
                                </div>
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">GC%</div>
                                    <div className="font-bold">{data.gc_content}%</div>
                                </div>
                                <div className="bg-slate-700 rounded-lg px-3 py-2 text-center">
                                    <div className="text-slate-400">Peso</div>
                                    <div className="font-bold">{data.molecular_weight_kda} kDa</div>
                                </div>
                            </div>
                        </div>

                        {/* NCBI Links */}
                        <div className="mt-4 flex flex-wrap gap-2">
                            {data.ncbi_protein_url && (
                                <a href={data.ncbi_protein_url} target="_blank" rel="noopener noreferrer"
                                    className="inline-flex items-center gap-1 px-3 py-1 bg-teal-600 hover:bg-teal-500 rounded-full text-xs font-medium transition-colors">
                                    🔗 NCBI Protein ({data.protein_id})
                                </a>
                            )}
                            {data.ncbi_gene_url && (
                                <a href={data.ncbi_gene_url} target="_blank" rel="noopener noreferrer"
                                    className="inline-flex items-center gap-1 px-3 py-1 bg-emerald-600 hover:bg-emerald-500 rounded-full text-xs font-medium transition-colors">
                                    🔗 NCBI Gene ({data.gene_name})
                                </a>
                            )}
                            {data.db_xref?.map((ref, i) => (
                                <span key={i} className="px-3 py-1 bg-slate-700 rounded-full text-xs text-slate-300">
                                    {ref}
                                </span>
                            ))}
                        </div>
                    </div>

                    {/* Central Dogma Flow Diagram */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <h3 className="text-sm font-semibold text-slate-700 mb-4">Flujo del Dogma Central</h3>

                        {/* Visual Pipeline */}
                        <div className="flex items-center justify-center gap-2 mb-6 flex-wrap">
                            <div className="bg-blue-50 border-2 border-blue-300 rounded-xl px-4 py-3 text-center">
                                <div className="text-xs text-blue-500 font-medium">DNA</div>
                                <div className="text-sm font-bold text-blue-700">Doble Hebra</div>
                                <div className="text-xs text-blue-400 mt-1">5'→3' / 3'→5'</div>
                            </div>
                            <div className="text-center">
                                <div className="text-2xl">→</div>
                                <div className="text-xs text-slate-500 font-medium mt-1">Transcripción</div>
                                <div className="text-xs text-slate-400">RNA Polimerasa</div>
                            </div>
                            <div className="bg-purple-50 border-2 border-purple-300 rounded-xl px-4 py-3 text-center">
                                <div className="text-xs text-purple-500 font-medium">mRNA</div>
                                <div className="text-sm font-bold text-purple-700">Hebra Simple</div>
                                <div className="text-xs text-purple-400 mt-1">5'→3'</div>
                            </div>
                            <div className="text-center">
                                <div className="text-2xl">→</div>
                                <div className="text-xs text-slate-500 font-medium mt-1">Traducción</div>
                                <div className="text-xs text-slate-400">Ribosoma</div>
                            </div>
                            <div className="bg-emerald-50 border-2 border-emerald-300 rounded-xl px-4 py-3 text-center">
                                <div className="text-xs text-emerald-500 font-medium">Proteína</div>
                                <div className="text-sm font-bold text-emerald-700">Cadena Péptido</div>
                                <div className="text-xs text-emerald-400 mt-1">N→C terminal</div>
                            </div>
                        </div>

                        {/* Sequence Display: DNA double strand */}
                        <div className="space-y-4">
                            {/* DNA Sense Strand (5'→3') */}
                            <div>
                                <div className="flex items-center gap-2 mb-1">
                                    <span className="text-xs font-semibold text-blue-600 bg-blue-50 px-2 py-0.5 rounded">DNA Hebra Codificante (Sense)</span>
                                    <span className="text-xs text-slate-400">5' → 3'</span>
                                </div>
                                <div className="font-mono text-xs bg-slate-50 p-3 rounded-lg overflow-x-auto whitespace-nowrap border border-slate-100">
                                    <span className="text-blue-500 font-bold mr-1">5'</span>
                                    {renderColoredDNA(data.dna_sense_5to3)}
                                    <span className="text-blue-500 font-bold ml-1">3'</span>
                                    {data.length_bp > 300 && <span className="text-slate-400 ml-2">...({data.length_bp - 300} más)</span>}
                                </div>
                            </div>

                            {/* Base pairing indicators */}
                            <div className="flex justify-center">
                                <div className="text-xs text-slate-300 font-mono tracking-widest">
                                    {'| '.repeat(Math.min(30, Math.floor(data.dna_sense_5to3?.length / 10 || 0)))}
                                </div>
                            </div>

                            {/* DNA Template Strand (3'→5') */}
                            <div>
                                <div className="flex items-center gap-2 mb-1">
                                    <span className="text-xs font-semibold text-red-500 bg-red-50 px-2 py-0.5 rounded">DNA Hebra Molde (Template)</span>
                                    <span className="text-xs text-slate-400">3' → 5'</span>
                                </div>
                                <div className="font-mono text-xs bg-slate-50 p-3 rounded-lg overflow-x-auto whitespace-nowrap border border-slate-100">
                                    <span className="text-red-500 font-bold mr-1">3'</span>
                                    {renderColoredDNA(data.dna_template_3to5)}
                                    <span className="text-red-500 font-bold ml-1">5'</span>
                                </div>
                            </div>

                            {/* Transcription Arrow */}
                            <div className="flex items-center justify-center gap-2 py-2">
                                <div className="h-px flex-1 bg-gradient-to-r from-transparent via-purple-300 to-transparent"></div>
                                <span className="text-xs font-semibold text-purple-600 bg-purple-50 px-3 py-1 rounded-full border border-purple-200">
                                    ↓ TRANSCRIPCIÓN (RNA Pol lee template 3'→5', sintetiza mRNA 5'→3')
                                </span>
                                <div className="h-px flex-1 bg-gradient-to-r from-transparent via-purple-300 to-transparent"></div>
                            </div>

                            {/* mRNA (5'→3') */}
                            <div>
                                <div className="flex items-center gap-2 mb-1">
                                    <span className="text-xs font-semibold text-purple-600 bg-purple-50 px-2 py-0.5 rounded">mRNA</span>
                                    <span className="text-xs text-slate-400">5' → 3' (igual a hebra codificante, con U en vez de T)</span>
                                </div>
                                <div className="font-mono text-xs bg-purple-50 p-3 rounded-lg overflow-x-auto whitespace-nowrap border border-purple-100">
                                    <span className="text-purple-500 font-bold mr-1">5'</span>
                                    {renderColoredDNA(data.mrna_5to3, true)}
                                    <span className="text-purple-500 font-bold ml-1">3'</span>
                                </div>
                            </div>

                            {/* Translation Arrow */}
                            <div className="flex items-center justify-center gap-2 py-2">
                                <div className="h-px flex-1 bg-gradient-to-r from-transparent via-emerald-300 to-transparent"></div>
                                <span className="text-xs font-semibold text-emerald-600 bg-emerald-50 px-3 py-1 rounded-full border border-emerald-200">
                                    ↓ TRADUCCIÓN (Ribosoma lee mRNA 5'→3', sintetiza proteína N→C)
                                </span>
                                <div className="h-px flex-1 bg-gradient-to-r from-transparent via-emerald-300 to-transparent"></div>
                            </div>

                            {/* Protein */}
                            <div>
                                <div className="flex items-center gap-2 mb-1">
                                    <span className="text-xs font-semibold text-emerald-600 bg-emerald-50 px-2 py-0.5 rounded">Proteína</span>
                                    <span className="text-xs text-slate-400">N-terminal → C-terminal ({data.length_aa} aminoácidos)</span>
                                </div>
                                <div className="font-mono text-xs bg-emerald-50 p-3 rounded-lg overflow-x-auto whitespace-nowrap border border-emerald-100">
                                    <span className="text-emerald-600 font-bold mr-1">NH₂</span>
                                    {renderColoredProtein(showFullProtein ? data.full_protein : data.protein)}
                                    <span className="text-emerald-600 font-bold ml-1">COOH</span>
                                    {!showFullProtein && data.length_aa > 100 && (
                                        <button
                                            onClick={() => setShowFullProtein(true)}
                                            className="ml-2 text-emerald-500 hover:text-emerald-700 underline"
                                        >
                                            ...ver todo ({data.length_aa} aa)
                                        </button>
                                    )}
                                </div>
                            </div>
                        </div>
                    </div>

                    {/* Codon-by-codon alignment table */}
                    <div className="bg-white rounded-xl border border-slate-200 p-5">
                        <div className="flex justify-between items-center mb-4">
                            <h3 className="text-sm font-semibold text-slate-700">
                                Alineamiento Codón por Codón
                                <span className="text-slate-400 font-normal ml-2">
                                    ({Math.min(visibleCodons, data.codons?.length || 0)} de {data.total_codons} codones)
                                </span>
                            </h3>
                            <div className="flex gap-2">
                                {visibleCodons < (data.codons?.length || 0) && (
                                    <button
                                        onClick={() => setVisibleCodons(prev => Math.min(prev + 10, data.codons.length))}
                                        className="text-xs px-3 py-1 bg-teal-50 text-teal-600 rounded-lg hover:bg-teal-100"
                                    >
                                        Mostrar más
                                    </button>
                                )}
                            </div>
                        </div>

                        <div className="overflow-x-auto">
                            <table className="w-full text-center text-xs font-mono">
                                <thead>
                                    <tr className="border-b border-slate-200">
                                        <th className="py-2 px-1 text-slate-500 text-left font-medium">#</th>
                                        <th className="py-2 px-1 text-blue-600 font-medium">DNA 5'→3'</th>
                                        <th className="py-2 px-1 text-red-500 font-medium">DNA 3'→5'</th>
                                        <th className="py-2 px-1 text-purple-600 font-medium">mRNA</th>
                                        <th className="py-2 px-1 text-emerald-600 font-medium">AA</th>
                                        <th className="py-2 px-1 text-slate-600 font-medium">Nombre</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {data.codons?.slice(0, visibleCodons).map((codon, idx) => (
                                        <tr
                                            key={idx}
                                            className={`border-b border-slate-50 ${codon.is_start ? 'bg-green-50' :
                                                    codon.is_stop ? 'bg-red-50' :
                                                        idx % 2 === 0 ? 'bg-white' : 'bg-slate-25'
                                                }`}
                                        >
                                            <td className="py-1.5 px-1 text-slate-400 text-left">{idx + 1}</td>
                                            <td className="py-1.5 px-1 font-bold">
                                                {codon.dna_codon.split('').map((b, i) => (
                                                    <span key={i} style={{ color: BASE_COLORS[b] }}>{b}</span>
                                                ))}
                                            </td>
                                            <td className="py-1.5 px-1">
                                                {codon.template_codon.split('').map((b, i) => (
                                                    <span key={i} style={{ color: BASE_COLORS[b] }}>{b}</span>
                                                ))}
                                            </td>
                                            <td className="py-1.5 px-1 font-bold">
                                                {codon.rna_codon.split('').map((b, i) => (
                                                    <span key={i} style={{ color: BASE_COLORS[b] }}>{b}</span>
                                                ))}
                                            </td>
                                            <td className="py-1.5 px-1">
                                                <span
                                                    style={{
                                                        backgroundColor: AA_COLORS[codon.amino_acid] || '#888',
                                                        color: '#fff',
                                                        borderRadius: '3px',
                                                        padding: '1px 5px',
                                                        fontWeight: 700,
                                                    }}
                                                >
                                                    {codon.amino_acid}
                                                </span>
                                            </td>
                                            <td className="py-1.5 px-1 text-slate-500 text-left">
                                                {codon.is_start && <span className="text-green-600 font-bold">🟢 Met (INICIO)</span>}
                                                {codon.is_stop && <span className="text-red-600 font-bold">🔴 STOP</span>}
                                                {!codon.is_start && !codon.is_stop && codon.amino_acid_name}
                                            </td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        </div>
                    </div>

                    {/* Legend */}
                    <div className="bg-slate-50 rounded-xl border border-slate-200 p-4">
                        <h4 className="text-xs font-semibold text-slate-600 mb-3">Leyenda de Colores</h4>
                        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-xs">
                            <div>
                                <p className="font-semibold text-slate-700 mb-1">Bases de DNA/RNA</p>
                                <div className="flex gap-2 flex-wrap">
                                    {Object.entries(BASE_COLORS).map(([base, color]) => (
                                        <span key={base} style={{ color }} className="font-mono font-bold">
                                            {base}
                                        </span>
                                    ))}
                                </div>
                            </div>
                            <div>
                                <p className="font-semibold text-slate-700 mb-1">Hidrofóbicos</p>
                                <span className="inline-block w-3 h-3 rounded" style={{ backgroundColor: '#4a9eff' }}></span>
                                <span className="ml-1 text-slate-500">A,L,I,V,M,F,W,P</span>
                            </div>
                            <div>
                                <p className="font-semibold text-slate-700 mb-1">Polares / Cargados+</p>
                                <span className="inline-block w-3 h-3 rounded" style={{ backgroundColor: '#50c878' }}></span>
                                <span className="ml-1 text-slate-500">S,T,N,Q,Y,C,G</span>
                                <span className="ml-2 inline-block w-3 h-3 rounded" style={{ backgroundColor: '#ff6b6b' }}></span>
                                <span className="ml-1 text-slate-500">K,R,H</span>
                            </div>
                            <div>
                                <p className="font-semibold text-slate-700 mb-1">Cargados−</p>
                                <span className="inline-block w-3 h-3 rounded" style={{ backgroundColor: '#ffd93d' }}></span>
                                <span className="ml-1 text-slate-500">D,E</span>
                            </div>
                        </div>
                    </div>
                </div>
            )}

            {/* No data state */}
            {!data && !loading && !error && (
                <div className="text-center py-16 bg-white rounded-xl border border-slate-200">
                    <div className="w-20 h-20 mx-auto bg-teal-50 rounded-2xl flex items-center justify-center mb-4">
                        <span className="text-3xl">🧬</span>
                    </div>
                    <h3 className="text-lg font-bold text-slate-700 mb-2">Explorar el Dogma Central</h3>
                    <p className="text-slate-500 text-sm max-w-md mx-auto">
                        Seleccione un gen para visualizar el proceso completo de
                        <strong> transcripción</strong> (DNA → mRNA) y
                        <strong> traducción</strong> (mRNA → Proteína)
                        con las hebras 5'→3' y 3'→5'.
                    </p>
                </div>
            )}
        </div>
    )
}
