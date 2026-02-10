/**
 * FileManager Component — Clean Laboratory Edition
 * Displays all genomic files grouped by accession with proper metadata.
 */
import { useState, useMemo } from 'react'
import { api } from '../services/api'
import toast from 'react-hot-toast'

// Format bytes to human readable sizes
function formatSize(bytes) {
  if (!bytes || bytes === 0) return '0 B'
  const units = ['B', 'KB', 'MB', 'GB']
  const i = Math.floor(Math.log(bytes) / Math.log(1024))
  return (bytes / Math.pow(1024, i)).toFixed(2) + ' ' + units[i]
}

// Extension to icon color mapping - Clean Laboratory Palette
const EXT_COLORS = {
  '.gbff': { bg: 'bg-blue-50', text: 'text-blue-600', border: 'border-blue-100', label: 'GBFF' },
  '.fna': { bg: 'bg-indigo-50', text: 'text-indigo-600', border: 'border-indigo-100', label: 'FNA' },
  '.gff': { bg: 'bg-slate-50', text: 'text-slate-600', border: 'border-slate-200', label: 'GFF' },
  '.faa': { bg: 'bg-blue-50', text: 'text-blue-500', border: 'border-blue-100', label: 'FAA' },
  '.jsonl': { bg: 'bg-slate-50', text: 'text-slate-500', border: 'border-slate-200', label: 'JSONL' },
  '.txt': { bg: 'bg-slate-50', text: 'text-slate-400', border: 'border-slate-200', label: 'TXT' },
}

function getExtStyle(ext) {
  return EXT_COLORS[ext] || { bg: 'bg-slate-50', text: 'text-slate-400', border: 'border-slate-200', label: ext?.replace('.', '').toUpperCase() || '?' }
}

export default function FileManager({ files, onRefresh, selectedGenomes }) {
  const [expandedGroup, setExpandedGroup] = useState(null)

  // Group files by accession
  const grouped = useMemo(() => {
    if (!files || files.length === 0) return {}
    const groups = {}
    for (const f of files) {
      const acc = f.accession || 'sin_accession'
      if (!groups[acc]) groups[acc] = []
      groups[acc].push(f)
    }
    // Sort each group: primary first, then by filename
    for (const acc of Object.keys(groups)) {
      groups[acc].sort((a, b) => {
        if (a.is_primary && !b.is_primary) return -1
        if (!a.is_primary && b.is_primary) return 1
        return (a.filename || '').localeCompare(b.filename || '')
      })
    }
    return groups
  }, [files])

  const accessions = Object.keys(grouped)

  // Summary counts
  const summary = useMemo(() => {
    if (!files) return { total: 0, genbank: 0, fasta: 0, gff: 0 }
    return {
      total: files.length,
      genbank: files.filter(f => f.extension === '.gbff').length,
      fasta: files.filter(f => f.extension === '.fna' || f.extension === '.faa').length,
      gff: files.filter(f => f.extension === '.gff' || f.extension === '.gtf').length,
    }
  }, [files])

  const toggleGroup = (acc) => {
    setExpandedGroup(prev => prev === acc ? null : acc)
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-700">

      {/* Header */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            Archivos de <span className="text-blue-600">{accessions.length} genomas</span>
          </h2>
          <p className="text-sm text-slate-600 font-medium leading-relaxed">
            Archivos genómicos descargados desde NCBI Datasets API v2
          </p>
          <div className="flex flex-wrap gap-3">
            {[
              { label: 'Formato', value: 'GenBank' },
              { label: 'API', value: 'NCBI v2' },
              { label: 'DB', value: 'nucleotide' },
              { label: 'Base', value: 'genome' },
            ].map(tag => (
              <span key={tag.label} className="inline-flex items-center gap-2 px-4 py-1.5 bg-slate-50 rounded-full border border-slate-100 text-[10px] font-black uppercase tracking-widest">
                <span className="text-slate-400">{tag.label}:</span>
                <span className="text-blue-600">{tag.value}</span>
              </span>
            ))}
          </div>
        </div>
        <div className="flex flex-col items-end gap-3 relative z-10">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Archivos Detectados</p>
          <button
            onClick={onRefresh}
            className="flex items-center gap-3 px-6 py-3 bg-slate-900 text-white rounded-2xl text-[10px] font-black uppercase tracking-widest shadow-xl shadow-slate-200 transition-all hover:bg-blue-600 hover:-translate-y-1 active:scale-95"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" /></svg>
            Actualizar
          </button>
        </div>
      </div>

      {/* Genome File Groups */}
      {accessions.length === 0 ? (
        <div className="bg-white rounded-[3rem] border-2 border-dashed border-slate-100 p-20 text-center space-y-6">
          <div className="w-20 h-20 bg-slate-50 rounded-3xl mx-auto flex items-center justify-center border border-slate-100 opacity-40">
            <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M9 13h6m-3-3v6m-9 1V7a2 2 0 012-2h6l2 2h6a2 2 0 012 2v8a2 2 0 01-2 2H5a2 2 0 01-2-2z" strokeWidth={2}/></svg>
          </div>
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">No se han detectado secuencias en el volumen local</p>
        </div>
      ) : (
        <div className="space-y-5">
          {accessions.map(acc => {
            const groupFiles = grouped[acc]
            const isExpanded = expandedGroup === acc || expandedGroup === null
            const totalSize = groupFiles.reduce((sum, f) => sum + (f.size_bytes || 0), 0)

            return (
              <div key={acc} className="bg-white rounded-[2.5rem] border-2 border-slate-100 overflow-hidden shadow-sm transition-all hover:shadow-md">
                {/* Accession Header */}
                <button
                  onClick={() => toggleGroup(acc)}
                  className="w-full flex items-center justify-between p-7 hover:bg-slate-50/50 transition-colors"
                >
                  <div className="flex items-center gap-5">
                    <div className="w-12 h-12 bg-blue-50 rounded-2xl flex items-center justify-center border border-blue-100 shadow-sm">
                      <span className="text-lg">🧬</span>
                    </div>
                    <div className="text-left">
                      <p className="text-sm font-black text-slate-900 uppercase tracking-tighter">{acc}</p>
                      <p className="text-[10px] font-bold text-slate-600 uppercase tracking-widest">
                        ({groupFiles.length} archivos)
                      </p>
                    </div>
                  </div>
                  <div className="flex items-center gap-6">
                    <span className="hidden sm:block text-[10px] font-black text-slate-400 uppercase tracking-widest">
                      {formatSize(totalSize)}
                    </span>
                    <svg className={`w-5 h-5 text-slate-400 transition-transform duration-300 ${isExpanded ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M19 9l-7 7-7-7" />
                    </svg>
                  </div>
                </button>

                {/* File List */}
                {isExpanded && (
                  <div className="border-t border-slate-50 divide-y divide-slate-50 animate-in fade-in slide-in-from-top-2 duration-300">
                    {groupFiles.map((file, idx) => {
                      const style = getExtStyle(file.extension)
                      return (
                        <div key={idx} className="flex items-center justify-between px-7 py-5 hover:bg-slate-50/40 transition-colors group">
                          <div className="flex items-center gap-5 flex-1 min-w-0">
                            <div className={`w-10 h-10 ${style.bg} rounded-xl flex items-center justify-center border ${style.border} shadow-sm flex-shrink-0`}>
                              <span className={`text-[9px] font-black ${style.text} uppercase`}>{style.label}</span>
                            </div>
                            <div className="flex-1 min-w-0">
                              <div className="flex items-center gap-3 flex-wrap">
                                <p className="font-mono text-xs font-black text-slate-900 truncate">{file.filename}</p>
                                {file.is_primary && (
                                  <span className="px-3 py-0.5 bg-blue-600 text-white text-[8px] font-black uppercase tracking-widest rounded-full shadow-sm flex-shrink-0">
                                    Principal
                                  </span>
                                )}
                              </div>
                              <p className="text-[10px] text-slate-500 font-bold uppercase tracking-widest mt-1">{file.file_type || 'Archivo'}</p>
                            </div>
                          </div>
                          <div className="flex items-center gap-6 flex-shrink-0 ml-4">
                            <p className="font-mono text-xs font-black text-slate-600 tabular-nums">
                              {formatSize(file.size_bytes)}
                            </p>
                          </div>
                        </div>
                      )
                    })}
                  </div>
                )}
              </div>
            )
          })}
        </div>
      )}

      {/* Summary Footer */}
      {accessions.length > 0 && (
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 shadow-sm p-8">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em] mb-5">Resumen de Archivos</p>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-6">
            <div className="text-center p-5 bg-slate-50 rounded-2xl border border-slate-100">
              <p className="text-3xl font-black text-slate-900 tracking-tighter">{summary.total}</p>
              <p className="text-[10px] font-black text-slate-500 uppercase tracking-widest mt-2">Total archivos</p>
            </div>
            <div className="text-center p-5 bg-blue-50/50 rounded-2xl border border-blue-100/50">
              <p className="text-3xl font-black text-blue-600 tracking-tighter">{summary.genbank}</p>
              <p className="text-[10px] font-black text-blue-500 uppercase tracking-widest mt-2">GenBank</p>
            </div>
            <div className="text-center p-5 bg-indigo-50/50 rounded-2xl border border-indigo-100/50">
              <p className="text-3xl font-black text-indigo-600 tracking-tighter">{summary.fasta}</p>
              <p className="text-[10px] font-black text-indigo-500 uppercase tracking-widest mt-2">FASTA</p>
            </div>
            <div className="text-center p-5 bg-slate-100 rounded-2xl border border-slate-200">
              <p className="text-3xl font-black text-slate-700 tracking-tighter">{summary.gff}</p>
              <p className="text-[10px] font-black text-slate-600 uppercase tracking-widest mt-2">GFF</p>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}