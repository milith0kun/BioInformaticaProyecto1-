/**
 * FileManager Component — Clean Laboratory Edition
 */
import { useState } from 'react'
import { api } from '../services/api'
import toast from 'react-hot-toast'

export default function FileManager({ files, onRefresh, selectedGenomes }) {
  const [isExtracting, setIsDownloading] = useState(false)

  const handleDownload = async (accession) => {
    try {
      toast.loading(`Iniciando descarga: ${accession}...`, { id: 'download' })
      await api.downloadGenome({
        accession,
        include_gbff: true,
        include_gff: true,
        include_fasta: true
      })
      toast.success(`Descarga de ${accession} iniciada. Verifica el estado en la pestaña de selección.`, { id: 'download', duration: 5000 })
      if (onRefresh) onRefresh()
    } catch (error) {
      toast.error(`Error: ${error.message}`, { id: 'download' })
    }
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-700">
      {/* Header Module */}
      <div className="bg-slate-900 rounded-[3rem] p-10 text-white relative overflow-hidden shadow-2xl shadow-blue-900/20">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32"></div>
        <div className="relative z-10 flex flex-col md:flex-row md:items-center justify-between gap-8">
          <div>
            <h2 className="text-3xl font-black italic uppercase tracking-tighter mb-2 text-blue-400">Gestión de Activos</h2>
            <p className="text-slate-400 text-[10px] font-bold uppercase tracking-[0.2em]">Repositorio de secuencias y metadatos locales</p>
          </div>
          <button 
            onClick={onRefresh}
            className="px-8 py-3 bg-white/10 hover:bg-white/20 rounded-2xl text-[10px] font-black uppercase tracking-widest transition-all border border-white/10 active:scale-95 shadow-xl"
          >
            Sincronizar Directorio
          </button>
        </div>
      </div>

      {/* File System Grid */}
      <div className="grid grid-cols-1 gap-6">
        {files.length === 0 ? (
          <div className="bg-white rounded-[3rem] border-2 border-dashed border-slate-100 p-20 text-center space-y-6">
            <div className="w-20 h-20 bg-slate-50 rounded-3xl mx-auto flex items-center justify-center border border-slate-100 opacity-40">
              <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M9 13h6m-3-3v6m-9 1V7a2 2 0 012-2h6l2 2h6a2 2 0 012 2v8a2 2 0 01-2 2H5a2 2 0 01-2-2z" strokeWidth={2}/></svg>
            </div>
            <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">No se han detectado secuencias en el volumen local</p>
          </div>
        ) : (
          <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 overflow-hidden shadow-sm">
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead className="bg-slate-50/50 border-b border-slate-100">
                  <tr>
                    <th className="px-8 py-5 text-[10px] font-black text-slate-400 uppercase tracking-widest">Identificador / Archivo</th>
                    <th className="px-8 py-5 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Tamaño</th>
                    <th className="px-8 py-5 text-[10px] font-black text-slate-400 uppercase tracking-widest text-right">Protocolo</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-50">
                  {files.map((file, i) => (
                    <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                      <td className="px-8 py-6">
                        <div className="flex items-center gap-4">
                          <div className="w-10 h-10 bg-slate-50 rounded-2xl flex items-center justify-center border border-slate-100 group-hover:bg-blue-50 group-hover:border-blue-100 transition-all text-slate-400 group-hover:text-blue-500 shadow-sm">
                            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M7 21h10a2 2 0 002-2V9.414a1 1 0 00-.293-.707l-5.414-5.414A1 1 0 0012.586 3H7a2 2 0 00-2 2v14a2 2 0 002 2z" strokeWidth={2}/></svg>
                          </div>
                          <div>
                            <p className="font-black text-slate-900 text-xs uppercase tracking-tight">{file.name}</p>
                            <p className="text-[9px] text-slate-400 font-bold uppercase tracking-widest">{file.type || 'Binario RAW'}</p>
                          </div>
                        </div>
                      </td>
                      <td className="px-8 py-6 text-right font-mono text-xs font-black text-slate-600">
                        {file.size || 'N/A'}
                      </td>
                      <td className="px-8 py-6 text-right">
                        <button className="text-[10px] font-black text-blue-600 uppercase tracking-widest hover:underline decoration-2 underline-offset-4">
                          Inspeccionar
                        </button>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}