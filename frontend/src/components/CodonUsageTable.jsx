/**
 * CodonUsageTable Component — Clean Laboratory Edition
 * Comprehensive codon usage analysis with filtering and comparative metrics
 */
import { useState, useEffect, useMemo } from 'react'
import api from '../services/api'

const COLORS = ['#2563eb', '#4f46e5', '#7c3aed', '#db2777', '#0f172a']

export default function CodonUsageTable() {
  const [data, setData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [searchTerm, setSearchTerm] = useState('')
  const [activeTab, setActiveTab] = useState('table')
  const [aaFilter, setAaFilter] = useState('All')

  useEffect(() => {
    loadData()
  }, [])

  const loadData = async () => {
    setLoading(true)
    try {
      const result = await api.getCompleteCodonUsage()
      setData(result)
    } catch (e) {
      console.error(e)
    } finally {
      setLoading(false)
    }
  }

  // Get unique amino acids for filter
  const aminoAcids = useMemo(() => {
    if (!data?.codon_table) return []
    const aas = new Set(data.codon_table.map(row => row.amino_acid_name))
    return ['All', ...Array.from(aas).sort()]
  }, [data])

  const filteredTable = useMemo(() => {
    if (!data?.codon_table) return []
    let table = data.codon_table

    if (aaFilter !== 'All') {
      table = table.filter(row => row.amino_acid_name === aaFilter)
    }

    if (searchTerm) {
      const s = searchTerm.toLowerCase()
      table = table.filter(row => 
        row.codon.toLowerCase().includes(s) || 
        row.amino_acid.toLowerCase().includes(s) ||
        row.amino_acid_name.toLowerCase().includes(s)
      )
    }
    return table
  }, [data, searchTerm, aaFilter])

  if (loading && !data) {
    return (
      <div className="flex flex-col items-center justify-center py-48">
        <div className="w-20 h-20 border-4 border-slate-100 border-t-blue-600 rounded-full animate-spin mb-8 shadow-inner"></div>
        <p className="text-[10px] font-black text-slate-900 uppercase tracking-[0.4em] animate-pulse text-center">Indexando Codones...</p>
      </div>
    )
  }

  if (!data) return null

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Metrics Header */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        <div className="bg-slate-900 rounded-[2.5rem] p-8 text-white shadow-2xl relative overflow-hidden group shadow-blue-900/20">
          <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 blur-3xl -mr-16 -mt-16 rounded-full group-hover:scale-150 transition-transform duration-1000"></div>
          <div className="relative z-10 space-y-2">
            <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest">Total Codones</p>
            <p className="text-3xl font-black tracking-tighter">{(data.total_codons || 0).toLocaleString()}</p>
            <p className="text-[9px] font-bold text-slate-400 uppercase tracking-[0.2em]">Secuencias Procesadas</p>
          </div>
        </div>
        
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-2">GC3 Content</p>
          <p className="text-3xl font-black text-blue-600 tracking-tighter">{data.gc3_content || 0}%</p>
          <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight mt-1">3ra posición del codón</p>
          <div className="w-full h-1 bg-slate-100 rounded-full mt-4 overflow-hidden">
            <div className="h-full bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.4)]" style={{ width: `${data.gc3_content || 0}%` }}></div>
          </div>
        </div>

        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-2">Nc (Codones Efectivos)</p>
          <p className="text-3xl font-black text-indigo-600 tracking-tighter">{data.effective_number_of_codons || '0.00'}</p>
          <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight mt-1">Wright (1990)</p>
          <div className="w-full h-1 bg-slate-100 rounded-full mt-4 overflow-hidden">
            <div className="h-full bg-indigo-500 shadow-[0_0_8px_rgba(99,102,241,0.4)]" style={{ width: `${((data.effective_number_of_codons || 0)/61)*100}%` }}></div>
          </div>
        </div>

        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-8 shadow-sm">
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-2">Aminoácidos Codificados</p>
          <p className="text-3xl font-black text-slate-900 tracking-tighter">20</p>
          <p className="text-[9px] font-bold text-slate-400 uppercase tracking-tight mt-1">+ 3 stop codons</p>
        </div>
      </div>

      {/* Tabs and Controls */}
      <div className="bg-white rounded-[3rem] border-2 border-slate-100 overflow-hidden shadow-sm">
        <div className="p-8 border-b border-slate-50 bg-slate-50/50 flex flex-col md:flex-row md:items-center justify-between gap-6">
          <div className="flex bg-white p-1 rounded-2xl border border-slate-100 shadow-sm">
            {[
              { id: 'table', label: 'Tabla de Codones', icon: '📋' },
              { id: 'amino', label: 'Por Aminoácido', icon: '🧪' },
              { id: 'viz', label: 'Visualizaciones', icon: '📊' }
            ].map(tab => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`flex items-center gap-2 px-6 py-2 rounded-xl text-[10px] font-black uppercase tracking-widest transition-all ${
                  activeTab === tab.id ? 'bg-slate-900 text-white shadow-lg' : 'text-slate-400 hover:bg-slate-50'
                }`}
              >
                <span>{tab.icon}</span>
                <span>{tab.label}</span>
              </button>
            ))}
          </div>

          <div className="flex flex-wrap gap-4">
            <div className="space-y-1">
              <label className="text-[8px] font-black text-slate-400 uppercase tracking-[0.2em] px-2">Filtrar por aminoácido:</label>
              <select 
                value={aaFilter}
                onChange={(e) => setAaFilter(e.target.value)}
                className="block w-full px-4 py-2 text-xs font-bold bg-white border-2 border-slate-100 rounded-xl focus:outline-none focus:border-blue-500 transition-all"
              >
                {aminoAcids.map(aa => (
                  <option key={aa} value={aa}>{aa === 'All' ? 'Todos (64 codones)' : aa}</option>
                ))}
              </select>
            </div>
            
            <div className="space-y-1">
              <label className="text-[8px] font-black text-slate-400 uppercase tracking-[0.2em] px-2">Búsqueda rápida:</label>
              <div className="relative">
                <input 
                  type="text" 
                  value={searchTerm}
                  onChange={(e) => setSearchTerm(e.target.value)}
                  placeholder="Buscar Codón o ID..."
                  className="pl-10 pr-4 py-2 text-xs font-bold bg-white border-2 border-slate-100 rounded-xl focus:outline-none focus:border-blue-500 transition-all"
                />
                <svg className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeWidth={2.5} /></svg>
              </div>
            </div>
          </div>
        </div>
        
        {activeTab === 'table' && (
          <>
            <div className="overflow-x-auto max-h-[600px] overflow-y-auto custom-scrollbar">
              <table className="w-full text-left">
                <thead className="bg-white sticky top-0 z-10 border-b border-slate-100">
                  <tr>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Codón</th>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Aminoácido</th>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Conteo</th>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Frecuencia</th>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">RSCU</th>
                    <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">RSCU Barra</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-50">
                  {filteredTable.map((row, i) => (
                    <tr key={i} className="hover:bg-slate-50/50 transition-colors group">
                      <td className="px-8 py-5 font-mono text-sm font-black text-blue-600 tracking-widest">{row.codon}</td>
                      <td className="px-8 py-5">
                        <div className="flex flex-col">
                          <span className="text-[10px] font-black text-slate-900 uppercase tracking-tighter">{row.amino_acid_name}</span>
                          <span className="text-[10px] font-bold text-slate-400">({row.letter})</span>
                        </div>
                      </td>
                      <td className="px-8 py-5 text-right font-mono text-xs font-black text-slate-700">{(row.count || 0).toLocaleString()}</td>
                      <td className="px-8 py-5 text-right font-mono text-xs font-black text-slate-600">{row.frequency.toFixed(2)}</td>
                      <td className="px-8 py-5 text-right font-mono text-xs font-black text-indigo-600">{row.rscu.toFixed(3)}</td>
                      <td className="px-8 py-5">
                        <div className="flex items-center gap-4">
                          <div className="w-24 h-1.5 bg-slate-100 rounded-full overflow-hidden">
                            <div className={`h-full ${
                              row.rscu >= 1.5 ? 'bg-blue-600 shadow-[0_0_8px_rgba(37,99,235,0.4)]' : 
                              row.rscu >= 1.0 ? 'bg-indigo-400' : 
                              row.rscu >= 0.5 ? 'bg-slate-400' : 'bg-slate-300'
                            }`} style={{ width: `${Math.min((row.rscu/3)*100, 100)}%` }}></div>
                          </div>
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            {/* RSCU Legend Footer */}
            <div className="p-8 bg-slate-50 border-t border-slate-100">
              <div className="flex flex-wrap items-center gap-8">
                <span className="text-[10px] font-black text-slate-400 uppercase tracking-widest">Escala RSCU:</span>
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 rounded bg-blue-600 shadow-sm"></div>
                  <span className="text-[10px] font-bold text-slate-600 uppercase">≥1.5 (Preferido)</span>
                </div>
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 rounded bg-indigo-400"></div>
                  <span className="text-[10px] font-bold text-slate-600 uppercase">1.0–1.5 (Normal)</span>
                </div>
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 rounded bg-slate-400"></div>
                  <span className="text-[10px] font-bold text-slate-600 uppercase">0.5–1.0 (Menos usado)</span>
                </div>
                <div className="flex items-center gap-2">
                  <div className="w-3 h-3 rounded bg-slate-300"></div>
                  <span className="text-[10px] font-bold text-slate-600 uppercase">&lt;0.5 (Poco usado)</span>
                </div>
              </div>
            </div>
          </>
        )}

        {activeTab === 'amino' && (
          <div className="p-20 text-center space-y-4">
            <div className="w-16 h-16 bg-blue-50 rounded-full flex items-center justify-center mx-auto text-2xl">🧪</div>
            <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter">Vista Agrupada por Aminoácido</h4>
            <p className="text-xs text-slate-500 max-w-sm mx-auto leading-relaxed">Seleccione un aminoácido en el filtro superior para ver sus codones sinónimos y su adaptabilidad relativa.</p>
          </div>
        )}

        {activeTab === 'viz' && (
          <div className="p-20 text-center space-y-4">
            <div className="w-16 h-16 bg-indigo-50 rounded-full flex items-center justify-center mx-auto text-2xl">📊</div>
            <h4 className="text-lg font-black text-slate-900 uppercase tracking-tighter">Visualización Estadística</h4>
            <p className="text-xs text-slate-500 max-w-sm mx-auto leading-relaxed">Mapas de calor y gráficos de dispersión de uso de codones optimizados para análisis de sesgo genómico.</p>
          </div>
        )}
      </div>
    </div>
  )
}