/**
 * AIValidation Component — Clean Laboratory Edition
 * Performs AI-driven audit and shows technical validation against reference standards
 */
import { api } from '../services/api'
import toast from 'react-hot-toast'

const STATUS_COLORS = {
  PASS: 'bg-emerald-50 text-emerald-700 border-emerald-200',
  WARNING: 'bg-amber-50 text-amber-700 border-amber-200',
  FAIL: 'bg-red-50 text-red-700 border-red-200',
}

// Simple Markdown Parser (Duplicated for independence)
const MarkdownText = ({ content }) => {
  if (!content) return null
  const lines = content.split('\n')
  return (
    <div className="space-y-3 font-sans text-slate-600 leading-relaxed">
      {lines.map((line, i) => {
        if (line.startsWith('### ')) return <h4 key={i} className="text-sm font-black text-blue-600 mt-6 mb-2 uppercase tracking-wide">{line.replace('### ', '')}</h4>
        if (line.startsWith('## ')) return <h3 key={i} className="text-base font-black text-slate-800 mt-8 mb-4 border-b border-blue-100 pb-2">{line.replace('## ', '')}</h3>
        if (line.startsWith('# ')) return <h2 key={i} className="text-xl font-black text-slate-900 mt-8 mb-4">{line.replace('# ', '')}</h2>
        if (line.trim().startsWith('- ')) return <div key={i} className="flex gap-3 ml-4"><span className="text-blue-400 font-bold">•</span><p className="flex-1" dangerouslySetInnerHTML={{ __html: parseInline(line.replace('- ', '')) }} /></div>
        if (!line.trim()) return <div key={i} className="h-2"></div>
        return <p key={i} dangerouslySetInnerHTML={{ __html: parseInline(line) }} />
      })}
    </div>
  )
}

const parseInline = (text) => {
  return text
    .replace(/\*\*(.*?)\*\*/g, '<strong class="font-bold text-slate-800">$1</strong>')
    .replace(/\*(.*?)\*/g, '<em class="text-slate-500">$1</em>')
    .replace(/`([^`]+)`/g, '<code class="bg-slate-100 px-1.5 py-0.5 rounded text-xs font-mono text-pink-600 font-bold">$1</code>')
}

export default function AIValidation({
  validationData,
  technicalValidation,
  isValidating,
  onValidate,
  hasAnalysis
}) {
  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* AI Orchestrator Header */}
      <div className="bg-slate-900 rounded-[3rem] p-10 text-white relative overflow-hidden shadow-2xl shadow-blue-900/20 group">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32 rounded-full transition-transform group-hover:scale-110 duration-1000"></div>
        <div className="relative z-10 flex flex-col md:flex-row md:items-center justify-between gap-8">
          <div className="space-y-2">
            <h2 className="text-3xl font-black italic uppercase tracking-tighter text-blue-400">Verificación <span className="text-white">IA Antrópica</span></h2>
            <p className="text-slate-400 text-[10px] font-bold uppercase tracking-[0.2em]">Auditoría bioinformática automatizada</p>
          </div>
          <button
            onClick={onValidate}
            disabled={isValidating || !hasAnalysis}
            className="px-10 py-4 bg-blue-600 text-white rounded-[1.5rem] font-black uppercase tracking-widest text-[10px] hover:bg-blue-500 transition-all shadow-xl shadow-blue-900/40 disabled:opacity-50 active:scale-95 flex items-center gap-4 group/btn"
          >
            {isValidating ? (
              <>
                <div className="w-4 h-4 border-2 border-white/20 border-t-white rounded-full animate-spin"></div>
                Procesando Auditoría...
              </>
            ) : (
              <>
                <span>Iniciar Validación</span>
                <span className="text-lg group-hover/btn:translate-x-1 transition-transform">🛡️</span>
              </>
            )}
          </button>
        </div>
      </div>

      {/* Technical Validation Table */}
      {technicalValidation && (
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 overflow-hidden shadow-sm animate-in slide-in-from-bottom-4 duration-500">
          <div className="p-8 border-b border-slate-50 bg-slate-50/50 flex items-center justify-between">
            <div>
              <h3 className="text-[10px] font-black text-slate-900 uppercase tracking-[0.3em]">Métricas de Referencia</h3>
              <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest mt-1">Comparativa vs Rangos Bacterianos</p>
            </div>
            <div className={`px-4 py-1.5 rounded-full text-[9px] font-black uppercase tracking-widest border ${STATUS_COLORS[technicalValidation.overall_status]}`}>
              {technicalValidation.overall_status === 'PASS' ? '✅ NORMAL' : '⚠️ ATENCIÓN'}
            </div>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-left">
              <thead>
                <tr className="bg-white border-b border-slate-100">
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest">Métrica</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Valor</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Referencia</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-right">Desviación</th>
                  <th className="px-8 py-5 text-[9px] font-black text-slate-400 uppercase tracking-widest text-center">Estado</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-50">
                {(technicalValidation.items || []).map((item, i) => (
                  <tr key={i} className="hover:bg-slate-50/50 transition-colors">
                    <td className="px-8 py-5">
                      <p className="text-[10px] font-black text-slate-900 uppercase tracking-tight">{item.metric.replace(/_/g, ' ')}</p>
                    </td>
                    <td className="px-8 py-5 text-right">
                      <p className="font-mono text-xs font-black text-slate-700">{typeof item.calculated === 'number' ? item.calculated.toLocaleString() : item.calculated}</p>
                    </td>
                    <td className="px-8 py-5 text-right">
                      <p className="font-mono text-xs font-bold text-slate-400">{typeof item.reference === 'number' ? item.reference.toLocaleString() : item.reference}</p>
                    </td>
                    <td className="px-8 py-5 text-right">
                      <p className={`font-mono text-xs font-black ${Math.abs(item.deviation_percent) > 20 ? 'text-rose-500' : 'text-slate-400'}`}>
                        {item.deviation_percent.toFixed(1)}%
                      </p>
                    </td>
                    <td className="px-8 py-5">
                      <div className="flex justify-center">
                        <span className={`w-2 h-2 rounded-full ${item.status === 'PASS' ? 'bg-emerald-500' : item.status === 'WARNING' ? 'bg-amber-500' : 'bg-red-500'}`}></span>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {!validationData && !isValidating && (
        <div className="bg-white rounded-[3rem] border-2 border-dashed border-slate-100 p-20 text-center space-y-6 opacity-60">
          <div className="w-20 h-20 bg-slate-50 rounded-3xl mx-auto flex items-center justify-center border border-slate-100">
            <span className="text-3xl grayscale">🛡️</span>
          </div>
          <p className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">Esperando inicialización del motor de auditoría</p>
        </div>
      )}

      {validationData && (
        <div className="bg-white rounded-[2.5rem] border-2 border-slate-100 p-10 shadow-sm animate-in slide-in-from-bottom-8 duration-700">
          <div className="flex items-center gap-4 mb-10 border-b border-slate-50 pb-8">
            <div className="w-12 h-12 bg-blue-50 rounded-2xl flex items-center justify-center text-blue-600 border border-blue-100 shadow-sm">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" strokeWidth={2.5} /></svg>
            </div>
            <div>
              <h3 className="text-[10px] font-black text-slate-900 uppercase tracking-[0.3em]">Reporte de Auditoría IA</h3>
              <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest mt-1">Generado por GenomicAI v2.0</p>
            </div>
          </div>

          <div className="prose prose-slate max-w-none">
            <div className="bg-slate-50 rounded-[2rem] p-10 border border-slate-100">
              <MarkdownText content={
                validationData.comprehensive_validation?.interpretation ||
                validationData.validation ||
                (typeof validationData === 'string' ? validationData : JSON.stringify(validationData, null, 2))
              } />
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
