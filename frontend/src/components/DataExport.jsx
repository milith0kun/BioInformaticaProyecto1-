/**
 * DataExport Component — Clean Laboratory Edition
 * Handles high-fidelity export of analysis results in professional formats.
 */
import { useState } from 'react'
import toast from 'react-hot-toast'
import { api } from '../services/api'

const EXPORT_OPTIONS = [
  {
    id: 'json',
    name: 'Dataset Estructurado (JSON)',
    description: 'Análisis completo: codones, genes y protocolos de validación.',
    color: 'blue',
    icon: (
      <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 7v10c0 2.21 3.582 4 8 4s8-1.79 8-4V7M4 7c0 2.21 3.582 4 8 4s8-1.79 8-4M4 7c0-2.21 3.582-4 8-4s8 1.79 8 4m0 5c0 2.21-3.582 4-8 4s-8-1.79-8-4" />
      </svg>
    ),
    action: () => api.exportJson()
  },
  {
    id: 'pdf',
    name: 'Reporte Ejecutivo (PDF)',
    description: 'Documento técnico con gráficos, tablas y métricas de alta resolución.',
    color: 'indigo',
    icon: (
      <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
      </svg>
    ),
    action: () => api.exportPdf()
  },
  {
    id: 'csv-genes',
    name: 'Matriz de Genes (CSV)',
    description: 'Dataset tabular optimizado para Excel, R o Python.',
    color: 'slate',
    icon: (
      <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 10h18M3 14h18m-9-4v8m-7 0h14a2 2 0 002-2V8a2 2 0 00-2-2H5a2 2 0 00-2 2v8a2 2 0 002 2z" />
      </svg>
    ),
    action: () => api.exportCsv('genes')
  }
]

const colorClasses = {
  blue: 'bg-blue-600',
  indigo: 'bg-indigo-600',
  slate: 'bg-slate-900',
}

export default function DataExport({ hasData, currentGenome, selectedGenomes }) {
  const [isExporting, setIsExporting] = useState({})

  const handleExport = async (option) => {
    if (!hasData) {
      toast.error('No hay datos disponibles. Ejecute el análisis genómico primero.')
      return
    }

    setIsExporting(prev => ({ ...prev, [option.id]: true }))

    try {
      toast.loading(`Generando ${option.name}...`, { id: option.id })
      option.action()
      setTimeout(() => {
        toast.success(`${option.name} generado correctamente`, { id: option.id })
        setIsExporting(prev => ({ ...prev, [option.id]: false }))
      }, 1500)
    } catch (error) {
      toast.error(`Fallo en exportación: ${error.message}`, { id: option.id })
      setIsExporting(prev => ({ ...prev, [option.id]: false }))
    }
  }

  return (
    <div className="space-y-10 animate-in fade-in duration-1000">
      {/* Topology Header */}
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-8 bg-white p-10 rounded-[3rem] border-2 border-slate-100 shadow-sm relative overflow-hidden">
        <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/5 blur-[100px] -mr-32 -mt-32"></div>
        <div className="space-y-4 relative z-10">
          <h2 className="text-3xl font-black text-slate-900 tracking-tighter uppercase italic">
            Exportación <span className="text-blue-600">de Activos</span>
          </h2>
          <div className="flex flex-wrap items-center gap-4">
            {hasData ? (
              <>
                <span className="px-4 py-1.5 bg-blue-50 text-blue-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-blue-100">
                  {currentGenome?.accession || 'MG1655 CORE'}
                </span>
                <p className="text-[10px] font-bold text-slate-500 uppercase tracking-widest italic">
                  {selectedGenomes?.length > 1 ? `${selectedGenomes.length} Genomas en entorno` : 'Dataset Individual Sincronizado'}
                </p>
              </>
            ) : (
              <span className="px-4 py-1.5 bg-amber-50 text-amber-600 text-[9px] font-black uppercase tracking-widest rounded-full border border-amber-100 animate-pulse">
                Esperando resultados del motor...
              </span>
            )}
          </div>
        </div>
      </div>

      {/* Main Export Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
        {EXPORT_OPTIONS.map((option) => (
          <button
            key={option.id}
            onClick={() => handleExport(option)}
            disabled={!hasData || isExporting[option.id]}
            className={`group bg-white rounded-[2.5rem] border-2 p-8 text-left transition-all duration-500 relative overflow-hidden ${
              !hasData 
                ? 'border-slate-50 opacity-40 cursor-not-allowed' 
                : 'border-slate-100 hover:border-blue-200 hover:shadow-2xl hover:shadow-blue-900/5 active:scale-95'
            }`}
          >
            <div className="relative z-10 flex flex-col h-full space-y-6">
              <div className={`w-14 h-14 ${colorClasses[option.color]} rounded-2xl flex items-center justify-center text-white shadow-xl transition-transform duration-700 group-hover:rotate-12`}>
                {isExporting[option.id] ? (
                  <svg className="w-6 h-6 animate-spin" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                ) : option.icon}
              </div>
              <div>
                <h3 className="text-sm font-black text-slate-900 uppercase tracking-tighter mb-2 group-hover:text-blue-600 transition-colors">{option.name}</h3>
                <p className="text-[11px] text-slate-500 font-medium leading-relaxed uppercase tracking-tight">{option.description}</p>
              </div>
              <div className="pt-4 border-t border-slate-50 flex items-center justify-between">
                <span className="text-[9px] font-black text-slate-400 uppercase tracking-widest">{isExporting[option.id] ? 'Procesando...' : 'Descargar Archivo'}</span>
                <svg className={`w-4 h-4 text-blue-600 transition-all duration-500 ${isExporting[option.id] ? 'opacity-0' : 'opacity-0 group-hover:opacity-100 translate-x-[-10px] group-hover:translate-x-0'}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M14 5l7 7m0 0l-7 7m7-7H3" />
                </svg>
              </div>
            </div>
          </button>
        ))}
      </div>

      {/* Technical Protocols Panel */}
      <div className="grid grid-cols-1 xl:grid-cols-2 gap-8">
        <div className="bg-slate-900 rounded-[3rem] p-12 text-white shadow-2xl shadow-blue-900/20 relative overflow-hidden">
          <div className="absolute top-0 right-0 w-64 h-64 bg-blue-500/10 blur-[100px] -mr-32 -mt-32"></div>
          <div className="relative z-10 space-y-8">
            <div className="flex items-center gap-6">
              <div className="w-12 h-12 bg-white/10 rounded-2xl flex items-center justify-center border border-white/10">
                <svg className="w-6 h-6 text-blue-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
              </div>
              <h3 className="text-xl font-black uppercase italic tracking-tighter">Protocolos de Exportación</h3>
            </div>
            <div className="space-y-6">
              <div className="p-6 bg-white/5 rounded-2xl border border-white/10 hover:bg-white/10 transition-all cursor-default group">
                <p className="text-[10px] font-black text-blue-400 uppercase tracking-widest mb-2">Formato JSON Estructural</p>
                <p className="text-xs text-slate-400 leading-relaxed font-medium">Contenedor de datos orientado a objetos. Incluye jerarquía completa de genes, distribuciones de codones y vectores de validación técnica.</p>
              </div>
              <div className="p-6 bg-white/5 rounded-2xl border border-white/10 hover:bg-white/10 transition-all cursor-default group">
                <p className="text-[10px] font-black text-indigo-400 uppercase tracking-widest mb-2">Matriz Tabular (CSV)</p>
                <p className="text-xs text-slate-400 leading-relaxed font-medium">Datos crudos optimizados para software de análisis estadístico. Incluye coordenadas de inicio/fin, productos y valores de GC.</p>
              </div>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-[3rem] border-2 border-slate-100 p-12 shadow-sm flex flex-col justify-between">
          <div className="space-y-6">
            <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-[0.3em]">Resumen de Contenido</h3>
            <div className="space-y-4">
              {[
                { label: 'Metadata Genómica', desc: 'Identificación de cepa, longitud y parámetros globales.' },
                { label: 'Perfiles de Codones', desc: 'RSCU, Nc y distribución espacial de terminadores.' },
                { label: 'Catálogo de Genes', description: 'Inventario completo de CDS con anotaciones RefSeq.' },
                { label: 'Auditoría de Validación', desc: 'Resultados de comparativa contra rangos bacterianos.' }
              ].map((item, i) => (
                <div key={i} className="flex items-start gap-4 group">
                  <div className="w-1.5 h-1.5 rounded-full bg-blue-600 mt-1.5 group-hover:scale-150 transition-transform"></div>
                  <div>
                    <p className="text-[11px] font-black text-slate-900 uppercase tracking-tight leading-none mb-1">{item.label}</p>
                    <p className="text-[10px] text-slate-500 font-medium leading-tight">{item.desc || item.description}</p>
                  </div>
                </div>
              ))}
            </div>
          </div>
          <div className="mt-10 p-6 bg-slate-50 rounded-2xl border border-slate-100 relative group overflow-hidden">
            <div className="absolute inset-0 bg-blue-600 translate-x-[-100%] group-hover:translate-x-0 transition-transform duration-1000 opacity-5"></div>
            <p className="text-[10px] text-slate-600 leading-relaxed font-bold relative z-10">
              <span className="text-blue-600 mr-2">NOTA TÉCNICA:</span> Los archivos generados son instantáneas del análisis actual. Los datos fuente GenBank se conservan íntegros en el directorio del servidor.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}