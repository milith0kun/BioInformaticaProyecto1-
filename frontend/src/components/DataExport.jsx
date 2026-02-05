/**
 * DataExport Component
 * Handles export of analysis results in various formats
 */
import { useState } from 'react'
import toast from 'react-hot-toast'
import { api } from '../services/api'

const EXPORT_OPTIONS = [
  {
    id: 'json',
    name: 'Análisis Completo (JSON)',
    description: 'Resultados de codones, genes y validación',
    color: 'teal',
    action: () => api.exportJson()
  },
  {
    id: 'csv-genes',
    name: 'Datos de Genes (CSV)',
    description: 'Tabla con locus_tag, posición, longitud',
    color: 'emerald',
    action: () => api.exportCsv('genes')
  },
  {
    id: 'csv-codons',
    name: 'Análisis de Codones (CSV)',
    description: 'ATG, TAA, TAG, TGA con conteos',
    color: 'slate',
    action: () => api.exportCsv('codons')
  },
  {
    id: 'csv-statistics',
    name: 'Estadísticas Genómicas (CSV)',
    description: 'Tamaño, GC%, densidad génica',
    color: 'teal',
    action: () => api.exportCsv('statistics')
  },
  {
    id: 'pdf',
    name: 'Informe Completo (PDF)',
    description: 'Documento con tablas y estadísticas',
    color: 'emerald',
    action: () => api.exportPdf()
  }
]

const colorClasses = {
  teal: 'bg-teal-500',
  emerald: 'bg-emerald-500',
  slate: 'bg-slate-600',
}

export default function DataExport({ hasData, comparisonData, currentGenome, selectedGenomes }) {
  const [isExporting, setIsExporting] = useState({})

  const handleExport = async (option) => {
    if (!hasData) {
      toast.error('No hay datos para exportar. Ejecute el análisis primero.')
      return
    }

    setIsExporting(prev => ({ ...prev, [option.id]: true }))

    try {
      toast.loading(`Exportando ${option.name}...`, { id: option.id })
      option.action()
      toast.success(`${option.name} descargado`, { id: option.id })
    } catch (error) {
      toast.error(`Error al exportar: ${error.message}`, { id: option.id })
    } finally {
      setIsExporting(prev => ({ ...prev, [option.id]: false }))
    }
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4">
        <div>
          <h2 className="text-xl font-bold text-slate-800">Exportar Datos</h2>
          <p className="text-slate-500 text-sm">
            {selectedGenomes && selectedGenomes.length > 1 
              ? `Análisis detallado: ${currentGenome?.accession || 'genoma activo'} • Comparación: ${selectedGenomes.length} genomas`
              : 'Descargue los resultados del análisis'}
          </p>
        </div>
        {!hasData && (
          <span className="px-3 py-1 bg-amber-100 text-amber-700 rounded-full text-sm font-medium">
            Sin datos - ejecute el análisis primero
          </span>
        )}
      </div>

      {/* Export Options */}
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
        {EXPORT_OPTIONS.map((option) => (
          <button
            key={option.id}
            onClick={() => handleExport(option)}
            disabled={!hasData || isExporting[option.id]}
            className="bg-white rounded-xl border border-slate-200 p-5 text-left transition-all hover:shadow-md hover:border-teal-200 disabled:opacity-50 disabled:cursor-not-allowed group"
          >
            <div className="flex items-start gap-4">
              <div className={`w-10 h-10 ${colorClasses[option.color]} rounded-xl flex items-center justify-center group-hover:scale-110 transition-transform flex-shrink-0`}>
                <svg className="w-5 h-5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                </svg>
              </div>
              <div className="flex-1 min-w-0">
                <h3 className="font-medium text-slate-800 text-sm">{option.name}</h3>
                <p className="text-xs text-slate-500 mt-0.5">{option.description}</p>
              </div>
              {isExporting[option.id] ? (
                <svg className="w-5 h-5 text-slate-400 animate-spin" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" fill="none" />
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                </svg>
              ) : (
                <svg className="w-5 h-5 text-slate-300 group-hover:text-teal-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                </svg>
              )}
            </div>
          </button>
        ))}
      </div>

      {/* Info Box */}
      <div className="bg-gradient-to-r from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
        <div className="flex items-start gap-4">
          <div className="w-10 h-10 bg-white/20 rounded-xl flex items-center justify-center flex-shrink-0">
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
          </div>
          <div>
            <h3 className="font-medium mb-2">Formatos de Exportación</h3>
            <ul className="text-sm text-teal-100 space-y-1">
              <li><strong>JSON:</strong> Todos los datos estructurados</li>
              <li><strong>CSV:</strong> Tablas para Excel, R o Python</li>
              <li><strong>PDF:</strong> Informe ejecutivo formateado</li>
            </ul>
          </div>
        </div>
      </div>

      {/* Format Details */}
      <div className="bg-white rounded-xl border border-slate-200 p-5">
        <h3 className="font-semibold text-slate-800 mb-4">¿Qué contiene cada formato?</h3>
        <div className="space-y-4">
          {[
            { icon: 'teal', title: 'JSON - Análisis Completo', desc: 'Metadata, codones, genes, estadísticas y validación.' },
            { icon: 'emerald', title: 'CSV - Datos Tabulares', desc: 'Archivos separados: genes, codones, estadísticas.' },
            { icon: 'slate', title: 'PDF - Informe Ejecutivo', desc: 'Resumen con tablas y top genes.' },
          ].map((item, i) => (
            <div key={i} className="flex items-start gap-3">
              <div className={`w-8 h-8 ${colorClasses[item.icon]} rounded-lg flex items-center justify-center flex-shrink-0`}>
                <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                </svg>
              </div>
              <div>
                <h4 className="font-medium text-slate-800 text-sm">{item.title}</h4>
                <p className="text-xs text-slate-500">{item.desc}</p>
              </div>
            </div>
          ))}
        </div>

        <div className="mt-5 p-4 bg-amber-50 border border-amber-100 rounded-lg">
          <p className="text-xs text-amber-800">
            <strong>Nota:</strong> Los archivos exportados contienen solo resultados del análisis.
            Los archivos GenBank originales están en <code className="bg-amber-100 px-1 rounded">ncbi_dataset/</code>.
          </p>
        </div>
      </div>
    </div>
  )
}
