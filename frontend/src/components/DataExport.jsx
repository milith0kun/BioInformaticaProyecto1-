/**
 * DataExport Component
 * Handles export of analysis results in various formats
 */
import { useState } from 'react'
import { 
  DocumentArrowDownIcon,
  DocumentTextIcon,
  TableCellsIcon,
  DocumentChartBarIcon,
  ArrowDownTrayIcon
} from '@heroicons/react/24/outline'
import toast from 'react-hot-toast'
import { api } from '../services/api'

const EXPORT_OPTIONS = [
  {
    id: 'json',
    name: 'Análisis Completo (JSON)',
    description: 'Resultados de codones, genes y validación en JSON',
    icon: DocumentTextIcon,
    color: 'bg-blue-500',
    action: () => api.exportJson()
  },
  {
    id: 'csv-genes',
    name: 'Datos de Genes (CSV)',
    description: 'Tabla de genes con locus_tag, posición, longitud y producto',
    icon: TableCellsIcon,
    color: 'bg-emerald-500',
    action: () => api.exportCsv('genes')
  },
  {
    id: 'csv-codons',
    name: 'Análisis de Codones (CSV)',
    description: 'Codones ATG, TAA, TAG, TGA con conteos y densidad',
    icon: TableCellsIcon,
    color: 'bg-purple-500',
    action: () => api.exportCsv('codons')
  },
  {
    id: 'csv-statistics',
    name: 'Estadísticas Genómicas (CSV)',
    description: 'Tamaño genoma, GC%, densidad génica, estadísticas de tamaño',
    icon: TableCellsIcon,
    color: 'bg-orange-500',
    action: () => api.exportCsv('statistics')
  },
  {
    id: 'pdf',
    name: 'Informe Completo (PDF)',
    description: 'Documento con todas las tablas y estadísticas del análisis',
    icon: DocumentChartBarIcon,
    color: 'bg-red-500',
    action: () => api.exportPdf()
  }
]

export default function DataExport({ hasData }) {
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
      <div className="flex justify-between items-center">
        <div>
          <h2 className="text-2xl font-bold text-gray-800">Exportar Datos</h2>
          <p className="text-gray-500 mt-1">
            Descargue los resultados del análisis en diferentes formatos
          </p>
        </div>
        {!hasData && (
          <span className="px-3 py-1 bg-yellow-100 text-yellow-800 rounded-full text-sm font-medium">
            Sin datos - ejecute el análisis primero
          </span>
        )}
      </div>

      {/* Export Options Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {EXPORT_OPTIONS.map((option) => {
          const Icon = option.icon
          return (
            <button
              key={option.id}
              onClick={() => handleExport(option)}
              disabled={!hasData || isExporting[option.id]}
              className={`relative bg-white rounded-xl shadow-sm p-6 text-left transition-all hover:shadow-md disabled:opacity-50 disabled:cursor-not-allowed group`}
            >
              <div className="flex items-start space-x-4">
                <div className={`p-3 rounded-xl ${option.color} group-hover:scale-110 transition-transform`}>
                  <Icon className="h-6 w-6 text-white" />
                </div>
                <div className="flex-1">
                  <h3 className="font-semibold text-gray-800 group-hover:text-gray-900">
                    {option.name}
                  </h3>
                  <p className="text-sm text-gray-500 mt-1">
                    {option.description}
                  </p>
                </div>
              </div>
              
              {isExporting[option.id] ? (
                <div className="absolute top-4 right-4">
                  <svg className="animate-spin h-5 w-5 text-gray-400" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                  </svg>
                </div>
              ) : (
                <ArrowDownTrayIcon className="absolute top-4 right-4 h-5 w-5 text-gray-300 group-hover:text-gray-500 transition-colors" />
              )}
            </button>
          )
        })}
      </div>

      {/* Info Box */}
      <div className="bg-gradient-to-r from-blue-600 to-cyan-600 rounded-xl p-6 text-white">
        <div className="flex items-start space-x-4">
          <DocumentArrowDownIcon className="h-8 w-8 flex-shrink-0" />
          <div>
            <h3 className="font-semibold text-lg mb-2">Formatos de Exportación</h3>
            <p className="text-gray-300 mt-1">
              Descargue todos los formatos disponibles de una vez
            </p>
          </div>
          <button
            onClick={() => {
            <p className="text-blue-100 mt-1 text-sm">
              Los formatos disponibles contienen únicamente datos del análisis genómico realizado:
            </p>
            <ul className="mt-3 space-y-1 text-sm text-blue-50">
              <li>• <strong>JSON:</strong> Codones, genes, estadísticas y validación estructurados</li>
              <li>• <strong>CSV:</strong> Tablas de datos específicos para análisis en Excel/R/Python</li>
              <li>• <strong>PDF:</strong> Informe ejecutivo con resumen y tablas principales</li>
            </ul>
            <p className="text-blue-100 mt-3 text-sm">
              No se incluyen archivos del genoma original (.fasta, .gbff) - solo resultados del análisis.
            </p>
          </div>
        </div>
      </div>

      {/* Format Information */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <h3 className="font-semibold text-gray-800 mb-4">¿Qué contiene cada formato?</h3>
        <div className="space-y-4">
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-blue-100 rounded-lg">
              <DocumentTextIcon className="h-5 w-5 text-blue-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">JSON - Análisis Completo</h4>
              <p className="text-sm text-gray-500">
                Archivo estructurado con todos los resultados: metadata, análisis de codones (ATG, TAA, TAG, TGA), 
                lista completa de genes con propiedades, estadísticas genómicas y validación.
              </p>
            </div>
          </div>
          
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-emerald-100 rounded-lg">
              <TableCellsIcon className="h-5 w-5 text-emerald-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">CSV - Datos Tabulares</h4>
              <p className="text-sm text-gray-500">
                Archivos separados por tipo: genes (locus_tag, posición, producto), codones (conteos y densidad), 
                estadísticas (tamaño, GC%, densidad génica). Compatible con Excel, R, Python.
              </p>
            </div>
          </div>
          
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-red-100 rounded-lg">
              <DocumentChartBarIcon className="h-5 w-5 text-red-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">PDF - Informe Ejecutivo</h4>
              <p className="text-sm text-gray-500">
                Documento formateado con resumen ejecutivo, tablas de estadísticas principales, 
                top 10 genes más largos/cortos. Ideal para reportes y presentaciones.
              </p>
            </div>
          </div>
        </div>
        
        <div className="mt-6 p-4 bg-amber-50 border border-amber-200 rounded-lg">
          <p className="text-sm text-amber-800">
            <strong>Nota:</strong> Los archivos exportados contienen únicamente los resultados del análisis realizado 
            por el sistema (codones, genes, estadísticas). Los archivos GenBank (.gbff) originales descargados 
            de NCBI se encuentran en la carpeta <code className="bg-amber-100 px-1 rounded">ncbi_dataset/</code> del proyecto.
          </p>
        </div>
      </div>
    </div>
  )
}
