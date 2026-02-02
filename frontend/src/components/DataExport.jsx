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
    name: 'JSON Completo',
    description: 'Todos los resultados del análisis en formato JSON estructurado',
    icon: DocumentTextIcon,
    color: 'bg-blue-500',
    action: () => api.exportJson()
  },
  {
    id: 'csv-genes',
    name: 'CSV - Genes',
    description: 'Lista completa de genes con sus propiedades en formato CSV',
    icon: TableCellsIcon,
    color: 'bg-emerald-500',
    action: () => api.exportCsv('genes')
  },
  {
    id: 'csv-codons',
    name: 'CSV - Codones',
    description: 'Resultados del análisis de codones en formato CSV',
    icon: TableCellsIcon,
    color: 'bg-purple-500',
    action: () => api.exportCsv('codons')
  },
  {
    id: 'csv-statistics',
    name: 'CSV - Estadísticas',
    description: 'Estadísticas generales del genoma en formato CSV',
    icon: TableCellsIcon,
    color: 'bg-orange-500',
    action: () => api.exportCsv('statistics')
  },
  {
    id: 'csv-validation',
    name: 'CSV - Validación',
    description: 'Resultados de validación contra valores de referencia',
    icon: TableCellsIcon,
    color: 'bg-cyan-500',
    action: () => api.exportCsv('validation')
  },
  {
    id: 'pdf',
    name: 'Informe PDF',
    description: 'Informe completo con tablas y gráficos en formato PDF',
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

      {/* Bulk Export */}
      <div className="bg-gradient-to-r from-gray-700 to-gray-800 rounded-xl p-6 text-white">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="font-semibold text-lg">Exportación Masiva</h3>
            <p className="text-gray-300 mt-1">
              Descargue todos los formatos disponibles de una vez
            </p>
          </div>
          <button
            onClick={() => {
              if (!hasData) {
                toast.error('No hay datos para exportar')
                return
              }
              EXPORT_OPTIONS.forEach(opt => opt.action())
              toast.success('Descargando todos los archivos...')
            }}
            disabled={!hasData}
            className="px-6 py-3 bg-white text-gray-800 rounded-lg font-semibold hover:bg-gray-100 transition-colors disabled:opacity-50 disabled:cursor-not-allowed flex items-center"
          >
            <DocumentArrowDownIcon className="h-5 w-5 mr-2" />
            Descargar Todo
          </button>
        </div>
      </div>

      {/* Format Information */}
      <div className="bg-white rounded-xl shadow-sm p-6">
        <h3 className="font-semibold text-gray-800 mb-4">Información de Formatos</h3>
        <div className="space-y-4">
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-blue-100 rounded-lg">
              <DocumentTextIcon className="h-5 w-5 text-blue-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">JSON</h4>
              <p className="text-sm text-gray-500">
                Formato estructurado ideal para procesamiento programático. Incluye todos los datos 
                del análisis con metadatos completos.
              </p>
            </div>
          </div>
          
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-emerald-100 rounded-lg">
              <TableCellsIcon className="h-5 w-5 text-emerald-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">CSV</h4>
              <p className="text-sm text-gray-500">
                Formato tabular compatible con Excel, Google Sheets y herramientas de análisis. 
                Perfecto para análisis estadísticos adicionales.
              </p>
            </div>
          </div>
          
          <div className="flex items-start space-x-3">
            <div className="p-2 bg-red-100 rounded-lg">
              <DocumentChartBarIcon className="h-5 w-5 text-red-600" />
            </div>
            <div>
              <h4 className="font-medium text-gray-800">PDF</h4>
              <p className="text-sm text-gray-500">
                Informe formateado con tablas listo para presentación o impresión. Incluye 
                las principales métricas y resultados de validación.
              </p>
            </div>
          </div>
        </div>
      </div>

      {/* Tips */}
      <div className="bg-blue-50 border border-blue-200 rounded-xl p-4">
        <div className="flex">
          <div className="flex-shrink-0">
            <svg className="h-5 w-5 text-blue-400" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="ml-3">
            <h4 className="text-sm font-medium text-blue-800">Consejo</h4>
            <p className="text-sm text-blue-700 mt-1">
              Para citar estos resultados en publicaciones, use el formato JSON que incluye 
              metadatos de fecha, versión del análisis y referencia del genoma (NC_000913.3).
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
