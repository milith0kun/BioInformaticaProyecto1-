/**
 * FileManager Component
 * Displays detected genomic files from programmatic download
 */
import { ArrowPathIcon } from '@heroicons/react/24/outline'
import { api } from '../services/api'

const FILE_TYPE_COLORS = {
  'GenBank Full Flat File': 'bg-emerald-100 text-emerald-800',
  'GenBank': 'bg-emerald-100 text-emerald-800',
  'Text File': 'bg-gray-100 text-gray-800',
}

function formatBytes(bytes) {
  if (bytes === 0) return '0 Bytes'
  const k = 1024
  const sizes = ['Bytes', 'KB', 'MB', 'GB']
  const i = Math.floor(Math.log(bytes) / Math.log(k))
  return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i]
}

export default function FileManager({ files, onRefresh }) {

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="bg-gradient-to-r from-emerald-50 to-cyan-50 border border-emerald-200 rounded-xl p-6 mb-6">
        <div className="flex items-start gap-4">
          <div className="flex-shrink-0">
            <svg className="w-12 h-12 text-emerald-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
          </div>
          <div className="flex-grow">
            <h3 className="text-xl font-bold text-gray-800 mb-2">✅ Descarga Programática con Bio.Entrez</h3>
            <p className="text-gray-700 leading-relaxed mb-3">
              El genoma de <strong>E. coli K-12 MG1655</strong> fue descargado automáticamente desde <strong>NCBI</strong> usando el módulo <code className="px-2 py-1 bg-white rounded border border-gray-300 text-sm font-mono">Bio.Entrez</code> de BioPython.
            </p>
            <div className="grid grid-cols-2 gap-4 text-sm">
              <div className="bg-white rounded-lg p-3 border border-gray-200">
                <span className="text-gray-600">Accession:</span>
                <strong className="ml-2 text-gray-800">NC_000913.3</strong>
              </div>
              <div className="bg-white rounded-lg p-3 border border-gray-200">
                <span className="text-gray-600">Método:</span>
                <strong className="ml-2 text-gray-800">Entrez.efetch()</strong>
              </div>
              <div className="bg-white rounded-lg p-3 border border-gray-200">
                <span className="text-gray-600">Base de datos:</span>
                <strong className="ml-2 text-gray-800">nucleotide</strong>
              </div>
              <div className="bg-white rounded-lg p-3 border border-gray-200">
                <span className="text-gray-600">Formato:</span>
                <strong className="ml-2 text-gray-800">GenBank</strong>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div className="flex justify-between items-center">
        <h2 className="text-2xl font-bold text-gray-800">Archivos Detectados</h2>
        <div className="flex space-x-3">
          <button
            onClick={onRefresh}
            className="flex items-center px-4 py-2 bg-gray-100 text-gray-700 rounded-lg hover:bg-gray-200 transition-colors"
          >
            <ArrowPathIcon className="h-5 w-5 mr-2" />
            Actualizar
          </button>
        </div>
      </div>

      {/* Files Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {files.length === 0 ? (
          <div className="col-span-full text-center py-12 bg-white rounded-xl shadow">
            <DocumentIcon className="h-16 w-16 mx-auto text-gray-300 mb-4" />
            <p className="text-gray-500">No se encontraron archivos genómicos</p>
            <p className="text-gray-400 text-sm mt-2">
              Los archivos se descargan automáticamente con Bio.Entrez
            </p>
          </div>
        ) : (
          files.map((file, index) => (
            <div
              key={index}
              className={`bg-white rounded-xl shadow-sm p-4 transition-all hover:shadow-md ${
                file.is_primary ? 'ring-2 ring-emerald-500' : ''
              }`}
            >
              <div className="flex items-start justify-between">
                <div className="flex-1 min-w-0">
                  <p className="font-medium text-gray-900 truncate" title={file.filename}>
                    {file.filename}
                  </p>
                  <p className="text-sm text-gray-500 mt-1">
                    {formatBytes(file.size_bytes)}
                  </p>
                </div>
                {file.is_primary && (
                  <span className="ml-2 px-2 py-1 text-xs font-medium bg-emerald-100 text-emerald-800 rounded-full">
                    Principal
                  </span>
                )}
              </div>
              <div className="mt-3">
                <span className={`inline-flex px-2 py-1 text-xs font-medium rounded-full ${
                  FILE_TYPE_COLORS[file.file_type] || 'bg-gray-100 text-gray-800'
                }`}>
                  {file.file_type}
                </span>
              </div>
            </div>
          ))
        )}
      </div>

      {/* Summary */}
      <div className="bg-gradient-to-r from-emerald-500 to-cyan-500 rounded-xl p-6 text-white">
        <h3 className="font-semibold mb-2">Resumen de Archivos</h3>
        <div className="grid grid-cols-2 md:grid-cols-2 gap-4">
          <div>
            <p className="text-3xl font-bold">{files.length}</p>
            <p className="text-emerald-100">Total archivos</p>
          </div>
          <div>
            <p className="text-3xl font-bold">
              {files.filter(f => f.extension === '.gbff').length}
            </p>
            <p className="text-emerald-100">GenBank</p>
          </div>
        </div>
      </div>
    </div>
  )
}
