/**
 * FileManager Component
 * Displays and manages detected genomic files
 */
import { useState } from 'react'
import { DocumentIcon, ArrowPathIcon, FolderOpenIcon } from '@heroicons/react/24/outline'
import toast from 'react-hot-toast'
import { api } from '../services/api'

const FILE_TYPE_COLORS = {
  'GenBank Full Flat File': 'bg-emerald-100 text-emerald-800',
  'GenBank': 'bg-emerald-100 text-emerald-800',
  'FASTA Nucleotide': 'bg-blue-100 text-blue-800',
  'FASTA': 'bg-blue-100 text-blue-800',
  'FASTA Amino Acid': 'bg-purple-100 text-purple-800',
  'General Feature Format': 'bg-orange-100 text-orange-800',
  'Gene Transfer Format': 'bg-yellow-100 text-yellow-800',
  'JSON Lines': 'bg-gray-100 text-gray-800',
  'Tab-Separated Values': 'bg-pink-100 text-pink-800',
}

function formatBytes(bytes) {
  if (bytes === 0) return '0 Bytes'
  const k = 1024
  const sizes = ['Bytes', 'KB', 'MB', 'GB']
  const i = Math.floor(Math.log(bytes) / Math.log(k))
  return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i]
}

export default function FileManager({ files, onRefresh }) {
  const [isExtracting, setIsExtracting] = useState(false)
  const [selectedFile, setSelectedFile] = useState(null)
  const [fileInfo, setFileInfo] = useState(null)

  const handleExtract = async () => {
    setIsExtracting(true)
    try {
      await api.extractZip()
      toast.success('Archivos extraídos correctamente')
      onRefresh()
    } catch (error) {
      toast.error('Error: ' + (error.response?.data?.detail || error.message))
    } finally {
      setIsExtracting(false)
    }
  }

  const handleFileClick = async (file) => {
    setSelectedFile(file)
    try {
      const info = await api.getFileInfo(file.filename)
      setFileInfo(info)
    } catch (error) {
      console.error('Error getting file info:', error)
    }
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex justify-between items-center">
        <h2 className="text-2xl font-bold text-gray-800">Archivos Genómicos</h2>
        <div className="flex space-x-3">
          <button
            onClick={handleExtract}
            disabled={isExtracting}
            className="flex items-center px-4 py-2 bg-emerald-600 text-white rounded-lg hover:bg-emerald-700 transition-colors disabled:opacity-50"
          >
            <FolderOpenIcon className="h-5 w-5 mr-2" />
            {isExtracting ? 'Extrayendo...' : 'Extraer ZIP'}
          </button>
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
              Haga clic en "Extraer ZIP" o coloque los archivos en la carpeta ncbi_dataset
            </p>
          </div>
        ) : (
          files.map((file, index) => (
            <div
              key={index}
              onClick={() => handleFileClick(file)}
              className={`bg-white rounded-xl shadow-sm p-4 cursor-pointer transition-all hover:shadow-md ${
                file.is_primary ? 'ring-2 ring-emerald-500' : ''
              } ${selectedFile?.filename === file.filename ? 'bg-emerald-50' : ''}`}
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

      {/* File Details Panel */}
      {selectedFile && fileInfo && (
        <div className="bg-white rounded-xl shadow-sm p-6">
          <h3 className="text-lg font-semibold text-gray-800 mb-4">
            Detalles del Archivo
          </h3>
          <dl className="grid grid-cols-2 gap-4">
            <div>
              <dt className="text-sm font-medium text-gray-500">Nombre</dt>
              <dd className="mt-1 text-sm text-gray-900">{fileInfo.filename}</dd>
            </div>
            <div>
              <dt className="text-sm font-medium text-gray-500">Extensión</dt>
              <dd className="mt-1 text-sm text-gray-900">{fileInfo.extension}</dd>
            </div>
            <div>
              <dt className="text-sm font-medium text-gray-500">Tamaño</dt>
              <dd className="mt-1 text-sm text-gray-900">
                {formatBytes(fileInfo.size_bytes)} ({fileInfo.size_mb} MB)
              </dd>
            </div>
            <div>
              <dt className="text-sm font-medium text-gray-500">Tipo</dt>
              <dd className="mt-1 text-sm text-gray-900">{fileInfo.file_type}</dd>
            </div>
            <div className="col-span-2">
              <dt className="text-sm font-medium text-gray-500">Hash MD5</dt>
              <dd className="mt-1 text-xs text-gray-600 font-mono break-all">
                {fileInfo.hash}
              </dd>
            </div>
            <div className="col-span-2">
              <dt className="text-sm font-medium text-gray-500">Ruta</dt>
              <dd className="mt-1 text-xs text-gray-600 break-all">
                {fileInfo.filepath}
              </dd>
            </div>
          </dl>
        </div>
      )}

      {/* Summary */}
      <div className="bg-gradient-to-r from-emerald-500 to-cyan-500 rounded-xl p-6 text-white">
        <h3 className="font-semibold mb-2">Resumen de Archivos</h3>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
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
          <div>
            <p className="text-3xl font-bold">
              {files.filter(f => ['.fna', '.fasta', '.faa'].includes(f.extension)).length}
            </p>
            <p className="text-emerald-100">FASTA</p>
          </div>
          <div>
            <p className="text-3xl font-bold">
              {files.filter(f => ['.gff', '.gtf'].includes(f.extension)).length}
            </p>
            <p className="text-emerald-100">GFF/GTF</p>
          </div>
        </div>
      </div>
    </div>
  )
}
