/**
 * FileManager Component
 * Displays detected genomic files from programmatic download
 */

const FILE_TYPE_COLORS = {
  'GenBank Full Flat File': 'bg-teal-100 text-teal-700',
  'GenBank': 'bg-teal-100 text-teal-700',
  'FASTA': 'bg-emerald-100 text-emerald-700',
  'GFF': 'bg-slate-100 text-slate-700',
  'Text File': 'bg-slate-100 text-slate-600',
}

function formatBytes(bytes) {
  if (bytes === 0) return '0 Bytes'
  const k = 1024
  const sizes = ['Bytes', 'KB', 'MB', 'GB']
  const i = Math.floor(Math.log(bytes) / Math.log(k))
  return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i]
}

export default function FileManager({ files, onRefresh, selectedGenomes }) {
  // Filtrar archivos solo de los genomas seleccionados/analizados
  const filteredFiles = selectedGenomes && selectedGenomes.length > 0
    ? files.filter(file => {
        const pathMatch = file.filepath.match(/(GCF_\d+\.\d+|GCA_\d+\.\d+)/)
        const accession = pathMatch ? pathMatch[1] : null
        return accession && selectedGenomes.includes(accession)
      })
    : files

  // Agrupar archivos por accession
  const filesByGenome = filteredFiles.reduce((acc, file) => {
    // Extraer accession del filepath
    const pathMatch = file.filepath.match(/(GCF_\d+\.\d+|GCA_\d+\.\d+)/)
    const accession = pathMatch ? pathMatch[1] : 'unknown'
    
    if (!acc[accession]) {
      acc[accession] = []
    }
    acc[accession].push(file)
    return acc
  }, {})

  const genomeCount = Object.keys(filesByGenome).length

  return (
    <div className="space-y-6">
      {/* Header Info */}
      <div className="bg-teal-50 border border-teal-100 rounded-xl p-5">
        <div className="flex items-start gap-4">
          <div className="w-12 h-12 bg-teal-500 rounded-xl flex items-center justify-center flex-shrink-0">
            <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
          </div>
          <div className="flex-1">
            <h3 className="font-semibold text-slate-800 mb-1">Archivos de {genomeCount} genoma{genomeCount > 1 ? 's' : ''}</h3>
            <p className="text-sm text-slate-600 mb-3">
              Archivos genómicos descargados desde NCBI Datasets API v2
            </p>
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-sm">
              <div className="bg-white rounded-lg p-2 border border-teal-100">
                <span className="text-slate-500">Formato:</span>
                <span className="ml-1 font-medium text-slate-700">GenBank</span>
              </div>
              <div className="bg-white rounded-lg p-2 border border-teal-100">
                <span className="text-slate-500">API:</span>
                <span className="ml-1 font-medium text-slate-700">NCBI v2</span>
              </div>
              <div className="bg-white rounded-lg p-2 border border-teal-100">
                <span className="text-slate-500">DB:</span>
                <span className="ml-1 font-medium text-slate-700">nucleotide</span>
              </div>
              <div className="bg-white rounded-lg p-2 border border-teal-100">
                <span className="text-slate-500">Base:</span>
                <span className="ml-1 font-medium text-slate-700">genome</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Actions */}
      <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4">
        <h2 className="text-xl font-bold text-slate-800">Archivos Detectados</h2>
        <button
          onClick={onRefresh}
          className="flex items-center gap-2 px-4 py-2 bg-slate-100 text-slate-700 rounded-lg hover:bg-slate-200 transition-colors text-sm font-medium"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
          </svg>
          Actualizar
        </button>
      </div>

      {/* Files by Genome */}
      {filteredFiles.length === 0 ? (
        <div className="text-center py-12 bg-white rounded-xl border border-slate-200">
          <div className="w-16 h-16 mx-auto bg-slate-100 rounded-xl flex items-center justify-center mb-4">
            <svg className="w-8 h-8 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
          </div>
          <p className="text-slate-500">No se encontraron archivos genómicos</p>
          <p className="text-sm text-slate-400 mt-1">
            {selectedGenomes && selectedGenomes.length > 0 
              ? 'Los genomas seleccionados no tienen archivos disponibles'
              : 'Analiza genomas primero para ver sus archivos aquí'}
          </p>
        </div>
      ) : (
        Object.entries(filesByGenome).map(([accession, genomeFiles]) => (
          <div key={accession} className="space-y-3">
            <h3 className="font-semibold text-slate-700 flex items-center gap-2">
              <span className="w-2 h-2 bg-teal-500 rounded-full"></span>
              {accession}
              <span className="text-sm text-slate-500 font-normal">({genomeFiles.length} archivos)</span>
            </h3>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
            {genomeFiles.map((file, index) => (
            <div
              key={index}
              className={`bg-white rounded-xl border p-4 transition-all hover:shadow-md ${file.is_primary ? 'border-teal-300 ring-1 ring-teal-200' : 'border-slate-200'
                }`}
            >
              <div className="flex items-start justify-between">
                <div className="flex-1 min-w-0">
                  <p className="font-medium text-slate-800 truncate text-sm" title={file.filename}>
                    {file.filename}
                  </p>
                  <p className="text-sm text-slate-500 mt-0.5">
                    {formatBytes(file.size_bytes)}
                  </p>
                </div>
                {file.is_primary && (
                  <span className="ml-2 px-2 py-0.5 text-xs font-medium bg-teal-100 text-teal-700 rounded-full">
                    Principal
                  </span>
                )}
              </div>
              <div className="mt-3">
                <span className={`inline-flex px-2 py-0.5 text-xs font-medium rounded ${FILE_TYPE_COLORS[file.file_type] || 'bg-slate-100 text-slate-600'
                  }`}>
                  {file.file_type}
                </span>
              </div>
            </div>
            ))}
            </div>
          </div>
        ))
      )}

      {/* Summary */}
      {filteredFiles.length > 0 && (
        <div className="bg-gradient-to-r from-teal-600 to-emerald-600 rounded-xl p-5 text-white">
          <h3 className="font-medium mb-3">Resumen de Archivos</h3>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-4">
            <div>
              <p className="text-2xl font-bold">{filteredFiles.length}</p>
              <p className="text-teal-100 text-sm">Total archivos</p>
            </div>
            <div>
              <p className="text-2xl font-bold">
                {filteredFiles.filter(f => f.extension === '.gbff').length}
              </p>
              <p className="text-teal-100 text-sm">GenBank</p>
            </div>
            <div>
              <p className="text-2xl font-bold">
                {filteredFiles.filter(f => f.extension === '.fna' || f.extension === '.fasta').length}
              </p>
              <p className="text-teal-100 text-sm">FASTA</p>
            </div>
            <div>
              <p className="text-2xl font-bold">
                {filteredFiles.filter(f => f.extension === '.gff').length}
              </p>
              <p className="text-teal-100 text-sm">GFF</p>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
