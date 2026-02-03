/**
 * AIValidation Component
 * Displays AI-powered scientific validation results with bioinformatics focus
 */

export default function AIValidation({ validationData, isValidating, onValidate, hasAnalysis }) {
  if (!hasAnalysis) {
    return (
      <div className="text-center py-20">
        <div className="max-w-xl mx-auto">
          <div className="w-20 h-20 mx-auto bg-slate-100 rounded-2xl flex items-center justify-center mb-6">
            <svg className="w-10 h-10 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
            </svg>
          </div>
          <h2 className="text-xl font-bold text-slate-800 mb-3">Validación Científica con IA</h2>
          <p className="text-slate-600">
            Primero ejecute el análisis del genoma para validar con IA.
          </p>
        </div>
      </div>
    )
  }

  if (isValidating) {
    return (
      <div className="text-center py-20">
        <div className="w-16 h-16 mx-auto border-4 border-teal-200 border-t-teal-600 rounded-full animate-spin mb-6"></div>
        <h2 className="text-xl font-bold text-slate-800 mb-2">Consultando Claude AI...</h2>
        <p className="text-slate-500">Validando resultados científicamente</p>
      </div>
    )
  }

  if (!validationData) {
    return (
      <div className="text-center py-20">
        <div className="max-w-xl mx-auto">
          <div className="w-20 h-20 mx-auto bg-teal-50 rounded-2xl flex items-center justify-center mb-6">
            <svg className="w-10 h-10 text-teal-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
            </svg>
          </div>
          <h2 className="text-xl font-bold text-slate-800 mb-3">Validación Científica con IA</h2>
          <p className="text-slate-600 mb-6">
            Use <strong>Claude AI (Anthropic)</strong> para validar los resultados contra conocimiento científico establecido de <em>E. coli</em> K-12 MG1655.
          </p>
          <button
            onClick={onValidate}
            className="px-6 py-3 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl font-semibold hover:from-teal-700 hover:to-emerald-700 transition-all shadow-lg"
          >
            Ejecutar Validación Científica
          </button>
        </div>
      </div>
    )
  }

  const { codon_validation, gene_validation, comprehensive_validation } = validationData

  const getConfidenceColor = (confidence) => {
    if (confidence >= 85) return 'bg-emerald-500'
    if (confidence >= 70) return 'bg-amber-500'
    return 'bg-red-500'
  }

  const getConfidenceBg = (confidence) => {
    if (confidence >= 85) return 'bg-emerald-50 border-emerald-200'
    if (confidence >= 70) return 'bg-amber-50 border-amber-200'
    return 'bg-red-50 border-red-200'
  }

  const ValidationCard = ({ title, subtitle, validation }) => (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-slate-50 to-slate-100 px-5 py-4 border-b border-slate-200">
        <div className="flex items-center justify-between flex-wrap gap-3">
          <div>
            <h3 className="font-semibold text-slate-800">{title}</h3>
            {subtitle && <p className="text-xs text-slate-500">{subtitle}</p>}
          </div>
          <div className="flex items-center gap-3">
            {/* Confidence Meter */}
            <div className="flex items-center gap-2">
              <div className="w-24 h-3 bg-slate-200 rounded-full overflow-hidden">
                <div
                  className={`h-full ${getConfidenceColor(validation.confidence)} transition-all`}
                  style={{ width: `${validation.confidence}%` }}
                />
              </div>
              <span className="text-sm font-bold text-slate-700">{validation.confidence}%</span>
            </div>
            {/* Status Badge */}
            <span className={`px-3 py-1 rounded-full text-xs font-semibold ${validation.is_valid
                ? 'bg-emerald-100 text-emerald-700'
                : 'bg-amber-100 text-amber-700'
              }`}>
              {validation.is_valid ? 'VÁLIDO' : 'REVISAR'}
            </span>
          </div>
        </div>
      </div>

      <div className="p-5 space-y-5">
        {/* Key Findings */}
        {validation.key_findings?.length > 0 && (
          <div className={`rounded-lg p-4 border ${getConfidenceBg(validation.confidence)}`}>
            <h4 className="text-xs uppercase tracking-wide text-slate-600 font-semibold mb-2">
              Hallazgos Clave
            </h4>
            <ul className="space-y-1">
              {validation.key_findings.map((finding, idx) => (
                <li key={idx} className="flex items-start gap-2 text-sm text-slate-700">
                  <span className="text-emerald-500 mt-1">•</span>
                  {finding}
                </li>
              ))}
            </ul>
          </div>
        )}

        {/* Scientific Context */}
        {validation.scientific_context && (
          <div>
            <h4 className="text-xs uppercase tracking-wide text-slate-500 font-medium mb-2">
              Contexto Científico
            </h4>
            <p className="text-sm text-slate-600 bg-slate-50 p-3 rounded-lg border border-slate-100">
              {validation.scientific_context}
            </p>
          </div>
        )}

        {/* Interpretation */}
        <div>
          <h4 className="text-xs uppercase tracking-wide text-slate-500 font-medium mb-2">
            Interpretación
          </h4>
          <p className="text-sm text-slate-700 leading-relaxed">
            {validation.interpretation}
          </p>
        </div>

        {/* Discrepancies */}
        {validation.discrepancies?.length > 0 && validation.discrepancies[0] !== '' && (
          <div>
            <h4 className="text-xs uppercase tracking-wide text-amber-600 font-medium mb-2">
              Discrepancias Detectadas
            </h4>
            <div className="space-y-2">
              {validation.discrepancies.map((disc, idx) => (
                <div key={idx} className="flex items-start gap-2 bg-amber-50 p-3 rounded-lg text-sm text-amber-800 border border-amber-100">
                  <svg className="w-4 h-4 text-amber-500 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
                  </svg>
                  <span>{disc}</span>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Recommendations */}
        {validation.recommendations?.length > 0 && (
          <div>
            <h4 className="text-xs uppercase tracking-wide text-teal-600 font-medium mb-2">
              Recomendaciones
            </h4>
            <div className="space-y-2">
              {validation.recommendations.map((rec, idx) => (
                <div key={idx} className="flex items-start gap-2 bg-teal-50 p-3 rounded-lg text-sm text-teal-800 border border-teal-100">
                  <svg className="w-4 h-4 text-teal-500 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
                  </svg>
                  <span>{rec}</span>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Footer */}
      <div className="px-5 py-3 bg-slate-50 border-t border-slate-100 text-xs text-slate-400">
        Validado: {new Date(validation.timestamp).toLocaleString('es-ES')}
      </div>
    </div>
  )

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center gap-4">
        <div>
          <h2 className="text-xl font-bold text-slate-800">Validación Científica</h2>
          <p className="text-slate-500 text-sm">Análisis por Claude 3.5 Haiku (Anthropic)</p>
        </div>
        <button
          onClick={onValidate}
          className="px-4 py-2 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-lg font-medium hover:from-teal-700 hover:to-emerald-700 transition-all text-sm flex items-center gap-2"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
          </svg>
          Revalidar
        </button>
      </div>

      {/* Comprehensive Validation - Main Card */}
      <ValidationCard
        title="Validación Global del Genoma"
        subtitle="Evaluación integral de E. coli K-12 MG1655"
        validation={comprehensive_validation}
      />

      {/* Specific Validations */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <ValidationCard
          title="Análisis de Codones"
          subtitle="ATG, TAA, TAG, TGA"
          validation={codon_validation}
        />
        <ValidationCard
          title="Análisis de Genes"
          subtitle="Anotación y estadísticas"
          validation={gene_validation}
        />
      </div>

      {/* Reference Info */}
      <div className="bg-slate-800 rounded-xl p-5 text-white">
        <h3 className="font-medium mb-3 flex items-center gap-2">
          <svg className="w-5 h-5 text-teal-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
          Referencia: NC_000913.3
        </h3>
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
          <div>
            <p className="text-slate-400">Organismo</p>
            <p className="font-medium">E. coli K-12 MG1655</p>
          </div>
          <div>
            <p className="text-slate-400">Genoma</p>
            <p className="font-medium">4,641,652 bp</p>
          </div>
          <div>
            <p className="text-slate-400">GC Content</p>
            <p className="font-medium">50.79%</p>
          </div>
          <div>
            <p className="text-slate-400">Genes</p>
            <p className="font-medium">4,651</p>
          </div>
        </div>
      </div>
    </div>
  )
}
