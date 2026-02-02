/**
 * AIValidation Component
 * Displays AI-powered scientific validation results from Google Gemini
 */

export default function AIValidation({ validationData, isValidating, onValidate, hasAnalysis }) {
  if (!hasAnalysis) {
    return (
      <div className="text-center py-20">
        <div className="max-w-2xl mx-auto">
          <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-gradient-to-br from-purple-100 to-pink-100 mb-8">
            <svg className="w-12 h-12 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
            </svg>
          </div>
          <h2 className="text-3xl font-bold text-gray-800 mb-4">
            Validación con Inteligencia Artificial
          </h2>
          <p className="text-gray-600 text-lg leading-relaxed mb-6">
            Primero debe ejecutar el análisis completo del genoma antes de validar con IA.
          </p>
        </div>
      </div>
    )
  }

  if (isValidating) {
    return (
      <div className="text-center py-20">
        <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-gradient-to-br from-purple-100 to-pink-100 mb-8 animate-pulse">
          <svg className="animate-spin h-12 w-12 text-purple-600" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
          </svg>
        </div>
        <h2 className="text-2xl font-bold text-gray-800 mb-2">Consultando a Google Gemini...</h2>
        <p className="text-gray-600">La IA está analizando los resultados científicos</p>
      </div>
    )
  }

  if (!validationData) {
    return (
      <div className="text-center py-20">
        <div className="max-w-2xl mx-auto">
          <div className="inline-flex items-center justify-center w-24 h-24 rounded-full bg-gradient-to-br from-purple-100 to-pink-100 mb-8">
            <svg className="w-12 h-12 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
            </svg>
          </div>
          <h2 className="text-3xl font-bold text-gray-800 mb-4">
            Validación Científica con IA
          </h2>
          <p className="text-gray-600 text-lg leading-relaxed mb-8">
            Utilice Google Gemini para validar científicamente los resultados del análisis genómico.
            La IA comparará los datos con el conocimiento científico establecido sobre E. coli K-12.
          </p>
          <button
            onClick={onValidate}
            className="px-8 py-4 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-xl font-bold hover:shadow-xl hover:scale-105 transition-all duration-200"
          >
            🤖 Validar con IA
          </button>
        </div>
      </div>
    )
  }

  const { codon_validation, gene_validation, comprehensive_validation } = validationData

  const getConfidenceColor = (confidence) => {
    if (confidence >= 90) return 'from-green-500 to-emerald-500'
    if (confidence >= 70) return 'from-yellow-500 to-orange-500'
    return 'from-red-500 to-rose-500'
  }

  const ValidationCard = ({ title, validation, icon }) => (
    <div className="bg-white rounded-xl shadow-lg border border-gray-200 overflow-hidden">
      <div className="bg-gradient-to-r from-purple-50 to-pink-50 px-6 py-4 border-b border-purple-100">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="text-2xl">{icon}</div>
            <h3 className="text-xl font-bold text-gray-800">{title}</h3>
          </div>
          <div className="flex items-center gap-3">
            <div className="flex items-center gap-2">
              <span className="text-sm font-semibold text-gray-600">Confianza:</span>
              <div className="relative w-32 h-3 bg-gray-200 rounded-full overflow-hidden">
                <div 
                  className={`absolute h-full bg-gradient-to-r ${getConfidenceColor(validation.confidence)} transition-all duration-500`}
                  style={{ width: `${validation.confidence}%` }}
                ></div>
              </div>
              <span className="text-sm font-bold text-gray-800">{validation.confidence}%</span>
            </div>
            <div className={`px-4 py-1.5 rounded-full text-sm font-bold ${
              validation.is_valid 
                ? 'bg-green-100 text-green-700 border border-green-300'
                : 'bg-red-100 text-red-700 border border-red-300'
            }`}>
              {validation.is_valid ? '✓ VÁLIDO' : '✗ REVISAR'}
            </div>
          </div>
        </div>
      </div>

      <div className="p-6 space-y-6">
        {/* Interpretación */}
        <div>
          <h4 className="text-sm font-bold text-gray-700 mb-3 flex items-center gap-2">
            <svg className="w-5 h-5 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
            Interpretación Científica
          </h4>
          <p className="text-gray-700 leading-relaxed bg-gray-50 p-4 rounded-lg border border-gray-200">
            {validation.interpretation}
          </p>
        </div>

        {/* Discrepancias */}
        {validation.discrepancies && validation.discrepancies.length > 0 && (
          <div>
            <h4 className="text-sm font-bold text-gray-700 mb-3 flex items-center gap-2">
              <svg className="w-5 h-5 text-orange-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
              </svg>
              Discrepancias Detectadas
            </h4>
            <ul className="space-y-2">
              {validation.discrepancies.map((disc, idx) => (
                <li key={idx} className="flex items-start gap-3 bg-orange-50 p-3 rounded-lg border border-orange-200">
                  <span className="text-orange-600 font-bold mt-0.5">⚠</span>
                  <span className="text-gray-700">{disc}</span>
                </li>
              ))}
            </ul>
          </div>
        )}

        {/* Recomendaciones */}
        {validation.recommendations && validation.recommendations.length > 0 && (
          <div>
            <h4 className="text-sm font-bold text-gray-700 mb-3 flex items-center gap-2">
              <svg className="w-5 h-5 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
              </svg>
              Recomendaciones
            </h4>
            <ul className="space-y-2">
              {validation.recommendations.map((rec, idx) => (
                <li key={idx} className="flex items-start gap-3 bg-blue-50 p-3 rounded-lg border border-blue-200">
                  <span className="text-blue-600 font-bold mt-0.5">💡</span>
                  <span className="text-gray-700">{rec}</span>
                </li>
              ))}
            </ul>
          </div>
        )}

        {/* Timestamp */}
        <div className="pt-4 border-t border-gray-200">
          <p className="text-xs text-gray-500">
            Validado: {new Date(validation.timestamp).toLocaleString('es-ES')}
          </p>
        </div>
      </div>
    </div>
  )

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-3xl font-bold text-gray-800 mb-2">Validación con IA</h2>
          <p className="text-gray-600">Análisis científico realizado por Google Gemini</p>
        </div>
        <button
          onClick={onValidate}
          className="px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-xl font-bold hover:shadow-xl hover:scale-105 transition-all duration-200"
        >
          🔄 Revalidar
        </button>
      </div>

      {/* Comprehensive Validation */}
      <ValidationCard 
        title="Validación Comprehensiva"
        validation={comprehensive_validation}
        icon="🔬"
      />

      {/* Grid de validaciones específicas */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <ValidationCard 
          title="Análisis de Codones"
          validation={codon_validation}
          icon="🧬"
        />
        
        <ValidationCard 
          title="Análisis de Genes"
          validation={gene_validation}
          icon="📊"
        />
      </div>
    </div>
  )
}
