/**
 * InteractiveChat Component — Clean Laboratory Edition
 * Expert Bioinformatics Assistant with NCBI Integration
 */
import { useState, useEffect, useRef } from 'react'
import { api } from '../services/api'
import toast from 'react-hot-toast'

// Simple Markdown Parser Component
const MarkdownText = ({ content }) => {
  if (!content) return null

  // Split by newlines to handle blocks
  const lines = content.split('\n')

  return (
    <div className="space-y-2">
      {lines.map((line, i) => {
        // Headers
        if (line.startsWith('### ')) return <h4 key={i} className="text-sm font-black text-blue-600 mt-4 mb-2 uppercase tracking-wide">{line.replace('### ', '')}</h4>
        if (line.startsWith('## ')) return <h3 key={i} className="text-base font-black text-slate-800 mt-6 mb-3 border-b border-blue-100 pb-1">{line.replace('## ', '')}</h3>
        if (line.startsWith('# ')) return <h2 key={i} className="text-lg font-black text-slate-900 mt-6 mb-4">{line.replace('# ', '')}</h2>

        // Lists
        if (line.trim().startsWith('- ')) {
          return (
            <div key={i} className="flex gap-3 ml-2">
              <span className="text-blue-400 font-bold">•</span>
              <p className="flex-1" dangerouslySetInnerHTML={{ __html: parseInline(line.replace('- ', '')) }} />
            </div>
          )
        }

        // Empty lines
        if (!line.trim()) return <div key={i} className="h-2"></div>

        // Regular paragraph with inline parsing
        return <p key={i} dangerouslySetInnerHTML={{ __html: parseInline(line) }} />
      })}
    </div>
  )
}

// Helper to parse bold, italic, links
const parseInline = (text) => {
  let parsed = text
    .replace(/\*\*(.*?)\*\*/g, '<strong class="font-black text-slate-900">$1</strong>')
    .replace(/\*(.*?)\*/g, '<em class="text-slate-600">$1</em>')
    .replace(/\[(.*?)\]\((.*?)\)/g, '<a href="$2" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline font-bold transition-colors">$1</a>')
    .replace(/`([^`]+)`/g, '<code class="bg-slate-100 px-1.5 py-0.5 rounded text-xs font-mono text-pink-600 font-bold">$1</code>')
  return parsed
}

export default function InteractiveChat({ hasAnalysis, currentGenome }) {
  const [messages, setMessages] = useState([])
  const [input, setInput] = useState('')
  const [isLoading, setIsLoading] = useState(false)
  const [suggestions, setSuggestions] = useState([])
  const messagesEndRef = useRef(null)

  useEffect(() => {
    loadSuggestions()
    // Welcome message
    if (messages.length === 0) {
      setMessages([{
        role: 'assistant',
        content: `Hola. Soy **GenomicAI**, tu asistente bioinformático experto.\n\nEstoy conectado a **NCBI GenBank** y tengo acceso completo al análisis del genoma activo (*${currentGenome?.accession || 'Ninguno seleccionado'}*).\n\nPuedes preguntarme sobre:\n- Interpretación biológica de los datos\n- Funciones de genes específicos\n- Comparaciones evolutivas\n- Referencias bibliográficas`
      }])
    }
  }, [currentGenome])

  const loadSuggestions = async () => {
    try {
      const data = await api.getChatSuggestions()
      if (data && data.suggestions) {
        setSuggestions(data.suggestions)
      }
    } catch (e) {
      console.error("Error loading suggestions", e)
    }
  }

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' })
  }

  useEffect(() => {
    scrollToBottom()
  }, [messages])

  const handleSend = async (e, textOverride = null) => {
    if (e) e.preventDefault()
    const textToSend = textOverride || input

    if (!textToSend.trim() || isLoading) return

    const userMessage = { role: 'user', content: textToSend }
    setMessages(prev => [...prev, userMessage])
    setInput('')
    setIsLoading(true)

    try {
      // Extract accession from currentGenome object
      const genomeAccession = currentGenome?.accession || null
      const response = await api.sendChatMessage(textToSend, 'default', true, genomeAccession)

      const aiMessage = {
        role: 'assistant',
        content: response.response,
        references: response.references,
        ncbi_links: response.ncbi_links
      }
      setMessages(prev => [...prev, aiMessage])
    } catch (error) {
      console.error(error)
      toast.error('Error de conexión con el servidor IA')
      setMessages(prev => [...prev, { role: 'assistant', content: '**Error de sistema**: No pude procesar tu solicitud. Verifica tu conexión o clave API.' }])
    } finally {
      setIsLoading(false)
    }
  }

  return (
    <div className="flex flex-col h-[800px] animate-in fade-in duration-1000 overflow-hidden rounded-[3rem] border-2 border-slate-100 shadow-sm bg-slate-50 relative">
      {/* Dynamic Background */}
      <div className="absolute inset-0 pointer-events-none overflow-hidden">
        <div className="absolute -top-20 -right-20 w-96 h-96 bg-blue-500/5 blur-[100px] rounded-full"></div>
        <div className="absolute top-40 -left-20 w-72 h-72 bg-purple-500/5 blur-[80px] rounded-full"></div>
      </div>

      {/* Header */}
      <div className="bg-white/80 backdrop-blur-md p-6 border-b border-slate-100 flex items-center justify-between z-10 sticky top-0">
        <div className="flex items-center gap-5">
          <div className="relative">
            <div className="w-12 h-12 bg-slate-900 rounded-2xl flex items-center justify-center shadow-lg shadow-slate-200">
              <span className="text-2xl">🧬</span>
            </div>
            <div className="absolute -bottom-1 -right-1 w-4 h-4 bg-emerald-500 border-2 border-white rounded-full animate-pulse shadow-sm" title="Sistema Online"></div>
          </div>
          <div>
            <h3 className="text-sm font-black uppercase tracking-widest text-slate-900">Dr. GenomicAI</h3>
            <div className="flex items-center gap-2 mt-0.5">
              <span className="text-[9px] font-bold text-blue-600 bg-blue-50 px-2 py-0.5 rounded-full uppercase tracking-wider border border-blue-100">
                {currentGenome ? `Contexto: ${currentGenome.accession || currentGenome}` : 'Esperando Genoma'}
              </span>
            </div>
          </div>
        </div>
        <div className="flex items-center gap-3">
          <span className="px-3 py-1 bg-slate-100 rounded-lg text-[9px] font-black text-slate-400 uppercase tracking-widest border border-slate-200">
            NCBI Connected
          </span>
        </div>
      </div>

      {/* Messages Area */}
      <div className="flex-1 overflow-y-auto p-8 space-y-8 custom-scrollbar relative z-0">
        {messages.map((msg, i) => (
          <div key={i} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'} animate-in fade-in slide-in-from-bottom-4 duration-500`}>

            {msg.role === 'assistant' && (
              <div className="w-8 h-8 rounded-full bg-blue-100 mr-4 flex-shrink-0 flex items-center justify-center border border-blue-200 mt-1">
                <span className="text-xs">🤖</span>
              </div>
            )}

            <div className={`max-w-[85%] space-y-3 ${msg.role === 'user' ? 'items-end flex flex-col' : ''}`}>
              <div className={`p-6 rounded-[2rem] text-sm leading-relaxed shadow-sm ${msg.role === 'user'
                  ? 'bg-slate-900 text-white rounded-tr-none shadow-xl shadow-slate-200'
                  : 'bg-white text-slate-600 border border-slate-100 rounded-tl-none'
                }`}>
                {msg.role === 'user' ? (
                  <p>{msg.content}</p>
                ) : (
                  <MarkdownText content={msg.content} />
                )}
              </div>

              {/* Citations & Links */}
              {msg.role === 'assistant' && (msg.references?.length > 0 || msg.ncbi_links?.length > 0) && (
                <div className="bg-slate-50 border border-slate-100 rounded-2xl p-4 text-xs space-y-3 ml-2 max-w-2xl animate-in fade-in duration-700">
                  {msg.references?.length > 0 && (
                    <div>
                      <h5 className="font-black text-[9px] text-slate-400 uppercase tracking-widest mb-2 flex items-center gap-2">
                        <span className="w-1 h-1 bg-blue-500 rounded-full"></span> Referencias PubMed
                      </h5>
                      <ul className="space-y-1.5">
                        {msg.references.map((ref, idx) => (
                          <li key={idx} className="bg-white p-2 rounded-lg border border-slate-100 hover:border-blue-200 transition-colors">
                            <a href={ref.url} target="_blank" rel="noopener noreferrer" className="flex items-start gap-2 group">
                              <span className="text-blue-500 font-bold group-hover:underline">[{idx + 1}]</span>
                              <span className="flex-1 text-slate-600 font-medium group-hover:text-slate-800 transition-colors">
                                {ref.title} <span className="text-slate-400 text-[10px] ml-1">({ref.year || 'N/A'})</span>
                              </span>
                            </a>
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  {msg.ncbi_links?.length > 0 && (
                    <div>
                      <h5 className="font-black text-[9px] text-slate-400 uppercase tracking-widest mb-2 flex items-center gap-2">
                        <span className="w-1 h-1 bg-emerald-500 rounded-full"></span> Enlaces NCBI
                      </h5>
                      <div className="flex flex-wrap gap-2">
                        {msg.ncbi_links.map((link, idx) => (
                          <a key={idx} href={link.url} target="_blank" rel="noopener noreferrer"
                            className="px-3 py-1.5 bg-white border border-slate-200 rounded-lg hover:border-emerald-400 hover:text-emerald-700 hover:shadow-sm transition-all text-[10px] font-bold text-slate-500 uppercase flex items-center gap-2">
                            <span>{link.type === 'gene' ? '🧬' : '📄'}</span>
                            {link.name}
                          </a>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </div>

            {msg.role === 'user' && (
              <div className="w-8 h-8 rounded-full bg-slate-200 ml-4 flex-shrink-0 flex items-center justify-center border border-slate-300 mt-1">
                <span className="text-xs font-black text-slate-500">YOU</span>
              </div>
            )}
          </div>
        ))}

        {isLoading && (
          <div className="flex justify-start animate-pulse padding-4">
            <div className="bg-white border border-slate-100 px-6 py-4 rounded-[2rem] rounded-tl-none flex items-center gap-3 shadow-sm ml-12">
              <span className="text-xs font-black text-slate-400 uppercase tracking-widest">Analizando</span>
              <div className="flex gap-1.5">
                <div className="w-1.5 h-1.5 bg-blue-500 rounded-full animate-bounce"></div>
                <div className="w-1.5 h-1.5 bg-blue-500 rounded-full animate-bounce delay-100"></div>
                <div className="w-1.5 h-1.5 bg-blue-500 rounded-full animate-bounce delay-200"></div>
              </div>
            </div>
          </div>
        )}
        <div ref={messagesEndRef} />
      </div>

      {/* Footer & Input */}
      <div className="p-6 bg-white border-t border-slate-100 z-10 relative">
        {/* Suggestions */}
        {messages.length < 3 && suggestions.length > 0 && !isLoading && (
          <div className="flex gap-3 overflow-x-auto pb-4 mb-2 custom-scrollbar mask-linear-fade">
            {suggestions.map((s, i) => (
              <button
                key={i}
                onClick={() => handleSend(null, s.text)}
                className="flex items-center gap-2 px-4 py-2 bg-slate-50 hover:bg-blue-50 border border-slate-100 hover:border-blue-200 rounded-xl whitespace-nowrap transition-all group active:scale-95"
              >
                <span className="text-lg grayscale group-hover:grayscale-0 transition-all">{s.icon}</span>
                <span className="text-xs font-bold text-slate-600 group-hover:text-blue-700">{s.text}</span>
              </button>
            ))}
          </div>
        )}

        <form onSubmit={handleSend} className="relative group">
          <input
            type="text"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            placeholder="Pregunta sobre genes, codones, o análisis evolutivo..."
            className="w-full pl-6 pr-32 py-5 bg-slate-50 border-2 border-slate-100 rounded-[2rem] text-sm font-medium text-slate-800 focus:outline-none focus:border-blue-500/50 focus:bg-white focus:ring-4 focus:ring-blue-500/5 transition-all shadow-inner placeholder-slate-400"
          />
          <button
            type="submit"
            disabled={isLoading || !input.trim()}
            className="absolute right-2 top-2 bottom-2 px-8 bg-slate-900 text-white rounded-[1.5rem] font-black uppercase tracking-widest text-[10px] hover:bg-blue-600 transition-all shadow-lg hover:shadow-blue-500/30 disabled:opacity-50 disabled:shadow-none active:scale-95"
          >
            {isLoading ? '...' : 'Enviar'}
          </button>
        </form>
        <p className="text-center mt-3 text-[9px] font-bold text-slate-300 uppercase tracking-widest">
          Alimentado por Modelos LLM Avanzados & NCBI Database
        </p>
      </div>
    </div>
  )
}
