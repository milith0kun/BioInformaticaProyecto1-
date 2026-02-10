/**
 * InteractiveChat Component
 * Full-screen conversational AI interface with a bioinformatics expert
 * Connected to GenBank, PubMed, and genome data for citations
 */
import { useState, useEffect, useRef } from 'react'
import { api } from '../services/api'

const TYPING_INDICATOR = '●●●'

export default function InteractiveChat({ hasAnalysis, currentGenome }) {
    const [messages, setMessages] = useState([])
    const [input, setInput] = useState('')
    const [isLoading, setIsLoading] = useState(false)
    const [suggestions, setSuggestions] = useState([])
    const [showReferences, setShowReferences] = useState(null)
    const messagesEndRef = useRef(null)
    const inputRef = useRef(null)

    useEffect(() => {
        loadSuggestions()
        loadHistory()
    }, [])

    useEffect(() => {
        scrollToBottom()
    }, [messages])

    const scrollToBottom = () => {
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' })
    }

    const loadSuggestions = async () => {
        try {
            const result = await api.getChatSuggestions()
            setSuggestions(result.suggestions || [])
        } catch (e) {
            setSuggestions([
                { text: '¿Qué puedes decirme sobre este genoma?', icon: '🧬' },
                { text: 'Explica el uso de codones en bacterias', icon: '📊' },
                { text: '¿Qué significa el contenido GC?', icon: '🔬' },
            ])
        }
    }

    const loadHistory = async () => {
        try {
            const result = await api.getChatHistory()
            if (result.messages && result.messages.length > 0) {
                setMessages(result.messages.map((m, i) => ({
                    id: i,
                    role: m.role,
                    content: m.content,
                    timestamp: new Date().toLocaleTimeString('es-ES', { hour: '2-digit', minute: '2-digit' })
                })))
            }
        } catch (e) {
            console.error('Error loading chat history:', e)
        }
    }

    const sendMessage = async (text) => {
        const messageText = text || input.trim()
        if (!messageText || isLoading) return

        const userMsg = {
            id: Date.now(),
            role: 'user',
            content: messageText,
            timestamp: new Date().toLocaleTimeString('es-ES', { hour: '2-digit', minute: '2-digit' })
        }

        setMessages(prev => [...prev, userMsg])
        setInput('')
        setIsLoading(true)

        try {
            const result = await api.sendChatMessage(messageText)

            const assistantMsg = {
                id: Date.now() + 1,
                role: 'assistant',
                content: result.response,
                timestamp: new Date().toLocaleTimeString('es-ES', { hour: '2-digit', minute: '2-digit' }),
                references: result.references || [],
                ncbi_links: result.ncbi_links || [],
                genome_context_used: result.genome_context_used
            }

            setMessages(prev => [...prev, assistantMsg])
        } catch (error) {
            const errorMsg = {
                id: Date.now() + 1,
                role: 'assistant',
                content: '❌ Error al comunicarse con el asistente. Verifica tu conexión y la API key.',
                timestamp: new Date().toLocaleTimeString('es-ES', { hour: '2-digit', minute: '2-digit' }),
                isError: true
            }
            setMessages(prev => [...prev, errorMsg])
        } finally {
            setIsLoading(false)
            inputRef.current?.focus()
        }
    }

    const clearChat = async () => {
        try {
            await api.clearChatHistory()
            setMessages([])
        } catch (e) {
            setMessages([])
        }
    }

    const handleKeyDown = (e) => {
        if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault()
            sendMessage()
        }
    }

    // Simple markdown renderer
    const renderMarkdown = (text) => {
        if (!text) return ''

        // Process markdown
        let html = text
            // Code blocks
            .replace(/```(\w+)?\n([\s\S]*?)```/g, '<pre class="bg-slate-900 text-green-400 rounded-lg p-4 overflow-x-auto text-sm my-3 font-mono"><code>$2</code></pre>')
            // Inline code
            .replace(/`([^`]+)`/g, '<code class="bg-slate-200 text-teal-800 px-1.5 py-0.5 rounded text-sm font-mono">$1</code>')
            // Links - markdown format [text](url)
            .replace(/\[([^\]]+)\]\((https?:\/\/[^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener noreferrer" class="text-teal-600 underline hover:text-teal-800 font-medium">$1</a>')
            // Bold
            .replace(/\*\*([^*]+)\*\*/g, '<strong class="font-semibold text-slate-900">$1</strong>')
            // Italic
            .replace(/\*([^*]+)\*/g, '<em class="italic">$1</em>')
            // Headers
            .replace(/^### (.+)$/gm, '<h4 class="text-base font-bold text-slate-800 mt-4 mb-2">$1</h4>')
            .replace(/^## (.+)$/gm, '<h3 class="text-lg font-bold text-slate-800 mt-4 mb-2">$1</h3>')
            // Lists
            .replace(/^- (.+)$/gm, '<li class="ml-4 list-disc text-slate-700">$1</li>')
            .replace(/^(\d+)\. (.+)$/gm, '<li class="ml-4 list-decimal text-slate-700">$2</li>')
            // Line breaks  
            .replace(/\n\n/g, '<br/><br/>')
            .replace(/\n/g, '<br/>')

        // Wrap consecutive <li> in <ul>
        html = html.replace(/((?:<li[^>]*>.*?<\/li><br\/?>)+)/g, '<ul class="my-2 space-y-1">$1</ul>')

        return html
    }

    return (
        <div className="flex flex-col h-[600px] sm:h-[700px]">
            {/* Chat Header */}
            <div className="flex items-center justify-between px-4 sm:px-6 py-3 bg-gradient-to-r from-slate-900 via-slate-800 to-teal-900 rounded-t-xl">
                <div className="flex items-center gap-3">
                    <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-teal-400 to-emerald-500 flex items-center justify-center text-white shadow-lg">
                        <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19.8 15.3M14.25 3.104c.251.023.501.05.75.082M19.8 15.3l-1.57.393A9.065 9.065 0 0112 15a9.065 9.065 0 00-6.23.693L5 14.5m14.8.8l1.402 1.402c1.232 1.232.65 3.318-1.067 3.611A48.309 48.309 0 0112 21c-2.773 0-5.491-.235-8.135-.687-1.718-.293-2.3-2.379-1.067-3.61L5 14.5" />
                        </svg>
                    </div>
                    <div>
                        <h3 className="text-white font-semibold text-sm sm:text-base">GenomicAI</h3>
                        <p className="text-teal-300 text-xs">Experto en Bioinformática • NCBI + PubMed</p>
                    </div>
                </div>
                <div className="flex items-center gap-2">
                    {currentGenome && (
                        <span className="hidden sm:inline-flex items-center gap-1.5 px-2.5 py-1 bg-teal-500/20 text-teal-300 rounded-lg text-xs">
                            <span className="w-1.5 h-1.5 bg-teal-400 rounded-full animate-pulse"></span>
                            {currentGenome.accession}
                        </span>
                    )}
                    <button
                        onClick={clearChat}
                        className="p-2 text-slate-400 hover:text-white hover:bg-white/10 rounded-lg transition-all"
                        title="Limpiar chat"
                    >
                        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
                        </svg>
                    </button>
                </div>
            </div>

            {/* Messages Area */}
            <div className="flex-1 overflow-y-auto p-4 sm:p-6 space-y-4 bg-gradient-to-b from-slate-50 to-white">
                {messages.length === 0 && (
                    <div className="flex flex-col items-center justify-center h-full text-center px-4">
                        <div className="w-20 h-20 rounded-2xl bg-gradient-to-br from-teal-100 to-emerald-100 flex items-center justify-center mb-6">
                            <span className="text-4xl">🧬</span>
                        </div>
                        <h3 className="text-lg font-bold text-slate-800 mb-2">¡Hola! Soy GenomicAI</h3>
                        <p className="text-slate-500 text-sm mb-6 max-w-md">
                            Soy tu asistente experto en bioinformática. Puedo analizar datos genómicos,
                            buscar referencias en PubMed, y ayudarte a interpretar resultados.
                        </p>

                        {/* Suggestions */}
                        <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 w-full max-w-lg">
                            {suggestions.slice(0, 6).map((s, i) => (
                                <button
                                    key={i}
                                    onClick={() => sendMessage(s.text)}
                                    className="flex items-center gap-2 px-4 py-3 bg-white border border-slate-200 rounded-xl text-left text-sm text-slate-700 hover:border-teal-300 hover:bg-teal-50/50 transition-all group"
                                >
                                    <span className="text-lg">{s.icon}</span>
                                    <span className="group-hover:text-teal-700 line-clamp-2">{s.text}</span>
                                </button>
                            ))}
                        </div>
                    </div>
                )}

                {messages.map((msg) => (
                    <div
                        key={msg.id}
                        className={`flex gap-3 ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}
                    >
                        {msg.role === 'assistant' && (
                            <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-teal-500 to-emerald-500 flex-shrink-0 flex items-center justify-center shadow-md">
                                <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19.8 15.3" />
                                </svg>
                            </div>
                        )}

                        <div className={`max-w-[80%] ${msg.role === 'user' ? 'order-1' : ''}`}>
                            <div
                                className={`rounded-2xl px-4 py-3 ${msg.role === 'user'
                                    ? 'bg-gradient-to-r from-teal-600 to-teal-700 text-white rounded-br-md'
                                    : msg.isError
                                        ? 'bg-red-50 text-red-700 border border-red-200 rounded-bl-md'
                                        : 'bg-white text-slate-700 border border-slate-200 shadow-sm rounded-bl-md'
                                    }`}
                            >
                                {msg.role === 'user' ? (
                                    <p className="text-sm leading-relaxed whitespace-pre-wrap">{msg.content}</p>
                                ) : (
                                    <div
                                        className="text-sm leading-relaxed prose-sm"
                                        dangerouslySetInnerHTML={{ __html: renderMarkdown(msg.content) }}
                                    />
                                )}
                            </div>

                            {/* References */}
                            {msg.references && msg.references.length > 0 && (
                                <div className="mt-2">
                                    <button
                                        onClick={() => setShowReferences(showReferences === msg.id ? null : msg.id)}
                                        className="flex items-center gap-1.5 text-xs text-teal-600 hover:text-teal-700 font-medium"
                                    >
                                        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
                                        </svg>
                                        {msg.references.length} referencia{msg.references.length > 1 ? 's' : ''} PubMed
                                        <svg className={`w-3 h-3 transition-transform ${showReferences === msg.id ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                                        </svg>
                                    </button>

                                    {showReferences === msg.id && (
                                        <div className="mt-2 space-y-2">
                                            {msg.references.map((ref, ri) => (
                                                <a
                                                    key={ri}
                                                    href={ref.url}
                                                    target="_blank"
                                                    rel="noopener noreferrer"
                                                    className="block p-3 bg-slate-50 rounded-lg border border-slate-200 hover:border-teal-300 hover:bg-teal-50/30 transition-all"
                                                >
                                                    <p className="text-xs font-medium text-slate-800 line-clamp-2">{ref.title}</p>
                                                    <p className="text-xs text-slate-500 mt-1">
                                                        {ref.authors} • {ref.journal} ({ref.year}) • PMID: {ref.pmid}
                                                    </p>
                                                </a>
                                            ))}
                                        </div>
                                    )}
                                </div>
                            )}

                            {/* NCBI Links */}
                            {msg.ncbi_links && msg.ncbi_links.length > 0 && (
                                <div className="mt-2 flex flex-wrap gap-1.5">
                                    {msg.ncbi_links.map((link, li) => (
                                        <a
                                            key={li}
                                            href={link.url}
                                            target="_blank"
                                            rel="noopener noreferrer"
                                            className="inline-flex items-center gap-1 px-2 py-1 bg-teal-50 border border-teal-200 rounded-lg text-xs text-teal-700 hover:bg-teal-100 transition-all"
                                        >
                                            {link.type === 'gene' ? '🧬' : '📄'}
                                            <span className="font-medium">{link.name}</span>
                                        </a>
                                    ))}
                                </div>
                            )}

                            {/* Genome context indicator */}
                            {msg.genome_context_used && msg.role === 'assistant' && (
                                <div className="flex items-center gap-1 mt-1.5 text-xs text-emerald-600">
                                    <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
                                        <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
                                    </svg>
                                    Datos del genoma + búsqueda NCBI
                                </div>
                            )}

                            <span className="text-[10px] text-slate-400 mt-1 block">{msg.timestamp}</span>
                        </div>

                        {msg.role === 'user' && (
                            <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-slate-600 to-slate-700 flex-shrink-0 flex items-center justify-center shadow-md order-2">
                                <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" />
                                </svg>
                            </div>
                        )}
                    </div>
                ))}

                {isLoading && (
                    <div className="flex gap-3 items-start">
                        <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-teal-500 to-emerald-500 flex-shrink-0 flex items-center justify-center shadow-md">
                            <svg className="w-4 h-4 text-white animate-pulse" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5" />
                            </svg>
                        </div>
                        <div className="bg-white rounded-2xl rounded-bl-md px-4 py-3 border border-slate-200 shadow-sm">
                            <div className="flex items-center gap-1">
                                <span className="w-2 h-2 bg-teal-500 rounded-full animate-bounce" style={{ animationDelay: '0ms' }}></span>
                                <span className="w-2 h-2 bg-teal-500 rounded-full animate-bounce" style={{ animationDelay: '150ms' }}></span>
                                <span className="w-2 h-2 bg-teal-500 rounded-full animate-bounce" style={{ animationDelay: '300ms' }}></span>
                            </div>
                        </div>
                    </div>
                )}

                <div ref={messagesEndRef} />
            </div>

            {/* Input Area */}
            <div className="px-4 sm:px-6 py-4 bg-white border-t border-slate-200 rounded-b-xl">
                <div className="flex items-end gap-3">
                    <div className="flex-1 relative">
                        <textarea
                            ref={inputRef}
                            value={input}
                            onChange={(e) => setInput(e.target.value)}
                            onKeyDown={handleKeyDown}
                            placeholder="Pregunta sobre el genoma, codones, genes, proteínas..."
                            rows={1}
                            className="w-full px-4 py-3 bg-slate-100 border border-slate-200 rounded-xl text-sm text-slate-800 placeholder-slate-400 resize-none focus:outline-none focus:ring-2 focus:ring-teal-400 focus:border-transparent transition-all"
                            style={{ minHeight: '44px', maxHeight: '120px' }}
                            onInput={(e) => {
                                e.target.style.height = 'auto'
                                e.target.style.height = Math.min(e.target.scrollHeight, 120) + 'px'
                            }}
                        />
                    </div>
                    <button
                        onClick={() => sendMessage()}
                        disabled={isLoading || !input.trim()}
                        className="w-11 h-11 flex-shrink-0 bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-xl flex items-center justify-center hover:from-teal-700 hover:to-emerald-700 disabled:opacity-40 transition-all shadow-md hover:shadow-lg"
                    >
                        <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8" />
                        </svg>
                    </button>
                </div>
                <p className="text-[10px] text-slate-400 mt-2 text-center">
                    GenomicAI usa Claude 3.5 + NCBI APIs. Conectado a PubMed y GenBank para referencias científicas.
                </p>
            </div>
        </div>
    )
}
