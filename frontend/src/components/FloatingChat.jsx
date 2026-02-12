import { useState } from 'react'
import InteractiveChat from './InteractiveChat'

export default function FloatingChat({ hasAnalysis, currentGenome }) {
  const [isOpen, setIsOpen] = useState(false)

  return (
    <div className="fixed bottom-8 right-8 z-[100] flex flex-col items-end">
      {/* Chat Window */}
      {isOpen && (
        <div className="mb-4 w-[450px] max-w-[calc(100vw-2rem)] h-[600px] shadow-2xl animate-in slide-in-from-bottom-8 fade-in duration-300 origin-bottom-right">
          <InteractiveChat 
            hasAnalysis={hasAnalysis} 
            currentGenome={currentGenome} 
            isFloating={true}
            onClose={() => setIsOpen(false)}
          />
        </div>
      )}

      {/* Toggle Button */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className={`w-16 h-16 rounded-[2rem] flex items-center justify-center shadow-2xl transition-all duration-500 hover:rotate-12 active:scale-90 border-4 border-white ${
          isOpen ? 'bg-rose-500 text-white rotate-45' : 'bg-slate-900 text-white hover:bg-blue-600'
        }`}
        title={isOpen ? "Cerrar Chat" : "Abrir Asistente IA"}
      >
        {isOpen ? (
          <svg className="w-8 h-8" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M6 18L18 6M6 6l12 12" />
          </svg>
        ) : (
          <div className="relative">
            <span className="text-2xl">🧬</span>
            <div className="absolute -top-1 -right-1 w-3 h-3 bg-emerald-500 border-2 border-slate-900 rounded-full animate-pulse"></div>
          </div>
        )}
      </button>
    </div>
  )
}
