import { useState } from 'react'

export default function FloatingMapButton({ onOpenConcept, onOpenNetwork, activeView }) {
  const [isOpen, setIsOpen] = useState(false)

  return (
    <div className="fixed bottom-24 right-8 z-[60] flex flex-col items-end gap-3">
      {/* Menu Options */}
      {isOpen && (
        <div className="flex flex-col gap-3 mb-2 animate-in slide-in-from-bottom-4 fade-in duration-300">
          <button
            onClick={() => {
              onOpenNetwork()
              setIsOpen(false)
            }}
            className={`flex items-center gap-3 px-6 py-3 rounded-2xl shadow-2xl transition-all hover:-translate-x-2 active:scale-95 ${
              activeView === 'network-map' 
              ? 'bg-emerald-600 text-white border-2 border-emerald-400' 
              : 'bg-slate-900 text-white hover:bg-emerald-600'
            }`}
          >
            <span className="text-[10px] font-black uppercase tracking-widest whitespace-nowrap">Red Genómica Pro</span>
            <div className="w-8 h-8 rounded-xl bg-white/10 flex items-center justify-center">
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
            </div>
          </button>

          <button
            onClick={() => {
              onOpenConcept()
              setIsOpen(false)
            }}
            className={`flex items-center gap-3 px-6 py-3 rounded-2xl shadow-2xl transition-all hover:-translate-x-2 active:scale-95 ${
              activeView === 'concept-map' 
              ? 'bg-blue-600 text-white border-2 border-blue-400' 
              : 'bg-slate-900 text-white hover:bg-blue-600'
            }`}
          >
            <span className="text-[10px] font-black uppercase tracking-widest whitespace-nowrap">Mapa Conceptual</span>
            <div className="w-8 h-8 rounded-xl bg-white/10 flex items-center justify-center">
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" /></svg>
            </div>
          </button>
        </div>
      )}

      {/* Main Toggle Button */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className={`w-16 h-16 rounded-[2rem] flex items-center justify-center shadow-2xl transition-all duration-500 hover:rotate-12 active:scale-90 border-4 border-white ${
          isOpen ? 'bg-rose-500 text-white rotate-45' : 'bg-blue-600 text-white'
        }`}
      >
        <svg className="w-8 h-8" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d={isOpen ? "M6 18L18 6M6 6l12 12" : "M9 20l-5.447-2.724A1 1 0 013 16.382V5.618a1 1 0 011.447-.894L9 7m0 13l6-3m-6 3V7m6 10l4.553 2.276A1 1 0 0021 18.382V7.618a1 1 0 00-.553-.894L15 4m0 13V4m0 0L9 7"} />
        </svg>
      </button>
    </div>
  )
}
