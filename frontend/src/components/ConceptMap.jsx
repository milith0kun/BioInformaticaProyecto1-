import React, { useEffect, useRef, useState } from 'react';
import mermaid from 'mermaid';
import { motion, AnimatePresence } from 'framer-motion';

const definitions = {
    "1. Célula": "Unidad estructural y funcional básica de los organismos, compuesta por componentes que permiten metabolismo, replicación y respuesta al entorno.",
    "2. Teoría celular": "Postulado biológico que establece que los organismos están formados por células y que la célula es la unidad fundamental de vida.",
    "3. Ciclo de vida celular": "Secuencia de estados por los que pasa una célula (crecimiento, replicación/división y muerte), regulada por procesos bioquímicos.",
    "4. Vía metabólica": "Red organizada de reacciones químicas catalizadas que sintetizan, degradan o señalizan procesos dentro de la célula.",
    "5. ADN": "Polímero de nucleótidos con bases A, T, G y C; normalmente bicatenario y químicamente estable; almacena información hereditaria.",
    "6. ARN": "Polímero de nucleótidos con bases A, U, G y C; usualmente monocatenario y más reactivo; participa en transferencia y procesamiento de información genética y en funciones catalíticas/regulatorias.",
    "7. Proteína": "Macromolécula formada por una o más cadenas polipeptídicas de aminoácidos; ejecuta funciones celulares (catálisis, estructura, transporte, señalización, regulación).",
    "8. Alfabeto molecular": "Representación secuencial de macromoléculas como \"cadenas\" (ADN/ARN: 4 letras; proteínas: 20 letras), base conceptual para análisis computacional de secuencias.",
    "9. Cromosoma": "Estructura que organiza y empaqueta ADN (y proteínas asociadas) y contiene genes; su número varía entre especies.",
    "10. Gen": "Segmento de ADN que contiene la información necesaria para producir un producto funcional (ARN o proteína) bajo un contexto regulatorio.",
    "11. Herencia mendeliana": "Transmisión de rasgos mediante unidades discretas (genes) que segregan y se combinan según reglas observables en la descendencia.",
    "12. Mutación": "Cambio en la secuencia de ADN (p. ej., sustitución de base) que puede alterar un rasgo o función biológica.",
    "13. Ligamiento genético": "Tendencia de genes cercanos en un cromosoma a heredarse juntos, debido a menor probabilidad de recombinación entre ellos.",
    "14. Mapa genético": "Ordenamiento relativo de genes en un cromosoma estimado a partir de frecuencias de recombinación (distancia genética).",
    "15. Un gen–una proteína": "Hypótesis histórica que asocia cada gen con la producción de una proteína; hoy se reconoce que un gen puede originar múltiples productos (p. ej., por splicing alternativo).",
    "16. Nucleótido": "Unidad básica de ADN/ARN compuesta por una base nitrogenada, un azúcar y un grupo fosfato.",
    "17. Bases nitrogenadas": "A (adenina), T (timina, solo ADN), U (uracilo, solo ARN), G (guanina), C (citosina).",
    "18. Regla de Chargaff": "Observación de que en ADN bicatenario la cantidad de A≈T y G≈C, consistente con apareamiento complementario.",
    "19. Doble hélice": "Estructura del ADN formada por dos hebras antiparalelas enrolladas, unidas por puentes de hidrógeno entre pares de bases.",
    "20. Complementariedad de bases": "Emparejamiento específico A–T (o A–U en ARN) y C–G; permite copia fiel en replicación/transcripción e hibridación.",
    "21. Replicación del ADN": "Proceso de copia del ADN donde cada hebra sirve como molde para sintetizar una hebra complementaria.",
    "22. Célula eucariota": "Célula con ADN encapsulado en un núcleo; en general, genes con exones e intrones y procesamiento de ARN.",
    "23. Célula procariota": "Célula sin núcleo; genes típicamente continuos; transcripción y traducción ocurren en el mismo compartimento.",
    "24. Exón": "Segmento de un gen eucariota que permanece en el ARNm maduro y contribuye a la secuencia codificante.",
    "25. Intrón": "Segmento interveniente transcrito en ARN pero eliminado durante el splicing antes de formar el ARNm maduro.",
    "26. Transcripción": "Síntesis de ARN a partir de un molde de ADN mediante una enzima; copia información de un gen a una molécula de ARN.",
    "27. ARN polimerasa": "Complejo enzimático que cataliza la transcripción, agregando ribonucleótidos complementarios al molde de ADN.",
    "28. ARN mensajero": "ARN que porta la información codificante de un gen hacia el ribosoma para síntesis proteica.",
    "29. Splicing": "Proceso eucariota de eliminación de intrones y unión de exones para generar ARNm funcional.",
    "30. Traducción": "Proceso por el cual el ribosoma lee codones del ARNm y ensambla una proteína incorporando aminoácidos en orden.",
    "31. Ribosoma": "Complejo ribonucleoproteico (ARN + proteínas) que cataliza la traducción y ensambla polipéptidos.",
    "32. Aminoácido": "Monómero de las proteínas; existen 20 tipos estándar codificados por el código genético.",
    "33. Polipéptido": "Cadena lineal de aminoácidos unida por enlaces peptídicos; puede plegarse para formar una proteína funcional.",
    "34. Codón": "Triplete de nucleótidos en ARNm que especifica un aminoácido o una señal de inicio/terminación.",
    "35. Código genético": "Correspondencia entre codones y aminoácidos (y señales start/stop) usada para traducir ARNm a proteína.",
    "36. Degeneración del código genético": "Propiedad por la cual múltiples codones distintos pueden codificar el mismo aminoácido.",
    "37. Codón de inicio": "Codón que inicia la traducción; típicamente AUG (Metionina) en el código estándar.",
    "38. Codón de terminación": "Codón que señala el fin de la traducción; típicamente UAA, UAG o UGA.",
    "39. ARNt": "ARN adaptador que transporta un aminoácido específico y reconoce un codón del ARNm mediante su anticodón.",
    "40. Anticodón": "Triplete en ARNt complementario al codón del ARNm; asegura incorporación del aminoácido correcto.",
    "41. Dogma central": "Principio que describe el flujo principal de información genética: ADN → ARN → proteína.",
    "42. PCR": "Técnica de amplificación de un fragmento específico de ADN mediante ciclos de desnaturalización, alineamiento de cebadores y extensión por ADN polimerasa.",
    "43. Cebador": "Oligonucleótido corto que se aparea con el molde y provee un extremo 3' para que la polimerasa inicie síntesis.",
    "44. Desnaturalización": "Separación de hebras de ADN por calentamiento para obtener moldes monocatenarios.",
    "45. Alineamiento": "Unión de cebadores a secuencias complementarias flanqueantes al enfriar la reacción.",
    "46. Extensión": "Síntesis de nuevas hebras complementarias a partir de cebadores mediante ADN polimerasa.",
    "47. Clonación molecular": "Inserción de un fragmento de ADN en un vector replicable para producir muchas copias en un organismo hospedero.",
    "48. Vector de clonación": "Molécula de ADN (p. ej., plasmidio/virus) capaz de replicarse y transportar un inserto de ADN.",
    "49. Biblioteca de clones": "Colección de clones que representan fragmentos (frecuentemente aleatorios) de un genoma.",
    "50. Enzima de restricción": "Proteína que corta ADN en secuencias específicas (sitios de reconocimiento), generando fragmentos.",
    "51. Sitio de reconocimiento": "Secuencia específica (a menudo palindrómica) donde una enzima de restricción se une y corta.",
    "52. Extremos romos": "Resultado de un corte que deja extremos sin salientes monocatenarios.",
    "53. Extremos pegajosos": "Extremos con salientes monocatenarios complementarios que facilitan unión por hibridación.",
    "54. Hibridación": "Unión de hebras complementarias de ácidos nucleicos por apareamiento de bases.",
    "55. Ligación": "Formación de enlaces covalentes para \"sellar\" el esqueleto azúcar-fosfato y unir fragmentos de ADN.",
    "56. Electroforesis en gel": "Técnica que separa fragmentos de ADN por tamaño al migrar en un gel bajo campo eléctrico.",
    "57. Sonda": "Oligonucleótido marcado de secuencia conocida que se une por hibridación para detectar una secuencia complementaria.",
    "58. Microarreglo": "Superficie con muchas sondas inmovilizadas usada para detectar presencia o abundancia de transcritos (expresión génica) por hibridación.",
    "59. Variación genética intraespecífica": "Diferencias de secuencia entre individuos de una misma especie (p. ej., cambios en un pequeño porcentaje de bases) que sustentan diversidad de rasgos.",
    "60. Genoma de referencia": "Secuencia representativa (\"maestra\") usada como patrón para una especie, aunque individuos difieran en posiciones específicas.",
    "61. Conservación genética": "Presencia de genes o secuencias similares entre especies, indicativa de origen común o función preservada.",
    "62. Evolución": "Proceso de cambio heredable en poblaciones a través del tiempo, reflejado en variaciones de secuencia genómica.",
    "63. Selección natural": "Mecanismo por el cual variantes genéticas que mejoran el éxito reproductivo tienden a aumentar su frecuencia.",
    "64. Adaptación": "Incremento de frecuencia de rasgos heredables que mejoran el ajuste de una población a su ambiente.",
    "65. Especiación": "Formación de nuevas especies cuando poblaciones divergen hasta volverse reproductivamente incompatibles.",
    "66. Genómica comparativa": "Análisis comparado de genomas para identificar similitudes/diferencias, inferir funciones y relaciones evolutivas.",
    "67. Alineamiento de secuencias": "Método computacional para comparar secuencias y localizar regiones homólogas o conservadas.",
    "68. BLAST": "Familia de algoritmos rápidos para buscar similitudes entre secuencias y evaluar su significancia estadística.",
    "69. Bioinformática": "Disciplina que utiliza algoritmos, estadística y modelos computacionales para analizar datos biológicos (especialmente secuencias) y \"decodificar\" patrones funcionales y evolutivos."
};

const ConceptMap = () => {
    const [selectedTerm, setSelectedTerm] = useState(null);
    const [searchTerm, setSearchTerm] = useState('');
    const [searchResults, setSearchResults] = useState([]);
    const [isLoading, setIsLoading] = useState(true);
    const [zoom, setZoom] = useState(0.7);
    const [pan, setPan] = useState({ x: 20, y: 20 });
    const [isDragging, setIsDragging] = useState(false);
    const [startPos, setStartPos] = useState({ x: 0, y: 0 });
    
    const chartRef = useRef(null);
    const containerRef = useRef(null);

    useEffect(() => {
        mermaid.initialize({
            startOnLoad: false,
            theme: 'default',
            securityLevel: 'loose',
            flowchart: { useMaxWidth: false, htmlLabels: true, curve: 'basis' }
        });
        renderDiagram();
    }, []);

    const renderDiagram = async () => {
        setIsLoading(true);
        setError(null);

        const graphDefinition = `
flowchart LR
    classDef default fill:#fff,stroke:#cbd5e1,stroke-width:2px,color:#1e293b,font-weight:bold,rx:10,ry:10;
    classDef highlight fill:#eff6ff,stroke:#3b82f6,stroke-width:2px,color:#1d4ed8;
    classDef process fill:#f0fdf4,stroke:#22c55e,stroke-width:2px,color:#15803d;

    C1["1. Célula"]:::highlight -->|sustenta| C2["2. Teoría celular"]
    C1 -->|presenta| C3["3. Ciclo de vida celular"]
    C1 -->|ejecuta| C4["4. Vía metabólica"]
    C1 -->|contiene| C5["5. ADN"]:::highlight
    C1 -->|contiene| C6["6. ARN"]:::highlight
    C1 -->|produce| C7["7. Proteína"]:::highlight
    C5 -->|se organiza en| C9["9. Cromosoma"]
    C5 -->|incluye| C10["10. Gen"]:::highlight
    C10 -->|explica| C11["11. Herencia mendeliana"]
    C5 -->|sufre| C12["12. Mutación"]
    C10 -->|presenta| C13["13. Ligamiento genético"]
    C10 -->|se ordena en| C14["14. Mapa genético"]
    C10 -->|se asocia con| C15["15. Un gen–una proteína"]
    C5 -->|está formado por| C16["16. Nucleótido"]
    C16 -->|incluye| C17["17. Bases nitrogenadas"]
    C5 -->|obedece| C18["18. Regla de Chargaff"]
    C5 -->|adopta| C19["19. Doble hélice"]
    C5 -->|usa| C20["20. Complementariedad de bases"]
    C5 -->|se copia por| C21["21. Replicación del ADN"]:::process
    C1 -->|puede ser| C22["22. Célula eucariota"]
    C1 -->|puede ser| C23["23. Célula procariota"]
    C10 -->|se transcribe mediante| C26["26. Transcripción"]:::process
    C26 -->|es catalizada por| C27["27. ARN polimerasa"]
    C26 -->|produce| C28["28. ARN mensajero"]
    C28 -->|es modificado por| C29["29. Splicing"]:::process
    C28 -->|contiene| C24["24. Exón"]
    C28 -->|elimina| C25["25. Intrón"]
    C28 -->|se traduce por| C30["30. Traducción"]:::process
    C30 -->|es realizada por| C31["31. Ribosoma"]
    C30 -->|usa| C32["32. Aminoácido"]
    C32 -->|forma| C33["33. Polipéptido"]
    C30 -->|produce| C7
    C28 -->|se lee en| C34["34. Codón"]
    C34 -->|pertenece al| C35["35. Código genético"]
    C35 -->|presenta| C36["36. Degeneración del código genético"]
    C34 -->|puede ser| C37["37. Codón de inicio"]
    C34 -->|puede ser| C38["38. Codón de terminación"]
    C30 -->|requiere| C39["39. ARNt"]
    C39 -->|usa| C40["40. Anticodón"]
    C5 -->|fluye según| C41["41. Dogma central"]:::highlight
    C6 -->|fluye según| C41
    C7 -->|fluye según| C41
    C5 -->|se representa en| C8["8. Alfabeto molecular"]
    C6 -->|se representa en| C8
    C7 -->|se representa en| C8
    C5 -->|se amplifica por| C42["42. PCR"]:::process
    C42 -->|requiere| C43["43. Cebador"]
    C42 -->|incluye| C44["44. Desnaturalización"]
    C42 -->|incluye| C45["45. Alineamiento"]
    C42 -->|incluye| C46["46. Extensión"]
    C5 -->|se inserta mediante| C47["47. Clonación molecular"]:::process
    C47 -->|utiliza| C48["48. Vector de clonación"]
    C47 -->|genera| C49["49. Biblioteca de clones"]
    C5 -->|es cortado por| C50["50. Enzima de restricción"]
    C50 -->|reconoce| C51["51. Sitio de reconocimiento"]
    C50 -->|produce| C52["52. Extremos romos"]
    C50 -->|produce| C53["53. Extremos pegajosos"]
    C5 -->|se une por| C54["54. Hibridación"]
    C5 -->|se sella por| C55["55. Ligación"]
    C5 -->|se separa por| C56["56. Electroforesis en gel"]
    C5 -->|se detecta con| C57["57. Sonda"]
    C10 -->|se analiza con| C58["58. Microarreglo"]
    C5 -->|presenta| C59["59. Variación genética intraespecífica"]
    C5 -->|se compara con| C60["60. Genoma de referencia"]
    C59 -->|contribuye a| C61["61. Conservación genética"]
    C61 -->|evidencia| C62["62. Evolución"]:::highlight
    C62 -->|opera por| C63["63. Selección natural"]
    C63 -->|favorece| C64["64. Adaptación"]
    C64 -->|origina| C65["65. Especiación"]
    C5 -->|se estudia con| C66["66. Genómica comparativa"]
    C5 -->|se compara mediante| C67["67. Alineamiento de secuencias"]
    C67 -->|se implementa en| C68["68. BLAST"]
    C66 -->|se apoya en| C69["69. Bioinformática"]:::highlight
`;
        try {
            const id = `mermaid-id-${Math.random().toString(36).substr(2, 9)}`;
            const { svg } = await mermaid.render(id, graphDefinition);
            
            if (chartRef.current) {
                chartRef.current.innerHTML = svg;
                const svgElement = chartRef.current.querySelector('svg');
                if (svgElement) {
                    svgElement.style.maxWidth = 'none';
                    svgElement.style.height = 'auto';
                    svgElement.style.width = 'auto';
                }
                attachListeners();
                setTimeout(autoFit, 100);
            }
        } catch (e) { 
            console.error("Mermaid Render Error:", e); 
            setError("Error al generar el mapa. Intente recargar.");
        } finally { setIsLoading(false); }
    };

    const autoFit = () => {
        if (!chartRef.current || !containerRef.current) return;
        const svgElement = chartRef.current.querySelector('svg');
        if (!svgElement) return;

        const cw = containerRef.current.clientWidth || 800;
        const ch = containerRef.current.clientHeight || 600;
        
        try {
            const bbox = svgElement.getBBox();
            if (bbox.width === 0 || bbox.height === 0) return;

            const scaleX = (cw - 60) / bbox.width;
            const scaleY = (ch - 60) / bbox.height;
            const scale = Math.min(scaleX, scaleY, 1);
            
            setZoom(scale);
            setPan({
                x: (cw - bbox.width * scale) / 2 - bbox.x * scale,
                y: (ch - bbox.height * scale) / 2 - bbox.y * scale
            });
        } catch (e) {
            console.warn("BBox calculation failed, using defaults");
            setZoom(0.8);
            setPan({ x: 20, y: 20 });
        }
    };

    const attachListeners = () => {
        if (!chartRef.current) return;
        const nodes = chartRef.current.querySelectorAll('.node');
        nodes.forEach(n => {
            n.style.cursor = 'pointer';
            n.onclick = (e) => {
                e.stopPropagation();
                const text = n.textContent.trim();
                const key = Object.keys(definitions).find(k => text.includes(k) || key?.includes(text));
                if (key) setSelectedTerm({ term: key, definition: definitions[key] });
            };
        });
    };

    const onMouseDown = (e) => {
        if (e.button !== 0) return;
        setIsDragging(true);
        setStartPos({ x: e.clientX - pan.x, y: e.clientY - pan.y });
    };
    const onMouseMove = (e) => {
        if (!isDragging) return;
        setPan({ x: e.clientX - startPos.x, y: e.clientY - startPos.y });
    };
    const onMouseUp = () => setIsDragging(false);

    const handleSearch = (e) => {
        const val = e.target.value;
        setSearchTerm(val);
        if (val.length < 2) { setSearchResults([]); return; }
        setSearchResults(Object.keys(definitions).filter(k => k.toLowerCase().includes(val.toLowerCase())));
    };

    const focusNode = (key) => {
        const nodes = Array.from(chartRef.current.querySelectorAll('.node'));
        const target = nodes.find(n => n.textContent.includes(key));
        if (target) {
            setZoom(1.1);
            setPan({ x: 50, y: 50 });
            target.style.filter = 'brightness(1.5) drop-shadow(0 0 15px #3b82f6)';
            setTimeout(() => target.style.filter = '', 3000);
        }
        setSearchTerm('');
        setSearchResults([]);
    };

    return (
        <div className="flex flex-col h-full bg-slate-50 relative overflow-hidden" ref={containerRef}>
            {/* Header */}
            <div className="p-6 bg-white border-b flex flex-col md:flex-row items-center justify-between gap-4 z-20 shrink-0 shadow-sm">
                <div className="flex items-center gap-4">
                    <div className="w-10 h-10 bg-blue-600 rounded-xl flex items-center justify-center text-white shadow-lg">
                        <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" /></svg>
                    </div>
                    <div>
                        <h2 className="text-sm font-black text-slate-900 uppercase tracking-tighter">Mapa Conceptual Pro</h2>
                        <p className="text-[9px] font-bold text-slate-400 uppercase tracking-widest">Navegación Dinámica</p>
                    </div>
                </div>
                <div className="flex items-center gap-3 relative z-30">
                    <div className="relative">
                        <input type="text" value={searchTerm} onChange={handleSearch} placeholder="Buscar..." className="w-48 pl-10 pr-4 py-2 bg-slate-50 border rounded-xl text-xs font-bold focus:ring-2 focus:ring-blue-500/20 outline-none" />
                        <svg className="w-4 h-4 text-slate-400 absolute left-3.5 top-1/2 -translate-y-1/2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" strokeWidth={3}/></svg>
                        {searchResults.length > 0 && (
                            <div className="absolute top-full left-0 right-0 mt-2 bg-white border rounded-xl shadow-2xl max-h-48 overflow-y-auto">
                                {searchResults.map(k => (
                                    <button key={k} onClick={() => focusNode(k)} className="w-full text-left px-4 py-2 hover:bg-blue-50 text-[10px] font-bold border-b last:border-0">{k}</button>
                                ))}
                            </div>
                        )}
                    </div>
                    <div className="flex bg-slate-100 p-1 rounded-xl">
                        <button onClick={() => setZoom(z => Math.max(z-0.1, 0.1))} className="p-1.5 hover:bg-white rounded-lg transition-all"><svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M20 12H4" strokeWidth={3}/></svg></button>
                        <button onClick={() => setZoom(z => Math.min(z+0.1, 3))} className="p-1.5 hover:bg-white rounded-lg transition-all"><svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path d="M12 4v16m8-8H4" strokeWidth={3}/></svg></button>
                    </div>
                    <button onClick={autoFit} className="px-3 py-2 bg-slate-900 text-white text-[9px] font-black uppercase rounded-xl">Reset</button>
                </div>
            </div>

            {/* Viewport */}
            <div 
                className={`flex-1 relative bg-slate-50 cursor-grab active:cursor-grabbing select-none overflow-hidden ${isDragging ? 'cursor-grabbing' : ''}`}
                onMouseDown={onMouseDown} onMouseMove={onMouseMove} onMouseUp={onMouseUp} onMouseLeave={onMouseUp}
            >
                {isLoading && (
                    <div className="absolute inset-0 flex items-center justify-center bg-white/80 z-30"><div className="w-8 h-8 border-4 border-blue-600 border-t-transparent rounded-full animate-spin"></div></div>
                )}
                <div 
                    style={{ 
                        transform: `translate(${pan.x}px, ${pan.y}px) scale(${zoom})`, 
                        transformOrigin: '0 0',
                        transition: isDragging ? 'none' : 'transform 0.3s cubic-bezier(0.2, 0, 0, 1)',
                        position: 'absolute',
                        top: 0,
                        left: 0,
                        willChange: 'transform'
                    }}
                    ref={chartRef}
                ></div>
            </div>

            {/* Modal */}
            <AnimatePresence>
                {selectedTerm && (
                    <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }} onClick={() => setSelectedTerm(null)} className="fixed inset-0 z-[100] flex items-center justify-center bg-slate-900/60 backdrop-blur-sm p-6">
                        <motion.div initial={{ scale: 0.9, y: 20 }} animate={{ scale: 1, y: 0 }} onClick={(e) => e.stopPropagation()} className="bg-white rounded-[2rem] shadow-2xl max-w-lg w-full overflow-hidden">
                            <div className="bg-blue-600 p-6 text-white"><h3 className="text-xl font-black italic tracking-tighter leading-none">{selectedTerm.term}</h3></div>
                            <div className="p-8">
                                <p className="text-slate-600 font-medium leading-relaxed italic text-lg">"{selectedTerm.definition}"</p>
                                <button onClick={() => setSelectedTerm(null)} className="mt-8 w-full py-4 bg-slate-900 text-white text-[10px] font-black uppercase rounded-xl">Cerrar</button>
                            </div>
                        </motion.div>
                    </motion.div>
                )}
            </AnimatePresence>
        </div>
    );
};

export default ConceptMap;
