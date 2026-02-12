import React, { useEffect, useRef, useState } from 'react';
import mermaid from 'mermaid';
import { motion, AnimatePresence } from 'framer-motion';

// Definiciones del Mapa (Hardcoded from Mapa.html)
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
    "15. Un gen–una proteína": "Hipótesis histórica que asocia cada gen con la producción de una proteína; hoy se reconoce que un gen puede originar múltiples productos (p. ej., por splicing alternativo).",
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
    const [zoomLevel, setZoomLevel] = useState(1);
    const [selectedTerm, setSelectedTerm] = useState(null);
    const chartRef = useRef(null);
    const containerRef = useRef(null);

    // Initial render logic
    useEffect(() => {
        mermaid.initialize({
            startOnLoad: false,
            theme: 'default', // Or 'base' for more customization
            securityLevel: 'loose',
            flowchart: {
                useMaxWidth: false,
                htmlLabels: true,
                curve: 'basis'
            }
        });

        const renderDiagram = async () => {
            if (chartRef.current) {
                const graphDefinition = `
flowchart LR
%% NÚCLEO
C1["1. Célula"] -->|sustenta| C2["2. Teoría celular"]
C1 -->|presenta| C3["3. Ciclo de vida celular"]
C1 -->|ejecuta| C4["4. Vía metabólica"]
C1 -->|contiene| C5["5. ADN"]
C1 -->|contiene| C6["6. ARN"]
C1 -->|produce| C7["7. Proteína"]
%% INFORMACIÓN GENÉTICA
C5 -->|se organiza en| C9["9. Cromosoma"]
C5 -->|incluye| C10["10. Gen"]
C10 -->|explica| C11["11. Herencia mendeliana"]
C5 -->|sufre| C12["12. Mutación"]
C10 -->|presenta| C13["13. Ligamiento genético"]
C10 -->|se ordena en| C14["14. Mapa genético"]
C10 -->|se asocia con| C15["15. Un gen–una proteína"]
%% ESTRUCTURA MOLECULAR
C5 -->|está formado por| C16["16. Nucleótido"]
C16 -->|incluye| C17["17. Bases nitrogenadas"]
C5 -->|obedece| C18["18. Regla de Chargaff"]
C5 -->|adopta| C19["19. Doble hélice"]
C5 -->|usa| C20["20. Complementariedad de bases"]
C5 -->|se copia por| C21["21. Replicación del ADN"]
%% TIPOS CELULARES
C1 -->|puede ser| C22["22. Célula eucariota"]
C1 -->|puede ser| C23["23. Célula procariota"]
%% EXPRESIÓN GÉNICA
C10 -->|se transcribe mediante| C26["26. Transcripción"]
C26 -->|es catalizada por| C27["27. ARN polimerasa"]
C26 -->|produce| C28["28. ARN mensajero"]
C28 -->|es modificado por| C29["29. Splicing"]
C28 -->|contiene| C24["24. Exón"]
C28 -->|elimina| C25["25. Intrón"]
%% TRADUCCIÓN
C28 -->|se traduce por| C30["30. Traducción"]
C30 -->|es realizada por| C31["31. Ribosoma"]
C30 -->|usa| C32["32. Aminoácido"]
C32 -->|forma| C33["33. Polipéptido"]
C30 -->|produce| C7
%% CÓDIGO GENÉTICO
C28 -->|se lee en| C34["34. Codón"]
C34 -->|pertenece al| C35["35. Código genético"]
C35 -->|presenta| C36["36. Degeneración del código genético"]
C34 -->|puede ser| C37["37. Codón de inicio"]
C34 -->|puede ser| C38["38. Codón de terminación"]
C30 -->|requiere| C39["39. ARNt"]
C39 -->|usa| C40["40. Anticodón"]
%% MARCO CONCEPTUAL
C5 -->|fluye según| C41["41. Dogma central"]
C6 -->|fluye según| C41
C7 -->|fluye según| C41
C5 -->|se representa en| C8["8. Alfabeto molecular"]
C6 -->|se representa en| C8
C7 -->|se representa en| C8
%% TÉCNICAS MOLECULARES
C5 -->|se amplifica por| C42["42. PCR"]
C42 -->|requiere| C43["43. Cebador"]
C42 -->|incluye| C44["44. Desnaturalización"]
C42 -->|incluye| C45["45. Alineamiento"]
C42 -->|incluye| C46["46. Extensión"]
C5 -->|se inserta mediante| C47["47. Clonación molecular"]
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
%% VARIACIÓN Y EVOLUCIÓN
C5 -->|presenta| C59["59. Variación genética intraespecífica"]
C5 -->|se compara con| C60["60. Genoma de referencia"]
C59 -->|contribuye a| C61["61. Conservación genética"]
C61 -->|evidencia| C62["62. Evolución"]
C62 -->|opera por| C63["63. Selección natural"]
C63 -->|favorece| C64["64. Adaptación"]
C64 -->|origina| C65["65. Especiación"]
%% ANÁLISIS COMPUTACIONAL
C5 -->|se estudia con| C66["66. Genómica comparativa"]
C5 -->|se compara mediante| C67["67. Alineamiento de secuencias"]
C67 -->|se implementa en| C68["68. BLAST"]
C66 -->|se apoya en| C69["69. Bioinformática"]
`;

                try {
                    const { svg } = await mermaid.render('mermaid-chart', graphDefinition);
                    chartRef.current.innerHTML = svg;
                    fitToScreen();
                } catch (error) {
                    console.error('Mermaid render error:', error);
                }
            }
        };

        renderDiagram();
    }, []);

    // Add click listeners to nodes
    useEffect(() => {
        const handleNodeClick = (e) => {
            const node = e.target.closest('.node');
            if (node) {
                // Try to find text content
                const textElement = node.querySelector('span, p, div') || node;
                const text = textElement.textContent.trim();

                // Match against keys in definitions
                // Sometimes mermaid adds extra whitespace or newlines
                const matchingKey = Object.keys(definitions).find(key => text.includes(key) || key.includes(text));

                if (matchingKey) {
                    setSelectedTerm({
                        term: matchingKey,
                        definition: definitions[matchingKey]
                    });
                }
            }
        };

        const chartElement = chartRef.current;
        if (chartElement) {
            chartElement.addEventListener('click', handleNodeClick);
        }

        return () => {
            if (chartElement) {
                chartElement.removeEventListener('click', handleNodeClick);
            }
        };
    }, []);

    const zoomIn = () => setZoomLevel(prev => Math.min(prev + 0.1, 3));
    const zoomOut = () => setZoomLevel(prev => Math.max(prev - 0.1, 0.3));
    const resetZoom = () => setZoomLevel(1);

    const fitToScreen = () => {
        if (containerRef.current && chartRef.current) {
            // Simple fit logic - adjust as needed based on actual dimensions
            setZoomLevel(0.8);
        }
    };

    return (
        <div className="flex flex-col h-full bg-slate-50 relative overflow-hidden">
            {/* Header / Title Area */}
            <div className="absolute top-0 left-0 right-0 z-10 p-6 pointer-events-none">
                <div className="max-w-4xl mx-auto text-center pointer-events-auto">
                    <motion.div
                        initial={{ opacity: 0, y: -20 }}
                        animate={{ opacity: 1, y: 0 }}
                        className="inline-block bg-white/90 backdrop-blur-md border border-slate-200 shadow-lg rounded-2xl px-8 py-4"
                    >
                        <h1 className="text-2xl font-black text-slate-900 tracking-tighter uppercase mb-1">
                            🧬 Mapa Conceptual Interactivo
                        </h1>
                        <p className="text-xs font-bold text-slate-500 uppercase tracking-widest">
                            Biología Molecular - 69 Términos Conectados
                        </p>
                    </motion.div>
                </div>
            </div>

            {/* Controls */}
            <div className="absolute bottom-6 left-1/2 -translate-x-1/2 z-10 flex items-center gap-2 bg-white/90 backdrop-blur-md p-2 rounded-2xl shadow-xl border border-slate-200">
                <ControlButton onClick={zoomOut} icon={
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M20 12H4" /></svg>
                } />
                <span className="w-16 text-center font-mono text-xs font-bold text-slate-600">
                    {Math.round(zoomLevel * 100)}%
                </span>
                <ControlButton onClick={zoomIn} icon={
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" /></svg>
                } />
                <div className="w-px h-6 bg-slate-200 mx-1"></div>
                <ControlButton onClick={resetZoom} label="1:1" />
                <ControlButton onClick={fitToScreen} icon={
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4" /></svg>
                } />
            </div>

            {/* Diagram Container */}
            <div
                ref={containerRef}
                className="flex-1 overflow-auto cursor-move bg-slate-50 dna-pattern relative"
            >
                <div
                    style={{
                        transform: `scale(${zoomLevel})`,
                        transformOrigin: 'top center',
                        transition: 'transform 0.3s ease-out'
                    }}
                    className="min-w-full min-h-full p-20 flex justify-center"
                >
                    <div ref={chartRef} className="mermaid-chart"></div>
                </div>
            </div>

            {/* Definition Modal */}
            <AnimatePresence>
                {selectedTerm && (
                    <motion.div
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        exit={{ opacity: 0 }}
                        onClick={() => setSelectedTerm(null)}
                        className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/40 backdrop-blur-sm p-4"
                    >
                        <motion.div
                            initial={{ scale: 0.9, y: 20 }}
                            animate={{ scale: 1, y: 0 }}
                            exit={{ scale: 0.9, y: 20 }}
                            onClick={(e) => e.stopPropagation()}
                            className="w-full max-w-lg bg-white rounded-3xl shadow-2xl overflow-hidden border border-slate-100"
                        >
                            <div className="bg-gradient-to-r from-blue-600 to-emerald-500 p-6 flex justify-between items-start">
                                <h2 className="text-xl font-black text-white pr-8 leading-tight">
                                    {selectedTerm.term}
                                </h2>
                                <button
                                    onClick={() => setSelectedTerm(null)}
                                    className="text-white/80 hover:text-white transition-colors bg-white/10 hover:bg-white/20 rounded-full p-1"
                                >
                                    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" /></svg>
                                </button>
                            </div>
                            <div className="p-8">
                                <p className="text-lg text-slate-600 leading-relaxed font-medium">
                                    {selectedTerm.definition}
                                </p>
                            </div>
                            <div className="px-8 pb-8 flex justify-end">
                                <button
                                    onClick={() => setSelectedTerm(null)}
                                    className="px-6 py-2 bg-slate-100 text-slate-600 rounded-xl font-bold text-sm hover:bg-slate-200 transition-colors"
                                >
                                    Cerrar
                                </button>
                            </div>
                        </motion.div>
                    </motion.div>
                )}
            </AnimatePresence>

            {/* Instructions */}
            <div className="absolute top-24 right-6 w-64 bg-white/80 backdrop-blur border border-slate-200 rounded-xl p-4 shadow-lg text-xs hidden lg:block">
                <h4 className="font-bold text-slate-800 uppercase mb-2 flex items-center gap-2">
                    <svg className="w-4 h-4 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>
                    Guía Rápida
                </h4>
                <ul className="space-y-1.5 text-slate-500">
                    <li className="flex gap-2">
                        <span className="text-blue-500 font-bold">•</span>
                        <span>Click en cualquier término para ver su definición.</span>
                    </li>
                    <li className="flex gap-2">
                        <span className="text-blue-500 font-bold">•</span>
                        <span>Usa los controles inferiores para hacer zoom.</span>
                    </li>
                    <li className="flex gap-2">
                        <span className="text-blue-500 font-bold">•</span>
                        <span>Arrastra y desplázate para navegar por el mapa.</span>
                    </li>
                </ul>
            </div>
        </div>
    );
};

const ControlButton = ({ onClick, icon, label }) => (
    <button
        onClick={onClick}
        className="p-2 hover:bg-slate-100 text-slate-600 rounded-xl transition-colors active:scale-95"
    >
        {icon || <span className="text-xs font-bold px-1">{label}</span>}
    </button>
);

export default ConceptMap;
