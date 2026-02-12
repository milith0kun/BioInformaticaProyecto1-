import React, { useEffect, useRef, useState, useMemo } from 'react';
import * as d3 from 'd3';
import { motion, AnimatePresence } from 'framer-motion';

// --- DATA ---
const graphData = {
    "nodes": [
        { "id": "Celula", "label": "1. Célula", "group": "fundamentales", "concepto": "Unidad estructural y funcional básica de los seres vivos." },
        { "id": "TeoriaCelular", "label": "2. Teoría celular", "group": "fundamentales", "concepto": "Postula que la célula es la unidad de vida, todas provienen de preexistentes y el organismo es un conjunto de ellas." },
        { "id": "CicloVida", "label": "3. Ciclo de vida celular", "group": "fundamentales", "concepto": "Serie de eventos que conducen al crecimiento y división de una célula (G1, S, G2, M)." },
        { "id": "Eucariota", "label": "22. Célula eucariota", "group": "fundamentales", "concepto": "Célula con núcleo definido y orgánulos membranosos (ej. células animales y vegetales)." },
        { "id": "Procariota", "label": "23. Célula procariota", "group": "fundamentales", "concepto": "Célula sin núcleo verdadero, el ADN está libre en el citoplasma (ej. bacterias)." },
        { "id": "ADN", "label": "5. ADN", "group": "estructurales", "concepto": "Ácido desoxirribonucleico; molécula que contiene la información genética hereditaria." },
        { "id": "ARN", "label": "6. ARN", "group": "estructurales", "concepto": "Ácido ribonucleico; molécula que transmite el código genético y participa en la síntesis de proteínas." },
        { "id": "Proteina", "label": "7. Proteína", "group": "estructurales", "concepto": "Macromolécula formada por aminoácidos que ejecuta la mayoría de las funciones celulares." },
        { "id": "Nucleotido", "label": "16. Nucleótido", "group": "estructurales", "concepto": "Monómero constituyente de los ácidos nucleicos, formado por un azúcar, un fosfato y una base nitrogenada." },
        { "id": "Bases", "label": "17. Bases nitrogenadas", "group": "estructurales", "concepto": "Componentes del ADN y ARN (Adenina, Guanina, Citosina, Timina/Uracilo) que codifican la información." },
        { "id": "Chargaff", "label": "18. Regla de Chargaff", "group": "estructurales", "concepto": "En el ADN, la cantidad de Adenina es igual a Timina, y la de Citosina igual a Guanina (A=T, C=G)." },
        { "id": "DobleHelice", "label": "19. Doble hélice", "group": "estructurales", "concepto": "Estructura tridimensional del ADN, formada por dos cadenas enrolladas en espiral." },
        { "id": "Complementariedad", "label": "20. Complementariedad", "group": "estructurales", "concepto": "Propiedad por la que las bases nitrogenadas se emparejan específicamente (A con T/U, C con G)." },
        { "id": "Cromosoma", "label": "9. Cromosoma", "group": "estructurales", "concepto": "Estructura organizada de ADN y proteínas que se encuentra en el núcleo de las células eucariotas." },
        { "id": "Gen", "label": "10. Gen", "group": "estructurales", "concepto": "Segmento de ADN que contiene la información para sintetizar una proteína o ARN funcional." },
        { "id": "Exon", "label": "24. Exón", "group": "estructurales", "concepto": "Región de un gen que se conserva en el ARN maduro y codifica para la proteína." },
        { "id": "Intron", "label": "25. Intrón", "group": "estructurales", "concepto": "Región no codificante de un gen que es eliminada durante el procesamiento del ARN (splicing)." },
        { "id": "HerenciaMendeliana", "label": "11. Herencia mendeliana", "group": "estructurales", "concepto": "Patrones de transmisión de rasgos de padres a hijos descritos por Gregor Mendel." },
        { "id": "Mutacion", "label": "12. Mutación", "group": "estructurales", "concepto": "Cambio permanente en la secuencia del ADN que puede alterar la información genética." },
        { "id": "Ligamiento", "label": "13. Ligamiento genético", "group": "estructurales", "concepto": "Tendencia de los genes cercanos en el mismo cromosoma a heredarse juntos." },
        { "id": "MapaGenetico", "label": "14. Mapa genético", "group": "estructurales", "concepto": "Representación lineal de la posición relativa de los genes en un cromosoma basada en la frecuencia de recombinación." },
        { "id": "UnGenUnaProteina", "label": "15. Un gen-una proteína", "group": "estructurales", "concepto": "Hipótesis que establece que cada gen codifica para una única proteína (o polipéptido)." },
        { "id": "Replicacion", "label": "21. Replicación del ADN", "group": "procesos", "concepto": "Proceso de duplicación del ADN para que las células hijas reciban una copia idéntica." },
        { "id": "Transcripcion", "label": "26. Transcripción", "group": "procesos", "concepto": "Síntesis de ARN a partir de una plantilla de ADN." },
        { "id": "Traduccion", "label": "30. Traducción", "group": "traduccion", "concepto": "Proceso en el que el ARNm es decodificado por el ribosoma para sintetizar una proteína." },
        { "id": "DogmaCentral", "label": "41. Dogma central", "group": "procesos", "concepto": "Flujo unidireccional de la información biológica: ADN -> ARN -> Proteína." },
        { "id": "ARNPolimerasa", "label": "27. ARN polimerasa", "group": "procesos", "concepto": "Enzima responsable de sintetizar ARN usando una hebra de ADN como molde." },
        { "id": "ARNm", "label": "28. ARN mensajero", "group": "procesos", "concepto": "Tipo de ARN que lleva la información genética del ADN al ribosoma para la síntesis de proteínas." },
        { "id": "Splicing", "label": "29. Splicing", "group": "procesos", "concepto": "Proceso de eliminación de intrones y unión de exones en el ARN pre-mensajero." },
        { "id": "Ribosoma", "label": "31. Ribosoma", "group": "traduccion", "concepto": "Complejo molecular (ARNr + proteínas) donde se lleva a cabo la traducción." },
        { "id": "Aminoacido", "label": "32. Aminoácido", "group": "traduccion", "concepto": "Monómero constituyente de las proteínas, unidos por enlaces peptídicos." },
        { "id": "Polipeptido", "label": "33. Polipéptido", "group": "traduccion", "concepto": "Cadena lineal de aminoácidos que se pliega para formar una proteína funcional." },
        { "id": "Codon", "label": "34. Codón", "group": "traduccion", "concepto": "Secuencia de tres nucleótidos en el ARNm que especifica un aminoácido o una señal de parada." },
        { "id": "CodigoGenetico", "label": "35. Código genético", "group": "traduccion", "concepto": "Conjunto de reglas que define cómo los codones se traducen en aminoácidos." },
        { "id": "Degeneracion", "label": "36. Degeneración del código", "group": "traduccion", "concepto": "Propiedad del código genético donde varios codones diferentes codifican el mismo aminoácido." },
        { "id": "CodonInicio", "label": "37. Codón de inicio", "group": "traduccion", "concepto": "Codón (AUG) que marca el sitio donde comienza la traducción de una proteína." },
        { "id": "CodonTerminacion", "label": "38. Codón de terminación", "group": "traduccion", "concepto": "Codones (UAA, UAG, UGA) que señalan el final de la síntesis proteica." },
        { "id": "ARNt", "label": "39. ARN de transferencia", "group": "traduccion", "concepto": "Tipo de ARN que transporta aminoácidos específicos al ribosoma durante la traducción." },
        { "id": "Anticodon", "label": "40. Anticodón", "group": "traduccion", "concepto": "Secuencia de tres bases en el ARNt que se aparean con el codón complementario del ARNm." },
        { "id": "ViaMetabolica", "label": "4. Vía metabólica", "group": "fundamentales", "concepto": "Serie de reacciones químicas enzimáticas consecutivas que transforman un sustrato en un producto." },
        { "id": "AlfabetoMolecular", "label": "8. Alfabeto molecular", "group": "estructurales", "concepto": "Conjunto limitado de símbolos (letras) que representan los bloques de construcción de la vida (A, T, C, G)." },
        { "id": "PCR", "label": "42. PCR", "group": "tecnicas", "concepto": "Reacción en Cadena de la Polimerasa; técnica para amplificar millones de copias de un fragmento de ADN." },
        { "id": "Cebador", "label": "43. Cebador/Primer", "group": "tecnicas", "concepto": "Oligonucleótido corto que sirve como punto de inicio para la síntesis de ADN en la PCR." },
        { "id": "Desnaturalizacion", "label": "44. Desnaturalización", "group": "tecnicas", "concepto": "Separación de las dos hebras de ADN mediante calor para exponer las secuencias." },
        { "id": "Alineamiento", "label": "45. Alineamiento/Priming", "group": "tecnicas", "concepto": "Unión específica de los cebadores a las secuencias complementarias del ADN molde." },
        { "id": "Extension", "label": "46. Extensión", "group": "tecnicas", "concepto": "Síntesis de nueva hebra de ADN por la ADN polimerasa a partir del cebador." },
        { "id": "ClonacionMolecular", "label": "47. Clonación molecular", "group": "tecnicas", "concepto": "Técnica para insertar un fragmento de ADN en un vector para su replicación en un huésped." },
        { "id": "Vector", "label": "48. Vector de clonación", "group": "tecnicas", "concepto": "Molécula de ADN (usualmente un plásmido) usada para transportar un gen foráneo." },
        { "id": "BibliotecaClones", "label": "49. Biblioteca de clones", "group": "tecnicas", "concepto": "Colección de células que portan diferentes fragmentos de ADN, representando un genoma completo." },
        { "id": "EnzimaRestriccion", "label": "50. Enzima de restricción", "group": "tecnicas", "concepto": "Proteína que corta el ADN en secuencias específicas (tijeras moleculares)." },
        { "id": "SitioReconocimiento", "label": "51. Sitio de reconocimiento", "group": "tecnicas", "concepto": "Secuencia de ADN específica donde una enzima de restricción corta." },
        { "id": "ExtremosRomos", "label": "52. Extremos romos", "group": "tecnicas", "concepto": "Cortes de ADN que terminan al mismo nivel, sin bases sobresalientes." },
        { "id": "ExtremosPegajosos", "label": "53. Extremos pegajosos", "group": "tecnicas", "concepto": "Cortes de ADN dejan cadenas simples sobresalientes que facilitan la unión con otros fragmentos." },
        { "id": "Hibridacion", "label": "54. Hibridación", "group": "tecnicas", "concepto": "Apareamiento de hebras de ADN o ARN complementarias para formar una molécula híbrida." },
        { "id": "Ligacion", "label": "55. Ligación", "group": "tecnicas", "concepto": "Unión enzimática de dos fragmentos de ADN mediante enlaces fosfodiéster." },
        { "id": "Electroforesis", "label": "56. Electroforesis en gel", "group": "tecnicas", "concepto": "Técnica para separar moléculas (ADN, proteínas) según su tamaño y carga mediante un campo eléctrico." },
        { "id": "Sonda", "label": "57. Sonda", "group": "tecnicas", "concepto": "Fragmento de ADN o ARN marcado (radioactivo o fluorescente) usado para detectar secuencias complementarias." },
        { "id": "Microarreglo", "label": "58. Microarreglo", "group": "tecnicas", "concepto": "Placa con miles de sondas fijas usada para medir la expresión de muchos genes a la vez." },
        { "id": "VariacionGenetica", "label": "59. Variación genética", "group": "genomica", "concepto": "Diferencias en las secuencias de ADN entre individuos de una misma especie." },
        { "id": "GenomaReferencia", "label": "60. Genoma de referencia", "group": "genomica", "concepto": "Secuencia genómica representativa de una especie usada como estándar para comparaciones." },
        { "id": "Conservacion", "label": "61. Conservación genética", "group": "genomica", "concepto": "Preservación de secuencias de ADN similares entre diferentes especies a lo largo de la evolución." },
        { "id": "Evolucion", "label": "62. Evolución", "group": "genomica", "concepto": "Cambio en las características heredables de las poblaciones biológicas a lo largo de las generaciones." },
        { "id": "SeleccionNatural", "label": "63. Selección natural", "group": "genomica", "concepto": "Proceso por el cual los organismos mejor adaptados tienden a sobrevivir y reproducirse más." },
        { "id": "Adaptacion", "label": "64. Adaptación", "group": "genomica", "concepto": "Rasgo heredable que mejora la supervivencia y reproducción de un organismo en su ambiente." },
        { "id": "Especiacion", "label": "65. Especiación", "group": "genomica", "concepto": "Formación de una nueva especie biológica a partir de una existente." },
        { "id": "GenomicaComparativa", "label": "66. Genómica comparativa", "group": "genomica", "concepto": "Estudio de las diferencias y similitudes en los genomas de diferentes especies." },
        { "id": "AlineamientoSecuencias", "label": "67. Alineamiento de secuencias", "group": "genomica", "concepto": "Disposición de dos o más secuencias de ADN/Proteínas para identificar regiones de similitud." },
        { "id": "BLAST", "label": "68. BLAST", "group": "genomica", "concepto": "Herramienta bioinformática para encontrar regiones de similitud local entre secuencias biológicas." },
        { "id": "Bioinformatica", "label": "69. Bioinformática", "group": "genomica", "concepto": "Disciplina que utiliza software y algoritmos para analizar y interpretar datos biológicos complejos." }
    ],
    "links": [
        { "source": "TeoriaCelular", "target": "Celula", "label": "postula que" },
        { "source": "Celula", "target": "CicloVida", "label": "pasa por" },
        { "source": "Celula", "target": "Eucariota", "label": "se clasifica en" },
        { "source": "Celula", "target": "Procariota", "label": "se clasifica en" },
        { "source": "Celula", "target": "ViaMetabolica", "label": "contiene" },
        { "source": "ADN", "target": "Nucleotido", "label": "está formado por" },
        { "source": "ARN", "target": "Nucleotido", "label": "está formado por" },
        { "source": "Nucleotido", "target": "Bases", "label": "contiene" },
        { "source": "Bases", "target": "Chargaff", "label": "cumplen" },
        { "source": "Bases", "target": "Complementariedad", "label": "permiten" },
        { "source": "ADN", "target": "DobleHelice", "label": "forma" },
        { "source": "DobleHelice", "target": "Complementariedad", "label": "se basa en" },
        { "source": "Complementariedad", "target": "Chargaff", "label": "explica" },
        { "source": "ADN", "target": "Cromosoma", "label": "se organiza en" },
        { "source": "Cromosoma", "target": "Gen", "label": "contiene" },
        { "source": "Gen", "target": "Exon", "label": "tiene regiones" },
        { "source": "Gen", "target": "Intron", "label": "tiene regiones" },
        { "source": "Eucariota", "target": "Exon", "label": "posee" },
        { "source": "Eucariota", "target": "Intron", "label": "posee" },
        { "source": "Procariota", "target": "Intron", "label": "carece de", "dashed": true },
        { "source": "Gen", "target": "HerenciaMendeliana", "label": "se transmite por" },
        { "source": "Cromosoma", "target": "Ligamiento", "label": "muestra" },
        { "source": "Ligamiento", "target": "MapaGenetico", "label": "permite construir" },
        { "source": "ADN", "target": "Mutacion", "label": "puede sufrir" },
        { "source": "Gen", "target": "UnGenUnaProteina", "label": "según hipótesis produce" },
        { "source": "UnGenUnaProteina", "target": "Proteina", "label": "codifica" },
        { "source": "DogmaCentral", "target": "ADN", "label": "describe flujo desde" },
        { "source": "DogmaCentral", "target": "ARN", "label": "describe flujo hacia" },
        { "source": "DogmaCentral", "target": "Proteina", "label": "describe flujo hasta" },
        { "source": "ADN", "target": "Replicacion", "label": "se duplica por" },
        { "source": "Replicacion", "target": "Complementariedad", "label": "usa" },
        { "source": "ADN", "target": "Transcripcion", "label": "se copia a" },
        { "source": "Transcripcion", "target": "ARN", "label": "produce" },
        { "source": "ARN", "target": "Traduccion", "label": "se decodifica en" },
        { "source": "Traduccion", "target": "Proteina", "label": "sintetiza" },
        { "source": "Transcripcion", "target": "ARNPolimerasa", "label": "es catalizada por" },
        { "source": "Transcripcion", "target": "ARNm", "label": "genera" },
        { "source": "ARNPolimerasa", "target": "ARNm", "label": "sintetiza" },
        { "source": "ARNm", "target": "Splicing", "label": "es procesado por" },
        { "source": "Splicing", "target": "Exon", "label": "une" },
        { "source": "Splicing", "target": "Intron", "label": "elimina" },
        { "source": "Eucariota", "target": "Splicing", "label": "requiere" },
        { "source": "ARNm", "target": "Traduccion", "label": "es leído en" },
        { "source": "Traduccion", "target": "Ribosoma", "label": "ocurre en" },
        { "source": "Traduccion", "target": "Aminoacido", "label": "incorpora" },
        { "source": "Ribosoma", "target": "Polipeptido", "label": "ensambla" },
        { "source": "Aminoacido", "target": "Polipeptido", "label": "forma" },
        { "source": "Polipeptido", "target": "Proteina", "label": "se pliega en" },
        { "source": "ARNm", "target": "Codon", "label": "se lee en" },
        { "source": "Codon", "target": "CodigoGenetico", "label": "está definido por" },
        { "source": "CodigoGenetico", "target": "Aminoacido", "label": "especifica" },
        { "source": "CodigoGenetico", "target": "Degeneracion", "label": "presenta" },
        { "source": "CodigoGenetico", "target": "CodonInicio", "label": "incluye" },
        { "source": "CodigoGenetico", "target": "CodonTerminacion", "label": "incluye" },
        { "source": "Traduccion", "target": "ARNt", "label": "requiere" },
        { "source": "ARNt", "target": "Anticodon", "label": "porta" },
        { "source": "ARNt", "target": "Aminoacido", "label": "transporta" },
        { "source": "Anticodon", "target": "Codon", "label": "reconoce" },
        { "source": "Complementariedad", "target": "Anticodon", "label": "permite" },
        { "source": "ADN", "target": "AlfabetoMolecular", "label": "se representa como" },
        { "source": "ARN", "target": "AlfabetoMolecular", "label": "se representa como" },
        { "source": "Proteina", "target": "AlfabetoMolecular", "label": "se representa como" },
        { "source": "PCR", "target": "ADN", "label": "amplifica" },
        { "source": "PCR", "target": "Cebador", "label": "usa" },
        { "source": "PCR", "target": "Desnaturalizacion", "label": "incluye fase" },
        { "source": "PCR", "target": "Alineamiento", "label": "incluye fase" },
        { "source": "PCR", "target": "Extension", "label": "incluye fase" },
        { "source": "Desnaturalizacion", "target": "DobleHelice", "label": "separa" },
        { "source": "Alineamiento", "target": "Cebador", "label": "posiciona" },
        { "source": "Alineamiento", "target": "Complementariedad", "label": "usa" },
        { "source": "Extension", "target": "Replicacion", "label": "replica mediante" },
        { "source": "Cebador", "target": "Hibridacion", "label": "se une por" },
        { "source": "ClonacionMolecular", "target": "ADN", "label": "inserta" },
        { "source": "ClonacionMolecular", "target": "Vector", "label": "utiliza" },
        { "source": "Vector", "target": "BibliotecaClones", "label": "genera" },
        { "source": "ADN", "target": "EnzimaRestriccion", "label": "se corta con" },
        { "source": "EnzimaRestriccion", "target": "SitioReconocimiento", "label": "reconoce" },
        { "source": "EnzimaRestriccion", "target": "ExtremosRomos", "label": "puede generar" },
        { "source": "EnzimaRestriccion", "target": "ExtremosPegajosos", "label": "puede generar" },
        { "source": "ExtremosPegajosos", "target": "Hibridacion", "label": "facilitan" },
        { "source": "Hibridacion", "target": "Complementariedad", "label": "se basa en" },
        { "source": "ClonacionMolecular", "target": "Ligacion", "label": "requiere" },
        { "source": "Ligacion", "target": "ADN", "label": "sella" },
        { "source": "Electroforesis", "target": "ADN", "label": "separa fragmentos de" },
        { "source": "Sonda", "target": "Hibridacion", "label": "detecta por" },
        { "source": "Sonda", "target": "Complementariedad", "label": "usa" },
        { "source": "Microarreglo", "target": "Sonda", "label": "contiene múltiples" },
        { "source": "Microarreglo", "target": "ARNm", "label": "detecta" },
        { "source": "Microarreglo", "target": "Hibridacion", "label": "funciona por" },
        { "source": "ADN", "target": "VariacionGenetica", "label": "presenta" },
        { "source": "VariacionGenetica", "target": "GenomaReferencia", "label": "se compara con" },
        { "source": "VariacionGenetica", "target": "Evolucion", "label": "impulsa" },
        { "source": "Mutacion", "target": "VariacionGenetica", "label": "genera" },
        { "source": "Gen", "target": "Conservacion", "label": "muestra" },
        { "source": "Evolucion", "target": "Conservacion", "label": "preserva por" },
        { "source": "Evolucion", "target": "SeleccionNatural", "label": "ocurre por" },
        { "source": "SeleccionNatural", "target": "Adaptacion", "label": "produce" },
        { "source": "Evolucion", "target": "Especiacion", "label": "puede llevar a" },
        { "source": "Bioinformatica", "target": "GenomicaComparativa", "label": "permite" },
        { "source": "Bioinformatica", "target": "AlineamientoSecuencias", "label": "utiliza" },
        { "source": "Bioinformatica", "target": "BLAST", "label": "implementa" },
        { "source": "GenomicaComparativa", "target": "AlineamientoSecuencias", "label": "requiere" },
        { "source": "AlineamientoSecuencias", "target": "BLAST", "label": "ejecuta" },
        { "source": "AlineamientoSecuencias", "target": "Conservacion", "label": "identifica" },
        { "source": "BLAST", "target": "Conservacion", "label": "detecta" },
        { "source": "Bioinformatica", "target": "AlfabetoMolecular", "label": "procesa" },
        { "source": "AlineamientoSecuencias", "target": "ADN", "label": "compara" },
        { "source": "AlineamientoSecuencias", "target": "Proteina", "label": "compara" },
        { "source": "ViaMetabolica", "target": "Proteina", "label": "es catalizada por" },
        { "source": "Gen", "target": "Transcripcion", "label": "es expresado por" },
        { "source": "Proteina", "target": "ViaMetabolica", "label": "ejecuta" }
    ]
};

const NetworkMap = () => {
    const svgRef = useRef(null);
    const containerRef = useRef(null);
    const [searchTerm, setSearchTerm] = useState('');
    const [searchResults, setSearchResults] = useState([]);
    const [selectedTooltip, setSelectedTooltip] = useState(null);
    const [zoomLevel, setZoomLevel] = useState(1);

    const [simulation, setSimulation] = useState(null);
    const [zoomTransform, setZoomTransform] = useState(d3.zoomIdentity);
    const nodesRef = useRef([]);

    // --- COLORS ---
    const groupColors = {
        fundamentales: { fill: "#e1f5ff", stroke: "#0277bd" },
        estructurales: { fill: "#fff3e0", stroke: "#ef6c00" },
        procesos: { fill: "#f3e5f5", stroke: "#6a1b9a" },
        tecnicas: { fill: "#e8f5e9", stroke: "#2e7d32" },
        genomica: { fill: "#fce4ec", stroke: "#c2185b" },
        traduccion: { fill: "#fff9c4", stroke: "#f9a825" }
    };
    const getColor = (group) => groupColors[group] || { fill: "#e3f2fd", stroke: "#1976d2" };

    useEffect(() => {
        if (!svgRef.current) return;

        const width = containerRef.current.clientWidth;
        const height = containerRef.current.clientHeight;

        const svg = d3.select(svgRef.current)
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", [0, 0, width, height]);

        svg.selectAll("*").remove(); // Clear previous render

        // --- DEFS ---
        const defs = svg.append("defs");
        const filter = defs.append("filter")
            .attr("id", "glow")
            .attr("x", "-50%")
            .attr("y", "-50%")
            .attr("width", "200%")
            .attr("height", "200%");
        filter.append("feGaussianBlur")
            .attr("stdDeviation", "3")
            .attr("result", "coloredBlur");
        const feMerge = filter.append("feMerge");
        feMerge.append("feMergeNode").attr("in", "coloredBlur");
        feMerge.append("feMergeNode").attr("in", "SourceGraphic");

        defs.append("marker")
            .attr("id", "arrow-normal")
            .attr("viewBox", "0 -5 10 10")
            .attr("refX", 20)
            .attr("refY", 0)
            .attr("markerWidth", 6)
            .attr("markerHeight", 6)
            .attr("orient", "auto")
            .append("path")
            .attr("d", "M0,-5L10,0L0,5")
            .attr("fill", "#999");

        // --- GROUPS ---
        const g = svg.append("g");
        const linkHitAreaGroup = g.append("g").attr("class", "link-hit-areas");
        const linkGroup = g.append("g").attr("class", "links");
        const linkLabelGroup = g.append("g").attr("class", "link-labels");
        const nodeGroup = g.append("g").attr("class", "nodes");
        const labelGroup = g.append("g").attr("class", "labels");

        // --- SIMULATION ---
        // Clone data to avoid mutation issues in React strict mode
        const nodes = graphData.nodes.map(d => ({ ...d }));
        const links = graphData.links.map(d => ({ ...d }));
        nodesRef.current = nodes;

        // Calculate radius
        const nodeConnections = {};
        nodes.forEach(n => nodeConnections[n.id] = 0);
        links.forEach(l => {
            const sourceId = typeof l.source === 'object' ? l.source.id : l.source;
            const targetId = typeof l.target === 'object' ? l.target.id : l.target;
            nodeConnections[sourceId] = (nodeConnections[sourceId] || 0) + 1;
            nodeConnections[targetId] = (nodeConnections[targetId] || 0) + 1;
        });
        const maxConnections = Math.max(...Object.values(nodeConnections), 1);
        const getNodeRadius = (nodeId) => {
            const connections = nodeConnections[nodeId] || 0;
            return 8 + (connections / maxConnections) * 8;
        };

        const sim = d3.forceSimulation(nodes)
            .force("link", d3.forceLink(links).id(d => d.id).distance(120))
            .force("charge", d3.forceManyBody().strength(-500))
            .force("center", d3.forceCenter(width / 2, height / 2))
            .force("collide", d3.forceCollide().radius(d => getNodeRadius(d.id) + 15));

        // --- DRAWING ---
        const link = linkGroup.selectAll("line")
            .data(links)
            .join("line")
            .attr("stroke", "#999")
            .attr("stroke-width", 2)
            .attr("stroke-opacity", 0.5)
            .attr("marker-end", "url(#arrow-normal)");

        const linkLabel = linkLabelGroup.selectAll("text")
            .data(links)
            .join("text")
            .attr("font-size", "9px")
            .attr("fill", "#666")
            .attr("text-anchor", "middle")
            .text(d => d.label);

        const node = nodeGroup.selectAll("circle")
            .data(nodes)
            .join("circle")
            .attr("r", d => getNodeRadius(d.id))
            .attr("fill", d => getColor(d.group).fill)
            .attr("stroke", d => getColor(d.group).stroke)
            .attr("stroke-width", 2)
            .attr("cursor", "pointer")
            .call(d3.drag()
                .on("start", (event, d) => {
                    if (!event.active) sim.alphaTarget(0.3).restart();
                    d.fx = d.x;
                    d.fy = d.y;
                })
                .on("drag", (event, d) => {
                    d.fx = event.x;
                    d.fy = event.y;
                })
                .on("end", (event, d) => {
                    if (!event.active) sim.alphaTarget(0);
                    d.fx = null;
                    d.fy = null;
                })
            );

        const label = labelGroup.selectAll("text")
            .data(nodes)
            .join("text")
            .attr("font-size", "10px")
            .attr("font-weight", "600")
            .attr("dx", d => getNodeRadius(d.id) + 4)
            .attr("dy", 4)
            .text(d => d.label || d.id)
            .style("pointer-events", "none")
            .style("text-shadow", "1px 1px 2px white");

        // --- UPDATES ---
        sim.on("tick", () => {
            link
                .attr("x1", d => d.source.x)
                .attr("y1", d => d.source.y)
                .attr("x2", d => d.target.x)
                .attr("y2", d => d.target.y);

            linkLabel
                .attr("x", d => (d.source.x + d.target.x) / 2)
                .attr("y", d => (d.source.y + d.target.y) / 2);

            node
                .attr("cx", d => d.x)
                .attr("cy", d => d.y);

            label
                .attr("x", d => d.x)
                .attr("y", d => d.y);
        });

        // --- ZOOM ---
        const zoomBehavior = d3.zoom()
            .scaleExtent([0.1, 4])
            .on("zoom", (event) => {
                g.attr("transform", event.transform);
                setZoomTransform(event.transform);
            });

        svg.call(zoomBehavior);

        // --- TOOLTIPS ---
        node.on("mouseover", (event, d) => {
            setSelectedTooltip({
                x: event.pageX,
                y: event.pageY,
                data: d
            });
        }).on("mouseout", () => {
            setSelectedTooltip(null);
        });

        // Save references
        setSimulation(sim);

        return () => {
            sim.stop();
        };
    }, []);

    // Focus mechanism
    const focusOnNode = (node) => {
        if (!svgRef.current || !node || !containerRef.current) return;

        const width = containerRef.current.clientWidth;
        const height = containerRef.current.clientHeight;
        const scale = 2.5; // Zoom in closer

        // Center the node: width/2 = x * k + tx  =>  tx = width/2 - x * k
        const tx = width / 2 - node.x * scale;
        const ty = height / 2 - node.y * scale;

        const svg = d3.select(svgRef.current);
        const transform = d3.zoomIdentity.translate(tx, ty).scale(scale);

        svg.transition()
            .duration(1000)
            .call(d3.zoom().on("zoom", (e) => {
                svg.select("g").attr("transform", e.transform);
                setZoomTransform(e.transform);
            }).transform, transform);

        // Highlight visual feedback
        setSelectedTooltip({
            x: width / 2 - 100, // Approximate center
            y: height / 2 - 100,
            data: node
        });
    };

    // Search functionality
    useEffect(() => {
        if (!searchTerm || searchTerm.length < 1) {
            setSearchResults([]);
            return;
        }
        const term = searchTerm.toLowerCase();
        const matches = graphData.nodes.filter(n =>
            (n.label || n.id).toLowerCase().includes(term) ||
            (n.concepto || "").toLowerCase().includes(term)
        ).sort((a, b) => {
            const aLabel = (a.label || "").toLowerCase();
            const bLabel = (b.label || "").toLowerCase();
            const aStarts = aLabel.startsWith(term);
            const bStarts = bLabel.startsWith(term);
            if (aStarts && !bStarts) return -1;
            if (!aStarts && bStarts) return 1;
            return 0;
        });
        setSearchResults(matches);
    }, [searchTerm]);

    const zoomIn = () => {
        const svg = d3.select(svgRef.current);
        svg.transition().call(d3.zoom().on("zoom", (e) => {
            svg.select("g").attr("transform", e.transform);
            setZoomTransform(e.transform);
        }).scaleBy, 1.2);
    };

    const zoomOut = () => {
        const svg = d3.select(svgRef.current);
        svg.transition().call(d3.zoom().on("zoom", (e) => {
            svg.select("g").attr("transform", e.transform);
            setZoomTransform(e.transform);
        }).scaleBy, 0.8);
    };

    const resetZoom = () => {
        const svg = d3.select(svgRef.current);
        svg.transition().call(d3.zoom().on("zoom", (e) => {
            svg.select("g").attr("transform", e.transform);
            setZoomTransform(e.transform);
        }).transform, d3.zoomIdentity);
    };

    // Auto-organize (re-heat simulation)
    const autoOrganize = () => {
        if (simulation) {
            simulation.alpha(1).restart();
        }
    };

    return (
        <div className="relative w-full h-[800px] bg-slate-50 overflow-hidden" ref={containerRef}>
            {/* Header Controls */}
            <div className="absolute top-0 left-0 right-0 z-10 p-4 flex justify-between items-start pointer-events-none">
                <div className="bg-white/90 backdrop-blur shadow-sm border border-slate-200 rounded-2xl p-4 pointer-events-auto">
                    <h1 className="text-xl font-black text-slate-800 uppercase tracking-tight">Red Genómica</h1>
                    <div className="mt-2 relative">
                        <input
                            type="text"
                            placeholder="Buscar N° o concepto..."
                            className="w-64 px-4 py-2 bg-slate-50 border border-slate-200 rounded-xl text-sm font-bold focus:outline-none focus:ring-2 focus:ring-blue-500"
                            value={searchTerm}
                            onChange={(e) => setSearchTerm(e.target.value)}
                        />
                        {searchResults.length > 0 && (
                            <div className="absolute top-full left-0 right-0 mt-2 bg-white border border-slate-200 rounded-xl shadow-xl max-h-60 overflow-y-auto z-20">
                                {searchResults.map((result) => (
                                    <div
                                        key={result.id}
                                        className="px-4 py-3 hover:bg-slate-50 cursor-pointer border-b border-slate-100 last:border-0"
                                        onClick={() => {
                                            const activeNode = nodesRef.current.find(n => n.id === result.id);
                                            if (activeNode) {
                                                focusOnNode(activeNode);
                                            }
                                            setSearchTerm(result.label);
                                            setSearchResults([]);
                                        }}
                                    >
                                        <div className="font-bold text-xs text-slate-800">{result.label}</div>
                                        <div className="text-[10px] text-slate-500 truncate">{result.concepto}</div>
                                    </div>
                                ))}
                            </div>
                        )}
                    </div>
                </div>

                <div className="flex gap-2 pointer-events-auto">
                    <button onClick={autoOrganize} className="px-4 py-2 bg-white text-slate-600 rounded-xl font-bold text-xs shadow-sm border border-slate-200 hover:bg-slate-50">
                        Auto-organizar
                    </button>
                    <button onClick={resetZoom} className="px-4 py-2 bg-blue-600 text-white rounded-xl font-bold text-xs shadow-lg shadow-blue-200 hover:bg-blue-700">
                        Centrar
                    </button>
                </div>
            </div>

            {/* SVG Canvas */}
            <svg ref={svgRef} className="w-full h-full cursor-move"></svg>

            {/* Tooltip Overlay */}
            <AnimatePresence>
                {selectedTooltip && (
                    <motion.div
                        initial={{ opacity: 0, scale: 0.9 }}
                        animate={{ opacity: 1, scale: 1 }}
                        exit={{ opacity: 0 }}
                        style={{ left: selectedTooltip.x + 20, top: selectedTooltip.y + 20 }}
                        className="fixed z-50 bg-white/95 backdrop-blur border-l-4 border-blue-500 shadow-2xl p-4 rounded-r-xl max-w-xs pointer-events-none"
                    >
                        <h3 className="font-bold text-slate-800 text-sm mb-1">{selectedTooltip.data.label}</h3>
                        <p className="text-xs text-slate-600 leading-relaxed font-medium">{selectedTooltip.data.concepto}</p>
                        <div className="mt-2 text-[10px] uppercase font-black text-slate-400 tracking-wider">
                            Grupo: {selectedTooltip.data.group}
                        </div>
                    </motion.div>
                )}
            </AnimatePresence>

            {/* Legend */}
            <div className="absolute bottom-6 left-6 bg-white/90 backdrop-blur p-4 rounded-2xl border border-slate-200 shadow-lg text-xs pointer-events-none">
                <h4 className="font-black text-slate-400 uppercase tracking-wider mb-2 text-[10px]">Categorías</h4>
                <div className="grid grid-cols-2 gap-x-4 gap-y-2">
                    {Object.entries(groupColors).map(([key, colors]) => (
                        <div key={key} className="flex items-center gap-2">
                            <div className="w-3 h-3 rounded-full border-2" style={{ backgroundColor: colors.fill, borderColor: colors.stroke }}></div>
                            <span className="capitalize text-slate-600 font-bold">{key}</span>
                        </div>
                    ))}
                </div>
            </div>

            {/* Zoom Controls */}
            <div className="absolute bottom-6 right-6 flex flex-col gap-2 pointer-events-auto">
                <button onClick={zoomIn} className="w-10 h-10 bg-white rounded-xl shadow-lg border border-slate-200 flex items-center justify-center text-slate-600 hover:bg-slate-50 font-bold text-lg">+</button>
                <button onClick={zoomOut} className="w-10 h-10 bg-white rounded-xl shadow-lg border border-slate-200 flex items-center justify-center text-slate-600 hover:bg-slate-50 font-bold text-lg">-</button>
            </div>
        </div>
    );
};

export default NetworkMap;
