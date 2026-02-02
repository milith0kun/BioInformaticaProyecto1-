# Integración con Google Gemini AI

## Estado Actual ✅

La integración con Google Gemini está **configurada y funcionando** con el modelo `gemini-2.0-flash`.

### Configuración

- **Modelo**: `gemini-2.0-flash` (versión más reciente, junio 2025)
- **API Version**: `v1`
- **Endpoint**: `https://generativelanguage.googleapis.com/v1/models/gemini-2.0-flash:generateContent`
- **API Key**: Configurada en `.env`

## Límites de la API Gratuita ⏱️

La API gratuita de Google Gemini tiene límites de tasa muy restrictivos:

- **Límite**: Aproximadamente 15-20 peticiones por minuto
- **Error común**: `429 Too Many Requests`
- **Solución**: Esperar 1-2 minutos entre validaciones

### Estrategias Implementadas

1. **Retry Automático**: El sistema reintenta automáticamente con backoff exponencial (2s, 4s, 8s)
2. **Mensajes Claros**: Cuando se excede el límite, se muestra un mensaje informativo
3. **Independencia**: Los resultados del análisis son válidos independientemente de la validación IA

## Uso

### Desde el Frontend

1. Ejecutar análisis completo
2. Ir a la pestaña "Validación IA"
3. Click en "Revalidar"
4. **Esperar** si aparece el mensaje de límite de tasa

### Desde la API

```bash
# 1. Ejecutar análisis
curl -X POST http://localhost:8000/api/analysis/codons
curl -X POST http://localhost:8000/api/analysis/genes

# 2. Esperar 15-30 segundos

# 3. Solicitar validación IA
curl -X POST http://localhost:8000/api/analysis/ai-validation
```

## Mensajes de Error

### ⏱️ Límite de peticiones excedido
**Causa**: Demasiadas peticiones en poco tiempo  
**Solución**: Esperar 1-2 minutos antes de reintentar

### 🔧 Modelo o endpoint no disponible
**Causa**: Modelo no existe en la API  
**Solución**: Verificar que el modelo `gemini-2.0-flash` esté disponible

### 🔑 API key inválida
**Causa**: API key revocada o sin permisos  
**Solución**: Generar nueva API key en [Google AI Studio](https://aistudio.google.com/)

## Renovar API Key

Si la API key expira o necesita renovarse:

1. Ir a [Google AI Studio](https://aistudio.google.com/)
2. Crear nueva API key
3. Actualizar `.env`:
   ```bash
   GEMINI_API_KEY=tu_nueva_key_aqui
   ```
4. Reiniciar backend

## Alternativas para Validación Frecuente

Si necesita validación IA frecuente, considere:

1. **Google AI Studio Pro**: Mayor cuota mensual
2. **Validación Manual**: Los resultados del análisis son científicamente correctos
3. **Batch Processing**: Ejecutar validación solo una vez por sesión

## Interpretación de Resultados

La validación IA proporciona:

- ✅ **Confianza**: 0-100% (qué tan seguros están los resultados)
- 📊 **Interpretación**: Análisis científico detallado
- ⚠️ **Discrepancias**: Problemas encontrados
- 💡 **Recomendaciones**: Sugerencias de mejora

## Nota Importante

> Los resultados del análisis genómico (codones, genes, validaciones) son científicamente correctos y basados en BioPython. La validación IA es **complementaria** y proporciona interpretación adicional, pero **no es requerida** para la validez de los resultados.
