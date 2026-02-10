/**
 * API Service for E. coli Genomic Analysis
 */
import axios from 'axios'

const API_BASE_URL = '/api'

// Create axios instance
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  timeout: 120000, // 2 minutes for large analyses
  headers: {
    'Content-Type': 'application/json',
  },
})

// Request interceptor
axiosInstance.interceptors.request.use(
  (config) => {
    console.log(`📤 ${config.method?.toUpperCase()} ${config.url}`)
    return config
  },
  (error) => {
    return Promise.reject(error)
  }
)

// Response interceptor
axiosInstance.interceptors.response.use(
  (response) => {
    console.log(`📥 ${response.status} ${response.config.url}`)
    return response
  },
  (error) => {
    console.error(`❌ Error: ${error.message}`)
    return Promise.reject(error)
  }
)

export const api = {
  // File endpoints
  detectFiles: async () => {
    const response = await axiosInstance.get('/files/detect')
    return response.data
  },

  extractZip: async () => {
    const response = await axiosInstance.post('/files/extract')
    return response.data
  },

  getFileInfo: async (filename) => {
    const response = await axiosInstance.get(`/files/${filename}/info`)
    return response.data
  },

  getPrimaryFile: async () => {
    const response = await axiosInstance.get('/files/primary')
    return response.data
  },

  // Analysis endpoints
  analyzeCodens: async () => {
    const response = await axiosInstance.post('/analysis/codons')
    return response.data
  },

  analyzeGenes: async () => {
    const response = await axiosInstance.post('/analysis/genes')
    return response.data
  },

  validateResults: async () => {
    const response = await axiosInstance.get('/analysis/validate')
    return response.data
  },

  runCompleteAnalysis: async () => {
    const response = await axiosInstance.get('/analysis/complete')
    return response.data
  },

  getAnalysisStatus: async () => {
    const response = await axiosInstance.get('/analysis/status')
    return response.data
  },

  // Results endpoints
  getCodonResults: async () => {
    const response = await axiosInstance.get('/analysis/results/codons')
    return response.data
  },

  getGeneResults: async (page = 1, pageSize = 50, search = '') => {
    const params = new URLSearchParams({ page, page_size: pageSize })
    if (search) params.append('search', search)
    const response = await axiosInstance.get(`/analysis/results/genes?${params}`)
    return response.data
  },

  getStatistics: async () => {
    const response = await axiosInstance.get('/analysis/results/statistics')
    return response.data
  },

  clearCache: async () => {
    const response = await axiosInstance.post('/analysis/clear-cache')
    return response.data
  },

  // AI Validation endpoint
  validateWithAI: async (apiKey = null) => {
    const params = apiKey ? { api_key: apiKey } : {}
    const response = await axiosInstance.post('/analysis/ai-validation', null, { params })
    return response.data
  },

  // Export endpoints
  exportJson: () => {
    window.open(`${API_BASE_URL}/export/json`, '_blank')
  },

  exportCsv: (type) => {
    window.open(`${API_BASE_URL}/export/csv/${type}`, '_blank')
  },

  exportPdf: () => {
    window.open(`${API_BASE_URL}/export/pdf`, '_blank')
  },

  // Genome endpoints (NCBI Datasets API Integration)
  searchGenomes: async (query, limit = 10) => {
    const response = await axiosInstance.get(`/genome/search/${encodeURIComponent(query)}`, {
      params: { limit }
    })
    return response.data
  },

  getGenomeInfo: async (accession) => {
    const response = await axiosInstance.get(`/genome/info/${accession}`)
    return response.data
  },

  downloadGenome: async (options) => {
    const response = await axiosInstance.post('/genome/download', options)
    return response.data
  },

  getGenomeDownloadStatus: async (accession) => {
    const response = await axiosInstance.get(`/genome/download-status/${accession}`)
    return response.data
  },

  getDownloadedGenomes: async () => {
    const response = await axiosInstance.get('/genome/downloaded')
    return response.data
  },

  activateGenome: async (accession) => {
    const response = await axiosInstance.post(`/genome/activate/${accession}`)
    return response.data
  },

  deleteGenome: async (accession) => {
    const response = await axiosInstance.delete(`/genome/${accession}`)
    return response.data
  },

  getCurrentGenome: async () => {
    const response = await axiosInstance.get('/genome/current')
    return response.data
  },

  getPopularGenomes: async () => {
    const response = await axiosInstance.get('/genome/popular')
    return response.data
  },

  // ================== COMPARACIÓN GENÓMICA ==================

  // Obtener cepas de E. coli relacionadas
  getRelatedStrains: async (category = null) => {
    const params = category ? { category } : {}
    const response = await axiosInstance.get('/genome/related-strains', { params })
    return response.data
  },

  // Obtener grupos funcionales disponibles
  getFunctionalGroups: async () => {
    const response = await axiosInstance.get('/genome/functional-groups')
    return response.data
  },

  // Comparar genomas descargados
  compareGenomes: async (accessions = null) => {
    const response = await axiosInstance.post('/genome/compare', null, {
      params: accessions ? { accessions: accessions.join(',') } : {}
    })
    return response.data
  },

  // Obtener genes por tamaño (mayor o menor)
  getGenesBySize: async (order = 'largest', count = 10) => {
    const response = await axiosInstance.get(`/genome/genes/by-size/${order}`, {
      params: { count }
    })
    return response.data
  },

  // Obtener genes por grupo funcional
  getGenesByGroup: async (groupId) => {
    const response = await axiosInstance.get(`/genome/genes/by-group/${groupId}`)
    return response.data
  },

  // Búsqueda avanzada de genes
  searchGenesAdvanced: async (params) => {
    const response = await axiosInstance.get('/genome/genes/search-advanced', { params })
    return response.data
  },

  // Obtener resumen de grupos funcionales
  getGeneGroupsSummary: async () => {
    const response = await axiosInstance.get('/genome/genes/groups-summary')
    return response.data
  },

  // ================== NCBI ENHANCED ==================

  // Proteins
  getProteins: async (page = 1, pageSize = 50, search = '') => {
    const params = { page, page_size: pageSize }
    if (search) params.search = search
    const response = await axiosInstance.get('/ncbi/proteins', { params })
    return response.data
  },

  getProteinDetail: async (proteinId) => {
    const response = await axiosInstance.get(`/ncbi/protein/${proteinId}`)
    return response.data
  },

  // Gene Location
  getGeneAtPosition: async (position) => {
    const response = await axiosInstance.get(`/ncbi/gene-at-position/${position}`)
    return response.data
  },

  getGeneDetail: async (locusTag) => {
    const response = await axiosInstance.get(`/ncbi/gene-detail/${locusTag}`)
    return response.data
  },

  // Sequence Viewer
  getSequenceSegment: async (start = 0, end = 1000) => {
    const response = await axiosInstance.get('/ncbi/sequence', { params: { start, end } })
    return response.data
  },

  // Complete Codon Usage
  getCompleteCodonUsage: async () => {
    const response = await axiosInstance.get('/ncbi/codon-usage')
    return response.data
  },

  // Genome Map
  getGenomeMapData: async () => {
    const response = await axiosInstance.get('/ncbi/genome-map')
    return response.data
  },

  // Literature
  searchLiterature: async (query, maxResults = 8) => {
    const response = await axiosInstance.get('/ncbi/literature/search', {
      params: { query, max_results: maxResults }
    })
    return response.data
  },

  getGenBankReference: async (accession) => {
    const response = await axiosInstance.get(`/ncbi/genbank-reference/${accession}`)
    return response.data
  },

  // ================== AI CHAT ==================

  sendChatMessage: async (message, sessionId = 'default', includeGenomeContext = true, genomeAccession = null) => {
    const response = await axiosInstance.post('/chat/message', {
      message,
      session_id: sessionId,
      include_genome_context: includeGenomeContext,
      genome_accession: genomeAccession
    })
    return response.data
  },

  getChatHistory: async (sessionId = 'default') => {
    const response = await axiosInstance.get(`/chat/history/${sessionId}`)
    return response.data
  },

  clearChatHistory: async (sessionId = 'default') => {
    const response = await axiosInstance.delete(`/chat/history/${sessionId}`)
    return response.data
  },

  getChatSuggestions: async () => {
    const response = await axiosInstance.get('/chat/suggestions')
    return response.data
  },

  // ================== NCBI SEARCH ==================

  searchNCBIGene: async (query, organism = '', maxResults = 5) => {
    const response = await axiosInstance.get('/ncbi/search-gene', {
      params: { query, organism, max_results: maxResults }
    })
    return response.data
  },

  searchNCBINucleotide: async (query, maxResults = 3) => {
    const response = await axiosInstance.get('/ncbi/search-nucleotide', {
      params: { query, max_results: maxResults }
    })
    return response.data
  },

  // ================== CENTRAL DOGMA ==================

  getCentralDogma: async (locusTag) => {
    const response = await axiosInstance.get(`/ncbi/central-dogma/${locusTag}`)
    return response.data
  },

  // ================== GC SLIDING WINDOW ==================

  getGCSlidingWindow: async (windowSize = 5000, step = 1000) => {
    const response = await axiosInstance.get('/ncbi/gc-window', {
      params: { window_size: windowSize, step }
    })
    return response.data
  },

  // ================== BLAST SEARCH ==================

  submitBLAST: async (sequence, program = 'blastn', database = 'nt', maxHits = 10) => {
    const response = await axiosInstance.post('/ncbi/blast/submit', {
      sequence, program, database, max_hits: maxHits
    })
    return response.data
  },

  getBLASTResults: async (rid) => {
    const response = await axiosInstance.get(`/ncbi/blast/results/${rid}`)
    return response.data
  },

  // ================== RNA ANALYSIS ==================

  getRNAAnalysis: async () => {
    const response = await axiosInstance.get('/ncbi/rna-analysis')
    return response.data
  },

  // ================== GENOME SUMMARY ==================

  getGenomeSummary: async () => {
    const response = await axiosInstance.get('/ncbi/genome-summary')
    return response.data
  },

  // ================== FUNCTIONAL CATEGORIES ==================

  getFunctionalCategories: async () => {
    const response = await axiosInstance.get('/ncbi/functional-categories')
    return response.data
  },

  // ================== CAI ANALYSIS ==================

  getCAIAnalysis: async (topN = 50) => {
    const response = await axiosInstance.get('/ncbi/cai-analysis', {
      params: { top_n: topN }
    })
    return response.data
  },

  // ================== PHYLOGENETIC TREE ==================

  getPhylogeneticTree: async (genbankPaths = []) => {
    const response = await axiosInstance.post('/ncbi/phylogenetic-tree', {
      genbank_paths: genbankPaths
    })
    return response.data
  },
}

export default api
