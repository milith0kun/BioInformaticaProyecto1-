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
}

export default api
