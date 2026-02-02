"""
Export API Routes
Handles data export in various formats
"""
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import StreamingResponse, FileResponse
import json
import csv
import io
from typing import Optional
from datetime import datetime

router = APIRouter()


def get_analysis_cache():
    """Import and get analysis cache from analysis module"""
    from app.api.routes.analysis import _analysis_cache
    return _analysis_cache


@router.get("/json")
async def export_json():
    """
    Export all analysis results as JSON
    """
    cache = get_analysis_cache()
    
    if cache["genes"] is None or cache["codons"] is None:
        raise HTTPException(
            status_code=404, 
            detail="No hay resultados para exportar. Ejecute el análisis primero."
        )
    
    # Build export data
    export_data = {
        "metadata": {
            "exported_at": datetime.now().isoformat(),
            "organism": "Escherichia coli K-12 MG1655",
            "analysis_version": "1.0.0"
        },
        "codon_analysis": {
            "genome_length": cache["codons"].genome_length,
            "atg_count": cache["codons"].atg_count,
            "atg_density": cache["codons"].atg_density,
            "stop_codons": {
                codon: {"count": info.count, "percentage": info.percentage}
                for codon, info in cache["codons"].stop_codons.items()
            },
            "gene_comparison": {
                "annotated_genes": cache["codons"].gene_comparison.annotated_genes,
                "atg_found": cache["codons"].gene_comparison.atg_found,
                "difference": cache["codons"].gene_comparison.difference
            }
        },
        "gene_analysis": {
            "total_genes": cache["genes"].total_genes,
            "total_cds": cache["genes"].total_cds,
            "genome_length": cache["genes"].genome_length,
            "gc_content": cache["genes"].gc_content,
            "gene_density": cache["genes"].gene_density,
            "size_statistics": {
                "mean": cache["genes"].size_statistics.mean,
                "median": cache["genes"].size_statistics.median,
                "min": cache["genes"].size_statistics.min,
                "max": cache["genes"].size_statistics.max,
                "std": cache["genes"].size_statistics.std
            },
            "genes": [
                {
                    "locus_tag": g.locus_tag,
                    "start": g.start,
                    "end": g.end,
                    "length": g.length,
                    "strand": g.strand,
                    "product": g.product,
                    "gc_content": g.gc_content
                }
                for g in cache["genes"].genes
            ]
        }
    }
    
    if cache["validation"]:
        export_data["validation"] = {
            "items": [
                {
                    "metric": item.metric,
                    "calculated": item.calculated,
                    "reference": item.reference,
                    "deviation_percent": item.deviation_percent,
                    "status": item.status.value
                }
                for item in cache["validation"].items
            ],
            "overall_status": cache["validation"].overall_status.value
        }
    
    # Create streaming response
    json_str = json.dumps(export_data, indent=2, ensure_ascii=False)
    
    return StreamingResponse(
        io.BytesIO(json_str.encode('utf-8')),
        media_type="application/json",
        headers={
            "Content-Disposition": "attachment; filename=ecoli_analysis.json"
        }
    )


@router.get("/csv/{type}")
async def export_csv(type: str):
    """
    Export analysis results as CSV
    Types: genes, codons, validation, statistics
    """
    cache = get_analysis_cache()
    
    if type == "genes":
        if cache["genes"] is None:
            raise HTTPException(status_code=404, detail="No hay datos de genes para exportar")
        
        output = io.StringIO()
        writer = csv.writer(output)
        
        # Header
        writer.writerow([
            "locus_tag", "start", "end", "length", "strand", "product", "gc_content"
        ])
        
        # Data
        for g in cache["genes"].genes:
            writer.writerow([
                g.locus_tag, g.start, g.end, g.length, g.strand, g.product or "", g.gc_content
            ])
        
        filename = "genes.csv"
        
    elif type == "codons":
        if cache["codons"] is None:
            raise HTTPException(status_code=404, detail="No hay datos de codones para exportar")
        
        output = io.StringIO()
        writer = csv.writer(output)
        
        # Codon summary
        writer.writerow(["Métrica", "Valor"])
        writer.writerow(["Longitud del genoma", cache["codons"].genome_length])
        writer.writerow(["Conteo ATG", cache["codons"].atg_count])
        writer.writerow(["Densidad ATG (por kb)", cache["codons"].atg_density])
        writer.writerow([])
        
        writer.writerow(["Codón de parada", "Conteo", "Porcentaje"])
        for codon, info in cache["codons"].stop_codons.items():
            writer.writerow([codon, info.count, info.percentage])
        
        filename = "codons.csv"
        
    elif type == "validation":
        if cache["validation"] is None:
            raise HTTPException(status_code=404, detail="No hay datos de validación para exportar")
        
        output = io.StringIO()
        writer = csv.writer(output)
        
        writer.writerow(["Métrica", "Calculado", "Referencia", "Desviación (%)", "Estado"])
        for item in cache["validation"].items:
            writer.writerow([
                item.metric, item.calculated, item.reference, 
                item.deviation_percent, item.status.value
            ])
        
        filename = "validation.csv"
        
    elif type == "statistics":
        if cache["genes"] is None:
            raise HTTPException(status_code=404, detail="No hay estadísticas para exportar")
        
        output = io.StringIO()
        writer = csv.writer(output)
        
        writer.writerow(["Estadística", "Valor"])
        writer.writerow(["Total de genes", cache["genes"].total_genes])
        writer.writerow(["Total de CDS", cache["genes"].total_cds])
        writer.writerow(["Longitud del genoma", cache["genes"].genome_length])
        writer.writerow(["Contenido GC (%)", cache["genes"].gc_content])
        writer.writerow(["Densidad génica (genes/Mb)", cache["genes"].gene_density])
        writer.writerow(["Longitud media de genes", cache["genes"].size_statistics.mean])
        writer.writerow(["Longitud mediana de genes", cache["genes"].size_statistics.median])
        writer.writerow(["Gen más corto", cache["genes"].size_statistics.min])
        writer.writerow(["Gen más largo", cache["genes"].size_statistics.max])
        
        filename = "statistics.csv"
        
    else:
        raise HTTPException(status_code=400, detail=f"Tipo de exportación no válido: {type}")
    
    output.seek(0)
    
    return StreamingResponse(
        io.BytesIO(output.getvalue().encode('utf-8')),
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename={filename}"
        }
    )


@router.get("/pdf")
async def export_pdf():
    """
    Generate PDF report with analysis results
    """
    cache = get_analysis_cache()
    
    if cache["genes"] is None or cache["codons"] is None:
        raise HTTPException(
            status_code=404,
            detail="No hay resultados para generar el informe. Ejecute el análisis primero."
        )
    
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
        
        # Create PDF in memory
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter)
        styles = getSampleStyleSheet()
        story = []
        
        # Title
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=18,
            spaceAfter=30
        )
        story.append(Paragraph("Análisis Genómico de E. coli K-12 MG1655", title_style))
        story.append(Spacer(1, 12))
        
        # Date
        story.append(Paragraph(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}", styles['Normal']))
        story.append(Spacer(1, 20))
        
        # Genome Statistics
        story.append(Paragraph("Estadísticas del Genoma", styles['Heading2']))
        story.append(Spacer(1, 12))
        
        genome_data = [
            ["Métrica", "Valor"],
            ["Longitud del genoma", f"{cache['genes'].genome_length:,} bp"],
            ["Contenido GC", f"{cache['genes'].gc_content}%"],
            ["Total de genes", f"{cache['genes'].total_genes:,}"],
            ["Total de CDS", f"{cache['genes'].total_cds:,}"],
            ["Densidad génica", f"{cache['genes'].gene_density} genes/Mb"],
        ]
        
        genome_table = Table(genome_data, colWidths=[3*inch, 2*inch])
        genome_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        story.append(genome_table)
        story.append(Spacer(1, 20))
        
        # Codon Analysis
        story.append(Paragraph("Análisis de Codones", styles['Heading2']))
        story.append(Spacer(1, 12))
        
        codon_data = [
            ["Métrica", "Valor"],
            ["Codones ATG encontrados", f"{cache['codons'].atg_count:,}"],
            ["Densidad ATG", f"{cache['codons'].atg_density} por kb"],
        ]
        
        for codon, info in cache['codons'].stop_codons.items():
            codon_data.append([f"Codón {codon}", f"{info.count:,} ({info.percentage}%)"])
        
        codon_table = Table(codon_data, colWidths=[3*inch, 2*inch])
        codon_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightblue),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        story.append(codon_table)
        story.append(Spacer(1, 20))
        
        # Gene Size Statistics
        story.append(Paragraph("Estadísticas de Tamaño de Genes", styles['Heading2']))
        story.append(Spacer(1, 12))
        
        size_data = [
            ["Estadística", "Valor (bp)"],
            ["Media", f"{cache['genes'].size_statistics.mean:,.1f}"],
            ["Mediana", f"{cache['genes'].size_statistics.median:,.1f}"],
            ["Mínimo", f"{cache['genes'].size_statistics.min:,}"],
            ["Máximo", f"{cache['genes'].size_statistics.max:,}"],
            ["Desviación estándar", f"{cache['genes'].size_statistics.std:,.1f}"],
        ]
        
        size_table = Table(size_data, colWidths=[3*inch, 2*inch])
        size_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.darkgreen),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgreen),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        story.append(size_table)
        
        # Validation if available
        if cache["validation"]:
            story.append(Spacer(1, 20))
            story.append(Paragraph("Validación contra Valores de Referencia", styles['Heading2']))
            story.append(Spacer(1, 12))
            
            val_data = [["Métrica", "Calculado", "Referencia", "Desviación", "Estado"]]
            for item in cache["validation"].items:
                val_data.append([
                    item.metric,
                    f"{item.calculated:,.2f}",
                    f"{item.reference:,.2f}",
                    f"{item.deviation_percent:.2f}%",
                    item.status.value
                ])
            
            val_table = Table(val_data, colWidths=[1.5*inch, 1.2*inch, 1.2*inch, 1*inch, 0.8*inch])
            val_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.purple),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            story.append(val_table)
        
        # Build PDF
        doc.build(story)
        buffer.seek(0)
        
        return StreamingResponse(
            buffer,
            media_type="application/pdf",
            headers={
                "Content-Disposition": "attachment; filename=ecoli_analysis_report.pdf"
            }
        )
        
    except ImportError:
        raise HTTPException(
            status_code=500,
            detail="ReportLab no está instalado. Instale con: pip install reportlab"
        )
