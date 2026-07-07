# Comparative genomic of cellular complexity in Thiotrichales 🦠

## Objetive ⭐

Building on the comparative framework introduced by Volland et al. (2022), this study characterizes cellular complexity and biosynthetic potential across an expanded set of Thiotrichales genomes, evaluating whether *Candidatus Thiomargarita magnifica*  a genomic outlier with an exceptionally large genome  represents an extreme case within a broader trend across the order. 

Cellular complexity is assessed via genome size, CDS count, and cell division/elongation genes, while biosynthetic potential is quantified through the genomic extent of biosynthetic gene clusters (BGCs).

## Workflow / Pipeline 👣

The pipeline processes raw *Thiotrichales* genome assemblies (`.fasta`) into a standardized table used to reproduce comparative figures.

1. **Quality control — CheckM**
   Assesses genome completeness and contamination. Genomes are filtered using the thresholds: Completeness > 50%, Contamination < 30%. Also provides assembly size (Mbp).

2. **Genome annotation — Prokka**
   Annotates coding sequences (CDS), tRNAs, and rRNAs. CDS counts are extracted and normalized (x1000) for filtered genomes.

3. **BGC detection — antiSMASH**
   Identifies biosynthetic gene clusters (BGCs) and calculates the genomic percentage occupied by secondary metabolite clusters (`bgc_extractor.py`).

4. **Ortholog and domain search — eggNOG + HMMER**
   Identifies orthologous gene groups and protein domain patterns associated with cell division (`ftsQ`, `ftsA`) and elongation (`mreD`, `mrdA`, `rodZ`) genes (`patterns.py`).

5. **Visualization — grafica.R**
   Integrates all filtered outputs (completeness, assembly size, CDS count, BGC percentage, gene patterns) into a single comparative figure.

> Detailed step-by-step commands for each stage are documented in **aun debo agregar el link | falta renderizar y corregir la guía ⚠️** .


## Repository structure 🌳

```
├── README.md
├── index.html                      # GitHub Pages site
├── bgc_extraccion/
│   └── bgc_extractor.py            # Extracts BGC genomic extent from antiSMASH output
├── Genome_filtering_guide/
│   └── Guia_grafica.ipynb          # Step-by-step genome filtering pipeline (CheckM, Prokka, antiSMASH, eggNOG)
├── grafica/
│   ├── grafica.R                   # Main script to reproduce comparative figures
│   └── grafica_con_comentarios.R   # Same script, annotated for reference
├── patterns_extraccion/
│   └── patterns.py                 # Extracts cell division/elongation gene patterns
```

## Cluster repository structure 💻




## Requirements / Dependencies 📦 

### Software

| Tool | Version | Use |
|:-------------|:---------:|:-----|
| R | 2026.01.1+403 | Visualization
| CheckM  | v1.2.2 | Genome quality assessment (completeness, contamination, assembly size)
| Prokka  | v1.14.6 | Genome annotation, CDS count
| antiSMASH  | v7.1.0 | Biosynthetic gene cluster (BGC) detection 

## Cluster dependencies 



## References 📚

- Volland, J. M., Gonzalez-Rizzo, S., Gros, O., Tyml, T., Ivanova, N., Schulz, F., ... & Girguis, P. R. (2022). A centimeter-long bacterium with DNA contained in metabolically active, membrane-bound organelles. *Science*, 376(6600), 1453-1458. https://doi.org/10.1126/science.abb3634




> 🚧 **Work in progress** — this README is still under construction.

<p align="center">
  <img src="https://github.com/ecoevolab/thiotrichales_genomics/blob/main/cartoon.gif" width="300"/>
</p>
