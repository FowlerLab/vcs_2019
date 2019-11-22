## load library
library(Seurat)

## Seurat cell cycle analysis
### Arguments: matrix file location; feature data; phenotypic data
### Output: phenotypic dataframe + cell cycle scores, cell cycle designations
### Based on https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html
seurat_cellcycle <- function(mat_file, pd, fd){

  # Read in the expression matrix The first row is a header row, the first column is rownames
  exp.mat <- readMM(file = mat_file)
  row.names(exp.mat) = fd$gene_short_name
  colnames(exp.mat) = pd$cell

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  # Create our Seurat object and complete the initalization steps
  marrow <- CreateSeuratObject(counts = exp.mat)
  marrow <- NormalizeData(marrow)
  marrow <- FindVariableFeatures(marrow, selection.method = "vst")
  marrow <- ScaleData(marrow, features = rownames(marrow))
  marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
  marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  # merg seurat cell cycle object with the phenodata frame
  orig_order    <- pd$cell # save the original order of cells
  pd            <- merge(pd, marrow[[]][,c('S.Score', 'G2M.Score', 'Phase')], by.x = 'cell', by.y = 'row.names') # add cell cycle info
  pd$sample     <- sapply(pd$cell, function(x) unlist(strsplit(as.character(x), '-'))[[2]])
  row.names(pd) <- pd$cell
  pd            <- pd[match(orig_order, row.names(pd)),]

  return(pd)

  }

## preprocess cds
### Inputs: matrix file location; feature data; phenotypic data
### Output: cds with seurat cell cycle annotations, processed with reduced_dimensions to plot UMAP
preprocess_seurat <- function(mat_file, pd, fd, annotate_seurat = TRUE){

  # add_seurat
  if(annotate_seurat) {
    print('Annotating with Seurat cell cycle analysis..')
    pd <- seurat_cellcycle(mat_file, pd, fd)}

  # load rep1 data
  print('Creating cds..')
  cds <- new_cell_data_set(readMM(mat_file),
                                cell_metadata = pd,
                                gene_metadata = fd)

  # preprocess
  print('Preprocessing cds..')
  cds  <- preprocess_cds(cds, num_dim = 100)

  # remove outlier cells
  print('Removing outlier cells..')
  low_size_factor_val <- mean(colData(cds)$Size_Factor) - (2*sd(colData(cds)$Size_Factor))
  hi_size_factor_val <- mean(colData(cds)$Size_Factor) + (2*sd(colData(cds)$Size_Factor))
  cds <- cds[,which(cds$Size_Factor > low_size_factor_val & cds$Size_Factor < hi_size_factor_val )]

  # remove genes detected < 10 times
  gene_total_umis <- rowSums(exprs(cds))
  to_keep         <- names(gene_total_umis)[which(gene_total_umis > 9)]
  cds             <- cds[to_keep,]

  # plot
  print('Preparing UMAP..')
  cds = reduce_dimension(cds, reduction_method = 'UMAP', umap.n_neighbors = 30)

  return(cds)}
