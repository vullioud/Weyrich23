#' load and shape the annotation table
#'
#' @param path_to_annotation_table path to a .gtf table containing the annotation.
#' @param keep the typoe of annotation to take into account.
#' @param output The format of the output: either "data.frame" or "GRanges"
#' @return A GRange of a data.frame  with the genomic coordinates, the geneID and the type.
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#'x <- load_and_shape_annotation_table(path_to_annotation_table, keep = "all", output = "data.frame")
#' }
#'
#' @export
load_and_shape_annotation_table <- function(path_to_annotation_table, keep = "all", output = "data.frame"){

  annotation_table <- rtracklayer::readGFF(path_to_annotation_table) %>%
    dplyr::select(.data$seqid, .data$start, .data$end, .data$type, .data$gene_id) %>%
    dplyr::rename(chr = .data$seqid,
                  stop = .data$end) %>%
    dplyr::distinct()

  if("all" %in% keep){
    out <- annotation_table
  } else {
    out <-  annotation_table %>%
      dplyr::filter(.data$type %in% keep)
  }

  if(output == "GRanges") {
    out <- GenomicRanges::makeGRangesFromDataFrame(out, keep.extra.columns = T)
  }
  out
}



#' Annotate tables
#'
#'Add a column to the table with the gene name.
#' @param path_to_annotation_table  path to a .gtf table containing the annotation.
#' @param table a table to annotate. Must contain the col. chr/start /stop
#' @param keep annotation to take into consideration.
#' @return the table to annotad with an extra list col: annotation
#'
#' @export
add_annotation <- function(path_to_annotation_table, table, keep){


  annotation_GR <- load_and_shape_annotation_table(path_to_annotation_table, keep = keep, output = "GRanges")

  table_GR <-  GenomicRanges::makeGRangesFromDataFrame(table, keep.extra.columns = T)

  annotated_table <- data.frame(annotatr::annotate_regions(table_GR, annotation = annotation_GR))

  out <-  as.data.frame(annotated_table) %>%
    dplyr::select(-.data$annot.start, -.data$annot.seqnames, -.data$annot.end,
                  -.data$width, -.data$strand, -.data$annot.width, -.data$annot.strand) %>%
    tidyr::nest(annotation = c(.data$annot.type, .data$annot.gene_id)) %>%
    dplyr::rename(chr = .data$seqnames,
                  stop = .data$end) %>%
    dplyr::mutate(chr = as.character(.data$chr),
                  start = as.character(.data$start),
                  stop = as.character(.data$stop))

  table %>% dplyr::mutate(chr = as.character(.data$chr),
                          start = as.character(.data$start),
                          stop = as.character(.data$stop)) %>%
    dplyr::left_join(out)
}

#' find annotation type
#' @param annotation_table A annotation data frame with "type" columns.
#' @param path_to_annotation_table the annotation table
#' @return the vector with the diffent unique annotation types
#'
#' @export
find_annotation_type <- function(annotation_table = NULL, path_to_annotation_table = NULL) {

  if(!is.null(annotation_table)){

    levels(as.factor(annotation_table$type))
  } else {
    annotation_table <- load_and_shape_annotation_table(path_to_annotation_table)
    levels(as.factor(annotation_table$type))
  }
}


#' Add annotation to loc genes
#' @param annotation_table A annotation data frame with "type" columns.
#' @param loc_table A table with the new gene names
#' @return A annotation tables with changed names for loc genes with new annotation
#'
.rename_loc <- function(annotation_table, loc_table) {

  if (is.null(annotation_table) | !any(annotation_table$annot.gene_id %in% loc_table$loc)) {

    out <- annotation_table

  } else  {
    out <- annotation_table %>%
      dplyr::mutate(annot.gene_id = ifelse(!(.data$annot.gene_id %in% loc_table$loc), .data$annot.gene_id, loc_table$annot.gene_id[match(.data$annot.gene_id, loc_table$loc)]))
  }
  out
}


