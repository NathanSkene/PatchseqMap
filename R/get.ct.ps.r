#' Get celltype p-values
#'
#' Evaluates for cell type 'i' which of the reference cell types it is most likely to be. The probability values are returned as a vector.
#'
#' @param cellTypesToConsider List of strings for cell type names from the reference dataset
#' @param specificity A specificity matrix, generally for level2 celltypes, obtained using EWCE (only the specificity matrix (not the ctd) should be passed)
#' @param hitGenes Genes which fulfull the following criteria: (1) Expressed in the reduced patchseq dataset, i.e. they have standard deviation over a given 
#' threshold; (2) They have at least one read in patchseq cell 'i'; (3) They are found in both patchseq and specificity datasets
#' @param patchCellI Index of which patchseq cell is being considered
#' @param patchSeqExp Gene by Cell matrix of the patchseq data
#' @param bg Vector of gene symbols. All of these must be in the specificity matrix.
#'
#' @return A vector of p-values corresponding to cellTypesToConsider, where lowest value indicates that celltype is the most likely type for cell i
#'
#' @examples
#' plot_crayons()
#'
#' @export
get.ct.ps <- function(cellTypesToConsider,specificity,hitGenes,patchCellI,patchSeqExp,bg){
    # Check that at least three cell types are being considered
    if(length(cellTypesToConsider)<2){
        stop("ERROR: at least three cell types must be considered")
    }
    
    # Check that all cellTypesToConsider are in the specificity data
    if(sum(!cellTypesToConsider %in% colnames(specificity))>0){
        stop("ERROR: all cellTypesToConsider must be a column name of specificity")
    }
    
    # Check that all bg and hitGenes are in the specificity matrix
    if(sum(!bg %in% rownames(specificity))>0){stop("ERROR: all background genes must be in the specificity matrix")}
    if(sum(!hitGenes %in% rownames(specificity))>0){stop("ERROR: all hit genes must be in the specificity matrix")}
    
    whichCT = rep(1,length(cellTypesToConsider))
    names(whichCT) = cellTypesToConsider
    
    # Check there are at least 20 hit genes
    if(length(hitGenes)<2){return(whichCT)}
    
    for(ct in cellTypesToConsider){
        print(ct)
        # Multiply the Patchseq reads by Celltype specificity
        actual = specificity[hitGenes,ct] * patchSeqExp[hitGenes,patchCellI]
        # Randomise the specificity values, multiply those with Patchseq data
        bootstrap = replicate(reps,sample(specificity[bg,ct],length(hitGenes),replace = FALSE)) * patchSeqExp[hitGenes,patchCellI]
        # What was the probability of getting a score this high (relative to the randomised values)?
        sumAc = sum(actual)
        sumBt = apply(bootstrap,2,sum)
        p = 1-sum(sumAc>sumBt)/length(sumBt)
        whichCT[ct] = p
        #print(whichCT)
    }    
    names(whichCT) = cellTypesToConsider
    return(whichCT)
}