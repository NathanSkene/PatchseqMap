
<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml">

<body>


<div class="container-fluid main-container">

<div class="fluid-row" id="header">



<h1 class="title toc-ignore">Mapping Patchseq Data onto Celltypes using bootstrapping</h1>

</div>


<div id="using-ewce-for-mapping-patchseq-data-onto-known-cell-types" class="section level1">
<h1>Using EWCE for mapping PatchSeq data onto known cell types</h1>
</div>
<div id="first-install-ewce-from-github-to-get-access-to-the-cortex_mrna-dataset" class="section level1">
<h1>First, install EWCE from github to get access to the cortex_mrna dataset</h1>
<pre class="r"><code># install.packages(&quot;devtools&quot;)
library(devtools)
install_github(&quot;nathanskene/ewce&quot;)</code></pre>
<pre><code>## Skipping install of 'EWCE' from a github remote, the SHA1 (d6666356) has not changed since last install.
##   Use `force = TRUE` to force installation</code></pre>
<pre class="r"><code>library(EWCE)
library(PatchseqMap)</code></pre>
<div id="simulating-patchseq-data" class="section level2">
<h2>Simulating PatchSeq data</h2>
<p>Patch-seq is a useful technique for obtaining both electrophysiology and mRNA from single cells. One of the primary purposes for collecting the mRNA is to be able to infer the underlying celltype. Generally fairly low amounts of mRNA are obtained. For demonstration purposes, let’s generate some simulated patchseq data based on the cortex/hippocampus single cell transcriptome data. We assume that there are 3x fewer reads through patchseq, sampled from the scRNA-Seq data with a poisson distribution.</p>
<pre class="r"><code># First reduce the cortex_mrna dataset to just interneurons
data(cortex_mrna)
int_annot = cortex_mrna$annot[cortex_mrna$annot$level1class==&quot;interneurons&quot;,]
int_exp = cortex_mrna$exp[,int_annot$cell_id]

# Use the cortex_mrna dataset to simulate patchseq data
downsample_expression &lt;- function(x){
    y = round(rpois(length(x),x/4))
    return(y)
}
simulated_patchdata = apply(int_exp,2,downsample_expression)
rownames(simulated_patchdata) = rownames(int_exp)
simulated_patchdata = simulated_patchdata[apply(simulated_patchdata,1,sum)&gt;0,]
MedianReadsPerCell = median(apply(simulated_patchdata,2,sum))
print(sprintf(&quot;Median reads per simulated cell: %s&quot;,MedianReadsPerCell))</code></pre>
<pre><code>## [1] &quot;Median reads per simulated cell: 4306.5&quot;</code></pre>
</div>
<div id="prepare-cell-type-specificity-data" class="section level2">
<h2>Prepare cell type specificity data</h2>
<pre class="r"><code># Generate celltype data for just the cortex/hippocampus data
exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=cortex_mrna$exp,level2annot = cortex_mrna$annot$level2class)
fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,level1class=cortex_mrna$annot$level1class,level2class=cortex_mrna$annot$level2class,groupName=&quot;kiCortexOnly&quot;,thresh=0,trim=0)
print(fNames_CortexOnly)</code></pre>
<pre><code>## [1] &quot;CellTypeData_kiCortexOnly.rda&quot;</code></pre>
<pre class="r"><code>fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
print(fNames_CortexOnly)</code></pre>
<pre><code>## [1] &quot;CellTypeData_kiCortexOnly.rda&quot;         
## [2] &quot;CellTypeData_kiCortexOnly_1to1only.rda&quot;</code></pre>
<pre class="r"><code>load(fNames_CortexOnly[1])

# Because the simulated patchseq data is just for interneurons, restrict the celltype specificity data to interneurons
ctd[[2]]$specificity = ctd[[2]]$specificity[,grep(&quot;^Int&quot;,colnames(ctd[[2]]$specificity))]</code></pre>
</div>
<div id="mapping-simulated-patchseq-data" class="section level2">
<h2>Mapping simulated PatchSeq data</h2>
<p>As ever, you should start by correcting the MGI symbols in your patchseq data:</p>
<pre class="r"><code>simulated_patchdata_corrected = fix.bad.mgi.symbols(simulated_patchdata)</code></pre>
<pre><code>## [1] &quot;845 rows do not have proper MGI symbols&quot;
##  [1] &quot;2310042E22Rik&quot; &quot;BC005764&quot;      &quot;C130030K03Rik&quot; &quot;Stmn1-rs1&quot;    
##  [5] &quot;Gm9846&quot;        &quot;E130309F12Rik&quot; &quot;Fam211b&quot;       &quot;AI848285&quot;     
##  [9] &quot;Acpl2&quot;         &quot;9630033F20Rik&quot; &quot;Adrbk2&quot;        &quot;Syne1_loc2&quot;   
## [13] &quot;Adc&quot;           &quot;Dlx1os&quot;        &quot;LOC106740&quot;     &quot;Pdxp&quot;         
## [17] &quot;Atp6v0c-ps2&quot;   &quot;2900056M20Rik&quot; &quot;Epb4.1l1&quot;      &quot;A330050F15Rik&quot;</code></pre>
<pre><code>## Warning in fix.bad.mgi.symbols(simulated_patchdata): Possible presence of
## excel corrupted date-like genes: Sepw1, Sepp1, Sept15, Sepn1</code></pre>
<pre><code>## [1] &quot;6 poorly annotated genes are replicates of existing genes. These are: &quot;
## [1] &quot;B3gnt2&quot;  &quot;Bhlhe40&quot; &quot;Hjurp&quot;   &quot;Naa38&quot;   &quot;Bcor&quot;    &quot;Ubl4a&quot;  
## [1] &quot;480 rows should have been corrected by checking synonms&quot;
## [1] &quot;373 rows STILL do not have proper MGI symbols&quot;</code></pre>
<pre class="r"><code>reps=100

allRes = map.patchseq(ctd,simulated_patchdata_corrected[,1:20],sdThresh=1,annotLevel=2,useSDthresh=TRUE,reps=reps)</code></pre>
<pre><code>## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int6&quot;
## [1] &quot;loop&quot;
## Int12 Int10  Int1 Int16  Int8 Int15 Int13  Int5  Int9 Int11  Int3  Int6 
##  0.11  0.13  0.37  0.38  0.50  0.77  0.82  0.86  0.89  0.91  0.93  0.95 
## Int14  Int2  Int4  Int7 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int12&quot;
## [1] 1
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## Int10 Int12 Int16  Int6  Int3  Int8  Int9 Int15 Int13  Int5 Int14  Int1 
##  0.02  0.14  0.15  0.16  0.18  0.23  0.29  0.44  0.66  0.82  0.94  0.95 
##  Int2 Int11  Int4  Int7 
##  0.97  0.98  0.99  0.99 
## [1] &quot;Int10&quot;
## [1] 2
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;loop&quot;
##  Int6 Int10 Int11 Int16 Int15 Int14  Int8 Int13  Int9  Int1  Int7  Int3 
##  0.09  0.15  0.22  0.43  0.62  0.69  0.71  0.77  0.84  0.91  0.93  0.94 
##  Int5 Int12  Int2  Int4 
##  0.96  0.99  1.00  1.00 
## [1] &quot;Int6&quot;
## [1] 3
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int6&quot;
## [1] &quot;loop&quot;
## Int10  Int5  Int8  Int3 Int15  Int9 Int13 Int14 Int11  Int7  Int1  Int2 
##  0.02  0.15  0.20  0.25  0.45  0.50  0.80  0.80  0.82  0.89  0.94  0.98 
##  Int6 Int12 Int16  Int4 
##  0.99  1.00  1.00  1.00 
## [1] &quot;Int10&quot;
## [1] 4
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
##  Int9 Int13  Int8 Int14 Int10  Int4 Int15  Int3 Int12 Int11  Int1  Int2 
##  0.00  0.49  0.54  0.55  0.62  0.70  0.80  0.80  0.87  0.92  0.98  0.99 
## Int16  Int5  Int6  Int7 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int9&quot;
## [1] 5
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
##  Int9  Int3  Int8 Int14 Int13 Int15  Int4  Int6  Int7 Int11  Int5 Int16 
##  0.03  0.16  0.43  0.44  0.45  0.48  0.78  0.78  0.78  0.92  0.93  0.96 
## Int10 Int12  Int2  Int1 
##  0.97  0.97  0.98  1.00 
## [1] &quot;Int9&quot;
## [1] 6
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## Int10 Int11  Int5 Int12  Int6 Int15 Int13  Int9  Int8  Int3  Int4 Int16 
##  0.00  0.11  0.14  0.15  0.16  0.40  0.42  0.49  0.59  0.79  0.80  0.94 
## Int14  Int7  Int2  Int1 
##  0.96  0.96  0.98  1.00 
## [1] &quot;Int10&quot;
## [1] 7
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
##  Int9  Int4 Int13 Int14  Int3  Int2 Int15  Int7 Int12  Int8  Int1  Int5 
##  0.06  0.11  0.14  0.17  0.24  0.37  0.44  0.63  0.64  0.82  0.87  0.88 
##  Int6 Int10 Int11 Int16 
##  0.89  0.92  0.94  0.96 
## [1] &quot;Int9&quot;
## [1] 8
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## Int10 Int11  Int8 Int13 Int16 Int12 Int14 Int15  Int9  Int3  Int6  Int4 
##  0.00  0.11  0.13  0.27  0.31  0.47  0.51  0.56  0.63  0.69  0.72  0.78 
##  Int1  Int7  Int2  Int5 
##  0.93  0.97  0.99  0.99 
## [1] &quot;Int10&quot;
## [1] 9
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;loop&quot;
## Int10  Int6 Int16  Int4  Int8 Int13  Int1 Int14  Int7 Int15  Int3  Int5 
##  0.01  0.19  0.21  0.27  0.27  0.34  0.35  0.40  0.40  0.46  0.47  0.59 
##  Int2 Int11  Int9 Int12 
##  0.82  0.88  0.97  1.00 
## [1] &quot;Int10&quot;
## [1] 10
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;loop&quot;
## Int10  Int9  Int6 Int13  Int8 Int15  Int5 Int11  Int3  Int4  Int1 Int16 
##  0.16  0.16  0.37  0.51  0.60  0.64  0.67  0.72  0.77  0.80  0.81  0.81 
## Int12  Int7 Int14  Int2 
##  0.85  0.89  0.93  0.94 
## [1] &quot;Int10&quot;
## [1] 11
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
##  Int1  Int2  Int3  Int8 Int13  Int5  Int7 Int15 Int14 Int16 Int11 Int12 
##  0.00  0.24  0.31  0.32  0.61  0.62  0.76  0.87  0.97  0.97  0.98  0.99 
##  Int6 Int10  Int4  Int9 
##  0.99  1.00  1.00  1.00 
## [1] &quot;Int1&quot;
## [1] 12
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
##  Int1  Int4 Int11 Int14  Int7 Int15 Int16  Int3  Int5 Int13  Int2  Int8 
##  0.21  0.26  0.76  0.76  0.81  0.85  0.89  0.89  0.90  0.92  0.93  0.93 
##  Int9 Int12 Int10  Int6 
##  0.98  0.99  1.00  1.00 
## [1] &quot;Int1&quot;
## [1] 13
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
##  Int1  Int2  Int3  Int7 Int13  Int8  Int5 Int15 Int11 Int16 Int14 Int10 
##  0.01  0.27  0.68  0.80  0.84  0.84  0.87  0.96  0.98  0.98  0.99  1.00 
## Int12  Int4  Int6  Int9 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int1&quot;
## [1] 14
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
##  Int1 Int14  Int7  Int3 Int16  Int5 Int11 Int15 Int12 Int13  Int2  Int8 
##  0.00  0.30  0.86  0.89  0.90  0.94  0.95  0.95  0.97  0.97  0.97  0.97 
## Int10  Int4  Int6  Int9 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int1&quot;
## [1] 15
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
##  Int1  Int7  Int3  Int2 Int13  Int8 Int11 Int14  Int5 Int10 Int15 Int16 
##  0.00  0.74  0.79  0.84  0.87  0.89  0.91  0.92  0.96  0.98  0.98  0.99 
##  Int6 Int12  Int4  Int9 
##  0.99  1.00  1.00  1.00 
## [1] &quot;Int1&quot;
## [1] 16
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;loop&quot;
##  Int2  Int5 Int13  Int7  Int8 Int14 Int15 Int16  Int6  Int1 Int11 Int10 
##  0.13  0.15  0.16  0.28  0.31  0.66  0.74  0.79  0.80  0.84  0.86  1.00 
## Int12  Int3  Int4  Int9 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int2&quot;
## [1] 17
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int1&quot;
## [1] &quot;Int2&quot;
## [1] &quot;loop&quot;
##  Int2 Int14 Int11 Int15 Int13 Int12  Int8  Int5  Int6  Int7  Int1 Int10 
##  0.04  0.11  0.41  0.59  0.77  0.78  0.91  0.96  0.96  0.97  0.99  1.00 
## Int16  Int3  Int4  Int9 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int2&quot;
## [1] 18
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
##  Int3 Int13  Int1  Int8  Int2  Int5 Int15 Int14 Int10 Int11 Int12 Int16 
##  0.05  0.36  0.38  0.45  0.54  0.77  0.91  0.98  1.00  1.00  1.00  1.00 
##  Int4  Int6  Int7  Int9 
##  1.00  1.00  1.00  1.00 
## [1] &quot;Int3&quot;
## [1] 19
## [1] &quot;Int1&quot;
## [1] &quot;Int10&quot;
## [1] &quot;Int11&quot;
## [1] &quot;Int12&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int2&quot;
## [1] &quot;Int3&quot;
## [1] &quot;Int4&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int6&quot;
## [1] &quot;Int7&quot;
## [1] &quot;Int8&quot;
## [1] &quot;Int9&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int14&quot;
## [1] &quot;Int15&quot;
## [1] &quot;Int16&quot;
## [1] &quot;Int5&quot;
## [1] &quot;Int8&quot;
## [1] &quot;loop&quot;
## [1] &quot;Int13&quot;
## [1] &quot;Int8&quot;
## [1] &quot;loop&quot;
## Int13  Int3  Int1  Int2 Int15  Int7  Int5  Int4  Int8 Int11 Int12 Int14 
##  0.14  0.30  0.37  0.42  0.52  0.56  0.73  0.84  0.87  0.88  0.98  0.98 
##  Int6 Int10 Int16  Int9 
##  0.98  1.00  1.00  1.00 
## [1] &quot;Int13&quot;
## [1] 20</code></pre>
<pre class="r"><code>print(data.frame(Actual=int_annot$level2class[1:20],Assigned=allRes$assigned))</code></pre>
<pre><code>##    Actual Assigned
## 1   Int10    Int12
## 2   Int10    Int10
## 3    Int6     Int6
## 4   Int10    Int10
## 5    Int9     Int9
## 6    Int9     Int9
## 7   Int10    Int10
## 8    Int9     Int9
## 9   Int10    Int10
## 10  Int10    Int10
## 11  Int10    Int10
## 12   Int2     Int1
## 13   Int4     Int1
## 14   Int2     Int1
## 15   Int1     Int1
## 16   Int1     Int1
## 17   Int2     Int2
## 18   Int2     Int2
## 19   Int3     Int3
## 20  Int13    Int13</code></pre>
</div>
</div>




</div>

<!-- dynamically load mathjax for compatibility with self-contained -->


</body>
</html>
