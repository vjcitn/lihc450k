#' use BiocFileCache discipline to acquire lihc450K SummarizedExperiment
#' @importFrom rhdf5client H5S_source
#' @importFrom rhdf5client setPath
#' @importFrom rhdf5client H5S_dataset2
#' @import restfulSE
#' @param remoteSEPath character(1) identifying remote RDS with assay-free SummarizedExperiment instance
#' @param remoteHSDSPath character(1) identifying HSDS server/port
#' @param hostPath character(1) the data 'host' -- looks strange with '/home' default but it is OK
#' @param cache defaults to BiocFileCache() value
#' @note The RDS for the SummarizedExperiment is in an AWS S3 bucket.
#' This function will check local cache for the data and will download
#' to cache if not found.  That download is a one-time operation for
#' any given value of \code{cache}.  Note further that the `remoteHSDSPath`
#' is provided courtesy of the HDF Group in early 2018 and may not be persistent.
#' @examples
#' lim = loadLIHC450k()
#' assay(lim)
#' assay(lim)[1,429]
#' # 'confirm' with restfulSE interface to ISB-CGC
#' # cgcConnection() %>% tbl("DNA_Methylation_chr16") %>% 
#' #    filter( ParticipantBarcode=="TCGA-ZS-A9CG", Probe_Id == "cg00000029") %>% 
#' #    glimpse()
#' # Observations: ??
#' # Running job \:  1s:18.9 gigabytes processed
#' # Variables: 8                                                                  
#' # $ ParticipantBarcode   <chr> "TCGA-ZS-A9CG"
#' # $ SampleBarcode        <chr> "TCGA-ZS-A9CG-01A"
#' # $ SampleTypeLetterCode <chr> "TP"
#' # $ AliquotBarcode       <chr> "TCGA-ZS-A9CG-01A-11D-A36Y-05"
#' # $ Platform             <chr> "HumanMethylation450"
#' # $ Study                <chr> "LIHC"
#' # $ Probe_Id             <chr> "cg00000029"
#' # $ Beta_Value           <dbl> 0.65
#' #
#' @export
loadLIHC450k = function(
  remoteSEPath = 
        "https://s3.amazonaws.com/biocfound-tcga/lihc450k.rds", 
  remoteHSDSPath =
        "http://52.4.181.237:5101",
  hostPath = "/home/stvjc/lihc450k.h5", 
  cache=BiocFileCache()) {
    if (!checkCache_lihc(cache)) message("adding RDS to local cache, future invocations will use local image")
    path = bfcrpath(cache, remoteSEPath)
    shell = readRDS(path)
  con = H5S_source(remoteHSDSPath)
  lic = setPath(con, hostPath)
  ds2 = H5S_dataset2(lic)
  bet = DelayedArray(new("H5S_ArraySeed", filepath="", domain="", host="",
       H5S_dataset=ds2))
  assays(shell) = S4Vectors::SimpleList(betas=bet)
  shell
}
    
checkCache_lihc = function(cache=BiocFileCache()) {
 allr = bfcinfo(cache)$rname
 "https://s3.amazonaws.com/biocfound-tcga/lihc450k.rds" %in% allr
}
