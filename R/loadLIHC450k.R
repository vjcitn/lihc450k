#' use BiocFileCache discipline to acquire lihc450K SummarizedExperiment
#' @importFrom rhdf5client H5S_source
#' @importFrom rhdf5client setPath
#' @importFrom rhdf5client H5S_dataset2
#' @importFrom HDF5Array HDF5ArraySeed
#' @importFrom DelayedArray DelayedArray
#' @import restfulSE
#' @param remoteSEPath character(1) identifying remote RDS with assay-free SummarizedExperiment instance
#' @param remoteHSDSPath character(1) identifying HSDS server/port
#' @param HDF5_in_S3_Path character(1) identifying file in S3 bucket with 1.2GB of HDF5 to be downloaded and cached for local use
#' @param hostPath character(1) the data 'host' -- looks strange with '/home' default but it is OK
#' @param cache defaults to BiocFileCache() value
#' @param useHSDS logical indicating, if TRUE, that 450k numerical data to be retrieved from HSDS; if FALSE, 1.2 GB of HDF5 will be retrieved from S3 and locally cached for future use
#' @note The RDS for the SummarizedExperiment is in an AWS S3 bucket.
#' This function will check local cache for the data and will download
#' to cache (this is a default behavior of bfcrpath())
#' if not found.  That download is a one-time operation for
#' any given value of \code{cache}.  Note further that the `remoteHSDSPath`
#' is provided courtesy of the HDF Group in early 2018 and may not be persistent.
#' @examples
#' lim = loadLIHC450k()
#' lim
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
  HDF5_in_S3_Path =
        "https://s3.amazonaws.com/biocfound-tcga/lihc450k_imp.h5",
  hostPath = "/home/stvjc/lihc450k.h5", 
  cache=BiocFileCache(),
  useHSDS = TRUE) {
    if (!checkCache_lihcSEshell(cache, remoteSEPath)) {
         message("adding RDS for SE 'shell' to local cache, future invocations will use local image")
         }
    else { message("found SE 'shell' in cache, will retrieve") }
    path = bfcrpath(cache, remoteSEPath) # downloads!
    shell = readRDS(path)
  if (useHSDS) {
    con = H5S_source(remoteHSDSPath)
    lic = setPath(con, hostPath)
    ds2 = H5S_dataset2(lic)
    bet = DelayedArray(new("H5S_ArraySeed", filepath="", domain="", host="",
       H5S_dataset=ds2))
    }
   else {
    if (!checkCache_lihc450kHDF5(cache=cache, remoteS3URL=HDF5_in_S3_Path))
        message("local HDF5 for LIHC450k not found in cache, one-time-download of 1.2GB into your cache for all future use occurring now")
    hdf5datpath = bfcrpath(cache, HDF5_in_S3_Path)
    bet = DelayedArray(HDF5ArraySeed(filepath=hdf5datpath, "assay001"))
    }
  assays(shell) = S4Vectors::SimpleList(betas=bet)
  shell
}

checkCache_lihcSEshell = function(cache=BiocFileCache(), remoteSEURL) {
 allr = bfcinfo(cache)$rname
 remoteSEURL %in% allr
}

checkCache_lihc450kHDF5 = function(cache=BiocFileCache(), remoteS3URL) {
 allr = bfcinfo(cache)$rname
 remoteS3URL %in% allr
}
