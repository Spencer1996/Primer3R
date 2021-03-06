#' Primer3_Hotset
#'
#'Primer3_Hotset is used to form the file, which is required by Primer3. Here tarquire sequence
#'do not need tp be sheared
#'
#' @param x is the is dataframe after extension，here, each taget sequence requires the names (NO.)
#' @param generic is the input of The PRIMER_TASK that tells primer3 which type of primers to pick.
#' You can select typical primers for PCR detection, primers for cloning or for sequencing. Here,
#' PRIMER_TASK alsoccalled 'pick_detection_primers' is default setting  to 'generic' while
#' retaining 'pick_detection_primers' as an alias for backward compatibility.
#' @param oligo is the input of PRIMER_PICK_INTERNAL_OLIGO, which default is 0 .If the associated
#'  value = 1 (non-0), then primer3 will attempt to pick an internal oligo (hybridization probe to
#'  detect the PCR product).
#' @param rangel is the left range  of PRIMER_PRODUCT_SIZE_RANGE input. The associated values specify the lengths
#' of the product that the user wants the primers to create, and is a space separated list of elements
#' of the form <x>-<y>,where an <x>-<y> pair is a legal range of lengths for the product. For example,
#' if one wants PCR products to be between 100 to 150 bases (inclusive) then one would set this parameter
#'  to 100-150. Here, default is 100-150
#' @param ranger is the lright  range  of PRIMER_PRODUCT_SIZE_RANGE input.
#' @param opt_size is the input of PRIMER_OPT_SIZE, which is the Optimum length (in bases) of a primer.
#'  Primer3 will attempt to pick primers close to this length. Here, default is 20.
#' @param min_size is the input of PRIMER_MIN_SIZE. The minimum acceptable length of a primer.
#'  Must be greater than 0 and less than or equal to PRIMER_MAX_SIZE. Here,default is 16
#' @param max_size is the input of PRIMER_MAX_SIZE, The maximum acceptable length of a primer.
#'  Must be greater than 0 and less than or equal to PRIMER_MIN_SIZE. Here,default is 22
#' @param explain_fiag is the input of If this flag is 1 (non-0), produce PRIMER_LEFT_EXPLAIN,
#' PRIMER_RIGHT_EXPLAIN, PRIMER_INTERNAL_EXPLAIN and/or PRIMER_PAIR_EXPLAIN output tags
#' as appropriate. These output tags are intended to provide information on the number
#'  of oligos and primer pairs that primer3 examined and counts of the number discarded
#'  for various reasons. If -format_output is set similar information is produced in the
#'   user-oriented output. Here,default is 1
#' @param return is the input of PRIMER_NUM_RETURN, which is the maximum number of primer (pairs) to return. Primer pairs
#' returned are sorted by their "quality", in other words by the value of the objective function
#' (where a lower number indicates a better primer pair). Caution: setting this parameter to a
#' large value will increase running time. Here,default is 5
#' @param min_tm Minimum acceptable melting temperature (Celsius) for a primer oligo. Here,default is 45
#' @param opt_tm Optimum melting temperature (Celsius) for a primer. Primer3 will try
#'to pick primers with melting temperatures are close to this temperature. The oligo
#'melting temperature formula used can be specified by user. Please see PRIMER_TM_FORMULA
#'for more information. Here,default is 60
#' @param max_tm Maximum acceptable melting temperature (Celsius) for a primer oligo.Here,default is 60
#'
#' @return dataframe
#' @export
#'
#' @examples
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' Hot_extension(hot,185,chromosome)
#' p.hot<-Primer3_Hotset(entend.hot,return=2)
#' View(p.hot)
#' write.table(x,"p.hot.txt",quote = FALSE,row.names = FALSE, col.names = FALSE)
Primer3_Hotset<-function(x,generic="generic", #x<-hot.ex
                               oligo=0,rangel=100,ranger=150,
                               opt_size=20,min_size=16,
                               max_size=22,explain_fiag=1,return=50,
                               min_tm=45,opt_tm=60,max_tm=60){
  wid_range<-function(x){nchar(as.character(x))}
  if(unique(x$strand)=="+"){right=1;left=0}else{right=0;left=1}
  NO<-x$NO.
  seq<-x$Extend.seq
  tseq<-x$Target.Seq
  IO <- character()
  for (i in 1:(length(seq))){
    IO[18*i-17] <- paste0("SEQUENCE_ID=",NO[i])
    IO[18*i-16] <- paste0("SEQUENCE_TEMPLATE=",seq[[i]])
    IO[18*i-15] <- paste0("PRIMER_TASK=",generic)
    IO[18*i-14] <- paste0("SEQUENCE_TARGET=",paste0(as.character(c(ifelse(x$strand[i]=="+",1,(wid_range(tseq[i])+1)),wid_range(tseq[i]))),collapse=","))
    IO[18*i-13] <- paste0("SEQUENCE_EXCLUDED_REGION=",paste0(as.character(c(ifelse(x$strand[i]=="+",1,(wid_range(tseq[i])+1)),wid_range(tseq[i]))),collapse=","))
    IO[18*i-12] <- paste0("PRIMER_PICK_RIGHT_PRIMER=",right)
    IO[18*i-11] <- paste0("PRIMER_PICK_LEFT_PRIMER=",left)
    IO[18*i-10] <- paste0("PRIMER_PICK_INTERNAL_OLIGO=",oligo)
    IO[18*i-9] <- paste0("PRIMER_PRODUCT_SIZE_RANGE=",paste(as.character(c(rangel,ranger)),collapse="-"))
    IO[18*i-8] <- paste0("PRIMER_OPT_SIZE=",opt_size)
    IO[18*i-7] <- paste0("PRIMER_MIN_SIZE=",min_size)
    IO[18*i-6] <- paste0("PRIMER_MAX_SIZE=",max_size)
    IO[18*i-5] <- paste0("PRIMER_EXPLAIN_FLAG=",explain_fiag)
    IO[18*i-4] <-paste0("PRIMER_NUM_RETURN=",return)
    IO[18*i-3] <-paste0("PRIMER_MIN_TM=",min_tm)
    IO[18*i-2] <-paste0("PRIMER_OPT_TM=",opt_tm)
    IO[18*i-1] <-paste0("PRIMER_MAX_TM=",max_tm)
    IO[18*i] <- c("=")
  }
  return(as.data.frame(IO))}
