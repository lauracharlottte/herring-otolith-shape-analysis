# ============================================================
# Outline extraction (Section 2.4 of the thesis)
# ------------------------------------------------------------
# Input: photographs of the right otoliths, already cropped, rotated
# and sorted into left/right by the Uploader application (MIR-PIB);
# these images are not part of this repository.
#
#   I.   binarise every photograph: greyscale, white square of
#        100 x 100 pixels over the centre, inversion, automatic
#        threshold (k-means on a pixel sample) x 1.4
#   II.  trace the outline in the binary image (Momocs::import_jpg);
#        extractions with fewer than 100 outline points count as failed
#   III. draw every extracted outline onto its photograph
#        (control images, Fig. 8d)
#
# Output: oto.RData (all outlines) and the folders Binary/ and
# Superimposed/ next to the photographs.
# ============================================================

## half width and half height (pixels) of the white centre square
xrange=50
yrange=50
adjust=0
threshold=0.5
cex=1

## ------------------------------------------------------------
## Which steps to run
## ------------------------------------------------------------
## All TRUE reproduces the original full pipeline: binarize every
## photograph, extract every contour, draw every control image,
## and overwrite oto.RData.
##
## The steps are separable because they are not always all
## needed. The common case after new photographs arrive is to
## draw the missing control images only, without touching the
## extraction the analysis already rests on:
##
##     RUN_BINARIZE <- FALSE
##     RUN_EXTRACT  <- FALSE
##     (SKIP_EXISTING stays TRUE)
##
## SKIP_EXISTING makes every step incremental -- it works on what
## is not already on disk. Set it to FALSE to redo everything.
RUN_BINARIZE    <- TRUE
RUN_EXTRACT     <- TRUE
RUN_SUPERIMPOSE <- TRUE
SKIP_EXISTING   <- TRUE

## Where step III takes its outlines from when step II is skipped.
## Results/oto_all.RData is the extraction the thesis
## results were computed on, so control images drawn from it show
## exactly the contours that were analysed. Falls back to F_OTO.
OUTLINE_SOURCE <- NULL   # NULL = decide below

library(Momocs)
library(tidyr)
library(foreach)
library(doParallel)
library(jpeg)
library(imager)

## ---- load project paths ------------------------------------
## Searches upwards from the working directory for config.R
## and sources it. Afterwards SHAPE_DIR,
## DIR_IMAGES, DIR_FOURIER, ensure_dir() etc. are available.
local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    hit <- file.path(d, "config.R")
    hit <- hit[file.exists(hit)]
    if (length(hit)) { source(hit[1]); break }
    if (identical(dirname(d), d)) {
      stop("config.R not found -- please setwd() into the repository folder.",
           call. = FALSE)
    }
    d <- dirname(d)
  }
})

## Image folder and output file are defined in config.R
fix.image.path <- DIR_IMAGES

if (is.null(OUTLINE_SOURCE)) {
  cand <- c(file.path(DIR_RESULTS, "oto_all.RData"), F_OTO)
  OUTLINE_SOURCE <- cand[file.exists(cand)][1]
}

cat("Steps: binarize", RUN_BINARIZE, "| extract", RUN_EXTRACT,
    "| superimpose", RUN_SUPERIMPOSE, "| skip existing", SKIP_EXISTING, "\n")

# Select pictures for shape extraction 
jpg.list.all<-list.files(path=fix.image.path, pattern = "Right.jpg")
jpg.list.all

oto<-list()

jpg.list<-jpg.list.all
### I. 1st step: Binarization 

cl <- makeCluster(detectCores())   #Creates parallel workers for faster computing. Useful if many images.
registerDoParallel(cl)
cl


### I. Binarisation with imager: each image is converted to greyscale, a white square
### is drawn over the centre (darker areas such as the nucleus would otherwise be taken
### for background), the grey values are inverted and an automatic threshold separates
### otolith and background. Binary images are saved in the sub-folder Binary/.

write.binary.path = paste(fix.image.path,"Binary/",sep="/")
if(!file.exists(write.binary.path)){
  suppressWarnings(dir.create(write.binary.path))
}

## SKIP_EXISTING leaves already-binarized photographs alone, so a
## new run after further images arrive only processes the new ones.
bin.todo <- jpg.list
if (SKIP_EXISTING) bin.todo <- setdiff(bin.todo, list.files(write.binary.path, pattern="\\.jpg$"))
if (!RUN_BINARIZE) bin.todo <- character(0)
cat("1. binarize:", length(bin.todo), "of", length(jpg.list), "photographs\n")

foreach(i=seq_along(bin.todo),.errorhandling="pass", .packages = c("doParallel",  "imager"))%dopar%{
  
  a <- paste(fix.image.path,"/", bin.todo[[i]], sep="") %>% load.image()
  im_width <- a %>% width()
  im_height <- a %>% height()
  
  a <- a%>%grayscale()%>%
    draw_rect( im_width/2-xrange, im_height/2-yrange, im_width/2+xrange, im_height/2+yrange, color = 1, opacity = 1, filled = TRUE) %>% 
    imchange(where = ~xs>0,fo = ~1-. )%>%#inverse color
    threshold("auto", approx = T, adjust = 1.4+adjust)    # automatic threshold x 1.4 (Section 2.4)
  a %>%as.cimg%>% imager::save.image( file = paste(write.binary.path,bin.todo[[i]], sep="")) 
}
message("1. Binary pictures are up to date!")


### II. 2nd step: Outline extraction 

if (RUN_EXTRACT) {

oto.new<-foreach(i=paste(fix.image.path,"/Binary/",jpg.list, sep=""),.errorhandling="pass", .packages = c("doParallel", "Momocs"))%dopar%{
  withCallingHandlers({ 
    setTimeLimit(10)
    import_jpg(i, threshold = threshold)  # binary image: 0.5 separates black and white
  }, 
  error = function(e) {
    print(i)###  error handling
  }
  )
}
message("2. All otolith outlines are extracted!")

### Verify that all outlines were properly extracted
oto.new <- do.call(c, oto.new)    


# Non-errors
noerror <- names(oto.new)[!names(oto.new)%in%c("message", "call")]
npoints <- unlist(lapply(oto.new, function(x) dim(x)[1]))
#npoints <- sapply(oto.new, function (x) dim(x%>%as.data.frame()))[1,]     
npoints <- names(npoints[npoints<100]) # shows number of points  along the otolith contour 
noerror <- noerror[!noerror%in%npoints]
noerror

# Errors:
errors <- jpg.list[!jpg.list%in%paste(noerror, ".jpg", sep="")]     # should show 0
errors

} else {
  ## Step II skipped: take the outlines from the stored extraction
  ## instead, so step III can still run on its own.
  cat("2. extraction skipped -- outlines from:\n   ", OUTLINE_SOURCE, "\n")
  load(OUTLINE_SOURCE)              # oto
  oto.new <- oto
  noerror <- names(oto.new)
  errors  <- character(0)
}


### III. 3rd step: Superimposition of outlines on the original picture to check their reliability

## Only specimens that have both a photograph and an extracted
## contour can be drawn; the extraction failures cannot.
jpgnames <- paste(intersect(names(oto.new), sub("\\.jpg$","",jpg.list)), ".jpg", sep="")

### This loop loads each otolith picture, extracts the coordinates of each outline generated before and  
### superimposes them on the photograph. Because the internal plot coordinates are different, we define
### X and Y coordinates for the outline so that it is aligned with the original picture.

write.superimposed.path = paste(fix.image.path,"Superimposed/",sep="/")
if(!file.exists(write.superimposed.path)){
  suppressWarnings(dir.create(write.superimposed.path))
}

sup.todo <- jpgnames
if (SKIP_EXISTING) sup.todo <- setdiff(sup.todo, list.files(write.superimposed.path, pattern="\\.jpg$"))
if (!RUN_SUPERIMPOSE) sup.todo <- character(0)
cat("3. superimpose:", length(sup.todo), "of", length(jpgnames), "specimens\n")
if (length(sup.todo)) {
  cat("   by year and subdivision:\n")
  print(table(sub(".*_(\\d{4})_SD(\\d+).*", "\\1 SD\\2", sup.todo)))
}

## Parallel: each worker writes its own file, nothing is shared.
sup.res <- foreach(nm=sup.todo, .errorhandling="pass",
                   .packages=c("magick","grDevices","graphics")) %dopar% {
  tryCatch({
    img     <- as.raster(magick::image_read(paste(fix.image.path, nm, sep="/")))
    otoname <- substr(nm, 1, nchar(nm)-4)
    pathi   <- paste(write.superimposed.path, otoname, ".jpg", sep="")

    coo <- as.data.frame(oto.new[[otoname]])
    x <- coo$V1 / (dim(img)[2]) + 1
    y <- coo$V2 / (dim(img)[1]) + 1

    jpeg(filename=pathi, width=dim(img)[2], height=dim(img)[1], units="px")
    plot(1:2, type='n')
    rasterImage(img, 1, 1, 2, 2)
    points(x, y, col="red", cex=cex)
    dev.off()
    TRUE
  }, error = function(e) conditionMessage(e))
}
sup.ok <- vapply(sup.res, isTRUE, logical(1))
if (length(sup.todo)) cat("   drawn:", sum(sup.ok), "of", length(sup.todo), "\n")
if (any(!sup.ok)) {
  print(data.frame(file=sup.todo[!sup.ok], message=unlist(sup.res[!sup.ok])))
}
message("3. Superimposed pictures are up to date! (",
        length(list.files(write.superimposed.path, pattern="\\.jpg$")), " in folder)")

stopCluster(cl)

if (length(errors) == 0){
  message("Otolith extraction finished! Any errors?\n--> No errors occurred!")
} else {
  message("Otolith extraction finished! Any errors?\n--> Yes, see file:")
  print(errors)
}
## Only a real extraction may overwrite oto.RData. A run that just
## redrew control images must not touch the file the analysis reads.
if (RUN_EXTRACT) {
  ind = names(oto.new)[!names(oto.new)%in%noerror]
  oto.new <- oto.new[!names(oto.new) %in% ind]
  oto <- c(oto,oto.new)
  save(oto, file = F_OTO)
  message("oto.RData written: ", length(oto), " outlines")
} else {
  message("Extraction was skipped -- oto.RData left untouched.")
}
