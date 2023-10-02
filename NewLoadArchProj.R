NewLoadArchProj <- function (path = "./", force = FALSE, showLogo = FALSE, orig.path = "/vf/users/CRGGH/heustonef/hpapdata/cellranger_snATAC/") 
{
  .newvalidInput(input = path, name = "path", valid = "character")
  .newvalidInput(input = force, name = "force", valid = "boolean")
  .newvalidInput(input = showLogo, name = "showLogo", valid = "boolean")
  path2Proj <- file.path(path, "Save-ArchR-Project.rds")
  if (!file.exists(path2Proj)) {
    stop("Could not find previously saved ArchRProject in the path specified!")
  }
  ArchRProj <- recoverArchRProject(readRDS(path2Proj))
  ArchRProj@projectMetadata$outputDirectory <- "./"
  outputDir <- getOutputDirectory(ArchRProj)
  outputDirNew <- outputDir
  ArrowFilesNew <- file.path(outputDirNew, "ArrowFiles", basename(ArchRProj@sampleColData$ArrowFiles))
  if (!all(file.exists(ArrowFilesNew))) {
    stop("ArrowFiles do not exist in saved ArchRProject!")
  }
  ArchRProj@sampleColData$ArrowFiles <- ArrowFilesNew
  if (length(ArchRProj@peakAnnotation) > 0) {
    keepAnno <- rep(TRUE, length(ArchRProj@peakAnnotation))
    for (i in seq_along(ArchRProj@peakAnnotation)) {
      if (!is.null(ArchRProj@peakAnnotation[[i]]$Positions)) {
        if (tolower(ArchRProj@peakAnnotation[[i]]$Positions) != 
            "none") {
          PositionsNew <- gsub(orig.path, outputDirNew, 
                               ArchRProj@peakAnnotation[[i]]$Positions, perl = FALSE)
          print(PositionsNew)
          if (!all(file.exists(PositionsNew))) {
            if (force) {
              keepAnno[i] <- FALSE
              message("Positions for peakAnnotation do not exist in saved ArchRProject!")
            }
            else {
              stop("Positions for peakAnnotation do not exist in saved ArchRProject!")
            }
          }
          ArchRProj@peakAnnotation[[i]]$Positions <- PositionsNew
        }
      }
      if (!is.null(ArchRProj@peakAnnotation[[i]]$Matches)) {
        MatchesNew <- gsub(orig.path, outputDirNew, 
                           ArchRProj@peakAnnotation[[i]]$Matches)
        if (!all(file.exists(MatchesNew))) {
          if (force) {
            message("Matches for peakAnnotation do not exist in saved ArchRProject!")
            keepAnno[i] <- FALSE
          }
          else {
            stop("Matches for peakAnnotation do not exist in saved ArchRProject!")
          }
        }
        ArchRProj@peakAnnotation[[i]]$Matches <- MatchesNew
      }
    }
    ArchRProj@peakAnnotation <- ArchRProj@peakAnnotation[keepAnno]
  }
  if (!is.null(getPeakSet(ArchRProj))) {
    if (!is.null(metadata(getPeakSet(ArchRProj))$bgdPeaks)) {
      bgdPeaksNew <- gsub(orig.path, outputDirNew, metadata(getPeakSet(ArchRProj))$bgdPeaks)
      if (!all(file.exists(bgdPeaksNew))) {
        if (force) {
          message("BackgroundPeaks do not exist in saved ArchRProject!")
          metadata(ArchRProj@peakSet)$bgdPeaks <- NULL
        }
        else {
          stop("BackgroundPeaks do not exist in saved ArchRProject!")
        }
      }
      else {
        metadata(ArchRProj@peakSet)$bgdPeaks <- bgdPeaksNew
      }
    }
  }
  ArchRProj@projectMetadata$outputDirectory <- outputDirNew
  message("Successfully loaded ArchRProject!")
  if (showLogo) {
    .ArchRLogo(ascii = "Logo")
  }
  ArchRProj
}
.newvalidInput <- function (input = NULL, name = NULL, valid = NULL) 
{
  valid <- unique(valid)
  if (is.character(valid)) {
    valid <- tolower(valid)
  }
  else {
    stop("Validator must be a character!")
  }
  if (!is.character(name)) {
    stop("name must be a character!")
  }
  if ("null" %in% tolower(valid)) {
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  av <- FALSE
  for (i in seq_along(valid)) {
    vi <- valid[i]
    if (vi == "integer" | vi == "wholenumber") {
      if (all(is.numeric(input))) {
        cv <- min(abs(c(input%%1, input%%1 - 1)), na.rm = TRUE) < 
          .Machine$double.eps^0.5
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "null") {
      cv <- is.null(input)
    }
    else if (vi == "bool" | vi == "boolean" | vi == "logical") {
      cv <- is.logical(input)
    }
    else if (vi == "numeric") {
      cv <- is.numeric(input)
    }
    else if (vi == "vector") {
      cv <- is.vector(input)
    }
    else if (vi == "matrix") {
      cv <- is.matrix(input)
    }
    else if (vi == "sparsematrix") {
      cv <- is(input, "dgCMatrix")
    }
    else if (vi == "character") {
      cv <- is.character(input)
    }
    else if (vi == "factor") {
      cv <- is.factor(input)
    }
    else if (vi == "rlecharacter") {
      cv1 <- is(input, "Rle")
      if (cv1) {
        cv <- is(input@values, "factor") || is(input@values, 
                                               "character")
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "palette") {
      cv <- all(.isColor(input))
    }
    else if (vi == "timestamp") {
      cv <- is(input, "POSIXct")
    }
    else if (vi == "dataframe" | vi == "data.frame" | vi == 
             "df") {
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
    }
    else if (vi == "fileexists") {
      cv <- all(file.exists(input))
    }
    else if (vi == "direxists") {
      cv <- all(dir.exists(input))
    }
    else if (vi == "granges" | vi == "gr") {
      cv <- is(input, "GRanges")
    }
    else if (vi == "grangeslist" | vi == "grlist") {
      cv <- .isGRList(input)
    }
    else if (vi == "list" | vi == "simplelist") {
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
    }
    else if (vi == "bsgenome") {
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text = input))
      }, error = function(e) {
        FALSE
      })
      cv <- any(cv1, cv2)
    }
    else if (vi == "se" | vi == "summarizedexperiment") {
      cv <- is(input, "SummarizedExperiment")
    }
    else if (vi == "seurat" | vi == "seuratobject") {
      cv <- is(input, "Seurat")
    }
    else if (vi == "txdb") {
      cv <- is(input, "TxDb")
    }
    else if (vi == "orgdb") {
      cv <- is(input, "OrgDb")
    }
    else if (vi == "bsgenome") {
      cv <- is(input, "BSgenome")
    }
    else if (vi == "parallelparam") {
      cv <- is(input, "BatchtoolsParam")
    }
    else if (vi == "archrproj" | vi == "archrproject") {
      cv <- is(input, "ArchRProject")
    }
    else {
      stop("Validator is not currently supported by ArchR!")
    }
    if (cv) {
      av <- TRUE
      break
    }
  }
  if (av) {
    return(invisible(TRUE))
  }
  else {
    stop("Input value for '", name, "' is not a ", paste(valid, 
                                                         collapse = ","), ", (", name, " = ", class(input), 
         ") please supply valid input!")
  }
}
