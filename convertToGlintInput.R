#!/usr/bin/env Rscript

# order is: datafile varname transpose
args <- commandArgs(trailingOnly = TRUE)

usagemsg = "USAGE:\nconvertToGlintInput.R <datafile> <varname>(optional) <transpose>(optional).\n<varname> should be provided in order to use <transpose>; alternatively, varname can be set to NULL and the script will try to find the variable name automatically.\n"

# 1. extract arguments
if (length(args) == 0){
    # cat is like print but with newline
    cat(usagemsg)
    quit()
}
if (length(args) == 1){
    datafile <- args[1]
    varname <- NULL
    transpose <- FALSE
}
if (length(args) == 2){
    datafile <- args[1]
    varname <- args[2]
    transpose <- FALSE
}
if (length(args) == 3){
    datafile <- args[1]
    varname <- args[2]
    transpose <- args[3]
}
if (length(args) > 3){
    cat(usagemsg)
    quit()
}

# 2. validate arguments

# datafile - check Rdata file exists
if (!file.exists(datafile)){
    print(paste("File", datafile, "does not exist."))
}

# transpose - if user specifyed transpose - check it is a boolean
if (typeof(transpose) == "character") { 

    if (toupper(transpose) %in% c('TRUE', 'FALSE')) {# all boolean options
        transpose <- type.convert(transpose)
    } else {
        print(paste("Not a boolean value:", transpose,"(booleans: true, false)"))
        quit()
    }
}
# varname - 
if(!is.null(varname)) {
    # if user specified NULL varname or numeric varname
    if (is.numeric(type.convert(varname))) {
        print(paste("varname should be a string:", varname))
        quit()
    } else if (!is.null(varname) && toupper(varname) %in% c('NULL')){
        varname <- NULL
    }
}


# 3. Start run - load Rdata file

#print(paste("Found datafile", datafile))

#if (!is.null(varname)){
#    print(paste("got argument name", varname))
#}


print(paste("Converting data file", datafile,'...'))
load(datafile)

# 4. find data argument
if (!is.null(varname)){
    if (varname %in% ls()){
        data <- get(varname)
        print(paste("Found variable", varname, "in data file", datafile))
    } else {
        print(paste("Cannot find variable", varname, "in data file", datafile))
        quit()
    }
} else {
    all_frame_args <- ls()[sapply(mget(ls()), is.data.frame)]
    all_matrix_args <- ls()[sapply(mget(ls()), is.matrix)]

    if(length(all_frame_args) != 1){
        if (length(all_matrix_args) != 1) { # there is no data frame and no matrix in datafile
            print("Cannot find only a single data frame or matrix in the data file. Please execute the script again and specify variable name.")
            quit()
        }
        else { # there is only data frame in datafile
            print(paste("Found matrix", all_matrix_args[1]))
            data <- get(all_matrix_args[1])
        }
    } else {
        if (length(all_matrix_args) == 1) { # there is both data frame and no matrix in datafile
            print(paste("Found data frame ", all_frame_args[1], " and matrix ", all_matrix_args[1], ".\nRdata must contain only one data variable. Otherwise, please run this script again and specify varname."))
            quit()
        }
        else { # there is only matrix in datafile
            print(paste("Found data frame ", all_frame_args[1]))
            data <- get(all_frame_args[1])
        }
    }
}

# 5. transpose if asked
if(toupper(transpose) == "TRUE"){
    print("transposing data...")
    data <-t(data)
}


# 6. save output
output_filename <- paste(datafile, ".txt", sep='') # do not add .glint extenstion since glint will think it's commpressed glint data file (which is not)
print(paste("Data file was saved into", output_filename))
cat("ID", file=output_filename, append=FALSE, sep = "")
write.table(data, output_filename, na = "NaN", sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE, append=TRUE) # sep must be something but space or tabs since there is no name for index [0][0] and glint won't be able to read it. col.name = NA is for allowing [0][0] to be "" (so there will be something there)

