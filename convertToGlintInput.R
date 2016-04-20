#!/usr/bin/env Rscript

# order is: datafile argname transpose
args <- commandArgs(trailingOnly = TRUE)

# 1. extract arguments
if (length(args) == 0){
    # cat is like print but with newline
    cat("USAGE IS:\nconvertToGlintInput.R <datafile> <argname(optional)> <transpose(optional)>.\n<argname> must be set in order to set <transpose>, you can set it to NULL and we will try to find the name automatically\n")
    quit()
}
if (length(args) == 1){
    datafile <- args[1]
    argname <- NULL
    transpose <- FALSE
}
if (length(args) == 2){
    datafile <- args[1]
    argname <- args[2]
    transpose <- FALSE
}
if (length(args) == 3){
    datafile <- args[1]
    argname <- args[2]
    transpose <- args[3]
}
if (length(args) > 3){
    cat("USAGE IS:\nconvertToGlintInput.R <datafile> <argname(optional)> <transpose(optional)>.\n<argname> must be set in order to set <transpose>, you can set it to  NULL and we will try to find the name automatically\n")
    quit()
}

# 2. validate arguments

# datafile - check Rdata file exists
if (!file.exists(datafile)){
    print(paste("file", datafile, "does not exists"))
}

# transpose - if user specifyed transpose - check it is a boolean
if (typeof(transpose) == "character") { 
    if (toupper(transpose) %in% c('TRUE', 'FALSE')) {# all boolean options
        transpose <- type.convert(transpose)
    } else {
        print(paste("not a boolean value", transpose,"(booleans: true, false)"))
        quit()
    }
} else if (is.numeric(type.convert(transpose))) {
    print(paste("not a boolean value", transpose,"(booleans: true, false)"))
     quit()
}

# argname - 
# if user specified NULL argname or numeric argname
if (is.numeric(type.convert(argname))) {
    print(paste("argname is not a string value", argname))
    quit()
} else if (!is.null(argname) && toupper(argname) %in% c('NULL')){
    argname <- NULL
}


# 3. Start run - load Rdata file

print(paste("got datafile", datafile))

if (!is.null(argname)){
    print(paste("got argument name", argname))
}


print(paste("converting datafile", datafile,'...'))
load(datafile)

# 4. find data argument
if (!is.null(argname)){
    if (argname %in% ls()){
        data <- get(argname)
    } else {
        print(paste("cant find argument ", argname, "in datafile", datafile))
        quit()
    }
} else {

    all_frame_args <- ls()[sapply(mget(ls()), is.data.frame)]
    # dfs <- ls()[sapply(get(ls()), is.list)]
    if(length(all_frame_args) != 1){
        print("cant find data in datafile, please execute script again and specify argument name")
        quit()
    } else {
        data <- get(all_frame_args[1])
    }
}

# 5. transpose if asked
if(transpose){
    print("transposing data...")
    data <-t(data)
}


# 6. save output
output_filename <- paste("output", datafile, sep='_') # do not add .glint extenstion since glint will think it's commpressed glint data file (which is not)
print(paste("data is saved to", output_filename, "as glint format"))
write.table(data, output_filename, na = "NaN", sep = "\t", col.names = FALSE, row.names = FALSE) # TODO: when adding names dont forget the total name (data[0,0])
