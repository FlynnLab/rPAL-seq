# Purpose: convert canis family style to match homo/mus
infile  <- "Family_TP_MDCK.csv"
outfile <- sub("\\.csv$", "_converted.csv", infile)

x <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("family" %in% names(x))

# Keep a copy (optional)
orig_family <- x$family

# SNORD -> "small nucleolar RNA, C/D box <rest>"
x$family <- sub("\\bSmall nucleolar RNA\\s+SNORD\\s*([[:alnum:]/._-]+)",
                "small nucleolar RNA, C/D box \\1",
                x$family, perl = TRUE)

# SNORA -> "small nucleolar RNA, H/ACA box <rest>"
x$family <- sub("\\bSmall nucleolar RNA\\s+SNORA\\s*([[:alnum:]/._-]+)",
                "small nucleolar RNA, H/ACA box \\1",
                x$family, perl = TRUE)

# Neaten whitespace
x$family <- gsub("\\s+", " ", trimws(x$family))

write.csv(x, outfile, row.names = FALSE)

# Quick report
n_cd  <- sum(grepl("\\bsmall nucleolar RNA, C/D box\\b", x$family))
n_hac <- sum(grepl("\\bsmall nucleolar RNA, H/ACA box\\b", x$family))
cat(sprintf("Wrote %s\nConverted to C/D box: %d\nConverted to H/ACA box: %d\n",
            outfile, n_cd, n_hac))

