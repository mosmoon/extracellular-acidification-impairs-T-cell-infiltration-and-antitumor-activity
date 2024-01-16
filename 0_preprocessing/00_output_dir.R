# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: make output dir
# ---------------------------------------------

loaddir <- function() {
    if(!dir.exists("../output")) dir.create("../output")
    if(!dir.exists("../output/00.RawCount")) dir.create("../output/00.RawCount")
}

