installed.pkgs <- rownames(installed.packages())
required.pkgs <- readLines("/tmp/packages.txt")
install.me <- setdiff(required.pkgs, installed.pkgs)

BiocManager::install(install.me, update = FALSE, ask = FALSE)

