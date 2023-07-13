args = commandArgs(trailingOnly = T)

rmarkdown::render(args[1], output_dir ="../../../../html_local")

# if (!file.exists("index.Rmd")){
#   file.create("index.Rmd")
# }
# 
# if (!file.exists("_site.yml")){
#   file.create("_site.yml")
# }
# 
# rmarkdown::render_site(args[1])