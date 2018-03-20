readme: 
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'bookdown::pdf_document2', output_file = 'README.pdf', encoding = 'UTF-8')"
