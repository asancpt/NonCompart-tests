pdf: 
	Rscript -e "rmarkdown::render('validation.Rmd', output_format = 'bookdown::pdf_document2', output_file = 'validation.pdf', encoding = 'UTF-8')"
