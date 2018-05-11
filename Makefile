all: validation.pdf validation.md validation.html

validation.pdf: 
	Rscript -e "rmarkdown::render('validation.Rmd', output_format = 'bookdown::pdf_document2', output_file = 'validation.pdf', encoding = 'UTF-8')"

validation.md: 
	Rscript -e "rmarkdown::render('validation.Rmd', output_format = 'github_document', output_file = 'validation.md', encoding = 'UTF-8')" ;\
	rm validation.html

validation.html: validation.pdf
	Rscript -e "rmarkdown::render('validation.Rmd', output_format = 'bookdown::tufte_html2', output_file = 'validation.html', encoding = 'UTF-8')" ;\
	cp validation.html docs/index.html

clean:
	rm validation.pdf validation.md validation.html
