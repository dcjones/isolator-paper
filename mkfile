
all:V: paper.pdf supplement.pdf

%.tex: %.md
    pandoc --standalone --from=markdown+simple_tables --to=latex -o $target $prereq

%.pdf: %.tex
    xelatex $prereq

