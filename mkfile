
all:V: paper.pdf supplement.pdf

%.tex: %.md
    pandoc \
        --standalone \
        --from=markdown+simple_tables \
        --template=template.tex \
        --to=latex \
        $prereq | ./fix-tables.pl > $target

%.pdf: %.tex
    xelatex $prereq

