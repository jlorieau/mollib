.PHONE: build clean test
.DEFAULT_GOAL := build

CLI_SRC := $(wildcard *.sh)
CLI_TXT := $(CLI_SRC:.sh=.txt)
CLI_HTML := $(CLI_SRC:.sh=.html)
CLI_TEX := $(CLI_SRC:.sh=.tex)
CLI_RST := $(CLI_SRC:.sh=.rst)
CLI := $(CLI_SRC:.sh=)

errors = 0

#%.rst: %.sh

%.txt: %.sh  ## Process the text output of commands
	@echo $<
	@exec ./$< > $@  # export the text file

%.html: %.sh
	@echo "  $$ `cat $<`" > tmp_$@
	@FORCE_COLOR=TRUE exec ./$< >> tmp_$@
	@pygmentize -l shell-session -f html tmp_$@ | sed 's/^/    /g' >> tmp2_$@

	@echo ".. only:: html\n\n.. raw:: html\n" > $@
	@cat -e tmp2_$@ | sed 's/\$$$$//g'>> $@                            # Copy command output, strip non-print chars
	@echo '\n' >> $@

	@sed -i '' 's/\^\[\[1m/<font style=\"font-weight:bold;\">/g' $@  # Replace bold codes
	@sed -i '' 's/\^\[\[22m/<\/font>/g' $@                           # Replace closing codes for bold
	@sed -i '' 's/\^\[\[36m/<font color=\"#008b8b\">/g' $@           # Replace cyan colors
	@sed -i '' 's/\^\[\[91m/<font color=\"red\">/g' $@               # Replace red colors
	@sed -i '' 's/\^\[\[92m/<font color=\"green\">/g' $@             # Replace green colors
	@sed -i '' 's/\^\[\[33m/<font color=\"#abb51f\">/g' $@
	@sed -i '' 's/\^\[\[94m/<font color=\"blue\">/g' $@              # Replace blue colors
	@sed -i '' 's/\^\[\[95m/<font color=\"magenta\">/g' $@           # Replace magenta colors
	@sed -i '' 's/\^\[\[96m/<font color=\"#008b8b\">/g' $@           # Replace cyan colors
	@sed -i '' 's/\^\[\[0m/<\/font>/g' $@                            # Replace closing codes for colors

	@rm tmp_$@ tmp2_$@  # Remove temporary files

%.tex: %.sh
	@FORCE_COLOR=TRUE exec ./$< >> tmp_$@
	@sed -i '' 's/^/  /g' tmp_$@                                # Add a 2-spacetab at the start of every line

	@echo ".. only:: latex\n\n.. raw:: latex\n" > $@
	@echo "  \\\\begin{sphinxVerbatim}[commandchars=\\\\\\\\\\\\{\\\\},fontsize=\\\\footnotesize]" >> $@
	@echo "  \\\\textcolor{darkorange}{$$} `cat $<`" >> $@
	@cat -e tmp_$@ | sed 's/\$$//g' >> $@                       # Copy command output, strip non-print chars
	@echo "  \\\\end{sphinxVerbatim}\n {} \n" >> $@

	@sed -i '' 's/--/-{-}/g' $@                             # Preserve --
	@sed -i '' 's/\^\[\[1m/\\textbf{/g' $@                 # Replace bold codes
	@sed -i '' 's/\^\[\[22m/}/g' $@                        # Replace closing bold codes
	@sed -i '' 's/\^\[\[36m/\\textcolor{cyan}{/g' $@       # Replace cyan colors
	@sed -i '' 's/\^\[\[91m/\\textcolor{red}{/g' $@        # Replace red colors
	@sed -i '' 's/\^\[\[92m/\\textcolor{olivegreen}{/g' $@ # Replace green colors
	@sed -i '' 's/\^\[\[33m/\\textcolor{darkyellow}{/g' $@ # Replace dark yelow colors
	@sed -i '' 's/\^\[\[94m/\\textcolor{blue}{/g' $@       # Replace blue colors
	@sed -i '' 's/\^\[\[95m/\\textcolor{magenta}{/g' $@    # Replace magenta colors
	@sed -i '' 's/\^\[\[96m/\\textcolor{cyan}{/g' $@       # Replace cyan colors
	@sed -i '' 's/\^\[\[0m/}/g' $@                         # Replace closing codes for colors

	@rm tmp_$@  # Remove temporary files

%.rst: %.html %.tex
	@cat $? > $@
	@rm $?

test: build  ## Test the CLI output against the reference .txt files
	@for i in $(CLI); do \
		echo "Testing $$i.sh" ; \
		exec ./$$i.sh | diff $$i.txt - ; \
		return_value=$$? ; \
	done

clean:
	@rm -f *.txt *.tex *.html *.rst

build: $(CLI_TXT) $(CLI_RST)  ## Build the CLI .txt outputs and documentation .rst
	@echo 'Building CLI files'
	@find . -size  0 -delete  # delete empty files

help:  ## Print this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'