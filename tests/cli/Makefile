.PHONE: build clean test
.DEFAULT_GOAL := build

CLI_SRC := $(wildcard *.sh)
CLI_TXT := $(CLI_SRC:.sh=.txt)
CLI_HTML := $(CLI_SRC:.sh=.html)
CLI_TEX := $(CLI_SRC:.sh=.tex)
CLI_RST := $(CLI_SRC:.sh=.rst)

errors = 0

#%.rst: %.sh

%.txt: %.sh  ## Process the text output of commands
	echo $<
	exec ./$< > $@  # export the text file
	FORCE_COLOR=TRUE exec ./$< > tmp_$@.txt

%.html: %.sh
	@echo "$ `cat $<`\n" > tmp_$@
	@FORCE_COLOR=TRUE exec ./$< >> tmp_$@
	@echo '.. only:: html\n\n.. raw:: html\n' > $@
	@pygmentize -l shell-session -f html tmp_$@ | sed 's/^/    /g' >> $@
	@echo '\n' >> $@

	@sed -i 's/\$$//g' $@                                         # Remove $ at the end of lines
	@sed -i 's/\^\[\[1m/<font style=\"font-weight:bold;\">/g' $@  # Replace bold codes
	@sed -i 's/\^\[\[22m/<\/font>/g' $@                           # Replace closing codes for bold
	@sed -i 's/\^\[\[36m/<font color=\"#008b8b\">/g' $@           # Replace cyan colors
	@sed -i 's/\^\[\[91m/<font color=\"red\">/g' $@               # Replace red colors
	@sed -i 's/\^\[\[92m/<font color=\"green\">/g' $@             # Replace green colors
	@sed -i 's/\^\[\[33m/<font color=\"#abb51f\">/g' $@
	@sed -i 's/\^\[\[94m/<font color=\"blue\">/g' $@              # Replace blue colors
	@sed -i 's/\^\[\[95m/<font color=\"magenta\">/g' $@           # Replace magenta colors
	@sed -i 's/\^\[\[96m/<font color=\"#008b8b\">/g' $@           # Replace cyan colors
	@sed -i 's/\^\[\[0m/<\/font>/g' $@                            # Replace closing codes for colors

	@rm tmp_$@  # Remove temporary files

%.tex: %.sh
	@echo "$ `cat $<`\n" > tmp_$@
	@FORCE_COLOR=TRUE exec ./$< >> tmp_$@
	@echo ".. only:: latex\n\n.. raw:: latex\n" > $@
	@echo "    \\\\begin{sphinxVerbatim}[commandchars=\\\\\\\\\\\\{\\\\},fontsize=\\\\footnotesize]" >> $@
	@cat -e tmp_$@ >> $@
	@echo "    \\\\end{sphinxVerbatim}\n {} \n" >> $@

#	@sed -i 's/^/    /g' $@                                  # Add a tab at the start of every line
#	@sed -i 's/^   //g' $@                                   # Strip leading spaces
	@sed -i 's/\$$//g' $@                                    # Remove $ at the end of lines
	@sed -i 's/--/-{-}/g' $@                                 # Preserve --
	@sed -i 's/\^\[\[1m/\\\\textbf{{/g' $@                   # Replace bold codes
	@sed -i 's/\^\[\[22m/}}/g' $@                            # Replace closing bold codes
	@sed -i 's/\^\[\[36m/\\\\textcolor{{cyan}}{{/g' $@       # Replace cyan colors
	@sed -i 's/\^\[\[91m/\\\\textcolor{{red}}{{/g' $@        # Replace red colors
	@sed -i 's/\^\[\[92m/\\\\textcolor{{olivegreen}}{{/g' $@ # Replace green colors
	@sed -i 's/\^\[\[33m/\\\\textcolor{{darkyellow}}{{/g' $@ # Replace dark yelow colors
	@sed -i 's/\^\[\[94m/\\\\textcolor{{blue}}{{/g' $@       # Replace blue colors
	@sed -i 's/\^\[\[95m/\\\\textcolor{{magenta}}{{/g' $@    # Replace magenta colors
	@sed -i 's/\^\[\[96m/\\\\textcolor{{cyan}}{{/g' $@       # Replace cyan colors
	@sed -i 's/\^\[\[0m/}}/g' $@                             # Replace closing codes for colors

	@rm tmp_$@  # Remove temporary files

%.rst: %.html %.tex
	@cat $? > $@
	@rm $?

test: $(CLI_SRC) $(CLI_TXT) ## Test the CLI output
	@echo "Testing $<"
	@exec ./$< | diff $(word 2,$^) -

	@$(eval return_value = $(shell echo $$? ))

	$(eval errors=$(shell echo $$(($(errors)+$(return_value)))))
	@echo "$(errors) mismatches found."

clean:
	@rm *.txt *.tex *.html *.rst

build: $(CLI_TXT) $(CLI_RST)
	@echo 'Building CLI files'
	@find . -size  0 -delete  # delete empty files
