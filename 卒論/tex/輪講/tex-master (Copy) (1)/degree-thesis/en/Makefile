LATEXMK   = latexmk
TEXINPUTS = TEXINPUTS='.//;'
OUTPUT    = ./out
TITLE     = thesis

.PHONY: default
default: build

.PHONY: build
build:
	mkdir -p $(OUTPUT)
	$(TEXINPUTS) $(LATEXMK) -jobname=$(OUTPUT)/$(TITLE) main.tex

.PHONY: watch
watch: init
	$(TEXINPUTS) $(LATEXMK) -pvc

.PHONY: clean
clean:
	$(TEXINPUTS) $(LATEXMK) -C
	rm -r $(OUT)
