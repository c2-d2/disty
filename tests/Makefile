.PHONY: all clean

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

DM=../disty
SEQ=$(wildcard *.fa)

all:
	@for seq in $(SEQ); do \
		pref=$${seq%.*}; \
		for s in $$(seq 0 4); do \
			CMD="$(DM) -s $$s $$seq"; \
			EXP="exp.$${pref}.$$s.txt"; \
			OBS="obs.$${pref}.$$s.txt"; \
			if [[ -f "$$EXP" ]]; then \
				echo ; \
				echo "$$CMD"; \
				echo "======================================"; \
				echo ; \
				echo $$EXP; \
				echo $$OBS; \
				echo ; \
				$$CMD > "$$OBS" 2> /dev/null; \
				diff "$$OBS" "$$EXP"; \
				echo "OK"; \
				echo; \
			fi; \
		done; \
	done;

clean:
	rm -f obs.*.txt
