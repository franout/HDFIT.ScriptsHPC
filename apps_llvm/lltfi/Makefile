SHELL = /bin/bash

LLVM_SRC_ROOT ?= $(HOME)/llvm-project
LLVM_DST_ROOT ?= $(HOME)/llvm-project/install

lltfi_REPO = https://github.com/DependableSystemsLab/LLTFI
lltfi_COMMIT = 2a2dc0e0d4747e4dc23c501a9701224e48331d15
hdfit_PATCH = 0001-Integrating-HDFIT.patch

.DEFAULT_GOAL := lltfi
.PRECIOUS: %/.cloned %/.patched %/.compiled
.PHONY: lltfi

all: lltfi clean veryclean

clean:
	-cd lltfi && rm -rf build

veryclean:
	rm -rf lltfi

lltfi/.cloned:
	@echo "Cloning $(@D)..."
	git clone $($(@D)_REPO) $(@D)
	cd $(@D) && git checkout $($(@D)_COMMIT)
	chmod -R u=rwX,go=rX $(@D)
	@touch $(@D)/.cloned

lltfi/.patched: lltfi/.cloned
	@echo "Patching $(@D)..."
	cd $(@D) && git am ../$(hdfit_PATCH)
	@touch $(@D)/.patched

lltfi: lltfi/.patched
	@echo "Compiling $(@)..."
	cd $(@) && rm -rf build
	cd $(@) && ./setup -LLFI_BUILD_ROOT $(CURDIR)/lltfi/build -LLVM_SRC_ROOT $(LLVM_SRC_ROOT) -LLVM_DST_ROOT $(LLVM_DST_ROOT)
