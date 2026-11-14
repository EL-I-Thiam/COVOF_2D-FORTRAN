PANDOC ?= $(shell which pandoc 2>/dev/null || echo /opt/homebrew/bin/pandoc)
PDF_ENGINE ?= xelatex
PANDOC_FLAGS = --from=markdown+tex_math_dollars+pipe_tables+grid_tables+fenced_code_attributes \
               --pdf-engine=$(PDF_ENGINE) \
               --number-sections \
               --toc

BUILD_DIR = build

all: user_guide dev_guide tutorials check_pages

user_guide:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume1_user_guide/00_metadata.yaml volume1_user_guide/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Guide_Utilisateur_CFD.pdf
	cat volume1_user_guide/[0-9]*.md > user_guide.md

dev_guide:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume2_dev_guide/00_metadata.yaml volume2_dev_guide/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Guide_Developpement_CFD.pdf
	cat volume2_dev_guide/[0-9]*.md > dev_guide.md

tutorials:
	mkdir -p $(BUILD_DIR)
	$(PANDOC) volume3_tutorials/00_metadata.yaml volume3_tutorials/[0-9]*.md \
		$(PANDOC_FLAGS) -o $(BUILD_DIR)/Recueil_Tutoriels_CFD.pdf
	cat volume3_tutorials/[0-9]*.md > tutorials_guide.md

check_pages:
	@echo "=== Vérification du nombre de pages (Minimum requis : 40 pages/volume) ==="
	@python3 ./check_pages.py $(BUILD_DIR)/Guide_Utilisateur_CFD.pdf 40
	@python3 ./check_pages.py $(BUILD_DIR)/Guide_Developpement_CFD.pdf 40
	@python3 ./check_pages.py $(BUILD_DIR)/Recueil_Tutoriels_CFD.pdf 40

clean:
	rm -rf $(BUILD_DIR) user_guide.md dev_guide.md tutorials_guide.md
