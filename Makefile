MS_PATH=../OpenMS/bin
XTANDEM=../tandem-linux-13-09-01-1/bin/64-bit/all_static/tandem.exe

LOGS=logs
DATA=data
INIS=inis

FASTA=$(DATA)/uniprot-mouse_taxonomy_10090_keyword_181_20140226_CON_REV.fasta
MZML=$(DATA)/20131122010914_20131111_CO_0418FrMe_QBIC_R19_02_QFMHS051X6_sorted.mzML

KGRN=\x1B[32m
KRED=\e[91m
RESET=\033[0m

# Run an OpenMS tool
#
# This macro takes 2 mandatory and one optional argument:
#   $(1) : The name of the OpenMS tool to run
#   $(2) : The arguments of the OpenMS tool
#   $(3) : An identifier for this particular instance of the tool to
#          differenciate between targets that use the same tool
#
# Do not specify -ini as an option. This is added by the macro.
#
# Usage:
# 	@$(call MSTOOL,FileMerger, -in $^ -out $@,targetid)
define MSTOOL
	$(MS_PATH)/$(1) \
		$(2) \
		-ini $(INIS)/$(1)$(3).ini \
	1> $(LOGS)/$(1)$(3).stdout \
	2> $(LOGS)/$(1)$(3).stderr
endef

# Check if the last command exited with error code 0.
#
# Print `done` if it did and `error` if it did not.
define CHECK
	if [ $$? -eq 0 ] ;\
	then \
		echo -e "$(KGRN)done$(RESET)" ;\
	else \
		echo -e "$(KRED)error$(RESET)" ;\
	fi
endef

.SECONDARY:

all: $(MZML:.mzML=.quantified.protein_groups.csv)

clean :
	@rm -f work/*
	@rm -f logs/*

%.merged.mzML : $(MZML)
	@printf %-70s "Merging files..."
	@$(call MSTOOL,FileMerger, -in $^ -out $@ )
	@$(CHECK)

%.protein_groups.idXML : %.peptides.probs.indexed.idXML
	@printf %-70s "Infer protein groups from peptides..."
	@$(call MSTOOL,FidoAdapter, -in $^ -out $@)
	@$(CHECK)

%.peptides.idXML : %.mzML
	@printf %-70s "Infer peptides from spectra..."
	@$(call MSTOOL,XTandemAdapter, -in $^ -out $@ -xtandem_executable $(XTANDEM) -database $(FASTA) -threads 3)
	@$(CHECK)

%.probs.idXML : %.idXML
	@printf %-70s "Convert scores to probabilities..."
	@$(call MSTOOL,IDPosteriorErrorProbability, -in $^ -out $@ -prob_correct)
	@$(CHECK)

%.indexed.idXML : %.idXML
	@printf %-70s "Find proteins that contain peptides..."
	@$(call MSTOOL,PeptideIndexer, -in $^ -out $@ -fasta $(FASTA) -decoy_string _rev -annotate_proteins -enzyme:specificity none)
	@$(CHECK)

%.quantified.protein_groups.csv : %.mapped_features_ids.featureXML %.protein_groups.idXML
	@printf %-70s "Quantify the protein groups"
	@$(call MSTOOL,ProteinQuantifier, -in $< -out $@ -protein_groups $(word 2,$^))
	@$(CHECK)

%.peaks.mzML : %.mzML
	@printf %-70s "Look for peaks..."
	@$(call MSTOOL,PeakPickerHiRes,-in $^ -out $@)
	@$(CHECK)

%.featureXML : %.peaks.mzML
	@printf %-70s "Look for features in peak list..."
	@$(call MSTOOL,FeatureFinderCentroided, -in $^ -out $@)
	@$(CHECK)

%.mapped_features_ids.featureXML : %.featureXML %.peptides.probs.indexed.idXML
	@printf %-70s "map features to peptide ids..."
	@$(call MSTOOL,IDMapper, -in $< -out $@ -id $(word 2,$^))
	@$(CHECK)
