MS_PATH=../OpenMS/bin
XTANDEM=../tandem-linux-13-09-01-1/bin/64-bit/all_static/tandem.exe

LOGS=logs
DATA=data
INIS=inis

FASTA=$(DATA)/uniprot-mouse_taxonomy_10090_keyword_181_20140226_CON_REV.fasta
MZML=$(DATA)/20131122010914_20131111_CO_0418FrMe_QBIC_R19_02_QFMHS051X6_sorted.mzML

KGRN="\x1B[32m"
RESET="\033[0m"

.SECONDARY:

all: $(MZML:.mzML=.quantified.protein_groups.csv)
#all: $(MZML:.mzML=.indexed.idXML)

clean :
	@rm -f *.merged.mzML
	@rm -f *.idXML

%.merged.mzML : $(MZML)
	@printf %-70s "Merging files $^..."
	@$(MS_PATH)/FileMerger -in $^ -out $@ > $(LOGS)/FileMerger.log
	@echo -e $(KGRN) done $(RESET)

%.protein_groups.idXML : %.peptides.probs.indexed.idXML
	@printf %-70s "Infer protein groups from peptides $^ ..."
	@$(MS_PATH)/FidoAdapter -in $^ -out $@ > FidoAdapter.log
	@echo -e $(KGRN) done $(RESET)

%.peptides.idXML : %.mzML
	@printf %-70s "Infer peptides from spectra $^ ..."
	@$(MS_PATH)/XTandemAdapter -xtandem_executable $(XTANDEM) -in $^ -out $@ -database $(FASTA) -threads 3 > $(LOGS)/XTandemAdapter.log
	@echo -e $(KGRN) done $(RESET)

%.probs.idXML : %.idXML
	@printf %-70s "Convert scores to probabilities in file $^ ..."
	@$(MS_PATH)/IDPosteriorErrorProbability -in $^ -out $@ -prob_correct > $(LOGS)/IDPosteriorErrorProbability.log
	@echo -e $(KGRN) done $(RESET)

%.indexed.idXML : %.idXML
	@printf %-70s "Find proteins that contain peptides $^ ..."
	@$(MS_PATH)/PeptideIndexer -in $^ -out $@ -fasta $(FASTA) -decoy_string _rev -annotate_proteins -enzyme:specificity none > $(LOGS)/PeptideIndexer.log
	@echo -e $(KGRN) done $(RESET)

%.quantified.protein_groups.csv : %.mapped_features_ids.featureXML %.protein_groups.idXML
	@printf %-70s "Use $< to quantify the protein groups"
	@$(MS_PATH)/ProteinQuantifier -out $@ -in $< -protein_groups $(word 2,$^) > $(LOGS)/ProteinQuantifier.log
	@echo -e $(KGRN) done $(RESET)

%.peaks.mzML : %.mzML
	@printf %-70s "Look for peaks in $^..."
	@$(MS_PATH)/PeakPickerHiRes -in $^ -out $@ -ini $(INIS)/PeakPickerHiRes.ini > $(LOGS)/PeakPickerHiRes.log
	@echo -e $(KGRN) done $(RESET)

%.featureXML : %.peaks.mzML
	@printf %-70s "Look for features in peak list..."
	@$(MS_PATH)/FeatureFinderCentroided -in $^ -out $@ -ini $(INIS)/FeatureFinderCentroided.ini > $(LOGS)/FeatureFinderCentroided.log
	@echo -e $(KGRN) done $(RESET)

%.mapped_features_ids.featureXML : %.featureXML %.peptides.probs.indexed.idXML
	@printf %-70s "map features to peptide ids..."
	@$(MS_PATH)/IDMapper -id $(word 2,$^) -in $< -out $@ -ini $(INIS)/IDMapper.ini > $(LOGS)/IDMapper.log
	@echo -e $(KGRN) done $(RESET)
