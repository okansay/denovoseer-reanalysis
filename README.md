# DeNovoSeer re-analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23590533.svg)](https://doi.org/10.5281/zenodo.23590533)

Code for the numbers in my letter on DeNovoSeer (Qiu et al., Briefings in Bioinformatics 2026, https://doi.org/10.1093/bib/bbag428).

I did not retrain anything. The script takes the variant-level test predictions the authors published (Supplementary File S1), joins them to the annotation files in their GitHub repository, and splits performance by where each label came from. In their study a label comes either from ClinVar or from a rule applied to Gene4Denovo observations. Those two groups behave very differently, and the pooled AUC hides that.

## What you need

Python 3.9 or newer, plus:

```
pip install -r requirements.txt
```

You also need Supplementary File S1. It is inside `supplementary_material_bbag428.zip` on the article page. I can't redistribute it, so download it yourself and point the script at the CSV.

The two Gene4Denovo files (annotation and ACMG evidence) are fetched from https://github.com/dadaxiaoxiao2333/DeNovoSeer the first time you run the script. They are about 175 MB together and end up in `data/`.

## Running it

```
python reanalysis.py --s1 Supplementary_File_S1_variant_predictions_Gene4Denovo_bbag428.csv
```

It takes a few minutes, mostly because of the bootstrap. Everything prints to the screen and is also written to `results/`:

- `summary.txt` has all the numbers quoted in the letter
- `auc_by_label_source.csv` has AUC (with 95% bootstrap CI) and AP for DeNovoSeer and a few single scores, by label group
- `rule_pathogenic_synonymous.csv` lists the synonymous test variants that got a pathogenic label from the Gene4Denovo rule
- `repeated_sample_ids.csv` lists rule-based pathogenic labels where every supporting record carries the same sample ID

## What the script checks

1. AUC and AP for all S1 variants, then separately for ClinVar-labelled and rule-labelled variants. Each variant is counted once, using the mean predicted probability S1 reports across the splits where it was in the test set. Pooled this way I get 0.876 and 0.881, the same as the paper.
2. How single inputs to the model do on their own in each group: the total clinical evidence score, PM2, BayesDel, AlphaMissense and CADD.
3. DeNovoSeer against AlphaMissense and CADD in the rule group, restricted to variants all three tools score, with paired bootstrap differences.
4. Synonymous variants labelled pathogenic by the rule, how often the model calls them pathogenic, and whether they had a SpliceAI score before imputation.
5. The PS2 and PM6 columns in the released ACMG evidence file.
6. Rule-based pathogenic labels where all the supporting records share one sample ID.

## Things to keep in mind

The S1 file only has labelled variants that showed up in at least one test split. Variants that appeared in several splits have their probabilities averaged in S1, which is not exactly the per-split mean ± SD reported in the paper. The pooled numbers still match, so I think the comparison is fair, but it is not the same procedure.

Bootstrap intervals use a fixed seed (2026) and 2,000 resamples. Your intervals may move in the third decimal if you change either.

Points 5 and 6 come from the full released Gene4Denovo files, not just the test set.

A repeated sample ID does not prove a label is wrong. Some of these variants may have other independent observations that were not recorded in the file. I report them because the labelling rule asks for two studies or two affected people, and these rows make it hard to tell which applies.

## Citation

Dolu KO. DeNovoSeer re-analysis code (v1.0). Zenodo; 2026. https://doi.org/10.5281/zenodo.23590533

## License

MIT. The DeNovoSeer data and predictions belong to their authors and are not included here.
