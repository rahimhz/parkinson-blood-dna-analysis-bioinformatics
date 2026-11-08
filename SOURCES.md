# Sources and data notes

## GEO accession

- GSE165081 record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165081
- Official supplementary directory: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165081/suppl/

## Files used

- `GSE165081_M.csv.gz` — processed methylation M-value matrix.
- `GSE165081_Sample_Sheet.PD-CN.csv.gz` — sample metadata.

The pipeline does not assume actual headers, sample counts, IDs, labels, or pair IDs. It inspects the data at run time and stops on ambiguous structure or invalid sample matching.

## Processed and raw data

The workflow uses the processed M-value matrix above. `GSE165081_RAW.tar` contains raw IDAT data and is deliberately not downloaded, extracted, or processed by this implementation.

## Design and interpretation

The GEO description mentions matched controls. This does not establish that the sample metadata contains verified one-to-one pairs. Paired analysis is permitted only with a supplied pair-ID column that passes validation. Independent analysis is explicit and exploratory. Neither analysis establishes causality or diagnostic utility.

## Citation

Henderson AR, Wang Q, Meechoovet B, et al. **DNA Methylation and Expression Profiles of Whole Blood in Parkinson's Disease.** *Frontiers in Genetics*. 2021;12:640266. doi:10.3389/fgene.2021.640266. PMID: 33981329.

- [PubMed record](https://pubmed.ncbi.nlm.nih.gov/33981329/)
- [GEO accession GSE165081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165081)

The publication describes a cohort design of 15 age-matched subjects per group, while the GSE165081 methylation subseries and its downloaded sample sheet contain 28 samples: 14 PD and 14 Control. This project reports the validated subseries sample count, not the broader study design count.
