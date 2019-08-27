# MedicalCorpusMaker
Create a custom medical corpus using PubMed API.

This has been created to enable development of a custom corpus for machine learning / natural language processing. 

This is based on searching a large medical index - [PubMed] (https://www.ncbi.nlm.nih.gov/pubmed) for custom search terms.

The abstracts are then collated to return either a Counter object with all terms and their frequency, or just the whole corpus. This is selected by setting count_list to TRUE (default is FALSE)

An email address must be provided to enable API searching of PubMed via. [BioPython](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html).

Multiple terms can be used, separated by commas.

```
make_corpus(email, terms*, count_list=FALSE)
```

Example use:

```
email = email@account.com
term_1 = 'renal'
term_2 = 'muscle'

corpus = make_corpus(email, term_1, term_2)
#returns corpus 

corpus_counted = make_corpus(email, term_1, term_2, count_list = TRUE)
#returns Counter object

```
All/any changes & pull requests welcomed. 

```
requirements:
Biopython

```
