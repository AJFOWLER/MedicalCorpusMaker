#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:00:36 2019

@author: AJF

module: make medical corpus
logic:
    -take terms
    - search pubmed for those terms
    - parse into single corpus
    - option to have a counted list for a p_distribution etc.
"""

#this is the main function and creates a corpus of words based on terms provided

def make_corpus(email, *terms, count_list = False):    
    #variable length of terms
    phrase = 'AND'.join(terms)
    search_result = pubmed_search(email, phrase)
    all_text = extract_records(search_result)
    if(count_list):
        corpus = counted(all_text)
    else:
        corpus = all_text
    return(corpus)

def pubmed_search(email, phrase):
    from Bio import Entrez
    '''login to pubmed'''
    Entrez.email = email
    #email login to enable API access
    handle = Entrez.esearch(db='pubmed', retmax = 100000, term = phrase)
    #find records from pubmed based on our phrase
    records = Entrez.read(handle)
    #parse the handle from esearch
    return(records)

def extract_records(records):
    from Bio import Entrez
    from Bio import Medline
    import string
    
    ids = records['IdList']
    #take ID list
    
    handle_full = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode = "text")
    #fetch full records from pubmed based on these.
    cleaned_list = Medline.parse(handle_full)
    cleaned_list = list(cleaned_list)

    all_text = []
        
    for entry in cleaned_list:
        abstract = entry.get('AB')
        all_text.append(abstract)
        #append abstract to the all_text list.
    cleaned = ''.join(item for item in all_text if item)
    #only join items that are actually items (i.e. drop none_types)
    cleaned = ' '.join(l for l in cleaned if l not in string.punctuation)
    return(cleaned)
        
    #may need to remove numbers/./formulae
    
def counted(all_text):
    '''
    return a counted list for probability distribution estimation
    
    '''
    from collections import Counter
    
    corpus_splitted = all_text.split()
    counted_list= Counter(corpus_splitted)

    return(counted_list)
