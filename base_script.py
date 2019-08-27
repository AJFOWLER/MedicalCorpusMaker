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
from Bio import Entrez
from Bio import Medline


def make_corpus(email, *terms, count_list = False):
    #this is the main function and creates a corpus of words based on terms provided
    #variable length of terms so join with AND (ignored if only 1 term)
    phrase = 'AND'.join(terms)
    #get search results
    search_result = pubmed_search(email, phrase)
    #extract these records
    all_text = extract_records(search_result)
    #if count_list then return counted version of the corpus.
    if(count_list):
        corpus = counted(all_text)
    else:
        corpus = all_text
    return(corpus)

def pubmed_search(email, phrase):
    '''login & search pubmed'''
    Entrez.email = email
    #email login to enable API access
    handle = Entrez.esearch(db='pubmed', retmax = 100000, term = phrase)
    #find records from pubmed based on our phrase
    records = Entrez.read(handle)
    #parse the handle from esearch
    return(records)

def extract_records(records):
    import string
    #take ID list
    ids = records['IdList']
    #return medline entries for these Ids
    handle_full = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode = "text")
    #fetch full records from pubmed based on these.
    cleaned_list = Medline.parse(handle_full)
    #as list.
    cleaned_list = list(cleaned_list)
    #set up all text
    all_text = []
        
    for entry in cleaned_list:
        #go over each entry
        abstract = entry.get('AB')
        all_text.append(abstract)
        #append abstract to the all_text list.
    cleaned = ''.join(item for item in all_text if item)
    #only join items that are actually items (i.e. drop none_types)
    cleaned = ' '.join(l for l in cleaned if l not in string.punctuation)
    #drop those items that aren't strings
    return(cleaned)
        
    #may need to remove numbers/./formulae
    
def counted(all_text):
    '''return a counted list for probability distribution estimation'''
    from collections import Counter
    #split corpus
    corpus_splitted = all_text.split()
    #return counter object of the corpus.
    counted_list= Counter(corpus_splitted)
    return(counted_list)
