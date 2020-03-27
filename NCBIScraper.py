#!/usr/bin/env python
# coding: utf-8

# In[58]:


from Bio import Entrez
from Bio.Entrez import esearch, esummary, read, efetch
from pathos.multiprocessing import Pool
from lxml import etree as ET
import pandas as pd
import gc
import numpy as np
from p_tqdm import p_map
import string
import re

# **What we need**
# 
# *Pubmed*
# - [x] publication date
# - [x] PubmedID
# - [x] Authors
# - [x] Affiliations
# - [x] Keywords
# - [x] MeSH Descriptors
# - [x] MeSH Qualifiers
# - [x] MeSH hierarchy
# - [x] Abstract available?
# - [x] Abstract
# - [ ] Free Full Text available? (not in esummary, check via PMC)
# - [ ] Abstract/Fulltext Scraper for Editorials and stuffs


# In[1236]:


class NCBITool(object):
    
  
    def __init__(self, db='pubmed'):
        self.db = db
        print('Setting Database: ', db)
        
        parsers = {
            'pubmed':NCBITool.parse_pubmed_ids,
            'pmc':NCBITool.parse_pmc_ids
        }
    
        self.parser = parsers[self.db]
        
    def search(self,
               search_term,
               retmax=2000):
        
        handle = esearch(db=self.db, term=search_term, retmax = retmax)
        
        self.search_result = read(handle)
        self.ids = [str(i) for i in self.search_result.get('IdList')]
        
        print('No. of results: ', len(self.ids))
        if len(self.ids)==retmax:
            print('###! There might be more results available !###')

	return self
            
    def parse(self, n_jobs=12, ids=None, keep_xml=False):
        
        if not hasattr(self, 'search_result') and ids==None:
            print('###! No ids Provided - Aborting !###')
            return
        
        elif ids==None:
            ids = self.ids
                        
        if n_jobs > len(self.ids):
            n_jobs = len(self.ids)
            
        with Pool(processes=n_jobs) as pool:
            data =pool.map(
                    self.parser,
                    list([list(chunk) for chunk in np.array_split(ids, n_jobs)]),
                )
            
#         data = p_map(self.parser,
#                      [list(chunk) for chunk in np.array_split(ids, n_jobs)],
#                      num_cpus=n_jobs)

        column_names = {
            'pubmed':[
                'pubmed_id',
               'title',
               'authors',
               'affiliations',
               'pub_date',
               'abstract',
               'doi',
               'pmcid',
               'journal',
               'pubmed_class',
               'pubmed_keywords',
               'mesh_descriptors',
               'mesh_qualifiers'
            ],
            'pmc':[
                'pmcid',
                'title',
                'authors',
                'affiliations',
                'pub_date',
                'pubmed_id',
                'doi',
                'abstract',
                'full_text',
                'journal',
                'pmcclass',
                'pmc_keywords'
            ]
        }

        self.data = pd.concat(
            [
                pd.DataFrame(
                    x[0],
                    columns = column_names[self.db]
                )
                for x in data
            ]
        ).reset_index()
        
        self.data.pub_date = pd.to_datetime(self.data.pub_date)
	
	return self
    
    @classmethod
    def parse_pmc_ids(self, pmcid, retmode = 'xml'):
        
        #get xml from pmc
        handle = efetch(db='pmc', id = pmcid, retmode = retmode)
        xml_string = handle.read()
        xml = ET.fromstring(xml_string)

        #check for keywords and MeshTerms
        keys = []
        
        for art in xml.getchildren():
            
            #title
            title = ''.join(art.xpath('.//article-meta//article-title/text()'))
            
            #authors
            auth = zip(art.findall('.//*[@contrib-type="author"]/name/given-names'),
                       art.findall('.//*[@contrib-type="author"]/name/surname'))
            auth = ';'.join([' '.join([i.text, j.text]) for (i,j) in auth])
            
            #affiliations
            aff = ';'.join(art.xpath('.//aff/text()'))
            
            #publication_date
            pub_date = '-'.join(art.xpath('.//article-meta/pub-date[@pub-type="epub"]/*/text()')[::-1])
            
            #ids
            ##pubmed_id
            pubmed_id = ''.join(art.xpath('.//article-meta/article-id[@pub-id-type="pmid"]/text()'))
            ##doi
            doi = ''.join(art.xpath('.//article-meta/article-id[@pub-id-type="doi"]/text()'))
            ##pmcid
            pmcid = ''.join(art.xpath('.//article-meta/article-id[@pub-id-type="pmc"]/text()'))
            
            #abstract
            abstract = ''.join(art.xpath('.//abstract//*/text()'))
            
            #fulltext
            full_text = ''.join(art.xpath('.//body//*/text()'))
            
            #journal
            journal = ''.join(art.xpath('//journal-meta/journal-id[@journal-id-type="iso-abbrev"]/text()'))
            
            #pmcclass
            pmcclass = ''.join(art.xpath('.//article-meta/article-categories//subject/text()'))

            #pmc_keywords
            pmc_keywords = art.xpath('.//kwd/text()')
        
            keys.append(
                (
                    pmcid,
                    title,
                    auth,
                    aff,
                    pub_date,
                    pubmed_id,
                    doi,
                    abstract,
                    full_text,
                    journal,
                    pmcclass,
                    pmc_keywords
                )
            )
            
        return keys, xml_string
        
    
    @classmethod
    def parse_pubmed_ids(cls, pub_id, retmode='xml'): 
    
        #get xml from pubmed
        handle = efetch(db='pubmed', id = pub_id, retmode = retmode)
        xml_string = handle.read()
        xml = ET.fromstring(xml_string)

        #check for keywords and MeshTerms
        keys = []

        for art in xml.getchildren():
            
            #authors
            auth = zip(art.findall('.//Author/ForeName'), art.findall('.//Author/LastName'))
            auth = ';'.join([' '.join([i.text, j.text]) for (i,j) in auth])
            
            #affiliations
            email_re = re.compile('Electronic address.*\.')
            aff = ';'.join(
                        np.unique(
                            [
                                email_re.sub('',aff.text).strip(string.punctuation+' ')
                                for aff in art.findall('.//Affiliation')
                            ]
                        )
                    )   
            
            #publication date
            pub_date = '-'.join(
                [
                    x.text for x in [
                        art.find('.//PubMedPubDate[@PubStatus="entrez"]/Year'),
                        art.find('.//PubMedPubDate[@PubStatus="entrez"]/Month'),
                        art.find('.//PubMedPubDate[@PubStatus="entrez"]/Day')
                    ]
                ]
            )
            
            if pub_date == '':
                pub_date = '1900-01-01'
                
            pubmed_id = ''.join(art.xpath('.//PubmedData/ArticleIdList/ArticleId[@IdType="pubmed"]/text()'))
            doi = ''.join(art.xpath('.//PubmedData/ArticleIdList/ArticleId[@IdType="doi"]/text()'))
            pmcid =  ''.join(art.xpath('.//PubmedData/ArticleIdList/ArticleId[@IdType="pmc"]/text()'))
            title = ''.join(art.xpath('.//ArticleTitle/text()')).strip()
            
            abstract = ''.join(art.xpath('.//AbstractText/text()'))
            journal = ''.join(art.xpath('.//Journal/ISOAbbreviation/text()'))
            pubmed_class = ','.join(art.xpath('.//PublicationType/text()'))
            keys.append((pubmed_id,
                         title,
                         auth,
                         aff,
                         pub_date,
                         abstract,
                         doi,
                         pmcid,
                         journal,
                         pubmed_class,
                         art.xpath('MedlineCitation/KeywordList/*/text()'),
                         art.xpath('MedlineCitation/MeshHeadingList/MeshHeading/*[@MajorTopicYN="Y"]/../DescriptorName/text()'),
                         art.xpath('MedlineCitation/MeshHeadingList/MeshHeading/*[@MajorTopicYN="Y"]/../QualifierName/text()')))


#         df = pd.DataFrame(keys, columns=['pubmed_id', 'pub_date', 'keywords', 'mesh_descriptors', 'mesh_qualifiers'])
#         df.pub_date = pd.to_datetime(df.pub_date)

        gc.collect()

        return keys, xml_string

