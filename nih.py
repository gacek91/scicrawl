import argparse
import sys
import os, os.path as op

pth = sys.argv[0]
os.chdir(op.dirname(pth))
parser = argparse.ArgumentParser(description='Search through NIH databases for articles matching your keywords and export them to excel.\
You must provide your e-mail.')
parser._action_groups.pop()

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-s', '--search', metavar='phrase',
                    type=str,required=True,action='append', nargs='*')
required.add_argument('-m', '--email', metavar='email', type=str, required=True, action='store', help='Your e-mail address')
optional.add_argument('-r','--records',metavar='number',help='Number of records to be analyzed; default is 100', type=int,
                    default=100,action='store')
optional.add_argument('-db', '--database', metavar='name',
                    type=str,help = 'Database: use pmc or pubmed. Default is pmc',default='pmc',action='store')

vars = parser.parse_args()
term = ' '.join([s for s in vars.search[0]])
print('Search phrase: {}'.format(term))
records = vars.records
print('Number of NIH records analyzed: {}'.format(records))
database = vars.database
email = vars.email
print('Database: {}'.format(database))

def nihsearch(term,email,records=100,database='pmc',return_as='df',fullOutput = False, save_excel=True):
    """
    types available = df, list

    todos - separate authors for connection networks
    remove keys from abstracts

    """
    from Bio import Entrez
    from Bio import Medline
    from tqdm import tqdm
    Entrez.email = email    # Always tell NCBI who you are
    handle = Entrez.esearch(db=database, term=term,retmax=records)
    record = Entrez.read(handle)
    handle = Entrez.efetch(db=database, id=record['IdList'], rettype="medline", retmode="text")
    record = Medline.parse(handle)
    export = []
    full_records = []
    i = 0
    indices = []
    for rec in tqdm(iterable=record, desc='Processed',unit='',total=records):
        currec = rec
        if fullOutput == True:
            full_records.append(currec)
        dickt = {}
        try:
            dickt['Authors'] = ', '.join(currec['AU'])
            doi = [d for d in currec['AID'] if 'doi' in d][0]
            doi = doi.split(' [doi]')[0]
            doi = 'https://doi.org/{}'.format(doi)
            dickt['DOI'] = doi
            dickt['Language'] = currec['LA'][0]
            if currec['LA'][0] == 'eng':
                dickt['Title'] = currec['TI']
            else:
                dickt['Title'] = currec['TI'].split('[')[1].split(']')[0]
            if len(currec['AB'].split(':')[0]) < 20:
                currec['AB'] = currec['AB'].split(currec['AB'][0:len(currec['AB'].split(':')[0])+2])[1]
            dickt['Abstract'] = currec['AB']
            dickt['Journal'] = currec['JT']
            try:
                dickt['Year'] = int(currec['DP'])
            except:
                dickt['Year'] = int(currec['DP'].split(' ')[0])
            if len(currec['PT']) > 1:
                dickt['Type'] = currec['PT'][1]
            else:
                dickt['Type'] = currec['PT'][0]
            export.append(dickt)
            indices.append(i)
        except:
            pass
        i += 1
    if return_as == 'df':
        import pandas as pd
        export = pd.DataFrame(export,index=indices)
        export = export[['Authors','Year','Title','Language','Type','Journal','DOI','Abstract']]
    elif return_as == list:
        pass
    else:
        print('Requested data type ({}) not available, exporting as a list of dictionaries'.format(return_as))
    if save_excel == True:
        export.to_excel('{}.xlsx'.format(term))
        print('Finished - {} articles in the final database'.format(export.shape[0]))
    if fullOutput == True:
        return export, full_records
    else:
        return export

nihsearch(term=term,records=records,database=database)
