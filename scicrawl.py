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

def excel2bibtex(file):
    """
    Create bibtex entries for all articles in excel file. DOI addresses needed!

    """
    import urllib.request
    from urllib.error import HTTPError
    import pandas as pd
    from copy import deepcopy
    from tqdm import tqdm
    term = op.split(file)[-1].split('.')[0]
    df = pd.read_excel(file)
    dois = list(df.DOI)
    bibtex = []
    for doi in tqdm(iterable=dois,desc='Processed',unit='',total=len(dois)):
        dickt = {}
        req = urllib.request.Request(doi)
        req.add_header('Accept', 'application/x-bibtex')
        try:
            with urllib.request.urlopen(req) as f:
                bibtex.append(f.read().decode())
        except HTTPError as e:
            pass
    bibtex_orig = deepcopy(bibtex)
    bibtex = '\n'.join([b for b in bibtex])
    print(bibtex)
    with open("{}.bib".format(term), "w+") as output:
        output.write(bibtex)
