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