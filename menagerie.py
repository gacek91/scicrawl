import os, os.path as op
from glob import glob
import warnings, git
from github import Github
import github


def new_project(name, path = os.getcwd(), folders = ['Data','Scripts','Images','Results','Articles','Presentations'],
               subfolders = None, rproj = False, pyproj = False,
               overwrite = False, gitinit = True):
    '''
    Create a new project directory together with pre-defined folders.
    
    Parameters
    ----------
    name: str
        project name
    path: str
        project path [default is the working directory]
    folders: list
        list of subfolders to be created in sthe project [default: Data, Scripts, Images, Results]
    subfolders: dict
        dictionary of folders to be created in subfolders [default is None]
    rproj: bool
        automatically create a new r project [default is False]
    pyproj: bool
        automatically create a new python project [default is False]
    overwrite: bool
        overwrite the project if already exists [default is False] - warning, you will lose all the files if you use this argument.
    gitinit: bool
        initialize a git repository at the given path [default is True]
    
    Returns
    ----------
    Information whether or not the project has successfully been created.
    
    '''
    
    import os, os.path as op
    from glob import glob
    import warnings, git
    from github import Github
    import github
    create = list()
    fullpath = op.join(path, name)
    if overwrite == True:
        from shutil import rmtree
        if op.exists(fullpath):
            rmtree(fullpath)
            print('Removed old project together with all files.')
    if not op.exists(fullpath):
        create.append(fullpath)
        folderspaths = [op.join(fullpath, folder) for folder in folders]
        create.extend(p for p in folderspaths)
        if subfolders != None:
            subfolderspaths = list()
            for folder in subfolders.keys():
                if type(subfolders[folder]) == list:
                    subfolderspaths.extend([op.join(fullpath,folder,p) for p in subfolders[folder]])
                else:
                    subfolderspaths.append(op.join(fullpath,folder,subfolders[folder]))
            create.extend(p for p in subfolderspaths)
        for c in create:
            os.mkdir(c)
        print('Successfully created the project <<{}>> with all the requested folders and subfolders.'.format(name))
    else:
        print('Project <<{}>> already exists. Doing nothing.'.format(name)) 
        ### dodaj tworzenie projektu w R
    if rproj == True:
        rprojcode = ['Version: 1.0', '','RestoreWorkspace: Default','SaveWorkspace: Default','AlwaysSaveHistory: Default','',
                     'EnableCodeIndexing: Yes','UseSpacesForTab: Yes','NumSpacesForTab: 2','Encoding: UTF-8','',
                     'RnwWeave: Sweave','LaTeX: pdfLaTeX']
        file = open(op.join(fullpath,"{}.Rproj").format(name),"w+")
        for line in rprojcode:
            file.write(line)
            file.write('\n')
        file.close()
        print('Created .Rproj file.')
    if gitinit:
        r = git.Repo.init(path = fullpath)
        print("Git repo initiated.")


def githubauth(token):
    '''
    Authenticate with a token on GitHub.
    
    Parameters
    ----------
    token: str
        GitHub token
    
    Returns
    ----------
    Instance: github.MainClass.Github
        Logged GitHub instance
    Link: str
        Partial GitHub remote link
    '''
    from github import Github
    from functools import partial
    g = Github(token)
    Instance = g.get_user()
    Link = partial('https://{userlogin}:{token}@github.com/{userlogin}/{repo}.git'.format)
    Link = Link.func(userlogin = Instance.login, token=token, repo='{repo}')
    return Instance, Link
    
def create_github_repo(Instance, name):
    '''
    Create a new repository on GitHub
    
    Parameters
    ----------
    Instance: github.MainClass.Github
        Logged GitHub instance
    name: str
        Repository name
        
    Returns
    ----------
    If repo exists, returns an error and does nothing.
    If repo does not exist, creates a new repo.
        
    '''
    from warnings import warn
    try:
        repo = Instance.create_repo(name)
        print('Repo {} created: https://github.com/{}/{}'.format(name, Instance.login, name))
    except:
        warn('Repo {} already exists (https://github.com/{}/{}). Doing nothing.'.format(name, Instance.login, name))
        
def overwrite_github_repo(Instance, name):
    '''
    Overwrrite an existing repository on GitHub
    
    Parameters
    ----------
    Instance: github.MainClass.Github
        Logged GitHub instance
    name: str
        Repository name
        
    Returns
    ----------
    If repo does not exist, returns an error and does nothing.
    If repo does exist, overwrites it creates a new repo.
        
    '''
    from warnings import warn
    x = ''
    i = 0
    while x != 'y' or x != 'n':
        if i == 0:
            x = str(input("You are about to delete an existing repository. There is no going back. Are you sure you want to proceed? (y [overwrite] /n [do nothing])"))
            i +=1
        elif i > 0:
            x = str(input("Nothing else is going to work - y or n please!"))
        if x == 'n':
            warn("Procedure halted. Doing nothing.")
            break
        if x == 'y':
            try:
                repo = Instance.get_repo(name)
                repo.delete()
                repo = Instance.create_repo(name)
                print('Repo {} overwritten: https://github.com/{}/{}'.format(name, Instance.login, name))
                break
            except:
                warn('\nRepo {} does not exist. Use create_github_repo function to create a repository.'.format(name))
                break

def git_commit_all(githubauth, name, path, commit_message = 'Just a commit'):
    '''
    Commit all changes to a GitHub repository
    
    Parameters
    ----------
    githubauth: tuple or str
        Tuple with a logged GitHub instance and a partial GitHub remote link (use githubauth function) or a string - a GitHub authentication token.
    name: str
        Repository name
    path: str
        Local path to the repo
    
    '''
    from warnings import warn
    from menagerie import githubauth as gthbauth
    from glob import glob
    import os.path as op
    import git
    
    fullpath = op.join(path, name)
    
    if type(githubauth) == tuple:
        try:
            Instance, Link = githubauth
        except:
            warn("\nNot a correct githubauth object. Use githubauth function or provide a github authentication token")
    elif type(githubauth) == str:
        try:
            Instance, Link = gthbauth(githubauth)
        except:
            warn("\nNot a correct githubauth object. Use githubauth function or provide a github authentication token")
    while type(githubauth) != tuple and type(githubauth) != str:
        warn("\nNot a correct githubauth object. Use githubauth function or provide a github authentication token")
        break
    while type(Instance) != github.AuthenticatedUser.AuthenticatedUser:
        warn('\nNo GitHub AuthenticatedUser instance! Use githubauth function or provide a correct instance!')
        break
    files =  glob(op.join(fullpath,'**'),recursive=True)
    try:
        r = git.Repo(fullpath,search_parent_directories=True)
    except:
        r = git.Repo.init(path = fullpath)
    
    [r.git.add(f) for f in file_list]
    r.index.commit(commit_message)
    try:
        r.git.push("--set-upstream", r.remotes.origin, r.head.ref)
        print('Changes committed to repo: https://github.com/{}/{}'.format(Instance.login, name))
    except:
        r.create_remote('origin', Link.format(repo=name))
        r.git.push("--set-upstream", r.remotes.origin, r.head.ref)
        print('Changes committed to repo: https://github.com/{}/{}'.format(Instance.login, name))