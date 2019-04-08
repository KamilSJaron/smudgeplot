# STANDARDS

- Classes/methods > functions
- verbose variable annotateSmudges
- snakecase variables and camel case for methods
- logging on error stream using `logging` module

# Repository

- `master` branch is a production ready branch. Only well tested code from `dev` branch (and no other) should be merged into master. Majority of the commits to `master` will be a tagged release too. Every merge should be done though pull request.
- `dev` branch is the branch that should be forked for implementation of additional features. Only in cases of one commit features the feature should be committed directly to `dev`. Use rebase before merging and don't use fast forward (I like to keep track of what was done where).
- other branches are feature specific

# Building of dependency list

Building up the list of required packagess

```
virtualenv -p python3 venv
source venv/bin/activate
pip3 install -r requirements.txt
pip3 install PySAIS==1.0.4
python3 setup.py install
# repeat following two lines until everything is fine
smudgeplot -v
pip install <package_name>
# get the list of installed packages without this one
pip freeze | grep -v "smudgeplot" > requirements.txt
deactivate
rm -r venv
```