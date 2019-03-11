# STANDARDS

- Classes/methods > functions
- verbose variable annotateSmudges
- snakecase variables and camel case for methods
- logging on error stream using `logging` module

# Repository

- `master` branch is a production ready branch. Only well tested code from `dev` branch (and no other) should be merged into master. Majority of the commits to `master` will be a tagged release too. Every merge should be done though pull request.
- `dev` branch is the branch that should be forked for implementation of additional features. Only in cases of one commit features the feature should be committed directly to `dev`. Use rebase before merging and don't use fast forward (I like to keep track of what was done where).
- other branches are feature specific