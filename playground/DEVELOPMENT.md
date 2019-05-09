# STANDARDS

- spaces around operators
- snake_case for variables and functions
- camelCase for classes and methods
- verbose naming is more important than a detailed documentation
- R code is tested using `testthat` and python code in `dev` branch using `unittest`

## versioning

- we try to keep `master` branch clean (i.e. production ready code).
- `dev` branch should be also a working code, however here mistakes are permitted. Most of the development is done in subbranches of `dev` branch, once a feature is implemented, it should be merged to `dev`. I like to use `--no-ff` to keep a record, what was developed where (this might be a bad practice and if you know why, let me know).
- One rule I would like to keep is this. There must be at least 72 hours incubation period between commits merged into `dev` and merging `dev` into `master`. The reason is simple, if there are any mistakes that were not spotted, there is a change to catch them. Also it takes while to run all the tests and stuff (the travis.ci testing is not working yet, but I hope it will quite soon).

## language

The future is `python3`, and potentially in a low level language if we will ever need a speed-up.

R will be abolished as some point. It's just that some of the R libraries (for kernel smoothing) seems to work way faster than the pythonic ones and I have not really figured why. The pythonic smudgpelot is not there yet (although a lot of code is already written, check `dev` branch).