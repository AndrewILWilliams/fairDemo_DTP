# fairDemo_DTP
Notebooks and binder instance for the Oxford DTP Physical Climate workshop 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AndrewWilliams3142/fairDemo_DTP/main)

Access this Binder by clicking the blue badge above or at the following URL:

https://mybinder.org/v2/gh/AndrewWilliams3142/fairDemo_DTP/main

## Notes
The `environment.yml` file should list all Python libraries on which your notebooks
depend, specified as though they were created using the following `conda` commands:

```
conda activate example-environment
conda env export --from-history -f environment.yml
```

Note that the only libraries available to you will be the ones specified in
the `environment.yml`, so be sure to include everything that you need! 

Also note that if you skip the `--from-history`, conda may include OS-specific
packages in `environment.yml`, which you would have to manually prune from
`environment.yml`.  For example, confirmed macOS-specific packages that should
be removed are:

* libcxxabi=4.0.1
* appnope=0.1.0
* libgfortran=3.0.1
* libcxx=4.0.1
