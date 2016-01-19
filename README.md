# A Biologically Constrained Basal Ganglia model

This code implements a mean-field model of the whole basal ganglia, whose parameterization was optimized to respect best a collection of numerous anatomical and electrophysiological data.

This model was used for the paper "A biologically constrained model of the whole basal ganglia addressing the paradoxes of connections and selection." available at [J Comput Neurosci](http://link.springer.com/article/10.1007/s10827-013-0476-2).

### Compiling

The code depends on the [Boost libraries](http://www.boost.org/) (1.37 or newer). 

To compile, type:

```sh
$ make
```

This creates an executable file called *single_run*.

### Reproduce the main results of the 2014 paper


Launch the script to run different parameterizations of the basal ganglia model and store the results in *model_output*:

```sh
$ bash tst_solutions.sh > model_output
```

The figures of the paper can be generated in R with the script *output_figures.R*:

```R
R> setwd('your-working-directory')
R> source('output_figures.R')
```

(don't forget to replace *your-working-directory* with the actual one)
