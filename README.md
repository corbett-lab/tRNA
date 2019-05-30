# tRNA gene family simulator

This program is a simulation of tRNA gene family evolution. The theory behind this is largely based on adaptation of ideas from this paper:
http://ped.fas.harvard.edu/files/ped/files/nature97_0.pdf
into a finite, diploid population with customizable parameters.

To compile, simply type 'make' while in the directory after downloading. This should work with most environments. If you want to compile using tcmalloc (which is not available on, for example, UCSC's hummingbird cluster), you can edit the makefile to comment out that line.

To run, simply type ./tRNA

## List of command line flags:

##### -n: number of diploid individuals in population, will be held constant throughout simulation (default: 10,000)
##### -g: number of generations for simulation to run (default: 10,000,000)
##### --start: number of tRNA genes in the genome initially (default: 1)
##### -m: length of genome in morgans (default: 30.0)

##### --print: generation printing interval (--print 500 will print every 500 generations, etc.)
##### -b: burn-in - output will start printing after a set number of generations pass (-b 0 prints results starting from beginning, -b 50000 prints only after 50000 generations have passed)
##### --path: the path to the directory where the functionDists directory is located (necessary for reading in mutation effects)
##### -s: seed, to ensure non-unique results on many simulations begun at the same time
##### --run: run name, to name folders for optional output of each run
##### --sample: whether or not you want to sample individuals from the population (default: false)
##### --sample-freq: sample every ____ generations (default: 10,000)
##### --sample-count: sample ____ individuals each time (default: 10)
##### --quiet: use if you are not interested in printing every tRNA's stats every time (default: false)
##### --output-lifespans: outputs lifespan log file, which contains lifespan data for all tRNAs created during the simulation (default: false)

##### --ug: baseline germline mutation rate (default: 1e-6)
##### --us: baseline somatic mutation rate (default: 1.96e-5)
##### --del: deletion rate (default: 3.63e-6)
##### --dup: duplication rate (default: 3.72e-6)
##### -c: fraction of duplications that are local (default: 0.6)
##### --segdup: fraction of local duplications that are segmental duplications (default: 0.25)
##### --geneconv: rate at which non-allelic gene conversion occurs (default: 2.5e-7)
##### --somatic-coefficient: for use in establishing relationship between expression and somatic mutation rate (default: 11.8898)
##### --max-mutations: number of mutations allowed before a tRNA gene is no longer considered a tRNA gene (default: 7)
##### --pseudo: begin with a tRNA pseudogene (0.0 sequence and 0.0 expression) in addition to the number of functional tRNA genes (default: false)

##### --function: fitness function to be used (default: "redundant"; other options: "model", "gaussian", "exponential")
##### --min-fitness: cumulative expression of tRNAs required for fitness 1.0 under "redundant" function (default: 4.0)
##### --fitmean: cumulative expression of tRNAs rqeuired for fitness 1.0 under "gaussian" function (default: 10.52110756)
##### --fitsd: standard deviation of gaussian fitness function (default: fitmean/4.0)
##### --fitlambda: lambda value for exponential fitness function (default: -15.0)

#### The following are based on the Nowak paper linked above:

##### --model-1:
Two tRNAs initially, with exactly the same mutation rate and function. All mutations are completely inactivating. No somatic mutations, duplications or deletions are possible.

##### --model-2:
Two tRNAs, one with function = 1 and mutation rate = --ug; one with function = 0.8 and mutation rate = --ug / 100. All mutations are completely inactivating. No somatic mutations, duplications or deletions are possible.

##### --model-4:
Any number of tRNAs, each with function = 1 and mutation rate == --ug, but somatic mutation rates are now possible. When invoking this model, individual fitness = 1 - ((developmental error rate)^(number of functional tRNA genes)). See Nowak paper for more detailed explanation.
###### --model-4-count: number of tRNA genes in genome initially in addition to --start value (only used with --model-4; default = 1).
###### --model-4-deverr: developmental error rate (only used with --model-4; default = 1e-4).

All questions should be directed to bthornlo@ucsc.edu.
