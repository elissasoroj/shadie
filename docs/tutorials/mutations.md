:octicons-alert-16: Under Construction

# Mutation Types in `shadie`

Mutation types are defined in `shadie` as `MutationType` class objects, created with the `shadie.mtype()` function. 

`shadie.mtype()` requires a minimum of 3 arguments:

* dominance (float): dominance coefficient of the mutation for diploid
* distribution (str): distribution that the fitness effect will be selected from
* params (float): additional arguments that define the distribution - number of arguments will vary based on the distribution 

`shadie` defaults:
```py
NEUT = MutationType(0.5, "f", 0.0)		#neutral mutation
DEL = MutationType(0.1, "g", -0.03, 0.2)  #deleterious
BEN = MutationType(0.8, "e", 0.1)         #beneficial
```

## Tutorial
*Click here to download this tutorial in Jupyter Notebook form*

Try creating your own mutation type:
```py
mut1 = shadie.mtype(0.5, 'n', 0.5, 0.1)
mut1.inspect();
```

What's happening here? The first float argument is the **dominance** coefficient, which determines the behavior of mutations in a heterozygote. The next string argument determines the **type of distribution** the fitness effects will be drawn from. The number of arguments that follow depend on the distribution chosen:

```py
# "f" is fixed fitness effect (no distribution), so takes a single argument
mut1 = shadie.mtype(0.5, "f", 0.1)

# "g" is a gamma distribution and takes 2 arguments: alpha, beta
mut2 = shadie.mtype(0.1, "g", 0.3, 0.2)

# "n" is a normal distribution and takes 2 arguments: mean, 
mut3 = shadie.mtype(0.5, "n", 0.05, 0.1)
```

Note that the variable name you assign to your mutation type (e.g. "mut1") will only be recognized by your python environment, not internally by shadie. shadie will only recognize mutation types using the `shadie`-assigned `idx` when creating your genomic element types. 

The repr lists the settings for the mutation type you just created:

```py
mut2
```
`# <MutationType: m10, 0.1, g, (0.3, 0.2)>`

Notice that your mutation type now has an `idx` (in the format m#) - the first value listed in the repr. This is a unique id that `shadie` uses to keep track of mutations. It may not start at `m1`, as some default mutation types have already been defined by `shadie`. 

If you'd like more information about your mutation, you can use the `inspect()` function to visualize your mutation distribution and see the parameters listed explicitly:

```py
mut2.inspect()
```
