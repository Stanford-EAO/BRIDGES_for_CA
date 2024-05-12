# Access and Contributing

## Access to BRIDGES

The BRIDGES code repository in not publicly available at the moment. However, if you are interested in collaborating, please do not hesitate to [reach out](contact.md). An early version of the BRIDGES model was published [here](https://www.sciencedirect.com/science/article/pii/S266679242200004X) and provides a good first insight into the underlying equations.


## Contributing and Coding conventions

### Getting started

Once you have access to the BRIDGES repository, we recommend reading through the sections "BRIDGES' structure" and "Running BRIDGES". The "Introduction to" sections offer more details and tutorials for the tools we are using.

### Gitflow

In large parts, we follow the Gitflow branching model as methodology to develop our code. Find more information [here](https://nvie.com/posts/a-successful-git-branching-model/).

!!! ToDo
    Define steps that need to be done as part of pull requests (run tests, update documentation, maybe have a change log file to keep track of changes over versions)


### Python and Julia

For Python and Julia code style conventions, check out the [Python Style Guide](https://peps.python.org/pep-0008/) and [Julia Style Guide](https://github.com/invenia/BlueStyle).

### Snakemake

Here, we want to collect some coding conventions for Snakemake syntax that we found helpful while developing BRIDGES. 

1. In this project, we prefer the use of scripts over shell commands in rules, as shell commands usually differ between Unix-like systems and Windows. In cases where the use of shell commands is absolutely necessary prefer bash syntax (c.f., section "A word on Windows and Mac").
2. Scripts used in rules should have the same name as the rule.
3. Use dot notation to specify rule inputs whenever possible (see `rule combine_datasets` in the example in "Introduction to Snakemake"). Providing plain file paths is prone to errors and makes later modification more cumbersome.