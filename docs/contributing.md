# Contributing to BRIDGES

!!! Warning
    **Just a draft.** On this page, mention the following:

    * Currently repo is not published and not open to external contributions.
    * Write a "How to get started" section. Point to Snakemake/Sherlock/MkDocs ... introductions.
    * Define steps that need to be done as part of pull requests (run tests, update documentation, maybe have a change log file to keep track of changes over versions)
    * Write something about using Git flow.
    * List coding conventions.

## Coding conventions

### Python and Julia

For Python and Julia code style conventions, check out the [Python Style Guide](https://peps.python.org/pep-0008/) and [Julia Style Guide](https://github.com/invenia/BlueStyle).

### Snakemake

Here, we want to collect some coding conventions for Snakemake syntax that we found helpful while developing BRIDGES. 

1. In this project, we prefer the use of scripts over shell commands in rules, as shell commands usually differ between Unix-like systems and Windows. In cases where the use of shell commands is absolutely necessary prefer bash syntax (c.f., section "A word on Windows and Mac").
2. Scripts used in rules should have the same name as the rule.
3. Use dot notation to specify rule inputs whenever possible (see `rule combine_datasets` in the example in "Introduction to Snakemake"). Providing plain file paths is prone to errors and makes later modification more cumbersome.