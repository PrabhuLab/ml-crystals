# Contributing to crystract

First off, thank you for considering contributing to `crystract`! We value and appreciate contributions from the community. Every little bit helps, and credit will always be given.

This document provides guidelines for contributing to the project. Please make sure to read the relevant section before making your contribution.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Pull Requests](#pull-requests)
- [Style Guides](#style-guides)
  - [R Code Style](#r-code-style)
  - [Commit Messages](#commit-messages)

## Code of Conduct

This project and everyone participating in it is governed by the [`crystract` Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainer.

## How Can I Contribute?

### Reporting Bugs

Bugs are tracked as [GitHub issues](https://github.com/PrabhuLab/ml-crystals/issues). Before opening a new issue, please perform a search to see if the problem has already been reported.

If you are reporting a new bug, please include as many details as possible to help us reproduce and fix it. A minimal reproducible example (MRE) is required. Your bug report should include:

1.  **A clear and descriptive title:** "Error when parsing CIF with multiple data blocks" is better than "CIF error".
2.  **A detailed description:** Explain the problem, what you expected to happen, and what actually happened.
3.  **A minimal reproducible example (MRE):** Provide a small piece of R code and a minimal example CIF file that reproduces the issue. This is the most important step!
4.  **Version information:** Include the output of `sessionInfo()`, which lists your R version, OS, and versions of attached packages.

### Suggesting Enhancements

Enhancement suggestions are also tracked as [GitHub issues](https://github.com/PrabhuLab/ml-crystals/issues). Before creating a new enhancement suggestion, please check the issue list to see if a similar idea has already been discussed.

When submitting an enhancement suggestion, please provide:

1.  **A clear and descriptive title.**
2.  **A detailed description of the proposed enhancement.** Explain the use case and why this feature would be valuable.
3.  **Example code (if possible):** Show how the new feature might work in practice.

### Pull Requests

We actively welcome your pull requests!

1.  **Fork the repository** and create your branch from `main`.
2.  **Set up your development environment.** It's recommended to use RStudio and have the `devtools` and `roxygen2` packages installed.
3.  **Make your changes.** Please follow the code style guide below.
4.  **Add tests.** If you add new functionality, please add a corresponding unit test in the `tests/testthat/` directory.
5.  **Update documentation.** If you change function arguments or behavior, please run `devtools::document()` to update the `NAMESPACE` and `.Rd` documentation files.
6.  **Ensure all checks pass.** Before submitting, run `devtools::check()` to make sure your changes haven't introduced any new problems.
7.  **Submit the Pull Request (PR)** against the `main` branch. Provide a clear description of the problem and your solution. Link to any relevant issues.

## Style Guides

### R Code Style

Please adhere to the [Tidyverse style guide](https://style.tidyverse.org/). You can use the `styler` package to automatically format your code.

-   Use `<-` for assignment, not `=`.
-   Use snake_case for variable and function names (e.g., `calculate_distances`).
-   Aim for a line length of 80 characters.

### Commit Messages

-   Use the present tense ("Add feature" not "Added feature").
-   Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
-   Limit the first line to 72 characters or less.
-   Reference issues and pull requests liberally in the body of the commit message.
