# Contributing to NuMAD

Thank you for considering contributing to NuMAD. Please follow these
guidelines for contributing.

## Reporting bugs

This section guides you through the process of [submitting a bug
report](https://github.com/sandialabs/NuMAD/issues/new) for NuMAD.
Following these guidelines will help maintainers understand your issue,
reproduce the behavior, and develop a fix in an expedient fashion. Before
submitting your bug report, please perform a cursory
search to see if the problem has been already reported. If it has been reported, and the
issue is still open, add a comment to the existing issue instead of opening a
new issue.

### Tips for effective bug reporting

- Use a clear descriptive title for the issue

- Describe the steps to reproduce the problem, the behavior you observed after
  following the steps, and the expected behavior

- Provide the SHA ID of the git commit that you are using

- For runtime errors, provide a function call stack

## Contributing code and documentation changes

Contributions can take the form of bug fixes, feature enhancements,
documentation updates. All updates to the repository are managed via [pull
requests](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests).
One of the easiest ways to get started is by looking at [open
issues](https://github.com/sandialabs/NuMAD/issues) and contributing fixes,
enhancements that address those issues. If your code contribution involves large
changes or additions to the codebase, we recommend opening an issue first and
discussing your proposed changes with the core development team to ensure that
your efforts are well directed, and so that your submission can be reviewed and
merged seamlessly by the maintenance team.

### Guidelines for preparing and submitting pull-requests

- Use a clear descriptive title for your pull-requests

- Describe if your submission is a bugfix, documentation update, or a feature
  enhancement. Provide a concise description of your proposed changes. 
  
- Provide references to open issues, if applicable, to provide the necessary
  context to understand your pull request
  
- Make sure that your pull-request merges cleanly with the `main` branch of
  NuMAD. When working on a feature, always create your feature branch off of
  the latest `main` commit
  
- Ensure that the code compiles without warnings, (leave for later? the unit tests and regression
  tests all pass without errors, and the documentation builds properly with your
  modifications)
  
- New physics models and code enhancements should be accompanied with relevant
  updates to the documentation, supported by necessary verification and
  validation, as well as unit tests and regression tests
  
  
Once a pull-request is submitted you will iterate with NuMAD maintainers
until your changes are in an acceptable state and can be merged in. You can push
addditional commits to the branch used to create the pull-request to reflect the
feedback from maintainers and users of the code.


### Coding conventions

Start reading our code and you'll get the hang of it. We optimize for readability:

- We indent using four spaces (soft tabs)
- We ALWAYS put spaces after list items and method parameters (`[1, 2, 3]`, not `[1,2,3]`), around operators (`x += 1`, not `x+=1`), and around hash arrows.
- This is open source software. Consider the people who will read your code, and make it look nice for them.


