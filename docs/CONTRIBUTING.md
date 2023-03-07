# Developer's Guidelines
This file contains guidelines and documentation for developers who want to contribute to this repository.

## Code Guidelines

This project uses the [_Black_ code style](https://github.com/psf/black) for consistency and clarity. Before committing,
format the project code by running `black .` in the project directory.  For docstrings and comments,
the [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) is used.

Some more guidelines to follow:
* Write a lot of [unit tests](https://docs.python.org/3/library/unittest.html). This catches errors close to the source and gives you confidence that your code works. If you're using PyCharm, create unit tests for a method by Right-click > Go To > Tests. From the console, all tests can be run using `python -m unittest` or using the script in `bin/tests.sh`.
* Use [type hints](https://docs.python.org/3/library/typing.html) and [docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for every method. This helps prevent errors and assists others with using your method properly.
* Ideally, only a single person works on a file at a time to prevent merge conflicts. This requires a certain file structure,
  avoiding long files and preferring small, specialized files.
* If you start using a package that is not yet in `requirements.txt`, add it with the specific version you are using, so your code works for everyone else as it does for you.
* If your code creates files that change with every execution or are user-specific, add them to the `.gitignore`.
* If you are writing a new feature extractor or clustering method, read the docstrings of the base classes carefully to implement your class correctly.



## Git Guidelines
[GitHub Flow](https://docs.github.com/en/get-started/quickstart/github-flow) is used as Git model for this repository. Please become familiar with it before committing. Follow these guidelines:
* The main branch should always be in a usable state.
* Never commit to the main branch directly (in fact, pushing to main is blocked). Always work on your feature branch, then use a pull request to merge your changes into main.
* Use pull requests as platform for discussions and questions. You can open a pull request even if your code is not done yet. Tagging people in pull requests to get feedback on work-in-progress code is explicitly encouraged. Once you're done, the pull request will be approved for merging into main.
* Always start new branch for new feature, do not reuse a branch for multiple features. A feature should be a single component of about a week's work, if more, split into smaller features. Always update main before starting a new feature branch, then start your feature branch from main.
* If it can be avoided, do not merge feature branches into another. This leads to messy pull requests with much manual labor. Instead, use [cherry-picking](https://gitbetter.substack.com/p/how-to-use-git-cherry-pick-effectively) if you need another branches' code before it has been merged into main.
* Commit often, after finishing a small part of a feature, even if the code does not work fully yet. Since you commit to your own feature branch, no one else is affected, and you can keep a history of your changes in case something goes wrong.
* Use descriptive commit messages (not just "fix bugs"). Follow [these guidelines](https://gist.github.com/robertpainsi/b632364184e70900af4ab688decf6f53). Use the imperative ("add" instead of "added") for verbs.
* No data should ever be committed to Git. The repository is for code only. Store any local data files in the `data/` folder which is ignored by commits.



## Data Structure
Since data is not committed but has to be downloaded by every developer, there should be a common structure of the `data/` folder:
* `measurements/`
  * `dataset_name/` (e.g., `ref_2023`)
    * `sheet.xlsx`: Post-flight data sheet
    * `measurements.mat`: Measurements from Flight Test Instrumentation System
