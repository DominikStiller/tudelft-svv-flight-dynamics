# Simulation, Verification and Validation: Flight dynamics assignment
> Flight dynamics of an aircraft

This is the repository for the TU Delft [AE3212-II](https://studiegids.tudelft.nl/a101_displayCourse.do?course_id=62266) project (second assignment) of group B24. The goal is to build a flight dynamics simulation, then verify and validate it with flight test data.

Links:
* See the [developer's guidelines](docs/CONTRIBUTING.md) if you want to contribute
* See the [assignment](docs/assignment.pdf) for a description of the tasks


## Setup
To get started with development:
1. Ensure that Python 3.11 installed.
2. Clone the GitHub repository by clicking on the green "Code" button above and follow the instructions.
3. Open the cloned folder in PyCharm (other IDEs can be used, adjust the following instructions accordingly).
4. Add a new interpreter in a [virtualenv environment](https://docs.python.org/3/tutorial/venv.html). This ensures isolation so that the packages for this project do not conflict with your preinstalled ones.
5. Install all required packages by opening `requirements.txt` and clicking "Install requirements". This ensures that everyone uses the same package versions, preventing bugs that might be hard to find.
6. Read the code and Git guidelines in this document so the code stays consistent and clear.


## Usage
To run the code from the console, make sure you are in the project directory and activated the virtual environment (`<venv>\Scripts\activate.bat` on Windows, `source <venv>/bin/activate` on Linux).

The main script can then be executed using `python -m fd`. Alternatively, the scripts in `bin` can be used which automatically activate the virtual environment and offer the same parameters.
