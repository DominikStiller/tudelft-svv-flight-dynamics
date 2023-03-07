from pathlib import Path

paths = []
paths.extend(Path("fd").glob("**/*.*"))
paths.extend(Path("tests").glob("**/*.*"))
paths.append(Path("bin/generate_code_for_appendix.py"))  # so meta

check_is_file = lambda f: f.is_file()
check_has_proper_extension = lambda f: f.suffix in [".py", ".yml"]
check_is_not_init = lambda f: f.name != "__init__.py"

paths = filter(
    lambda f: check_is_file(f) and check_has_proper_extension(f) and check_is_not_init(f), paths
)

with Path("data/appendix_code_generated.tex").open("w") as f:
    for code_file in paths:
        f.write("\\paragraph{" + str(code_file).replace("_", "\\_") + "}\n")
        f.write("\\begin{pythoncode}\n")
        f.write(code_file.read_text())
        f.write("\\end{pythoncode}\n")
        f.write("\n")
