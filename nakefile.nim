import nake, std/strformat

task "docs", "Generate documentation":
  # https://nim-lang.github.io/Nim/docgen.html
  let
    src = "manu.nim"
    dir = "docs/"
    doc = dir / src.changeFileExt(".html")
    url = "https://github.com/planetis-m/manu"
  if doc.needsRefresh(src):
    echo "Generating the docs..."
    direShell(nimExe,
        &"doc --project --verbosity:0 --git.url:{url} --git.devel:master --git.commit:master --out:{dir} {src}")
    withDir(dir):
      moveFile("theindex.html", "index.html")
  else:
    echo "Skipped generating the docs."

task "test", "Run the tests":
  withDir("tests/"):
    for f in walkFiles("t*.nim"):
      direShell(nimExe, &"c -r --verbosity:0 --path:../ {f}")
