# Package

packageName   = "manu"
version       = "2.3.1"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

# Deps

requires "nim >= 1.5.0"

import os

const
  ProjectUrl = "https://github.com/planetis-m/manu"
  PkgDir = thisDir().quoteShell
  DocsDir = PkgDir / "docs"

task docs, "Generate documentation":
  # https://nim-lang.github.io/Nim/docgen.html
  withDir(PkgDir):
    let src = "manu.nim"
    # Generate the docs for {src}
    exec("nim doc --project --verbosity:0 --git.url:" & ProjectUrl &
        " --git.devel:master --git.commit:master --out:" & DocsDir & " " & src)
    mvFile(DocsDir / "theindex.html", DocsDir / "index.html")

# task test, "Run the tests":
#   withDir(PkgDir):
#     for f in listFiles("tests"):
#       if f.endsWith(".nim"):
#         echo "Running ", f, "..."
#         exec("nim c -r --hints:off -w:off " & quoteShell(f))
