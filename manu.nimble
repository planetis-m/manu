# mode = ScriptMode.Verbose

packageName   = "manu"
version       = "1.2"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

requires "nim >= 1.0.0"

switch "forceBuild"

proc configForTests() =
   switch "hints", "off"
   switch "linedir", "on"
   switch "stacktrace", "on"
   switch "linetrace", "on"
   switch "debuginfo"
   switch "path", "."
   switch "run"

task tests, "run tests":
   configForTests()
   #setCommand "c", "tests/testLu.nim"
   #setCommand "c", "tests/testCholesky.nim"
   setCommand "c", "tests/testMatrix.nim"

task docs, "generate documentation":
   switch "project"
   switch "index"
   switch "docSeeSrcUrl", "https://github.com/b3liever/manu/master"
   setCommand "doc", "manu.nim"

after docs:
   exec("mv htmldocs/ docs/")
   switch "out", "htmldocs/index.html"
   setCommand "buildIndex", "docs/"
   exec("mv docs/theindex.html docs/index.html")
