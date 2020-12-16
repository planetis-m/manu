# mode = ScriptMode.Verbose

packageName   = "manu"
version       = "2.2.0"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

requires "nim >= 1.4.0"

switch "forceBuild"

proc configForTests() =
   switch "hints", "off"
   switch "linedir", "on"
   switch "stacktrace", "on"
   switch "linetrace", "on"
   switch "debuginfo"
   switch "path", "."
   switch "gc", "arc"
   switch "run"

task test, "run tests":
   configForTests()
   setCommand "c", "tests/testMatrix.nim"

task doc, "generate documentation":
   switch "project"
   #switch "index"
   switch "out", "docs/"
   switch "git.url", "https://github.com/planetis-m/manu"
   #switch "docSeeSrcUrl", "https://github.com/planetis-m/manu/tree/master"
   setCommand "doc", "manu.nim"

after doc:
   switch "out", "docs/index.html"
   setCommand "buildIndex", "docs/"
