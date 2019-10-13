# mode = ScriptMode.Verbose

packageName   = "manu"
version       = "1.0"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "htmldocs", "examples", "experiments"]

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

task test, "run tests":
  configForTests()
  setCommand "c", "tests/testMatrix.nim"

task docs, "generate documentation":
  # exec("mkdir -p htmldocs/manu")
  switch "project"
  switch "docSeeSrcUrl", "https://github.com/b3liever/manu/blob/master"
  setCommand "doc", "manu.nim"
