# mode = ScriptMode.Verbose

packageName   = "manu"
version       = "2.1.1"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

requires "nim >= 1.4.0"

switch "forceBuild"

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
