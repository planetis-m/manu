# mode = ScriptMode.Verbose

packageName   = "manu"
version       = "2.0"
author        = "Antonis Geralis"
description   = "Nim Matrix library"
license       = "MIT"
skipDirs = @["tests", "docs", "examples", "experiments"]

requires "nim >= 1.2.0"

switch "forceBuild"

task doc, "generate documentation":
   switch "project"
   switch "index"
   switch "out", "docs/"
   switch "git.url", "https://github.com/b3liever/manu"
   #switch "docSeeSrcUrl", "https://github.com/b3liever/manu/tree/master"
   setCommand "doc", "manu.nim"

after doc:
   switch "out", "docs/index.html"
   setCommand "buildIndex", "docs/"
