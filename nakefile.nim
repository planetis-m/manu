import nake, std/strformat

task "docs", "Generate documentation":
  # https://nim-lang.github.io/Nim/docgen.html
  let
    manu = "manu"
    src = [
      manu.addFileExt(".nim"),
      manu / "matrix.nim", manu / "cholesky.nim",
      manu / "qr.nim", manu / "lu.nim",
      manu / "eigen.nim", manu / "svd.nim"
    ]
    dir = "docs/"
    doc = dir / manu.addFileExt(".html")
    url = "https://github.com/planetis-m/manu"
  if doc.needsRefresh(src):
    echo "Generating the docs..."
    direShell(nimExe,
        &"doc --project --verbosity:0 --git.url:{url} --git.devel:master --git.commit:master --out:{dir} {src[0]}")
    withDir(dir):
      moveFile("theindex.html", "index.html")
  else:
    echo "Skipped generating the docs."

task "test", "Run the tests":
  withDir("tests/"):
    for f in walkFiles("t*.nim"):
      direShell(nimExe, &"c -r --hints:off -w:off --path:../ {f}")
