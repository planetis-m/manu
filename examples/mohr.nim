# From: https://csmbrannon.net/2011/04/25/stress-state-analysis-python-script/
import std/[algorithm, os, math, strutils, strformat], manu

# ------------------------------------------------------------
# Helper functions for uniform printing throughout the script.
# ------------------------------------------------------------

proc headerprint(s: string) =
  # Prints a centered string to divide output sections.
  const width = 64
  const ch = '='
  let numspaces = width - len(s)
  let before = int(ceil(numspaces / 2))
  let after = int(floor(numspaces / 2))
  echo("\n", ch.repeat(before), s, ch.repeat(after), "\n")

proc valprint(s: string, value: float) =
  # Ensure uniform formatting of scalar value outputs.
  echo(s.align(30), ": ", value.formatEng)

proc matprint(s: string, value: Matrix) =
  # Ensure uniform formatting of matrix value outputs.
  echo(s, ":\n", value)

proc usage() =
  # When the user needs help, print the script usage.
  headerprint(" Analyze Stress State ")
  let appname = getAppFilename().extractFilename()
  let s = &""" For a given stress state this script computes many
 useful quantities that help to analyze the stress state.

 Currently, the following values are output:
   Isotropic Matrix
   Deviatoric Matrix
   Principal Stresses
   Maximum Shear
   Mean Stress
   Equivalent Stress
   Invariant I1
   Invariant J2
   Invariant J3
   Lode Coordinates
   Triaxiality

 Command line syntax option 1:

     > ./{appname} sig11 sig22 sig33

 Command line syntax option 2:

     > ./{appname} sig11 sig22 sig33 sig12 sig13 sig23"""
  quit(s)

# ------------------------------
# The main section of the script
# ------------------------------

proc main() =
  # If the script is run with an incorrect number of arguments
  # or if the user is asking for help, print the usage information.
  let params = commandLineParams()
  if "--help" in params or "-h" in params or
      params.len != 3 and params.len != 6:
    usage()

  # load stress components from the command line in a temporary
  # container
  var dum = newSeq[float](6)
  for idx in 0 ..< params.len:
    try:
      dum[idx] = parseFloat(params[idx])
    except ValueError:
      quit(&"Argument '{params[idx]}' is not a valid float")

  # load the stresses into our matrix and compute the
  # deviatoric and isotropic stress matricies
  let
    sigma = matrix(@[
      @[dum[0], dum[3], dum[4]],
      @[dum[3], dum[1], dum[5]],
      @[dum[4], dum[5], dum[2]]
    ])
    sigmaIso = 1.0/3.0*trace(sigma)*eye64(sigma.m)
    sigmaDev = sigma - sigmaIso

  # compute principal stresses
  var eigVals = eig(sigma).getRealEigenvalues
  sort(eigVals)
  reverse(eigVals)

  # compute max shear stress
  let maxShear = (eigVals[0]-eigVals[2])/2.0

  # compute the stress invariants
  let
    I1 = trace(sigma)
    J2 = 1.0/2.0*trace(sigmaDev*sigmaDev)
    J3 = 1.0/3.0*trace(sigmaDev*sigmaDev*sigmaDev)

  # compute other common stress measures
  let
    meanStress = 1.0/3.0*I1
    eqvStress = sqrt(3.0*J2)

  # compute lode coordinates
  let
    lodeR = sqrt(2.0*J2)
    lodeZ = I1/sqrt(3.0)

    temp = 3.0*sqrt(6.0)*det(sigmaDev/lodeR)
    lodeTheta = 1.0/3.0*arcsin(temp)

  # compute the stress triaxiality
  let triaxiality = meanStress/eqvStress

  # Print out what we've found
  headerprint(" Stress State Analysis ")
  matprint("Input Stress", sigma)
  headerprint(" Component Matricies ")
  matprint("Isotropic Stress", sigmaIso)
  matprint("Deviatoric Stress", sigmaDev)
  headerprint(" Scalar Values ")
  valprint("P1", eigVals[0])
  valprint("P2", eigVals[1])
  valprint("P3", eigVals[2])
  valprint("Max Shear", maxShear)
  valprint("Mean Stress", meanStress)
  valprint("Equivalent Stress", eqvStress)
  valprint("I1", I1)
  valprint("J2", J2)
  valprint("J3", J3)
  valprint("Lode z", lodeZ)
  valprint("Lode r", lodeR)
  valprint("Lode theta (rad)", lodeTheta)
  valprint("Lode theta (deg)", radToDeg(lodeTheta))
  valprint("Triaxiality", triaxiality)
  headerprint(" End Output ")


main()

# End of script
