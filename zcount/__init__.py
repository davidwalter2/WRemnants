import ROOT
import narf
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

narf.clingutils.Declare('#include "definitions.h"')
narf.clingutils.Declare('#include "utilities.h"')
narf.clingutils.Declare('#include "pairing.h"')
