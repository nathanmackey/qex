import qex
import mcmc/mcmc
import gauge/gaugeAction
import strutils

qexInit()

let 
  # Get command line information
  cmdInfo = readCMD()
  ntraj = cmdInfo["ntraj"].getInt()
  saveInterval = cmdInfo["saveInterval"].getInt()
  ioPath = cmdInfo["ioPath"].getStr()
  fn = cmdInfo["fn"].getStr()
  jsonPath = cmdInfo["jsonPath"].getStr()
  flowString = cmdInfo["flow"].getStr()
  flow = parseBool(flowString)
  flowInterval = case flow
    of true: cmdInfo["flowInterval"].getInt()
    of false: 0
  
  # Get information from JSON file
  jsonInfo = readJSON(jsonPath)
  tau = jsonInfo["hmc"]["trajectory-length"].getFloat()
  nspv = jsonInfo["staggered-pauli-villars"]["species"].getInt()
  nsf = jsonInfo["staggered-fermions"]["species"].getInt()
  nrsf = jsonInfo["rooted-staggered-fermions"]["species"].getInt()
var 
  cfg = cmdInfo["cfg"].getInt()
  flowInfo = case flow 
    of true: jsonInfo["flow"]
    of false: parseJson("{}")

# Construct lattice field theory
var 
  hmc = newLatticeFieldTheory(jsonInfo["hmc"]):
    fieldTheory.addGaugeMatterAction(jsonInfo["action"]):
      action.addGaugeField(jsonInfo["gauge"])
# Read in gauge field
if cfg != 0: hmc.read(ioPath & fn & "_" & $(cfg), onlyGauge = false)

# Run HMC
let g = hmc.u[]
let lo = g[0].l
unit(g)
let coords = @[0,0,0,0]
let (rk, idx) = lo.rankIndex(coords)

var f = newOneOf g

var gc = GaugeActionCoeffs(plaq: 1.0)

var s = gc.gaugeActionBP(g)
echo "s"

let eps = 1e-4

var r = lo.newRNGField(MRG32k3a, 987654321)
echo "r"

var rDir = lo.newGauge

rDir.randomTAH r 

echo "rDir"

var x: type(g[0]{1})

x = rDir[0]{idx}

echo "x"

var pert: type(x*1)

pert :=  x
echo "Pert: ", pert

var u: type(g[0]{1})
u = g[0]{idx} 
g[0]{idx} := exp(eps*x) * u
var Up: type(g[0]{1} * 1)
Up := g[0]{idx}
echo "perturbed link: ", Up
var sp = gc.gaugeActionBP(g)
echo "Starting Action: ", s
echo "Perturbed Action: ", sp
gc.gaugeForceBP(g,f)
echo "Force: ", (f[0]{idx}*1)
echo "delta S (force): ", redot(x,f[0]{idx})
echo "delta S (action): ", s-sp
qexFinalize()


