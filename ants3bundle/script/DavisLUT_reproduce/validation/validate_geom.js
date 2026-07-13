// Literal geometry validation of the DavisLUT rule, run headless:  ants3 -j validate_geom.js
// (run from this validation/ directory).
//
// Setup: a collimated photon source inside a LYSO block fires at a fixed incidence angle
// (10, 45, 80 deg) at the LYSO->air interface (top face of the block), where the DavisLUT
// rule acts. A thin photon monitor just below the interface catches the reflected photons
// and records their polar angle theta_out (angle from the surface normal). The recorded
// theta_out distribution is written to data/geom_<angle>.txt for comparison with the LUT
// (see plot_geom.py). This exercises the full tracer -> rule -> sensor path.

var FWD = "../luts/pooled_28um_fwd.lut";   // LYSO -> air
var REV = "../luts/pooled_28um_rev.lut";   // air  -> LYSO (reverse direction)
var OUT = "geom_out";
var angles = [10, 45, 80];

core.createDir(OUT);
core.createDir(OUT + "/sim");
config.load("../config/config-28Yuan-airGaps-04-LUT.json");

// collimated single-point source just below the interface; 300k photons; save monitors
config.replace("PhotonSim.PhotonBombs.GenerationMode", "single");
config.replace("PhotonSim.PhotonBombs.Single.Position", [0, 0, -0.05]);
config.replace("PhotonSim.PhotonBombs.PhotonsPerBomb.Mode", "constant");
config.replace("PhotonSim.PhotonBombs.PhotonsPerBomb.FixedNumber", 300000);
config.replace("PhotonSim.PhGenOverrides.Direction.DirectionMode", "Fixed");
config.replace("PhotonSim.Run.OutputDirectory", OUT + "/sim");
config.replace("PhotonSim.Run.SaveMonitors", true);

for (var a = 0; a < angles.length; a++)
{
    var th = angles[a], rad = th * Math.PI / 180.0;
    config.replace("PhotonSim.PhGenOverrides.Direction.DirectionVector",
                   [Math.sin(rad), 0, Math.cos(rad)]);   // incidence plane = x-z
    config.updateConfig();

    // geometry: air world box, LYSO block (top face at z=0), monitor 0.1 mm below interface
    geo.clearWorld();
    geo.box("AirBox", [200, 200, 80], 0, "World", [0, 0, 0], [0, 0, 0]);      // mat 0 = air
    geo.box("LYSO",   [100, 100, 20], 1, "AirBox", [0, 0, -10], [0, 0, 0]);   // mat 1 = LYSO, z in [-20,0]
    geo.monitor("Mon", 0, 48, 48, "LYSO", [0, 0, 9.9], [0, 0, 0], true, true, true); // global z = -0.1
    geo.configurePhotonMonitor("Mon", [10, 10], [], [90, 0, 90], []);          // angle: 90 bins, 0..90 deg
    geo.updateGeometry();

    rules.setLutMaterialRule("LYSO", "air", FWD);
    rules.setLutMaterialRule("air", "LYSO", REV);

    lsim.simulate();
    lsim.loadMonitorData(OUT + "/sim/PhotonMonitors.txt");

    var ang = lsim.getMonitorAngle(0), tot = 0, s = "";
    for (var i = 0; i < ang.length; i++) { s += ang[i][0] + " " + ang[i][1] + "\n"; tot += ang[i][1]; }
    s = "# angle=" + th + " total_hits=" + tot + "\n" + s;
    core.saveText(s, "data/geom_" + th + ".txt");
}
