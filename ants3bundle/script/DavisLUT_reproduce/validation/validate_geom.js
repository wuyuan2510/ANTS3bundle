// Full geometry validation of DavisLUT, run headless from this directory:
//     ../../../bin/ants3 -j validate_geom.js
//
// Deprecated: repeated dispatcher simulations in one GUI JavaScript event loop can reuse a
// later source configuration for several nominally different cases. Use validate_geom.py,
// which runs one static configuration in a fresh lsim process per direction/angle.
core.abort("validate_geom.js is deprecated because its dispatcher loop is state-racy; run python3 validate_geom.py");

// The setup uses ideal, non-absorbing transport and two very close photon monitors so that
// absolute R/T probabilities as well as conditional theta_out distributions can be compared
// with the LUT. Both forward (LYSO->air) and reverse (air->LYSO) LUTs are exercised.

var FWD = "../luts/pooled_28um_fwd.lut";
var REV = "../luts/pooled_28um_rev.lut";
var OUT = "geom_out";
var DATA = "data/geometry";
var PHOTONS = 300000;
var ANGLES = [0, 10, 32, 33.3, 34, 45, 80, 89];
var cases = [
    { tag: "fwd", reverse: false },
    { tag: "rev", reverse: true  }
];

core.createDir(OUT);
core.createDir(OUT + "/sim");
core.createDir(DATA);
for (var staleFile of [DATA + "/summary.json", DATA + "/metrics.json"])
    if (core.isFileExist(staleFile)) core.deleteFile(staleFile);
config.load("../config/config-28Yuan-airGaps-04-LUT.json");

// Remove unrelated bulk-transport effects: this is an interface-rule conformance experiment.
config.replace("Materials[0].AbsCoeff", 0);
config.replace("Materials[0].BulkAbsorptionWave", []);
config.replace("Materials[0].RayleighMFP", 0);
config.replace("Materials[1].AbsCoeff", 0);
config.replace("Materials[1].BulkAbsorptionWave", []);
config.replace("Materials[1].RayleighMFP", 0);
config.replace("PhotonSim.PhGenOverrides.Wave.Enabled", true);
config.replace("PhotonSim.PhGenOverrides.Wave.FixedWavelength", 420);
config.replace("PhotonSim.PhotonBombs.GenerationMode", "single");
config.replace("PhotonSim.PhotonBombs.PhotonsPerBomb.Mode", "constant");
config.replace("PhotonSim.PhotonBombs.PhotonsPerBomb.FixedNumber", PHOTONS);
config.replace("PhotonSim.PhGenOverrides.Direction.DirectionMode", "Fixed");
config.replace("PhotonSim.Run.OutputDirectory", OUT + "/sim");
config.replace("PhotonSim.Run.SaveMonitors", true);
config.replace("PhotonSim.Run.SaveTracks", false);

var summaries = [];
for (var ic = 0; ic < cases.length; ic++)
{
    var c = cases[ic];
    for (var ia = 0; ia < ANGLES.length; ia++)
    {
        var angle = ANGLES[ia];
        var rad = angle*Math.PI/180.0;
        var depth = 0.002; // mm: source is between the interface and its reflection monitor
        var zSource = c.reverse ? depth : -depth;
        var dz = c.reverse ? -Math.cos(rad) : Math.cos(rad);
        var sourceX = -depth*Math.tan(rad); // pencil beam hits the interface at x=0

        config.replace("PhotonSim.PhotonBombs.Single.Position", [sourceX, 0, zSource]);
        config.replace("PhotonSim.PhGenOverrides.Direction.DirectionVector", [Math.sin(rad), 0, dz]);
        config.updateConfig();

        // Interface at z=0. Monitors are only 4 um away, avoiding side loss even at 89.9 deg.
        geo.clearWorld();
        geo.box("AirBox", [120, 120, 40], 0, "World", [0, 0, 0], [0, 0, 0]);
        geo.box("LYSO",   [100, 100, 20], 1, "AirBox", [0, 0, -10], [0, 0, 0]);
        geo.monitor("LowerMonitor", 0, 40, 40, "LYSO",   [0, 0, 9.996], [0, 0, 0], true, true, true);
        geo.monitor("UpperMonitor", 0, 40, 40, "AirBox", [0, 0, 0.004], [0, 0, 0], true, true, true);
        geo.configurePhotonMonitor("LowerMonitor", [80, 80], [], [90, 0, 90], []);
        geo.configurePhotonMonitor("UpperMonitor", [80, 80], [], [90, 0, 90], []);
        geo.updateGeometry();

        rules.setLutMaterialRule("LYSO", "air", FWD);
        rules.setLutMaterialRule("air", "LYSO", REV);

        lsim.setSeed(20260714 + 1000*ic + ia);
        lsim.simulate();
        lsim.loadMonitorData(OUT + "/sim/PhotonMonitors.txt");

        // Monitors are registered in creation order: lower=0, upper=1.
        var reflectedIndex = c.reverse ? 1 : 0;
        var transmittedIndex = c.reverse ? 0 : 1;
        var reflected = lsim.getMonitorAngle(reflectedIndex);
        var transmitted = lsim.getMonitorAngle(transmittedIndex);
        var hits = lsim.getMonitorHitsAll();
        var reflectedHits = hits[reflectedIndex];
        var transmittedHits = hits[transmittedIndex];
        var tag = angleTag(angle);

        saveDistribution(reflected,
                         DATA + "/geom_" + c.tag + "_reflected_" + tag + ".txt",
                         c.tag, angle, "reflected", reflectedHits);
        saveDistribution(transmitted,
                         DATA + "/geom_" + c.tag + "_transmitted_" + tag + ".txt",
                         c.tag, angle, "transmitted", transmittedHits);

        summaries.push({direction:c.tag, angle:angle,
                        reflected:reflectedHits, transmitted:transmittedHits,
                        launched:PHOTONS, missing:PHOTONS-reflectedHits-transmittedHits});
        core.print(c.tag, "angle", angle,
                   "R", reflectedHits/PHOTONS,
                   "T", transmittedHits/PHOTONS,
                   "missing", PHOTONS-reflectedHits-transmittedHits);
    }
}

core.saveObject({photonsPerRun:PHOTONS, wavelength:420, runs:summaries}, DATA + "/summary.json");
core.print("Geometry validation completed; results are in", DATA);

function angleTag(angle)
{
    return String(angle).replace(".", "p");
}

function saveDistribution(data, fileName, direction, angle, outcome, total)
{
    var text = "# direction=" + direction + " angle=" + angle +
               " outcome=" + outcome + " total_hits=" + total + "\n";
    for (var i = 0; i < data.length; i++) text += data[i][0] + " " + data[i][1] + "\n";
    core.saveText(text, fileName);
}
