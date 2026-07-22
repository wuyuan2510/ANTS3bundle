// Direct runtime validation of DavisLUT through ALutInterfaceRule::calculate().
// Run headless from this directory:
//     ../../../bin/ants3 -j validate_rule.js
//
// This validates both LUT directions, reflected and transmitted 2D distributions,
// probability conservation, normal-incidence handling, rotated-normal covariance,
// output-vector normalization/signs, and a synthetic non-zero absorption branch.

var cases = [
    { tag: "fwd", lut: "../luts/pooled_28um_fwd.lut" },
    { tag: "rev", lut: "../luts/pooled_28um_rev.lut" }
];

var params = {
    angles: [0, 10, 32, 33.3, 34, 45, 80, 89],
    photonsPerAngle: 500000,
    seed: 20260714,
    // The full 2D distributions are always checked. Set this true only when the flattened
    // measured/expected maps are also wanted in the JSON output (about 4 MB for both files).
    include2D: false
};

core.createDir("data/runtime");
for (var oldTag of ["fwd", "rev"])
{
    var oldResult = "data/runtime/" + oldTag + ".json";
    if (core.isFileExist(oldResult)) core.deleteFile(oldResult);
}
var allPassed = true;

for (var ic = 0; ic < cases.length; ic++)
{
    var c = cases[ic];
    core.print("Validating", c.tag, c.lut);
    var result = rules.validateSurfaceLut(c.lut, params);
    core.saveObject(result, "data/runtime/" + c.tag + ".json");

    core.print("  pass:", result.passed,
               "max |v|-1:", result.maximumNormError,
               "direction errors:", result.directionErrors,
               "status errors:", result.statusErrors);

    for (var ia = 0; ia < result.angles.length; ia++)
    {
        var a = result.angles[ia];
        core.print("  angle", a.angle,
                   "pass", a.passed,
                   "R", a.measuredR, a.expectedR,
                   "T", a.measuredT, a.expectedT);
    }
    allPassed = allPassed && result.passed;
}

if (!allPassed) core.abort("DavisLUT direct runtime validation FAILED; inspect data/runtime/*.json");
core.print("DavisLUT direct runtime validation PASSED");
