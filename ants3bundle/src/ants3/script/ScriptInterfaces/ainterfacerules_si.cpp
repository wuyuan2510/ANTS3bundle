#include "ainterfacerules_si.h"
#include "ainterfacerulehub.h"
#include "alutinterfacerule.h"
#include "alutsurfacedata.h"
#include "alutsurfacegenerator.h"
#include "amaterialhub.h"
#include "ageometryhub.h"
#include "ageoobject.h"
#include "aphoton.h"
#include "arandomhub.h"
#include "ascripthub.h"
#include "ajsontools.h"

#include <QFileInfo>
#include <QJsonObject>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

AInterfaceRules_SI::AInterfaceRules_SI() :
    RuleHub(AInterfaceRuleHub::getInstance())
{
    Description = "Optical interface rules: generation of surface LUTs (\"DavisLUT\" rule) and rule assignment";

    Help["generateSurfaceLut"] = "Generates a surface LUT for the 'DavisLUT' interface rule by ray tracing over a 3D surface\n"
                                 "topography Z(x,y) (e.g. measured with AFM), see Roncali & Cherry, Phys.Med.Biol. 58 (2013) 2185.\n"
                                 "Arguments: heightmap file name, output LUT file name, object with parameters.\n"
                                 "Parameter keys (* = mandatory):\n"
                                 "  n1*: refractive index of the medium the photons arrive from (e.g. the crystal)\n"
                                 "  n2*: refractive index of the medium behind the surface\n"
                                 "  format: 'matrix' (default) = text file with a matrix of heights (rows <-> y, columns <-> x),\n"
                                 "          requires pixelSizeX and pixelSizeY; 'xyz' = text file with 3 columns (x y z) on a regular grid\n"
                                 "  pixelSizeX, pixelSizeY: pixel size for the 'matrix' format (same length units as the heights)\n"
                                 "  wavelength: [nm] recorded in the LUT metadata (default 420)\n"
                                 "  thetaBins (default 40), thetaOutBins (45), phiOutBins (36): LUT binning\n"
                                 "  photonsPerBin: photons per incidence angle bin (default 80000)\n"
                                 "  phiSteps: 0 (default) = random incidence azimuth, >0 = discrete azimuth steps\n"
                                 "  maxBounces (default 100), seed (0 = do not reseed), comment\n"
                                 "  reverseGeometry: true if the photons arrive from the medium ABOVE the heightmap surface (default false)\n"
                                 "  alsoReverse: true = also generate the LUT for the reverse direction (swapped n1/n2) and save it\n"
                                 "               to a file with '_reverse' appended to the name\n"
                                 "Opposite heightmap edges must match (C0-periodic); mirror-tile generic AFM scans first.\n"
                                 "Returns an object with generation statistics, including wraps and\n"
                                 "direction/medium-inconsistent escape counters for periodic-seam diagnostics.\n"
                                 "Note: the LUT is direction-specific, assign it only for the n1->n2 material pair!";
    Help["getLutInfo"]         = "Returns the metadata and binning info of the given surface LUT file as an object";
    Help["validateSurfaceLut"] = "Validates a surface LUT through the actual ALutInterfaceRule::calculate() runtime path.\n"
                                 "Arguments: LUT file name and an optional parameters object. Parameters:\n"
                                 "  angles: incidence angles in degrees (default [0,10,32,33.3,34,45,80,89])\n"
                                 "  photonsPerAngle: photons sampled at each angle (default 500000)\n"
                                 "  seed: deterministic random seed (default 12345)\n"
                                 "  include2D: include flattened measured/expected theta-phi arrays in the result (default false)\n"
                                 "Returns pass/fail thresholds and measured/expected R/T/A probabilities and angular distributions.\n"
                                 "It also checks rotated-normal covariance, unit output vectors, direction signs and a synthetic absorption branch.";
    Help["setLutMaterialRule"] = "Creates a 'DavisLUT' interface rule with the LUT loaded from the given file\n"
                                 "and sets it for the material pair (from, to). Materials are given by name.\n"
                                 "The LUT data are embedded in the config, the file is not needed afterwards.";
    Help["setLutVolumeRule"]   = "Creates a 'DavisLUT' interface rule with the LUT loaded from the given file\n"
                                 "and sets it for the volume pair (from, to). Volumes are given by name.";
    Help["clearMaterialRule"]  = "Removes the interface rule (of any type) defined for the material pair (from, to)";
    Help["clearVolumeRule"]    = "Removes the interface rule (of any type) defined for the volume pair (from, to)";
}

QVariantMap AInterfaceRules_SI::generateSurfaceLut(QString heightmapFile, QString outputLutFile, QVariantMap params)
{
    QVariantMap result;

    if (!params.contains("n1") || !params.contains("n2"))
    {
        abort("generateSurfaceLut: parameters 'n1' and 'n2' are mandatory");
        return result;
    }

    ALutSurfaceGenerator generator;
    generator.n1 = params["n1"].toDouble();
    generator.n2 = params["n2"].toDouble();
    if (params.contains("wavelength"))    generator.wavelength         = params["wavelength"].toDouble();
    if (params.contains("thetaBins"))     generator.thetaIncBins       = params["thetaBins"].toInt();
    if (params.contains("thetaOutBins"))  generator.thetaOutBins       = params["thetaOutBins"].toInt();
    if (params.contains("phiOutBins"))    generator.phiOutBins         = params["phiOutBins"].toInt();
    if (params.contains("photonsPerBin")) generator.photonsPerThetaBin = params["photonsPerBin"].toInt();
    if (params.contains("phiSteps"))      generator.phiSteps           = params["phiSteps"].toInt();
    if (params.contains("maxBounces"))    generator.maxBounces         = params["maxBounces"].toInt();
    if (params.contains("seed"))          generator.seed               = params["seed"].toInt();
    if (params.contains("comment"))       generator.comment            = params["comment"].toString();
    if (params.contains("reverseGeometry")) generator.reverseGeometry  = params["reverseGeometry"].toBool();

    const QString format = (params.contains("format") ? params["format"].toString() : "matrix");
    QString err;
    if (format == "matrix")
    {
        const double pixelSizeX = params["pixelSizeX"].toDouble();
        const double pixelSizeY = params["pixelSizeY"].toDouble();
        if (pixelSizeX <= 0 || pixelSizeY <= 0)
        {
            abort("generateSurfaceLut: 'matrix' format requires positive 'pixelSizeX' and 'pixelSizeY' parameters");
            return result;
        }
        err = generator.loadHeightmapMatrix(heightmapFile, pixelSizeX, pixelSizeY);
    }
    else if (format == "xyz") err = generator.loadHeightmapXYZ(heightmapFile);
    else
    {
        abort("generateSurfaceLut: 'format' should be 'matrix' or 'xyz'");
        return result;
    }
    if (!err.isEmpty())
    {
        abort("generateSurfaceLut: " + err);
        return result;
    }

    generator.progressCallback = [this](int percent) -> bool
    {
        AScriptHub::getInstance().reportProgress(percent, Lang);
        return !AScriptHub::isAborted(Lang);
    };

    auto runAndSave = [&generator, this](const QString & fileName) -> QVariantMap
    {
        QVariantMap summary;
        ALutSurfaceData data;
        QString err = generator.generate(data);
        if (!err.isEmpty())
        {
            abort("generateSurfaceLut: " + err);
            return summary;
        }
        QJsonObject json;
        data.writeToJson(json);
        if (!jstools::saveJsonToFile(json, fileName))
        {
            abort("generateSurfaceLut: cannot save LUT to file " + fileName);
            return summary;
        }
        summary["lutFile"]     = fileName;
        summary["n1"]          = generator.n1;
        summary["n2"]          = generator.n2;
        summary["meanBounces"] = generator.meanBounces();
        summary["anomalies"]   = (qlonglong)generator.anomalies();
        summary["wraps"]       = (qlonglong)generator.wraps();
        summary["upEscapeReclassified"]   = (qlonglong)generator.upEscapeReclassified();
        summary["downEscapeReclassified"] = (qlonglong)generator.downEscapeReclassified();
        summary["upEscapeAfterSeam"]      = (qlonglong)generator.upEscapeAfterSeam();
        summary["downEscapeAfterSeam"]    = (qlonglong)generator.downEscapeAfterSeam();
        summary["degenerateDiscarded"]    = (qlonglong)generator.degenerateDiscarded();
        summary["gridSizeX"]   = data.GridSizeX;
        summary["gridSizeY"]   = data.GridSizeY;
        return summary;
    };

    result = runAndSave(outputLutFile);
    if (result.isEmpty()) return result;   // aborted

    if (params.contains("alsoReverse") && params["alsoReverse"].toBool())
    {
        std::swap(generator.n1, generator.n2);
        generator.reverseGeometry = !generator.reverseGeometry;

        const QFileInfo fileInfo(outputLutFile);
        const QString suffix = fileInfo.suffix();
        QString reverseName = fileInfo.path() + "/" + fileInfo.completeBaseName() + "_reverse";
        if (!suffix.isEmpty()) reverseName += "." + suffix;

        const QVariantMap reverseSummary = runAndSave(reverseName);
        if (reverseSummary.isEmpty()) return result;
        result["reverse"] = reverseSummary;
    }

    return result;
}

QVariantMap AInterfaceRules_SI::getLutInfo(QString lutFile)
{
    QVariantMap result;

    QJsonObject json;
    if (!jstools::loadJsonFromFile(json, lutFile))
    {
        abort("Cannot open or parse LUT file " + lutFile);
        return result;
    }
    ALutSurfaceData data;
    const QString err = data.readFromJson(json);
    if (!err.isEmpty())
    {
        abort(err);
        return result;
    }

    result["n1"]              = data.n1;
    result["n2"]              = data.n2;
    result["wavelength"]      = data.Wavelength;
    result["sourceHeightmap"] = data.SourceHeightmap;
    result["gridSizeX"]       = data.GridSizeX;
    result["gridSizeY"]       = data.GridSizeY;
    result["pixelSizeX"]      = data.PixelSizeX;
    result["pixelSizeY"]      = data.PixelSizeY;
    result["thetaIncBins"]    = data.ThetaIncBins;
    result["thetaOutBins"]    = data.ThetaOutBins;
    result["phiOutBins"]      = data.PhiOutBins;
    result["photonsPerBin"]   = data.PhotonsPerThetaBin;
    result["meanBounces"]     = data.MeanBounces;
    result["anomalyCount"]    = data.AnomalyCount;
    result["generationDate"]  = data.GenerationDate;
    result["comment"]         = data.Comment;
    return result;
}

QVariantMap AInterfaceRules_SI::validateSurfaceLut(QString lutFile, QVariantMap params)
{
    QVariantMap result;

    ALutInterfaceRule rule(-1, -1);
    QString err = rule.loadLUT(lutFile);
    if (err.isEmpty()) err = rule.Data.buildRuntime();
    if (!err.isEmpty())
    {
        abort("validateSurfaceLut: " + err);
        return result;
    }

    QVariantList angleVariants;
    if (params.contains("angles")) angleVariants = params["angles"].toList();
    else angleVariants = QVariantList{0.0, 10.0, 32.0, 33.3, 34.0, 45.0, 80.0, 89.0};
    if (angleVariants.isEmpty())
    {
        abort("validateSurfaceLut: 'angles' should not be empty");
        return result;
    }

    const int photonsPerAngle = params.contains("photonsPerAngle") ? params["photonsPerAngle"].toInt() : 500000;
    const int seed = params.contains("seed") ? params["seed"].toInt() : 12345;
    const bool include2D = params.contains("include2D") && params["include2D"].toBool();
    if (photonsPerAngle < 1000)
    {
        abort("validateSurfaceLut: 'photonsPerAngle' should be at least 1000");
        return result;
    }

    const ALutSurfaceData & data = rule.Data;
    const int nTheta = data.ThetaOutBins;
    const int nPhi = data.PhiOutBins;
    const int nCells = nTheta * nPhi;

    auto vectorToVariants = [](const std::vector<double> & vec)
    {
        QVariantList list;
        list.reserve(vec.size());
        for (double v : vec) list.push_back(v);
        return list;
    };

    auto thetaMarginal = [nTheta, nPhi](const std::vector<double> & distribution)
    {
        std::vector<double> marginal(nTheta, 0.0);
        for (int it = 0; it < nTheta; it++)
            for (int ip = 0; ip < nPhi; ip++)
                marginal[it] += distribution[it*nPhi + ip];
        return marginal;
    };

    auto totalVariation = [](const std::vector<double> & a, const std::vector<double> & b)
    {
        double sum = 0;
        for (size_t i = 0; i < a.size(); i++) sum += std::abs(a[i] - b[i]);
        return 0.5 * sum;
    };

    auto normalizeCounts = [](const std::vector<qlonglong> & counts, qlonglong total)
    {
        std::vector<double> distribution(counts.size(), 0.0);
        if (total > 0)
            for (size_t i = 0; i < counts.size(); i++) distribution[i] = counts[i] / (double)total;
        return distribution;
    };

    auto stochasticNoiseTvd = [](const std::vector<double> & expected, qlonglong samples)
    {
        if (samples < 1) return 1.0;
        double meanTvd = 0;
        for (double p : expected)
            meanTvd += std::sqrt(std::max(0.0, 2.0*p*(1.0-p)/(M_PI*samples)));
        return 0.5 * meanTvd;
    };

    auto expectedDistribution = [&data, nCells](double angle, bool reflected,
                                                 std::vector<double> & distribution, double & probability)
    {
        const double binWidth = 90.0 / data.ThetaIncBins;
        const double x = angle/binWidth - 0.5;
        int k0 = (int)std::floor(x);
        double f = x - k0;
        if (k0 < 0) { k0 = 0; f = 0; }
        if (k0 >= data.ThetaIncBins-1) { k0 = data.ThetaIncBins-1; f = 0; }
        const int k1 = (f > 0 ? k0+1 : k0);
        const double w0 = 1.0-f;
        const double w1 = (k1 == k0 ? 0.0 : f);

        auto count = [&data, reflected](int k) -> double
        {
            return reflected ? data.ReflectedCounts[k] : data.TransmittedCounts[k];
        };
        auto hist = [&data, reflected](int k) -> const std::vector<int> &
        {
            return reflected ? data.ReflectedHist[k] : data.TransmittedHist[k];
        };

        probability = w0*count(k0)/data.Launched[k0] + w1*count(k1)/data.Launched[k1];
        distribution.assign(nCells, 0.0);
        if (probability <= 0) return;

        const std::vector<int> & h0 = hist(k0);
        const std::vector<int> & h1 = hist(k1);
        for (int i = 0; i < nCells; i++)
            distribution[i] = (w0*h0[i]/(double)data.Launched[k0] +
                               w1*h1[i]/(double)data.Launched[k1]) / probability;
    };

    auto expectedAbsorption = [&data](double angle)
    {
        const double binWidth = 90.0 / data.ThetaIncBins;
        const double x = angle/binWidth - 0.5;
        int k0 = (int)std::floor(x);
        double f = x-k0;
        if (k0 < 0) { k0 = 0; f = 0; }
        if (k0 >= data.ThetaIncBins-1) { k0 = data.ThetaIncBins-1; f = 0; }
        const int k1 = (f > 0 ? k0+1 : k0);
        const double w1 = (k1 == k0 ? 0.0 : f);
        return (1.0-f)*data.AbsorbedCounts[k0]/(double)data.Launched[k0] +
               w1*data.AbsorbedCounts[k1]/(double)data.Launched[k1];
    };

    auto probabilityTolerance = [photonsPerAngle](double probability)
    {
        const double sigma = std::sqrt(std::max(probability*(1.0-probability), 1.0/photonsPerAngle) / photonsPerAngle);
        return std::max(0.001, 6.0*sigma);
    };

    ARandomHub & randomHub = ARandomHub::getInstance();
    QVariantList angleResults;
    bool allPassed = true;
    double maximumNormError = 0;
    qlonglong totalDirectionErrors = 0;
    qlonglong totalStatusErrors = 0;

    for (int ia = 0; ia < angleVariants.size(); ia++)
    {
        const double angle = angleVariants[ia].toDouble();
        if (angle < 0 || angle >= 90.0)
        {
            abort("validateSurfaceLut: every incidence angle should be in [0, 90) degrees");
            return QVariantMap();
        }

        randomHub.setSeed(seed + 104729*ia);
        const double theta = angle*M_PI/180.0;
        const double normal[3] = {0, 0, 1};
        const double ex[3] = {1, 0, 0};
        const double ey[3] = {0, 1, 0};

        std::vector<qlonglong> measuredReflection(nCells, 0);
        std::vector<qlonglong> measuredTransmission(nCells, 0);
        qlonglong numReflection = 0;
        qlonglong numTransmission = 0;
        qlonglong numAbsorption = 0;
        qlonglong numErrors = 0;
        qlonglong directionErrors = 0;
        qlonglong statusErrors = 0;

        for (int i = 0; i < photonsPerAngle; i++)
        {
            APhoton photon;
            photon.v[0] = std::sin(theta);
            photon.v[1] = 0;
            photon.v[2] = std::cos(theta);

            const AInterfaceRule::EInterfaceRuleResult outcome = rule.calculate(&photon, normal);
            if (outcome == AInterfaceRule::Absorbed)
            {
                numAbsorption++;
                if (rule.Status != AInterfaceRule::Absorption) statusErrors++;
                continue;
            }
            if (outcome != AInterfaceRule::Back && outcome != AInterfaceRule::Forward)
            {
                numErrors++;
                continue;
            }

            const bool reflected = (outcome == AInterfaceRule::Back);
            if (reflected && rule.Status != AInterfaceRule::LobeReflection) statusErrors++;
            if (!reflected && rule.Status != AInterfaceRule::Transmission) statusErrors++;

            const double norm = std::sqrt(photon.v[0]*photon.v[0] + photon.v[1]*photon.v[1] + photon.v[2]*photon.v[2]);
            maximumNormError = std::max(maximumNormError, std::abs(norm-1.0));
            const double dotN = photon.v[0]*normal[0] + photon.v[1]*normal[1] + photon.v[2]*normal[2];
            if ((reflected && dotN >= 0) || (!reflected && dotN <= 0)) directionErrors++;

            double cosOut = reflected ? -dotN : dotN;
            cosOut = std::clamp(cosOut, 0.0, 1.0);
            const double thetaOut = std::acos(cosOut)*180.0/M_PI;
            double phiOut = std::atan2(photon.v[0]*ey[0] + photon.v[1]*ey[1] + photon.v[2]*ey[2],
                                       photon.v[0]*ex[0] + photon.v[1]*ex[1] + photon.v[2]*ex[2])*180.0/M_PI;
            if (phiOut < 0) phiOut += 360.0;
            const int it = std::min(nTheta-1, (int)(thetaOut/90.0*nTheta));
            const int ip = std::min(nPhi-1, (int)(phiOut/360.0*nPhi));
            if (reflected)
            {
                numReflection++;
                measuredReflection[it*nPhi + ip]++;
            }
            else
            {
                numTransmission++;
                measuredTransmission[it*nPhi + ip]++;
            }
        }

        std::vector<double> expectedReflection, expectedTransmission;
        double expectedR, expectedT;
        expectedDistribution(angle, true,  expectedReflection,   expectedR);
        expectedDistribution(angle, false, expectedTransmission, expectedT);
        const double expectedA = expectedAbsorption(angle);
        const std::vector<double> measuredR = normalizeCounts(measuredReflection, numReflection);
        const std::vector<double> measuredT = normalizeCounts(measuredTransmission, numTransmission);
        const std::vector<double> expectedRTheta = thetaMarginal(expectedReflection);
        const std::vector<double> expectedTTheta = thetaMarginal(expectedTransmission);
        const std::vector<double> measuredRTheta = thetaMarginal(measuredR);
        const std::vector<double> measuredTTheta = thetaMarginal(measuredT);

        const double measuredRProb = numReflection/(double)photonsPerAngle;
        const double measuredTProb = numTransmission/(double)photonsPerAngle;
        const double measuredAProb = numAbsorption/(double)photonsPerAngle;
        const double rThetaTvd = expectedR > 0 ? totalVariation(measuredRTheta, expectedRTheta) : 0;
        const double tThetaTvd = expectedT > 0 ? totalVariation(measuredTTheta, expectedTTheta) : 0;
        // At normal incidence the incidence plane (and therefore absolute phi=0) is undefined.
        // R/T and theta remain testable there; full theta/phi comparison starts at non-zero angle.
        const bool full2dApplicable = std::sin(theta) >= 1.0e-6;
        const double rFullTvd = expectedR > 0 && full2dApplicable ? totalVariation(measuredR, expectedReflection) : 0;
        const double tFullTvd = expectedT > 0 && full2dApplicable ? totalVariation(measuredT, expectedTransmission) : 0;
        const double rThetaLimit = expectedR > 0 ? std::max(0.004, 3.0*stochasticNoiseTvd(expectedRTheta, numReflection)) : 0;
        const double tThetaLimit = expectedT > 0 ? std::max(0.004, 3.0*stochasticNoiseTvd(expectedTTheta, numTransmission)) : 0;
        const double rFullLimit = expectedR > 0 && full2dApplicable ? std::max(0.01, 3.0*stochasticNoiseTvd(expectedReflection, numReflection)) : 0;
        const double tFullLimit = expectedT > 0 && full2dApplicable ? std::max(0.01, 3.0*stochasticNoiseTvd(expectedTransmission, numTransmission)) : 0;

        const bool passed = std::abs(measuredRProb-expectedR) <= probabilityTolerance(expectedR) &&
                            std::abs(measuredTProb-expectedT) <= probabilityTolerance(expectedT) &&
                            std::abs(measuredAProb-expectedA) <= probabilityTolerance(expectedA) &&
                            rThetaTvd <= rThetaLimit && tThetaTvd <= tThetaLimit &&
                            rFullTvd <= rFullLimit && tFullTvd <= tFullLimit &&
                            numErrors == 0 && directionErrors == 0 && statusErrors == 0;
        allPassed = allPassed && passed;
        totalDirectionErrors += directionErrors;
        totalStatusErrors += statusErrors;

        QVariantMap angleResult;
        angleResult["angle"] = angle;
        angleResult["passed"] = passed;
        angleResult["reflected"] = numReflection;
        angleResult["transmitted"] = numTransmission;
        angleResult["absorbed"] = numAbsorption;
        angleResult["errors"] = numErrors;
        angleResult["directionErrors"] = directionErrors;
        angleResult["statusErrors"] = statusErrors;
        angleResult["measuredR"] = measuredRProb;
        angleResult["expectedR"] = expectedR;
        angleResult["measuredT"] = measuredTProb;
        angleResult["expectedT"] = expectedT;
        angleResult["measuredA"] = measuredAProb;
        angleResult["expectedA"] = expectedA;
        angleResult["reflectedThetaTvd"] = rThetaTvd;
        angleResult["reflectedThetaTvdLimit"] = rThetaLimit;
        angleResult["transmittedThetaTvd"] = tThetaTvd;
        angleResult["transmittedThetaTvdLimit"] = tThetaLimit;
        angleResult["full2dApplicable"] = full2dApplicable;
        angleResult["reflected2dTvd"] = rFullTvd;
        angleResult["reflected2dTvdLimit"] = rFullLimit;
        angleResult["transmitted2dTvd"] = tFullTvd;
        angleResult["transmitted2dTvdLimit"] = tFullLimit;
        angleResult["measuredReflectedTheta"] = vectorToVariants(measuredRTheta);
        angleResult["expectedReflectedTheta"] = vectorToVariants(expectedRTheta);
        angleResult["measuredTransmittedTheta"] = vectorToVariants(measuredTTheta);
        angleResult["expectedTransmittedTheta"] = vectorToVariants(expectedTTheta);
        if (include2D)
        {
            angleResult["measuredReflected2D"] = vectorToVariants(measuredR);
            angleResult["expectedReflected2D"] = vectorToVariants(expectedReflection);
            angleResult["measuredTransmitted2D"] = vectorToVariants(measuredT);
            angleResult["expectedTransmitted2D"] = vectorToVariants(expectedTransmission);
        }
        angleResults.push_back(angleResult);
    }

    // Rotational covariance: with the same random sequence, rotating the normal and incidence
    // frame must leave the recovered (theta,phi,outcome) bin sequence unchanged.
    const int covariancePhotons = std::min(50000, photonsPerAngle);
    auto makeCovarianceSignature = [&rule, &randomHub, nTheta, nPhi, seed](const double * normal, const double * tangent,
                                                                          int numPhotons)
    {
        std::vector<int> signature;
        signature.reserve(numPhotons);
        const double angle = 45.0*M_PI/180.0;
        double ey[3] = {normal[1]*tangent[2]-normal[2]*tangent[1],
                        normal[2]*tangent[0]-normal[0]*tangent[2],
                        normal[0]*tangent[1]-normal[1]*tangent[0]};
        randomHub.setSeed(seed + 7919);
        for (int i = 0; i < numPhotons; i++)
        {
            APhoton photon;
            for (int d = 0; d < 3; d++) photon.v[d] = std::sin(angle)*tangent[d] + std::cos(angle)*normal[d];
            const AInterfaceRule::EInterfaceRuleResult outcome = rule.calculate(&photon, normal);
            if (outcome == AInterfaceRule::Absorbed) { signature.push_back(-1); continue; }
            if (outcome != AInterfaceRule::Back && outcome != AInterfaceRule::Forward) { signature.push_back(-2); continue; }
            const bool reflected = outcome == AInterfaceRule::Back;
            double dotN = 0, dotX = 0, dotY = 0;
            for (int d = 0; d < 3; d++)
            {
                dotN += photon.v[d]*normal[d];
                dotX += photon.v[d]*tangent[d];
                dotY += photon.v[d]*ey[d];
            }
            const double thetaOut = std::acos(std::clamp(reflected ? -dotN : dotN, 0.0, 1.0))*180.0/M_PI;
            double phiOut = std::atan2(dotY, dotX)*180.0/M_PI;
            if (phiOut < 0) phiOut += 360.0;
            const int it = std::min(nTheta-1, (int)(thetaOut/90.0*nTheta));
            const int ip = std::min(nPhi-1, (int)(phiOut/360.0*nPhi));
            signature.push_back((reflected ? 0 : nTheta*nPhi) + it*nPhi + ip);
        }
        return signature;
    };

    const double baseNormal[3] = {0, 0, 1};
    const double baseTangent[3] = {1, 0, 0};
    double rotatedNormal[3] = {0.2, -0.3, 0.9327379053088815};
    const double normalLength = std::sqrt(rotatedNormal[0]*rotatedNormal[0] + rotatedNormal[1]*rotatedNormal[1] + rotatedNormal[2]*rotatedNormal[2]);
    for (double & v : rotatedNormal) v /= normalLength;
    double rotatedTangent[3] = {1.0-rotatedNormal[0]*rotatedNormal[0],
                                -rotatedNormal[0]*rotatedNormal[1],
                                -rotatedNormal[0]*rotatedNormal[2]};
    const double tangentLength = std::sqrt(rotatedTangent[0]*rotatedTangent[0] + rotatedTangent[1]*rotatedTangent[1] + rotatedTangent[2]*rotatedTangent[2]);
    for (double & v : rotatedTangent) v /= tangentLength;
    const std::vector<int> baseSignature = makeCovarianceSignature(baseNormal, baseTangent, covariancePhotons);
    const std::vector<int> rotatedSignature = makeCovarianceSignature(rotatedNormal, rotatedTangent, covariancePhotons);
    qlonglong covarianceMismatches = 0;
    for (int i = 0; i < covariancePhotons; i++)
        if (baseSignature[i] != rotatedSignature[i]) covarianceMismatches++;
    const bool covariancePassed = covarianceMismatches <= std::max(2, covariancePhotons/100000);
    allPassed = allPassed && covariancePassed;

    // A minimal synthetic LUT exercises a non-zero absorption branch, which generated v1 LUTs
    // normally do not contain.
    ALutInterfaceRule absorptionRule(-1, -1);
    ALutSurfaceData & synthetic = absorptionRule.Data;
    synthetic.ThetaIncBins = 1;
    synthetic.ThetaOutBins = 1;
    synthetic.PhiOutBins = 1;
    synthetic.Launched = {10};
    synthetic.ReflectedCounts = {2};
    synthetic.TransmittedCounts = {3};
    synthetic.AbsorbedCounts = {5};
    synthetic.ReflectedHist = {{2}};
    synthetic.TransmittedHist = {{3}};
    err = synthetic.buildRuntime();
    if (!err.isEmpty())
    {
        abort("validateSurfaceLut: internal absorption test setup failed: " + err);
        return QVariantMap();
    }
    const int absorptionPhotons = std::min(100000, photonsPerAngle);
    qlonglong syntheticR = 0, syntheticT = 0, syntheticA = 0;
    qlonglong syntheticErrors = 0, syntheticStatusErrors = 0;
    randomHub.setSeed(seed + 15485863);
    const double normal[3] = {0, 0, 1};
    for (int i = 0; i < absorptionPhotons; i++)
    {
        APhoton photon;
        photon.v[0] = 0; photon.v[1] = 0; photon.v[2] = 1;
        const AInterfaceRule::EInterfaceRuleResult outcome = absorptionRule.calculate(&photon, normal);
        if (outcome == AInterfaceRule::Back)
        {
            syntheticR++;
            if (absorptionRule.Status != AInterfaceRule::LobeReflection) syntheticStatusErrors++;
        }
        else if (outcome == AInterfaceRule::Forward)
        {
            syntheticT++;
            if (absorptionRule.Status != AInterfaceRule::Transmission) syntheticStatusErrors++;
        }
        else if (outcome == AInterfaceRule::Absorbed)
        {
            syntheticA++;
            if (absorptionRule.Status != AInterfaceRule::Absorption) syntheticStatusErrors++;
        }
        else syntheticErrors++;
    }
    auto syntheticTolerance = [absorptionPhotons](double p)
    {
        return std::max(0.002, 6.0*std::sqrt(p*(1.0-p)/absorptionPhotons));
    };
    const bool absorptionPassed = syntheticErrors == 0 && syntheticStatusErrors == 0 &&
                                  std::abs(syntheticR/(double)absorptionPhotons-0.2) <= syntheticTolerance(0.2) &&
                                  std::abs(syntheticT/(double)absorptionPhotons-0.3) <= syntheticTolerance(0.3) &&
                                  std::abs(syntheticA/(double)absorptionPhotons-0.5) <= syntheticTolerance(0.5);
    allPassed = allPassed && absorptionPassed;

    QVariantMap covarianceResult;
    covarianceResult["passed"] = covariancePassed;
    covarianceResult["photons"] = covariancePhotons;
    covarianceResult["mismatches"] = covarianceMismatches;
    QVariantMap absorptionResult;
    absorptionResult["passed"] = absorptionPassed;
    absorptionResult["photons"] = absorptionPhotons;
    absorptionResult["measuredR"] = syntheticR/(double)absorptionPhotons;
    absorptionResult["measuredT"] = syntheticT/(double)absorptionPhotons;
    absorptionResult["measuredA"] = syntheticA/(double)absorptionPhotons;
    absorptionResult["errors"] = syntheticErrors;
    absorptionResult["statusErrors"] = syntheticStatusErrors;

    result["passed"] = allPassed;
    result["lutFile"] = lutFile;
    result["photonsPerAngle"] = photonsPerAngle;
    result["seed"] = seed;
    result["thetaOutBins"] = nTheta;
    result["phiOutBins"] = nPhi;
    result["maximumNormError"] = maximumNormError;
    result["directionErrors"] = totalDirectionErrors;
    result["statusErrors"] = totalStatusErrors;
    result["angles"] = angleResults;
    result["rotatedNormalCovariance"] = covarianceResult;
    result["syntheticAbsorption"] = absorptionResult;
    return result;
}

ALutInterfaceRule * AInterfaceRules_SI::makeLutRule(int matFrom, int matTo, const QString & lutFile)
{
    ALutInterfaceRule * rule = new ALutInterfaceRule(matFrom, matTo);

    QString err = rule->loadLUT(lutFile);
    if (err.isEmpty()) err = rule->checkOverrideData();
    if (!err.isEmpty())
    {
        delete rule;
        abort(err);
        return nullptr;
    }

    const QString warning = rule->getMaterialConsistencyWarning();
    if (!warning.isEmpty())
        AScriptHub::getInstance().outputText("WARNING - " + warning, Lang);

    return rule;
}

void AInterfaceRules_SI::setLutMaterialRule(QString matFrom, QString matTo, QString lutFile)
{
    const AMaterialHub & MatHub = AMaterialHub::getConstInstance();
    const int iFrom = MatHub.findMaterial(matFrom);
    if (iFrom == -1)
    {
        abort("Material not found: " + matFrom);
        return;
    }
    const int iTo = MatHub.findMaterial(matTo);
    if (iTo == -1)
    {
        abort("Material not found: " + matTo);
        return;
    }
    if (iFrom == iTo)
    {
        abort("Materials 'from' and 'to' should be different");
        return;
    }

    ALutInterfaceRule * rule = makeLutRule(iFrom, iTo, lutFile);
    if (!rule) return;

    RuleHub.setMaterialRule(iFrom, iTo, rule);
    RuleHub.announceRulesChanged();
}

void AInterfaceRules_SI::setLutVolumeRule(QString volFrom, QString volTo, QString lutFile)
{
    // material indices are needed for the consistency check; taken from the current geometry
    // if the volumes exist (the rule itself is applied by volume names at simulation time)
    int iFrom = 0;
    int iTo   = 0;
    AGeoObject * world = AGeometryHub::getConstInstance().World;
    const AGeoObject * objFrom = world->findObjectByName(volFrom);
    if (objFrom) iFrom = objFrom->Material;
    const AGeoObject * objTo = world->findObjectByName(volTo);
    if (objTo) iTo = objTo->Material;

    ALutInterfaceRule * rule = makeLutRule(iFrom, iTo, lutFile);
    if (!rule) return;

    RuleHub.setVolumeRule(TString(volFrom.toLatin1().data()), TString(volTo.toLatin1().data()), rule);
    RuleHub.announceRulesChanged();
}

void AInterfaceRules_SI::clearMaterialRule(QString matFrom, QString matTo)
{
    const AMaterialHub & MatHub = AMaterialHub::getConstInstance();
    const int iFrom = MatHub.findMaterial(matFrom);
    if (iFrom == -1)
    {
        abort("Material not found: " + matFrom);
        return;
    }
    const int iTo = MatHub.findMaterial(matTo);
    if (iTo == -1)
    {
        abort("Material not found: " + matTo);
        return;
    }

    RuleHub.setMaterialRule(iFrom, iTo, nullptr);
    RuleHub.announceRulesChanged();
}

void AInterfaceRules_SI::clearVolumeRule(QString volFrom, QString volTo)
{
    RuleHub.removeVolumeRule(TString(volFrom.toLatin1().data()), TString(volTo.toLatin1().data()));
    RuleHub.announceRulesChanged();
}
