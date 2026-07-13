// MC validation of the DavisLUT runtime sampling. Fires photons at fixed incidence angles onto
// a LYSO->air normal, runs the ACTUAL runtime sampling routines (ALutSurfaceData::selectThetaBin
// + sampleOutgoing -- the code ALutInterfaceRule::calculate() calls, incl. the 90deg clamp),
// reconstructs the outgoing direction with the SAME formula as calculate(), then "detects" the
// reflected photons by recovering (theta_out, phi_out) and histogramming them -- the round trip a
// real reflected photon undergoes. Writes measured vs LUT-predicted distributions for comparison.
#include "alutsurfacedata.h"
#include "ajsontools.h"
#include <QJsonObject>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

int main(int argc, char** argv)
{
    const char* lutFile = argv[1];
    const char* outDir  = argv[2];
    const long  Nphot   = (argc > 3) ? atol(argv[3]) : 2000000;
    const double angles[3] = {10.0, 45.0, 80.0};

    ALutSurfaceData D;
    QJsonObject js;
    if (!jstools::loadJsonFromFile(js, lutFile)) { printf("cannot open %s\n", lutFile); return 1; }
    QString err = D.readFromJson(js);   if (!err.isEmpty()) { printf("read: %s\n", err.toLatin1().data()); return 1; }
    err = D.buildRuntime();             if (!err.isEmpty()) { printf("build: %s\n", err.toLatin1().data()); return 1; }

    const int nTo = D.ThetaOutBins, nPh = D.PhiOutBins, nTi = D.ThetaIncBins;
    const double N[3] = {0, 0, 1};
    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    for (int a = 0; a < 3; a++)
    {
        const double th = angles[a] * M_PI / 180.0;
        const double cosTi = cos(th), sinTi = sin(th);
        const double vin[3] = {sinTi, 0, cosTi};
        // incidence-plane frame, identical to ALutInterfaceRule::calculate(): ex = tangential of vin
        double ex[3] = {(vin[0]-cosTi*N[0])/sinTi, (vin[1]-cosTi*N[1])/sinTi, (vin[2]-cosTi*N[2])/sinTi};
        double ey[3] = {N[1]*ex[2]-N[2]*ex[1], N[2]*ex[0]-N[0]*ex[2], N[0]*ex[1]-N[1]*ex[0]};

        std::vector<long> meas(nTo * nPh, 0);
        long reflected = 0, transmitted = 0;

        for (long i = 0; i < Nphot; i++)
        {
            const int k = D.selectThetaBin(angles[a], U(gen));
            const double r = U(gen);
            const double R = D.getReflectionProbability(k), T = D.getTransmissionProbability(k);
            bool refl;
            if (r < R) refl = true; else if (r < R + T) refl = false; else continue;   // absorbed (v1: none)

            double thetaOut, phiOut;
            if (!D.sampleOutgoing(refl, k, U(gen), U(gen), U(gen), thetaOut, phiOut)) continue;

            const double tr = thetaOut * M_PI/180.0, pr = phiOut * M_PI/180.0;
            const double s = sin(tr), tn = cos(tr) * (refl ? -1.0 : 1.0);
            double v[3];
            for (int d = 0; d < 3; d++) v[d] = s*(cos(pr)*ex[d] + sin(pr)*ey[d]) + tn*N[d];

            if (refl)
            {
                reflected++;
                double cosOut = -(v[0]*N[0]+v[1]*N[1]+v[2]*N[2]);      // detector: angle from normal
                if (cosOut > 1) cosOut = 1; else if (cosOut < 0) cosOut = 0;
                double thOut = acos(cosOut) * 180.0/M_PI;
                double phOut = atan2(v[0]*ey[0]+v[1]*ey[1]+v[2]*ey[2],
                                     v[0]*ex[0]+v[1]*ex[1]+v[2]*ex[2]) * 180.0/M_PI;
                if (phOut < 0) phOut += 360.0;
                int it = std::min(nTo-1, (int)(thOut/90.0*nTo));
                int ip = std::min(nPh-1, (int)(phOut/360.0*nPh));
                meas[it*nPh + ip]++;
            }
            else transmitted++;
        }

        // LUT prediction = interpolation-weighted mix of the two adjacent incidence bins
        const double x = angles[a] / (90.0/nTi) - 0.5;
        int k0 = (int)floor(x); double f = x - k0;
        if (k0 < 0) { k0 = 0; f = 0; } if (k0 >= nTi-1) { k0 = nTi-1; f = 0; }
        const int k1 = (f > 0) ? k0+1 : k0;
        auto Rof = [&](int k){ return D.ReflectedCounts[k]/(double)D.Launched[k]; };
        const double R_lut = (1-f)*Rof(k0) + f*Rof(k1);
        std::vector<double> pred(nTo*nPh, 0);
        auto addbin = [&](int k, double w){
            const std::vector<int>& h = D.ReflectedHist[k]; double sm=0; for (int v:h) sm+=v;
            if (sm>0) for (size_t j=0;j<h.size();j++) pred[j] += w*h[j]/sm; };
        addbin(k0, 1-f); if (k1!=k0) addbin(k1, f);

        char fn[512]; snprintf(fn, sizeof(fn), "%s/valid_%02.0f.txt", outDir, angles[a]);
        FILE* fp = fopen(fn, "w");
        fprintf(fp, "# angle=%.1f nThetaOut=%d nPhiOut=%d R_meas=%.4f R_lut=%.4f reflected=%ld transmitted=%ld\n",
                angles[a], nTo, nPh, reflected/(double)Nphot, R_lut, reflected, transmitted);
        for (int it=0; it<nTo; it++) for (int ip=0; ip<nPh; ip++)
            fprintf(fp, "%d %d %.6e %.6e\n", it, ip, meas[it*nPh+ip]/(double)std::max(1L,reflected), pred[it*nPh+ip]);
        fclose(fp);
        printf("angle %2.0f: R_meas=%.4f R_lut=%.4f reflected=%ld/%ld -> %s\n",
               angles[a], reflected/(double)Nphot, R_lut, reflected, Nphot, fn);
    }
    return 0;
}
