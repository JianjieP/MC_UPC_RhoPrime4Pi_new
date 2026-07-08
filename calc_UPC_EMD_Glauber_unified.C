/*
 * calc_UPC_EMD_Glauber_unified.C
 * Unified Pb+Pb / Au+Au UPC EMD probability table macro.
 */

#include "TMath.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double alpha = 1.0 / 137.0359895; 
const double hbarc = 0.197327053;       // GeV*fm
const double m_N   = 0.938;             // GeV
const double pi    = TMath::Pi();

struct SystemParams {
    TString name;
    int Z;
    int A;
    double R_WS;
    double a_WS;
    double neutronSepThr;
    double gdrSigmaMax;
    double gdrMean;
    double gdrWidth;
    double gdrFactorXn;
    double b_min_lin;
    bool usePbPartialGraphs;
};

SystemParams GetSystemParams(const char* system) {
    TString sys(system);
    sys.ToLower();
    if (sys == "auau" || sys == "au") {
        return {"AuAu", 79, 197, 6.38, 0.535, 8.1, 540.0, 13.70, 4.75, 1.0, 0.0, false};
    }
    if (sys != "pbpb" && sys != "pb") {
        cout << "WARNING: Unknown system '" << system << "', using PbPb." << endl;
    }
    return {"PbPb", 82, 208, 6.62, 0.546, 7.4, 640.0, 13.42, 4.05, 0.93, 12.0, true};
}

// =========================================================================
// Glauber 计算器：处理 Woods-Saxon 分布与无强子碰撞生存概率
// =========================================================================
class GlauberCalc {
public:
    double R_WS, a_WS, sigma_NN_fm2;
    int A;
    bool usePbGlauberGrid;
    vector<double> s_grid, TA_grid;
    vector<double> b_grid, Psurv_grid;

    GlauberCalc(const SystemParams& params, double sig_mb = 68.0) {
        A = params.A; R_WS = params.R_WS; a_WS = params.a_WS;
        usePbGlauberGrid = params.usePbPartialGraphs;
        sigma_NN_fm2 = sig_mb * 0.1; // 1 mb = 0.1 fm^2
        Init();
    }

    double rho(double r) { return 1.0 / (1.0 + exp((r - R_WS) / a_WS)); }

    void Init() {
        double norm = 0, dr = 0.01;
        for (double r = 0; r < 25.0; r += dr) norm += rho(r) * 4.0 * pi * r * r * dr;
        double C = 1.0 / norm;

        double ds = usePbGlauberGrid ? 0.05 : 0.02;
        for (double s = 0; s <= 25.0; s += ds) {
            double dz = usePbGlauberGrid ? 0.02 : 0.01, val = 0;
            for (double z = 0; z < 25.0; z += dz) {
                double r = sqrt(s*s + z*z);
                val += C * rho(r) * dz;
            }
            s_grid.push_back(s);
            TA_grid.push_back(2.0 * val);
        }

        auto get_TA = [&](double s) {
            if (s >= 25.0) return 0.0;
            int idx = int(s / ds);
            if (static_cast<size_t>(idx) >= TA_grid.size() - 1) return 0.0;
            double w = (s - s_grid[idx]) / ds;
            return TA_grid[idx] * (1.0 - w) + TA_grid[idx+1] * w;
        };

        double db = usePbGlauberGrid ? 0.1 : 0.02;
        for (double b = 0; b <= 25.0; b += db) {
            double val_TAA = 0;
            double s_step = usePbGlauberGrid ? 0.1 : 0.05;
            double phi_step = usePbGlauberGrid ? pi / 50.0 : pi / 100.0; 
            for (double s = 0; s < 25.0; s += s_step) {
                double ta1 = get_TA(s);
                if (ta1 == 0) continue;
                for (double phi = 0; phi < 2.0 * pi; phi += phi_step) {
                    double s2 = sqrt(s*s + b*b - 2.0*s*b*cos(phi));
                    val_TAA += ta1 * get_TA(s2) * s * s_step * phi_step;
                }
            }
            b_grid.push_back(b);
            Psurv_grid.push_back(exp(-A * A * sigma_NN_fm2 * val_TAA));
        }
    }

    double GetSurvivalProb(double b) {
        if (b >= 20.0) return 1.0; 
        double db = usePbGlauberGrid ? 0.1 : 0.02;
        int idx = int(b / db);
        if (static_cast<size_t>(idx) >= Psurv_grid.size() - 1) return 1.0;
        double w = (b - b_grid[idx]) / db;
        return Psurv_grid[idx] * (1.0 - w) + Psurv_grid[idx+1] * w;
    }
};

// =========================================================================
// EMD 光核反应截面与通量计算器
// =========================================================================
class EMDCalculator {
public:
    int Znu, Anu;
    SystemParams params;
    vector<double> E_gamma_centers, E_gamma_widths;
    vector<double> sigma_Xn, sigma_1n;

    EMDCalculator(const SystemParams& params_in) : params(params_in) {
        Znu = params.Z; Anu = params.A;
        InitCrossSections();
    }

    void InitCrossSections() {
        double neutronSepThr = params.neutronSepThr; 
        double eStepLowE = 0.05;    
        double BrtWgnr_xsection = params.gdrSigmaMax; 
        double BrtWgnr_width = params.gdrWidth;
        double BrtWgnr_mean = params.gdrMean;
        double BrtWgnr_Con = 0.1 * BrtWgnr_width * BrtWgnr_width * BrtWgnr_xsection; 

        vector<double> temp_E_edge, temp_sig_Xn; 

        // 1. GDR & Lepretre & Carlos
        int nStepsLowE = int((25.0 - neutronSepThr) / eStepLowE) + 1;
        for (int i = 0; i < nStepsLowE; i++) {
            double E = i * eStepLowE + neutronSepThr;
            temp_E_edge.push_back(E);
            temp_sig_Xn.push_back(params.gdrFactorXn * BrtWgnr_Con * E * E / (pow(BrtWgnr_mean * BrtWgnr_mean - E * E, 2) + E * E * BrtWgnr_width * BrtWgnr_width));
        }

        double eLepretre[28] = {26.,28.,30.,32.,34.,36.,38.,40.,44.,46.,48.,50.,52.,55.,57.,62.,64.,66.,69.,72.,74.,76.,79.,82.,86.,92.,98.,103.};
        double sigLepretre[28] = {30.,21.5,22.5,18.5,17.5,15.,14.5,19.,17.5,16.,14.,20.,16.5,17.5,17.,15.5,18.,15.5,15.5,15.,13.5,18.,14.5,15.5,12.5,13.,13.,12.};
        for(int j=0; j<27; j++) { temp_E_edge.push_back(eLepretre[j]); temp_sig_Xn.push_back(0.1 * Anu * sigLepretre[j] / 208.0); }

        double eCarlos[22] = {103.,106.,112.,119.,127.,132.,145.,171.,199.,230.,235.,254.,280.,300.,320.,330.,333.,373.,390.,420.,426.,440.};
        double sigCarlos[22] = {12.0,11.5,12.0,12.0,12.0,15.0,17.0,28.0,33.0,52.0,60.0,70.0,76.0,85.0,86.0,89.0,89.0,75.0,76.0,69.0,59.0,61.0};
        for(int j=0; j<22; j++) { temp_E_edge.push_back(eCarlos[j]); temp_sig_Xn.push_back(0.1 * Anu * sigCarlos[j] / 208.0); }

        // 2. Armstrong & MichaDwell & Regge
        double xSection_p_Armstrong[160]={0.,.4245,.4870,.5269,.4778,.4066,.3341,.2444,.2245,.2005,.1783,.1769,.1869,.1940,.2117,.2226,.2327,.2395,.2646,.2790,.2756,.2607,.2447,.2211,.2063,.2137,.2088,.2017,.2050,.2015,.2121,.2175,.2152,.1917,.1911,.1747,.1650,.1587,.1622,.1496,.1486,.1438,.1556,.1468,.1536,.1544,.1536,.1468,.1535,.1442,.1515,.1559,.1541,.1461,.1388,.1565,.1502,.1503,.1454,.1389,.1445,.1425,.1415,.1424,.1432,.1486,.1539,.1354,.1480,.1443,.1435,.1491,.1435,.1380,.1317,.1445,.1375,.1449,.1359,.1383,.1390,.1361,.1286,.1359,.1395,.1327,.1387,.1431,.1403,.1404,.1389,.1410,.1304,.1363,.1241,.1284,.1299,.1325,.1343,.1387,.1328,.1444,.1334,.1362,.1302,.1338,.1339,.1304,.1314,.1287,.1404,.1383,.1292,.1436,.1280,.1326,.1321,.1268,.1278,.1243,.1239,.1271,.1213,.1338,.1287,.1343,.1231,.1317,.1214,.1370,.1232,.1301,.1348,.1294,.1278,.1227,.1218,.1198,.1193,.1342,.1323,.1248,.1220,.1139,.1271,.1224,.1347,.1249,.1163,.1362,.1236,.1462,.1356,.1198,.1419,.1324,.1288,.1336,.1335,.1266};
        double xSection_n_Armstrong[160]={0.,.3125,.3930,.4401,.4582,.3774,.3329,.2996,.2715,.2165,.2297,.1861,.1551,.2020,.2073,.2064,.2193,.2275,.2384,.2150,.2494,.2133,.2023,.1969,.1797,.1693,.1642,.1463,.1280,.1555,.1489,.1435,.1398,.1573,.1479,.1493,.1417,.1403,.1258,.1354,.1394,.1420,.1364,.1325,.1455,.1326,.1397,.1286,.1260,.1314,.1378,.1353,.1264,.1471,.1650,.1311,.1261,.1348,.1277,.1518,.1297,.1452,.1453,.1598,.1323,.1234,.1212,.1333,.1434,.1380,.1330,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12};
        for (int j = 9; j <= 70; j++) { temp_E_edge.push_back(temp_E_edge.back() + 25.0); temp_sig_Xn.push_back(0.1 * (Znu * xSection_p_Armstrong[j] + (Anu - Znu) * xSection_n_Armstrong[j])); }

        double eMichaDwell[11]={2000.0,3270.0,4100.0,4810.0,6210.0,6600.0,7790.0,8400.0,9510.0,13600.0,16400.0};
        double sigMichaDwell[11]={0.1266,0.1080,0.0805,0.1017,0.0942,0.0844,0.0841,0.0755,0.0827,0.0626,0.0740};
        for (int j = 0; j < 11; j++) { temp_E_edge.push_back(eMichaDwell[j]); temp_sig_Xn.push_back(0.1 * Anu * sigMichaDwell[j]); }

        double x=0.0677, y=0.129, eps=0.0808, eta=0.4525, em=0.94, exx=pow(10,0.05);
        double s = 0.002 * em * temp_E_edge.back();
        for (int j = 1; j <= 150; j++) {
            s = s * exx;
            temp_E_edge.push_back(1000.0 * 0.5 * (s - em*em) / em);
            temp_sig_Xn.push_back(0.1 * 0.65 * Anu * (x * pow(s, eps) + y * pow(s, -eta)));
        }

        if (params.usePbPartialGraphs) {
        // =========================================================================
        // ★ 严格全量引入 Pb 专属的 1n-10n 分部截面图！
        // =========================================================================
        TGraph* gSection_Nn[10];
        
        double e1n[76]={7.5,7.57,7.64,7.71,7.79,7.88,7.95,8.05,8.12,8.19,8.26,8.33,8.47,8.56,8.65,8.74,8.87,8.98,9.12,9.23,9.35,9.49,9.6,9.74,9.84,9.96,10.08,10.22,10.37,10.49,10.61,10.73,10.85,11.0,11.13,11.24,11.38,11.51,11.64,11.76,11.89,12.0,12.13,12.26,12.53,12.68,12.95,13.12,13.5,13.79,14.05,14.28,14.58,14.8,15.08,15.37,15.66,15.95,16.21,16.53,16.79,17.05,17.33,17.6,17.86,18.13,18.4,18.66,18.9,19.2,19.46,19.74,20.03,20.28,20.55,20.85};
        double s1n[76]={11.0,21.0,34.0,31.0,21.0,19.0,32.0,36.0,29.0,25.0,26.0,35.0,31.0,34.0,32.0,38.0,38.0,46.0,54.0,60.0,71.0,73.0,67.0,80.0,115.0,134.0,113.0,101.0,132.0,164.0,176.0,183.0,208.0,241.0,306.0,311.0,302.0,321.0,336.0,386.0,407.0,418.0,432.0,459.0,512.0,534.0,589.0,641.0,645.0,625.0,598.0,560.0,470.0,418.0,336.0,287.0,230.0,184.0,137.0,116.0,95.0,82.0,71.0,59.0,53.0,46.0,42.0,35.0,31.0,27.0,26.0,25.0,25.0,23.0,23.0,19.0};
        gSection_Nn[0] = new TGraph(76, e1n, s1n);
        for(int i=0; i<76; i++) gSection_Nn[0]->GetY()[i] *= 0.1 * 0.93;

        double e2nV[50]={14.587,14.758,14.915,15.007,15.092,15.281,15.386,15.602,15.836,16.077,16.316,16.590,16.801,17.005,17.226,17.447,17.682,17.933,18.115,18.303,18.504,18.822,19.067,19.308,19.545,19.776,20.116,20.344,20.575,20.736,20.951,21.187,21.471,21.685,21.930,22.270,22.500,22.790,23.030,23.255,23.486,23.711,23.919,24.151,24.375,24.600,24.808,25.040,25.475,25.708};
        double s2nV[50]={0.904,1.809,2.932,3.551,4.267,5.114,6.017,6.850,7.613,7.988,8.537,8.664,8.811,8.637,8.493,8.210,8.139,7.693,7.574,7.471,7.128,6.798,6.290,6.124,5.603,5.313,5.029,5.037,4.724,4.521,4.602,4.532,4.507,4.435,4.270,4.247,4.023,3.897,3.887,3.793,3.632,3.402,3.262,3.122,3.041,2.842,2.535,2.244,1.778,1.589};
        double e2nL[40]={26.001,26.141,26.149,26.435,27.091,27.098,27.885,29.022,29.927,31.763,34.364,37.141,42.241,45.281,48.289,51.868,55.422,58.917,62.668,64.688,68.408,73.211,76.919,84.075,87.430,91.308,92.320,95.786,99.933,102.943,106.520,107.693,111.102,115.628,120.750,123.759,127.338,128.881,131.890,135.306};
        double s2nL[40]={30.110,29.326,27.974,26.828,25.854,24.851,23.790,22.571,21.159,19.817,18.783,17.806,16.979,16.985,16.244,16.008,15.672,15.403,15.350,14.835,14.566,14.384,14.240,13.961,13.602,13.489,13.044,13.429,13.338,13.054,13.032,12.714,12.784,12.773,12.817,12.589,12.285,12.596,12.424,12.303};
        gSection_Nn[1] = new TGraph(90);
        for(int i=0; i<50; i++) gSection_Nn[1]->SetPoint(i, e2nV[i], s2nV[i] * 0.93);
        for(int i=0; i<40; i++) gSection_Nn[1]->SetPoint(50+i, e2nL[i], s2nL[i] * 0.1 / 2.0);

        double e3n[44]={22.515,22.994,23.148,23.305,23.623,23.834,23.935,24.255,24.733,25.308,26.009,26.813,28.595,30.444,32.659,34.611,37.150,37.875,40.211,43.416,46.641,49.378,53.185,56.672,60.200,63.925,67.504,70.432,74.660,77.751,82.468,85.720,90.274,93.852,98.730,101.750,106.210,109.535,114.341,117.593,122.471,125.561,129.951,133.366};
        double s3n[44]={3.124,4.637,5.964,6.951,8.114,9.123,10.271,11.182,12.646,14.889,16.752,18.192,19.338,20.086,19.558,18.536,17.781,17.352,16.573,15.312,14.790,13.715,13.306,12.812,12.410,11.949,11.669,11.510,11.398,11.123,10.971,10.946,10.807,10.716,10.696,10.655,10.583,10.581,10.459,10.508,10.456,10.483,10.460,10.457};
        gSection_Nn[2] = new TGraph(44, e3n, s3n);
        for(int i=0; i<44; i++) gSection_Nn[2]->GetY()[i] *= 0.1 / 3.0;

        double e4n[35]={31.473,32.280,33.087,34.054,35.344,35.826,36.469,37.276,38.830,40.193,45.066,48.583,52.221,55.706,59.377,62.954,66.393,69.840,73.689,80.844,84.421,88.091,91.576,94.502,98.556,102.261,105.979,110.765,114.927,118.281,121.124,123.813,127.902,131.101,134.206};
        double s4n[35]={0.838,1.718,2.775,4.111,5.857,6.894,8.109,9.103,10.259,10.645,11.446,11.615,11.435,11.333,11.200,11.156,11.056,10.838,10.695,10.645,10.732,10.687,10.580,10.654,10.450,10.129,9.925,10.244,10.346,10.204,9.874,10.085,9.906,9.961,9.851};
        gSection_Nn[3] = new TGraph(35, e4n, s4n);
        for(int i=0; i<35; i++) gSection_Nn[3]->GetY()[i] *= 0.1 / 4.0;

        double e5n[30]={39.137,40.408,41.704,43.161,45.105,46.920,48.867,52.326,55.655,60.369,63.782,68.496,71.423,75.975,79.551,82.153,86.705,90.282,95.160,98.737,102.856,105.891,109.288,113.046,115.973,120.525,123.858,127.678,131.256,134.345};
        double s5n[30]={0.699,2.029,2.950,3.920,5.131,6.049,6.819,7.618,8.208,8.513,8.683,8.941,9.056,9.201,9.292,9.331,9.456,9.495,9.542,9.559,9.546,9.547,9.563,9.551,9.553,9.710,9.725,9.845,9.779,9.814};
        gSection_Nn[4] = new TGraph(30, e5n, s5n);
        for(int i=0; i<30; i++) gSection_Nn[4]->GetY()[i] *= 0.1 / 5.0;

        double e6n[26]={46.139,47.999,50.162,54.058,57.885,61.042,65.158,68.191,71.674,75.342,78.268,83.144,86.558,91.922,95.498,99.660,102.651,105.414,110.129,113.380,118.258,121.834,125.690,128.988,132.686,135.166};
        double s6n[26]={0.363,1.617,2.388,3.561,4.416,4.891,5.348,5.642,5.942,6.245,6.460,6.792,6.952,7.209,7.404,7.573,7.668,7.809,7.977,8.074,8.197,8.329,8.371,8.440,8.559,8.605};
        gSection_Nn[5] = new TGraph(26, e6n, s6n);
        for(int i=0; i<26; i++) gSection_Nn[5]->GetY()[i] *= 0.1 / 6.0;

        double e7n[20]={58.466,62.292,65.342,69.187,73.225,76.173,80.885,84.094,88.036,91.449,96.487,100.063,105.264,108.352,113.066,116.154,120.998,124.281,129.321,132.247};
        double s7n[20]={0.599,1.326,1.871,2.467,3.056,3.389,3.957,4.275,4.530,4.836,5.259,5.484,5.818,6.036,6.301,6.485,6.762,6.947,7.157,7.282};
        gSection_Nn[6] = new TGraph(20, e7n, s7n);
        for(int i=0; i<20; i++) gSection_Nn[6]->GetY()[i] *= 0.1 / 7.0;

        double e8n[19]={66.435,71.147,75.264,77.972,82.522,85.512,89.673,93.443,96.824,102.350,105.113,109.248,112.264,116.652,119.578,124.551,127.380,131.736,133.718};
        double s8n[19]={0.369,0.959,1.399,1.639,2.095,2.379,2.716,3.038,3.315,3.730,3.971,4.309,4.445,4.779,5.003,5.356,5.515,5.802,5.962};
        gSection_Nn[7] = new TGraph(19, e8n, s8n);
        for(int i=0; i<19; i++) gSection_Nn[7]->GetY()[i] *= 0.1 / 8.0;

        double e9n[18]={74.077,78.116,80.740,85.127,88.432,92.279,95.854,99.431,103.192,106.420,111.135,114.061,118.288,121.214,125.441,128.652,132.594,134.871};
        double s9n[18]={0.388,0.909,1.121,1.579,1.788,2.070,2.303,2.490,2.698,2.873,3.054,3.130,3.281,3.356,3.489,3.536,3.653,3.659};
        gSection_Nn[8] = new TGraph(18, e9n, s9n);
        for(int i=0; i<18; i++) gSection_Nn[8]->GetY()[i] *= 0.1 / 9.0;

        double e10n[16]={81.558,85.134,88.709,92.285,95.861,99.437,102.986,106.589,110.165,113.742,117.318,120.895,124.471,128.047,131.624,134.551};
        double s10n[16]={0.186,0.418,0.760,1.048,1.295,1.447,1.642,1.851,2.036,2.162,2.269,2.431,2.607,2.691,2.775,2.831};
        gSection_Nn[9] = new TGraph(16, e10n, s10n);
        for(int i=0; i<16; i++) gSection_Nn[9]->GetY()[i] *= 0.01;

        // 生成网格并提取真实的 ratio_1n
        for(size_t i=0; i<temp_E_edge.size()-1; i++) {
            double E_center = (temp_E_edge[i] + temp_E_edge[i+1]) / 2.0;
            E_gamma_widths.push_back((temp_E_edge[i+1] - temp_E_edge[i]) / 1000.0); 
            E_gamma_centers.push_back(E_center / 1000.0); 
            sigma_Xn.push_back(temp_sig_Xn[i]);

            // 1. 求出 tcross (1n 到 10n 的物理总和)
            double tcross = 0;
            for(int j=0; j<10; j++) {
                double e_min = gSection_Nn[j]->GetX()[0];
                double e_max = gSection_Nn[j]->GetX()[gSection_Nn[j]->GetN()-1];
                if(E_center > e_min && E_center < e_max) {
                    tcross += gSection_Nn[j]->Eval(E_center);
                }
            }

            // 2. 判定分母 tm (还原 drawcross.C 中对 hXsectionXn_T 的覆盖逻辑)
            double tm = tcross;
            if (E_center > 20.85 && E_center < 24.0) tm = temp_sig_Xn[i];
            if (E_center > 130.0) tm = temp_sig_Xn[i];

            // 3. 计算 1n 占比 ratio_1n
            double ratio_1n = 0.0;
            double e_min_1n = gSection_Nn[0]->GetX()[0];
            double e_max_1n = gSection_Nn[0]->GetX()[gSection_Nn[0]->GetN()-1];
            
            if (E_center > e_min_1n && E_center < e_max_1n) {
                ratio_1n = gSection_Nn[0]->Eval(E_center) / tm;
            } else if (E_center > 20.85 && E_center < 24.0) {
                // 原码中在此区间 ratio_1n = 1.0 - ratio_2n
                ratio_1n = 1.0 - (gSection_Nn[1]->Eval(E_center) / tm);
            }
            
            // 安全保护
            if (ratio_1n < 0) ratio_1n = 0.0;
            if (ratio_1n > 1.0) ratio_1n = 1.0;

            // 存入绝对 1n 截面
            sigma_1n.push_back(temp_sig_Xn[i] * ratio_1n); 
        }

        for(int j=0; j<10; j++) delete gSection_Nn[j];
        } else {
            double sa[160]={0.,.004,.008,.013,.017,.021,.025,.029,.034,.038,.042,.046,.051,.055,.059,.063,.067,.072,.076,.08,.085,.09,.095,.1,.108,.116,
            .124,.132,.14,.152,.164,.176,.188,.2,.22,.24,.26,.28,.3,.32,.34,.36,.38,.4,.417,.433,.450,.467,.483,.5,.51,.516,.52,.523,.5245,.525,.5242,
            .5214,.518,.512,.505,.495,.482,.469,.456,.442,.428,.414,.4,.386,.370,.355,.34,.325,.310,.295,.280,.265,.25,.236,.222,.208,.194,.180,.166,
            .152,.138,.124,.11,.101,.095,.09,.085,.08,.076,.072,.069,.066,.063,.06,.0575,.055,.0525,.05,.04875,.0475,.04625,.045,.04375,.0425,.04125,.04,.03875,.0375,.03625,.035,.03375,.0325,.03125,.03,
            .02925,.0285,.02775,.027,.02625,.0255,.02475,.024,.02325,.0225,.02175,.021,.02025,.0195,.01875,.018,.01725,.0165,.01575,.015,
            .01425,.0135,.01275,.012,.01125,.0105,.00975,.009,.00825,.0075,.00675,.006,.00525,.0045,.00375,.003,.00225,.0015,.00075,0.};
            TGraph* g_1n_orig = new TGraph();
            for(int i=0; i<160; i++) { g_1n_orig->SetPoint(i, 8.1 + i*0.1, sa[i] * 100.0); }

            for(size_t i=0; i<temp_E_edge.size()-1; i++) {
                double E_center = (temp_E_edge[i] + temp_E_edge[i+1]) / 2.0;
                E_gamma_widths.push_back((temp_E_edge[i+1] - temp_E_edge[i]) / 1000.0);
                E_gamma_centers.push_back(E_center / 1000.0);
                sigma_Xn.push_back(temp_sig_Xn[i]);

                if(E_center > 24.0) {
                    sigma_1n.push_back(0.0);
                } else {
                    double sig1n_val = g_1n_orig->Eval(E_center);
                    if(sig1n_val > temp_sig_Xn[i]) sig1n_val = temp_sig_Xn[i];
                    sigma_1n.push_back(sig1n_val);
                }
            }
            delete g_1n_orig;
        }
    }

    double GetPhotonFlux(double b_fm, double E_gamma_GeV, double Gamma_target_frame) {
        if (b_fm < 1e-3) return 0.0;
        double X = (b_fm * E_gamma_GeV) / (Gamma_target_frame * hbarc);
        if (X > 50) return 0.0; 
        double factor1 = (double(Znu * Znu) * alpha) / (pi * pi);  
        double bessel  = TMath::BesselK1(X);
        return factor1 * (1.0 / (E_gamma_GeV * b_fm * b_fm)) * X * X * bessel * bessel;
    }

    void GetProbabilities(double b, double Gamma_target_frame, double &P_1n, double &P_Xn) {
        double mu_Xn = 0, mu_1n = 0;
        for(size_t i=0; i<E_gamma_centers.size(); i++) {
            double flux = GetPhotonFlux(b, E_gamma_centers[i], Gamma_target_frame);
            mu_Xn += flux * sigma_Xn[i] * E_gamma_widths[i];
            mu_1n += flux * sigma_1n[i] * E_gamma_widths[i];
        }
        double P_0n = exp(-mu_Xn);
        P_Xn = 1.0 - P_0n;                  
        P_1n = mu_1n * P_0n;                
    }
};

void calc_UPC_EMD_Glauber_unified(
    const char* system = "PbPb",
    double sqrt_s_NN_GeV = 5360.0,
    double sigma_NN_mb = 68.0,
    const char* outputFileName = ""
) {
    SystemParams params = GetSystemParams(system);
    if (params.name == "AuAu" && sigma_NN_mb == 68.0) {
        sigma_NN_mb = 42.0;
    }

    double gamma_c = sqrt_s_NN_GeV / (2.0 * m_N);
    double Gamma_target_frame = 2.0 * gamma_c * gamma_c - 1.0; 
    
    EMDCalculator emdCalc(params);
    GlauberCalc glauber(params, sigma_NN_mb);

    cout << "=========================================================" << endl;
    cout << " UPC EMD Calculator - " << params.name << " @ " << sqrt_s_NN_GeV << " GeV [UNIFIED]" << endl;
    cout << " Woods-Saxon: R=" << params.R_WS << " fm, a=" << params.a_WS << " fm" << endl;
    cout << " sigma_NN = " << sigma_NN_mb << " mb" << endl;
    cout << "=========================================================" << endl;

    int n_bins_b_lin = 20000;
    double b_min_lin = params.b_min_lin, b_max_lin = 1000.0; 
    double db_lin = (b_max_lin - b_min_lin) / n_bins_b_lin;

    TString outName = outputFileName;
    if (outName.IsNull()) {
        outName = Form("UPC_Probabilities_%s_%.0fGeV_Glauber_Final.root", params.name.Data(), sqrt_s_NN_GeV);
    }
    TFile *fOut = new TFile(outName, "RECREATE");
    
    // --- 独立基础概率直方图 ---
    TH1D *hProb_Surv     = new TH1D("hProb_Surv", "Probability of no hadronic collision; b (fm); P_{surv}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_1n   = new TH1D("hProb_EMD_1n", "Pure EMD SCD 1n Prob; b (fm); P_{1n}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_Xn   = new TH1D("hProb_EMD_Xn", "Pure EMD SCD Xn Prob; b (fm); P_{Xn}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_1n1n = new TH1D("hProb_EMD_1n1n", "Pure EMD MCD 1n1n Prob; b (fm); P_{1n}*P_{1n}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_XnXn = new TH1D("hProb_EMD_XnXn", "Pure EMD MCD XnXn Prob; b (fm); P_{Xn}*P_{Xn}", n_bins_b_lin, b_min_lin, b_max_lin);

    // --- 物理乘积直方图 (兼容旧有干涉代码) ---
    TH1D *hProb_MCD_1n1n = new TH1D("hProb_MCD_1n1n", "UPC MCD 1n1n Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_MCD_XnXn = new TH1D("hProb_MCD_XnXn", "UPC MCD XnXn Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_SCD_1n = new TH1D("hProb_SCD_1n", "UPC SCD 1n Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_SCD_Xn = new TH1D("hProb_SCD_Xn", "UPC SCD Xn Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);

    double total_sig_XnXn = 0, total_sig_1n1n = 0;
    double total_sig_SCD_Xn = 0, total_sig_SCD_1n = 0;

    // === 阶段 1：高精度制表 (供 ROOT 输出) ===
    for(int i=0; i<n_bins_b_lin; i++) {
        double b = b_min_lin + i*db_lin + db_lin/2.0; 
        double P_1n, P_Xn;
        emdCalc.GetProbabilities(b, Gamma_target_frame, P_1n, P_Xn);
        
        double P_surv = glauber.GetSurvivalProb(b);

        // 填充独立的成分直方图
        hProb_Surv->SetBinContent(i+1, P_surv);
        hProb_EMD_1n->SetBinContent(i+1, P_1n);
        hProb_EMD_Xn->SetBinContent(i+1, P_Xn);
        hProb_EMD_1n1n->SetBinContent(i+1, P_1n * P_1n);
        hProb_EMD_XnXn->SetBinContent(i+1, P_Xn * P_Xn);

        // 填充物理乘积直方图
        hProb_SCD_1n->SetBinContent(i+1, P_1n * P_surv);
        hProb_SCD_Xn->SetBinContent(i+1, P_Xn * P_surv);
        hProb_MCD_1n1n->SetBinContent(i+1, (P_1n * P_1n) * P_surv);
        hProb_MCD_XnXn->SetBinContent(i+1, (P_Xn * P_Xn) * P_surv);

        double area_element = 2.0 * pi * b * db_lin; 
        total_sig_XnXn += area_element * (P_Xn * P_Xn) * P_surv; 
        total_sig_1n1n += area_element * (P_1n * P_1n) * P_surv; 
        total_sig_SCD_Xn += area_element * P_Xn * P_surv;
        total_sig_SCD_1n += area_element * P_1n * P_surv;
    }

    // === 阶段 2：长程对数积分 (补齐 SCD 的宏观长程拖尾截面) ===
    int n_bins_b_log = 500;
    double b_min_log = 1000.0, b_max_log = 1.0e9; 
    double log_step = (log(b_max_log) - log(b_min_log)) / n_bins_b_log;

    for(int i=0; i<n_bins_b_log; i++) {
        double log_b_left = log(b_min_log) + i * log_step;
        double b_left = exp(log_b_left), b_right = exp(log_b_left + log_step);
        double b = (b_left + b_right) / 2.0; 
        
        double P_1n, P_Xn;
        emdCalc.GetProbabilities(b, Gamma_target_frame, P_1n, P_Xn);
        
        double area_element = 2.0 * pi * b * (b_right - b_left); 
        total_sig_XnXn += area_element * (P_Xn * P_Xn); 
        total_sig_1n1n += area_element * (P_1n * P_1n); 
        total_sig_SCD_Xn += area_element * P_Xn;
        total_sig_SCD_1n += area_element * P_1n;
    }

    // 写入独立成分
    hProb_Surv->Write();
    hProb_EMD_1n->Write(); hProb_EMD_Xn->Write();
    hProb_EMD_1n1n->Write(); hProb_EMD_XnXn->Write();
    
    // 写入乘积结果
    hProb_SCD_1n->Write(); hProb_SCD_Xn->Write();
    hProb_MCD_1n1n->Write(); hProb_MCD_XnXn->Write();
    fOut->Close();

    cout << "✅ 成功: UPC Probability 数组(包含基础概率拆分与乘积)已保存 -> " << outName << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "=> Mutual Coulomb Dissociation (双向 ZDC 巧合, MCD):" << endl;
    cout << "    Sigma(1n1n MCD) = " << total_sig_1n1n * 0.01 << " barns" << endl;
    cout << "    Sigma(XnXn MCD) = " << total_sig_XnXn * 0.01 << " barns" << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "=> Single Coulomb Dissociation (单侧 ZDC 触发, 本征 SCD 截面):" << endl;
    cout << "    Sigma(1n SCD)   = " << total_sig_SCD_1n * 0.01 << " barns" << endl;
    cout << "    Sigma(Xn SCD)   = " << total_sig_SCD_Xn * 0.01 << " barns" << endl;
    cout << "=========================================================" << endl;
}
