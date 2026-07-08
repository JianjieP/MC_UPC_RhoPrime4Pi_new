/*
 * calc_UPC_EMD_AuAu_Glauber.C
 * Au+Au (RHIC 200 GeV) 高保真平滑完善版：
 * 1. 包含 Au 专属的光学 Glauber 模型 (R=6.38, a=0.535, sigNN=42mb)
 * 2. 直方图 Bin 密度提升 20 倍 (起点 0.0 fm，db = 0.05 fm)。
 * 3. 独立保存 P_surv、纯 EMD 概率以及它们的物理乘积。
 * 4. 修正 SCD 定义（不乘2），并保留宏观长程对数积分确保截面完整。
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

// 物理常数
const double alpha = 1.0 / 137.0359895; 
const double hbarc = 0.197327053;       // GeV*fm
const double m_N   = 0.938;             // GeV
const double pi    = TMath::Pi();

// =========================================================================
// Glauber 计算器：高精度数值积分版本
// =========================================================================
class GlauberCalc {
public:
    double R_WS, a_WS, sigma_NN_fm2;
    int A;
    vector<double> s_grid, TA_grid;
    vector<double> b_grid, Psurv_grid;

    // 默认参数设为 Au+Au @ 200 GeV
    GlauberCalc(int A_ = 197, double R_ = 6.38, double a_ = 0.535, double sig_mb = 42.0) {
        A = A_; R_WS = R_; a_WS = a_;
        sigma_NN_fm2 = sig_mb * 0.1; // 1 mb = 0.1 fm^2
        Init();
    }

    double rho(double r) { return 1.0 / (1.0 + exp((r - R_WS) / a_WS)); }

    void Init() {
        double norm = 0, dr = 0.01;
        for (double r = 0; r < 25.0; r += dr) norm += rho(r) * 4.0 * pi * r * r * dr;
        double C = 1.0 / norm;

        // 🌟 提升厚度函数积分精度
        double ds = 0.02; 
        for (double s = 0; s <= 25.0; s += ds) {
            double dz = 0.01, val = 0;
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
            if (idx >= TA_grid.size() - 1) return 0.0;
            double w = (s - s_grid[idx]) / ds;
            return TA_grid[idx] * (1.0 - w) + TA_grid[idx+1] * w;
        };

        // 🌟 提升生存概率的网格精度
        double db = 0.02; 
        for (double b = 0; b <= 25.0; b += db) {
            double val_TAA = 0;
            double s_step = 0.05, phi_step = pi / 100.0; 
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
        double db = 0.02; // 与 Init 中的 db 保持一致
        int idx = int(b / db);
        if (idx >= Psurv_grid.size() - 1) return 1.0;
        double w = (b - b_grid[idx]) / db;
        return Psurv_grid[idx] * (1.0 - w) + Psurv_grid[idx+1] * w;
    }
};

// =========================================================================
// EMD 光核反应截面与通量计算器 (Au 专属)
// =========================================================================
class EMDCalculator {
public:
    int Znu, Anu;
    vector<double> E_gamma_edges, E_gamma_centers, E_gamma_widths;
    vector<double> sigma_Xn, sigma_1n;

    EMDCalculator(int Z_in = 79, int A_in = 197) {
        Znu = Z_in; Anu = A_in;
        InitCrossSections();
    }

    void InitCrossSections() {
        double neutronSepThr = 8.1; 
        double eStepLowE = 0.05;    
        
        // Au GDR Breit-Wigner parameters
        double BrtWgnr_xsection = 540.0; 
        double BrtWgnr_width = 4.75;
        double BrtWgnr_mean = 13.70;
        double BrtWgnr_Con = 0.1 * BrtWgnr_width * BrtWgnr_width * BrtWgnr_xsection; 

        vector<double> temp_E_edge, temp_sig_Xn; 

        // 1. GDR
        int nStepsLowE = int((25.0 - neutronSepThr) / eStepLowE) + 1;
        for (int i = 0; i < nStepsLowE; i++) {
            double E = i * eStepLowE + neutronSepThr;
            temp_E_edge.push_back(E);
            temp_sig_Xn.push_back(BrtWgnr_Con * E * E / (pow(BrtWgnr_mean * BrtWgnr_mean - E * E, 2) + E * E * BrtWgnr_width * BrtWgnr_width));
        }

        // 2. Lepretre et al.
        double eLepretre[28] = {26.,28.,30.,32.,34.,36.,38.,40.,44.,46.,48.,50.,52.,55.,57.,62.,64.,66.,69.,72.,74.,76.,79.,82.,86.,92.,98.,103.};
        double sigLepretre[28] = {30.,21.5,22.5,18.5,17.5,15.,14.5,19.,17.5,16.,14.,20.,16.5,17.5,17.,15.5,18.,15.5,15.5,15.,13.5,18.,14.5,15.5,12.5,13.,13.,12.};
        for(int j=0; j<27; j++) {
            temp_E_edge.push_back(eLepretre[j]);
            temp_sig_Xn.push_back(0.1 * Anu * sigLepretre[j] / 208.0);
        }

        // 3. Carlos et al.
        double eCarlos[22] = {103.,106.,112.,119.,127.,132.,145.,171.,199.,230.,235.,254.,280.,300.,320.,330.,333.,373.,390.,420.,426.,440.};
        double sigCarlos[22] = {12.0,11.5,12.0,12.0,12.0,15.0,17.0,28.0,33.0,52.0,60.0,70.0,76.0,85.0,86.0,89.0,89.0,75.0,76.0,69.0,59.0,61.0};
        for(int j=0; j<22; j++) {
            temp_E_edge.push_back(eCarlos[j]);
            temp_sig_Xn.push_back(0.1 * Anu * sigCarlos[j] / 208.0);
        }

        // 4. Armstrong et al.
        double xSection_p_Armstrong[160]={0.,.4245,.4870,.5269,.4778,.4066,.3341,.2444,.2245,.2005,
		.1783,.1769,.1869,.1940,.2117,.2226,.2327,.2395,.2646,.2790,.2756,.2607,.2447,.2211,.2063,.2137,.2088,.2017,.2050,.2015,.2121,.2175,
		.2152,.1917,.1911,.1747,.1650,.1587,.1622,.1496,.1486,.1438,.1556,.1468,.1536,.1544,.1536,.1468,.1535,.1442,.1515,.1559,.1541,.1461,
		.1388,.1565,.1502,.1503,.1454,.1389,.1445,.1425,.1415,.1424,.1432,.1486,.1539,.1354,.1480,.1443,.1435,.1491,.1435,.1380,.1317,.1445,
		.1375,.1449,.1359,.1383,.1390,.1361,.1286,.1359,.1395,.1327,.1387,.1431,.1403,.1404,.1389,.1410,.1304,.1363,.1241,.1284,.1299,.1325,
		.1343,.1387,.1328,.1444,.1334,.1362,.1302,.1338,.1339,.1304,.1314,.1287,.1404,.1383,.1292,.1436,.1280,.1326,.1321,.1268,.1278,.1243,
		.1239,.1271,.1213,.1338,.1287,.1343,.1231,.1317,.1214,.1370,.1232,.1301,.1348,.1294,.1278,.1227,.1218,.1198,.1193,.1342,.1323,.1248,
		.1220,.1139,.1271,.1224,.1347,.1249,.1163,.1362,.1236,.1462,.1356,.1198,.1419,.1324,.1288,.1336,.1335,.1266};
        double xSection_n_Armstrong[160]={0.,.3125,.3930,.4401,.4582,.3774,.3329,.2996,.2715,.2165,
		.2297,.1861,.1551,.2020,.2073,.2064,.2193,.2275,.2384,.2150,.2494,.2133,.2023,.1969,.1797,.1693,.1642,.1463,.1280,.1555,.1489,.1435,
		.1398,.1573,.1479,.1493,.1417,.1403,.1258,.1354,.1394,.1420,.1364,.1325,.1455,.1326,.1397,.1286,.1260,.1314,.1378,.1353,.1264,.1471,
		.1650,.1311,.1261,.1348,.1277,.1518,.1297,.1452,.1453,.1598,.1323,.1234,.1212,.1333,.1434,.1380,.1330,.12,.12,.12,.12,.12,.12,.12,.12,
		.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
		.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
		.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12};
        for (int j = 9; j <= 70; j++) {
            temp_E_edge.push_back(temp_E_edge.back() + 25.0);
            temp_sig_Xn.push_back(0.1 * (Znu * xSection_p_Armstrong[j] + (Anu - Znu) * xSection_n_Armstrong[j]));
        }

        // 5. MichaDwell
        double eMichaDwell[11]={2000.0,3270.0,4100.0,4810.0,6210.0,6600.0,7790.0,8400.0,9510.0,13600.0,16400.0};
        double sigMichaDwell[11]={0.1266,0.1080,0.0805,0.1017,0.0942,0.0844,0.0841,0.0755,0.0827,0.0626,0.0740};
        for (int j = 0; j < 11; j++) {
            temp_E_edge.push_back(eMichaDwell[j]);
            temp_sig_Xn.push_back(0.1 * Anu * sigMichaDwell[j]);
        }

        // 6. Regge Model
        double x=0.0677, y=0.129, eps=0.0808, eta=0.4525, em=0.94, exx=pow(10,0.05);
        double s = 0.002 * em * temp_E_edge.back();
        for (int j = 1; j <= 150; j++) {
            s = s * exx;
            temp_E_edge.push_back(1000.0 * 0.5 * (s - em*em) / em);
            temp_sig_Xn.push_back(0.1 * 0.65 * Anu * (x * pow(s, eps) + y * pow(s, -eta)));
        }

        // 7. 处理 Au 专属单中子截面 1n
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
                if(sig1n_val > temp_sig_Xn[i]) sig1n_val = temp_sig_Xn[i]; // 保护
                sigma_1n.push_back(sig1n_val);
            }
        }
        delete g_1n_orig;
    }

    double GetPhotonFlux(double b_fm, double E_gamma_GeV, double Gamma_target_frame) {
        if (b_fm < 1e-3) return 0.0; // 🌟 安全保护防止中心发散
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

void calc_UPC_EMD_AuAu_Glauber_new() {
    int Z_Au = 79; int A_Au = 197;
    // Au+Au 在 200 GeV 时的光学 Glauber 参数
    double R_Au = 6.38, a_Au = 0.535, sigNN_200 = 42.0; 
    
    double sqrt_s_NN_c = 200.0; 
    double gamma_c = sqrt_s_NN_c / (2.0 * m_N);
    double Gamma_target_frame = 2.0 * gamma_c * gamma_c - 1.0; 
    
    EMDCalculator emdCalc(Z_Au, A_Au);
    GlauberCalc glauber(A_Au, R_Au, a_Au, sigNN_200);

    cout << "=========================================================" << endl;
    cout << " UPC EMD Calculator - Au+Au @ 200 GeV [HIGH-RES FINAL]" << endl;
    cout << " Woods-Saxon: R=" << R_Au << " fm, a=" << a_Au << " fm" << endl;
    cout << " Sigma_NN(200GeV) = " << sigNN_200 << " mb" << endl;
    cout << "=========================================================" << endl;

    // 🌟 核心升级：20000 个 Bin，范围 0-1000，步长 0.05 fm！
    int n_bins_b_lin = 20000;
    double b_min_lin = 0.0, b_max_lin = 1000.0; 
    double db_lin = (b_max_lin - b_min_lin) / n_bins_b_lin;

    TFile *fOut = new TFile("UPC_Probabilities_AuAu_200GeV_Glauber_Final.root", "RECREATE");
    
    // --- 独立基础概率直方图 ---
    TH1D *hProb_Surv     = new TH1D("hProb_Surv", "Probability of no hadronic collision; b (fm); P_{surv}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_1n   = new TH1D("hProb_EMD_1n", "Pure EMD SCD 1n Prob; b (fm); P_{1n}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_Xn   = new TH1D("hProb_EMD_Xn", "Pure EMD SCD Xn Prob; b (fm); P_{Xn}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_1n1n = new TH1D("hProb_EMD_1n1n", "Pure EMD MCD 1n1n Prob; b (fm); P_{1n}*P_{1n}", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_EMD_XnXn = new TH1D("hProb_EMD_XnXn", "Pure EMD MCD XnXn Prob; b (fm); P_{Xn}*P_{Xn}", n_bins_b_lin, b_min_lin, b_max_lin);

    // --- 物理乘积直方图 (兼容下游干涉代码) ---
    TH1D *hProb_MCD_1n1n = new TH1D("hProb_MCD_1n1n", "Au+Au MCD 1n1n Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_MCD_XnXn = new TH1D("hProb_MCD_XnXn", "Au+Au MCD XnXn Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_SCD_1n = new TH1D("hProb_SCD_1n", "Au+Au SCD 1n Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);
    TH1D *hProb_SCD_Xn = new TH1D("hProb_SCD_Xn", "Au+Au SCD Xn Prob; b (fm); P(b)*P_{surv}(b)", n_bins_b_lin, b_min_lin, b_max_lin);

    double total_sig_XnXn = 0, total_sig_1n1n = 0;
    double total_sig_SCD_Xn = 0, total_sig_SCD_1n = 0;

    // === 阶段 1：超高精度制表区 (0 ~ 1000 fm) ===
    for(int i=0; i<n_bins_b_lin; i++) {
        double b = b_min_lin + i*db_lin + db_lin/2.0; 
        double P_1n, P_Xn;
        emdCalc.GetProbabilities(b, Gamma_target_frame, P_1n, P_Xn);
        
        // 调用光学 Glauber 计算此 b 下无强子碰撞的生存概率
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

    // === 阶段 2：长程对数积分区 (1000 ~ 10^9 fm) ===
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

    cout << "✅ 成功: Au+Au 200GeV 概率数组(0.05 fm超高精度与成分拆分)已保存 -> UPC_Probabilities_AuAu_200GeV_Glauber_Final.root" << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "=> Mutual Coulomb Dissociation (双向 ZDC 巧合, 互激发):" << endl;
    cout << "    Sigma(1n1n MCD) = " << total_sig_1n1n * 0.01 << " barns" << endl;
    cout << "    Sigma(XnXn MCD) = " << total_sig_XnXn * 0.01 << " barns" << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "=> Single Coulomb Dissociation (单侧 ZDC 触发, 本征截面不乘2):" << endl;
    cout << "    Sigma(1n SCD)   = " << total_sig_SCD_1n * 0.01 << " barns" << endl;
    cout << "    Sigma(Xn SCD)   = " << total_sig_SCD_Xn * 0.01 << " barns" << endl;
    cout << "=========================================================" << endl;
}
