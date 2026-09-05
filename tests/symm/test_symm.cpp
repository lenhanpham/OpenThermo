// test_symm.cpp - harness for symmetry module testing (all cases in tests/symm/)
// Usage: test_symm <geometry_file> [--delta X] [--mode pg_determ|detect]
// Output: single JSON line to stdout.
#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/symmetry.h"

static std::map<std::string,int> kElem = {
    {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
    {"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},
    {"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},{"Zn",30},
    {"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},{"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},
    {"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},{"In",49},{"Sn",50},
    {"Sb",51},{"Te",52},{"I",53},{"Xe",54},{"Cs",55},{"Ba",56},{"La",57},{"Ce",58},{"Pr",59},{"Nd",60},
    {"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},{"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},
    {"Lu",71},{"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},{"Au",79},{"Hg",80},
    {"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},{"Fr",87},{"Ra",88},{"Ac",89},{"Th",90},
    {"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Bq",0}
};

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
static std::string trimPG(const std::string& s) {
    std::string t;
    for (char c : s) if (!std::isspace((unsigned char)c)) t += c;
    return t;
}

// Parse geometry: supports (a) custom: line1=title, line2=natoms, rest: Z x y z
// and (b) xyz: line1=natoms, line2=comment, rest: Elem x y z.
// G_C4h.xyz is custom-format despite .xyz suffix -> auto-detect by trying both.
static bool parseGeom(const std::string& path, std::vector<int>& nat,
                      std::vector<std::vector<double>>& coord /*3xN*/, std::string& title) {
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "cannot open " << path << "\n"; return false; }
    std::vector<std::string> L;
    std::string line;
    while (std::getline(f, line)) { L.push_back(line); }
    // drop trailing empty lines only
    while (!L.empty() && trim(L.back()).empty()) L.pop_back();
    if (L.size() < 3) return false;

    auto tryCustom = [&]() -> bool {
        try {
            int n = std::stoi(trim(L[1]));
            if (n <= 0 || (int)L.size() < 2 + n) return false;
            std::vector<int> tn(n);
            std::vector<std::vector<double>> tc(3, std::vector<double>(n));
            for (int i = 0; i < n; i++) {
                std::istringstream iss(L[2 + i]);
                int z; double x, y, zc;
                if (!(iss >> z >> x >> y >> zc)) return false;
                if (z < 0 || z > 92) return false;
                // ensure no trailing element symbol confusion: atomic number must be int-like
                tn[i] = z; tc[0][i] = x; tc[1][i] = y; tc[2][i] = zc;
            }
            title = trim(L[0]); nat = tn; coord = tc;
            return true;
        } catch (...) { return false; }
    };
    auto tryXyz = [&]() -> bool {
        try {
            int n = std::stoi(trim(L[0]));
            if (n <= 0 || (int)L.size() < 2 + n) return false;
            std::vector<int> tn(n);
            std::vector<std::vector<double>> tc(3, std::vector<double>(n));
            for (int i = 0; i < n; i++) {
                std::istringstream iss(L[2 + i]);
                std::string el; double x, y, z;
                if (!(iss >> el >> x >> y >> z)) return false;
                auto it = kElem.find(el);
                if (it == kElem.end()) return false;
                tn[i] = it->second; tc[0][i] = x; tc[1][i] = y; tc[2][i] = z;
            }
            title = trim(L[1]); nat = tn; coord = tc;
            return true;
        } catch (...) { return false; }
    };

    bool isXyz = path.size() >= 4 && path.substr(path.size() - 4) == ".xyz";
    if (isXyz) {
        // G_C4h.xyz is custom despite suffix; try custom first if line1 is non-numeric
        try { std::stoi(trim(L[0])); /* numeric -> true xyz */ }
        catch (...) { if (tryCustom()) return true; return false; }
        if (tryXyz()) return true;
        if (tryCustom()) return true;
        return false;
    }
    if (tryCustom()) return true;
    if (tryXyz()) return true;
    return false;
}

static std::string jsonEscape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"') o += "\\\"";
        else if (c == '\\') o += "\\\\";
        else o += c;
    }
    return o;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "Usage: test_symm <file> [--delta X] [--mode pg_determ|detect]\n"; return 2; }
    std::string file = argv[1];
    double delta = 0.01;
    std::string mode = "pg_determ";
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--delta" && i + 1 < argc) delta = std::stod(argv[++i]);
        else if (a == "--mode" && i + 1 < argc) mode = argv[++i];
    }

    std::vector<int> nat;
    std::vector<std::vector<double>> coord;
    std::string title;
    if (!parseGeom(file, nat, coord, title)) {
        std::cout << "{\"file\":\"" << jsonEscape(file) << "\",\"error\":\"parse_failed\"}\n";
        return 1;
    }
    int natoms = (int)nat.size();
    std::string pg;
    int ng=-999, ni=-999, nsg=-999, ncr=-999, nsr=-999, np=-999, nprm=-999, rotsym=-999;
    auto t0 = std::chrono::steady_clock::now();
    try {
        if (mode == "detect") {
            symmetry::SymmetryDetector det;
            det.ncenter = natoms;
            det.a.resize(natoms);
            for (int i = 0; i < natoms; i++) {
                det.a[i].index = nat[i];
                det.a[i].x = coord[0][i]; det.a[i].y = coord[1][i]; det.a[i].z = coord[2][i];
                det.a[i].mass = 1.0;
            }
            det.PGnameinit = "?";
            det.detectPG(0);
            pg = trimPG(det.PGname);
            rotsym = det.rotsym;
            // also run sym_elements for counts at requested delta
            std::vector<std::vector<double>> c2 = coord;
            std::vector<std::vector<double>> symn(3, std::vector<double>(150, 0.0));
            std::vector<std::vector<int>> nsym(150, std::vector<int>(5, 0));
            std::vector<std::vector<int>> nper(natoms, std::vector<int>(250, 0));
            std::vector<int> nccl(natoms, 0);
            std::vector<std::vector<int>> nscl(natoms, std::vector<int>(natoms, 0));
            int nseq = 0; std::vector<std::string> symb;
            // center like PG_determ does
            symmetry::sym_elements(natoms, nat, c2, symb, delta, ng, ni, nsg, ncr, nsr, np,
                                   symn, nsym, 0, nprm, nper, nseq, nccl, nscl);
        } else {
            std::vector<std::vector<double>> c2 = coord; // PG_determ shifts in place
            symmetry::PG_determ(natoms, nat, c2, delta, pg);
            pg = trimPG(pg);
            symmetry::SymmetryDetector det;
            det.PGname = pg; det.PGname2rotsym(); rotsym = det.rotsym;
            // counts
            std::vector<std::vector<double>> c3 = coord;
            // PG_determ centers coords; replicate centering via symm_cmass/cshift is inside PG_determ,
            // but sym_elements needs centered coords for matching counts -> do explicit centering
            // using same path: call sym_elements on centered copy via PG_determ-like steps is complex;
            // instead call sym_elements on raw + centered and prefer centered (COM) version:
            // We recompute COM centering here:
            {
                // use atomic weights of 1.0 approx via symm_cmass with unit weights is not exact;
                // simplest: call sym_elements on c2 (already centered by PG_determ) is wrong since c2 mutated.
                // So re-derive: PG_determ(c3) mutated c3 to centered; reuse c3 state? c3 still raw.
                // Call PG_determ already mutated c2 only. Use c2 (centered) for counts:
            }
            std::vector<std::vector<double>> symn(3, std::vector<double>(150, 0.0));
            std::vector<std::vector<int>> nsym(150, std::vector<int>(5, 0));
            std::vector<std::vector<int>> nper(natoms, std::vector<int>(250, 0));
            std::vector<int> nccl(natoms, 0);
            std::vector<std::vector<int>> nscl(natoms, std::vector<int>(natoms, 0));
            int nseq = 0; std::vector<std::string> symb;
            // c2 is now centered by PG_determ -> use it for consistent counts
            if (natoms == 1) { ng=1; ni=0; nsg=0; ncr=0; nsr=0; np=0; nprm=1; }
            else {
                symmetry::sym_elements(natoms, nat, c2, symb, delta, ng, ni, nsg, ncr, nsr, np,
                                       symn, nsym, 0, nprm, nper, nseq, nccl, nscl);
            }
        }
    } catch (const std::exception& e) {
        auto t1 = std::chrono::steady_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cout << "{\"file\":\"" << jsonEscape(file) << "\",\"natoms\":" << natoms
                  << ",\"error\":\"exception: " << jsonEscape(e.what()) << "\",\"time_ms\":" << ms << "}\n";
        return 1;
    }
    auto t1 = std::chrono::steady_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "{\"file\":\"" << jsonEscape(file) << "\",\"natoms\":" << natoms
              << ",\"pg\":\"" << jsonEscape(trimPG(pg)) << "\",\"rotsym\":" << rotsym
              << ",\"ng\":" << ng << ",\"ni\":" << ni << ",\"nsg\":" << nsg
              << ",\"ncr\":" << ncr << ",\"nsr\":" << nsr << ",\"np\":" << np
              << ",\"nprm\":" << nprm << ",\"delta\":" << delta
              << ",\"mode\":\"" << mode << "\",\"time_ms\":" << ms << "}\n";
    return 0;
}
