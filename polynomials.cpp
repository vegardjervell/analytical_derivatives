/*
Author : Vegard Gjeldvik Jervell
Contains : See Factorial.h for documentation and comments.
*/

#include "polynomials.h"
#include <cmath>
#include <shared_mutex>

using namespace analytical;

double binom(int n, int k) {
    if (k > n) return 0;
    if ( (k == 0) || (k == n) ) return 1;
    if ( (k == 1) || (k == n - 1) ) return n;
    return (factorial::Fac(n) / (factorial::Fac(k) * factorial::Fac(n - k))).eval();
}

void fill_partitions(const int N, int m, std::vector<std::vector<int>>& partitions){
    if (N == 0) {
        partitions.resize(1); partitions[0].clear(); return; // There is exactly one partition of zero: The empty partition.
    }
    if (N == 1){
        partitions.resize(1); partitions[0].resize(1); partitions[0][0] = 1; return;
    }
    if (m < 0) m = N;
    if (m == N){
        switch (N){
        case 2: partitions = {{2}, {1, 1}}; return;
        case 3: partitions = {{3}, {2, 1}, {1, 1, 1}}; return;
        case 4: partitions = {{4}, {3, 1}, {2, 2}, {2, 1, 1}, {1, 1, 1, 1}}; return;
        case 5: partitions = {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}; return;
        }
    }
    partitions.clear();
    std::vector<std::vector<int>> sub_partitions(N);
    // for (auto& sub_inner : sub_partitions){
    //     sub_inner.reserve(N);
    // }
    int min_k = (N - m < 0) ? 0 : N - m;
    for (int k = min_k; k < N; k++){
        int n = N - k;
        fill_partitions(k, n, sub_partitions);
        partitions.reserve(partitions.size() + sub_partitions.size());
        for (std::vector<int>& part : sub_partitions){
            part.insert(part.begin(), n);
            partitions.emplace_back(std::move(part));
        }
    }
}

void extend_partitions(int N, std::vector<std::vector<std::vector<int>>>& partitions){
    if (partitions.size() >= N + 1) return;
    int k0 = partitions.size();
    partitions.resize(N + 1);
    for (int k = k0; k <= N; k++){
        for (size_t m = k; m > 0; m--){
            for (const std::vector<int>& prev_part : partitions[k - m]){
                bool include_part = true;
                for (const int p : prev_part){
                    if (p > m) {include_part=false; break;}
                }
                if (!include_part) continue;
                partitions[k].push_back(prev_part);
                partitions[k].back().insert(partitions[k].back().begin(), m);
            }
        }
    }
}

std::vector<std::vector<int>> get_partitions(int N, int m){
    // Find all unique integer partitions of the number N, with largest value smaller than or equal to m
    // See: Memo on derivatives
    if (m < 0) m = N;
    std::vector<std::vector<int>> partitions;
    fill_partitions(N, m, partitions);
    return partitions;

    int min_k = (N - m < 0) ? 0 : N - m;
    for (int k = min_k; k < N; k++){
        int n = N - k;
        std::vector<std::vector<int>> sub_partitions = get_partitions(k, n);
        for (std::vector<int>& part : sub_partitions){
            part.insert(part.begin(), n);
            partitions.push_back(std::move(part));
        }
    }
    return partitions;
}

const std::vector<std::vector<int>>& get_partitions(int N) {
    static std::vector<std::vector<std::vector<int>>> all_partitions = {{{}}};
    static std::shared_mutex mtx;

    {
        std::shared_lock lock(mtx);
        if (all_partitions.size() > static_cast<size_t>(N)) {
            return all_partitions[N];
        }
    }

    {
        std::unique_lock lock(mtx);
        if (all_partitions.size() <= static_cast<size_t>(N)) {
            extend_partitions(N, all_partitions);
        }
        return all_partitions[N];
    }
}


std::vector<std::vector<std::vector<int>>> build_partitions(int N, int maxval){
    if (maxval < 0) maxval = N;
    std::vector<std::vector<std::vector<int>>> partitions(1, std::vector<std::vector<int>>());
    partitions.reserve(N + 1);
    partitions[0] = {{}};
    extend_partitions(N, partitions);
    return partitions;
}

inline double factorial_d(int n){
    double v = 1;
    for (; n > 1; n--) {
        v *= n;
    }
    return v;
}

double partition_multiplicity(const std::vector<int>& partition){
    // The "Multiplicity" of a partition is the number of ways a set of N elements can be subdivided into subsets of k_1, k_2, ... elements
    // For N = k_1 + k_2 + ...
    // See: Memo on derivatives
    int n = 0;
    double denom = 1;
    for (const int p : partition){
        n += p;
        denom *= factorial_d(p);
    }
    double num = factorial_d(n);

    int prev_counted = -1;
    // counted.reserve(partition.size());
    for (const int p : partition){
        if (p == prev_counted) continue;
        // bool counted = false;
        // for (const int c : counted){
        //     if (c == p) {counted = true; break;}
        // }
        // if (counted) continue;
        denom *= factorial_d(std::count(partition.begin(), partition.end(), p));
        // counted.push_back(p);
        prev_counted = p;
    }
    return num / denom; // + 0.5 To ensure correct flooring in case of funky floating point error.
}

double Analytical::numerical_derivative(double x, int n, double eps) const {
    if (n == 0) return (*this)(x);
    double v0 = derivative(x - eps, n - 1);
    double v1 = derivative(x + eps, n - 1);
    return (v1 - v0) / (2 * eps);
}

double pow::derivative(double x, int n) const {
    if (expo >= 0 && n > expo) return 0;
    double p = 1;
    double tmp = expo;
    for (int k = 0; k < n; k++){
        p *= tmp;
        tmp -= 1;
    }
    return p * std::pow(x, expo - n);
}

double log::derivative(double x, int n) const {
    if (n == 0) return std::log(x);
    if (n == 1) return 1 / x;
    return pow(-1).derivative(x, n - 1);
}

double exp::derivative(double x, int n) const {
    return std::exp(x);
}

double sin::derivative(double x, int n) const {
    const int p = (n % 4 < 2) ? 1 : -1;
    if (n % 2 == 0){
        return p * std::sin(x);
    }
    return p * std::cos(x);
}

double cos::derivative(double x, int n) const {
    const int p = ((n % 4 == 0) || (n % 4 == 3)) ? 1 : -1;
    if (n % 2 == 0){
        return p * std::cos(x);
    }
    return p * std::sin(x);
}

double Constant::derivative(double x, int n) const {
    if (n == 0) return value;
    return 0;
}

Composed::Composed(const Analytical& outer, const Analytical& inner)
    : inner{inner.make_unique()},
      outer{outer.make_unique()}
    {}

Composed::Composed(const Composed& other)
    : inner{other.inner->make_unique()},
      outer{other.outer->make_unique()}
    {}

Composed& Composed::operator=(const Composed& other){
    if (this == &other) return *this;
    inner = other.inner->make_unique();
    outer = other.outer->make_unique();
    return *this;
}

double Composed::derivative(double x, int n) const {
    const std::vector<std::vector<int>>& partitions = get_partitions(n);
    double value = 0;
    for (const std::vector<int>& part : partitions){
        double tmp = 1;
        for (const int l : part){
            tmp *= inner->derivative(x, l);
            if (tmp == 0) break;
        }
        if (tmp == 0) continue;
        double Z = partition_multiplicity(part);
        value += Z * outer->derivative((*inner)(x), part.size()) * tmp;
    }
    return value;
}

Product::Product(const Analytical& f, const Analytical& g)
    : f{f.make_unique()},
      g{g.make_unique()}
    {}

Product::Product(const Product& other)
    : f{other.f->make_unique()},
      g{other.g->make_unique()}
    {}
    
Product& Product::operator=(const Product& other){
    if (this == &other) return *this;
    f = other.f->make_unique();
    g = other.g->make_unique();
    return *this;
}

double Product::derivative(double x, int n) const {
    double val = 0;
    for (int k = 0; k <= n; k++){
        val += binom(n, k) * f->derivative(x, n - k) * g->derivative(x, k);
    }
    return val;
}

Product& Product::operator*=(const Product& other){
    g = std::make_unique<Product>(*g, other);
    return *this;
}

Sum::Sum(const Analytical& t1, const Analytical& t2){
    terms.push_back(t1.make_unique());
    terms.push_back(t2.make_unique());
}

Sum::Sum(const Sum& other) {
    terms.reserve(other.terms.size());
    for (const auto& term : other.terms) {
        terms.push_back(term->make_unique());
    }
}

Sum& Sum::operator=(const Sum& other) {
    if (this == &other) return *this;
    std::vector<std::unique_ptr<Analytical>> new_terms;
    new_terms.reserve(other.terms.size());
    for (const auto& term : other.terms) {
        new_terms.push_back(term->make_unique());
    }
    terms = std::move(new_terms);  // Strong exception safety
    return *this;
}

double Sum::derivative(double x, int n) const {
    double val = 0;
    for (const auto& term : terms){
        val += term->derivative(x, n);
    }
    return val;
}

Sum& Sum::operator+=(const Sum& other) {
    terms.reserve(terms.size() + other.terms.size());
    for (const auto& term : other.terms) {
        terms.push_back(term->make_unique());
    }
    return *this;
}

Sum& Sum::operator-=(const Sum& other){
    terms.reserve(terms.size() + other.terms.size());
    for (const auto& term : other.terms) {
        terms.push_back( (-(*term)).make_unique() );
    }
    return *this;
}

namespace analytical{
Sum operator+(const Analytical& lhs, const Analytical& rhs){
    return Sum(lhs, rhs);
}

Sum operator-(const Analytical& lhs, const Analytical& rhs){
    return Sum(lhs, -rhs);
}

Product operator*(const Analytical& lhs, const Analytical& rhs){
    return Product(lhs, rhs);
}

Product operator-(const Analytical& f){
    return Product(Constant(-1), f);
}

Product operator/(const Analytical& lhs, const Analytical& rhs){
    return Product(lhs, Composed(analytical::pow(-1), rhs));
}
} // namespace analytical

Polynomial::Polynomial(int k_min, int k_max, std::vector<double> coeff, int k_step) 
    : k_min{k_min}, k_max{k_max}, k_step{k_step},
    _is_linear{((k_min == 0) || (k_min == 1)) && (k_max == 1)}, 
    _is_constant{(k_min == 0) && (k_max == 0)}, 
    coeff(coeff)
    {
        if ((k_max - k_min) % k_step != 0) {
            throw std::out_of_range("Polynomial: k_min - k_max is not a multiple of k_step : (" + std::to_string(k_min) + ", " + std::to_string(k_max) + ", " + std::to_string(k_step) + ")");
        }
        if (((k_max - k_min) / k_step) + 1 != coeff.size()) {
            throw std::out_of_range("Polynomial: Wrong number of coefficients (Expected " + std::to_string(((k_max - k_min) / k_step) + 1) + ", got " + std::to_string(coeff.size())
                                    + ")\n(min, max, step) = (" + std::to_string(k_min) + ", " + std::to_string(k_max) + ", " + std::to_string(k_step) + ")");
        }
    }

double Polynomial::operator()(double x) const {
    if (_is_constant) return coeff[0];
    double p = 0;
    size_t C_idx = 0;
    for (int k = k_min; k <= k_max; k += k_step){
        p += coeff[C_idx++] * std::pow(x, k);
    }
    return p;
}

double Polynomial::derivative(double x, int n) const {
    if (n == 0) return this->operator()(x);
    if (_is_constant) return 0;
    if (_is_linear && (n > 1)) return 0;
    if (k_min >= 0 && n > k_max) return 0;

    double p = 0;
    int prefactor = ((n % 2 == 0) ? 1 : -1);
    int max_k_negative = (k_max < 0) ? k_max : -1;
    size_t C_idx = 0;
    int k = k_min;
    double px = std::pow(x, k - n);
    double px_step = std::pow(x, k_step);
    for (; k <= max_k_negative; k += k_step){
        p += prefactor * factorial::partialfactorial(- k - 1, n - k - 1) * coeff[C_idx++] * px; // pow(x, k - n);
        px *= px_step;
    }
    for (; k < n; k += k_step) C_idx++;
    px = std::pow(x, k - n);
    for (; k <= k_max; k += k_step){
        // if (k < n){C_idx++; continue;}
        p += factorial::partialfactorial(k - n, k) * coeff[C_idx++] * px; // pow(x, k - n);
        px *= px_step; 
    }
    return p;
}

PolyExp::PolyExp(Polynomial pref, Polynomial expo)
    : pref{pref}, expo{expo}
    {}

PolyExp::PolyExp(Polynomial pref, double expo)
    : pref{pref}, expo{Polynomial::constant(expo)}
    {}

PolyExp::PolyExp(double pref, Polynomial expo)
    : pref{Polynomial::constant(pref)}, expo{expo}
    {}

double PolyExp::operator()(double x) const {
    return pref(x) * std::exp(expo(x));
}

double PolyExp::derivative(double x, int n) const {
    if (expo._is_constant){
        if (expo.coeff[0] == 0) return pref.derivative(x, n);
        return pref.derivative(x, n) * std::exp(expo.coeff[0]);
    }
    if (n == 0) {
        double p = pref(x);
        double E = std::exp(expo(x));
        return p * E;
    }

    vector1d df(n + 1);
    vector1d dg(n + 1);
    for (int k = 0; k <= n; k++){
        df[k] = pref.derivative(x, k);
        dg[k] = expo.derivative(x, k);
    }

    double val = 0;
    double binom_val = 1;
    int k_start = ((pref.k_min >= 0) && (pref.k_max < n)) ? n - pref.k_max : 0;
    int k = 0;
    for (; k < k_start; k++){
        binom_val *= static_cast<double>(n - k) / (k + 1);
    }
    for (; k <= n; k++){
        if (df[n - k] == 0){
            binom_val *= static_cast<double>(n - k) / (k + 1);
            continue;
        }
        double Gk = get_Gk(x, k, dg);
        val += binom_val * Gk * df[n - k];
        binom_val *= static_cast<double>(n - k) / (k + 1);
    }
    return val * std::exp(dg[0]);
}

double PolyExp::get_Gk(double x, int k, const vector1d& dg, const std::vector<std::vector<int>>& partitions) const {
    if (k == 0) return 1;
    if (expo._is_linear) return std::pow(expo.derivative(x, 1), k);
    int max_dg_order = ((expo.k_min >= 0) && (expo.k_max < k)) ? expo.k_max : k;
    double G = 0;
    for (const std::vector<int>& partition : partitions){
        double tmp = 1;
        for (int l : partition){
            if (l > max_dg_order) break;
            tmp *= dg[l];
        }
        if (tmp != 0){
            double Z = partition_multiplicity(partition);
            G += Z * tmp;
        }
    }
    return G;
}

double PolyExp::get_Gk(double x, int k, const vector1d& dg) const {
    const std::vector<std::vector<int>>& partitions = get_partitions(k);
    // int max_dg_order = ((expo.k_min >= 0) && (expo.k_max < k)) ? expo.k_max : k;
    // const std::vector<std::vector<int>> partitions = get_partitions(k, max_dg_order);
    return get_Gk(x, k, dg, partitions);
}

PolyFrac::PolyFrac(Polynomial num, Polynomial denom)
    : num{num}, denom{denom}
    {}

double PolyFrac::operator()(double x) const {
    return num(x) / denom(x);
}

double PolyFrac::derivative(double x, int n) const {
    vector1d df(n + 1);
    vector1d dg(n + 1);
    for (int k = 0; k <= n; k++){
        df[k] = num.derivative(x, k);
        dg[k] = denom.derivative(x, k);
    }

    double val = 0;
    double binom_val = 1;
    int k_start = ((num.k_min >= 0) && (num.k_max < n)) ? n - num.k_max : 0;
    int k = 0;
    for (; k < k_start; k++){
        binom_val *= static_cast<double>(n - k) / (k + 1);
    }

    for (; k <= n; k++){
        if (df[n - k] == 0){
            binom_val *= static_cast<double>(n - k) / (k + 1);
            continue;
        }
        double Hk = get_Hk(x, k, dg);
        val += binom_val * Hk * df[n - k];
        binom_val *= static_cast<double>(n - k) / (k + 1);
    }
    return val / dg[0];

}

double PolyFrac::get_Hk(double x, int k, const vector1d& dg) const {
    const std::vector<std::vector<int>>& partitions = get_partitions(k);
    int max_dg_order = ((denom.k_min >= 0) && (denom.k_max < k)) ? denom.k_max : k;
    double H = 0;
    for (const std::vector<int>& partition : partitions){
        double tmp = 1;
        for (int l : partition){
            if (l > max_dg_order) break;
            tmp *= dg[l];
        }
        if (tmp != 0){
            double Z = partition_multiplicity(partition);
            H += Z * tmp * factorial_d(static_cast<double>(partition.size())) * std::pow(-1. / dg[0], partition.size());
        }
    }
    return H;
}

template<typename T>
std::unique_ptr<Analytical> __make_unique(const T& obj) {
    return std::make_unique<T>(obj);
}

#define DEFINE_MAKE_UNIQUE(func) std::unique_ptr<Analytical> func::make_unique() const {return __make_unique(*this);}
DEFINE_MAKE_UNIQUE(pow)
DEFINE_MAKE_UNIQUE(log)
DEFINE_MAKE_UNIQUE(exp)
DEFINE_MAKE_UNIQUE(sin)
DEFINE_MAKE_UNIQUE(cos)
DEFINE_MAKE_UNIQUE(Constant)
DEFINE_MAKE_UNIQUE(Composed)
DEFINE_MAKE_UNIQUE(Product)
DEFINE_MAKE_UNIQUE(Sum)
DEFINE_MAKE_UNIQUE(Polynomial)
DEFINE_MAKE_UNIQUE(PolyExp)
DEFINE_MAKE_UNIQUE(PolyFrac)
#undef DEFINE_MAKE_UNIQUE

namespace analytical {

std::ostream& operator<<(std::ostream& strm, const Polynomial& p){
    size_t C_idx = 0;
    for (int k = p.k_min; k < p.k_max; k += p.k_step){
        strm << "(" << p.coeff[C_idx++] << ", " << k << ") + ";
    }
    strm << "(" << p.coeff[C_idx++] << ", " << p.k_max << ")";
    return strm;
}

std::ostream& operator<<(std::ostream& strm, const PolyExp& p){
    strm << "[ " << p.pref << " ] * exp[ " << p.expo << " ]";
    return strm;
}

std::ostream& operator<<(std::ostream& strm, const PolyFrac& p){
    strm << "( " << p.num << " ) / ( " << p.denom << " )";
    return strm;
}

}