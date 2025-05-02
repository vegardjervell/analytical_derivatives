#include "polynomials.h"

using namespace analytical; // NOTE: Contains "sin", "cos", "log" and "exp". Be careful when combining with <cmath>.

int main(){

{
    std::cout << "------------ PARTITION DEMO START ------------" << std::endl;
    for (int n = 1; n < 6; n++){
        const std::vector<std::vector<int>>& partitions = get_partitions(n);
        std::cout << "n = " << n << std::endl;
        for (const auto& part : partitions){
            std::cout << "\t" << partition_multiplicity(part) << " : {";
            for (int l : part){
                std::cout << l << ", ";
            }
            std::cout << "}" << std::endl;
        }
    }
    std::cout << "------------ PARTITION DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ BASE FUNCTION DEMO START ------------" << std::endl;
    double x = 1.3456;
    sin s; // sin(x)
    cos c; // cos(x)
    log l; // ln(x)
    exp e; // exp(x)
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << s.derivative(x, n) << " / " << s.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << c.derivative(x, n) << " / " << c.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << l.derivative(x, n) << " / " << l.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << e.derivative(x, n) << " / " << e.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ BASE FUNCTION DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ COMPOSED DEMO START ------------" << std::endl;
    Sum func1 = 0 + pow(6); // 0 + x^6
    Composed func2(pow(2), pow(3)); // (x^3)^2
    double x = 1.3925;
    std::cout << func1(x) << std::endl;
    std::cout << func2(x) << std::endl;
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << func1.derivative(x, n) << " / " << func1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << func2.derivative(x, n) << " / " << func2.numerical_derivative(x, n, 1e-8) << std::endl;
    }

    Composed func3{pow(4), sin()}; // sin^4(x)
    Composed func4{exp(), sin()}; // exp(sin(x))
    Composed func5{exp(), pow(-2)}; // exp(x^{-2})
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << func3.derivative(x, n) << " / " << func3.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << func4.derivative(x, n) << " / " << func4.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << func5.derivative(x, n) << " / " << func5.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ COMPOSED DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ PRODUCT DEMO START ------------" << std::endl;
    double x = 1.5342;
    Sum func1 = 2 * pow(-2) + 3 * pow(-1) + 4 - 2 * pow(1) - 3 * pow(2);
    Sum func2 = -2 * pow(-3) + 1.5 * pow(-1) + 1 + 3 * pow(2) + 4 * pow(3);
    Product func3 = func1 / func2;
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << func1.derivative(x, n) << " / " << func1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << func2.derivative(x, n) << " / " << func2.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << func3.derivative(x, n) << " / " << func3.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ PRODUCT DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ Polynomial DEMO START ------------" << std::endl;
    double x = 2.3578;

    Polynomial p1(-2, 2, {1, 2, 0.5, 2.3, -3});
    Sum s1 = 1 * pow(-2) + 2 * pow(-1) + 0.5 * pow(0) + 2.3 * pow(1) - 3 * pow(2);

    Polynomial p2(-8, -4, {-1, 2, -3}, 2);
    Sum s2 = -1 * pow(-8) + 2 * pow(-6) -3 * pow(-4);

    Product frac1 = p1 / p2;
    Product frac2 = s1 / s2;
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << p1.derivative(x, n) << " / " << p1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << s1.derivative(x, n) << " / " << s1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << p2.derivative(x, n) << " / " << p2.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << s2.derivative(x, n) << " / " << s2.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << frac1.derivative(x, n) << " / " << frac1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << frac2.derivative(x, n) << " / " << frac2.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ Polynomial DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ PolyExp DEMO START ------------" << std::endl;
    double x = 5.27381;
    Polynomial p1(-2, 2, {1, 2, 0.5, 2.3, -3});
    Polynomial p2(-8, -4, {-1, 2, -3}, 2);
    PolyExp pe1(p1, p2); // p1(x) * exp(p2(x))
    Product pe2 = p1 * Composed(exp(), p2); // p1(x) * exp(p2(x))
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << pe1.derivative(x, n) << " / " << pe1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << pe2.derivative(x, n) << " / " << pe2.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ PolyExp DEMO END ------------\n" << std::endl;
}{
    std::cout << "------------ PolyFrac DEMO START ------------" << std::endl;
    double x = 2.1457;
    Polynomial p1(0, 2, {1.0, 2.0, 3.0});
    Polynomial p2(-2, 2, {-1.0, 2.0, -3.0}, 2);
    PolyFrac pf1(p1, p2); // p1(x) / p2(x)
    Product pf2 = p1 / p2;
    for (int n = 0; n < 8; n++){
        std::cout << "n = " << n << std::endl;
        std::cout << "\t" << pf1.derivative(x, n) << " / " << pf1.numerical_derivative(x, n, 1e-8) << std::endl;
        std::cout << "\t" << pf2.derivative(x, n) << " / " << pf2.numerical_derivative(x, n, 1e-8) << std::endl;
    }
    std::cout << "------------ PolyFrac DEMO END ------------\n" << std::endl;
}
    return 0;
}