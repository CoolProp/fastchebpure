// C++ standard library
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <optional>
#include <numeric>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <list>
#include <vector>
#include <Eigen/Dense>

#include "teqp/models/multifluid.hpp"
#include "teqp/models/multifluid_ancillaries.hpp"
#include "teqp/algorithms/VLE_pure.hpp"
#include "teqp/algorithms/critical_pure.hpp"
#include "teqp/derivs.hpp"

// Imports from boost
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "boost/functional/hash.hpp"
using namespace boost::multiprecision;
using my_float_mp = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100>>;

#include "ChebTools/ChebTools.h"

extern const std::filesystem::path teqp_datapath{ "../externals/CoolProp" };
extern const std::filesystem::path output_prefix{"../output/"};
extern const std::filesystem::path check_destination{"../outputcheck/"};

using namespace ChebTools;

/**
* @brief This class stores sets of L matrices (because they are a function only of the degree of the expansion)
*
* The L matrix is used to convert from functional values to coefficients, as in \f[ \vec{c} = \mathbf{L}\vec{f} \f]
*/
class LMatrixLibrary {
private:
    std::map<std::size_t, Eigen::MatrixXd> matrices;
    void build(std::size_t N) {
        Eigen::MatrixXd L(N + 1, N + 1); ///< Matrix of coefficients
        for (int j = 0; j <= N; ++j) {
            for (int k = j; k <= N; ++k) {
                double p_j = (j == 0 || j == N) ? 2 : 1;
                double p_k = (k == 0 || k == N) ? 2 : 1;
                L(j, k) = 2.0 / (p_j*p_k*N)*cos((j*EIGEN_PI*k) / N);
                // Exploit symmetry to fill in the symmetric elements in the matrix
                L(k, j) = L(j, k);
            }
        }
        matrices[N] = L;
    }
public:
    /// Get the \f$\mathbf{L}\f$ matrix of degree N
    const Eigen::MatrixXd & get(std::size_t N) {
        auto it = matrices.find(N);
        if (it != matrices.end()) {
            return it->second;
        }
        else {
            build(N);
            return matrices.find(N)->second;
        }
    }
};
static LMatrixLibrary l_matrix_library;

using VectorDyadicSplittingFunction = std::function<std::vector<double>(double)>;

auto vector_factory(const std::size_t N, const VectorDyadicSplittingFunction& func, const double xmin, const double xmax) -> std::vector<ChebyshevExpansion>{
    
    // Get the precalculated Chebyshev-Lobatto nodes
    const Eigen::VectorXd & x_nodes_n11 = get_CLnodes(N);

    // Step 1&2: Grid points functional values (function evaluated at the
    // extrema of the Chebyshev polynomial of order N - there are N+1 of them)
    Eigen::VectorXd fL(N + 1), fR(N+1), fp(N+1);
    for (int k = 0; k <= N; ++k) {
        // The extrema in [-1,1] scaled to real-world coordinates
        double x_k = ((xmax - xmin)*x_nodes_n11(k) + (xmax + xmin)) / 2.0;
        auto funcvals = func(x_k);
        fL(k) = funcvals[0];
        fR(k) = funcvals[1];
        fp(k) = funcvals[2];
    }

    // Step 3: Get coefficients for the L matrix from the library of coefficients
    const Eigen::MatrixXd &L = l_matrix_library.get(N);
    // Step 4: Obtain coefficients from vector - matrix product
    return {ChebyshevExpansion(L*fL, xmin, xmax), ChebyshevExpansion(L*fR, xmin, xmax), ChebyshevExpansion(L*fp, xmin, xmax)};
}

// A convenience function that returns true if the paired expansions have converged
// according to the convergence criterion
template<typename Container = std::vector<std::vector<ChebTools::ChebyshevExpansion>>>
bool are_converged(int Msplit, double tol, const Container& ce, std::size_t i){
    // Convenience function to get the M-element norm ratio, which is our convergence criterion
    auto get_err = [Msplit](const ChebTools::ChebyshevExpansion& ce) { return ce.coef().tail(Msplit).norm() / ce.coef().head(Msplit).norm(); };
    
    for (auto iexpansion = 0; iexpansion < ce.size(); ++iexpansion){
        if (get_err(ce[iexpansion][i]) > tol){
            return false;
        }
    }
    return true;
};

using Container = std::vector<std::vector<ChebyshevExpansion>>;
using VectorDyadicSplittingCallback = std::function<void(int, const Container&)>;

template<typename Container = std::vector<std::vector<ChebyshevExpansion>>>
auto vectored_dyadic_splitting(const std::size_t N, const VectorDyadicSplittingFunction& func, const double xmin, const double xmax,
    const int M, const double tol, const int max_refine_passes = 8,
    const VectorDyadicSplittingCallback& callback = nullptr) -> Container
{
    
    // Start off with the full domain from xmin to xmax
    auto vector_expansions = vector_factory(N, func, xmin, xmax);
    Container expansions(vector_expansions.size());
    for (auto i = 0; i < vector_expansions.size(); ++i){
        expansions[i].emplace_back(vector_expansions[i]);
    }

    // Now enter into refinement passes
    for (int refine_pass = 0; refine_pass < max_refine_passes; ++refine_pass) {
        bool all_converged = true;
        // Start at the right and move left because insertions will make the length increase
        for (int iexpansion = static_cast<int>(expansions[0].size())-1; iexpansion >= 0; --iexpansion) {
            if (!are_converged(M, tol, expansions, iexpansion)) {
                double xmin = expansions[0][iexpansion].xmin();
                double xmax = expansions[0][iexpansion].xmax();
                // Splitting is required, do a dyadic split
                auto xmid = xmin*0.25 + xmax*0.75;
                std::cout << "s" << std::endl;
                auto newleft = vector_factory(N, func, xmin, xmid);
                auto newright = vector_factory(N, func, xmid, xmax);
                using ArrayType = decltype(newleft[0].coef());

                // Function to check if any coefficients are invalid (evidence of a bad function value)
                auto all_coeffs_ok = [](const ArrayType& v) {
                    for (auto i = 0; i < v.size(); ++i) {
                        if (!std::isfinite(v[i])) { return false; }
                    }
                    return true;
                };
                // Check if any coefficients are invalid, stop if so
                for (auto& expansion: newleft){
                    if (!all_coeffs_ok(expansion.coef())) {
                        throw std::invalid_argument("At least one coefficient is non-finite in new left");
                    }
                }
                for (auto& expansion: newright){
                    if (!all_coeffs_ok(expansion.coef())) {
                        throw std::invalid_argument("At least one coefficient is non-finite in new right");
                    }
                }
                
                for (auto iancillary = 0; iancillary < expansions.size(); ++iancillary){
                    std::swap(expansions[iancillary][iexpansion], newleft[iancillary]);
                    expansions[iancillary].insert(expansions[iancillary].begin() + iexpansion+1, newright[iancillary]);
                }
                
                all_converged = false;
            }
        }
        if (callback != nullptr) {
//            const PairedDyadicSplittingCallback func = callback.value();
            callback(refine_pass, expansions);
        }
        
        if (all_converged) { break; }
    }
    return expansions;
};

/// A custom exception class that holds onto the temperature
struct FailedIteration : public std::exception {
    std::string msg;
    double T;
    FailedIteration(double T, const std::string& msg) : T(T), msg(msg) {};
    const char* what() const noexcept override {
        return msg.c_str();
    }
};

/**
 With a range of starting points covering a range of points around the ceritical point given
 by the EOS developers, search for the numerical critical critical point satisfying
 dp/drho|T = 0 and d2p/drho^2|T = 0. There are in some cases multiple solutions to
 these constraints, so apply additional screening to rule out numerical wiggles that
 satisfy the derivative constraints but should be rejected
 */
template<class Model>
auto aggressively_solve_pure_critical(const Model& model, double Tcrit0, double rhocrit0, double fracerr_rho_tol=1e-7){
    
    using tdx = teqp::TDXDerivatives<Model, double, Eigen::ArrayXd>;
    auto molefrac = (Eigen::ArrayXd(1) << 1.0).finished();
    auto R = model.R(molefrac);
    auto Trange = Eigen::ArrayXd::LinSpaced(50, 0.9*Tcrit0, 2*Tcrit0);
    auto rhorange = Eigen::ArrayXd::LinSpaced(50, 0.5*rhocrit0, 2*rhocrit0);
    
    std::vector<double>Tsolns, rhosolns;
    for (auto Tstart: Trange){
        for (auto rhostart: rhorange){
            double Tsoln, rhosoln;
            std::tie(Tsoln, rhosoln) = solve_pure_critical(model, Tstart, rhostart);
            
            if ((Tsoln < Trange.minCoeff()) || (Tsoln > Trange.maxCoeff())){ continue;  }
            if ((rhosoln < rhorange.minCoeff()) || (rhosoln > rhorange.maxCoeff())){ continue;  }
            if (!std::isfinite(Tsoln) || !std::isfinite(rhosoln)){ continue; }
            
            auto get_derivs = [&](double T, double rho){
                auto Ar0n = tdx::template get_Ar0n<3>(model, T, rho, molefrac);
                double Ar00 = Ar0n[0], Ar01 = Ar0n[1], Ar02 = Ar0n[2], Ar03 = Ar0n[3];
                
                double dpdrho = R*T*(1 + 2*Ar01 + Ar02);
                double d2pdrho2 = R*T/rho*(2*Ar01 + 4*Ar02 + Ar03);
                return std::make_tuple(dpdrho, d2pdrho2);
            };
            auto [dpdrho, d2pdrho2] = get_derivs(Tsoln, rhosoln);
            if (std::abs(dpdrho) > 1e-9){ continue; }
            if (std::abs(d2pdrho2) > 1e-9){ continue; }
            
            double drho = 1e-2*rhosoln;
            // dp/drho|T should be positive to the left and right of solution along the isotherm
            // and d2p/drho2|T should change sign on either side of the solution
            auto [dpdrhoR, d2pdrho2R] = get_derivs(Tsoln, rhosoln+drho);
            auto [dpdrhoL, d2pdrho2L] = get_derivs(Tsoln, rhosoln-drho);
            if (dpdrhoL < 0 || dpdrhoR < 0){ continue; }
            if (d2pdrho2L*d2pdrho2R > 0) { continue; }
            
//            std::cout << Tsoln << "," << rhosoln << std::endl;
            Tsolns.push_back(Tsoln);
            rhosolns.push_back(rhosoln);
        }
    }
    Eigen::ArrayXd Tsolns_ = Eigen::Map<Eigen::ArrayXd>(&Tsolns[0], Tsolns.size());
    Eigen::ArrayXd rhosolns_ = Eigen::Map<Eigen::ArrayXd>(&rhosolns[0], rhosolns.size());
    auto stddev = [](Eigen::ArrayXd& vec){ return std::sqrt((vec - vec.mean()).square().sum()/(vec.size()-1)); };
    
    double fracerr_T = stddev(Tsolns_)/Tsolns_.mean();
    double fracerr_rho = stddev(rhosolns_)/rhosolns_.mean();
//    std::cout << stddev(Tsolns_) << "," << Tsolns_.mean() << std::endl;
//    std::cout << stddev(rhosolns_) << "," << rhosolns_.mean() << std::endl;
    
    if (fracerr_T < 1e-10 && fracerr_rho < fracerr_rho_tol){
        return std::make_tuple(Tsolns_.mean(), rhosolns_.mean());
    }
    else{
        throw std::invalid_argument("Could not solve for critical point. There were "+std::to_string(Tsolns_.size())+" solutions and relative log10(std(y)/mean(y)) for T and rho were: "+std::to_string(log10(fracerr_T)) + "," + std::to_string(log10(fracerr_rho)));
    }
}

/** The function to  fit a superancillary function for a given fluid
 
 \note It is thread-safe, so can be run in parallel
 */
void build_superancillaries(const std::string &fluid, const std::filesystem::path &ofpath){
    const std::filesystem::path fluid_json_path = teqp_datapath / "dev" / "fluids" / (fluid + ".json");
    auto model = teqp::build_multifluid_model({ fluid_json_path.string()}, teqp_datapath.string());

    // Build conventional ancillaries
    auto build_ancillaries = [](const auto& c, double Tctrue, double rhoctrue) {
        if (c.redfunc.Tc.size() != 1) {
            throw teqp::InvalidArgument("Can only build ancillaries for pure fluids");
        }
        auto jancillaries = nlohmann::json::parse(c.get_meta()).at("pures")[0].at("ANCILLARIES");
        // Hack the ancillaries to have the true critical point as their reducing point
        jancillaries["rhoL"]["T_r"] = Tctrue;
        jancillaries["rhoL"]["Tmax"] = Tctrue;
        jancillaries["rhoL"]["reducing_value"] = rhoctrue;
        jancillaries["rhoV"]["T_r"] = Tctrue;
        jancillaries["rhoV"]["Tmax"] = Tctrue;
        jancillaries["rhoV"]["reducing_value"] = rhoctrue;
        return teqp::MultiFluidVLEAncillaries(jancillaries);
    };
    
    // Convenience function to get the density derivatives
    auto getdrhodTs = [](const auto& model, double T, double rhoL, double rhoV){
        auto molefrac = (Eigen::ArrayXd(1) << 1.0).finished();
        double R = model.R(molefrac);
        double dpsatdT = dpsatdT_pure(model, T, rhoL, rhoV);
        using tdx = teqp::TDXDerivatives<decltype(model)>;

        auto get_drhodT = [&](double T, double rho){
            double dpdrho = R*T*(1 + 2*tdx::get_Ar01(model, T, rho, molefrac) + tdx::get_Ar02(model, T, rho, molefrac));
            double dpdT = R*rho*(1 + tdx::get_Ar01(model, T, rho, molefrac) - tdx::get_Ar11(model, T, rho, molefrac));
            return -dpdT/dpdrho + dpsatdT/dpdrho;
        };

        return std::make_tuple(get_drhodT(T, rhoL), get_drhodT(T, rhoV));
    };
    
    // The calculation cache for results of calculations
    struct DensitiesType { const my_float_mp rhoL, rhoV, pL, pV, DeltarhocritL, DeltarhocritV; };
    std::unordered_map<double, DensitiesType, boost::hash<double>> densitydb;

    auto j = teqp::load_a_JSON_file(fluid_json_path.string());
    double Tcrit = j.at("STATES").at("critical").at("T"); // Critical temperature, according to the EOS developers
    double rhomolarcrit = j.at("STATES").at("critical").at("rhomolar"); // Critical density, according to the EOS developers
    double Treducing = j.at("EOS")[0].at("STATES").at("reducing").at("T"); // Reducing temperature
//    double rhomolar_reducing = j.at("EOS")[0].at("STATES").at("reducing").at("rhomolar"); // Reducing molar density
    double Ttriple = j.at("STATES").at("triple_liquid").at("T"); // Triple-point (SLV) temperature
    double Tmin_sat = j.at("EOS")[0].at("STATES").at("sat_min_liquid").at("T"); // Minimum saturation temperature
    bool pseudo_pure = j.at("EOS")[0].at("pseudo_pure");
    double R = j.at("EOS")[0].at("gas_constant"); // Gas constant being used

    double Tcrittrue, rhocrittrue;
    if (fluid == "NITROGEN" || fluid == "Nitrogen"){
        std::tie(Tcrittrue, rhocrittrue) = solve_pure_critical(model, Tcrit, rhomolarcrit);
    }
    else if (fluid == "DMC" || fluid == "MXYLENE" || fluid == "DimethylCarbonate" || fluid == "m-Xylene"){
        // Need to relax the tolerance a bit because these ones have multiple critical points
        // The temperature is good, but the densities are not. Maybe this is catastrophic truncation?
        std::tie(Tcrittrue, rhocrittrue) = aggressively_solve_pure_critical(model, Tcrit, rhomolarcrit, 1e-5);
    }
    else{
        std::tie(Tcrittrue, rhocrittrue) = aggressively_solve_pure_critical(model, Tcrit, rhomolarcrit);
    }
    
    using tdxflt = teqp::TDXDerivatives<decltype(model), my_float_mp, Eigen::ArrayX<my_float_mp>>;
    const auto z = (Eigen::ArrayX<my_float_mp>(1) << 1.0).finished();
    auto pcrittrue = rhocrittrue*R*Tcrittrue*(1.0 + tdxflt::get_Ar01<teqp::ADBackends::multicomplex>(model, Tcrittrue, rhocrittrue, z));
    
    auto anc = build_ancillaries(model, Tcrittrue, rhocrittrue);
    
    // Build critical region polynomial for each phase
    // to be used in place of conventional ancillary equation
    Eigen::ArrayXd cLarray, cVarray;
    double critical_polynomial_Theta = 0.01;
    {
        std::vector<double> Thetas, rhoLs, rhoVs;
        double T = Tcrittrue*(1-critical_polynomial_Theta), dT = 0.0, rhoL=anc.rhoL(T), rhoV=anc.rhoV(T);
        double numsteps = 1000, acceleration = 0.5;
        for (auto counter = 0; counter < numsteps; ++counter){
            // Set up the residual function
            teqp::IsothermPureVLEResiduals<decltype(model), my_float_mp, teqp::ADBackends::multicomplex> residual(model, T);
            auto rhovec = teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, rhoL, rhoV, 10).cast<double>();
            rhoL = rhovec[0]; rhoV = rhovec[1];
            auto [drhodTL, drhodTV] = getdrhodTs(model, T, rhoL, rhoV);
            rhoL += drhodTL*dT;
            rhoV += drhodTV*dT;
            T += dT;
            auto Tnew = acceleration*Tcrittrue + (1-acceleration)*T;
            dT = Tnew-T;
            auto Theta = (Tcrittrue-T)/Tcrittrue;
            if (Theta < 1e-10){
                break;
            }
//            std::cout << Theta << "," << rhoL << "," << rhoV << std::endl;
            if (!std::isfinite(rhoL) || !std::isfinite(rhoV) || !std::isfinite(Theta) || rhoL <= rhoV){
                break;
            }
            
            Thetas.push_back(Theta);
            rhoLs.push_back(rhoL);
            rhoVs.push_back(rhoV);
        }
        
        // Solve the least-squares problem for the polynomial coefficients
        auto N = Thetas.size();
        Eigen::MatrixXd A(N,7);
        Eigen::VectorXd bL(N), bV(N);
        for (auto i = 0; i < 7; ++i){
            auto view = (Eigen::Map<Eigen::ArrayXd>(&(Thetas[0]), N)).log();
            A.col(i) = view.pow(i);
        }
        bL = (Eigen::Map<Eigen::ArrayXd>(&(rhoLs[0]), N)-rhocrittrue).log();
        bV = (rhocrittrue-Eigen::Map<Eigen::ArrayXd>(&(rhoVs[0]), N)).log();
        cLarray = A.colPivHouseholderQr().solve(bL).array();
        cVarray = A.colPivHouseholderQr().solve(bV).array();
//        std::cout << cLarray << std::endl;
//        std::cout << cVarray << std::endl;
        
        /////// Code to check the results of the fit
        auto rhoLcheck = ((A*cLarray.matrix()).array()).exp() + rhocrittrue;
        auto rhoVcheck = rhocrittrue - ((A*cVarray.matrix()).array()).exp();
        auto rhoLdev = rhoLcheck/Eigen::Map<Eigen::ArrayXd>(&(rhoLs[0]), N)-1;
        auto rhoVdev = rhoVcheck/Eigen::Map<Eigen::ArrayXd>(&(rhoVs[0]), N)-1;
        std::cout << rhoLdev.abs().mean()*100 << "% (liq) for critical ancillary for " << fluid << std::endl;
        std::cout << rhoVdev.abs().mean()*100 << "% (vap) for critical ancillary for " << fluid << std::endl;
    }
    
    // Get the Brho values for each phase at a specified value of T far enough
    // from the critical point to allow VLE calcs to succeed, but close enough to get
    // the right scaling behavior
    auto calc_Brhos = [&](double T){
        auto Theta = (Tcrittrue-T)/Tcrittrue;
        auto rhoLrhoV = pure_VLE_T(model, T, anc.rhoL(T), anc.rhoV(T), 10);
        return std::make_tuple((rhoLrhoV[0]-rhocrittrue)/sqrt(Theta), (rhoLrhoV[1]-rhocrittrue)/sqrt(Theta));
    };
    double BrhoL, BrhoV;
    std::tie(BrhoL, BrhoV) = calc_Brhos(Tcrittrue*(1-0.001));

//    if (pseudo_pure) {
//        return std::make_tuple(std::vector<ChebTools::ChebyshevExpansion>{}, nlohmann::json{});
//    }
    
    struct CriticalEstimation {
        double Brho, beta, Tcrittrue, rhocrittrue;
        double operator ()(double T){ return rhocrittrue + Brho*pow((Tcrittrue-T)/Tcrittrue, beta); }
    };
    std::optional<CriticalEstimation> last_estimationL, last_estimationV;
    
    // A function to get the co-existing densities and pressure for a given value of temperature
    VectorDyadicSplittingFunction get_VLE_values = [&](double T) -> std::vector<double>{
        if (std::abs(T / Tcrittrue - 1) < 1e-14) {
            return {rhocrittrue, rhocrittrue, static_cast<double>(pcrittrue)};
        }
        else if (densitydb.count(T) == 0) { // If not in cache...
            // Do the calculation and store in cache
            
            auto Theta = (Tcrittrue-T)/Tcrittrue;
            
            // Set up the residual function
            teqp::IsothermPureVLEResiduals<decltype(model), my_float_mp, teqp::ADBackends::multicomplex> residual(model, T);
            decltype(teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, 1.0, 1.0, 10)) rhovec;
            
            auto x = log(Theta);
            auto rhoLpoly = rhocrittrue + exp(cLarray[0] + cLarray[1]*x + cLarray[2]*pow(x, 2) + cLarray[3]*pow(x, 3) + cLarray[4]*pow(x, 4) + cLarray[5]*pow(x, 5) + cLarray[6]*pow(x, 6));
            auto rhoVpoly = rhocrittrue - exp(cVarray[0] + cVarray[1]*x + cVarray[2]*pow(x, 2) + cVarray[3]*pow(x, 3) + cVarray[4]*pow(x, 4) + cVarray[5]*pow(x, 5) + cVarray[6]*pow(x, 6));
            
            // Try to just do the iteration, let's hope this will work
            if (Theta < critical_polynomial_Theta){
                // Use a polynomial fit to the density in the critical region as starting densities
                rhovec = teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, rhoLpoly, rhoVpoly, 10);
            }
            else{
                // Use the conventional ancillary functions as the starting densities
                rhovec = teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, anc.rhoL(T), anc.rhoV(T), 10);
            }
            
            bool bad_solution = false;
            if (!std::isfinite(static_cast<double>(rhovec[0])) || rhovec[1] >= rhovec[0] || rhovec[0] < 0){
                bad_solution = true;
            }
            
            // Now we see if an error has occurred
            if (bad_solution) {
//                double rhoLanc = anc.rhoL(T), rhoVanc = anc.rhoV(T);
//                std::cout << rhoLanc << "," << rhoLextrap << "," << rhoVanc << "," << rhoVextrap << std::endl;
                
                if (Theta < critical_polynomial_Theta){
                    // The first fallback method is an extrapolation based on the closest converged expansion to the critical point
                    if (!std::isfinite(static_cast<double>(rhovec[0])) || rhovec[1] >= rhovec[0] || rhovec[0] < 0){
                        // And if that doesn't work, we use the critical extrapolation formula based on the expansion closest
                        // to the critical point that is fully converged
                        if (last_estimationL && last_estimationV){
                            double rhoLextrap = last_estimationL.value()(T), rhoVextrap = last_estimationV.value()(T);
                            rhovec = teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, rhoLextrap, rhoVextrap, 10);
                        }
                        else{
                            throw FailedIteration(T, "last_estimation not available @T="+std::to_string(T)+". Tcrittrue is "+std::to_string(Tcrittrue)+" K");
                        }
                    }
                }
                else{
                    throw FailedIteration(T, "Iteration failed below the polynomial @T="+std::to_string(T)+". Tcrittrue is "+std::to_string(Tcrittrue)+" K");
                }
                
                if (static_cast<double>(rhovec[0]) < 0 || (static_cast<double>(rhovec[0]) < static_cast<double>(rhovec[1]))){
                    throw FailedIteration(T, "Iteration failed @T="+std::to_string(T)+". Tcrittrue is "+std::to_string(Tcrittrue)+" K");
                }
                if (!std::isfinite(static_cast<double>(rhovec[0]))){
                    throw FailedIteration(T, "Iteration invalid liquid density @T="+std::to_string(T)+". Tcrittrue is "+std::to_string(Tcrittrue)+" K");
                }
            }
            // Calculate the pressures in each phase
            using tdxflt = teqp::TDXDerivatives<decltype(model), my_float_mp, Eigen::ArrayX<my_float_mp>>;
            const auto z = (Eigen::ArrayX<my_float_mp>(1) << 1.0).finished();
            auto pL = rhovec[0]*R*T*(1.0 + tdxflt::get_Ar01<teqp::ADBackends::multicomplex>(model, T, rhovec[0], z));
            auto pV = rhovec[1]*R*T*(1.0 + tdxflt::get_Ar01<teqp::ADBackends::multicomplex>(model, T, rhovec[1], z));
            // We have a good solution, store it in the database
            densitydb.insert(std::make_pair(T, DensitiesType{ rhovec[0], rhovec[1], pL, pV, rhocrittrue+BrhoL*pow(Theta, 0.5), rhocrittrue + BrhoV*pow(Theta, 0.5) }));
        }
        // Now we obtain values from the database
        auto d = densitydb.at(T); // Retrieve from cache
        return {static_cast<double>(d.rhoL), static_cast<double>(d.rhoV), static_cast<double>(d.pL)};
    };

    double Tmin = std::max(Ttriple, Tmin_sat), Tmax = Tcrittrue, tol = 1e-12;
    int N = 12, Msplit = 3, max_refine_passes = 12;
    
    std::vector<Container> last_good_exsL, last_good_exsV;
    
    // This callback function is called after each pass of refinement, to allow you to monitor the process,
    // or in this case, store diagnostic information about the last good expansion
    VectorDyadicSplittingCallback callback = [&last_good_exsL, &last_good_exsV, &last_estimationL, &last_estimationV, Msplit, tol, &BrhoL, &BrhoV, &Tcrittrue, &rhocrittrue, &model, &getdrhodTs](
      int num_pass, const Container& expansions)
    {
        std::cout << ".";
//        last_good_exsL = exsA;

        double T, rhoL, rhoV;

        // Work backwards since we start at the critical point
        for (int k = static_cast<int>(expansions[0].size())-1; k >= 0; --k) {
            // If is converged, stop, this is the one we seek
            auto ceL = expansions[0][k];
            auto ceV = expansions[1][k];
            if (are_converged(Msplit, tol, expansions, k)){
                T = ceL.xmax();
                auto Theta = (Tcrittrue-T)/Tcrittrue;
                rhoL = ceL.y_Clenshaw(T);
                rhoV = ceV.y_Clenshaw(T);
                // Find the curve for the critical region to be used for estimation
                double deltarhoL = rhoL-rhocrittrue;
                double deltarhoV = rhoV-rhocrittrue;
                auto [drhoLdT, drhoVdT] = getdrhodTs(model, T, rhoL, rhoV);
                double betaL = Theta/(-1/Tcrittrue)*drhoLdT/deltarhoL;
                double betaV = Theta/(-1/Tcrittrue)*drhoVdT/deltarhoV;
                auto BrhoL = deltarhoL/pow(Theta, betaL);
                auto BrhoV = deltarhoV/pow(Theta, betaV);
//                std::cout << T << ";" << rhoL << ";" << rhoV << ";" << deltarhoL << ";" << betaL << ";" << BrhoL << std::endl;
                last_estimationL = CriticalEstimation{BrhoL, betaL, Tcrittrue, rhocrittrue};
                last_estimationV = CriticalEstimation{BrhoV, betaV, Tcrittrue, rhocrittrue};
                break;
            }
        }
    };
    
    // Here we drive the fitting, using the custom dyadic splitting function
    // developed in this work
    Container exps;
    try {
        exps = vectored_dyadic_splitting(
            N,
            get_VLE_values,
            Tmin, Tmax, Msplit, tol, max_refine_passes, callback
        );
    }
    catch (FailedIteration&f) {
        if (f.T > (1-1e-9)*Tmax){
//            exps = last_good_exs;
        }
        else {
            throw;
        }
    }
    std::cout << std::endl;
    
    // Write out the expansions, metadata, and
    // anything else to the output file
    // in JSON format
    auto tovec = [](const Eigen::ArrayXd& a) {
        std::vector<double> z(a.size());
        for (auto i = 0; i < a.size(); ++i) { z[i] = a[i]; }
        return z;
    };
    nlohmann::json jcrit_anc = {
        { "cL", tovec(cLarray) },
        { "cV", tovec(cVarray) },
        { "Tc / K", Tcrittrue },
        { "rhoc / mol/m^3", rhocrittrue },
        { "Theta_min", critical_polynomial_Theta},
        { "_note", R"(coefficients are for the function like ln(|rho^A-rhoc|) = sum_i c_i ln(Theta)^i with Theta=(Tc-T)/Tc)" }
    };
    nlohmann::json meta = {
        { "Tcrit / K", Tcrit },
        { "Tcrittrue / K", Tcrittrue },
        { "Treducing / K", Treducing },
        { "Ttriple / K", Ttriple },
        { "rhocrittrue / mol/m^3", rhocrittrue },
        { "BrhoL / mol/m^3", BrhoL },
        { "BrhoV / mol/m^3", BrhoV },
        { "gas_constant / J/mol/K", R }
    };
    
    // Collect all the expansions
    nlohmann::json jexpansionsrhoL = nlohmann::json::array(),
                   jexpansionsrhoV = nlohmann::json::array(),
                   jexpansionsp = nlohmann::json::array();
    
    for (auto j = 0; j < exps[0].size(); ++j) {
        auto& exrhoL = exps[0][j];
        auto& exrhoV = exps[1][j];
        auto& expp = exps[2][j];
        jexpansionsrhoL.push_back({
            {"coef", tovec(exrhoL.coef())},
            {"xmin", exrhoL.xmin()},
            {"xmax", exrhoL.xmax()},
        });
        jexpansionsrhoV.push_back({
            {"coef", tovec(exrhoV.coef())},
            {"xmin", exrhoV.xmin()},
            {"xmax", exrhoV.xmax()},
        });
        jexpansionsp.push_back({
            {"coef", tovec(expp.coef())},
            {"xmin", expp.xmin()},
            {"xmax", expp.xmax()},
        });
    }
    nlohmann::json jo = {
        {"meta", meta},
        {"crit_anc", jcrit_anc},
        {"jexpansions_rhoL", jexpansionsrhoL},
        {"jexpansions_rhoV", jexpansionsrhoV},
        {"jexpansions_p", jexpansionsp},
    };
    // Stream the output into the file you specified
    std::ofstream ofs(ofpath); ofs << jo.dump(2);
    
    std::cout << exps[0].size() << " expansions" << std::endl;
}

/**
\brief Check the superancillaries
 
 This obtains the degree-doubled nodes for each expansions in the superancillary, and carries out a VLE calculation at the given point.  An output JSON data structure is written to file at the specified location
 
 \param fluid The FLD name coming from REFPROP
 \param input_file_path The path to the superancillaries to be loaded from, in JSON format
 \param outfile The path to the file to be written by this file
 */
void check_superancillaries(const std::string& fluid, const std::filesystem::path& input_file_path, const std::filesystem::path& outfile) {
    const std::filesystem::path fluid_json_path = teqp_datapath / "dev" / "fluids" / (fluid + ".json");
    auto model = teqp::build_multifluid_model({ fluid_json_path.string()}, teqp_datapath.string());

    auto get_collection = [](const std::filesystem::path & expansion_file) {
        const nlohmann::json jfile = teqp::load_a_JSON_file(expansion_file.string());
        auto parser_ = [&](const auto& key){
            std::vector<ChebTools::ChebyshevExpansion> o;
            for (const auto& ex : jfile.at(key)) {
                o.emplace_back(ex.at("coef").template get<std::vector<double>>(), ex.at("xmin"), ex.at("xmax"));
            }
            return ChebTools::ChebyshevCollection(o);
        };
        
        return std::make_tuple(parser_("jexpansions_rhoL"), parser_("jexpansions_rhoV"), parser_("jexpansions_p"));
    };

    auto db = nlohmann::json::array();

    // Load expansions from file for liquid and vapor
    auto [ccL, ccV, ccp] = get_collection(input_file_path);
    auto& ccL_ = ccL, ccV_ = ccV;
    auto meta = teqp::load_a_JSON_file(input_file_path.string())["meta"];
    double Tcrittrue = meta.at("Tcrittrue / K");
    double R = meta.at("gas_constant / J/mol/K");
    double rhocrittrue = meta.at("rhocrittrue / mol/m^3");
    
    auto get_degreedoubled_nodes = [&]() {
        std::vector<double> x;
        for (auto cc : { ccL_ }) {
            for (auto& ex : cc.get_exps()) {
                auto N = ex.coef().size() - 1;
                auto nodes_doubled = ChebTools::ChebyshevExpansion::factory(2*N, [](double x) { return x; }, ex.xmin(), ex.xmax()).get_nodes_realworld();
                for (auto& n : nodes_doubled) {
                    x.push_back(n);
                }
            }
        }
        std::sort(x.begin(), x.end());
        return x;
    };

    // Collect all the nodes from the expansions, and their degree-doubled in-between nodes
    std::vector<double> Tnodes = get_degreedoubled_nodes();
    
    for (auto T : Tnodes) {
        try {
            double Theta = (Tcrittrue-T)/Tcrittrue;
            
            double rhoSAL = ccL(T);
            double rhoSAV = ccV(T);
            double pSA = ccp(T);
            
            Eigen::ArrayXd rhovec;
            if (std::abs(T - Tcrittrue) > 1e-14 && T < ccL.get_exps().back().xmax() && T < ccV.get_exps().back().xmax()) {
                teqp::IsothermPureVLEResiduals<decltype(model), my_float_mp, teqp::ADBackends::multicomplex> residual(model, T);
                rhovec = teqp::do_pure_VLE_T<decltype(residual), my_float_mp>(residual, rhoSAL*(1+1e-5), rhoSAV*(1-1e-5), 10).cast<double>();
                if (!std::isfinite(rhovec[0])) {
                    throw FailedIteration(T, "Iteration failed @T=" + std::to_string(T) + " K for " + fluid + ". Tcrittrue is " + std::to_string(Tcrittrue) + " K.");
                }
            }
            else {
                rhovec = (Eigen::Array2d() << rhocrittrue, rhocrittrue).finished();
            }
            if (rhovec.size() != 2) {
                std::cout << "rhovec is not 2 elements in length" << std::endl;
            }
            
            // Calculate the pressures in each phase
            using tdxflt = teqp::TDXDerivatives<decltype(model), my_float_mp, Eigen::ArrayX<my_float_mp>>;
            const auto z = (Eigen::ArrayX<my_float_mp>(1) << 1.0).finished();
            auto pmp = static_cast<double>(rhovec[0]*R*T*(1.0 + tdxflt::get_Ar01<teqp::ADBackends::multicomplex>(model, T, rhovec[0], z)));
            
            db.push_back(nlohmann::json{
                {"T / K", T},
                {"rho'(SA) / mol/m^3", rhoSAL},
                {"rho'(mp) / mol/m^3", rhovec[0]},
                {"rho'(SA)/rho'(mp)", rhoSAL/rhovec[0]},
                {"rho''(SA) / mol/m^3", rhoSAV},
                {"rho''(mp) / mol/m^3", rhovec[1]},
                {"rho''(SA)/rho''(mp)", rhoSAV/rhovec[1]},
                {"p(mp) / Pa", pmp},
                {"p(SA) / Pa", pSA}
            });
        }
        catch (FailedIteration& f) {
            if (f.T > 0.9999 * Tcrittrue) {

            }
            else {
                std::cout << f.msg << std::endl;
            }
            db.push_back(nlohmann::json{
                {"T / K", T},
                {"errmsg", f.msg}
            });
        }
        catch(std::exception& e){
            std::cout << e.what() << std::endl;
            db.push_back(nlohmann::json{
                {"T / K", T},
                {"errmsg", e.what()}
            });
        }
    }
    
    // Return results
    nlohmann::json jo = {
        {"meta", meta},
        {"data", db}
    };
    std::ofstream ofs(outfile); ofs << jo.dump(2);
}
