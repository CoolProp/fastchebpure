// C++ standard library
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <set>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include <CLI/CLI.hpp>

// Paths are defined in the actual code (fastcheb.cpp)
extern const std::filesystem::path teqp_datapath;
extern const std::filesystem::path output_prefix;
extern const std::filesystem::path check_destination;

// Prototypes for builder, checker, injector
void build_superancillaries(const std::string&, const std::filesystem::path&, const std::filesystem::path&);
void check_superancillaries(const std::string&, const std::filesystem::path&, const std::filesystem::path&, const std::filesystem::path&);
void inject_superancillary(const std::string&, const std::filesystem::path&, const std::filesystem::path&, const std::filesystem::path&, bool, const std::vector<double>&);

int main(int argc, char** argv){

    CLI::App app{
        "fastchebpure - build, validate, and inject Chebyshev superancillary\n"
        "equations for the saturation properties of pure fluids.\n\n"
        "Subcommands (run `fitcheb <cmd> --help` for details):\n"
        "  fit      Fit expansions and write {output}/{fluid}_exps.json\n"
        "  check    Validate expansions and write {checkdir}/{fluid}_check.json\n"
        "  inject   Embed expansions + check points into CoolProp dev/fluids/{fluid}.json\n\n"
        "With no subcommand, `fit` and `check` are run for each selected fluid.",
        "fitcheb"
    };
    app.set_version_flag("--version", "1.0.0");
    app.require_subcommand(0, 1);

    // --- Shared options (defined once on the app; visible to all subcommands) ---
    std::vector<std::string> fluids_arg;
    app.add_option("-f,--fluid", fluids_arg,
        "Fluid(s) to process (stem name without .json, e.g. 'Water').\n"
        "May be specified multiple times. Default: all fluids in the data path.");

    std::vector<std::string> skip_arg = {"Air", "SES36"};
    app.add_option("-s,--skip", skip_arg,
        "Fluid(s) to skip. Default: Air SES36.")
        ->capture_default_str();

    std::string data_path_arg = teqp_datapath.string();
    app.add_option("-d,--datapath", data_path_arg,
        "Path to the CoolProp data directory (must contain dev/fluids/).")
        ->capture_default_str();

    std::string output_arg = output_prefix.string();
    app.add_option("-o,--output", output_arg,
        "Directory for the generated expansion JSON files.")
        ->capture_default_str();

    std::string check_arg = check_destination.string();
    app.add_option("-c,--checkdir", check_arg,
        "Directory for the validation check JSON files.")
        ->capture_default_str();

    int nthreads = 6;
    app.add_option("-j,--jobs", nthreads,
        "Number of parallel worker threads (0 = serial).")
        ->capture_default_str();

    bool force = false;
    app.add_flag("--force", force,
        "fit/check: overwrite existing outputs instead of skipping them.\n"
        "inject: inject even when source_eos_hash disagrees with the destination EOS hash.");

    std::vector<double> thetas_vec = {0.5, 0.3, 0.1};
    app.add_option("--thetas", thetas_vec,
        "Theta = (T_c - T)/T_c values at which to sample check_points (inject only).\n"
        "Default: 0.5 0.3 0.1  (three points spanning the SA usable range).")
        ->capture_default_str();

    // --- Subcommands (no per-command options — all options are shared) ---
    auto* cmd_fit    = app.add_subcommand("fit",    "Fit superancillary expansions and write {fluid}_exps.json.");
    auto* cmd_check  = app.add_subcommand("check",  "Validate expansions against the EOS and write {fluid}_check.json.");
    auto* cmd_inject = app.add_subcommand("inject", "Inject expansions + check points into CoolProp dev/fluids/{fluid}.json.");
    cmd_fit->fallthrough();
    cmd_check->fallthrough();
    cmd_inject->fallthrough();

    CLI11_PARSE(app, argc, argv);

    // Determine which phases to run. With no subcommand, default to fit+check
    // (the historical "both" mode); inject is never run implicitly.
    const bool do_fit    = cmd_fit->parsed()    || (!cmd_check->parsed() && !cmd_inject->parsed());
    const bool do_check  = cmd_check->parsed()  || (!cmd_fit->parsed()   && !cmd_inject->parsed());
    const bool do_inject = cmd_inject->parsed();

    // Resolve effective paths
    std::filesystem::path eff_datapath{data_path_arg};
    std::filesystem::path eff_output{output_arg};
    std::filesystem::path eff_check{check_arg};

    for (auto [label, path] : {
        std::pair{"data path", eff_datapath},
        std::pair{"output directory", eff_output},
        std::pair{"check directory", eff_check}
    }) {
        if (!std::filesystem::exists(path)) {
            std::cerr << "Error: " << label << " does not exist: " << path << "\n";
            return EXIT_FAILURE;
        }
    }

    std::filesystem::path fluids_dir = eff_datapath / "dev" / "fluids";
    if (!std::filesystem::exists(fluids_dir)) {
        std::cerr << "Error: fluids directory does not exist: " << fluids_dir << "\n";
        return EXIT_FAILURE;
    }

    // Build the set of fluids to skip
    std::set<std::string> skip_set(skip_arg.begin(), skip_arg.end());

    // Collect fluid names to process
    std::vector<std::string> fluids_to_run;
    if (!fluids_arg.empty()) {
        for (auto& f : fluids_arg) {
            if (skip_set.count(f)) {
                std::cout << "Skipping (skip list): " << f << "\n";
                continue;
            }
            fluids_to_run.push_back(f);
        }
    } else {
        for (auto const& dir_entry : std::filesystem::directory_iterator{fluids_dir}) {
            if (!dir_entry.is_regular_file()) { continue; }
            if (dir_entry.path().extension() != ".json") { continue; }
            auto stem = dir_entry.path().stem().string();
            if (skip_set.count(stem)) { continue; }
            fluids_to_run.push_back(stem);
        }
    }

    if (fluids_to_run.empty()) {
        std::cerr << "No fluids to process.\n";
        return EXIT_FAILURE;
    }

    auto make_job = [&](const std::string& fluid) {
        return [fluid, &eff_output, &eff_check, &eff_datapath, do_fit, do_check, do_inject, force, &thetas_vec]() {
            auto outfile   = eff_output / (fluid + "_exps.json");
            auto checkfile = eff_check  / (fluid + "_check.json");
            auto fluidfile = eff_datapath / "dev" / "fluids" / (fluid + ".json");
            try {
                if (do_fit) {
                    if (force || !std::filesystem::exists(outfile)) {
                        std::cout << "Building -> " << outfile.filename().string() << "\n";
                        build_superancillaries(fluid, outfile, eff_datapath);
                    } else {
                        std::cout << "Skipping (exists): " << outfile.filename().string() << "\n";
                    }
                }
                if (do_check) {
                    if (force || !std::filesystem::exists(checkfile)) {
                        std::cout << "Checking -> " << checkfile.filename().string() << "\n";
                        check_superancillaries(fluid, outfile, checkfile, eff_datapath);
                    } else {
                        std::cout << "Skipping (exists): " << checkfile.filename().string() << "\n";
                    }
                }
                if (do_inject) {
                    inject_superancillary(fluid, outfile, checkfile, fluidfile, force, thetas_vec);
                }
            } catch (const std::exception& e) {
                std::cerr << "[" << fluid << "]: " << e.what() << "\n";
            }
        };
    };

    if (nthreads > 0) {
        boost::asio::thread_pool pool(static_cast<std::size_t>(nthreads));
        for (auto& fluid : fluids_to_run) {
            std::cout << "Submitting: " << fluid << "\n";
            boost::asio::post(pool, make_job(fluid));
        }
        pool.join();
    } else {
        for (auto& fluid : fluids_to_run) {
            make_job(fluid)();
        }
    }

    return EXIT_SUCCESS;
}
