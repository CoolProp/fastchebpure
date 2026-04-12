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

// Prototype for builder and checker
void build_superancillaries(const std::string &, const std::filesystem::path &);
void check_superancillaries(const std::string &, const std::filesystem::path&, const std::filesystem::path&);

int main(int argc, char** argv){

    CLI::App app{
        "fastchebpure - build and validate Chebyshev superancillary equations\n"
        "for the saturation properties of pure fluids.",
        "fitcheb"
    };
    app.set_version_flag("--version", "1.0.0");

    // --- Fluid selection ---
    std::vector<std::string> fluids_arg;
    app.add_option("-f,--fluid", fluids_arg,
        "Fluid(s) to process (stem name without .json, e.g. 'Water').\n"
        "May be specified multiple times. Default: all fluids in the data path.");

    // --- Skip list ---
    std::vector<std::string> skip_arg = {"Air", "SES36"};
    app.add_option("-s,--skip", skip_arg,
        "Fluid(s) to skip. Default: Air SES36.")
        ->capture_default_str();

    // --- Paths ---
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

    // --- Thread count ---
    int nthreads = 6;
    app.add_option("-j,--jobs", nthreads,
        "Number of parallel worker threads (0 = serial).")
        ->capture_default_str();

    // --- Mode ---
    std::string mode = "both";
    app.add_option("-m,--mode", mode,
        "Operation mode: 'build', 'check', or 'both' (default).")
        ->capture_default_str()
        ->check(CLI::IsMember({"build", "check", "both"}));

    // --- Force rebuild ---
    bool force = false;
    app.add_flag("--force", force,
        "Overwrite existing output files instead of skipping them.");

    CLI11_PARSE(app, argc, argv);

    // Resolve effective paths (command-line overrides compile-time defaults)
    std::filesystem::path eff_datapath{data_path_arg};
    std::filesystem::path eff_output{output_arg};
    std::filesystem::path eff_check{check_arg};

    // Validate directories
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
        // Explicit list provided on the command line
        for (auto& f : fluids_arg) {
            if (skip_set.count(f)) {
                std::cout << "Skipping (skip list): " << f << "\n";
                continue;
            }
            fluids_to_run.push_back(f);
        }
    } else {
        // Scan the data directory
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

    // Build a job lambda for each fluid
    auto make_job = [&](const std::string& fluid) {
        return [fluid, &eff_output, &eff_check, &mode, force]() {
            auto outfile   = eff_output / (fluid + "_exps.json");
            auto checkfile = eff_check  / (fluid + "_check.json");
            try {
                if (mode == "build" || mode == "both") {
                    if (force || !std::filesystem::exists(outfile)) {
                        std::cout << "Building -> " << outfile.filename().string() << "\n";
                        build_superancillaries(fluid, outfile);
                    } else {
                        std::cout << "Skipping (exists): " << outfile.filename().string() << "\n";
                    }
                }
                if (mode == "check" || mode == "both") {
                    if (force || !std::filesystem::exists(checkfile)) {
                        std::cout << "Checking -> " << checkfile.filename().string() << "\n";
                        check_superancillaries(fluid, outfile, checkfile);
                    } else {
                        std::cout << "Skipping (exists): " << checkfile.filename().string() << "\n";
                    }
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
        // Serial execution
        for (auto& fluid : fluids_to_run) {
            make_job(fluid)();
        }
    }

    return EXIT_SUCCESS;
}
