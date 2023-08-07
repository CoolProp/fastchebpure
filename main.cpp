// C++ standard library
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

// Paths are defined in the actual code
extern const std::string teqp_datapath;
extern const std::string output_prefix;
extern const std::string check_destination;

const int FASTCHEB_PROCESSORS = 6;

// Prototype for builder and checker
void build_superancillaries(const std::string &, const std::string &);
void check_superancillaries(const std::string &, const std::string &, const std::string&);

int main(){
    
    // Check for existence of output locations, they must be present
    if (!std::filesystem::exists(output_prefix)){
        std::cout << "output prefix doesn't exist: " << output_prefix << std::endl;
        return EXIT_FAILURE;
    }
    if (!std::filesystem::exists(check_destination)) {
        std::cout << "output checkfile destination doesn't exist: " << check_destination << std::endl;
        return EXIT_FAILURE;
    }

    // Launch the pool with desired number of threads.
    boost::asio::thread_pool pool(FASTCHEB_PROCESSORS);
    
    for (auto const& dir_entry : std::filesystem::directory_iterator{ teqp_datapath + "/dev/fluids" }){
        if (dir_entry.is_regular_file()) {
            auto fluid = dir_entry.path().stem().string();
            
//            if (!(fluid == "WATER")){ continue; } // Rudimentary filtering for testing purposes
            
            // Skip the .DS_Store file on OSX
            if (fluid == ".DS_Store"){ continue; }
            
            auto job = [fluid]() {
                auto outfile_path = output_prefix + fluid + "_exps.json";
                auto checkfile_path = output_prefix + "/check/" + fluid + "_check.json";
                
                try {
                    if (!std::filesystem::exists(outfile_path)) {
                        std::cout << "Building ->: " << outfile_path << std::endl;
                        build_superancillaries(fluid, outfile_path);
                    }
                    if (!std::filesystem::exists(checkfile_path)){
                        std::cout << "Checking ->: " << checkfile_path << std::endl;
                        check_superancillaries(fluid, outfile_path, checkfile_path);
                    }
                }
                catch (const std::exception& e) {
                    std::cout << "[" << fluid << "]: " << e.what() << std::endl;
                }
            };

            if (FASTCHEB_PROCESSORS > 0) {
                // Submit a lambda object to the pool.
                std::cout << "Submitting: " << fluid << std::endl;
                boost::asio::post(pool, job);
            }
            else {
                // Run the job serially
                job();
            }
        }
    }
    pool.join();
    return EXIT_SUCCESS;
}
