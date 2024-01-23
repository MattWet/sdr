#include <Rcpp.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace Rcpp;
using namespace std;

// Small function which returns our magic number
int my_magic_number() {
    return 74143467; // Reto's credit card number
}

// Writes magic number into current ifstream
void write_magic_number(std::ofstream& ofile, int bytes = 4) {
    int m = my_magic_number();
    ofile.write(reinterpret_cast<const char*>(&m), bytes);
}
    
// Reading magic number from current position of ifstream
int read_magic_number(std::ifstream& ifile, int bytes = 4) {
    int res;
    ifile.read(reinterpret_cast<char*>(&res), sizeof(res));
    return res;
}

void write_string(std::ofstream& ofile, const std::string& val, const int bytes) {
    ofile.write(val.c_str(), bytes);
}

// Reading std::string of length nchar from ifstream
std::string read_string(std::ifstream& ifile, const int nchar) {
    std::string result(nchar, '\0'); // Create string with 'nchar' chars
    ifile.read(&result[0], nchar);
    return result;
}

// Parameters
// ----------
// ofile: Binary file stream (out | binary)
// value: double, value to write to binary stream.
// bytes: 4 (float) or 8 (double)
//
// Return
// ------
// Always returns double.
void write_double(std::ofstream& ofile, double value, const int &bytes) {
    // Encode the double value and write to the file
    if (bytes == 4) {
        float enc_value = static_cast<float>(value);
        ofile.write(reinterpret_cast<char*>(&enc_value), bytes);
    } else {
        ofile.write(reinterpret_cast<char*>(&value), bytes);
    }
}

// Parameters
// ----------
// ifile: Binary file stream (in | binary)
// bytes: 4 (float) or 8 (double)
//
// Return
// ------
// Always returns double.
double read_double(std::ifstream& ifile, const int &bytes) {
    double result;
    if (bytes == 4) {
        float enc_value;
        ifile.read(reinterpret_cast<char*>(&enc_value), bytes);
        result = static_cast<double>(enc_value);
    } else {
        ifile.read(reinterpret_cast<char*>(&result), bytes);
    }
    return result;
}


void write_int(std::ofstream& ofile, int val) {
    ofile.write(reinterpret_cast<const char*>(&val), sizeof(int));
}

int read_int(std::ifstream& ifile) {
    int res;
    ifile.read(reinterpret_cast<char*>(&res), sizeof(res));
    return res;
}

int get_bytes(std::string& type) {

    int res = -999;
    if (type == "float") {
        res = 4;
    } else if (type == "double") {
        res = 8;
    } else {
        std::cerr << "Unrecognized type!";
    }
    return res;
}


// Small helper function; checks if the maximum absolute value
// found in the data exceeds what we can store with X bytes
// (controlled by the type argument).
// If exceeding, throw an error.
void check_data_range(std::string& type, double& maxabs) {

    bool error = false;
    if (type == "float") {
        float limit = std::numeric_limits<float>::max();
        if (maxabs > limit) error = true;
    } else if (type == "double") {
        double limit = std::numeric_limits<double>::max();
        if (maxabs > limit) error = true;
    }

    if (error) {
        std::string errormsg = "Max absolute value (" + to_string(maxabs) + ") exceeds what can be stored as " + type + "; exit (change type)";
        std::cerr << errormsg;
    }
}


// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// Parameters
// ----------
// file : Path/name to the csv file to be read.
// binfile: Name of the binary output file created by this function.
// type: Can be 'float' or 'double'. Will write the
//       data in 4 bytes or 8 bytes respectively.
// skip: Number of lines to skip at the start of the CSV file.
// header: true if there is a header line with variable description, else false.
// sep: Single character, value separator.
// verbose: true/false, increases verbosity if set true.

// [[Rcpp::export]]
Rcpp::List create_binmm(std::string file,
                        std::string binfile,
                        std::string type = "double",
                        int skip = 0,
                        bool header = true,
                        char sep = ',', bool verbose = false) {

    int ncol = 0, nrow = 0;
    int nchar_filename = file.size(); // Number of character of input file name
    std::string line;
    std::string val;

    // Getting bytes to write the data later on (float, double)
    int bytes = get_bytes(type);
    if (verbose) {
        Rcout << "[cpp] Bytes for numerics: " << bytes << " (" << type << ")\n";

        Rcout << "[cpp] Range limit for " << type << ": " << std::numeric_limits<double>::max() << "  " << std::numeric_limits<double>::min() << "\n";
    }


    // Open input file connection
    if (verbose) Rcout << "[cpp] Reading file " << file.c_str() << " (filename len = " << nchar_filename << ")\n";
    std::ifstream ifile(file.c_str());
    if (!ifile) std::cerr << "Whoops, input file not found/problem to open stream.";

    // Skipping lines
    if (skip > 0) {
        for (int i = 0; i < skip; i++) {
            // Will be 'wasted'; do we need to read the full line?
            // Btw: '\n' is crucial ("\n" would not work).
            std::getline(ifile, line, '\n');
        }
    }

    std::stringstream iss(""); // Must be 'globally' available

    // Reading first line to count number of columns.
    // This can be a header line (if header = true we will extract
    // the column names from it) or a data line. To be able to catch
    // both cases we need to store the file pointer to revert in case
    // header = false to re-read the same line as 'data'.
    std::streampos line1 = ifile.tellg(); // Store position of first line
    std::getline(ifile, line, '\n'); // Reading line
    iss.clear(); iss.str(line);

    //std::stringstream iss(line);
    // First run: count number of columns
    while (std::getline(iss, val, sep)) {
        ncol++;
    }
    if (verbose) Rcout << "[cpp] Number of columns: " << ncol << "\n";

    // Open output stream for binary data
    std::ofstream ofile(binfile.c_str(), std::ios::binary);
    if (!ofile) std::cerr << "Problems opening output file";

    // Writing binary file 'meta'
    // - magic_number: 4 byte int
    // - nrow:            4 byte int
    // - ncol:            4 byte int
    // - file_name_bytes: 4 byte int
    // - file name:       8 * file_name_bytes char
    write_magic_number(ofile);    // 4 byte integer
    write_int(ofile, -999);       // nrow, will be replaced at the end as soon as we know
    write_int(ofile, ncol);       // ncol
    write_int(ofile, bytes);      // Number of bytes for data
    write_int(ofile, nchar_filename); // Length of original file name
    write_string(ofile, file.c_str(), nchar_filename); // Original file name
    
    // Creating vector for column names as well as for the length
    // of the colnames (number of characters of each name).
    CharacterVector colnames(ncol);

    // Extracting header information
    if (header) {
        // Create new CharacterVector and store the result
        // Clear and reset (processing the same line once again)
        iss.clear();
        iss.str(line);
        // First run: count number of columns
        for (int i = 0; i < ncol; i++) {
            std::getline(iss, val, ',');
            // Remove double quotes if there are any
            val.erase(remove(val.begin(), val.end(), '"'), val.end());
            // Remove single quotes if there are any
            val.erase(remove(val.begin(), val.end(), '\''), val.end());
            colnames[i] = val.c_str();
        }
    } else {
        for (int i = 0; i < ncol; i++) {
            colnames[i] = "V" + to_string(i);
        }
        // Revert to first line as the first line (we already used
        // once to count number of columns) as a data line.
        if (verbose) Rcout << "[cpp] Revert line\n";
        ifile.seekg(line1);
    }

    // Store number of characters for longest colname; used
    // for byte length for binary writing/reading later.
    if (verbose) Rcout << "[cpp] Writing length and content of column names\n";
    for (int i = 0; i < ncol; i++)
        write_int(ofile, colnames[i].size());
    for (int i = 0; i < ncol; i++)
        write_string(ofile, Rcpp::as<std::string>(colnames[i]), colnames[i].size());


    // Keep track of current std::string pos; insert
    // 2 * ncol missing values (-999) to be filled with
    // ncol * means (row-wise mean) and ncol * sd (row-wise standard
    // devation) at the end.
    std::streampos scale_pos = ofile.tellp();
    for (int i = 0; i < 2 * ncol; i ++) write_double(ofile, -999, bytes); // Dummies

    // Creating numeric vectors to store mean and
    NumericVector mean = rep(0.0, ncol);
    NumericVector m2   = rep(0.0, ncol); // Used to calculate sd at the end
    NumericVector sd   = rep(0.0, ncol);

    // Reading data; read line by line until we find EOF.
    // At the same time we calculate the row sums and
    // squared row sums for mean/sd.
    // TODO(R): Test what happens if we find non-numeric
    //          elements and if the length of the row (number
    //          ov values) is not equal to ncol!

    // Keeping track of the max(abs(value)) to do a range check
    // to see if we can accurately store this given the type
    // used.

    double x, xabs, maxabs = 0.0;
    double d1, d2;

    if (verbose) Rcout << "[cpp] Processing data ...\n";
    while (std::getline(ifile, line, '\n')) {
        nrow++;
        iss.clear(); iss.str(line);
        for (int j = 0; j < ncol; j++) {
            std::getline(iss, val, ',');
            x    = std::stod(val);
            xabs = std::abs(x);
            if (xabs > maxabs) maxabs = xabs;
            write_double(ofile, x, bytes);

            // Updating mean/m2 for mean/sd; in-line algorithm
            // to avoid calculating sums to not get outside limits of double
            d1       = x - mean[j];  mean[j] += d1 / nrow;
            d2       = x - mean[j];  m2[j]   += d1 * d2;
        }
    }

    // Calculating standard deviation
    for (int i = 0; i < ncol; i++) sd[i] = sqrt(m2[i] / (nrow - 1));

    // Revert to scale_pos and write all means, then all standard deviations.
    ofile.seekp(scale_pos);
    for (int j = 0; j < ncol; j++)  write_double(ofile, mean[j], bytes);
    for (int j = 0; j < ncol; j++)  write_double(ofile, sd[j], bytes);

    if (verbose) Rcout << "[cpp] Maximum absolute value: " << maxabs << "\n";
    check_data_range(type, maxabs);

    if (verbose) Rcout << "[cpp] Dimension found/read: " << nrow << " x " << ncol << "\n";
    if (verbose) Rcout << "[cpp] Closing file connections\n";

    ifile.close();

    // Setting number of rows (after 4 bytes; after magic number)
    ofile.seekp(4); 
    write_int(ofile, nrow);

    ofile.close();

    if (verbose) Rcout << "[cpp] Number of rows: " << nrow << "\n";


    if (verbose) Rcout << "[cpp] All done, creating return object ...\n";

    mean.attr("names") = colnames;
    sd.attr("names") = colnames;
    List dim = Rcpp::List::create(Named("nrow") = nrow, Named("ncol") = ncol);

    List rval = Rcpp::List::create(Named("original_file") = file,
                                   Named("binfile")       = binfile,
                                   Named("dim")           = dim,
                                   Named("colnames")      = colnames,
                                   // mean and standard deviation for scaling
                                   Named("scale") = Rcpp::List::create(Named("mean") = mean,
                                                                       Named("sd") = sd)
                                  );

    // Dummy class; typically not used!
    rval.attr("class") = "binmm_create";
    return rval;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List meta_binmm(std::string file, bool verbose = false) {


    if (verbose) Rcout << "[cpp] Reading file " << file.c_str() << "\n";

    int j, nrow, ncol, bytes, nchar_filename;

    // Open input file connection (binary)
    std::ifstream ifile(file.c_str(), ios::in | std::ios::binary);
    if (!ifile)
        std::cerr << "Whoops, input file not found/problem to open stream.";

    // Reading binary file meta
    // - magic_number: 4 byte int
    // - nrow: 4 byte int
    // - ncol: 4 byte int
    // - bytes: 4 byte int, bytes used for writing the data
    // - file_name_bytes: 4 byte int
    int magic_number = read_magic_number(ifile);
    if (magic_number != my_magic_number())
        std::cerr << "Wrong magic number; content of binary file not what is expected";
    if (verbose) Rcout << "[cpp] Got magic number " << magic_number << "\n";

    nrow = read_int(ifile);
    ncol = read_int(ifile);
    if (verbose) Rcout << "[cpp] Total dim: " << nrow << " x " << ncol << "\n";

    // Bytes used to write the data
    bytes    = read_int(ifile);
    if (verbose) Rcout << "[cpp] Data size (bytes): " << bytes << "\n";

    nchar_filename = read_int(ifile);
    std::string original_file = read_string(ifile, nchar_filename);

    if (verbose) Rcout << "[cpp] Original file: " << original_file << " (" << nchar_filename << ")\n";

    // Reading column names
    // First nchar of each of the colum names, then string (name)
    IntegerVector nchar_colnames = rep(-9, ncol);
    CharacterVector colnames(ncol);
    for (j = 0; j < ncol; j++)  nchar_colnames[j] = read_int(ifile);
    for (j = 0; j < ncol; j++)  colnames[j] = read_string(ifile, nchar_colnames[j]);

    // Reading scaling parameters; mean and standard deviation;
    // two numeric vectors following the column names.
    NumericVector mean = rep(-999., ncol);
    NumericVector sd   = rep(-999., ncol);
    for (j = 0; j < ncol; j++)  mean[j] = read_double(ifile, bytes);
    for (j = 0; j < ncol; j++)  sd[j]   = read_double(ifile, bytes);

    // Adding names
    mean.attr("names") = colnames;
    sd.attr("names")   = colnames;
    
    // Pointer position start of data
    std::streampos pointer_pos = ifile.tellg();

    Rcpp::List dim   = Rcpp::List::create(Named("nrow") = nrow, Named("ncol") = ncol);
    Rcpp::List scale = Rcpp::List::create(Named("mean") = mean, Named("sd") = sd);


    List rval = Rcpp::List::create(
                        Named("original_file") = original_file,
                        Named("binfile")       = file,
                        Named("dim")           = dim,
                        Named("colnames")      = colnames,
                        Named("scale")         = scale,
                        Named("bytes")         = bytes,
                        Named("pointer_pos") = static_cast<int>(pointer_pos));
    rval.attr("class") = "binmm";
    return rval;
}


// -------------------------------------------------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix subset_binmm(std::string file,
                        Rcpp::IntegerVector& ii, Rcpp::IntegerVector& jj,
                        bool standardize = false, bool verbose = false) {


    if (verbose) Rcout << "[cpp] Loading meta information ...\n";

    Rcpp::List meta = meta_binmm(file, verbose);
    Rcpp::List dim = meta["dim"];
    Rcpp::CharacterVector colnames = meta["colnames"];

    Rcpp::List          scale = meta["scale"];
    Rcpp::NumericVector mean  = scale["mean"];
    Rcpp::NumericVector sd    = scale["sd"];

    int ncol = dim["ncol"];
    int bytes = meta["bytes"], pointer_pos = meta["pointer_pos"];
    int i, j;

    // Need to open the file again unfortunately ...
    // TODO(R): This is a drawback of this method separating meta and data
    std::ifstream ifile(file.c_str(), ios::in | std::ios::binary);
    if (!ifile)
        std::cerr << "Whoops, input file not found/problem to open stream.";

    // Create numeric matrix of the dimension requested
    // by the user; loop over i/j indices (0 based) and
    // read the data.
    // Sets pointers to jump to the value of interest
    // (long stream of doubles).
    NumericMatrix rmat(ii.size(), jj.size());

    // Appending dimension names
    CharacterVector rmat_colnames(jj.size());
    for (j = 0; j < jj.size(); ++j) rmat_colnames[j] = colnames[jj[j]];
    rmat.attr("dimnames") = Rcpp::List::create(R_NilValue, rmat_colnames);

    if (verbose & !standardize) {
        Rcout << "[cpp] Reading data\n";
    } else if (verbose) {
        Rcout << "[cpp] Reading data and standardize columns\n";
    }

    int pos; // Pointer position in bytes
    for (i = 0; i < ii.size(); ++i) {
        for (j = 0; j < jj.size(); ++j) {
            // Calculate pointer position
            pos = pointer_pos + ((ii[i] * ncol) + jj[j]) * bytes; //sizeof(val);
            // Setting pointer; read and store double
            ifile.seekg(pos);
            rmat(i, j) = read_double(ifile, bytes);
            if (standardize) rmat(i, j) = (rmat(i, j) - mean[jj[j]]) / sd[jj[j]];
        }
    }

    if (verbose) Rcout << "[cpp] Closing file connection\n";
    ifile.close();

    // Dummy return
    return rmat;
}


