#include <Rcpp.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace Rcpp;
using namespace std;


// retoMat(file, skip = 0, header = true, verbose = true, sep = ',')
// 
// Parameters
// ----------
// file : std::string
//      Path/name of the CSV file to be read.
// skip : int
//      Positive integer, number of lines to be skipped in the CSV.
// header : bool
//      If true it is expected that the first line we find is the header line.
// sep : char
//      Separator of the data, defaults to ','.
// verbose : bool
//      Creates some output if true.
//
// Return
// ------
// List : Returns a list with information about the dimension of the data set,
// column names, and row-wise mean and standard deviation.

// [[Rcpp::export]]
Rcpp::List retoMat(std::string file, int skip = 0,
                    bool header = true,
                    char sep = ',', bool verbose = true) {

    if (verbose)
        Rcout << "[cpp] Reading file " << file.c_str() << "\n";

    int ncol = 0, nrow = 0;

    // Open file connection
    std::ifstream myfile(file.c_str());
    std::string line;
    std::string val;

    // Skipping lines
    if (skip > 0) {
        for (int i = 0; i < skip; i++) {
            // Will be 'wasted'; do we need to read the full line?
            // Btw: '\n' is crucial ("\n" would not work).
            std::getline(myfile, line, '\n');
        }
    }

    std::stringstream iss(""); // Must be 'globally' available

    // Header?
    std::getline(myfile, line, '\n');
    iss.clear(); iss.str(line);
    //std::stringstream iss(line);
    // First run: count number of columns
    while (std::getline(iss, val, sep)) {
        ncol++;
        //Rcout << " | " << val << "\n";
    }
    if (verbose)
        Rcout << "[cpp] Number of columns: " << ncol << "\n";
    
    // Creating colnames vector (length ncol) and fill.
    CharacterVector colnames(ncol);

    // Creating numeric vectors to store row sums and
    // squared row sums to calculate mean and standard deviation (sd).
    // We could re-use some of them, but they are fairly small
    // (numeric vectors of length ncol).
    NumericVector xsum        = rep(0.0, ncol);
    NumericVector xsumsquared = rep(0.0, ncol);
    NumericVector mean        = rep(0.0, ncol);
    NumericVector sd          = rep(0.0, ncol);

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
            Rcout << i;
            colnames[i] = "V" + to_string(i);
        }
    }
    //Rcout << "[cpp] Colnames: " << colnames << "\n";

    // Reading data; read line by line until we find EOF.
    // At the same time we calculate the row sums and
    // squared row sums for mean/sd.
    // TODO(R): Test what happens if we find non-numeric
    //          elements and if the length of the row (number
    //          ov values) is not equal to ncol!
    double x = 0.;
    if (verbose)
        Rcout << "[cpp] Processing data ...\n";
    while (std::getline(myfile, line, '\n')) {
        nrow++;
        iss.clear(); iss.str(line);
        for (int i = 0; i < ncol; i++) {
            std::getline(iss, val, ',');
            x = std::stod(val);
            xsum[i]        = xsum[i] + x;
            xsumsquared[i] = xsumsquared[i] + x * x;
        }
    }

    if (verbose)
        Rcout << "[cpp] Number of rows: " << nrow << "\n";

    // Calculate mean and standard deviation (sd)
    if (verbose)
        Rcout << "[cpp] Calculating mean and standard deviation (row-wise)\n";
    for (int i = 0; i < ncol; i++) {
        sd[i]    = sqrt((xsumsquared[i] - xsum[i] * xsum[i] / nrow) / nrow);
        mean[i]  = xsum[i] / nrow;
    }

    if (verbose)
        Rcout << "[cpp] All done, creating return object ...\n";

    mean.attr("names") = colnames;
    sd.attr("names") = colnames;
    List dim = Rcpp::List::create(Named("nrow") = nrow, Named("ncol") = ncol);

    //Rcpp::List res = Rcpp::List::create(colnames);
    List rval = Rcpp::List::create(Named("file") = file,
                                   Named("header") = header,
                                   Named("skip") = skip,
                                   Named("dim")  = dim,
                                   Named("colnames") = colnames,
                                   // mean and standard deviation for scaling
                                   Named("scale") = Rcpp::List::create(Named("mean") = mean,
                                                                       Named("sd") = sd)
                                  );
    rval.attr("class") = "retoMat";
    return rval;
}

// Helper function to check if integer 'b' is in the
// IntegerVector 'a'. Returns true if found, false otherwise.
bool contains(IntegerVector a, int b) {
    return std::find(a.begin(), a.end(), b)!=a.end();
}



// retoMat_subset(x, i, j, standardize = false, sep = ',', verbose = false)
// 
// Parameters
// ----------
// x : List
//      Object of class retoMat with data set meta information.
// i : IntegerVector
//      Index of rows to be read.
// j : IntegerVector
//      Index of columns to be read.
// standardize : bool
//      If set true, the values will be standardized using the row-wise
//      mean and sd stored in x.
// sep : char
//      Separator of the data, defaults to ','.
// verbose : bool
//      Creates some output if true.
//
// Return
// ------
// Numericmatrix : Matrix of dimension i.size() x j.size() with
// column names (row names empty).


// [[Rcpp::export]]
Rcpp::NumericMatrix retoMat_subset(const List & x, IntegerVector i, IntegerVector j,
        bool standardize = false, char sep = ',', bool verbose = false) {

    int ni = i.size(), nj = j.size();
    Rcpp::List dim = x["dim"];
    int ncol = as<int>(dim["ncol"]);
    int nrow = as<int>(dim["nrow"]);
    int skip = as<int>(x["skip"]);
    bool header = as<bool>(x["header"]);
    std::string file = x["file"];

    // Mean and standard deviation for standardization
    Rcpp::List scale = x["scale"];
    NumericVector mean = scale["mean"];
    NumericVector sd   = scale["sd"];

    // Index counters for csv (ci, cj) file and our results matrix (mi, mj)
    int ci = 0, cj = 0, mi = 0, mj = 0;

    // Get original column names (colnames), set up new
    // vector (mcname) for the 'matrix column names' added
    // later as dimnames.
    CharacterVector colnames = x["colnames"];
    CharacterVector mcname(nj);
    mj = 0; mj = 0;
    for (cj = 0; cj < ncol; cj++) {
        if (contains(j, cj + 1)) {
            mcname[mj] = colnames[cj];
            mj++;
        }
    }

    // Setting up resulting matrix; adding dimension names
    NumericMatrix res(ni, nj);
    res.attr("dimnames") = Rcpp::List::create(R_NilValue, mcname);

    // Open file connection
    std::ifstream myfile(file);
    std::string line;
    std::string val;

    // Skipping non-data lines
    if (verbose) Rcout << "[cpp] Skipping unrequired lines\n";
    if (skip > 0) {
        for (int i = 0; i < skip; i++) std::getline(myfile, line, '\n');
    }
    if (header) {
        std::getline(myfile, line, '\n'); // Skipping header
    }
    std::stringstream iss(""); // Must be 'globally' available
                               //
    // Start extracting the data ...
    if (verbose) Rcout << "[cpp] Start extracting required data\n";
    ci = 0; // count from zero
    mi = 0;
    while (std::getline(myfile, line, '\n')) {
        // End search
        if (ci > max(i)) break;
        // Check if current row (index) is needed
        if (!contains(i, ci + 1)) {
            ci++; // Increasing csv line counter
            continue;
        }

        // Else parse the line and extract the required values
        iss.clear(); iss.str(line);
        mj = 0;
        for (cj = 0; cj < ncol; cj++) {
            std::getline(iss, val, sep);
            if (contains(j, cj + 1)) {
                if (!standardize) {
                    res(mi, mj) = std::stod(val);
                } else {
                    res(mi, mj) = (std::stod(val) - mean[cj]) / sd[cj];
                }
                mj++;
            }
        }
        mi++; // Increasing matrix row counter
        ci++; // Increasing csv line counter
    }

    return res;

}


