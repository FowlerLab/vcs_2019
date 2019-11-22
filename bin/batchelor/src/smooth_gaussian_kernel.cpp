#include "batchelor.h"

#include "utils.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"

#include <vector>
#include <set>
#include <stdexcept>
#include <cmath>

/* Perform smoothing with the Gaussian kernel */

SEXP smooth_gaussian_kernel(SEXP vect, SEXP index, SEXP data, SEXP sigma) {
    BEGIN_RCPP
    auto correction_vectors=beachmat::create_numeric_matrix(vect);
    const size_t npairs=correction_vectors->get_nrow();
    const size_t ngenes=correction_vectors->get_ncol();
    const Rcpp::IntegerVector _index(index);
    if (npairs!=_index.size()) { 
        throw std::runtime_error("number of rows in 'vect' should be equal to length of 'index'");
    }

    auto mat=beachmat::create_numeric_matrix(data);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes_for_dist=mat->get_nrow();
    const double s2=check_numeric_scalar(sigma, "sigma");
    
    // Constructing the average vector for each MNN cell.
    std::vector<Rcpp::NumericVector> averages(ncells);
    std::set<int> mnncell;
    {
        std::vector<int> number(ncells);
        Rcpp::NumericVector currow(ngenes);
        for (size_t row=0; row<npairs; ++row) {
            correction_vectors->get_row(row, currow.begin());

            int pair_dex=_index[row];
            auto& target=averages[pair_dex];
            ++(number[pair_dex]);
            auto mIt=mnncell.lower_bound(pair_dex);

            if (mIt!=mnncell.end() && *mIt==pair_dex) {
                auto cIt=currow.begin();
                for (auto& t : target){
                    t += *cIt;
                    ++cIt;
                }

            } else {
                target=Rcpp::clone(currow);
                mnncell.insert(mIt, pair_dex);
            }
        }
    
        for (size_t i : mnncell) {
            auto& target=averages[i];
            const auto num=number[i];
            for (auto& t : target) {
                t/=num;
            }
        }
    }

    // Setting up output constructs.
    Rcpp::NumericMatrix output(ngenes, ncells); // yes, this is 'ngenes' not 'ngenes_for_dist'.
    Rcpp::NumericMatrix exponent(ngenes, ncells);
    std::fill(exponent.begin(), exponent.end(), R_NegInf);

    bool starting_prob=true;
    std::vector<double> distances2(ncells), totalprob(ncells, R_NaReal); 
    beachmat::const_column<beachmat::numeric_matrix> mnn_incoming(mat.get(), false), // no support for sparsity here.
        other_incoming(mat.get(), false);

    // Using distances between cells and MNN-involved cells to smooth the correction vector per cell.
    for (size_t mnn : mnncell) {
        mnn_incoming.fill(mnn);
        auto mnn_iIt=mnn_incoming.get_values();

        for (size_t other=0; other<ncells; ++other) {
            double& curdist2=(distances2[other]=0);
            other_incoming.fill(other);
            auto other_iIt=other_incoming.get_values();
            auto iIt_copy=mnn_iIt;

            for (size_t g=0; g<ngenes_for_dist; ++g) {
                const double tmp=(*iIt_copy  - *other_iIt);
                curdist2+=tmp*tmp;
                ++other_iIt;
                ++iIt_copy;
            }

            // Compute log-probabilities using a Gaussian kernel based on the squared distances.
            // We keep things logged to avoid float underflow, and just ignore the constant bit at the front.
            curdist2/=-s2;
        }
        
        // We sum the probabilities to get the relative MNN density. 
        // This requires some care as the probabilities are still logged at this point.
        double density=0;
        bool starting=true;
        for (const auto& other_mnn : mnncell) {
            if (starting) {
                density=distances2[other_mnn];
                starting=false;
            } else {
                density=R::logspace_add(density, distances2[other_mnn]);
            }
        }

        // Each correction vector is weighted by the Gaussian probability (to account for distance)
        // and density (to avoid being dominated by high-density regions).
        // Summation (and then division, see below) yields smoothed correction vectors.
        const auto& correction = averages[mnn];
        auto oIt=output.begin();
        auto eIt=exponent.begin();

        for (size_t other=0; other<ncells; ++other) {
            const double logmult=distances2[other] - density;
            double& logtotal=totalprob[other];

            if (!starting_prob) {
                logtotal=R::logspace_add(logtotal, logmult);
            } else {
                logtotal=logmult;
            }

            // Again, a **lot** of care is required to avoid underflow!
            // The current sum for the 'i'th entry is determined by 'output[i] * exp(exponent[i])'. 
            // This avoids underflow during summation from very large negative 'exponent'.
            // We try to keep 'exponent' as large as possible to avoid overflow in 'output'.
            for (const auto& corval : correction) {
                if (logmult > *eIt) {
                    *oIt *= std::exp(*eIt - logmult);
                    *oIt += corval;
                    *eIt = logmult;
                } else {
                    *oIt += corval * std::exp(logmult - *eIt);
                }
                ++oIt;
                ++eIt;
            }
        }

        starting_prob=false;
    }

    // Dividing by the total probability, and un-logging.
    for (size_t other=0; other<ncells; ++other) {
        auto curcol=output.column(other);
        auto curexp=exponent.column(other);
        const double logtotal=totalprob[other];

        auto eIt=curexp.begin();
        for (auto& val : curcol) {
            val *= std::exp(*eIt - logtotal);
            ++eIt;
        }
    }

    return output;
    END_RCPP
}
