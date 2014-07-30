
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP filee8c1c7bd55f( SEXP wtmin, SEXP wtmex, SEXP hist, SEXP stimu, SEXP clampT, SEXP t, SEXP d, SEXP a, SEXP n) ;
}

// definition

SEXP filee8c1c7bd55f( SEXP wtmin, SEXP wtmex, SEXP hist, SEXP stimu, SEXP clampT, SEXP t, SEXP d, SEXP a, SEXP n ){
BEGIN_RCPP


arma::mat wtmatrixinhib = Rcpp::as<arma::mat>(wtmin);
arma::mat wtmatrixexhib = Rcpp::as<arma::mat>(wtmex);
arma::mat history       = Rcpp::as<arma::mat>(hist);
arma::mat stimuli       = Rcpp::as<arma::mat>(stimu);
arma::mat clampTask     = Rcpp::as<arma::mat>(clampT);
arma::colvec tau        = Rcpp::as<arma::colvec>(t);
arma::colvec decay      = Rcpp::as<arma::colvec>(d);
arma::colvec act        = Rcpp::as<arma::colvec>(a);
int N                   = Rcpp::as<int>(n);

arma::vec::iterator act_it_b = act.begin();
arma::vec::iterator act_it_e = act.end();

arma::vec::const_iterator tau_it_b = tau.begin();
arma::vec::const_iterator tau_it_e = tau.end();

arma::vec::const_iterator decay_it_b = decay.begin();
arma::vec::const_iterator decay_it_e = decay.end();

arma::mat::row_iterator hist_it = history.begin_row(0);

double gstr=0.6;
double estr=0.4;
double gamma=0.2;
double beta=0.1;
double actrest = -0.1;

double actmax = 1.0;
double actmin = -1.0;

int CTI = 10;
int CSI = 20;
int ITI = 50;
int cycles = 0;
int trial = 0;

int index = 0;

bool last_resp = false;

for (int t = 0; t < N; t++) {
    
    int stim = std::rand()%4;
    int task = std::rand()%2;
    cycles = 0;
    bool response = false;
    int responsetime = 0;
    while (!response || (cycles < (responsetime + 50))) {
        

        
        arma::colvec::fixed<12> posact;
        for (int i = 0; i < 12; i++) {
            if ((*act_it_b) > 0) 
                posact(i) = ((*act_it_b));
            else
                posact(i) = 0; 
            ++act_it_b;      
        }
        
        arma::colvec excitation = wtmatrixexhib*posact;
        
        arma::colvec inhibition = wtmatrixinhib*posact;
    
        arma::colvec netinput;
                
        if(!response) {
            if (cycles < 10) {
                netinput = gstr*arma::trans(clampTask.row(task)) + beta*excitation + gamma*inhibition;
            }
            if ( (cycles > 9) & (cycles < 20) ) {
                netinput = estr*arma::trans(stimuli.row(stim))+ gstr*arma::trans(clampTask.row(task)) + beta*excitation + gamma*inhibition;
    }
    if ( cycles > 19 ) {
                netinput = estr*arma::trans(stimuli.row(stim)) + beta*excitation + gamma*inhibition; 
           
           
            }        
        }
        else {
            netinput = beta*excitation + gamma*inhibition;        
        } 
        
        arma::vec::const_iterator ni_it_b = netinput.begin();
        arma::vec::const_iterator ni_it_e = netinput.end();
        act_it_b = act.begin();
        tau_it_b = tau.begin();
decay_it_b = decay.begin();
        for (arma::vec::iterator i = act_it_b; i != act_it_e; ++i) {

            double deltact;
            if ((*ni_it_b) > 0 ) { 
                deltact = actmax - (*i);
            }
            else {
                deltact = (*i) - actmin;
            }
            deltact = deltact * (*ni_it_b) * (*tau_it_b);
            deltact = deltact - (*decay_it_b * ((*i) - actrest));
            (*i) = (*i)+deltact;
            if ((*i) > actmax) (*i) = actmax;
            if ((*i) < actmin) (*i) = actmin;
            ++tau_it_b;
     ++decay_it_b;
            ++ni_it_b;
        }
        
        if (act(2) > 0.2 || act(3) > 0.2) {
            response = true;        
        }
        
        if (index > 0) {
            if (response && last_resp==0) {
                responsetime = cycles;            
            }        
        }
             
        for (arma::vec::iterator i = act_it_b; i != act_it_e; ++i) {
            (*hist_it) = (*i);
            ++hist_it;        
        }
        
        (*hist_it) = cycles+1;++hist_it;
        (*hist_it) = task+1;++hist_it;
        (*hist_it) = stim+1;++hist_it;
        (*hist_it) = t+1;++hist_it;
        (*hist_it) = response;++hist_it;
        last_resp = response;
        cycles = cycles + 1;
        index = index + 1;
        if(cycles > 200) break;
    }


}

return wrap(history);

END_RCPP
}



