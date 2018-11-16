//Model.h

//     Copyright 2008 Marylyn Ritchie

//     This file is part of MDR.
// 
//     MDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     MDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with MDR.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>

using namespace std;

namespace mdr{

struct model_stats{
	int classhigh;
	int classlow;
	int misclasshigh;
	int misclasslow;
	int correctclass;
	int incorrectclass;
	int totaltiecells;
	int totalnotrecognized;
	int totalmissing;
	float error;
	float balanced_error;
	vector<int> affected;
	vector<int> unaffected;
};

class Model{

  public:
    
    model_stats training, testing;
    vector<short> cell_status;
    vector<unsigned int> combination;
    
    void reset_stats();
    void resize_cell_status_vector(unsigned int vector_size);
    
    
    void set_train_stats(model_stats& stats){training=stats;}
    void set_test_stats(model_stats& stats){testing=stats;}
    model_stats& get_train_stats(){return training;}
    model_stats& get_test_stats(){return testing;}
       
    float get_combined_fitness(){return combined_fitness;}
    float get_trainerr(){return trainerror;}
    float get_predicterror(){return predicterror;}
    float get_predictavg(){return predictavg;}
    float get_balpredavg() const{return balpredavg;}
    float get_balpredicterror(){return balpredicterror;}
    float get_totalnotclassified(){return totalnotclassified;}
    float get_threshold(){return threshold;}
    float get_total(){return total;}
    float get_totalmissdata(){return totalmissdata;}
    float get_missingPredValue(){return missingPredValue;}
    float get_pvalue(){return pvalue;}
    int get_cvc() const{return cvc;}
    float get_interact_llr(){return interact_llr;}
    float get_lr_pvalue(){return lr_pvalue;}
    
    void set_cvc(int cv_total){cvc = cv_total;}
    void set_combined_fitness(float cb){ combined_fitness=cb;}
    void set_trainerr(float te){ trainerror=te;}
    void set_predicterror(float pe){ predicterror=pe;}
    void set_predictavg(float peavg){ predictavg=peavg;}
    void set_balpredavg(float bpavg){ balpredavg=bpavg;}
    void set_balpredicterror(float bpe){ balpredicterror=bpe;}
    void set_totalnotclassified(float tnc){ totalnotclassified=tnc;}
    void set_threshold(float thresh){ threshold=thresh;}
    void set_total(float tot){ total=tot;}
    void set_totalmissdata(float totmiss){ totalmissdata=totmiss;}
    void set_missingPredValue(float missPval){ missingPredValue=missPval;}
    void set_pvalue(float pval){ pvalue=pval;}
    void set_interact_llr(float llr){interact_llr = llr;}
    void set_lr_pvalue(float pval){ lr_pvalue=pval;}
        
  private:
    float combined_fitness, trainerror, predicterror, predictavg, balpredavg, 
      balpredicterror, totalnotclassified, threshold, total, totalmissdata,
      missingPredValue, pvalue, interact_llr, lr_pvalue;
    int cvc;
    
};

}

#endif

