//código da primeira tentativa, implementação em classe da população
#include <cmath>
#include "indv.hpp"
#include "Utils.hpp"
#include <forward_list>

#define normCoef 1.1
#define norm

class Pop{
private:
    /* data */
    int lkindex = 0;
    double lkacc = 0.0;
public:
    indv* indvs;
    unsigned int tam;
    double coRate;
    double mutRate;
    double acc = 0.0;
    double maxApt = -INFINITY;
    Pop(unsigned int tam, double coRate = 0.75, double mutRate = 0.02,
        bool random = true);
    ~Pop();
    inline void normalizar();
    indv** select();
    inline indv* lookupAcc(double value);
    std::forward_list<unsigned int>* operatorsIndex();
    void proxGen(indv*** selecao,std::forward_list<unsigned int>** opLists);
    //void loop();
};

Pop::Pop(unsigned int _tam, double _coRate, double _mutRate, bool random){
    indvs = new indv[_tam];
    this->tam = _tam;
    this->coRate = _coRate;
    this->mutRate = _mutRate;
    //indv* ph;
    if (!random) {return;}
    
    for (int i = 0; i<tam; i++){
        indv* ph = indvs +i;
        //indv ph = indvs[i];
        //ph = indvs +i;
        //ph->genCrom();
        genCrom(*ph);
        //ph->calcApt();
        calcApt(*ph);
        if ((*ph).aptidao > maxApt) {maxApt = (*ph).aptidao;}
    }   
}

Pop::~Pop(){
}

inline void Pop::normalizar(){
    for (int i=0; i<tam; i++){
        indvs[i].aptidao = normCoef*maxApt - indvs[i].aptidao;
    }
}

indv** Pop::select(){
    acc = 0;
    for (int i=0; i < tam; i++){
        #ifdef norm
            indvs[i].aptidao = normCoef*maxApt - indvs[i].aptidao;
        #endif
        acc += indvs[i].aptidao;
    }
    
    lkindex = 0;
    lkacc = indvs[0].aptidao;
    //indv* ires[tam];
    indv** ires = new indv* [tam];
    setrParam(acc);
    
    for (int i=0; i < tam; i++){
        ires[i] = lookupAcc(getrRNG());
    }
    return ires;
}

std::forward_list<unsigned int>* Pop::operatorsIndex(){

    #define mutOffset 1
    unsigned int coIndex = 0;
    unsigned int mutIndex = 0;

    auto res = new std::forward_list<unsigned int> [tam+1];

    res[0].assign(tam/2, 0);
    setgParam(coRate);
    unsigned int nStep = getgRNG(); 
    auto coIt = res[0].begin();
    while (true){
        if (coIndex + nStep >= tam/2){
            break;
        }
        std::advance(coIt, nStep);
        *coIt = 1;
        //std::cout << "Crossover em " << coIndex+nStep << ", " << coIndex+nStep+1 \n";
        if (coIndex + nStep +1 >= tam/2){
            break;
        }
        coIndex += nStep+1;
        std::advance(coIt, 1);
        nStep = getgRNG();
    }

    setgParam(mutRate);
    nStep = getgRNG(); 
    int bitot = tam*sbits;
    while (true){
        if (mutIndex + nStep >= bitot){
            break;
        }
        res[mutOffset+((mutIndex+nStep)/sbits)].push_front((mutIndex + nStep)%sbits);
        //std::cout << "Mutacao ocorreu! " << ((mutIndex+nStep)/sbits) << " e " << (mutIndex + nStep)%sbits << "\n";
        if (mutIndex + nStep +1 >= bitot){
            break;
        }
        mutIndex += nStep+1;
        nStep = getgRNG();
    }

    return res;
    #undef coOffset
}

inline indv* Pop::lookupAcc(double value){
    if (value <= lkacc){
        while (value < lkacc - indvs[lkindex].aptidao){
            lkacc -= indvs[lkindex].aptidao;
            lkindex--;
        }
        return indvs +lkindex;
    }

    while (value >= lkacc + indvs[lkindex+1].aptidao){
        lkindex++;
        lkacc+=indvs[lkindex].aptidao;
    }
    return indvs +lkindex +1;
}

static void crossOver(indv* orig1, indv* orig2, indv* dest1, indv* dest2){
    setiParam(sbits-1);
    auto pos = getiRNG();
    unsigned long long mask = ((1 << (pos+1))-1);
    unsigned long long ksam = reprMax - mask;
    dest1->cromossomo = dna((ksam & orig1->cromossomo.to_ullong()) + (mask & orig1->cromossomo.to_ullong()));
    dest2->cromossomo = dna((ksam & orig2->cromossomo.to_ullong()) + (mask & orig2->cromossomo.to_ullong()));
}

void Pop::proxGen(indv*** _selecao, std::forward_list<unsigned int>** opLists){
    indv* newGen = new indv[tam];
    maxApt = -INFINITY;
    auto coIt = (*opLists)[0].begin();
    for (int i = 0; i < tam; i+=2){
        auto pS = (*_selecao) + i;
        auto pG = newGen + i;
        if (*coIt == 1){
            crossOver(pS[0],pS[1], pG, pG+1);
        } else {
            *(pG) = *((*_selecao)[i]);
            *(pG+1) = *((*_selecao)[i+1]);
        }
        ++coIt;
        
        if (!((*opLists)[i+1].empty())){
        for (unsigned int& pos : (*opLists)[i+1]){
            newGen[i].cromossomo.flip(pos);
        }}
        if (!((*opLists)[i+2].empty())){
        for (unsigned int& pos : (*opLists)[i+2]){
            newGen[i+1].cromossomo.flip(pos);
        }}
        //newGen[i].calcApt();
        calcApt(newGen[i]);
        //newGen[i+1].calcApt();
        calcApt(newGen[i+1]);
        if (newGen[i].aptidao > maxApt) {maxApt = newGen[i].aptidao;}
        if (newGen[i+1].aptidao > maxApt) {maxApt = newGen[i+1].aptidao;}
    }
    //delete [] indvs;
    delete[] newGen;
    delete[] (*_selecao); //memory leak???
    for (int i = 0; i<tam+1; i++){
    (*opLists)[i].~forward_list();
    }
    delete[] (*opLists);
    //indvs = newGen;
}
