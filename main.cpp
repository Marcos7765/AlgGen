#include "indv.hpp"
#include "Funcao.hpp"
#include <forward_list>
#include <iostream>
#include <algorithm>
#include <cstring>
//#include <bitset>

//template <typename T, std::size_t N> void parr(const T(&a)[N], std::ostream& o = std::cout){o << "{";for (std::size_t i = 0; i < N-1; ++i){o << a[i] << ", ";} o << a[N-1] << "}\n";}

/*
fancy lookup that tries to find the index of a given accumullated aptitude
starting from wherever it stopped before (or from the start in initialization)
consider changing aptitude to accumullated aptitude and sorting both arrays
with it so binary search gets viable
*/

inline unsigned short lookupAcc(long double value, long double* aptidoes, 
    bool reset = false){
    static unsigned short lkindex = 0;
    static long double lkacc = aptidoes[0];
    if (reset){lkindex = 0; lkacc = aptidoes[0];}
    if (value <= lkacc){
        while (value < lkacc - aptidoes[lkindex]){
            lkacc -= aptidoes[lkindex];
            lkindex--;
        }
        return lkindex;
    }

    while (value >= lkacc + aptidoes[lkindex+1]){
        lkindex++;
        lkacc+=aptidoes[lkindex];
    }
    return lkindex +1;
}
void selectIndex(unsigned short* dest, long double* aptidoes, size_t tam){
    long double acc = 0;
    for (int i=0; i < tam; i++){
        acc += aptidoes[i];
    }
    
    setrParam(acc);
    lookupAcc(getrRNG(), aptidoes, true);
    for (int i=0; i < tam; i++){
        dest[i] = lookupAcc(getrRNG(), aptidoes);
    }
}

/*
Using an array with half the size of the population is more memory efficient
than using a forward_list for any rate greater than 1/19 (given a big enough
population).
As the rate is usually big enough for it to stop at every item, it could be par-
ralelized rolling each element against a random number from realGer.
*/
void coIndexes(bool* dest, double rate, size_t tam){    
    unsigned int coIndex = 0;
    setgParam(rate);
    unsigned int nStep = getgRNG(); 
    while (true){
        if (coIndex + nStep >= tam/2){
            break;
        }
        dest[coIndex+nStep] = 1;
        //std::cout << "Crossover em " << coIndex+nStep << ", " << coIndex+nStep+1 \n";
        if (coIndex + nStep +1 >= tam/2){
            break;
        }
        coIndex += nStep+1;
        nStep = getgRNG();
    }
}

void mutIndexes(std::forward_list<unsigned char>* dest, double rate, size_t bitsize,
    size_t tam){
    unsigned int mutIndex = 0;
    setgParam(rate);
    unsigned int nStep = getgRNG(); 
    //int bitot = tam*bitsize;
    int bitot = tam*bitsize;
    while (true){
        if (mutIndex + nStep >= bitot){
            break;
        }
        dest[((mutIndex+nStep)/bitsize)].push_front((mutIndex + nStep)%bitsize);
        //std::cout << "Mutacao ocorreu! " << ((mutIndex+nStep)/bitsize) << " e " << (mutIndex + nStep)%bitsize << "\n";
        if (mutIndex + nStep +1 >= bitot){
            break;
        }
        mutIndex += nStep+1;
        nStep = getgRNG();
    }
}

template<size_t N> void nextGen(indv<N>* dest, indv<N>* orig, unsigned short* select, bool* coIndex,
    std::forward_list<unsigned char>* mutIndex, size_t tam){
    for (int i = 0; i < tam; i+=2){
        if (coIndex[i/2] == true){
            crossOver(orig[select[i]], orig[select[i+1]], dest[i], dest[i+1]);
        } else {
            dest[i] = orig[select[i]];
            dest[i+1] = orig[select[i+1]];
        }
        
        if (!(mutIndex[i].empty())){
        for (unsigned char& pos : mutIndex[i]){
            dest[i].flip((unsigned short) pos);
        }}
        if (!(mutIndex[i+1].empty())){
        for (unsigned char& pos : mutIndex[i+1]){
            dest[i+1].flip((unsigned short) pos);
        }}
    }
}

//std::string printai(std::forward_list<unsigned char>& arg){std::string ph;for (auto it = arg.begin(); it != arg.end(); ++it){ph += std::to_string((unsigned short) (*it)) += ", ";}return ph;}
/*
void indv_print(indv<8>& alvo){
    std::cout << std::bitset<8>( 
        ((unsigned long long)alvo.dna[1] << 4) +
        ((unsigned long long)alvo.dna[0])) << "\n";
}
*/
//inline long double aptFunc(long double x, long double y){
//    return w27(x,y);
//}
inline long double aptFunc(long double x, long double y){
    //return x*x - 3*x + 4;
    return w27(x, y);
}
inline long double normApt(long double apt, long double maxX, long double minX, long double media, double temp){
    //return 1.1*maxX - apt;
    //return (1.1*maxX - apt)*(apt <= media);
    //return (media - 1.1*apt)*(apt <= media); //melhor resultado
    //return (media - 1.07*apt)*(media - 1.07*apt)*(apt <= media);
    //return std::abs((media - 1.1*apt)*(apt <= media));
    //return 1.0/std::exp((media - 1.1*apt));
    //return 1.0/std::exp(std::log10(apt + 500.0));
    //return 1.0/std::exp(std::log10(apt + 500.0));
    //return std::exp(-std::log10(apt + 500.0))*(apt <= media);
    return std::exp(-std::log10(apt + std::abs(minX) + 10e-1)/(26-(temp/4)));//*(apt <= media);
    //return std::exp(std::log10((media - apt +1)))*(apt <= media);
    //return std::exp(-apt)*(apt <= media);
    //return 1/(2+std::exp(apt));
    //return std::exp(-apt/(temp/2))*(apt <= ((media+minX)/2));
    //return std::exp(std::exp(-apt)/(500-temp))*(apt <= ((media+2*minX)/2));
}

//template<size_t N> void iter(indv<N>* source, long double* aptitude, size_t tam, 
//    int agmax, int agmin, double coRate, double mutRate){
template<size_t N> void iter(indv<N>* source, long double* aptitude, size_t tam, 
    int agmaxX, int agminX, int agmaxY, int agminY, double coRate, double mutRate){
    
    static auto selecionados = new unsigned short[tam];
    static auto cruzada = new bool[tam/2];
    static auto mutacoes = new std::forward_list<unsigned char>[tam];
    static auto temp = 0.0;
    //if coIndexes went through all elements, it wouldn't be necessary
    std::fill(cruzada, cruzada+(tam/2), false);
    long double maxApt = -INFINITY;
    long double minApt = INFINITY;
    long double media = 0;
    for (int i =0; i<tam; i++){
        //aptitude[i] = aptFunc(source[i].recVal(agmax, agmin));
        //aptitude[i] = aptFunc(source[i].recValX(agmax, agmin), source[i].recValY(agmax, agmin));
        aptitude[i] = aptFunc(source[i].recValX(agmaxX, agminX), source[i].recValY(agmaxY, agminY));
        media += aptitude[i];
        if (aptitude[i] > maxApt){maxApt = aptitude[i];}
        if (aptitude[i] < minApt){minApt = aptitude[i];}
        mutacoes[i].clear();
    }
    for (int i =0; i<tam; i++){
        aptitude[i] = normApt(aptitude[i], maxApt, minApt, media/tam, temp);
    } 
    temp++;
    selectIndex(selecionados, aptitude, tam);
    coIndexes(cruzada, coRate, tam/2);
    mutIndexes(mutacoes, mutRate, N, tam);
    
    static auto novaGer = new indv<N>[tam];
    nextGen(novaGer, source, selecionados, cruzada, mutacoes, tam);
    //std::cout << "Velha geração:\n";
    //for (int i = 0; i<tam;i++){
    //    indv_print(source[i]);
    //}
    //std::cout << "Nova geração:\n";
    //for (int i = 0; i<tam;i++){
    //    indv_print(novaGer[i]);
    //}
    std::memcpy(source, novaGer, tam*(sizeof(indv<N>)));
}

int main(){
    /* melhores parâmetros
    //const double maxD = 500, minD = -500;
    //const double mutRate = 0.0002, coRate = 0.75;
    //const size_t popSize = 1000;
    //const int numGen = 500;
    */
    //const int bitsize = 32;
    //const double maxD = 500, minD = -500;
    ////const double mutRate = 0.01, coRate = 0.7;
    //const double mutRate = 0.6/(bitsize*2), coRate = 0.9;
    //const size_t popSize = 4000;
    //const int numGen = 50;
    const int bitsize = 32;
    const double maxD = 500, minD = -500;
    //const double mutRate = 0.01, coRate = 0.7;
    const double mutRate = 0.6/(bitsize*2), coRate = 0.9;
    const size_t popSize = 1000;
    const int numGen = 100;

    auto teste = new indv<bitsize>[popSize];
    auto aptidoes = new long double[popSize];
    
    setiParam(teste[0].reprMax());
    for (int i=0; i<popSize;i++){
        teste[i].genCrom();
    }

    double janelaMinX = minD, janelaMaxX = maxD;
    double janelaMinY = minD, janelaMaxY = maxD;
    long double totalApt = 0;
    long double gMinInfo[3] = {INFINITY, 0.0, 0.0};
    long double gMaxInfo[3] = {-INFINITY, 0.0, 0.0};

    for (int j=0; j<popSize;j++){
            totalApt+=aptFunc(teste[j].recValX(janelaMaxX, janelaMinX), \
                teste[j].recValY(janelaMaxY, janelaMinY));
            
            if (aptFunc(teste[j].recValX(janelaMaxX,janelaMinX), teste[j].recValY(janelaMaxY,janelaMinY)) < gMinInfo[0]){
                gMinInfo[2] = teste[j].recValY(janelaMaxY,janelaMinY);
                gMinInfo[1] = teste[j].recValX(janelaMaxX,janelaMinX);
                gMinInfo[0] = aptFunc(gMinInfo[1], gMinInfo[2]);
            } else{
                if (aptFunc(teste[j].recValX(janelaMaxX,janelaMinX), teste[j].recValY(janelaMaxY,janelaMinY)) > gMaxInfo[0]){
                    gMaxInfo[2] = teste[j].recValY(janelaMaxY,janelaMinY);
                    gMaxInfo[1] = teste[j].recValX(janelaMaxX,janelaMinX);
                    gMaxInfo[0] = aptFunc(gMinInfo[1], gMinInfo[2]);
                }
            }
        }
    std::cout << totalApt/popSize <<",\t" << gMinInfo[0] <<",\t" << gMinInfo[1] <<",\t"<< gMinInfo[1] <<"\n";

    for (int i=1; i<=numGen;i++){
        //iter(teste, aptidoes, popSize, maxD, minD, coRate, mutRate);
        iter(teste, aptidoes, popSize, janelaMaxX, janelaMinX, janelaMaxY, janelaMinY, coRate, mutRate);
        //iter(teste, aptidoes, 10, 10, -10, 0.90, 0.02);
        totalApt = 0;
        gMinInfo[0] = INFINITY; gMinInfo[1] = 0.0; gMinInfo[2] = 0.0;
        gMaxInfo[0] = INFINITY; gMaxInfo[1] = 0.0; gMaxInfo[2] = 0.0;
        for (int j=0; j<popSize;j++){//std::cout << aptidoes[j] << "\t";}
            //totalApt+=aptFunc(teste[j].recValX(maxD,minD), teste[j].recValY(maxD,minD));
            totalApt+=aptFunc(teste[j].recValX(janelaMaxX, janelaMinX), \
                teste[j].recValY(janelaMaxY, janelaMinY));
            
            if (aptFunc(teste[j].recValX(janelaMaxX,janelaMinX), teste[j].recValY(janelaMaxY,janelaMinY)) < gMinInfo[0]){
                gMinInfo[2] = teste[j].recValY(janelaMaxY,janelaMinY);
                gMinInfo[1] = teste[j].recValX(janelaMaxX,janelaMinX);
                gMinInfo[0] = aptFunc(gMinInfo[1], gMinInfo[2]);
            } else{
                if (aptFunc(teste[j].recValX(janelaMaxX,janelaMinX), teste[j].recValY(janelaMaxY,janelaMinY)) > gMaxInfo[0]){
                    gMaxInfo[2] = teste[j].recValY(janelaMaxY,janelaMinY);
                    gMaxInfo[1] = teste[j].recValX(janelaMaxX,janelaMinX);
                    gMaxInfo[0] = aptFunc(gMinInfo[1], gMinInfo[2]);
                }
            }
        }
        //int ph = i%10;
        //if (i == 100){
        //    //if (
        //    //    std::pow((janelaMaxX - gMaxInfo[i],2)) > std::std::pow((janelaMaxX - gMaxInfo[i],2))
        //    //)
        //    std::pow(
        //        (janelaMaxX - gMinInfo[1]),2) < std::pow((janelaMinX - gMinInfo[1]),2) ? \
        //            janelaMaxX = (80*gMinInfo[1]+20*maxD)/100  \
        //            : janelaMinX = (80*gMinInfo[1]+20*minD)/100;
        //    std::pow(
        //        (janelaMaxY - gMinInfo[2]),2) < std::pow((janelaMinY - gMinInfo[2]),2) ? \
        //            janelaMaxY = (80*gMinInfo[2]+20*maxD)/100 \
        //            : janelaMinY = (80*gMinInfo[2]+20*minD)/100;
        //    setiParam(teste[0].reprMax());
        //    for (int i=0; i<popSize;i++){
        //        teste[i].genCrom();
        //    }
        //}
        //std::cout << "x_medio " << i <<": "<< totalApt/popSize <<"\t";
        //std::cout << "x_minimo " << i <<": " << gMinInfo[1] <<"\t";
        //std::cout << "valor_minimo " << i <<": " << gMinInfo[0] <<"\n";
        std::cout << totalApt/popSize <<",\t" << gMinInfo[0] <<",\t" << gMinInfo[1] <<",\t"<< gMinInfo[1] <<"\n";
    }
    //std::cout << "valor teste: " << aptFunc(1.57895847e-07,-1.64704659e-07) << "\n";
}
/*
int main(){
    #define tam 10
    #define nbits 8
    indv<nbits> teste[tam];
    double aptidao[tam];
    unsigned short selecionados[tam];
    bool cruzada[tam/2];
    std::fill(cruzada, cruzada+(tam/2), false);
    std::forward_list<unsigned char> mutacoes[tam];
    
    auto selecionades = new unsigned short[tam];
    int agmax = 10000;
    int agmin = 0;
    setiParam(teste[0].reprMax());
    for (int i =0; i<tam; i++){
        teste[i].genCrom();
        aptidao[i] = teste[i].recVal(agmax, agmin);
    }
    //parr(aptidao);

    selectIndex(selecionados, aptidao, tam);
    //parr(selecionados);
    coIndexes(cruzada, 0.5, tam/2);
    parr(cruzada);
    mutIndexes(mutacoes, 0.1, nbits, tam);
    
    for (int i = 0; i < tam; i++){
        std::cout << "Operados: " << printai(mutacoes[i]) << "\n";
    }
    
    indv<nbits> novaGer[tam];
    //for (int i = 0;i<tam;i+=2){
    //    if (cruzada[i/2] == true){
    //        std::cout << "Velha geração:\n";
    //        indv_print(teste[i]);
    //        indv_print(teste[i+1]);
    //        crossOver(teste[i], teste[i+1], novaGer[i], novaGer[i+1]);
    //        std::cout << "Nova geração:\n";
    //        indv_print(novaGer[i]);
    //        indv_print(novaGer[i+1]);
    //    }
    //}
    for (int i = 0; i<tam;i++){
        indv_print(teste[i]);
    }
    nextGen(novaGer, teste, selecionados, cruzada, mutacoes, tam);
    std::cout << "Velha geração:\n";
    std::cout << "Nova geração:\n";
    for (int i = 0; i<tam;i++){
        indv_print(novaGer[i]);
    }

    #undef tam
    return 0;
}
*/

//void exIndv(){
//    indv<8> ckc;
//    
//    setiParam(ckc.reprMax());
//    ckc.genCrom();
//    
//    std::cout << std::bitset<8>( 
//        ((unsigned long long)ckc.dna[1] << 4) +
//        ((unsigned long long)ckc.dna[0])) << "\n";
//    ckc = (unsigned long long) 0b11111111;
//    std::cout << std::bitset<8>( 
//        ((unsigned long long)ckc.dna[1] << 4) +
//        ((unsigned long long)ckc.dna[0])) << "\n";
//
//    std::cout << ckc.recVal(2,1) << "\n";
//}