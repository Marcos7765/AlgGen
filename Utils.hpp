#ifndef Utils
#define Utils
#include <random>
//#include <cstring>

std::random_device rdRng;
std::seed_seq seeder{rdRng()};
std::mt19937 prRng(seeder);
#define bytemax 15
std::uniform_int_distribution<unsigned long long> intGer(0, bytemax);
#undef bytemax
std::uniform_real_distribution<double> realGer(0.0,1.0);
std::geometric_distribution<unsigned int> geoGer(1.0);

void setrParam(double max){
    std::uniform_real_distribution<double>::param_type ph{0.0, max};
    realGer.param(ph);
}

void setgParam(double rate){
    std::geometric_distribution<unsigned int>::param_type ph{rate};
    geoGer.param(ph);
}

void setiParam(unsigned long long max){
    std::uniform_int_distribution<unsigned long long>::param_type ph{0, max};
    intGer.param(ph);
}

inline double getrRNG(){return realGer(prRng);}

inline double getgRNG(){return geoGer(prRng);}

inline unsigned long long getiRNG(){return intGer(prRng);}
#endif
/*
int main(){
    #define t 10000000
    double teste = 0;
    for (int i = 0; i < t; i++){
        unsigned int val = geoGer(prRng);
        //printf("%d\n", val);
        teste += val;
    }
    printf("Média: %f\n", teste/t);
   return 0; 
}*/