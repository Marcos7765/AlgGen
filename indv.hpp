#ifndef indvClass
#define indvClass
#include "Utils.hpp"
#include <cstring>

template <size_t Nbits> struct indv{
    unsigned char dnaX[Nbits/__CHAR_BIT__ + (Nbits%__CHAR_BIT__ != 0)];
    unsigned char dnaY[Nbits/__CHAR_BIT__ + (Nbits%__CHAR_BIT__ != 0)];
    double recValX(double max, double min);
    double recValY(double max, double min);
    size_t bytesize();
    size_t reprMax();
    void genCrom();
    void flip(unsigned short pos);
    //void operator=(const unsigned long long& val);
};

template <size_t Nbits> inline size_t indv<Nbits>::bytesize(){
    return Nbits/__CHAR_BIT__ + (Nbits%__CHAR_BIT__ != 0);
}
template <size_t Nbits> inline size_t indv<Nbits>::reprMax(){
    return (std::pow(2,Nbits) -1);
}

//template <size_t Nbits> void indv<Nbits>::operator=(const unsigned long long& val){
//    std::memcpy(&dna, &val, bytesize());
//}

template <size_t Nbits> void indv<Nbits>::genCrom(){
    auto val = getiRNG();
    std::memcpy(&dnaX, &val, bytesize());
    val = getiRNG();
    std::memcpy(&dnaY, &val, bytesize());
    //std::memcpy(&dna, &(getiRNG()), bytesize());
}

template <size_t Nbits> inline double indv<Nbits>::recValX(double max, double min){
    unsigned long long ph = 0;
    std::memcpy(&ph, &dnaX, bytesize());
    //std::cout << ph << "\n";
    return min + (max - min)*((ph/((double)reprMax())));
}
template <size_t Nbits> inline double indv<Nbits>::recValY(double max, double min){
    unsigned long long ph = 0;
    std::memcpy(&ph, &dnaY, bytesize());
    //std::cout << ph << "\n";
    return min + (max - min)*((ph/((double)reprMax())));
}

//#include <iostream>
//#include <bitset>
template <size_t Nbits> void indv<Nbits>::flip(unsigned short pos){
    unsigned long long ph = 0;
    if (pos >= Nbits){
        std::memcpy(&ph, &dnaY, bytesize());
        ph = ph ^ (1 << pos%Nbits);
        std::memcpy(&dnaY, &ph, bytesize());
        return;
    }
    std::memcpy(&ph, &dnaX, bytesize());
    ph = ph ^ (1 << pos%Nbits);
    std::memcpy(&dnaX, &ph, bytesize());
    //std::cout << "orig: " << std::bitset<Nbits>(ph) << "\n";
    //std::cout << "pos:" << pos <<" mask: " << std::bitset<Nbits>((1<<pos)) << "\n";
    //*(this) = ph ^ (1 << pos);
}

//#include <iostream>
//#include <bitset>
template <size_t N> void crossOver(indv<N>& orig1, indv<N>& orig2, indv<N>& dest1, indv<N>& dest2){
    setiParam(N-1);
    auto posX = getiRNG();
    auto posY = getiRNG();
    unsigned long long mask = ((1 << (posX+1))-1);
    unsigned long long ksam = orig1.reprMax() - mask;
    unsigned long long ph1 = 0;
    unsigned long long ph2 = 0;
    unsigned long long val;
    //std::cout << "i = " << pos << "\n";
    //std::cout << "mask = " << std::bitset<N>(mask) << "\n";
    std::memcpy(&ph1, &(orig1.dnaX), orig1.bytesize());
    std::memcpy(&ph2, &(orig2.dnaX), orig1.bytesize());
    
    val = ((ksam & ph1) + (mask & ph2));
    std::memcpy(&(dest1.dnaX), &val, orig1.bytesize());
    val = ((ksam & ph2) + (mask & ph1));
    std::memcpy(&(dest2.dnaX), &val, orig1.bytesize());
    //dest1 = ((ksam & ph1) + (mask & ph2));
    //dest2 = ((ksam & ph2) + (mask & ph1));
    std::memcpy(&ph1, &(orig1.dnaY), orig1.bytesize());
    std::memcpy(&ph2, &(orig2.dnaY), orig1.bytesize());

    mask = ((1 << (posY+1))-1);
    ksam = orig1.reprMax() - mask;

    val = ((ksam & ph1) + (mask & ph2));
    std::memcpy(&(dest1.dnaY), &val, orig1.bytesize());
    val = ((ksam & ph2) + (mask & ph1));
    std::memcpy(&(dest2.dnaY), &val, orig1.bytesize());
}

#endif