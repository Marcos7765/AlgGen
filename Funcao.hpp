#include<cmath>

static inline const double a = 500.0;
static inline const double b = 0.1;
static inline const double c = M_PI/2;
inline long double x1(long double x){return 25*x;}
inline long double x2(long double y){return 25*y;}

inline long double rd(long double x, long double y){
    return 1 + (y - x*x)*(y - x*x) + (1-x)*(1-x);
}
inline long double z(long double x, long double y){
    return x*std::sin(std::sqrt(std::abs(x))) - y*std::sin( \
        std::sqrt(std::abs(y)));
}
inline long double w23(long double x, long double y){
    return z(x,y)/rd(x,y);
}
inline long double zsh(long double x, long double y){
    return 0.5 - (( std::pow(std::sin(std::sqrt(x*x + y*y)),2) ) -0.5) / \
    ( std::pow((1 + 0.1* (x*x + y*y) ),2));
}
inline long double F10(long double x, long double y){
    return -a * std::exp(-b * std::sqrt( (x1(x)*x1(x) + x2(y)*x2(y)) /2) ) - \
        std::exp((std::cos(c*x1(x)) + std::cos(c*x2(y))) /2) + std::exp(1);
}
inline long double Fobj(long double x, long double y){
    return F10(x,y)*zsh(x,y);
}
inline long double r(long double x, long double y){
    return 100*(y - x*x)*(y - x*x) + (1 - x)*(1 - x);
}
inline long double w4(long double x, long double y){
    return std::sqrt(r(x,y)*r(x,y) + z(x,y)*z(x,y)) + Fobj(x,y);
}
inline long double w27(long double x, long double y){
    return w4(x,y) + w23(x,y);
}