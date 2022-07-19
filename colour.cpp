#include "colour.h"
#include <cmath>
#include <algorithm>

#define BLACK Colour::from_rgb(0.0, 0.0, 0.0)

template<typename T>
static T clamp(T d, T min, T max) {
  const T t = d < min ? min : d;
  return t > max ? max : t;
}

template<typename T>
static T mix(T x, T y, T a)  // lerp
{
    return x * (1 - a) + y * a;
}

template<typename T>
Colour Colour::from_rgb(T r_, T g_, T b_)
{
    Colour out;
    out.r = static_cast<unsigned int>(r_ * 255.99);
    out.g = static_cast<unsigned int>(g_ * 255.99);
    out.b = static_cast<unsigned int>(b_ * 255.99);
    return out;
}

template<typename T>
Colour Colour::from_brightness(T grey)
{
    Colour out;
    unsigned int greyval = static_cast<unsigned int>(grey * 255.99);
    out.r = greyval;
    out.g = greyval;
    out.b = greyval;
    return out;
}

template<typename T>
Colour Colour::from_hsv(T h, T s, T v)
{
    if (h == -1.0)
        return Colour::from_brightness(0.0);

    T H = h * 360.0;
    T C = s*v;
    T X = C*(1-fabs(fmod(H/60.0, 2)-1));
    T m = v-C;
    T r,g,b;
    if(H >= 0 && H < 60){
        r = C,g = X,b = 0;
    }
    else if(H >= 60 && H < 120){
        r = X,g = C,b = 0;
    }
    else if(H >= 120 && H < 180){
        r = 0,g = C,b = X;
    }
    else if(H >= 180 && H < 240){
        r = 0,g = X,b = C;
    }
    else if(H >= 240 && H < 300){
        r = X,g = 0,b = C;
    }
    else{
        r = C,g = 0,b = X;
    }
    r += m;
    g += m;
    b += m;

    return Colour::from_rgb(r, g, b);
}

