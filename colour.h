#ifndef COLOUR_H
#define COLOUR_H

class Colour
{
public:
    unsigned int r;
    unsigned int g;
    unsigned int b;

    template<typename T>
    static Colour from_rgb(T r_, T g_, T b_);

    template<typename T>
    static Colour from_brightness(T grey);
    template<typename T>
    static Colour from_rgb_float(T r_, T g_, T b_);

    template<typename T>
    static Colour from_hsv(T h, T s, T v);
};


#endif // COLOUR_H
