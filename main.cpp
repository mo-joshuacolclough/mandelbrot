#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <thread>
#include <string>
#include <chrono>
#include <functional>
#include <gmp.h>

#include "complex.h"
#include "complex.cpp"
#include "colour.h"
#include "colour.cpp"

#define CPU_NUM 42

template<typename T>
T julia_set(Complex<T> z, const Complex<T>& c, unsigned int max_iter)
{
    unsigned int i = 0;

    while (i < max_iter && z.magnitude_squared() <= 4.0)
    {
        T rtemp = z.r*z.r - z.i*z.i + c.r;
        z.i = 2 * z.r * z.i + c.i;
        z.r = rtemp;

        i++;
    }
    i--; // Take away 1 - if the loop breaks after the first go then the result should be 0.0;

    return (static_cast<T>(i) / static_cast<T>(max_iter));
}

template<typename T>
inline void complex_square(Complex<T>& z)
{
    T rtemp = z.r*z.r - z.i*z.i;
    z.i = 2 * z.r * z.i;
    z.r = rtemp;
}

template<typename T>
T mandelbrot(Complex<T> in, const unsigned int max_iter, const unsigned int palette_lim, const T& log2)
{
    // Returns a value between 0.0 and 1.0. This is a fraction of how 'deep' the number is in the set. 1.0 = max iter reached
    unsigned int i = 0;

    Complex<T> c = in;  // c = z0
    Complex<T> z = Complex<T>(static_cast<T>(0), static_cast<T>(0));

    while (z.magnitude_squared() <= 4.0 && i < max_iter)
    {
        complex_square(z);
        z += c;
        i++;
    }
    if (i > max_iter-2)
        return 0.0;

    complex_square(z);
    z += c;
    complex_square(z);
    z += c;

    i = (i + 1) % palette_lim;

    T inew = static_cast<T>(i) + 1.0 - (std::log(std::log(z.magnitude())) / log2);
    T depth = inew/static_cast<T>(palette_lim);
    depth = depth - std::floor(depth);

    return depth;
}

template<typename T>
T burning_ship(Complex<T> in, unsigned int max_iter)
{
    // Returns a value between 0.0 and 1.0. This is a fraction of how 'deep' the number is in the set. 1.0 = max iter reached
    unsigned int i = 0;

    Complex<T> c = in;  // c = z0
    Complex<T> z = Complex<T>(static_cast<T>(0), static_cast<T>(0));

    while (z.magnitude_squared() <= 4.0 && i < max_iter)
    {
        z.r = std::abs(z.r);
        z.i = std::abs(z.i);

        T rtemp = z.r*z.r - z.i*z.i + c.r;
        z.i = 2 * z.r * z.i + c.i;
        z.r = rtemp;

        i++;
    }

    T inew = static_cast<T>(i);
    T depth = inew / static_cast<T>(max_iter);

    return depth;
}

template<typename T>
T testfunc(Complex<T> in, unsigned int max_iter, unsigned int palette_lim)
{
    // Returns a value between 0.0 and 1.0. This is a fraction of how 'deep' the number is in the set. 1.0 = max iter reached
    unsigned int i = 0;

    Complex<T> c = in;  // c = z0
    Complex<T> z = Complex<T>(static_cast<T>(0), static_cast<T>(0));

    while (z.magnitude_squared() < 4.0 && i < max_iter)
    {
        z.r = z.r * c.i + c.r * 2.0;
        z.i = z.i * c.r + c.i * 2.0;
        complex_square(z);


        i++;
    }

    T inew = static_cast<T>(i % palette_lim);
    T depth = inew / static_cast<T>(max_iter);

    return depth;
}



template<typename T>
inline Complex<T> julia_scale(unsigned int px, unsigned int py, unsigned int width, unsigned int height, T scale)
{
    Complex<T> scaled = Complex<T>( static_cast<T>(px), static_cast<T>(py) );

    scaled.r /= width;
    scaled.r = (scaled.r * 3.0 - 1.5) * scale;
    scaled.i /= height;
    scaled.i = (scaled.i * 2.0 - 1.0) * scale;

    return scaled;
}

template<typename T>
Complex<T> mandelbrot_scale(unsigned int px, unsigned int py, unsigned int width, unsigned int height, T scale)
{
    Complex<T> scaled = Complex<T>( static_cast<T>(px), static_cast<T>(py) );

    scaled.r /= width;
    scaled.r = (scaled.r * 3.0 - 2.0) * scale;
    scaled.i /= height;
    scaled.i = (scaled.i * 2.0 - 1.0) * scale;

    return scaled;
}

template<typename T>
Complex<T> scale_with_focus_centred(unsigned int px, unsigned int py, unsigned int width, unsigned int height, const Complex<T>& focus, T scalex, T scaley)
{
    Complex<T> scaled = Complex<T>( static_cast<T>(px), static_cast<T>(py) );

    scaled.r /= width;
    scaled.r -= 0.5;
    scaled.r = scaled.r / scalex + focus.r;

    scaled.i /= height;
    scaled.i -= 0.5;
    scaled.i = scaled.i / scaley + focus.i;

    return scaled;
}

template<typename T>
inline void write_next_pixel_bw(std::ostream& imgfile, const T& depth)
{
    Colour c = Colour::from_brightness(depth);
    imgfile << c.r << ' ' << c.g << ' ' << c.b << '\n';
}

/*
template<typename T>
inline void write_depth(std::ostream& outfile, const T& depth)
{
    outfile << depth << std::endl;
}
*/


template<typename ... Args>
std::string string_format( const std::string& format, Args ... args ) // From stack overflow (but CC)
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
//    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

template<typename T>
void runth(std::shared_ptr<bool> busyflag, const unsigned int index, const Complex<T> focus, T sc, const unsigned int img_width,
           const unsigned int img_height, const unsigned int palette_lim, const unsigned int max_iter)
{
    *busyflag = true;

    T log2 = std::log(static_cast<T>(2.0));

    std::ofstream outfile(string_format("/scratch/jcolclou/frames/mandelbrot%d.ppm", index));

    // Writes PPM header regardless of if it is a txt or PPM file - contains useful information for both.
    outfile << "P3\n" << img_width << ' ' << img_height << "\n255\n";

    Complex<long double> z0;
    Colour c;
    long double depth;



    for (unsigned int py = 0; py < img_height; py++)
    {
//        std::cout << "Scanlines remaining: " << img_height-py << '\r' << std::flush;
        for (unsigned int px = 0; px < img_width; px++)
        {
            z0 = scale_with_focus_centred<long double>(px, py, img_width, img_height, focus, sc/3.0, sc/2.0);
            depth = mandelbrot<long double>(z0, max_iter, palette_lim, log2);
            c = Colour::from_brightness<long double>(depth);
            ppm_write_next_pixel(outfile, c);
        }
    }

    outfile.close();

    std::cout << "THREAD " << index << " DONE." << std::endl;
    *busyflag = false;
}

template<size_t N>
struct MandThreadpool
{
    std::shared_ptr<bool> busy[N];
    bool initialised[N];
    std::thread workers[N];

    MandThreadpool()
    {
        for (int i = 0; i < N; i++)
        {
            busy[i] = std::make_shared<bool>(false);
            initialised[i] = false;
        }
    }

    template<typename T>
    void new_job(const unsigned int job_index, const Complex<T> focus, T sc, const unsigned int img_width,
                 const unsigned int img_height, const unsigned int palette_lim, const unsigned int max_iter) // WARNING: Waits for free thread
    {
        std::cout << "THREAD STATUS" << std::endl;
        std::cout << "\tB\t|\tI\t" << std::endl;
        for (size_t th = 0; th < N; th++)
        {
            std::cout << '\t' << *busy[th] << "\t|\t" << initialised[th] << '\t' << std::endl;
        }

        bool found = false;
        size_t thread_to_use = 0;
        while (!found)
        {
            for (size_t th = 0; th < N; th++)
            {
                if (!*busy[th])
                {
                    std::cout << "Found thread: " << th << std::endl;
                    thread_to_use = th;
                    found = true;
                    break;
                }
            }

            if (!found)
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }

        // Make sure to join finished thread if needed
        if (initialised[thread_to_use]) {
            std::cout << "JOINING" << std::endl;
            workers[thread_to_use].join();
        }
        else
            initialised[thread_to_use] = true;

        // Free thread is found :) , deploy the job
        workers[thread_to_use] = std::thread(runth<T>, busy[thread_to_use], job_index, focus,
                                             sc, img_width, img_height, palette_lim, max_iter);

        std::this_thread::sleep_for(std::chrono::milliseconds(10)); // Allow some time for the busy flag to be set
        std::cout << "Thread deployed on " << thread_to_use << "." << std::endl;
    }

    void wait_on_close() // joins all threads at the end
    {
        for (int i = 0; i < N; i++)
            workers[i].join();
    }
};

int main()
{
    const unsigned int img_width = 1920;
    const unsigned int img_height = 1080;

    /* --- Zoom ---
    */
    Complex<long double> focus = Complex<long double>(
        0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781,
        -0.641313061064803174860375015179302066579494952282305259556177543064448574172753690255637023068968116237074056553707214
    );

    const unsigned int max_frames = CPU_NUM * 100;
    unsigned int i = 0;
    long double sc = 1.0;
    long double max_iter;
    long double palette_lim = 100;

    MandThreadpool<CPU_NUM> thpool = MandThreadpool<CPU_NUM>();

    while (i < max_frames)
    {
        max_iter = 50 + static_cast<unsigned int>(std::round(100 * std::sqrt(std::sqrt(sc))));
        std::cout << "MAX ITER: " << max_iter << std::endl;

        thpool.new_job(i, focus, sc, img_width, img_height, static_cast<unsigned int>(std::round(palette_lim)), max_iter);
        std::cout << "New job started. " << i << std::endl << std::endl;
        i++;
        sc *= 1.05;
    }

    thpool.wait_on_close();

    // ----------------------

    // --- Single image ---
    /*
    std::ofstream imgfile("/scratch/jcolclou/test_fractal.ppm");

    const unsigned int max_iter = 100;
    const unsigned int palette_lim = 50;
    Complex<long double> focus = Complex<long double>(0.0, 0.0);
    const long double sc = 1.0;

    imgfile << "P3\n" << img_width << ' ' << img_height << "\n255\n";

    Complex<long double> z0;
    Colour c;

    for (unsigned int py = 0; py < img_height; py++)
    {
//        std::cout << "Scanlines remaining: " << img_height-py << '\r' << std::flush;
        for (unsigned int px = 0; px < img_width; px++)
        {
            z0 = scale_with_focus_centred<long double>(px, py, img_width, img_height, focus, sc/3.0, sc/2.0);
            long double depth = testfunc<long double>(z0, max_iter, palette_lim);
            if (depth < 0.001)
                c = BLACK;
            else
                c = Colour::from_hsv<long double>(depth, 0.9, 1.0);

            ppm_write_next_pixel(imgfile, c);
        }
    }

    imgfile.close();
    */
    // --------------------------------

    return 0;
}
