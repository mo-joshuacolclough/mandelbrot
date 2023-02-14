#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <thread>
#include <string>
#include <chrono>
#include <functional>
#include <iomanip>

//#include <gmp.h>
#include <mpi.h>

#include "complex.h"
#include "complex.cpp"
#include "colour.h"
#include "colour.cpp"

using BaseType = long double;

inline void complex_square(Complex<BaseType>& z) {
  const BaseType temp = z.r*z.r - z.i*z.i;
  z.i = 2 * z.r * z.i;
  z.r = temp;
}

BaseType mandelbrot(Complex<BaseType> in, const unsigned int max_iter, const unsigned int palette_lim, const BaseType& log2) {
  // Returns a value between 0.0 and 1.0. This is a fraction of how 'deep' the number is in the set. 1.0 = max iter reached
  unsigned int i = 0;

  Complex<BaseType> c = in;  // c = z0
  Complex<BaseType> z = Complex<BaseType>(0, 0);

  while (z.magnitude_squared() <= 4.0 && i < max_iter) {
    complex_square(z);
    z += c;
    i++;
  }
  if (i > max_iter-2) {
    return 0.0;
  }

  complex_square(z);
  z += c;
  complex_square(z);
  z += c;

  i = (i + 1) % palette_lim;

  const BaseType inew = static_cast<BaseType>(i) + 1.0 - (std::log(std::log(z.magnitude())) / log2);
  BaseType depth = inew/static_cast<BaseType>(palette_lim);
  depth = depth - std::floor(depth);

  return depth;
}

Complex<BaseType> mandelbrot_scale(unsigned int px, unsigned int py, unsigned int width, unsigned int height, const BaseType scale) {
  Complex<BaseType> scaled(static_cast<BaseType>(px), static_cast<BaseType>(py));

  scaled.r /= width;
  scaled.r = (scaled.r * 3.0 - 2.0) * scale;
  scaled.i /= height;
  scaled.i = (scaled.i * 2.0 - 1.0) * scale;

  return scaled;
}

Complex<BaseType> scale_with_focus_centred(unsigned int px, unsigned int py,
                                           unsigned int width, unsigned int height,
                                           const Complex<BaseType>& focus,
                                           const BaseType scalex, const BaseType scaley) {
  Complex<BaseType> scaled( static_cast<BaseType>(px), static_cast<BaseType>(py) );

  scaled.r /= width;
  scaled.r -= 0.5;
  scaled.r = scaled.r / scalex + focus.r;

  scaled.i /= height;
  scaled.i -= 0.5;
  scaled.i = scaled.i / scaley + focus.i;

  return scaled;
}

inline void write_next_pixel_bw(std::ostream& imgfile, const BaseType& depth) {
  Colour c = Colour::from_brightness(depth);
  imgfile << c.r << ' ' << c.g << ' ' << c.b << '\n';
}

inline void write_next_pixel_viridis(std::ostream& imgfile, const BaseType depth) {
  Colour c = Colour::viridis(depth);
  imgfile << c.r << ' ' << c.g << ' ' << c.b << '\n';
}

inline void write_depth(std::ostream& outfile, const BaseType depth[], const unsigned int length) {
  outfile.write(reinterpret_cast<const char*>(depth), std::streamsize(length * sizeof(BaseType)));
}

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args ) {
  int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  //    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

void runframe(const unsigned int index, const Complex<BaseType> focus,
              BaseType sc, const unsigned int img_width, const unsigned int img_height,
              const unsigned int palette_lim, const unsigned int max_iter) {
  const BaseType log2 = std::log(2.0);

  //    std::ofstream outfile(string_format("/scratch/jcolclou/frames/mandelbrot%d.ppm", index));

  Complex<BaseType> z0;
  const unsigned int total_pix = img_width * img_height;
  BaseType* depthbuffer = new BaseType[total_pix]; // depth for every pixel


#pragma omp parallel for collapse(2)
  for (unsigned int py = 0; py < img_height; py++) {
    for (unsigned int px = 0; px < img_width; px++) {
      z0 = scale_with_focus_centred(px, py, img_width, img_height, focus, sc/3.0, sc/2.0);
      depthbuffer[py * img_width + px] = mandelbrot(z0, max_iter, palette_lim, log2);
    }
  }

  // WRITE BUFFER TO FILE
  // --- txt --
  //    std::ofstream outfile(string_format("/scratch/jcolclou/frames/mandelbrot%d.txt", index), std::ios::binary | std::ios::out);


  // --- ppm ---
  std::ofstream outfile(string_format("/scratch/jcolclou/frames/mandelbrot%d.ppm", index), std::ios::binary | std::ios::out);

  outfile << "P3\n" << img_width << ' ' << img_height << "\n255\n";

  for (unsigned int p = 0; p < total_pix; p++) {
    write_next_pixel_viridis(outfile, depthbuffer[p]);
  }

  // --- txt ---
  // write_depth(outfile, depthbuffer, total_pix);

  outfile.close();

  delete [] depthbuffer;
}


int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  const unsigned int img_width = 1920;
  const unsigned int img_height = 1080;

  /* --- Zoom ---
    */
  Complex<BaseType> focus(
    0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781,
    -0.641313061064803174860375015179302066579494952282305259556177543064448574172753690255637023068968116237074056553707214
  );

  const unsigned int max_frames = 4;
  const BaseType sc0 = 1e10;
  BaseType sc = sc0;
  BaseType palette_lim = 256;
  BaseType max_iter;

  // Split frames to different ranks
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int frames_per_rank = max_frames / world_size;

  int f0 = frames_per_rank * my_rank;
  std::cout << "My frame offset: " << f0 << std::endl;

  for (int i = f0; i < f0 + frames_per_rank; ++i) {
    sc = sc0 * pow(1.01, i);
    max_iter = 50 + static_cast<unsigned int>(std::round(100 * std::sqrt(std::sqrt(sc))));
    std::cout << "MAX ITER: " << max_iter << std::endl;

    runframe(i, focus, sc, img_width, img_height, static_cast<unsigned int>(std::round(palette_lim)), max_iter);
  }

  // ----------------------

  // --- Single image ---
  /*
    std::ofstream imgfile("/scratch/jcolclou/test_fractal.ppm");

    const unsigned int max_iter = 100;
    const unsigned int palette_lim = 50;
    Complex<BaseType> focus = Complex<BaseType>(0.0, 0.0);
    const BaseType sc = 1.0;

    imgfile << "P3\n" << img_width << ' ' << img_height << "\n255\n";

    Complex<BaseType> z0;
    Colour c;

    for (unsigned int py = 0; py < img_height; py++)
    {
//        std::cout << "Scanlines remaining: " << img_height-py << '\r' << std::flush;
        for (unsigned int px = 0; px < img_width; px++)
        {
            z0 = scale_with_focus_centred(px, py, img_width, img_height, focus, sc/3.0, sc/2.0);
            BaseType depth = testfunc(z0, max_iter, palette_lim);
            if (depth < 0.001)
                c = BLACK;
            else
                c = Colour::from_hsv<BaseType>(depth, 0.9, 1.0);

            ppm_write_next_pixel(imgfile, c);
        }
    }

    imgfile.close();
    */
  // --------------------------------


  MPI_Finalize();
  return 0;
}
