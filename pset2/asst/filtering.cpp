/* -----------------------------------------------------------------
 * File:    filtering.cpp
 * Created: 2015-09-22
 * -----------------------------------------------------------------
 *
 * Image convolution and filtering
 *
 * ---------------------------------------------------------------*/


#include "filtering.h"
#include <cmath>
#include <cassert>

using namespace std;

Image boxBlur(const Image &im, int k, bool clamp) {
    // --------- HANDOUT  PS02 ------------------------------
    // convolve an image with a box filter of size k by k
    // assuming k is odd
    int k_length = k*k;
    int k_start = -k/2;
    int k_end = k/2;
    float k_value = 1.0/k_length;
    vector<float> kernel;
    for (int i=0; i < k_length; i++) {
        kernel.push_back(k_value);
    }

    Image output(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int z=0; z < im.channels(); z++) {
                float computed_val = 0; // Output of the kernel

                // Kernel loops
                for (int k_i = k_start; k_i <= k_end; k_i++) {
                    for (int k_j = k_start; k_j <= k_end; k_j++) {
                        // k_i and k_j make sense for indexing into the image
                        // but not into the kernel so we make an adjustment to
                        // make it start from 0
                        int k_idx = k*(k_i - k_start) + (k_j - k_start);
                        float kernel_val = kernel[k_idx];
                        computed_val += kernel_val * im.smartAccessor(i+k_i, j+k_j, z, clamp);
                    }
                }
                output(i,j,z) = computed_val;
            }
        }
    }
    return output;
}

Image Filter::convolve(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Write a convolution function for the filter class
    int ki_start = -width/2;
    int ki_end = width/2;
    int kj_start = -height/2;
    int kj_end = height/2;

    Image output(im.width(), im.height(), im.channels());

    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int z=0; z < im.channels(); z++) {
                float computed_val = 0; // Output of the kernel

                // Kernel loops, kernel is row major
                // k_i indexes into rows, k_j into columns
                for (int k_i = ki_start; k_i <= ki_end; k_i++) {
                    for (int k_j = kj_start; k_j <= kj_end; k_j++) {
                        // k_i and k_j make sense for indexing into the image
                        // but not into the kernel so we make an adjustment to
                        // make it start from 0
                        int k_idx = (-k_i - ki_start) + width*(-k_j - kj_start);
                        float kernel_val = kernel[k_idx];
                        computed_val += kernel_val * im.smartAccessor(i+k_i, j+k_j, z, clamp);
                    }
                }
                output(i,j,z) = computed_val;
            }
        }
    }
    return output;
}



Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
    // --------- HANDOUT  PS02 ------------------------------
    // Reimplement the box filter using the filter class.
    // check that your results match those in the previous function "boxBlur"
    int k_length = k*k;
    float k_value = 1.0/k_length;
    vector<float> kernel;
    for (int i=0; i < k_length; i++) {
        kernel.push_back(k_value);
    }

    Filter boxblur(kernel, k, k);
    Image output = boxblur.convolve(im, clamp);

    return output;
}


Image gradientMagnitude(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Uses a Sobel kernel to compute the horizontal and vertical
    // components of the gradient of an image and returns the gradient magnitude.

    vector<float> h_sobel_kernel {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    vector<float> v_sobel_kernel {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    Filter h_sobel(h_sobel_kernel, 3, 3);
    Filter v_sobel(v_sobel_kernel, 3, 3);

    Image h_sobel_im = h_sobel.convolve(im, clamp);
    Image v_sobel_im = v_sobel.convolve(im, clamp);

    Image output(im.width(), im.height(), im.channels());

    // Loop through to find the gradient magnitude
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                float h_val = h_sobel_im(i,j,k);
                float v_val = v_sobel_im(i,j,k);
                output(i,j,k) = sqrt(pow(h_val, 2) + pow(v_val, 2));
            }
        }
    }

    return output;

}

vector<float> gauss1DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // Create a vector containing the normalized values in a 1D Gaussian filter
    // Truncate the gaussian at truncate*sigma.
    vector<float> kernel;
    float sum = 0;
    for (int x=-truncate*sigma; x <= truncate*sigma; x++) {
        float val = exp(-pow(x, 2)/(2*pow(sigma,2)));
        sum += val;
        kernel.push_back(val);
    }

    for (unsigned i=0; i < kernel.size(); i++) {
        kernel[i] /= sum;
    }

    return kernel;
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Gaussian blur across the rows of an image
    vector<float> kernel = gauss1DFilterValues(sigma, truncate);
    int k_width = 1 + 2*ceil(sigma * truncate);
    Filter horizontal_gaussian(kernel, k_width, 1);
    Image output = horizontal_gaussian.convolve(im, clamp);
    return output;
}

Image gaussianBlur_vertical(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Gaussian blur across the rows of an image
    vector<float> kernel = gauss1DFilterValues(sigma, truncate);
    int k_height = 1 + 2*ceil(sigma * truncate);
    Filter horizontal_gaussian(kernel, 1, k_height);
    Image output = horizontal_gaussian.convolve(im, clamp);
    return output;
}

vector<float> gauss2DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // create a vector containing the normalized values in a 2D Gaussian
    // filter. Truncate the gaussian at truncate*sigma.
    vector<float> kernel;
    float sum = 0;
    for (int i=-sigma*truncate; i <= sigma*truncate; i++) {
        for (int j=-sigma*truncate; j <= sigma*truncate; j++) {
            int rad = round(sqrt(pow(i,2) + pow(j,2)));
            float val = exp(-pow(rad, 2)/(2*pow(sigma,2)));
            sum += val;
            kernel.push_back(val);
        }
    }

    for (unsigned i=0; i < kernel.size(); i++) {
        kernel[i] /= sum;
    }

    return kernel;
}


Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    //  Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
    vector<float> kernel = gauss2DFilterValues(sigma, truncate);
    int k_dim = 1 + 2*ceil(sigma * truncate);
    Filter gauss2D(kernel, k_dim, k_dim);
    Image output = gauss2D.convolve(im, clamp);
    return output;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
    Image horizontal = gaussianBlur_horizontal(im, sigma, truncate, clamp);
    Image vertical = gaussianBlur_vertical(horizontal, sigma, truncate, clamp);
    return vertical;
}


Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // sharpen an image
    Image lowpass = gaussianBlur_separable(im, sigma, truncate, clamp);
    Image highpass = im - lowpass;
    Image output = im + highpass*strength;
    return output;
}


Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Denoise an image using the bilateral filter
    return im;
}


Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // 6.865 only
    // Bilaterial Filter an image seperatly for
    // the Y and UV components of an image
    return im;
}




/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(int k){
    // initlize a kxkx1 image of all 0's
    Image impulse(k, k, 1);

    // set the center pixel to have intensity 1
    int center = floor(k/2);
    impulse(center,center,0) = 1.0f;

    return impulse;
}


// ------------- FILTER CLASS -----------------------
Filter::Filter(const vector<float> &fData, int fWidth, int fHeight)
    : kernel(fData), width(fWidth), height(fHeight)
{
        assert(fWidth*fHeight == fData.size());
}


Filter::Filter(int fWidth, int fHeight)
    : kernel(std::vector<float>(fWidth*fHeight,0)), width(fWidth), height(fHeight) {}


Filter::~Filter() {}


const float & Filter::operator()(int x, int y) const {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();

    return kernel[x + y*width];
}


float & Filter::operator()(int x, int y) {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();

    return kernel[x +y*width];
}
// --------- END FILTER CLASS -----------------------
