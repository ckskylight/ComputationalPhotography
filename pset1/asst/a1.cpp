/* -----------------------------------------------------------------
 * File:    a1.cpp
 * Created: 2015-09-15
 * -----------------------------------------------------------------
 *
 * Assignment 01
 *
 * ---------------------------------------------------------------*/


#include "a1.h"
using namespace std;

// Change the brightness of the image
// const Image & means a reference to im will get passed to the function,
// but the compiler won't let you modify it within the function.
// So you will return a new image
Image brightness(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
	Image output(im.width(), im.height(), im.channels());
	// Modify image brightness
    for (int i=0; i < im.number_of_elements(); i++) {
        output(i) = im(i) * factor;
    }
	// return output;
	return output;
}


Image contrast(const Image &im, float factor, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    Image output(im.width(), im.height(), im.channels());
    // Modify image contrast
    for (int i=0; i < im.number_of_elements(); i++) {
        output(i) = factor * (im(i) - midpoint) + midpoint;
    }
    return output;
}


Image color2gray(const Image &im, const std::vector<float> &weights) {
    // --------- HANDOUT  PS01 ------------------------------
    Image output(im.width(), im.height(), 1);
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            float average = 0;
            for (int k=0; k < im.channels(); k++) {
                average += weights[k] * im(i,j,k);
            }
            output(i,j) = average;
        }
    }
	return output;
}


// For this function, we want two outputs, a single channel luminance image
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create the luminance image
    Image luminance = color2gray(im);
    // Create the chrominance image
    Image chrominance(im.width(), im.height(), im.channels());
    for (int k=0; k < im.channels(); k++) {
        for (int i=0; i < im.width(); i++) {
            for (int j=0; j < im.height(); j++) {
                chrominance(i,j,k) = im(i,j,k) / luminance(i,j);
            }
        }
    }

    // Create the output vector as (luminance, chrominance)
	vector<Image> output;
    output.push_back(luminance);
    output.push_back(chrominance);

    return output;
}


// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, float brightF, float contrastF, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    vector<Image> lumichromi = lumiChromi(im);
    Image luminance = lumichromi[0];
    Image chrominance = lumichromi[1];

    // Modify brightness, then contrast of luminance image
    Image brightness_adj = brightness(luminance, brightF);
    Image contrast_adj = contrast(brightness_adj, contrastF, midpoint);

    Image output(im.width(), im.height(), im.channels());
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            for (int k=0; k < im.channels(); k++) {
                output(i,j,k) = chrominance(i,j,k) * contrast_adj(i,j);
            }
        }
    }
    return output;
}


vector<float> matrix_multiply(vector<vector<float> > matrix, vector<float> vec) {
    // Init result vector
    vector<float> result;

    // Make sure matrix and vector height are the same
    // Assumes matrix is square
    if (matrix.size() == vec.size()) {
        for (unsigned i=0; i < matrix.size(); i++) {
            float value = 0;
            for (unsigned j=0; j < matrix.size(); j++) {
                value += matrix[i][j] * vec[j];
            }
            result.push_back(value);
        }
    }

    return result;
}

Image rgb2yuv(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    vector<vector<float> > matrix;
    float row1[3] = {0.299, 0.587, 0.114};
    float row2[3] = {-0.147, -0.289, 0.436};
    float row3[3] = {0.615, -0.515, -0.100};
    matrix.push_back(vector<float>(row1, row1+3));
    matrix.push_back(vector<float>(row2, row2+3));
    matrix.push_back(vector<float>(row3, row3+3));

    // Create output image of appropriate size
    Image output(im.width(), im.height(), im.channels());

    // Change colorspace
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            vector<float> rgb;
            rgb.push_back(im(i,j,0));
            rgb.push_back(im(i,j,1));
            rgb.push_back(im(i,j,2));
            vector<float> new_space = matrix_multiply(matrix, rgb);

            output(i,j,0) = new_space[0];
            output(i,j,1) = new_space[1];
            output(i,j,2) = new_space[2];
        }
    }
    return output;
}


Image yuv2rgb(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    vector<vector<float> > matrix;
    float row1[3] = {1, 0, 1.14};
    float row2[3] = {1, -0.395, -0.581};
    float row3[3] = {1, 2.032, 0};
    matrix.push_back(vector<float>(row1, row1+3));
    matrix.push_back(vector<float>(row2, row2+3));
    matrix.push_back(vector<float>(row3, row3+3));

    // Create output image of appropriate size
    Image output(im.width(), im.height(), im.channels());

    // Change colorspace
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            vector<float> yuv;
            yuv.push_back(im(i,j,0));
            yuv.push_back(im(i,j,1));
            yuv.push_back(im(i,j,2));
            vector<float> new_space = matrix_multiply(matrix, yuv);

            output(i,j,0) = new_space[0];
            output(i,j,1) = new_space[1];
            output(i,j,2) = new_space[2];
        }
    }
    return output;
}


Image saturate(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    Image yuv_output(im.width(), im.height(), im.channels());
    Image yuv = rgb2yuv(im);
    // Saturate image
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            yuv_output(i,j,0) = yuv(i,j,0); // y
            yuv_output(i,j,1) = yuv(i,j,1) * factor; // u
            yuv_output(i,j,2) = yuv(i,j,2) * factor; // v
        }
    }

    Image rgb_output = yuv2rgb(yuv_output);
    // return output;
    return rgb_output; // Change this
}


// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Remember to create the output images and the output vector
    vector<Image> output;

    // Dot position
    int dot_x = im.width()/2;
    int dot_y = im.height()/2;

    // Image 2 (gray)
    Image gray = color2gray(im);

    // Image 1 (inverted colors)
    Image yuv_output(im.width(), im.height(), im.channels());
    Image yuv = rgb2yuv(im);
    // Saturate image
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            yuv_output(i,j,0) = 0.5; // y
            yuv_output(i,j,1) = yuv(i,j,1) * -1; // u
            yuv_output(i,j,2) = yuv(i,j,2) * -1; // v

            // Black dot
            if (i == dot_x && j == dot_y) {
                yuv_output(i,j,0) = 0; // y
                yuv_output(i,j,1) = 0; // u
                yuv_output(i,j,2) = 0; // v
                gray(i,j) = 0; // set the gray image to 0 as well
            }
        }
    }
    Image rgb_output = yuv2rgb(yuv_output);

    // Push to vector
    output.push_back(rgb_output);
    output.push_back(gray);

	return output; //Change this
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Implement automatic white balance by multiplying each channel
    // of the input by a factor such that the three channel of the output image
    // have the same mean value. The mean value of the green channel
    // is taken as reference.

    float red_ch_mean = 0;
    float green_ch_mean = 0; // reference
    float blue_ch_mean = 0;
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            red_ch_mean += im(i,j,0);
            green_ch_mean += im(i,j,1); // take the green channel value
            blue_ch_mean += im(i,j,2);
        }
    }
    // Average
    red_ch_mean /= im.width() * im.height();
    green_ch_mean /= im.width() * im.height();
    blue_ch_mean /= im.width() * im.height();

    float red_factor = green_ch_mean / red_ch_mean;
    float blue_factor = green_ch_mean / blue_ch_mean;

    Image output(im.width(), im.height(), im.channels());
    for (int i=0; i < im.width(); i++) {
        for (int j=0; j < im.height(); j++) {
            // Modify red and blue channels
            output(i,j,0) = im(i,j,0) * red_factor;
            output(i,j,1) = im(i,j,1);
            output(i,j,2) = im(i,j,2) * blue_factor;
        }
    }

    return output;
}
