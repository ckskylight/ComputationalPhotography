#include "a1.h"
#include <iostream>

using namespace std;


// This is a way for you to test your functions.
// We will only grade the contents of a1.cpp and Image.cpp
int main() {

    Image im("./Input/castle_small.png");
    cout << "Image has " << im.number_of_elements() << " elements" << "\n";

    Image brightness_out = brightness(im, 2.0);
    brightness_out.write("./Output/brightness.png");

    Image contrast_out = contrast(im, 2.0, 0.5);
    contrast_out.write("./Output/contrast.png");

    vector<float> weights;
    weights.assign(3, 1./3);
    Image gray_out = color2gray(im, weights);
    gray_out.write("./Output/gray.png");

    weights[2] = 0.0;
    Image gray_out_no_blue = color2gray(im, weights);
    gray_out_no_blue.write("./Output/gray_no_blue.png");

    Image gray_default = color2gray(im);
    gray_default.write("./Output/gray_default.png");

    Image bc_lumi = brightnessContrastLumi(im, 1.5, 1.5, 0.5);
    bc_lumi.write("./Output/bc_lumi.png");

    // Test going to yuv and back to see if we end up with the same image
    Image yuv_image = rgb2yuv(im);
    Image rgb_image = yuv2rgb(yuv_image);
    rgb_image.write("./Output/color_space.png");

    // Test saturate
    Image saturated = saturate(im, 1.5);
    saturated.write("./Output/saturate.png");

    // Test spanish
    vector<Image> spanish_pair = spanish(im);
    spanish_pair[0].write("./Output/spanish_colors.png");
    spanish_pair[1].write("./Output/spanish_gray.png");

    // Test zebra spanish
    Image zebra("./Input/zebra.png");
    vector<Image> zebra_spanish = spanish(zebra);
    zebra_spanish[0].write("./Output/zebra_spanish_colors.png");
    zebra_spanish[1].write("./Output/zebra_spanish_gray.png");

    // Test white balance with flower
    Image flower("./Input/wb.png");
    Image wb = grayworld(flower);
    wb.write("./Output/white_balance_flower.png");
    wb = grayworld(zebra);
    wb.write("./Output/white_balance_zebra.png");

    // std::vector<Image> LC = lumiChromi(im);
    // LC[0].write("./Output/castle_luminance.png");
    // LC[1].write("./Output/castle_chrominance.png");
}
