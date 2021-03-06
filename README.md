# OUT OF DATE

This is an archived and NOT UPDATED repository containing a previous generation of the SMALL-LABS code. Please go to https://github.com/BiteenMatlab/SMALL-LABS to get the updated version


# BGSUBFIT

BGSUBFit is code to do real background subtraction of single molecule imaging movies so that the molecules can be fit and their intensity accurately measured even in the presence of an arbitrarily complex low frequency fluorescent background.

The program can also do fitting without doing background subtraction. See the User Guide for more details.

Written by Benjamin P Isaacoff at the University of Michigan.

## Installation

Download the entire folder and unzip if you downloaded the .zip folder. Move the @TIFFStack directory to your MATLAB home directory. Change the working directory in Matlab to this folder and simply call the functions in the command window as described in the User Guide.

## Usage

See the Quick Start Guide for a quick introduction to using BGSUBFit. 

See the User Guide for the details. Briefly, the function *BGSUBFit* is a wrapper for the other code to perform all of the steps in the correct order. Simply run *BGSUBFit* by specifying the directory containing your .tif stack movies, specify the three required parameters, and any optional parameters, then run it and click to choose the movies you want to fit. Or run the various programs independently.

## Contributing

1. Please inform me before making any changes, then follow the directions below: 
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

Written & first uploaded by BPI on 6/12/16.

Last update 10/20/16

## Credits

All individual programs should have all their individual attributions still in place including authors. 

The development of this code is greatly indebted to the work of David J Rowland (often referred to as DJR in the code), in addition to containing some functions written by him, I’ve borrowed a lot of code snippets from his programs.

The code not written by BPI:

*TiffStack* by DR Muir and BM Kampa  
*saveastiff* by YoonOh Tak  
*bpass* by John C. Crocker and David G. Grier  
*MLEwG* by KI Mortensen, LS Churchman, JA Spudich, H Flyvbjerg  
*gaussFit* by David J Rowland  
*Track_3D2* by David J Rowland   
*hungarian* by Yi Cao

## License

                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

  See LICENSE.txt
