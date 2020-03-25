# DFT-based watermarking method with Gray Component Replacement masking
The demand for security of shared digital and printed images is increasing year over year. There is a need for the robust watermarking scheme which is capable of offering high detection rates even for very aggressive attacks such as Print-Scan process. Yet the perceptibility of the watermark is usually restricting robustness of the watermarking scheme. In this paper a novel DFT based watermarking method with Gray Component Replacement (GCR) masking is proposed. The watermark is additively embedded in the magnitude coefficients of the image within the Fourier domain and then masked using GCR masking to hide the artefacts introduced by the embedding. Experimental results show that GCR masking is capable of completely hide introduced artefacts and increase perceived quality of the watermarked image using both objective and subjective measures. It also outperforms similar state-of-the-art methods in respect to robustness against image attacks. The method is especially suited for printing since the watermark is hidden in visual part of the electromagnetic spectrum but is detectable in infra-red (IR) part of spectrum using IR or hyper-spectral camera. Finally, proposed watermarking method with GCR masking enables the embedding of strong and robust watermark and avoids the introduction of visible artefacts.

## Dependencies
Following libraries need to be on the system:
1. larmadillo (libarmadillo-dev)
2. llcms (liblcms2-dev)

Following python modules are needed:
1. scikit-image
2. pandas
3. outlier-utils

## Installation

This repository contains large files in bin directory so make sure to have [git-lfs](https://git-lfs.github.com) installed. Optionally you can download tarball of the bin directory separately from [here](http://petox-design.com/grf/wmcode-bin.tar.gz) and replace files in bin directory. 

Run [install.sh](https://github.com/Call1st0/dft-based-watermarking-method-with-GCR-masking/blob/master/install.sh) script to compile code from Cpp directory, install dependencies and create link library *libwmgcr* neaded to run the code. Make sure to have pkg-config available on the system since it is used to check for the existence of the dependencies.
```Console
foo@bar:~$ sudo bash install.sh
```

Create conda environment using [environment.yml](https://github.com/Call1st0/dft-based-watermarking-method-with-GCR-masking/blob/master/environments.yml) file and activate it. If you prefer pip just install three requiered modules mentioned in [dependencies](#dependencies).

```Console
foo@bar:~$ conda env create -f environment.yml
foo@bar:~$ conda activate wmgcr
```

Finally run [demo.py](https://github.com/Call1st0/dft-based-watermarking-method-with-GCR-masking/blob/master/python/demo.py) to get the idea how the method works.

```Console
foo@bar:~$ python demo.py
```
## Uninstall
Run [uninstall.sh](https://github.com/Call1st0/dft-based-watermarking-method-with-GCR-masking/blob/master/uninstall.sh) to remove installed library and optionally remove dependencies.

```Console
foo@bar:~$ sudo bash uninstall.sh
```