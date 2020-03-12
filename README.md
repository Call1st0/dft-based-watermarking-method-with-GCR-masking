# DFT-based watermarking method with Gray Component Replacement masking
The demand for security of shared digital and printed images is increasing year over year. There is a need for the robust watermarking scheme which is capable of offering high detection rates even for very aggressive attacks such as Print-Scan process. Yet the perceptibility of the watermark is usually restricting robustness of the watermarking scheme. In this paper a novel DFT based watermarking method with Gray Component Replacement (GCR) masking is proposed. The watermark is additively embedded in the magnitude coefficients of the image within the Fourier domain and then masked using GCR masking to hide the artefacts introduced by the embedding. Experimental results show that GCR masking is capable of completely hide introduced artefacts and increase perceived quality of the watermarked image using both objective and subjective measures. It also outperforms similar state-of-the-art methods in respect to robustness against image attacks. The method is especially suited for printing since the watermark is hidden in visual part of the electromagnetic spectrum but is detectable in infra-red (IR) part of spectrum using IR or hyper-spectral camera. Finally, proposed watermarking method with GCR masking enables the embedding of strong and robust watermark and avoids the introduction of visible artefacts.

## Installation

Run install script to compile code from Cpp directory, install dependencies and create link library *libwmgcr* neaded to run the code. Make sure to have pkg-config available on the system since it used to check for the existence of the libraries neaded.

```Console
foo@bar:~$ sudo bash install.sh
```

Create conda environment using [environment.yml](https://github.com/Call1st0/dft-based-watermarking-method/blob/master/environment.yml) file and activate it.

```Console
foo@bar:~$ conda env create -f environment.yml
foo@bar:~$ conda activate wmgcr
```

Finally run [demo.py](https://github.com/Call1st0/dft-based-watermarking-method/blob/master/python/demo.py) to get the idea how the method works.

```Console
foo@bar:~$ python demo.py
```
## Dependencies
Following libraries need to be on the system:
1. blas
2. lapack
3. armadillo
4. lcms

Following python modules are needed:
1. scikit-image
2. pandas
3. outlier-utils
