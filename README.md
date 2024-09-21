## 3D deblending for astronomical data cubes

Blending of neighboring HI sources is a common issue in upcoming HI surveys, especially when pushing the utilization of data to the limits of spatial resolution. In this repository, we present an algorithm that effectively separates the fluxes in blended HI detections by leveraging information from optical images as prior knowledge. Furthermore, this algorithm has the potential to be utilized with datasets in other wavelengths (e.g., CO, optical IFU, ...).

### Demo and usage
This repository contains a demo showcasing the deblending of HI in galaxy pairs using the ```watershed``` algorithm from ```scikit-image```. The documentation for the watershed algorithm can be found [here](https://scikit-image.org/docs/stable/api/skimage.segmentation.html#skimage.segmentation.watershed).

To see how the code works, simply run the Jupyter notebook [deblending.ipynb](https://github.com/BetaGem/wallaby-galaxy-pair/deblending.ipynb). For batch processing, an additional script is required. Adjust the free parameters according to your data before running the script.

#### Dependences
- ```numpy```, ```matplotlib```, ```astropy```, ```photutils```, ```scikit-image```, ```astroquery```(optional), ```ipywidgets```(optional)

### Citation
For more comprehensive details on the algorithm, please refer to our paper: [Huang et al. (2024)](). If you find the code useful and plan to utilize it for your own research, we kindly request that you cite our paper as a reference.

Thanks for your interest in our work!
