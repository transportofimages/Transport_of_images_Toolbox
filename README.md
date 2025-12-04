# Transport of Images

The **Transport of images Toolbox** is a Matlab package for computing all zeros of harmonic mappings by continuation.

---

## Installation

To clone the **Transport of images Toolbox** repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/transportofimages/Transport_of_images_Toolbox.git
```
To use our Toolbox in Matlab, you will need to add the **Transport_of_images_Toolbox** directory to the Matlab path.

---

## Basic usage

To compute all zeros of the harmonic mapping

$$f(z) = z^2 + \overline{\left(\frac{2z+1}{z^2 + z}\right)} + 2 \log |z|$$

type:

```matlab
fun = harmonicRat([1 0 0], [1], [2 1], [1 1 0], [2], [0]);
zer = tiroots(fun)
res = max(abs(fun.f(zer)))
```
Output:

```matlab
zer =
  -0.8775 - 0.9278i
  -0.8775 + 0.9278i
  -1.4569 + 0.0000i
  -0.5872 - 0.0000i
res =
   8.8818e-16
```

For more examples we refer to the m-files ex_*.m.


## Citing Transport of Images

If you find the **Transport of images Toolbox** useful in your work, we kindly request that you cite the [following preprint](https://arxiv.org/abs/2011.00079):

```latex
@article{SeteZur2022,
	Author = {S\`{e}te, Olivier \and Zur, Jan},
	 Title = {{The transport of images method: computing all zeros of harmonic mappings by continuation}},
   Journal = {IMA J. Numer. Anal.},
  fjournal = {IMA Journal of Numerical Analysis},
    volume = {42},
      year = {2022},
    number = {3},
     pages = {2403--2428},
      ISSN = {0272-4979},
       DOI = {10.1093/imanum/drab040},
       URL = {https://doi.org/10.1093/imanum/drab040},
}
```
