Here is a brief guide on how to use the current version of the pyther code:

1 - Data Preparation

Start by saving your gene expression signature (GES) and interactome object as .tsv files. For the interactome, you should use the `InteractomeToTable` function in the `r-helper.R` script.

2 - Prepare your environment.

Set your path variables and source in the code as follows:

```python
path_str = "PATH_TO_PYTHER"
os.chdir(path_str)

from pyther_classes import *
from pyther_fn import *
```

3 - Load Data

Load your data as follows (presuming they are in the same directory as pyther code, modify your script as necessary):

```python
gesObj = anndata.read_csv('ges.tsv', delimiter = '\t')
intObj = interactome_from_tsv('network.tsv', 'network_name')

```

4 - Run aREA Function

Run aREA. This is equivalent to running Viper with `eset.filter = FALSE` and `method = 'none'`. In future updates, we will include the sample-shuffling null model that is also available in VIPER for group-vs-group comparisons.

```python
nesMat = aREA(gesObj, intObj)

```
