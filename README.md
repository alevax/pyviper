# Pyther - beta version

This package enabling network-based protein activity estimation on python. Functions are partly transplanted from R package [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html
).

Here is a brief guide on how to use the current version of the pyther code:

---

### Dependencies

- `scanpy` for single cell pipeline
- `pandas` and `anndata` for data computing and storage. 
- `numpy` and `scipy`  for scientific computation.
- `joblib` for parallel computing
- `tqdm` show progress bar


### Path and Environment

Set your path variables and source in the code as follows:

```python
path_str = "PATH_TO_PYTHER"
os.chdir(path_str)

from pyther_classes import *
from pyther_fn import *
from pyther_narnea import *
```


### Data Preparation

Start by saving your gene expression signature (GES) and interactome object as .tsv files. For the interactome, you should use the `InteractomeToTable` function in the `r-helper.R` script.

taking the example from test.ipynb
```
HMP_16_pruned = pruneRegulon(HMP_16_filtered)

```


### Load data

Load your data as follows (presuming they are in the same directory as pyther code, modify your script as necessary):

```python
gesObj = anndata.read_csv('ges.tsv', delimiter = '\t')
intObj = interactome_from_tsv('network.tsv', 'network_name')

```

### Pyther Functions

-   **pyther**

    ```
    def pyther(
              gex_data,  # anndata obj
              interactome, # list of interactomes
              layer = None,
              njobs = 3,
              eset_filter = True, 
              mvws=1, # an integer or a list of 2 numerical values
              method=[None, "scale", "rank", "mad", "ttest"], 
              verbose= True,
              output_type  = ['anndata', 'ndarray'])

    ```

    python version of metaViper, equivalent to running Viper with `minsize = 0` and `method = 'none'`

    ```
    pyther(testset,[hmp_16, hmp_19])
    ```

- **aREA**
  
    ```
    aREA(
        gex_data, # anndata obj
        interactome, # a single interactome
        layer = None)
    
    ```
  
    Run aREA. This is equivalent to running Viper with `eset.filter = FALSE` and `method = 'none'`. In future updates, we will include the sample-shuffling null model that is also available in VIPER for group-vs-group comparisons.

    ```
    nesMat = aREA(gesObj, intObj)

    ```

- **matrix_narnea**
  
  ```
  matrix_narnea(
                gesObj, # anndata obj
                int_table, # dataframe
                intermediate = False, # when true output is dataframe
                min_targets = 30)
  
  ```

  matrix narnea in still under development. Current version takes a dataframe of regulon network as input `int_table`. Output is a list contains two dataframes: `nes` and `pes` score,  ordered alphabetically. For input and output details, please refer to test_narnea.ipynb

  ```
  test_narnea = matrix_narnea(testset, int_table)
  ```

- **meta_narnea**
  
    ```
    meta_narnea(
                gesObj, # anndata obj
                intList, # list of dataframes
                njobs = 3)
    
    ```

  meta_narnea enables matrix narnea to be performed on a list of networks, with results being intergraded by weighted average. Output is an anndata object of NES score with an extra layer of `.pes` of PES score.


- **slice_concat**
  
  ```
  slice_concat(
              inner_function, 
              gex_data ,
              bins = 10, 
              write_local = True, 
              **kwargs)
  ```

  this function works for large dataset, which allows you to slice the data, run the function and concat the results.In current version, it only support functions with output of dataframe or anndata object. For pyther function , the output type should be 'ndarray'. For detailed example, please refer to test.ipynb

  ```
  slice_concat(pyther,gex_data, bins = 2, interactome = intList, output_type = 'ndarray',verbose = False)
  ```