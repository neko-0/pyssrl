# Post-processing packages for SSRL beam test. 

This pacakge is extended based on the the `collinearw` pacakge.

## Implementing event analysis routine

A custom analysis routine can be build from the `HistMaker` class.
The core function to implement is `plevel_process`, which has the following signature:

```python
def plevel_process(self, process, file_name):
  '''
  process: collinearw.Process object
  file_name: str, path to the ROOT file being used.
  '''
  pass # implement your own routine on how to interact with the Process object.
```

Then you can process the `ConfigMgr` as the following:

```python
from pyssrl import SSRLHistMaker

histmaker = SSRLHistMaker() # your own custom routine
histmaker.process(config) # config is the collinearw.ConfigMgr object
config.save("filled.pkl")
```