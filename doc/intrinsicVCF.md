## Intrinsic VCF
This documentation describes the steps to create intrinsic VCF for variant selection towards the last step in the workflow. 

The functionality is implemented in `src/makeIntrinsicVCF.sh`.

### 1. Masked alignment
We merge the coordinates of all principal components (in `principal_components/`), then perform masked alignment against the all principal components.
Based on a user-defined read length (>100), we extract regions with length >= 100 and less than the specified read length. 

We merge 
