GaMRed is an algorithm or adaptive filtering of insignificant features in high-throughput data, based on Gaussian mixture decomposition. Here, we provide a MATLAB to run the algorithm.

Main scripts:
- GaMRed.m - Estimating noise components using Gaussian mixture model on summarized data (average or variance of gene expression)
- gmm_init_dp - Compute initial conditions for GMM by using dynamic programming based approach
- filter_type.m - Choose method of data summarization (filter type) based on LR slope coefficient
- run_GaMRed.m - Run GaMRed.m function using sample data.

Citations:
1. Marczyk M, Jaksik R, Polanski A, Polanska J. GaMRed - adaptive filtering of high-throughput biological data. IEEE/ACM Transactions on Computational Biology and Bioinformatics 2018. doi: 10.1109/TCBB.2018.2858825. [Epub ahead of print]
2. Polanski A, Marczyk M, Pietrowska M, Widlak P, Polanska J. Initializing the EM Algorithm for Univariate Gaussian, Multi-Component, Heteroscedastic Mixture Models by Dynamic Programming Partitions. International Journal of Computational Methods. 2017;15(03):1850012.
3. Marczyk M, Jaksik R, Polanski A, Polanska J. Adaptive filtering of microarray gene expression data based on Gaussian mixture decomposition. BMC Bioinformatics. 2013;14:101.