# dimension-detector
Please see the original research associated with dimension-detector "Detecting the ultra low dimensionality of real
networks" at XXX.

Pedro Almagro, Marián Boguñá, M. Ángeles Serrano

Universitat de Barcelona | Universidad de Sevilla

September 7, 2022

Questions related to code: palmagro@us.es

This repository contains three folders:

1. create_SD (Bash Script)
2. cyclesmap (Bash Script)
3. dimension-detector (Python Library)

The first folder contains the code to obtain $\mathbb{S}^D$ surrogates of a given network, the second folder contains code to calculate chordless cycles in a given netework, the third folder contains code to detect the optimal dimension of a given network.

The workflow to detect the optimal dimension of a given network is:

1. Generate the folder with the surrogates of the network:

    ./create_SD/create_SD.sh *network* *resolution* *n_poll* *wsize* *nrealizations* *maxD*

2. To obtain the feature maps of the surrogates, execute (from the FMC-UB cluster) from bash using the surrogates generated in the previous step:

    ./crete_feats.sh SDnets/network SDfeats/network

3. Execute using Python and the feature maps obtained in the previous step:

    python dimension.py "SDfeats/network" network_features predictors maxk

The following is a brief description of the three blocks (more information can be found in the corresponding sh and py files)

# create_SD.sh 

Script to generate $\mathbb{S}^D$  surrogates of a given network. We assume clustering (T) as proportion of triangles calculated in cyclesmap fortran script. 

The script requires an edgelist of the given network and a file containing a description of its features obtained with cyclesmap fortran script. 

When T is not enough even for D=1, the script generate only a folder D* (it will be interpreted by dimension python script).

Parameters:

- network: name of the network (edgelist file must be located in RealNets folder with edge extension and features file in RealFeats with csv extension)
- resolution: number of surrogates per dimension
- n_poll: number of points used to infer the relation T vs. Beta
- wsize: size of the clustering interval in which to create the surrogates
- nrealizations: number of realizations per random Beta value (default=1)  
- maxD: the script will generate surrogates from D=1 to D=maxD

The resulting surrogates will be placed in RealSD folder. RealFeats folder is used for calculation purposes.

# create_feats.sh 

Script to generate feature maps of a set of surrogates using cyclesmap fortran script.

Parameters:

- surrogates folder: name of the folder containing surrogates (they should be organized by dimension, as obtained with create_SD.sh script)
- results folder: name of the destination folder to store the results

# dimension-detector

Python library to infer the dimension of a network given a set of feature maps of its surrogates (generated with create_feats.sh script). 
The code requires a folder with surrogates organized by dimension, as obtained with create_feats.sh script. 

This main function of this library is dimension(surrogate_set,network_features,predictors,maxk) that recieves a set of feature maps of surrogates of a network, the feature map of the network, a set of predictors (a subset of {'triangles', 'squares','pentagons'}) and a maximum value of k to explore (maxk) and returns the infered dimension for the given network, the value of k and the accuracy for the kNN method. surrogate_set indicates the name of te folder with surrogates organized by dimension, as obtained with create_feats.sh script. network_features indicates the name of a file containing the featuremaps of the network obtained with cyclesmap fortran script.


