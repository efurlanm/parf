# Parallel Random Forest Algorithm

The [Random Forests](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm) algorithm is one of the best among the known classification algorithms, able to classify big quantities of data with great accuracy. Also, this algorithm is inherently parallelizable.

Originally, the algorithm was written in the programming language Fortran 77, which is obsolete and does not provide many of the capabilities of modern programming languages; also, the original code is not an example of "clear" programming, so it is very hard to employ in education. Within this project the program is adapted to Fortran 90. In contrast to Fortran 77, Fortran 90 is a structured programming language, legible — to researchers as well as to students.

The creator of the algorithm [[4]](#1), Berkeley professor emeritus [Leo Breiman](https://www.stat.berkeley.edu/users/breiman/), expressed a big interest in this idea in our correspondence. He has confirmed that noone has yet worked on a parallel implementation of his algorithm, and promised his support and help. Leo Breiman is one of the pioneers in the fields of machine learning and data mining, and a co-author of the first significant programs (CART – Classification and Regression Trees) in that field.

The most up-to-date version of PARF's source code can be checked out by SVN (go to the *Source* tab) [[1]](#1) . Snapshots of the SVN repository are created occasionaly, and available from the *Downloads* tab [[2]](#2) for convenience.

To make an executable, a Fortran 90 compiler is required. The currently supported compilers are: [Intel Fortran](http://www.intel.com/cd/software/products/asmo-na/eng/282048.htm) (free on Linux for academic use), [Portland Group Fortran](http://www.pgroup.com/products/workpgi.htm) (commercial) and [GNU g95](http://www.g95.org/) (free) [[3]](#3) .

For an efficient Java implementation of the Random Forest algorithm that integrates into the [Weka](http://www.cs.waikato.ac.nz/ml/weka/) environment, see [FastRandomForest](https://code.google.com/archive/p/fast-random-forest/).

*RF* and *Random Forests* are registered trademarks of Leo Breiman and Adele Cutler.

PARF was developed in the Centre for Informatics and Computing (http://www.irb.hr/en/cir/) and Division of Electronics (http://www.irb.hr/en/str/zel/) of Rudjer Boskovic Institute (http://www.irb.hr/), with the financial support of Ministry of Science, Technology and Sports of Croatia, i-Project 2004-111.

Authors: Goran Topić and Tomislav Šmuc.

---
<span id="1"></span>[1] See: https://github.com/efurlanm/ml/tree/master/parf2016 . Source: https://code.google.com/archive/p/parf/source/default/source . </br>
<span id="2"></span>[2] See: https://github.com/efurlanm/ml/tree/master/parf2008 . Source: https://code.google.com/archive/p/parf/downloads . </br>
<span id="3"></span>[3] Currently, there are many free compilers like GNU Fortran, LLVM Flang, Intel Fortran, Nvidia Fortran, and LFortran. The code compiled correctly on the Intel Fortran 2021.2 . </br>
<span id="4"></span>[4] See: https://www.stat.berkeley.edu/users/breiman/RandomForests/



# References

* https://waikato.github.io/weka-wiki/documentation/#general-documentation
* https://www.cs.waikato.ac.nz/ml/weka/Witten_et_al_2016_appendix.pdf
* https://ufpr.dl.sourceforge.net/project/weka/documentation/3.9.x/WekaManual-3-9-3.pdf
* https://www.cs.waikato.ac.nz/ml/weka/mooc/dataminingwithweka/
