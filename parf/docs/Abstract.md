# PARF - Parallel RF Algorithm

Source: https://web.archive.org/web/20060925120144/http://www.irb.hr/en/research/projects/it/2004/2004-111/

The random forests algorithm is one of the best among the known classification algorithms, able to classify big quantities of data with great accuracy. Also, this algorithm is inherently parallelisable. Originally, the algorithm was written in the programming language Fortran 77, which is obsolete and does not provide many of the capabilities of modern programming languages; also, the original code is not an example of "clear" programming, so it is very hard to employ in education. Within this project the program is adapted to Fortran 90. In contrast to Fortran 77, Fortran 90 is a structured programming language, legible — to researchers as well as to students. Besides, recently several free implementations were released, which is significant, since one of the reasons for the original choice of Fortran 77 was the lack of free Fortran 90 compilers, and the fear that the software would not be well accepted in the academic community. The creator of the algorithm, [Berkeley](https://web.archive.org/web/20060925120144/http://stat-www.berkeley.edu/users/breiman/) professor emeritus Leo Breiman, expressed a big interest in this idea in our correspondence. He has confirmed that noone has yet worked on a parallel implementation of his algorithm, and promised his support and help. Leo Breiman is one of the pioneers in the fields of machine learning and data mining, and a co-author of the first significant programs (CART – Classification and Regression Trees) in that field.

PARF's source code can be found [here](https://web.archive.org/web/20060925120144/http://www.irb.hr/en/research/projects/it/2004/2004-111/parf.tgz). To make an executable, a Fortran 90 compiler is required. The currently supported compilers are: Intel Fortran, Portland Group Fortran and GNU g95.

RF and Random Forests are registered trademarks of Leo Breiman and Adele Cutler.

PARF was developed in [Centre for Informatics and Computing](https://web.archive.org/web/20060925120144/http://www.irb.hr/hr/cir/) of Ruđer Bošković Institute, with the financial support of Ministry of Science, Technology and Sports, i-Project 2004-111.

Contact: [Goran Topić](https://web.archive.org/web/20060925120144/http://www.irb.hr/en/home/amadan/)
