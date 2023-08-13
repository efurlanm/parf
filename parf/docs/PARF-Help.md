./parf -h

    PARF (C) 2005 Rudjer Boskovic Institute
    Goran Topic, Tomislav Smuc; algorithm by Leo Breiman and Adele Cutler
    Licensed under GNU GPL 2.0

    Usage: rf [OPTION...]
    -h | --help   show this message
    -t file       file to use as training set
    -a file       file to analyse and classify
    -tv [file]    training set votes output file
    -tc [file]    training set confusion matrix output file
    -av [file]    test set votes output file
    -ac [file]    test set confusion matrix output file
    -ar [file]    test set classification results output file
    -aa [file]    test set ARFF output file
    -ta [file]    train + test set ARFF output file
    -c class      the class attribute, or NEW, or LAST (default)
    -cq [n[%]]    quantity of generated class instances (only with -c NEW)
    -cp category  positive category
    -n trees      the number of trees to grow
    -f n          the fill method: 0=none, 1=rough, 2+=# of passes
    -v n          redo the forest with n most important variables
    -vs n         redo the forest with variables more significant than n
    -p [n[%]]     number of proximate cases to take into consideration
    -st [file]    training set scaling output file
    -sa [file]    test set scaling output file
    -sy [file]    prototype scaling output file
    -s n          number of scaling coordinates (default 2)
    -sd [n]       max allowed divergent iterations in scaling calculation
    -y [file]     prototypes output file
    -ya [file]    detailed prototypes analysis output file
    -yn n         number of prototypes per class (if available)
    -yp n         number of proximates to look at for prototypes
    -m mvar       the number of variables to split on
    -mv mvar      the number of split variables with most important variables
    -xs size[%]   the minimum node size to split
    -xr ratio[%]  the split cutoff ratio limit (0-1)
    -b categories the number of categories that use the fast split algorithm
    -bi iters     the number of fast split iterations per category
    -u(u) var,... comma-separated list of used or unused attributes
    -ri var       row identifier variable
    -r seed       random number generator seed
    -fd [file]    dump the forest as a text
    -fs file      save the forest
    -fl file      load the forest
    -i [file]     variable importances output file
    -ic [file]    case-by-case variable importances output file
    -ii [file]    variable interaction output file
    -if [file]    fast variable importances output file
    -ot [file]    training set outlier measure output file
    -oa [file]    testing set outlier measure output file
    -w wt,...     class weight override
    -g [file]     generate gnuplot graphics script
    -gt type      gnuplot graphics terminal type
    --verbose     print what is done