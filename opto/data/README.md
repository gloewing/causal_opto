# Data
'optoDA_data.csv' is the dataset analyzed with the excursion effect methods. It is in a compressed format.
To uncompress use the following R code:

`R.utils::gunzip("optoDA_data.csv.xz", 
                ext = "xz", FUN = xzfile, 
                remove = FALSE)`
