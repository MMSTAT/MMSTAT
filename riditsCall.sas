%let mtitle=%str(GLOMERULONEPHROPATHY SEVERITY);
title1 &mtitle;
%let noGroup=4;
%let codeno=5;
%let one=Group1;
data DataFile;
input severity $ Group1-Group&noGroup  @@;
one=&one;
cards;
none      05  05  07  04
minimum   20  25  16  13
mild      21  13  18  13
moderate  16  14  12  14
severe    10  15  19  28
;
%ridits(DataFile,&noGroup,&codeno,.05,No,&mtitle,&one);

%let mtitle=%str(MONONUCLEAR CELL LEUKEMIA SEVERITY);
title1 &mtitle;
%let noGroup=3;
%let codeno=4;
%let one=sum(of Group1-Group&noGroup);
data DataFile;
input severity $ Group1-Group&noGroup  @@;
one=&one;
cards;
none     39 30 29
mild     04 05 02
moderate 02 05 09
severe   05 10 10
;
%ridits(DataFile,&noGroup,&codeno,.05,No,&mtitle,&one);


%let mtitle=%str(1974 Transportation Study -Population);
title1 &mtitle;
%let noGroup=2;
%let codeno=4;
%let one=sum(of Group1-Group&noGroup);
data DataFile;
input severity $ Group1-Group&noGroup  @@;
one=&one;
cards;
none      357 417
minor     540 330
moderate   53  33
serious    35  17
;
%ridits(DataFile,&noGroup,&codeno,.05,No,&mtitle,&one);


%let mtitle=%str(1974 Transportation Study -Group1);
title1 &mtitle;
%let noGroup=2;
%let codeno=4;
%let one=Group1;
data DataFile;
input severity $ Group1-Group&noGroup  @@;
one=&one;
cards;
none      357 417
minor     540 330
moderate   53  33
serious    35  17
;
%ridits(DataFile,&noGroup,&codeno,.05, No ,&mtitle,&one);
