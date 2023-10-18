%macro compGroup(compGroup);
/* TIE Group and Fij Statistic Calculations */
/* Fleiss, Design and Analysis of Clinical Experiments
           page 77 (3.37) */
data fij; set ridit;
keep dum tie ndot fnum1 fdenom1 sum1-sum&noGroup
     ndot MnRidit1-MnRidit&noGroup;
tie=Group1+Group&compGroup;
ndot=sum1+sum&compGroup;
fnum1=tie*(tie-1)*(tie+1);
fdenom1=ndot*(ndot-1)*(ndot+1);
%mmeans2(fij,fnum1,sum);
data fij; merge fij meansout; by dum; drop sum;
fnum=sum;
%mmeans2(fij,fdenom1,max);
data fij; merge fij meansout; by dum; drop max;
fdenom=max;
Fij=1-(fnum/fdenom);

data scheffe; set fij; format Chi 10.5;
if _N_ > 1 then delete;
/* Fleiss, Design and Analysis of Clinical Experiments
           page 82 (3.43) */
above=12 * sum1 * sum&compGroup * (MnRidit1-MnRidit&compGroup)**2;
below = ( sum1 + sum&compGroup + 1) * Fij;
Chi=above/below;
df=&noGroup-1;
p_value=1-probchi(Chi,df);

/* OUTPUT of Scheffe ChiSquare */
title;
proc print data=scheffe split='*' noobs; var Chi df p_value;
label Chi="Scheffe`*ChiSquare" df='Degrees*of*Freedom'
p_value='p_value' ;
title "Group1 vs Group&compGroup"; run;
title;
%mend compGroup;


%macro dosums(codeno);
%local i ii stop;
data td; set td; keep sum1-sum&codeno;
sum1=0;
sum2=sum(col1);
%let stop=&codeno+1;
%do i=3 %to &stop;
   %do ii=&i-2 %to &i-2; %end;
   sum&i=sum(of col1-col&ii);
%end;
%mend dosums;


%macro equalmns(noGroup);
/* TIE Group and F Statistic Calculations */
/* Fleiss, Design and Analysis of Clinical Experiments
      page 77 (3.37) and page 82 (3.42) */
data f; set ridit;
keep dum tie ndot fnum1 fdenom1 sum1-sum&noGroup
     ndot MnRidit1-MnRidit&noGroup PopulationRidit tot;
tie=sum(of Group1-Group&noGroup);
fnum1=tie*(tie-1)*(tie+1);
fdenom1=ndot*(ndot-1)*(ndot+1);
%do j=1 %to &noGroup;
comp&j = sum&j * (MnRidit&j-PopulationRidit)**2;
%end;
tot=sum(of comp1-comp&noGroup);

%mmeans2(f,fnum1,sum);
data f; merge f meansout; by dum; drop sum;
fnum=sum;

%mmeans2(f,fdenom1,max);
data f; merge f meansout; by dum; drop max;
fdenom=max;
F=1-(fnum/fdenom);

data equalmns; set f; format Chi 10.5;
if _N_ > 1 then delete;
/* Fleiss, Design and Analysis of Clinical Experiments page 82 (3.42)
*/
above=12 * ndot * tot;
below = ( ndot + 1) * F;
Chi=above/below;
df=&noGroup-1;
p_value=1-probchi(Chi,df);
run; quit;

/* OUTPUT of Scheffe` Test of Equal Means ChiSquare */
proc print data= equalmns split='*' noobs; var Chi df p_value;
   label Chi="ChiSquare" df='Degrees*of*Freedom' p_value='p_value' ;
   title1 "Scheffe` Analysis";
   title2 "Test of Equal Mean Ridits"; run; quit; title;
%mend equalmns;


%macro grafit(noGroup,alpha);
data grafit; set interval; keep MnRiditL MnRidit MnRiditU group;
%local i;
%do i=1 %to &noGroup;
/*  MnRiditL=usualL&i; MnRiditL=SchefL&i; */
/*  MnRiditU=usualU&i; MnRiditU=SchefU&i; */
    MnRiditL=BonL&i;
    MnRidit=MnRidit&i;
    MnRiditU=BonU&i;
    Group=&i;
    output;
%end;
data grafit; set grafit;
   confid1=1-&alpha;
   confid2=put(confid1,3.2);
   confid3=substr(confid2,2);
   call symput('conlim',trim(confid3));
   run;
proc print data=grafit noobs;
  title "grafit";  run cancel;

title1 h=1.65 f='Arial' "&mtitle";                
proc sgplot data=grafit;
series x=group y=MnRiditL / markers;
series x=group y=MnRidit / markers;
series x=group y=MnRiditU / markers;
xaxis grid;
yaxis grid;
run; quit;
%mend grafit;

%macro group;
%local j;
%do j=1 %to &noGroup;
%mmeans2(ridit,Group&j,sum);
   data ridit; merge ridit meansout; by dum;
   data ridit; set ridit; drop sum;
   sum&j=sum;
   product&j=Group&j*ridit;
%mmeans2(ridit,product&j,sum);
   data ridit; merge ridit meansout; by dum;
   data ridit; set ridit; drop sum; 
   sumprod&j=sum;
   MnRidit&j=sumprod&j/sum&j;
   Std&j=1/(2*sqrt(3*sum&j));
   wtMnRidit&j=sum&j*MnRidit&j;
%end;
data ridit; set ridit;
  ndot=sum(of sum1-sum&noGroup);
  PopulationRidit=sum(of wtMnRidit1-wtMnRidit&noGroup)/ndot;
%mend group;


%macro interval(noGroup,alpha);
data interval; set ridit;
%local i;
quanparm=1-(&alpha/2);
quantile=probit(quanparm);
%do i=1 %to &noGroup;
   diff&i   = 1-MnRidit&i;
   sqval1&i = (MnRidit&i * diff&i) / ( sum&i-1 );
   sqval&i  = sqrt(sqval1&i);
   usual&i  = quantile * sqval&i;
   usualL&i = MnRidit&i - usual&i;
   usualU&i = MnRidit&i + usual&i;
   roughL&i = MnRidit&i - 1 / sqrt(3*sum&i);
   roughU&i = MnRidit&i + 1 / sqrt(3*sum&i);
   odds&i   = ( 0.5 + ( MnRidit&i-MnRidit1 ) ) /
        ( 1 - ( 0.5 + (MnRidit&i-MnRidit1) ) );
%end;

data interval; set interval; if _n_ > 1 then delete;
keep MnRidit1-MnRidit&noGroup diff1-diff&noGroup
     usualL1-usualL&noGroup usualU1-usualU&noGroup
     roughL1-roughL&noGroup roughU1-roughU&noGroup
     odds1-odds&noGroup quantile sum1-sum&noGroup;

data interval; set interval; drop quanparm;
%local i;
%do i=1 %to &noGroup;
   num&i = sqrt(sum&i+sum1);
   den&i = 2*sqrt(3*sum&i*sum1);
   se&i  =  num&i/den&i;
   paramter = 1- (&alpha / (2* (&noGroup-1) ) );
   B = probit(paramter);
   paramter =  &noGroup-1;
   quanparm = 1- (&alpha/2);
   S = cinv(quanparm, paramter);
   BonL&i   = MnRidit&i-MnRidit1+.5 - B*se&i;
   BonU&i   = MnRidit&i-MnRidit1+.5 + B*se&i;
   SchefL&i = MnRidit&i-MnRidit1+.5 - S*se&i;
   SchefU&i = MnRidit&i-MnRidit1+.5 + S*se&i;
%end;

data interval; set interval;
   %global conlim;
   confid1=1-&alpha;
   confid2=put(confid1,3.2);
   confid3=substr(confid2,2);
   call symput('conlim',trim(confid3));
   run;

data interval; set interval;
 file print;
    put "         "
        "    &conlim% Confidence Intervals on Individual Mean Ridits"
/;
    put "         " "Group" "   Mean Ridit"
                     "        ROUGH     " "      USUAL"
        "      ODDS" / ;
    %do i= 1 %to &noGroup;
    ii=&i;
    put "         " ii 5.0
        MnRidit&i   10.3
        roughL&i 09.3 "," roughU&i 06.3
        usualL&i 09.3 "," usualU&i 06.3
        odds&i   08.2 ":1";
    %end;
 run;
%mend interval;


%macro mmeans2(dsname,varlst,stat);
proc datasets library=work nolist; delete meansout; run cancel;
proc delete data=meansout; run; quit;
proc means data=&dsname noprint;
  var &varlst;
  output out=meansout
     n=n nmiss=nmiss mean=mean Std=Std min=min max=max range=range
     sum=sum var=var uss=uss css=css Stderr=Stderr cv=cv
/*   skewness=skewness kurtosis=kurtosis sumwgt=sumwgt */
     t=t prt=prt;
run;
data meansout; set meansout;
dum=1;
keep &stat dum; run;
%mend mmeans2;


%macro myprints;
options formdlim='';
title1;
data _null_; file print; put _page_; run;
title1 &mtitle;
proc print data=ridit noobs;
var severity Group1-Group&noGroup one two three four ridit
    product1-product&noGroup sum1-sum&noGroup sumprod1-sumprod&noGroup
    MnRidit1-MnRidit&noGroup Std1-Std&noGroup wtMnRidit1-wtMnRidit&noGroup ndot
    PopulationRidit;
    title2 "Calculations -&noGroup Groups"; run;
options formdlim='';
title1 &mtitle;
proc print data=fij; title2 'fij'; run;
proc print data=scheffe noobs; title2 "scheffe`"; run ;

options formdlim='';
title1 &mtitle;
proc print data=f noobs;        title2 'f'       ; run ;
proc print data=equalmns noobs; title2 "equalmns"; run ;
proc print data=interval noobs label;
var  MnRidit1-MnRidit&noGroup diff1-diff&noGroup
     usualL1-usualL&noGroup usualU1-usualU&noGroup
     roughL1-roughL&noGroup roughU1-roughU&noGroup
     odds1-odds&noGroup quantile sum1-sum&noGroup
     num2-num&noGroup den2-den&noGroup se2-se&noGroup;
title2 "interval"; run;
%mend myprints;


%macro ridits(DataFile,noGroup,codeno,alpha,diagnose,mtitle,one);
/***************************************************************************
*  RIDIT ANALYSIS written  July 22,1995 and                                *
*                 modified 1996,1997,1998,2005                             *
*                                                                          *
*  (C) Copyright Mary A. Marion, Jan 12, 1998                              *
*                                                                          *
*  Inputs:                                                                 *
*                                                                          *
*  NOGroup  = Number of groups                                             *
*  CODENO = No of severity codes (levels)                                  *
*  ALPHA  = Significance Level such as .05                                 *
*  DIAGNOSE = Diagnostic print indicator (Yes or No)                       *
*  MTITLE = title of the experiment                                        *
*  DSNAME = Input data matrix of severity codes by group                   *
*  ONE    = reference population                                           *
*                                                                          *
*  Constraints:                                                            *
*                                                                          *
*  Scheffe'-type comparisons between groups always compare to Group1.      *
*  Thus always enter the control group as Group1 when not combining        *
*  across all groups to form the reference population.                     *                                            *
*  all groups to form the reference population                             *
*                                                                          *
*  Macros called:  compGroup, dosums, equalmns, group, interval,           *
*                  mmeans2, myprints                                       *
*                                                                          *
************************************************************************* */
data timetrak;
time1=time();

data DataFile; set &DataFile;
one=&one;
two=one/2;

/* THREE (Column 3) Computation */
proc transpose data=DataFile out=td; var one; run;
%dosums(&codeno);
proc transpose data=td out=td2; var sum1-sum&codeno; run;
data ridit; merge td2 DataFile; rename col1=three;
keep severity one two col1 dum Group1-Group&noGroup;
dum=1;
run;

/* RIDIT Calculations */
%mmeans2(ridit,one,sum);
   data ridit; merge ridit meansout; by dum;
   data ridit; set ridit; drop sum;
   sum0=sum;
four=two+three;
ridit=four/sum0;

/* GROUP CALCULATIONS */
%group;

/* OUTPUT of table of dose group X severity levels +
   ridits for the severity categories                */
data _null_; file print; put _page_; run;
options formdlim='';
options nonumber;
proc print data=ridit noobs label;
   label severity='Severity' ridit='Ridit';
   var severity Group1-Group&noGroup /* one */ ridit;
   sum Group1-Group&noGroup /* one */;
   title &mtitle; run;

/* OUTPUT of MnRidits, Population Mean Ridit and
          standard errors of MnRidits */
data riditout; set ridit; if _n_ > 1 then delete;
format Std1-Std&noGroup 7.5;
options formdlim=' ';
title;
proc print data=riditout noobs label;
  var MnRidit1-MnRidit&noGroup PopulationRidit;
  %do i=1 %to &noGroup; label MnRidit&i="Mean Ridit&i"; %end;;
  run;
proc print data=riditout noobs;
  var Std1-Std&noGroup;
  run;

/* CONFIDENCE INTERVALS on the RIDIT MEANS */
%interval(&noGroup,&alpha);

/* TESTING the HYPOTHESIS of EQUAL MEAN RIDITS */
%equalmns(&noGroup);

/* (G-1) SCHEFFE`-type GROUP COMPARISONS to the control group (Group1)
*/
     %macro generate(noGroup);
        %do i=2 %to &noGroup;
        %compGroup(&i); %end;
        %mend;
     %generate(&noGroup);

/* Table of Confidence Intervals and Odds */
data interval; set interval;
file print;
put "           &conlim% Simultaneous Confidence Intervals on Mean
Ridits" /;
put "             "
    "Group"  "       Bonferonni"  "         Scheffe`"  "      Odds";
%do i=1 %to &noGroup;
   ii=&i;
   put "             "     ii 3.0  "  "
       BonL&i   10.3  ","  BonU&i   6.3
       SchefL&i 10.3  ","  SchefU&i 6.3
       odds&i 08.2   ":1";
%end;
run;

/* Graphical Analysis */
%grafit(&noGroup,&alpha);

/* Output of Intermediate Calculations */
%if &DIAGNOSE=No %then %goto TRAK;
%else
%if &DIAGNOSE=Yes %then %do;
%myprints;
%end;

%TRAK:
title1; title2; title3;
data timetrak; set timetrak;
time2=time();
Xtime=(time2-time1)/60;
file print;
put _page_ ;
put // "Total Execution Time is " xtime 5.3 " Minutes";
run; quit;

proc datasets nolist;
delete DataFile td td2 meansout ridit riditout fij scheffe
       f equalmns interval grafit timetrak;
run; quit;
%mend ridits;
