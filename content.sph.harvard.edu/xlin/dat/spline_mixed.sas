*Fit a smoothing spline using PROC MIXED; 

options ls=72;
data diabetes;
  infile 'diabetes.dat';
  input age basedef cpep;
lncpep=log(cpep);


* Identify distinct values of age;
data knots;
 set diabetes;

proc sort data=knots;
by age;

data knots;
 set knots (keep=age);
 by age;
 if first.age;


proc iml;

use diabetes;
read all into w;

use knots;
read all into t0;

t       = w[,1];
y       = w[,4];
n       = nrow(w);
nk      = nrow(t0);
h       = j(nk-1,1,0);
q       = j(nk,nk-2,0);
rr      = j(nk-2,nk-2,0);
g       = j(nk-2,nk-2,0);
l       = j(nk,nk-2,0);
b       = j(nk,nk-2,0);
z       = j(n,nk-2,0);

*calculate h;
do i=1 to nk-1;
 h[i]=t0[i+1]-t0[i];
end;

*calculate Q;
do j=1 to nk-2;
 q[j,j]=1/h[j];
 q[j+1,j]=-1/h[j+1]-1/h[j];
 q[j+2,j]=1/h[j+1];
end;

*calculate R;
do j=1 to nk-2;
 rr[j,j]=(h[j]+h[j+1])/3;
end;
rr[1,2]=h[2]/6;
do j=2 to nk-3;
 rr[j,j-1]=h[j]/6;
 rr[j,j+1]=h[j+1]/6;
end;
rr[nk-2,nk-3]=h[nk-2]/6;

*calculate G;
g=t(half(rr));

*calculate L;
l=q*t(inv(g));

*calculate B;
b=l*inv(t(l)*l);

*find index of element e in vector t0, i.e, calculate the incidence matrix N;
start index(e,t0);
 do i=1 to nrow(t0);
  if e=t0[i] then j=i;
 end;
 return(j);
finish;

*define design matrix Z=NB for introduced random effects;
do i=1 to n;
 do j=1 to nk-2;
  z[i,j]=b[index(t[i],t0),j];
 end;
end;

*create SAS data set for PROC MIXED;
datmtr=y||t||z;
create diabetes2 from datmtr;
append from datmtr;

quit;

data diabetes2;
 set diabetes2;
 rename col1=lncpep col2=age col3-col37=z1-z35;

proc mixed data=diabetes2 covtest;
 model lncpep=age/s outpred=fit;
 random z1-z35/type=toep(1) s;


proc sort data=fit;
by age;

data fit;
 set fit (keep=age Pred StdErrPred);
 by age;
 if first.age;
file "spline_mixed.out";
put age Pred StdErrPred;

proc print data=fit;
var age Pred StdErrPred;

*Fit a smoothing spline via mixed models using macro spmm.mac;
%include '$HOME/gamm/Daowen_program/spmm.mac';

%spmm(data=diabetes,dep=lncpep,smthvar=age,outband=ci, print=y);


data one;
set ci;
file 'spmm_fit.out';
put age lncpep F_Se F_L95 F_U95  B_Se  B_L95 B_U95;




