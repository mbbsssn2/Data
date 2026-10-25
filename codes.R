############################################
#PREAMPLE
############################################

install.packages("moments")
install.packages("evd")
install.packages("KSgeneral")
install.packages("DescTools")
install.packages("extRemes")
library(moments)
library(evd)
library(KSgeneral)
library(DescTools)
library(extRemes)

############################################
#READ THE DATA
############################################

x2=matrix(0,132,142)
z=read.csv("c:\\Work\\bfaso\\rainfall2.csv",header=T)
year=z[,1]
month=z[,2]
for (k in seq(1,142))
{ii=0
x=z[,(k+3)]
for (i in seq(2011,2021)) for (j in seq(1,12))
{ii=ii+1
x2[ii,k]=max(x[year==i&month==j&!is.na(x)])}}

x2[abs(x2)==Inf]=NA
label=names(z[seq(4,145)])

############################################
#COMPUTE PROPORTIONS OF ZEROS
############################################

pp=rep(0,142)
for (i in seq(1,142))
{xx=x2[,i]
xx2=xx[!is.na(xx)]
pp[i]=length(xx2[xx2==0])/length(xx2)}

############################################
#PERFORM RUNS TEST AND
#COMPUTE EXTREMAL INDEX
############################################

pv=rep(0,142)
ei=rep(0,142)
pv2=rep(0,132)
for (k in seq(1,142))
{xx=x2[,k]
th=quantile(xx[xx>0&!is.na(xx)],probs=c(0.95))
for (j in seq(1,132))
{xx2=xx[-j]
yy=xx2[xx2>0&!is.na(xx2)]
pv2[j]=RunsTest(yy)$p.value}
ei[k]=extremalindex(xx[xx>0&!is.na(xx)],threshold=th)[1]
pv[k]=median(pv2)}

############################################
#COMPUTE BASIC STATISTICS
############################################

for (i in seq(1,142))
{xx=x2[,i]
cat(label[i]," ",min(xx[!is.na(xx)])," ",max(xx[!is.na(xx)])," ",
mean(xx[!is.na(xx)])," ",median(xx[!is.na(xx)])," ",
sd(xx[!is.na(xx)])," ",skewness(xx[!is.na(xx)])," ",
kurtosis(xx[!is.na(xx)]),"\n")}

############################################
#DETERMINE THE BEST MODEL
############################################

q1=qchisq(0.95,df=1)
q2=qchisq(0.95,df=2)

ibest=rep(0,142)
for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]),
scale=exp(p[2]+p[4]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est3=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est4=optim(f,par=c(est2$par[1],est2$par[2],est2$par[3],est2$par[4],0),
control=list(maxit=1000000))

if (2*(est1$value-est2$value)<q1&2*(est1$value-est3$value)<q1
&2*(est1$value-est4$value)<q2) ibest[ii]=1
if (2*(est1$value-est2$value)>q1&est2$value<est3$value
&2*(est2$value-est4$value)<q1) ibest[ii]=2
if (2*(est1$value-est3$value)>q1&est3$value<est2$value
&2*(est3$value-est4$value)<q1) ibest[ii]=3
if (2*(est1$value-est4$value)>q2&2*(est2$value-est4$value)>q1
&2*(est3$value-est4$value)>q1) ibest[ii]=4}


############################################
#PARAMETER ESTIMATES, SES AND KS P-VALUE WHEN
#THE BEST MODEL IS CONSTANT
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==1)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
ee=est1$value
p1=est1$par[1]
p2=est1$par[2]
p3=est1$par[3]

pv=rep(0,n)
pp1=pv
pp2=pv
pp3=pv
for (j in seq(1,n))
{s=fgev(x[-j],std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x[-j],loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
pp1[j]=est1$par[1]
pp2[j]=est1$par[2]
pp3[j]=est1$par[3]
pv[j]=ks.test(x[-j],"pgev",exp(pp1[j]),exp(pp2[j]),pp3[j])$p.value}

cat(label[ii]," & ",formatC(ee,format="f",digits=3),
" & ",formatC(p1,format="f",digits=3)," (",
formatC((n-1)*sd(pp1)/sqrt(n),format="f",digits=3),"), ",
formatC(p2,format="f",digits=3)," (",
formatC((n-1)*sd(pp2)/sqrt(n),format="f",digits=3),"), ",
formatC(p3,format="f",digits=3)," (",
formatC((n-1)*sd(pp3)/sqrt(n),format="f",digits=3),") ",
" & ",formatC(median(pv),format="f",digits=3),"\\\\","\n")}}

############################################
#PARAMETER ESTIMATES, SES AND KS P-VALUE WHEN
#THE BEST MODEL IS LINEAR LOCATION
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==2)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
ee=est2$value
p1=est2$par[1]
p2=est2$par[2]
p3=est2$par[3]
p4=est2$par[4]

pv=rep(0,n)
pp1=pv
pp2=pv
pp3=pv
pp4=pv
pobs=pv
for (j in seq(1,n))
{j1=j-1
j2=j+1
f=function (p)
{tt=0
if (j1>=1) {for (i in seq(1,j1))
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),scale=exp(p[2]),shape=p[3],log=TRUE)}
if (j2<=n) {for (i in seq(j2,n))
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),scale=exp(p[2]),shape=p[3],log=TRUE)}
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(p1,p2,p3,p4),control=list(maxit=1000000))
pp1[j]=est2$par[1]
pp2[j]=est2$par[2]
pp3[j]=est2$par[3]
pp4[j]=est2$par[4]

if (j==1) pexp2=seq(2,n)/(n+1)
if (j==n) pexp2=seq(1,(n-1))/(n+1)
if (j1>=1&j2<=n) pexp2=c(seq(1,j1),seq(j2,n))/(n+1)
for (i in seq(1,n))
pobs[i]=pgev(x[i],loc=exp(pp1[j]+pp4[j]*time2[i]),scale=exp(pp2[j]),shape=pp3[j])
if (j==1) pobs2=pobs[seq(2,n)]
if (j==n) pobs2=pobs[seq(1,(n-1))]
if (j1>=1&j2<=n) pobs2=c(pobs[seq(1,j1)],pobs[seq(j2,n)])
ss=sort(pobs2)
pv[j]=cont_ks_c_cdf(max(abs(ss-pexp2)),n-1)}

cat(label[ii]," & ",formatC(ee,format="f",digits=3),
" & ",formatC(p1,format="f",digits=3)," (",
formatC((n-1)*sd(pp1)/sqrt(n),format="f",digits=3),"), ",
formatC(p2,format="f",digits=3)," (",
formatC((n-1)*sd(pp2)/sqrt(n),format="f",digits=3),"), ",
formatC(p3,format="f",digits=3)," (",
formatC((n-1)*sd(pp3)/sqrt(n),format="f",digits=3),"), ",
formatC(p4,format="f",digits=3)," (",
formatC((n-1)*sd(pp4)/sqrt(n),format="f",digits=3),") ",
" & ",formatC(median(pv),format="f",digits=3),"\\\\","\n")}}

############################################
#PARAMETER ESTIMATES, SES AND KS P-VALUE WHEN
#THE BEST MODEL IS LINEAR SCALE
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==3)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=0
for (i in 1:n)
tt=tt-dgev(x[i],loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n)
tt=tt-dgev(x[i],loc=exp(p[1]),scale=exp(p[2]+p[4]*time2[i]),
shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
ee=est2$value
p1=est2$par[1]
p2=est2$par[2]
p3=est2$par[3]
p4=est2$par[4]

pv=rep(0,n)
pp1=pv
pp2=pv
pp3=pv
pp4=pv
pobs=pv
for (j in seq(1,n))
{j1=j-1
j2=j+1
f=function (p)
{tt=0
if (j1>=1) {for (i in seq(1,j1))
tt=tt-dgev(x[i],loc=exp(p[1]),scale=exp(p[2]+p[4]*time2[i]),shape=p[3],log=TRUE)}
if (j2<=n) {for (i in seq(j2,n))
tt=tt-dgev(x[i],loc=exp(p[1]),scale=exp(p[2]+p[4]*time2[i]),shape=p[3],log=TRUE)}
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(p1,p2,p3,p4),control=list(maxit=1000000))
pp1[j]=est2$par[1]
pp2[j]=est2$par[2]
pp3[j]=est2$par[3]
pp4[j]=est2$par[4]

if (j==1) pexp2=seq(2,n)/(n+1)
if (j==n) pexp2=seq(1,(n-1))/(n+1)
if (j1>=1&j2<=n) pexp2=c(seq(1,j1),seq(j2,n))/(n+1)
for (i in seq(1,n))
pobs[i]=pgev(x[i],loc=exp(pp1[j]),scale=exp(pp2[j]+pp4[j]*time2[i]),shape=pp3[j])
if (j==1) pobs2=pobs[seq(2,n)]
if (j==n) pobs2=pobs[seq(1,(n-1))]
if (j1>=1&j2<=n) pobs2=c(pobs[seq(1,j1)],pobs[seq(j2,n)])
ss=sort(pobs2)
pv[j]=cont_ks_c_cdf(max(abs(ss-pexp2)),n-1)}

cat(label[ii]," & ",formatC(ee,format="f",digits=3),
" & ",formatC(p1,format="f",digits=3)," (",
formatC((n-1)*sd(pp1)/sqrt(n),format="f",digits=3),"), ",
formatC(p2,format="f",digits=3)," (",
formatC((n-1)*sd(pp2)/sqrt(n),format="f",digits=3),"), ",
formatC(p3,format="f",digits=3)," (",
formatC((n-1)*sd(pp3)/sqrt(n),format="f",digits=3),"), ",
formatC(p4,format="f",digits=3)," (",
formatC((n-1)*sd(pp4)/sqrt(n),format="f",digits=3),") ",
" & ",formatC(median(pv),format="f",digits=3),"\\\\","\n")}}

############################################
#PARAMETER ESTIMATES, SES AND KS P-VALUE WHEN
#THE BEST MODEL IS LINEAR LOCATION AND SCALE
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==4)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=0
for (i in 1:n)
tt=tt-dgev(x[i],loc=exp(p[1]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n)
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))

f=function (p)
{tt=0
for (i in 1:n)
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est3=optim(f,par=c(est2$par[1],est2$par[2],est2$par[3],est2$par[4],0),
control=list(maxit=1000000))
ee=est3$value
p1=est3$par[1]
p2=est3$par[2]
p3=est3$par[3]
p4=est3$par[4]
p5=est3$par[5]

pv=rep(0,n)
pp1=pv
pp2=pv
pp3=pv
pp4=pv
pp5=pv
pobs=pv
for (j in seq(1,n))
{j1=j-1
j2=j+1
f=function (p)
{tt=0
if (j1>=1) {for (i in seq(1,j1))
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)}
if (j2<=n) {for (i in seq(j2,n))
tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)}
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(p1,p2,p3,p4,p5),control=list(maxit=1000000))
pp1[j]=est2$par[1]
pp2[j]=est2$par[2]
pp3[j]=est2$par[3]
pp4[j]=est2$par[4]
pp5[j]=est2$par[5]

if (j==1) pexp2=seq(2,n)/(n+1)
if (j==n) pexp2=seq(1,(n-1))/(n+1)
if (j1>=1&j2<=n) pexp2=c(seq(1,j1),seq(j2,n))/(n+1)
for (i in seq(1,n))
pobs[i]=pgev(x[i],loc=exp(pp1[j]+pp4[j]*time2[i]),
scale=exp(pp2[j]+pp5[j]*time2[i]),shape=pp3[j])
if (j==1) pobs2=pobs[seq(2,n)]
if (j==n) pobs2=pobs[seq(1,(n-1))]
if (j1>=1&j2<=n) pobs2=c(pobs[seq(1,j1)],pobs[seq(j2,n)])
ss=sort(pobs2)
pv[j]=cont_ks_c_cdf(max(abs(ss-pexp2)),n-1)}

cat(label[ii]," & ",formatC(ee,format="f",digits=3),
" & ",formatC(p1,format="f",digits=3)," (",
formatC((n-1)*sd(pp1)/sqrt(n),format="f",digits=3),"), ",
formatC(p2,format="f",digits=3)," (",
formatC((n-1)*sd(pp2)/sqrt(n),format="f",digits=3),"), ",
formatC(p3,format="f",digits=3)," (",
formatC((n-1)*sd(pp3)/sqrt(n),format="f",digits=3),"), ",
formatC(p4,format="f",digits=3)," (",
formatC((n-1)*sd(pp4)/sqrt(n),format="f",digits=3),") ",
formatC(p5,format="f",digits=3)," (",
formatC((n-1)*sd(pp5)/sqrt(n),format="f",digits=3),") ",
" & ",formatC(median(pv),format="f",digits=3),"\\\\","\n")}}

############################################
#COMPUTE THE UPPER BOUNDS AND SES
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==1)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
if (est1$par[3]<0) up=exp(est1$par[1])-exp(est1$par[2])/est1$par[3]

up2=rep(0,n)
for (j in seq(1,n))
{s=fgev(x[-j],std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x[-j],loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
ee=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
up2[j]=NA
if (ee$par[3]<0) up2[j]=exp(ee$par[1])-exp(ee$par[2])/ee$par[3]}

if (est1$par[3]<0) cat(label[ii]," & ",formatC(up,format="f",digits=3),
" (",formatC((n-1)*sd(up2[!is.na(up2)])/sqrt(n),format="f",digits=3),")",
"\\\\","\n")}}


############################################
#COMPUTE THE RETURN LEVELS AND SES
############################################

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
n=length(x)

s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])

rr2=rep(0,n)
rr3=rr2
rr5=rr2
rr10=rr2
for (j in seq(1,n))
{s=fgev(x[-j],std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x[-j],loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
rr2[j]=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
rr3[j]=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
rr5[j]=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
rr10[j]=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])}

cat(label[ii]," & ",formatC(r2,format="f",digits=3)," (",
formatC((n-1)*sd(rr2)/sqrt(n),format="f",digits=3),"), & ",
formatC(r3,format="f",digits=3)," (",
formatC((n-1)*sd(rr3)/sqrt(n),format="f",digits=3),"), & ",
formatC(r5,format="f",digits=3)," (",
formatC((n-1)*sd(rr5)/sqrt(n),format="f",digits=3),") & ",
formatC(r10,format="f",digits=3)," (",
formatC((n-1)*sd(rr10)/sqrt(n),format="f",digits=3),") ","\\\\","\n")}

############################################
#COMPUTE RETURN LEVELS FOR STATIONARY &
#NON-STATIONARY MODELS ONE YEAR AHEAD
############################################

q1=qchisq(0.95,df=1)
q2=qchisq(0.95,df=2)

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==1)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==2)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+12)),
scale=exp(est2$par[2]),shape=est2$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+12)),
scale=exp(est2$par[2]),shape=est2$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+12)),
scale=exp(est2$par[2]),shape=est2$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+12)),
scale=exp(est2$par[2]),shape=est2$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==3)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]),
scale=exp(p[2]+p[4]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est3=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+12)),shape=est3$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+12)),shape=est3$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+12)),shape=est3$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+12)),shape=est3$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==4)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est4=optim(f,par=c(est2$par[1],est2$par[2],est2$par[3],est2$par[4],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+12)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+12)),shape=est4$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+12)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+12)),shape=est4$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+12)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+12)),shape=est4$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+12)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+12)),shape=est4$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}}

############################################
#COMPUTE RETURN LEVELS FOR STATIONARY &
#NON-STATIONARY MODELS TWO YEARS AHEAD
############################################

q1=qchisq(0.95,df=1)
q2=qchisq(0.95,df=2)

for (ii in seq(1,142))
{xx=x2[,ii]
time=seq(1,132)
x=xx[!is.na(xx)&xx>0]
time2=time[!is.na(xx)&xx>0]
n=length(x)

if (ibest[ii]==1)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est1$par[1]),
scale=exp(est1$par[2]),shape=est1$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==2)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+24)),
scale=exp(est2$par[2]),shape=est2$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+24)),
scale=exp(est2$par[2]),shape=est2$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+24)),
scale=exp(est2$par[2]),shape=est2$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est2$par[1]+est2$par[4]*(time2[n]+24)),
scale=exp(est2$par[2]),shape=est2$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==3)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]),
scale=exp(p[2]+p[4]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est3=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+24)),shape=est3$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+24)),shape=est3$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+24)),shape=est3$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est3$par[1]),
scale=exp(est3$par[2]+est3$par[4]*(time2[n]+24)),shape=est3$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}

if (ibest[ii]==4)
{s=fgev(x,std.err=FALSE)
f=function (p)
{tt=-sum(dgev(x,loc=exp(p[1]),scale=exp(p[2]),shape=p[3],log=TRUE))
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est1=optim(f,par=c(log(s$estimate[1]),log(s$estimate[2]),s$estimate[3]),control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est2=optim(f,par=c(est1$par[1],est1$par[2],est1$par[3],0),
control=list(maxit=1000000))
#
f=function (p)
{tt=0
for (i in 1:n) tt=tt-dgev(x[i],loc=exp(p[1]+p[4]*time2[i]),
scale=exp(p[2]+p[5]*time2[i]),shape=p[3],log=TRUE)
if (is.na(tt)) tt=1.0e20
if (abs(tt)>1.0e20) tt=1.0e20
return(tt)}
est4=optim(f,par=c(est2$par[1],est2$par[2],est2$par[3],est2$par[4],0),
control=list(maxit=1000000))
r2=qgev(1-1/(12*2*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+24)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+24)),shape=est4$par[3])
r3=qgev(1-1/(12*3*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+24)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+24)),shape=est4$par[3])
r5=qgev(1-1/(12*5*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+24)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+24)),shape=est4$par[3])
r10=qgev(1-1/(12*10*(1-pp[ii])),loc=exp(est4$par[1]+est4$par[4]*(time2[n]+24)),
scale=exp(est4$par[2]+est4$par[5]*(time2[n]+24)),shape=est4$par[3])
cat(label[ii]," & ",formatC(r2,format="f",digits=3),
" & ",formatC(r3,format="f",digits=3),
" & ",formatC(r5,format="f",digits=3)," & ",
formatC(r10,format="f",digits=3),"\n")}}