data(nelson)
# Max-likekihood estimate
addt(kv~week|11605/(temp+273.15),xref=150,data=nelson,transf=c(log10,function(x) 10^x))->addtNelson
addtNelson
# Two-pass estimation method used to initialize 
addtTwoPass(kv~week|11605/(temp+273.15),xref=150,data=nelson,transf=c(log10,function(x) 10^x))->addtNelson2
addtNelson2
# degradation vs time
plot(addtNelson,1)
# degradition vs time with free acceleration fit
plot(addtNelson,1,fitFreeAccel=TRUE)
# g-degradation vs time
plot(addtNelson,2)
# g-degradition vs time with free acceleration fit
plot(addtNelson,2,fitFreeAccel=TRUE)