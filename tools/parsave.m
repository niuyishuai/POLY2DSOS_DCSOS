function parsave(filename,p,N,d,density,PSOS,NSOS,time)
    if isa(PSOS,'struct')
        maxdeg=2*PSOS.b.degree;
        basislen=length(PSOS.b);
        save(filename,'p','N','d','density','maxdeg','basislen','time');
    else
        maxdeg=PSOS.degree;
        save(filename,'p','N','d','density','maxdeg','time');
    end
end