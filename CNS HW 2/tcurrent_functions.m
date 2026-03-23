function out=tcurrent_functions(t,Y)
    %These are right-hand-side of equations for ODE simulation with Ca
    %T-current
    global i0 i1 ton tdur
    
    %set up the variables
    out=zeros(4,1); %makes output a 4d column vector
    v=Y(1);
    h=Y(2);
    hna=Y(3);
    n=Y(4);
    
    
    %parameters
    
    if(t<ton|t>ton+tdur)
        i=i0;
    else
        i=i1;
    end
    
    faraday=96520;
    rgas=8.3134;
    celsius=25;temp=273.15+celsius;
    
    gl=.05;el=-60;
    cao=2;pcat=.15;cai=1e-4;
    %for blocking spikes, othersise uncomment next line
    gna=0;gk=0;
    %gna=8;gk=4;
    ena=55;ek=-90;
    
    %functions
    xi=v*faraday*2/(rgas*1000*temp);
    cfedrive=.002*faraday*xi.*(cai-cao*exp(-xi))./(1-exp(-xi));
    minf=1./(1+exp(-(v+59)/6.2));
    m=minf;
    i_cat=pcat*m.^2.*h.*cfedrive;
    
    amna=.091*(v+38)./(1-exp(-(v+38)/5));
    bmna=-.062*(v+38)./(1-exp((v+38)/5));
    mna=amna./(amna+bmna);
    
    hinf=1./(1+exp((v+83)/4));
    tauh=22.7+.27./(exp((v+48)/4)+exp(-(v+407)/50));
    
    ahna=.016*exp((-55-v)/15);
    bhna=2.07./(1+exp((17-v)/21));
    
    ank=.01*(-45-v)./(exp((-45-v)/5)-1);
    bnk=.17*exp((-50-v)/40);
    
    %equations
    out(1)=-gl*(v-el)-i_cat+i-gna*mna^3.*hna.*(v-ena)-gk.*n^4.*(v-ek);
    out(2)=(hinf-h)./tauh;
    out(3)=ahna.*(1-hna)-bhna.*hna;
    out(4)=ank.*(1-n)-bnk.*n;