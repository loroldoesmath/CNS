function out=hh_functions(t,Y)
    %These are right-hand-side of HH equations for ODE simulation
    
    % LW test 
    disp('TEST')
    %set up the variables
    out=zeros(4,1); %makes output a 2d column vector
    v=Y(1);
    m=Y(2);
    h=Y(3);
    n=Y(4);
    
    %define parameters of the model
    i0=0;
    vna=50;vk=-77;vl=-54.4;
    gna=120;gk=36;gl=0.3;
    c=1;phi=1;
    
    %these are parameters of the applied pulse of current
    ip=10; %strenth
    pon=50;% time on
    poff=300; % time off 
    
    %define functions to be used in the equations
    am=phi*.1*(v+40)./(1-exp(-(v+40)/10));
    bm=phi*4*exp(-(v+65)/18);
    ah=phi*.07*exp(-(v+65)/20);
    bh=phi*1./(1+exp(-(v+35)/10));
    an=phi*.01*(v+55)./(1-exp(-(v+55)/10));
    bn=phi*.125*exp(-(v+65)/80);
    
    %input pulse
    if(t>pon & t<poff)
        is=ip;
    else
        is=0;
    end
    
    
    out(1)=(i0+is-gna.*h.*(v-vna).*m.^3-gk.*(v-vk).*n.^4-gl*(v-vl))/c;
    out(2)=am.*(1-m)-bm.*m;
    out(3)=ah.*(1-h)-bh.*h;
    out(4)=an.*(1-n)-bn.*n;