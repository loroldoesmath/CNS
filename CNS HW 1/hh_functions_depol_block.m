function out = hh_functions_depol_block(t, Y)
    % Modified HH equations for depolarization block (very strong current)
    
    out = zeros(4,1);
    v = Y(1);
    m = Y(2);
    h = Y(3);
    n = Y(4);
    
    % Parameters
    i0 = 0;
    vna = 50;
    vk = -77;
    vl = -54.4;
    gna = 120;
    gk = 36;
    gl = 0.3;
    c = 1;
    phi = 1;
    
    % Very strong constant current to cause depolarization block
    is = 100;  % Much stronger than needed for normal spiking
    
    % Rate functions
    am = phi * 0.1 * (v + 40) ./ (1 - exp(-(v + 40) / 10));
    bm = phi * 4 * exp(-(v + 65) / 18);
    ah = phi * 0.07 * exp(-(v + 65) / 20);
    bh = phi * 1 ./ (1 + exp(-(v + 35) / 10));
    an = phi * 0.01 * (v + 55) ./ (1 - exp(-(v + 55) / 10));
    bn = phi * 0.125 * exp(-(v + 65) / 80);
    
    % ODEs
    out(1) = (i0 + is - gna .* h .* (v - vna) .* m.^3 - ...
              gk .* (v - vk) .* n.^4 - gl * (v - vl)) / c;
    out(2) = am .* (1 - m) - bm .* m;
    out(3) = ah .* (1 - h) - bh .* h;
    out(4) = an .* (1 - n) - bn .* n;
    end