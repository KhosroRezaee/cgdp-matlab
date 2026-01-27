function desc = compute_physchem(seq)
%COMPUTE_PHYSCHEM Compute (Length, NetCharge@pH7, GRAVY, BomanIndex)

tbl = aa_tables();
seq = char(seq);
L = length(seq);

% Charge at pH7 via Henderson-Hasselbalch with standard pKa
pH = 7.0;
pKa = tbl.pKa;

% N-term positive
charge = 1 / (1 + 10^(pH - pKa.Nterm));
% C-term negative
charge = charge - 1 / (1 + 10^(pKa.Cterm - pH));

% Side chains
for i=1:L
    aa = seq(i);
    switch aa
        case 'D'
            charge = charge - 1 / (1 + 10^(pKa.D - pH));
        case 'E'
            charge = charge - 1 / (1 + 10^(pKa.E - pH));
        case 'C'
            charge = charge - 1 / (1 + 10^(pKa.C - pH));
        case 'Y'
            charge = charge - 1 / (1 + 10^(pKa.Y - pH));
        case 'H'
            charge = charge + 1 / (1 + 10^(pH - pKa.H));
        case 'K'
            charge = charge + 1 / (1 + 10^(pH - pKa.K));
        case 'R'
            charge = charge + 1 / (1 + 10^(pH - pKa.R));
    end
end

% GRAVY (mean KD hydropathy)
kd_sum = 0;
bom_sum = 0;
for i=1:L
    aa = seq(i);
    kd_sum = kd_sum + tbl.kd(aa);
    bom_sum = bom_sum + tbl.boman(aa);
end
gravy = kd_sum / L;
boman = bom_sum / L;

desc = struct('length', L, 'charge', charge, 'gravy', gravy, 'boman', boman);
end
