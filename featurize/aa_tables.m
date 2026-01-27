function tbl = aa_tables()
%AA_TABLES returns per-residue contributions used throughout the pipeline.
%
% - KD: Kyte-Doolittle hydropathy for GRAVY
% - Boman: 'solubility' table used by modlAMP (Boman-like index)
% - Charge: side chain charges at pH7 approximated via pKa
%
% NOTE (Windows/MATLAB quirk): Use CHAR vectors for containers.Map keys so
% cellstr() yields 20 single-letter keys (not 1 string key).

aa = 'ACDEFGHIKLMNPQRSTVWY';   % char row vector (20)
tbl.aa = aa;

keys = cellstr(aa(:));         % 20x1 cell array of 1-char strings

% Kyte-Doolittle hydropathy (standard values)
kd_vals = [ ...
    1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, ...
    1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3 ];
tbl.kd = containers.Map(keys, num2cell(kd_vals));

% Boman "solubility" values (used by modlAMP BomanIndex)
boman_vals = [ ...
    0.17, -0.02, 1.23, 1.23, -0.38, 0.01, 0.96, -0.31, 0.99, -0.56, ...
    -0.37, 0.95, 0.60, 0.89, 1.01, 0.13, 0.14, -0.22, 0.50, 0.96 ];
tbl.boman = containers.Map(keys, num2cell(boman_vals));

% Side-chain pKa (approx) for charge computation at pH7
tbl.pKa = struct();
tbl.pKa.Nterm = 9.69;
tbl.pKa.Cterm = 2.34;
tbl.pKa.D = 3.90;
tbl.pKa.E = 4.07;
tbl.pKa.C = 8.37;
tbl.pKa.Y = 10.07;
tbl.pKa.H = 6.04;
tbl.pKa.K = 10.54;
tbl.pKa.R = 12.48;

end
