function [bin_centers,pdf_ngrc,pdf_exact] = plot_MG_pdf(Xtot)
% PDF comparison for exact and NGRC systems
% Exact system: (Tpr1, U1)
% NGRC system:  (Tpr1, XP1)
% Xtot=[Tpr1 U1 XP1]; % source array
%Tpr1=Xtot(:,1);
U1=Xtot(:,2);
XP1=Xtot(:,3);

u_exact = U1(:);
u_ngrc  = XP1(:);

% Remove NaN/Inf values if present
u_exact = u_exact(isfinite(u_exact));
u_ngrc  = u_ngrc(isfinite(u_ngrc));

% Use common bin limits so both PDFs are directly comparable
all_vals = [u_exact; u_ngrc];
nBins = 80;
bin_edges = linspace(min(all_vals), max(all_vals), nBins + 1);

[pdf_exact, edges] = histcounts(u_exact, bin_edges, 'Normalization', 'pdf');
[pdf_ngrc,  ~]     = histcounts(u_ngrc,  bin_edges, 'Normalization', 'pdf');

bin_centers = 0.5 * (edges(1:end-1) + edges(2:end));



end