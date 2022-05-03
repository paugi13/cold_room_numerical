function [coord_total] = node_coord(nod_reforc_1,nod_poli, nod_reforc_2, nod_poma, ax_reforc, ax_poli, ax_poma, l_air, total_nod)
% Sets coordinates for every registered node in the mesh.

coord_reforc_1 = zeros(1,size(nod_reforc_1,2));
coord_poli = zeros(1,size(nod_poli,2));
coord_reforc_2 = zeros(1,size(nod_reforc_2,2));
coord_poma = zeros(1,size(nod_poma,2));
coord_total = zeros(1,total_nod);

for i=2:1:size(coord_reforc_1,2)
    coord_reforc_1(1,i) = coord_reforc_1(1,i-1) + ax_reforc;
end

coord_poli(1,1) = coord_reforc_1(end) + ax_poli;

for i=2:1:size(coord_poli,2)
    coord_poli(1,i) = coord_poli(1,i-1) + ax_poli;
end

coord_reforc_2(1,1) = coord_poli(end) + ax_reforc;

for i=2:1:size(coord_reforc_2,2)
    coord_reforc_2(1,i) = coord_reforc_2(1,i-1) + ax_reforc;
end

coord_poma(1,1) = coord_reforc_2(end) + l_air;

for i=2:1:size(coord_poma,2)
    coord_poma(1,i) = coord_poma(1,i-1) + ax_poma;
end

coord_total(1, nod_reforc_1) = coord_reforc_1;
coord_total(1, nod_poli) = coord_poli;
coord_total(1, nod_reforc_2) = coord_reforc_2;
coord_total(1, nod_poma) = coord_poma;



