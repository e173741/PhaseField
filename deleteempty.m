file1=['coord_93.mat'] ;
A2x=load(file1) ;
coord1=A2x.coord ;
coord2=coord1(~any(cellfun('isempty',coord1), 2), :);
save('coord_93.mat','coord2','-v7.3') ;
