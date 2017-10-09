function score = GRassymetry(GRs, GRd)

%%
GRcoord = [ToColumn(GRs-1) ToColumn(GRd) zeros(length(GRd),1)];
d = sum(GRcoord.^2,2);

extremum = GRcoord(argmax(d),:);

score = 0;
for i=2:length(GRd)
    s = cross(GRcoord(i-1,:), GRcoord(i,:));
    score = score + s(3);
end

