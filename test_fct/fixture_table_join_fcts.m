t_A = table({'a','b','c','d','e'}', (1:5)', 'variablenames',{'lowletter' 'number'});
t_Arev = flipud(t_A);

t_B = table({'A','B','C','D','E'}', (1:5)', 'variablenames',{'upletter' 'number'});
Brand = [3 2 5 1 4]';
t_Brand = t_B(Brand,:);

t_C = table({'A','B','B','H','I' 'BB'}', [1 2 7:9 2]', 'variablenames',{'upletter' 'number'});

%%

[t_out, Ridx] = leftjoin(t_Arev,t_Brand);

assert(all(Ridx == flipud(sortidx(Brand))))
assert(all(strcmp(t_out.lowletter, t_Arev.lowletter)))
assert(all(strcmpi(t_out.lowletter, t_out.upletter)))

[t_out, Ridx] = leftjoin(t_Arev,t_C);
assert(all(strcmpi(t_C.upletter(Ridx(Ridx>0)), t_out.upletter(Ridx>0))))

%%

[t_out, Lidx, Ridx] = leftinnerjoin(t_Arev,t_Brand)

assert(all(Ridx == flipud(sortidx(Brand))))
assert(all(strcmp(t_out.lowletter, t_Arev.lowletter)))
assert(all(strcmpi(t_out.lowletter, t_out.upletter)))

[t_out, Lidx, Ridx] = leftinnerjoin(t_Arev,t_C)
assert(all(strcmpi(t_C.upletter(Ridx), t_out.upletter)))


[t_out, Lidx, Ridx] = leftinnerjoin(t_C, t_Arev)
