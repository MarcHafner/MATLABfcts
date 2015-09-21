function mIdx = memberidx(A, B)
% mIdx = memberidx(A, B)
%
% If A and B are numeric arrays, logical arrays, character arrays,
% categorical arrays, datetime arrays, duration arrays, or cell arrays of
% strings, then Locb contains the lowest index in B for each value in A
% that is a member of B. The output array, Locb, contains 0 wherever A is
% not a member of B.
%
% If A and B are tables, then Locb contains the lowest index in B for each
% row in A that is also a row in B. The output vector, Locb, contains 0
% whenever A is not a row of B.

[~,mIdx] = ismember(A, B);
