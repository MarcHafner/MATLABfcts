function set_RandSeed(Seed)
% function set_RandSeed(Seed)
%   set the random stream to mt 19937ar with seed  Seed

s = RandStream('mt19937ar','Seed',Seed);
RandStream.setGlobalStream(s);
