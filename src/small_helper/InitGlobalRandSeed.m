function InitGlobalRandSeed(s)
% InitGlobalRandSeed(s)
%   initialize the global stream for random numbers (using algorithm
%   mt19937ar with seed  s

    s = RandStream('mt19937ar','Seed',s);
    RandStream.setGlobalStream(s);
