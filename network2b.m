% leveling network
% data: coefficient, heights, observations for the leveling network

% coefficient matrix
Coeff=[1 0 0 0 -1 0 0
    1 0 0 0 0 -1 0
    1 -1 0 0 0 0 0
    0 -1 0 0 0 1 0
    0 -1 0 0 0 0 1
    0 1 -1 0 0 0 0
    0 0 -1 0 0 0 1
    0 0 -1 1 0 0 0
    0 0 0 1 0 0 -1
    0 0 0 1 0 -1 0
    0 0 0 1 -1 0 0];

% heights
H=[104
    102.4
    100.72
    105.21
    101.15
    103.3
    104.13];

[nn uu]=size(Coeff);

Olculer=Coeff*H; % observations

% Leveling lines in km
ghetero=[0.920
    1.254
    1.602
    0.987
    1.210
    1.320
    0.895
    0.959
    1.302
    0.842
    0.801];

