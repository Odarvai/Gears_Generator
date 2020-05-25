%% main code for minimizing the fitness function using GA (Genetic Algorithm)
ObjFcn = @gearFitness;
nvars = 5;
LB = [1 1 -0.5 0.2 7.25];
UB = [32 16 1 0.8 15];
ConsFcn = @Constraints;

[x, fval] = ga(ObjFcn, nvars, [], [], [], [], LB, UB, ConsFcn);