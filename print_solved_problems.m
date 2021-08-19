
PS = load('indices_solved_problems_49_starting point.mat');

vector_with_solved_problems = PS.vector_with_solved_problems;

cd ..
cd Dataset;
a = dir;
[p , l] = size(a);


for i  = 3: p
    NS = load(a(i).name);
    if vector_with_solved_problems(i-2)==1
        a(i).name
    end
end