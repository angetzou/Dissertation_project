# Dissertation_project
Dissertation_project

Author: Angelis Tzouchas

Format of the file

create a folder such as : dissertation folder and paste Dataset and matrix-free-ipm files inside

dissertation_folder
	Dataset
	matrix-free-ipm


In matrix-free-ipm we have all the codes.
In Dataset we have the Netlib collection.


In order to use the interior-point method for the regularized primal-dual problem.

We should run the loaddata.m file

In the loaddata.m we can set the follwing parameters


maxit_ipm =250;

maxit_kryl =400;

gamma = 1e-4; 

delta = 1e-3;

Rd_max = 1e-4; % is usefull for ipm with dynamic regularization

Rd_min = 1e-6; % is usefull for ipm with dynamic regularization

k = 50; % indicates the number of pivots for the Partial Cholesky Decomposition

on the line 63-65 we can select which interior point method we want to use

line: 63. The first one is with uniform regularization
line: 64. The second one is with dynamic regularization
These models uses the starting point x=1,y=0,s=1

line: 65. The third interior point method, uses the starting_point functions that we choose


where the third starting point is the function starting_point while the second is the starting_point2

Depending which driver we choose on the specific driver file we can set the tolerances
in lines 56-59.


So suppose that we need to use the third starting point we use the command
in line 69, or if we want to use the starting point 2 we use the line 68

Once the model has finished, with the command: sum(vector_with_solved_problems)
we find the number of problems that are solved, out of 82 problems with x=>0.

In order to see which problems are solved, we save the 'vector_with_solved_problems'
with the command save('name.mat','vector_with_solved_problems')

Up to this poin the save command save the name.mat file on the Dataset file. 
We navigate into Dataset file cut the name.mat file from the dataset and paste it on
matrix-free-ipm file. Then by running the print_solved_problems.m 
We can see the names of the solved problems from Netlib collection.

By running the load_qap_problems we use IPM  for solving the QAP problems

The LP_Convert_to_Standard_Form.m was extracted by the following Github link: https://github.com/spougkakiotis/IP_PMM/blob/master/Matlab_code/LP_Convert_to_Standard_Form.m
The author fo the above code is Spyridon Pougkakiotis.
