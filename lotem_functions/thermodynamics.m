function [outputArg1,outputArg2] = thermodynamics(inputArg1,inputArg2)

n = 216;
FILE *ptr;
char name[FILENAME_MAX];
    
for 1:1:n
    %fid1 = fopen('C:\Users\ldrul\OneDrive\Documents\GitHub\Genetic-Supression-predictor\lotem_functions\txt_files\file1.txt','W+');
    snprintf(name, sizeof(name), "%d.txt", i);
    fopen_s(&ptr, name, "w");
    %operations to fill data into file i.txt;
    fclose(ptr);
end
%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
end

