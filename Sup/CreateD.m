%% CIV 2118 - Método dos Elementos Finitos - 2022.2
% Trabalho Final - Parte 1
% Aluno: João Guilherme M. Alvarez & Camila Alves
% Matricula: 2220784 & 
%
% Função para o cálculo da matrix D.
%
%%
function [D] = CreateD(E,Ni)
    m = [1 Ni 0; Ni 1 0; 0 0 ((1 - Ni) / 2)];
    c = E / (1 - (Ni^2));

    D = c .* m;
end