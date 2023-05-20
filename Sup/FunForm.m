%% CIV 2118 - Método dos Elementos Finitos - 2022.2
% Trabalho Final - Parte 1
% Aluno: João Guilherme M. Alvarez & Camila Alves
% Matricula: 2220784 & 
%
% Função para o cálculo das funções de formas e suas derivadas.
%
%%
function [FF, dFF] = FunForm(gauss, i, j)
    xi = gauss(i, 1);
    eta = gauss(j, 1);

    %Funções de forma
    N1 = 0.25 * (1 - xi) * (eta - 1) * (eta + xi + 1);
    N2 = 0.25 * (1 + xi) * (eta - 1) * (eta - xi + 1);
    N3 = 0.25 * (1 + xi) * (eta + 1) * (eta + xi - 1);
    N4 = 0.25 * (xi - 1) * (eta + 1) * (xi - eta + 1);
    N5 = 0.5 * (1 - eta) * (1 - xi^2);
    N6 = 0.5 * (1 + xi) * (1 - eta^2);
    N7 = 0.5 * (1 + eta) * (1 - xi^2);
    N8 = 0.5 * (1 - xi) * (1 - eta^2);

    %Derivadas das Funções de Forma
    dN1e = 0.25 * (-((xi - 1) * ((2 * eta) + xi)));
    dN1x = 0.25 * (-((eta - 1) * ((2 * xi) + eta)));
    dN2e = 0.25 * ((xi + 1) * ((2 * eta) - xi));
    dN2x = 0.25 * (-((eta - 1) * ((2 * xi) - eta)));
    dN3e = 0.25 * ((xi + 1) * ((2 * eta) + xi));
    dN3x = 0.25 * ((eta + 1) * ((2 * xi) + eta));
    dN4e = 0.25 * (-((xi - 1) * ((2 * eta) - xi)));
    dN4x = 0.25 * ((eta + 1) * ((2 * xi) - eta));
    dN5e = 0.5 * (xi^2 - 1);
    dN5x = (eta * xi) - xi;
    dN6e = (-(eta * xi)) - eta;
    dN6x = 0.5 * (1 - (eta^2));
    dN7e = 0.5 * (1 - (xi^2));
    dN7x = (-(eta * xi)) - xi;
    dN8e = (eta * xi) - eta;
    dN8x = 0.5 * ((eta^2) - 1);
    
    FF = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0;...
          0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];

%     dFF = [diff(N1,eta) diff(N2,eta) diff(N3,eta) diff(N4,eta) diff(N5,eta) diff(N6,eta) diff(N7,eta) diff(N8,eta);...
%            diff(N1,xi) diff(N2,xi) diff(N3,xi) diff(N4,xi) diff(N5,xi) diff(N6,xi) diff(N7,xi) diff(N8,xi)];
% 
    dFF = [dN1x dN2x dN3x dN4x dN5x dN6x dN7x dN8x;...
           dN1e dN2e dN3e dN4e dN5e dN6e dN7e dN8e];
end