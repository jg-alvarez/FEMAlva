%% CIV 2118 - Método dos Elementos Finitos - 2022.2
% Trabalho Final - Parte 1
% Aluno: João Guilherme M. Alvarez & Camila Alves
% Matricula: 2220784 & 
%
% Função para o cálculo da matrix B a partir da matriz jacobiana e a derivada das funções de forma.
%
%%
function [B] = CreateB(J, dFF)
     J1 = inv(J);

     tB = J1 * dFF; %Matriz que contem termos de B


     B = [tB(1,1) 0 tB(1,2) 0 tB(1,3) 0 tB(1,4) 0 tB(1,5) 0 tB(1,6) 0 tB(1,7) 0 tB(1,8) 0;...
          0 tB(2,1) 0 tB(2,2) 0 tB(2,3) 0 tB(2,4) 0 tB(2,5) 0 tB(2,6) 0 tB(2,7) 0 tB(2,8);...
          tB(2,1) tB(1,1) tB(2,2) tB(1,2) tB(2,3) tB(1,3) tB(2,4) tB(1,4) tB(2,5) tB(1,5) tB(2,6) tB(1,6) tB(2,7) tB(1,7) tB(2,8) tB(1,8)];

end