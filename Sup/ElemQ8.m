%% CIV 2118 - Método dos Elementos Finitos - 2022.2
% Trabalho Final - Parte 1
% Aluno: João Guilherme M. Alvarez & Camila Alves
% Matricula: 2220784 & 
%
% Objeto do Elemento isoparamétrico quadrilateral Q8 em estado plano de tensões.
%
%% Class definition
classdef ElemQ8 < handle
    %% Public, tunable properties
    properties (Access = public)
        E = 0;          %MPa
        ni = 0;
        h = 0;          %mm
        
        Gauss = [];     %Matriz com posições e pesos dos pontos de gauss
        coord = [];
        g = [];         %Condições de contorno
        
        Nf = [];        %Funções de forma
        D = [];
        Jac = [];       %Matriz Jacobiana
        J = 0;          %Jacobiano
        B = [];         %Matriz de transformação

        ke = [];        %Matriz de rigidez do elemento
    end

    %% Pre-computed constants
    properties (Constant)
        n = 8;          %Numero de nós no elemento
        nGL = 2;        %Numero de graus de liberdade por nó
    end

    %% Constructor method
    methods
        function Q8 = ElemQ8(E, ni, h, nGauss, conect, geom)
            if nargin == 0
                Q8.E = 200 * 10^6;
                Q8.ni = 0.3;
                Q8.h = 0.001;
                nGauss = 3;
                conect = [1 2 3 4 5 6 7 8];
                geom = [0 0; 1 0; 1 1; 0 1; 0.5 0; 1 0.5; 0.5 1; 0 0.5];
            else
                Q8.E = E;
                Q8.ni = ni;
                Q8.h = h;
            end

            Q8.defGauss(nGauss);
            Q8.StartupGL(geom, conect);
            Q8.formD();
            Q8.formK(nGauss);
        end
    end

    %% Private methods
    methods (Access = private)
        %Inicializa o objeto com suas coordenadas reais e graus de liberdade
        function StartupGL(Q8, geom, conect)
%             i = 0;
            Q8.coord = zeros(Q8.n, Q8.nGL);
            
            for k = 1: Q8.n
                for j = 1: Q8.nGL
                    Q8.coord(k, j) = geom(conect(k), j);

%                     i = i + 1;
%                     Q8.g(i) = nf(conect(1, k), j);
                end
            end
        end

        %Gera a matriz com os valores resultantes do metodo de Gauss
        function defGauss(Q8, ngp)
            Q8.Gauss = gauss(ngp);
        end

        %Calcula a matriz D
        function formD(Q8)
            Q8.D = CreateD(Q8.E, Q8.ni);
        end

        %Calcula a matriz de rigidez do elemento Q8
        function formK(Q8, nGauss)
            k = zeros((Q8.n * Q8.nGL), (Q8.n * Q8.nGL));

            for i = 1: nGauss
                for j = 1: nGauss
                    [Q8.Nf, dff] = FunForm(Q8.Gauss, i, j);
                    Q8.Jac = dff * Q8.coord;
                    Q8.J = det(Q8.Jac);
                    
                    Q8.B = CreateB(Q8.Jac, dff);
    
                    wi = Q8.Gauss(i, 2);
                    wj = Q8.Gauss(j, 2);
                    
                    k = k + ((Q8.B' * Q8.D * Q8.B) * wi * wj * Q8.J * Q8.h);
                end
            end

            Q8.ke = k;
        end

    end

    %% Public methods
    methods (Access = public)
        function ClearElement(Q8)
            Q8.E = 0;
            Q8.ni = 0;
            Q8.h = 0;
            Q8.Gauss = [];
            Q8.coord = [];
            Q8.g = [];
            Q8.Nf = [];
            Q8.D = [];
            Q8.Jac = [];
            Q8.J = 0;
            Q8.B = [];
            Q8.ke = [];
        end

    end
end
