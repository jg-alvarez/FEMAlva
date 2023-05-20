%% CIV 2118 - Método dos Elementos Finitos - 2022.2
% Trabalho Final - Parte 1
% Aluno: João Guilherme M. Alvarez & 
% Matricula: 2220784 & 
%
% Objeto do Elemento isoparamétrico quadrilateral Q8 em estado plano de tensões.
%
%% Class definition
classdef Barra < handle
    %% Public, tunable properties
    properties (Access = public)
        E = 0;
        ni = 0;

        e = 0;          %Epessura
        b = 0;          %Base
        h = 0;          %Altura

        nNE = 0;        %Numero de nós por elementos
        nElem = 0;      %Numero de elementos
        Elem = [];      %Vetor de elementos
        GLE = [];       %Matriz com graus de liberdade de cada elemento

        GLn = [];       %Matriz com Graus de Liberdade por nó
%         nNos = 0;       %Numero de nós total
        NosElem = [];   %Vetor de nós por elemento

        k = [];         %Matriz de rigidez da barra
        d = [];         %Vetor de deslocamento nodal

        tested = [];    %Vetor para verificação do calculo dos deslocamentos, pode ser apagado posteriormente
        defelem = [];   %Vetor de matriz de deformações dos elementos
        tenelem = [];   %Vetor de matriz de tensões dos elementos
        sig = [];       %Matriz de tensões xx, yy e xy

        CC = [];        %Condição de contorno
        P = [];         %Condição de carregamento
    end

    %% Pre-computed constants
    properties (Constant)
        nGL = 2;        %Numero de graus de liberdade por nó
    end

    %% Constructor method
    methods
        function B = Barra(E, ni, e, b, h, nNE, Elem, GLElem, GLNos, NosElem, CC, P)
            B.E = E;
            B.ni = ni;

            B.e = e;
            B.b = b;
            B.h = h;

            B.nNE = nNE;
            B.Elem = Elem;
            B.nElem = length(Elem);
            B.GLE = GLElem;

            B.GLn = GLNos;
            B.NosElem = NosElem;

            B.CC = CC;
            B.P = P;

            B.OrderK();
            B.Delta();
            B.DefTenElem();
            B.TenNos();
        end
    end

    %% Private methods
    methods (Access = private)
%--------------------------------------------------------------------------%
% Organiza a matriz k global a partir das matrizes k de cada elemento
        function OrderK(B)
            B.k = zeros(2 * height(B.GLn), 2 * height(B.GLn));

            for n = 1: B.nElem
                for i = 1: length(B.GLE(n,:))
                    for j = 1: length(B.GLE(n,:))
                        B.k(B.GLE(n, i), B.GLE(n, j)) = B.k(B.GLE(n, i), B.GLE(n, j)) + B.Elem(n).ke(i, j);
                    end
                end
            end

        end
        
%--------------------------------------------------------------------------%
% Calcula os deslocamentos de cada nó e as reações nos apoios
        function Delta(B)
            for i = 1: height(B.P)
                GLsup(i, 1) = i;
            end
            for n = 1: height(B.CC)
                for t = 2: 3
                    if B.CC(n, t) == 0
                        GLsup(B.GLn(B.CC(n, 1), t - 1), 1) = 0;
                    end
                end
            end
            GLsup2 = find(GLsup == 0);
            ksup = B.k;
            usup = zeros(height(B.P), 1);
            fsup = B.P;
            for i = height(GLsup2): -1: 1
                usup(GLsup2(i, 1)) = [];
                fsup(GLsup2(i, 1)) = [];
                ksup(GLsup2(i, 1), :) = [];
                ksup(:, GLsup2(i, 1)) = [];
            end
            usup = ksup \ fsup;

            usup2 = ones(height(B.P), 1);            
            for i = 1: height(GLsup2)
                usup2(GLsup2(i)) = 0;
            end
            loc = find(usup2);
            usup2(loc) = usup;

            B.d = usup2;

            %Calculo das reações nos nós de apoio
            for i = 1: height(B.CC)
                for j = 2: 3
                    if B.CC(i, j) == 0
                        B.P(B.GLn(B.CC(i, 1), (j - 1))) = B.k(B.GLn(B.CC(i, 1), (j - 1)), :) * B.d;
                    end
                end
            end

        end

%--------------------------------------------------------------------------%
% Calcula as deformações e tensões em cada nó individualmente
        function DefTenElem(B)
            for i = 1: length(B.Elem)
                deselem(:,:,i) = zeros((B.nNE * B.Elem(i).nGL), 1);
                for n = 1: length(B.GLE(i, :))
                    deselem(n, 1, i) = B.d(B.GLE(i, n), 1);
                end

                B.tested(:,:,i) = deselem(:,:,i);
                B.defelem(:,:,i) = B.Elem(i).B * deselem(:,:,i);
                B.tenelem(:,:,i) = B.Elem(i).D * B.defelem(:,:,i);
            end

        end
        
%--------------------------------------------------------------------------%
% Calcula as tensões totais em cada
        function TenNos(B)
            B.sig = zeros(height(B.GLn), 3);

            for i = 1: height(B.GLn)
                ne = 0;                                                     %Número de elementos que contem o nó

                for n = 1: B.nElem
                    for t = 1: B.nNE
                        if B.NosElem(n, t) == i
                            B.sig(i, 1) = B.sig(i, 1) + B.tenelem(1, 1, n); %Tensão em xx
                            B.sig(i, 2) = B.sig(i, 2) + B.tenelem(2, 1, n); %Tensão em yy
                            B.sig(i, 3) = B.sig(i, 3) + B.tenelem(3, 1, n); %Tensão em xy

                            ne = ne + 1;
                        end
                    end
                end
                B.sig(i, 1) = (B.sig(i, 1) / ne) / 1000;
                B.sig(i, 2) = (B.sig(i, 2) / ne) / 1000;
                B.sig(i, 3) = (B.sig(i, 3) / ne) / 1000;

            end
        end

    end

    %% Public methods
    methods (Access = public)
%--------------------------------------------------------------------------%
        function ClearElement(B)
            B.E = 0;
            B.ni = 0;
            B.e = 0;
            B.b = 0;
            B.h = 0;
            B.nNE = 0;
            B.nElem = 0;
            B.Elem = [];
            B.GLE = [];
            B.GLn = [];
%             B.nNos = 0;
            B.NosElem = [];
            B.k = [];
            B.d = [];
            B.CC = [];
            B.P = [];
        end

    end
end
