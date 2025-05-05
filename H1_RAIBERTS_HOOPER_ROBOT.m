%% Script for the kinematic model of the Raibert's hopper
clear; clc;

%% Definition of symbolic variables
syms m l d             % System parameters
syms theta psi I       % Generalized coordinates

%% Definition of the Pfaffian constraint A(q)*qdot = 0
A = [ I + m*(l + d)^2,  m*(l + d)^2,  0 ];
disp('matrix A');
disp(A);

%% Calculation of the null space of A
% The null space of A, i.e., g1 and g2
G = null(A);
disp('null space of A   G(q): ');
disp(G);
% G is a 3x matrix

%% Definition of inputs and kinematic model
% We introduce two inputs u1 and u2.
syms u1 u2 
u = [u1; u2];
disp('input u: ');
disp(u);

qdot = G * u;
disp('qdot=G(q)*u');
disp(qdot);

g1 = G(:,1);
g2 = G(:,2);

% Calculation of the Lie brackets [g1, g2]
J_g1 = jacobian(g1, [theta, psi, I]);  % Jacobian of g1 with respect to [theta, psi, I]
J_g2 = jacobian(g2, [theta, psi, I]);  % Jacobian of g2 with respect to [theta, psi, I]

% Lie bracket
Lie_g1_g2 = J_g2 * g1 - J_g1 * g2;
F = [g1, g2, Lie_g1_g2];
disp('F=span(g1,g2,[g1,g2])');
disp(F)
disp('rank(F)=dim(deltaA(q))=v=');
disp(rank(F))
%% Check if the system is holonomic

disp('qdot=G(q)u is controllable');
disp('The system is subject to nonholonomic constraints only');
disp('dim(deltaA(q))=3=v=n');
disp('A(q)*qdot=0 is not integrable, and it is a set of nonholonomic constraints in the Pfaffian form');