function [S] = vdWDFinitialize_InputKernel(S)
% @file    vdWDFinitialize_InputKernel.m
% @brief   This file contains the functions for directly loading the needed model 
%          kernel functions and the spline functions at model energy ratios.
%          The input kernel and spline function are generated by
%          vdWDF_initial_Genkernel. For more information, please see the
%          manual.
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Dion, Max, Henrik Rydberg, Elsebeth Schröder, David C. Langreth, and Bengt I. Lundqvist. 
% "Van der Waals density functional for general geometries." 
% Physical review letters 92, no. 24 (2004): 246401.
% Román-Pérez, Guillermo, and José M. Soler. 
% "Efficient implementation of a van der Waals density functional: application to double-wall carbon nanotubes." 
% Physical review letters 103, no. 9 (2009): 096102.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
%% Initialization, set parameters and grids
    S.vdWDF_Nrpoints = 1024; %% radial points for composing Phi functions in real and reciprocal space
    S.vdWDF_rMax = 100.0; %% max radius in real space and minimum k point 2*pi/r_max in reciprocal space
    S.vdWDF_dr = S.vdWDF_rMax/S.vdWDF_Nrpoints;
    S.vdWDF_dk = 2.0*pi/S.vdWDF_rMax;
    S.vdWDF_Nqs = 20;
    S.vdWDF_qmesh = [
    1.0e-5            0.0449420825586261 0.0975593700991365 0.159162633466142
    0.231286496836006 0.315727667369529  0.414589693721418  0.530335368404141
    0.665848079422965 0.824503639537924  1.010254382520950  1.227727621364570
    1.482340921174910 1.780437058359530  2.129442028133640  2.538050036534580
    3.016440085356680 3.576529545442460  4.232271035198720  5.0];
    S.vdWDF_qmesh = reshape(S.vdWDF_qmesh', [], 1);
    S.vdWDF_kernel = zeros(1 + S.vdWDF_Nrpoints, S.vdWDF_Nqs, S.vdWDF_Nqs); %% kernal Phi, index 0, reciprocal
    S.vdWDF_d2Phidk2 = zeros(1 + S.vdWDF_Nrpoints, S.vdWDF_Nqs, S.vdWDF_Nqs); %% 2nd derivative of kernal
%% input vdWDF_kernel and vdWDF_d2Phidk2

% load MATLAB computed data
vdWDF_kernel = load('vdWDF_kernel.mat');
vdWDF_d2Phidk2 = load('vdWDF_d2Phidk2.mat');
vdWDF_D2yDx2 = load('vdWDF_D2yDx2.mat');
S.vdWDF_kernel = vdWDF_kernel.vdWDF_kernel;
S.vdWDF_d2Phidk2 = vdWDF_d2Phidk2.vdWDF_d2Phidk2;
S.vdWDF_D2yDx2 = vdWDF_D2yDx2.vdWDF_D2yDx2;