%%%
%%% Defines commonly-used constant parameters as global variables.
%%%

%%% Physical constants
m1km = 1000;
t1day = 86400;
t1year = 365*t1day;
Omega = 2*pi*366/365/t1day;    
g = 9.81;

%%% Output file name definitions
OUTN_PV = 'PV';
OUTN_PSI = 'PSI';
OUTN_TRACER = 'TRACER';
OUTN_TFILE = 'time.txt';

%%% File system
exec_name = 'DIYsimulate.exe';
model_code_dir = 'model_code';
