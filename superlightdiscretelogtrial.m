%trial discrete logarithm
clc;

%for ii=1:100
%dimension
%d= randi([400 500]) 
dim=5

%choose the power
%m=64


%M= [[-317, 937, -963, -705, 260, -82, 143, -74, 265, -135]; [589, 837, 522, 589, -817, 135, 146, -783, -79, 318]; [398, 455, -379, 399, 280, 82, -800, 802, 824, -95]; [-709, 743, 971, -454, -332, 406, 938, 331, -178, -965]; [207, -126, 844, -275, 412, 820, 763, 843, 19, 885]; [-323, -220, -199, 903, -891, -272, -428, -949, -586, 894]; [-490, -914, -520, 422, 829, 262, 981, -1, -115, -29]; [383, 641, -718, -940, 370, -657, 144, -232, -466, 27]; [-20, 445, 223, 469, 686, 614, 653, 978, 7, -789]; [-369, -850, 957, -60, 628, 904, 641, 864, -8, 433]]
%H= [[-166, -542, 599, 559, -520, 7, 355, 589, 906, 399]; [976, -414, 508, -36, 847, -24, 75, 308, -405, -535]; [-261, 653, 100, -932, 849, 36, 331, 998, 304, -417]; [433, 303, 517, 38, -547, 192, 315, -167, -852, -967]; [269, -619, -155, 996, -150, -361, 861, -897, -248, -168]; [-393, 700, 905, -208, 477, -700, -736, -945, 605, 476]; [151, -166, 114, -466, 464, -505, 521, 620, 159, -665]; [-51, 725, -307, -470, -858, 271, 199, 580, 527, 926]; [-904, 594, 243, -328, 562, 938, 377, -3, -251, -545]; [-871, 707, 463, 228, 73, 218, -37, -602, -403, -471]]
M=GenerateRandomMatrix(5,5,-1000,1000)
H=GenerateRandomMatrix(5,5,-1000,1000)


m=randi([20 200])
n= randi([20 200])

%choose matrices M and H

I=maxeins(dim)
V=add(otimes(M,add(I,H)),H)
%=[2 10 0 -2;2 -3 4 4;-8 -4 -6 12;-2 2 -10 -5]
%H=generatespecialmatrixbalik(d);
%H=generatespecialmatrix(d);
%H= [9 -10 -9 -7 -8;-9 -12 1 18 11;7 7 3 -30 -50;3 7 7 16 10;-20 9 3 18 -100]
%H=maxplussprand(d,d,100)
%H=generatespecialmatrixnocheck(d);
%H=GenerateRandomMatrix(d,d,1,1000);


% V= [259   501   778   233   922  -892  -131  -344   473   -28;
%   -316   737   166   -21   -81   -81    62  -139   689    79;
%   -490   130  -517   187  -450  -611   969   564   565    79;
%   -323  -505   930   471   747   156   288  -995  -281   690;
%   -867   700  -938   302   324  -688  -908   -35  -957   353;
%   -889   980  -460  -616  -488  -648  -771   575   249    44;
%    212   368   280  -965  -522  -101  -775  -988   101  -655;
%    421   296  -607   131   375  -252   710   729   -36  -205;
%   -338    17   152  -849   142  -538   749  -667   571  -781;
%    866  -994   590  -408  -937   974  -459    75  -936  -745];
% 
% H= [378  -578  -354   463   978    17  -636   205  -563  -916;
%     63   547   743  -851  -975   240  -268   562  -248  -369;
%    584  -801   891   923   954   117  -354   280  -297   122;
%   -729   306   642  -798  -716    76    98  -261   425    24;
%   -140   174  -617  -715   914  -636  -690   796  -745  -535;
%     67    97  -506    74   500  -715    14  -673  -293   434;
%   -979   591  -233   558  -345   240   454  -668   447  -895;
%   -568  -268    64  -734   719   649  -653   260   217  -201;
%    172  -948  -432     1  -534  -265   348   957  -613   645;
%   -557  -824   377   115  -748   284  -792   311  -265   653]

%disp('Now generating an instance to which discrete log will be applied');
B=otimes(V,powmaxplus(H,n))
A=otimes(V,powmaxplus(H,m))
%disp('Now applying the discrete log');
nlog=superlightdiscretelogavfk(B,V,H)
mlog=superlightdiscretelogavfk(A,V,H)
d;




%end





