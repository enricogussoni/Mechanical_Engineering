%% STRAIN
clear all;
close all;

prova1_ax=-[23 16 35 51 60 63 62 57 58 56 27];
prova1_cir=-[-115 -150 -165 -170 -171 -170 -168 -169 -170 -171 -169];
prova2_ax=-[23 22 37 52 62 63 63 59 59 55 27];
prova2_cir=-[-116 -150 -164 -164 -170 -168 -168 -171 -173 -174 -172];
distance=[6.5 11.5 16.5 21.5 26.5 31.5 36.5 41.5 46.5 51.5 67.5];

EPSexp_ax=(prova1_ax+prova2_ax)/2;
EPSexp_cir=(prova1_cir+prova2_cir)/2;

xdis_shell=[0
2.0125
3.01873
4.02498
5.03117
6.03692
7.04142
8.04382
9.04354
10.0422
11.037
12.0298
13.0224
14.0158
16.0069
17.0041
19.0066
21.013
23.0126
25.0101
27.0088
28.0086
29.0084
30.0083
31.0082
32.0081
33.008
34.0079
35.0078
36.0077
37.0076
38.0075
39.0074
40.0073
41.0072
42.0071
43.007
44.0069
45.0068
46.0067
47.0066
48.0065
49.0064
50.0063
51.0062
52.0061
53.006
54.0059
55.0058
56.0057
57.0055
58.0054
59.0053
60.0052
61.0051
62.005
63.0049
64.0048
65.0047
66.0046
67.0045
68.0044
69.0043
70.0042
71.0041
72.004
73.0039
74.0038
75.0037
76.0036
77.0035
78.0034
79.0033
80.0032
81.0031
82.003
83.0029
84.0028
85.0027]';

EPSshell_cir=10^6*[4.24516E-006
2.01288E-005
3.30231E-005
4.7168E-005
6.19992E-005
7.70496E-005
9.19392E-005
0.000106369
0.000120123
0.00013304
0.000145003
0.000155967
0.000165924
0.000174876
0.000189844
0.000195938
0.000205662
0.000212465
0.000216888
0.000219484
0.000220733
0.000220982
0.000221042
0.00022095
0.000220741
0.000220446
0.00022009
0.000219695
0.00021928
0.000218862
0.000218453
0.000218065
0.000217706
0.000217382
0.000217101
0.000216865
0.000216677
0.000216541
0.000216457
0.000216427
0.00021645
0.000216527
0.000216656
0.000216835
0.000217063
0.000217336
0.00021765
0.000218001
0.000218381
0.000218782
0.000219196
0.000219609
0.000220008
0.000220376
0.000220692
0.000220935
0.000221075
0.000221085
0.000220928
0.000220566
0.000219956
0.000219052
0.000217801
0.000216148
0.000214036
0.000211403
0.000208187
0.000204324
0.000199753
0.000194414
0.000188254
0.000181228
0.000173297
0.000164453
0.000154669
0.000144018
0.000132418
0.000120475
0.000107261]';

EPSshell_ax=10^6*[-0.000274459
-0.000192458
-0.000147469
-0.000110002
-7.93214E-005
-5.47024E-005
-3.54391E-005
-2.08465E-005
-1.02656E-005
-3.08584E-006
1.25055E-006
3.24295E-006
3.32714E-006
1.89166E-006
-4.23819E-006
-8.39093E-006
-1.78215E-005
-2.7631E-005
-3.6931E-005
-4.52369E-005
-5.2326E-005
-5.5388E-005
-5.81311E-005
-6.05669E-005
-6.27111E-005
-6.45828E-005
-6.62034E-005
-6.75951E-005
-6.87809E-005
-6.97834E-005
-7.06245E-005
-7.13247E-005
-7.19034E-005
-7.23779E-005
-7.27638E-005
-7.30746E-005
-7.33213E-005
-7.35131E-005
-7.36565E-005
-7.37559E-005
-7.38131E-005
-7.38278E-005
-7.37971E-005
-7.3716E-005
-7.3577E-005
-7.33705E-005
-7.30846E-005
-7.27056E-005
-7.22176E-005
-7.1603E-005
-7.0843E-005
-6.99172E-005
-6.88045E-005
-6.74836E-005
-6.59332E-005
-6.41331E-005
-6.20645E-005
-5.97114E-005
-5.70618E-005
-5.41085E-005
-5.08509E-005
-4.72969E-005
-4.34642E-005
-3.93831E-005
-3.5098E-005
-3.06707E-005
-2.61825E-005
-2.17374E-005
-1.74645E-005
-1.35222E-005
-1.00989E-005
-7.4203E-006
-5.7435E-006
-5.37534E-006
-6.64026E-006
-9.94766E-006
-1.56299E-005
-2.43895E-005
-3.63355E-005]';

xdis_solid=[0
2
4
6
8
10
11.9767
13.9535
15.9302
17.907
19.8837
21.8605
24.4461
27.0316
29.0084
30.9851
32.9619
34.9386
36.9154
38.8921
40.8689
42.8456
44.8223
46.7991
48.7758
50.7526
52.7293
54.7061
56.6828
58.6596
60.6363
62.613
64.5898
66.5665
68.5433
70.52
72.4968
74.4735
76.4503
78.427
80.4037
82.3805
84.3572
86.334
88.3107
90.2875
92.2642
94.241]'-10;

xdis_solid_ax=[0
2
4
6
8
10
11.9767
13.9535
15.9302
17.907
19.8837
21.8605
23.8372
25.8139
27.7907
29.7674
31.7442
33.7209
35.6977
37.6744
39.6512
41.6279
43.6046
45.5814
47.5581
49.5349
51.5116
53.4884
55.4651
57.4419
59.4186
61.3953
63.3721
65.3488
67.3256
69.3023
71.2791
73.2558
75.2326
77.2093
79.186
81.1628
83.1395
85.1163
87.093
89.0698
91.0465
93.0233
95
96.9767
98.9535
100.93
102.907]'-10;

EPSsolid_ax=10^6*[-5.86354E-005
-5.80254E-005
-6.04915E-005
-7.45254E-005
-8.86844E-005
-8.44247E-005
-6.57525E-005
-3.62044E-005
-1.49273E-005
-9.17198E-006
-7.43939E-006
-1.12233E-005
-1.65249E-005
-2.35448E-005
-3.036E-005
-3.68212E-005
-4.24997E-005
-4.70833E-005
-5.08624E-005
-5.35286E-005
-5.56121E-005
-5.68617E-005
-5.77945E-005
-5.82232E-005
-5.85493E-005
-5.86358E-005
-5.87492E-005
-5.8784E-005
-5.88872E-005
-5.89655E-005
-5.90702E-005
-5.91015E-005
-5.90352E-005
-5.87499E-005
-5.81641E-005
-5.71289E-005
-5.55336E-005
-5.32159E-005
-5.00923E-005
-4.60409E-005
-4.11132E-005
-3.53379E-005
-2.9075E-005
-2.27173E-005
-1.71547E-005
-1.32048E-005
-1.30363E-005
-1.91493E-005
-2.93584E-005
-3.95575E-005
-4.56532E-005
-4.54478E-005
-4.14545E-005]';

EPSsolid_cir=10^6*[1.54862E-006
2.91109E-006
5.74954E-006
1.00844E-005
1.80411E-005
3.15848E-005
5.1224E-005
7.45788E-005
9.85426E-005
0.000120517
0.000139294
0.000154727
0.000174676
0.000175673
0.000181917
0.000185978
0.000188332
0.000189412
0.000189605
0.000189216
0.000188503
0.00018764
0.00018677
0.000185974
0.000185312
0.000184811
0.00018449
0.000184352
0.000184397
0.000184622
0.00018502
0.000185582
0.000186289
0.000187113
0.000188006
0.000188893
0.000189665
0.000190165
0.000190187
0.000189464
0.000187668
0.000184415
0.000179283
0.000171832
0.000161677
0.000148508
0.000132216
0.00011321]';

figure(1)
hold on;
xlim([-10 100]);
plot(distance,EPSexp_ax,'r-o');
plot(xdis_shell,EPSshell_ax);
plot(xdis_solid_ax,EPSsolid_ax,'g');
legend('Experimental','Shell','Solid')
xlabel('Distance (mm)');
ylabel('Strain (microeps)');
title('Axial Strain')

figure(2)
hold on;
plot(distance,EPSexp_cir,'-ro')
plot(xdis_shell,EPSshell_cir)
plot(xdis_solid,EPSsolid_cir,'g')
legend('Experimental','Shell','Solid')
xlabel('Distance (mm)');
ylabel('Strain (microeps)');
title('Circumferential Strain')

%% STRESS
E=206000;
poss=0.3;
den=(poss+1)*(2*poss-1);

SIGexp_cir=E/den*((poss-1)*EPSexp_cir-poss*EPSexp_ax)/10^6;
SIGexp_ax=E/den*((poss-1)*EPSexp_ax-poss*EPSexp_cir)/10^6;

SIGshell_cir=[-19.9564
-10.0289
-3.61598
2.51225
8.28141
13.6388
18.5505
22.9995
26.9863
30.5201
33.6165
36.303
38.6113
40.5713
43.5639
44.6558
46.1844
47.0088
47.3218
47.2881
47.036
46.8602
46.6647
46.4576
46.2458
46.0347
45.8287
45.6312
45.445
45.2718
45.1132
44.97
44.8428
44.7318
44.6372
44.5589
44.4968
44.4507
44.4204
44.4059
44.4069
44.4234
44.4555
44.5032
44.5664
44.6453
44.7398
44.8497
44.9748
45.1146
45.2681
45.4343
45.6115
45.7972
45.9887
46.182
46.3724
46.554
46.7199
46.8615
46.9691
47.0311
47.0343
46.9636
46.8022
46.5313
46.1302
45.5764
44.8455
43.9122
42.7491
41.3293
39.6239
37.6073
35.2489
32.5345
29.4193
25.9903
22.0241]';

SIGshell_ax=[-63.1241
-42.9559
-31.5719
-21.8314
-13.6073
-6.76789
-1.1788
3.29545
6.79077
9.43596
11.3511
12.648
13.4271
13.7782
13.503
13.0079
11.5696
9.8209
8.00838
6.28629
4.74272
4.05394
3.42432
2.85423
2.34262
1.88738
1.48557
1.13371
0.827967
0.564316
0.338716
0.147208
-0.0139784
-0.148346
-0.259064
-0.348911
-0.420245
-0.474967
-0.514504
-0.539784
-0.551235
-0.548795
-0.53189
-0.49945
-0.449938
-0.381363
-0.291303
-0.176943
-0.035127
0.137583
0.344835
0.590399
0.878055
1.21146
1.59401
2.02864
2.5176
3.06227
3.66283
4.31796
5.02452
5.7771
6.56768
7.38508
8.21455
9.03717
9.82937
10.5623
11.2013
11.7054
12.0268
12.1101
11.8927
11.3031
10.2643
8.68716
6.48862
3.55257
-0.21717]';

SIGsolid_cir=[-5.37207
-4.38462
-3.56041
-3.31648
-2.02737
1.08844
7.0577
14.0618
20.9286
26.8068
31.0439
34.4569
36.5527
38.3351
39.2439
39.6962
39.8217
39.7032
39.4833
39.1651
38.8633
38.5448
38.2883
38.0515
37.8832
37.7459
37.6682
37.6232
37.6281
37.6672
37.7519
37.8748
38.0433
38.2522
38.5027
38.7833
39.0824
39.3703
39.6131
39.7482
39.7096
39.3896
38.6844
37.4377
35.5309
32.7791
29.0009
24.2755]';

SIGsolid_ax=[-14.6658
-13.7973
-13.7942
-16.3907
-18.5234
-16.8463
-11.037
-3.08677
3.35231
6.70851
8.29761
8.72507
8.32731
7.43429
6.31462
5.13429
4.01609
3.02211
2.1889
1.52008
1.0086
0.633009
0.370117
0.194002
0.0813959
0.0119418
-0.030177
-0.0557892
-0.0708303
-0.0757802
-0.0663129
-0.0329025
0.0383309
0.165597
0.370075
0.674976
1.10289
1.6732
2.39712
3.2725
4.27588
5.35469
6.41484
7.30697
7.83481
7.73451
6.50596
3.77754
0.0157607
-3.74464
-6.46886
-7.69042
-7.781];

figure(3)
hold on;
xlim([-10 100])
plot(distance,SIGexp_ax,'r-o');
plot(xdis_shell,SIGshell_ax,'y');
plot(xdis_solid_ax,SIGsolid_ax,'g');
line([0 85],[22 22])
line([0 85],[0 0],'color',[.1 .0 .0])
legend('Experimental','Shell','Solid','Mariotte','Expected value')
xlabel('Distance (mm)');
ylabel('Stress (MPa)');
title('Axial Stress')

figure(4)
hold on;
plot(distance,SIGexp_cir,'-ro')
plot(xdis_shell,SIGshell_cir,'y')
plot(xdis_solid,SIGsolid_cir,'g')
line([0 85],[44 44])
legend('Experimental','Shell','Solid','Mariotte')
xlabel('Distance (mm)');
ylabel('Stress (MPa)');
title('Circumferential Stress')

%% MISES
solid_mis=[-5.37207
-4.38462
-3.56041
-3.31648
-2.02737
1.08844
7.0577
14.0618
20.9286
26.8068
31.0439
34.4569
36.5527
38.3351
39.2439
39.6962
39.8217
39.7032
39.4833
39.1651
38.8633
38.5448
38.2883
38.0515
37.8832
37.7459
37.6682
37.6232
37.6281
37.6672
37.7519
37.8748
38.0433
38.2522
38.5027
38.7833
39.0824
39.3703
39.6131
39.7482
39.7096
39.3896
38.6844
37.4377
35.5309
32.7791
29.0009
24.2755];