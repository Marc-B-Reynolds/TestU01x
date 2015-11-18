xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                 Starting SmallCrush
                 Version: TestU01 1.2.3
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

smarsa_BirthdaySpacings test:
-----------------------------------------------
   N =  1,  n = 5000000,  r =  0,    d = 1073741824,    t = 2,    p = 1

      Number of cells = d^t = 1152921504606846976
      Lambda = Poisson mean =      27.1051
----------------------------------------------------
Total expected number = N*Lambda      :      27.11
Total observed number                 :     301
p-value of test                       :4.5e-198    *****
-----------------------------------------------

***********************************************************
Test sknuth_Collision calling smultin_Multinomial

***********************************************************

smultin_Multinomial test:
-----------------------------------------------
   N =  1,  n = 5000000,  r =  0,   d = 65536,   t =  2,
       Sparse =   TRUE

       GenerCell = smultin_GenerCellSerial
       Number of cells = d^t =         4294967296
       Expected number per cell =  1 /  858.99346
       EColl = n^2 / (2k) =  2910.383046
       Hashing =   TRUE

       Collision test,    Mu =      2909.2534,    Sigma =     53.896

-----------------------------------------------
Test Results for Collisions

Expected number of collisions = Mu    :     2909.25
Observed number of collisions         :     2833
p-value of test                       :    0.92

-----------------------------
Total number of cells containing j balls

  j =  0                              :       4289970129
  j =  1                              :          4994335
  j =  2                              :             2831
  j =  3                              :                1
  j =  4                              :                0
  j =  5                              :                0

***********************************************************

sknuth_Gap test:
-----------------------------------------------
   N =  1,  n = 200000,  r = 22,   Alpha =        0,   Beta  = 0.00390625
-----------------------------------------------
Number of degrees of freedom          : 1114
Chi-square statistic                  : 1.78e+5
p-value of test                       :   eps      *****

***********************************************************

sknuth_SimpPoker test:
-----------------------------------------------
   N =  1,  n = 400000,  r = 24,   d =   64,   k =   64
-----------------------------------------------
Number of degrees of freedom          :   19
Chi-square statistic                  :  988.96
p-value of test                       :   eps      *****

***********************************************************

sknuth_CouponCollector test:
-----------------------------------------------
   N =  1,  n = 500000,  r = 26,   d =   16
-----------------------------------------------
Number of degrees of freedom          :   44
Chi-square statistic                  :   40.20
p-value of test                       :    0.64

***********************************************************

sknuth_MaxOft test:
-----------------------------------------------
   N =  1,  n = 2000000,  r =  0,   d = 100000,   t =  6

      Number of categories = 100000
      Expected number per category  = 20.00

-----------------------------------------------
Number of degrees of freedom          : 99999
Chi-square statistic                  : 1.04e+5
p-value of test                       :   eps      *****
-----------------------------------------------
Anderson-Darling statistic            : 7.39e-5
p-value of test                       : 1 -  7.4e-5    *****
-----------------------------------------------

***********************************************************

svaria_WeightDistrib test:
-----------------------------------------------
   N =  1,  n = 200000,  r = 27,  k = 256,  Alpha =      0,  Beta =  0.125
-----------------------------------------------
Number of degrees of freedom          :   41
Chi-square statistic                  :   40.42
p-value of test                       :    0.50

-----------------------------------------------

***********************************************************
smarsa_MatrixRank test:
-----------------------------------------------
   N =  1,  n = 20000,  r = 20,    s = 10,    L = 60,    k = 60
-----------------------------------------------
Number of degrees of freedom          :    3
Chi-square statistic                  : 3.76e+6
p-value of test                       :   eps      *****

***********************************************************

sstring_HammingIndep test:
-----------------------------------------------
   N =  1,  n = 500000,  r = 20,   s = 10,   L = 300,   d = 0

Counters with expected numbers >= 10
-----------------------------------------------
Number of degrees of freedom          : 2209
Chi-square statistic                  : 2253.75
p-value of test                       :    0.25

***********************************************************

swalk_RandomWalk1 test:
-----------------------------------------------
   N =  1,  n = 1000000,  r =  0,   s = 30,   L0 =  150,   L1 =  150
-----------------------------------------------
Test on the values of the Statistic H

Number of degrees of freedom          :   52
ChiSquare statistic                   :  221.12
p-value of test                       :   eps      *****
-----------------------------------------------
Test on the values of the Statistic M

Number of degrees of freedom          :   52
ChiSquare statistic                   :  110.69
p-value of test                       :  4.0e-6    *****
-----------------------------------------------
Test on the values of the Statistic J

Number of degrees of freedom          :   75
ChiSquare statistic                   :  103.38
p-value of test                       :    0.02
-----------------------------------------------
Test on the values of the Statistic R

Number of degrees of freedom          :   44
ChiSquare statistic                   :   36.06
p-value of test                       :    0.80
-----------------------------------------------
Test on the values of the Statistic C

Number of degrees of freedom          :   26
ChiSquare statistic                   :   33.03
p-value of test                       :    0.16
-----------------------------------------------


========= Summary results of SmallCrush =========

 Version:          TestU01 1.2.3
 Generator:        xorshift
 Number of statistics:  15
 Total CPU time:   00:00:20.77
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  BirthdaySpacings              4.5e-198
  3  Gap                              eps  
  4  SimpPoker                        eps  
  6  MaxOft                           eps  
  6  MaxOft AD                      1 -  7.4e-5
  8  MatrixRank                       eps  
 10  RandomWalk1 H                    eps  
 10  RandomWalk1 M                   4.0e-6
 ----------------------------------------------
 All other tests were passed
