# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Dictionary.ipynb ,
#   this module will construct a Python dictionary
#   of Butcher tables for a large number of explicit
#   RK methods.

# Authors: Brandon Clark
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#	   Gabriel M Steward

# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python

# Step 2a: Generating a Dictionary of Butcher Tables for Explicit Runge Kutta Techniques

# Initialize the dictionary Butcher_dict
Butcher_dict = {}

# Step 2.a.i: Euler's Method

Butcher_dict['Euler'] = (
[[sp.sympify(0)],
 ["", sp.sympify(1)]]
, 1)

# Step 2.a.ii: RK2 Heun's Method

Butcher_dict['RK2 Heun'] = (
[[sp.sympify(0)],
[sp.sympify(1), sp.sympify(1)],
["", sp.Rational(1,2), sp.Rational(1,2)]]
, 2)

# Step 2.a.iii: RK2 Midpoint (MP) Method

Butcher_dict['RK2 MP'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
["", sp.sympify(0), sp.sympify(1)]]
, 2)

# Step 2.a.iv: RK2 Ralston's Method

Butcher_dict['RK2 Ralston'] = (
[[sp.sympify(0)],
[sp.Rational(2,3), sp.Rational(2,3)],
["", sp.Rational(1,4), sp.Rational(3,4)]]
, 2)

# Step 2.a.v: Kutta's  Third-order Method

Butcher_dict['RK3'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.sympify(1), sp.sympify(-1), sp.sympify(2)],
["", sp.Rational(1,6), sp.Rational(2,3), sp.Rational(1,6)]]
, 3)

# Step 2.a.vi: RK3 Heun's Method

Butcher_dict['RK3 Heun'] = (
[[sp.sympify(0)],
[sp.Rational(1,3), sp.Rational(1,3)],
[sp.Rational(2,3), sp.sympify(0), sp.Rational(2,3)],
["", sp.Rational(1,4), sp.sympify(0), sp.Rational(3,4)]]
, 3)

# Step 2.a.vii: RK3 Ralton's Method

Butcher_dict['RK3 Ralston'] = ([[0],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(3,4), sp.sympify(0), sp.Rational(3,4)],
["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)]]
, 3)

# Step 2.a.viii: Strong Stability Preserving Runge-Kutta (SSPRK3) Method
Butcher_dict['SSPRK3'] = (
[[0],
[sp.sympify(1), sp.sympify(1)],
[sp.Rational(1,2), sp.Rational(1,4), sp.Rational(1,4)],
["", sp.Rational(1,6), sp.Rational(1,6), sp.Rational(2,3)]]
, 3)

# Step 2.a.ix: Classic RK4 Method

Butcher_dict['RK4'] = ([[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(1,2), sp.sympify(0), sp.Rational(1,2)],
[sp.sympify(1), sp.sympify(0), sp.sympify(0), sp.sympify(1)],
["", sp.Rational(1,6), sp.Rational(1,3), sp.Rational(1,3), sp.Rational(1,6)]]
, 4)

# Step 2.a.x:  RK5 Dormand-Prince Method

Butcher_dict['DP5'] = ([[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
[sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
[sp.sympify(1), sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
[sp.sympify(1), sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
["", sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), sp.sympify(0)]]
, 5)

# Step 2.a.xi:  RK5 Dormand-Prince Method Alternative

Butcher_dict['DP5alt'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
["", sp.Rational(821, 10800), sp.sympify(0), sp.Rational(19683, 71825), sp.Rational(175273, 912600), sp.Rational(395, 3672), sp.Rational(785, 2704), sp.Rational(3, 50)]]
, 5)

# Step 2.a.xii:  RK5 Dormand-Prince Method Alternative
Butcher_dict['DP5alt'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
["", sp.Rational(821, 10800), sp.sympify(0), sp.Rational(19683, 71825), sp.Rational(175273, 912600), sp.Rational(395, 3672), sp.Rational(785, 2704), sp.Rational(3, 50)]]
, 5)

# Step 2.a.xiii:  RK5 Cash-Karp Method

Butcher_dict['CK5'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
[sp.sympify(1), sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
[sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
["",sp.Rational(37,378), sp.sympify(0), sp.Rational(250,621), sp.Rational(125,594), sp.sympify(0), sp.Rational(512,1771)]]
, 5)

# Step 2.a.xiv:  RK6 Dormand-Prince Method

Butcher_dict['DP6'] = (
[[0],
[sp.Rational(1,10), sp.Rational(1,10)],
[sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
[sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
[sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
[sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
[sp.sympify(1), sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
[sp.sympify(1), sp.Rational(465467, 266112), sp.Rational(-2945, 1232), sp.Rational(-5610201, 14158144), sp.Rational(10513573, 3212352), sp.Rational(-424325, 205632), sp.Rational(376225, 454272), sp.sympify(0)],
["", sp.Rational(61, 864), sp.sympify(0), sp.Rational(98415, 321776), sp.Rational(16807, 146016), sp.Rational(1375, 7344), sp.Rational(1375, 5408), sp.Rational(-37, 1120), sp.Rational(1,10)]]
, 6)

#  Step 2.a.xv:  RK6 Luther's Method

q = sp.sqrt(21)
Butcher_dict['L6'] = (
[[0],
[sp.sympify(1), sp.sympify(1)],
[sp.Rational(1,2), sp.Rational(3,8), sp.Rational(1,8)],
[sp.Rational(2,3), sp.Rational(8,27), sp.Rational(2,27), sp.Rational(8,27)],
[(7 - q)/14, (-21 + 9*q)/392, (-56 + 8*q)/392, (336 -48*q)/392, (-63 + 3*q)/392],
[(7 + q)/14, (-1155 - 255*q)/1960, (-280 -  40*q)/1960, (-320*q)/1960, (63 + 363*q)/1960, (2352 + 392*q)/1960],
[sp.sympify(1), ( 330 + 105*q)/180, sp.Rational(2,3), (-200 + 280*q)/180, (126 - 189*q)/180, (-686 - 126*q)/180, (490 -  70*q)/180],
["", sp.Rational(1, 20), sp.sympify(0), sp.Rational(16, 45), sp.sympify(0), sp.Rational(49, 180), sp.Rational(49, 180), sp.Rational(1, 20)]]
, 6)

# Step 2.a.xvi: RK8 Dormand-Prince Method

Butcher_dict['DP8']=(
[[0],
[sp.Rational(1, 18), sp.Rational(1, 18)],
[sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
[sp.Rational(1, 8), sp.Rational(1, 32), sp.sympify(0), sp.Rational(3, 32)],
[sp.Rational(5, 16), sp.Rational(5, 16), sp.sympify(0), sp.Rational(-75, 64), sp.Rational(75, 64)],
[sp.Rational(3, 8), sp.Rational(3, 80), sp.sympify(0), sp.sympify(0), sp.Rational(3, 16), sp.Rational(3, 20)],
[sp.Rational(59, 400), sp.Rational(29443841, 614563906), sp.sympify(0), sp.sympify(0), sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
[sp.Rational(93, 200), sp.Rational(16016141, 946692911), sp.sympify(0), sp.sympify(0), sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
[sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), sp.sympify(0), sp.sympify(0), sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
[sp.Rational(13, 20), sp.Rational(246121993, 1340847787), sp.sympify(0), sp.sympify(0), sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
[sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), sp.sympify(0), sp.sympify(0), sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
[sp.sympify(1), sp.Rational(185892177, 718116043), sp.sympify(0), sp.sympify(0), sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
[sp.sympify(1), sp.Rational(403863854, 491063109), sp.sympify(0), sp.sympify(0), sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), sp.sympify(0)],
["", sp.Rational(14005451, 335480064), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)]]
, 8)

# Step 3.a:  Generating a Dictionary of Butcher Tables for Explicit Runge Kutta Techniques

# Step 3.a.i: Adaptive Heun-Euler Method

Butcher_dict['AHE'] = (
[[sp.sympify(0)],
[sp.sympify(1), sp.sympify(1)],
["", sp.Rational(1,2), sp.Rational(1,2)],
["", sp.sympify(1), sp.sympify(0)]]
, 2)

# Step 3.a.ii: Adaptive Bogacki-Shampine Method

Butcher_dict['ABS'] = (
[[sp.sympify(0)],
[sp.Rational(1,2), sp.Rational(1,2)],
[sp.Rational(3,4), sp.sympify(0), sp.Rational(3,4)],
[sp.sympify(1), sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)],
["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9), sp.sympify(0)],
["", sp.Rational(7,24), sp.Rational(1,4), sp.Rational(1,3), sp.Rational(1,8)]]
, 3)

# Step 3.a.iii: Adaptive Runge-Kutta-Fehlberg

Butcher_dict['ARKF'] = (
[[sp.sympify(0)],
[sp.Rational(1,4), sp.Rational(1,4)],
[sp.Rational(3,8), sp.Rational(3,32), sp.Rational(9,32)],
[sp.Rational(12,13), sp.Rational(1932,2197), sp.Rational(-7200,2197), sp.Rational(7296,2197)],
[sp.sympify(1), sp.Rational(439,216), sp.sympify(-8), sp.Rational(3680,513), sp.Rational(-845,4104)],
[sp.Rational(1,2), sp.Rational(-8,27), sp.sympify(2), sp.Rational(-3544,2565), sp.Rational(1859,4104), sp.Rational(-11,40)],
["", sp.Rational(16,135), sp.sympify(0), sp.Rational(6656,12825), sp.Rational(28561,56430), sp.Rational(-9,50), sp.Rational(2,55)],
["", sp.Rational(25,216), sp.sympify(0), sp.Rational(1408,2565), sp.Rational(2197,4104), sp.Rational(-1,5), sp.sympify(0)]]
, 5)

# Step 3.a.iv: Adaptive Cash-Karp

Butcher_dict['ACK'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
[sp.sympify(1), sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
[sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
["",sp.Rational(37,378), sp.sympify(0), sp.Rational(250,621), sp.Rational(125,594), sp.sympify(0), sp.Rational(512,1771)],
["",sp.Rational(2825,27648), sp.sympify(0), sp.Rational(18575,48384), sp.Rational(13525,55296), sp.Rational(277,14336), sp.Rational(1,4)]]
, 5)

# Step 3.a.v: Adaptive Dormand-Prince 5(4)

Butcher_dict['ADP5'] = (
[[0],
[sp.Rational(1,5), sp.Rational(1,5)],
[sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
[sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
[sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
[sp.sympify(1), sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
[sp.sympify(1), sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
["", sp.Rational(35,384), sp.sympify(0), sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), sp.sympify(0)],
["", sp.Rational(5179,57600), sp.sympify(0), sp.Rational(7571,16695), sp.Rational(393,640), sp.Rational(-92097,339200), sp.Rational(187,2100), sp.Rational(1,40)]]
, 5)

# Step 3.a.vi: Adaptive Dormand-Prince 8(7)

Butcher_dict['ADP8']=(
[[0],
[sp.Rational(1, 18), sp.Rational(1, 18)],
[sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
[sp.Rational(1, 8), sp.Rational(1, 32), sp.sympify(0), sp.Rational(3, 32)],
[sp.Rational(5, 16), sp.Rational(5, 16), sp.sympify(0), sp.Rational(-75, 64), sp.Rational(75, 64)],
[sp.Rational(3, 8), sp.Rational(3, 80), sp.sympify(0), sp.sympify(0), sp.Rational(3, 16), sp.Rational(3, 20)],
[sp.Rational(59, 400), sp.Rational(29443841, 614563906), sp.sympify(0), sp.sympify(0), sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
[sp.Rational(93, 200), sp.Rational(16016141, 946692911), sp.sympify(0), sp.sympify(0), sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
[sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), sp.sympify(0), sp.sympify(0), sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
[sp.Rational(13, 20), sp.Rational(246121993, 1340847787), sp.sympify(0), sp.sympify(0), sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
[sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), sp.sympify(0), sp.sympify(0), sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
[sp.sympify(1), sp.Rational(185892177, 718116043), sp.sympify(0), sp.sympify(0), sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
[sp.sympify(1), sp.Rational(403863854, 491063109), sp.sympify(0), sp.sympify(0), sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), sp.sympify(0)],
["", sp.Rational(14005451, 335480064), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)],
["", sp.Rational(13451932, 455176623), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.sympify(0), sp.Rational(-808719846, 976000145), sp.Rational(1757004468, 5645159321), sp.Rational(656045339, 265891186), sp.Rational(-3867574721, 1518517206), sp.Rational(465885868, 322736535), sp.Rational(53011238, 667516719), sp.Rational(2, 45), sp.sympify(0)]]
, 8)

# Step 4: The Adams-Bashforth Method

import numpy as np
from sympy import symbols
from sympy import integrate

# okay so our goal here is to calculate the coefficients to arbitrary precision for
# Adams-Bashforth methods.
# Careful! AB methods are NOT structured the same as RK methods, the same code
# cannot read them! 

# this does not work to generate an order 1 method. But why would you need to generate
# that, it's just a 1x1 matrix with a 1 in it.
order = 19 # change this if you want different orders. Default 19, which takes several seconds to generate.
pythonButcher = [[sp.Rational(0.0/1.0),sp.Rational(0.0/1.0)],
                 [sp.Rational(0.0/1.0),sp.Rational(0.0/1.0)]]
# just to initialize it, we'll set it in a minute

n = 1
while n < order+1: 
    j = 0
    while j < n: # the number of coefficients in each method is equal to the order.

        # set up the product
        x = symbols('x')
        expr = x
        i = 0
        while i < n:
            if i == j:
                i = i+1
            else:
                expr = expr*(x + i)
                i = i+1
        expr = expr/x
        expr2 = integrate(expr,(x,0,1))
        expr2 = expr2 * sp.Rational((-1.0)**j,(sp.factorial(j)* sp.factorial(n - j - 1)))
        if (len(pythonButcher) < n):
            pythonButcher.append([sp.Rational(0.0/1.0),sp.Rational(0.0/1.0)])
            # add another row if we need it.
        if (len(pythonButcher[n-1]) < j+1):
            pythonButcher[n-1].append(expr2)
            i = 0
            while i < n:
                # Uncomment this section to fill the matrix with zeroes
                # if(len(pythonButcher[i]) < len(pythonButcher[n-1])):
                    # pythonButcher[i].append(sp.sympify(0))
                i = i+1
        else:
            pythonButcher[n-1][j] = expr2
        j = j+1
    n = n+1

# Due to the way we set this up to make sure everything was explicitly an array,
# We need to remove a single trailing zero. 
# Comment out if you want the zeroes.
pythonButcher[0] = [1]

Butcher_dict['AB']=(
pythonButcher
, order)
