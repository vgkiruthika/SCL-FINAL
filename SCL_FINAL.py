# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 21:32:19 2023

@author: KIRUTHIKA GOPAL
"""

#WS6

#NEWTON FORWARD INTERPOLATION


import matplotlib.pyplot as plt
from sympy import *
import numpy as np

def create_table(x,y,n,table):
  for col in range(1,n):
    for row in range(n-col):
      table[row][col]=round(table[row+1][col-1]-table[row][col-1],2)

  return table

def calculate(u,n):
  temp = u
  for i in range(1, n):
      temp = temp * (u - i)
  print(temp)
  return temp

def fact(n):
    f = 1;
    for i in range(2, n + 1):
        f *= i
    return f

x_val=list(map(float,input("Enter x values:").split()))
y_val=list(map(float,input("Enter y values:").split()))
x,y=symbols('x y')

n=len(x_val)

table=[[0.0 for i in range(n)]for i in range(n)]
for i in range(n):
    table[i][0]=y_val[i]

#creating forward difference
table=create_table(x_val,y_val,n,table)

for row in table:
  print(row)


inp=float(input("Enter the value to be interpolated at:"))

y=table[0][0]
u=(x-x_val[0])/(x_val[1]-x_val[0])

for i in range(1,n):
  y = y + (((calculate(u, i)) * (table[0][i])) / fact(i));


pol=simplify(y)
res=pol.subs(x,inp)
print("Result: ",res)

x_pts=np.linspace(min(x_val),max(x_val),100)
y_pts=[]


for pt in x_pts:
  y_pts.append(pol.subs(x,pt))

x_val.append(inp)
y_val.append(res)

ax,fig=plt.subplots()
plt.plot(x_pts,y_pts)
plt.scatter(x_val,y_val,color='red')
first_diff=diff(pol)
res1=first_diff.subs(x,inp)
print("Res_1: ",res1)

second_diff=diff(first_diff)
res2=second_diff.subs(x,inp)
print("Res_2: ",res2)

#40 50 60 70 80
#31 73 124 159 190

#100 150 200 250 300 350 400
#10.63 13.03 15.04 16.81 18.42 19.90 21.27
#160

#0.96 0.98 1.00 1.02 1.04
#0.7825 0.7739 0.7651 0.7563 0.7473
#0.96

#45 50 55 60
#0.7071 0.7660 0.8192 0.8660
#52

#--------------------------------------------------------------------------------

#Newton's backward interpolation

import matplotlib.pyplot as plt
from sympy import *
import numpy as np

def create_table(x, y, n, table):
    for i in range(1, n):
      for j in range(n - 1, i - 1, -1):
        table[j][i] =table[j][i - 1] - table[j - 1][i - 1]
    return table

def calculate(u, n):
    temp = u
    for i in range(1, n):
        temp = temp * (u + i)
    return temp

def fact(n):
    f = 1
    for i in range(2, n + 1):
        f *= i
    return f

x_val = list(map(float, input("Enter x values: ").split()))
y_val = list(map(float, input("Enter y values: ").split()))
x, y = symbols('x y')

n = len(x_val)

table = [[0.0 for i in range(n)] for i in range(n)]
for i in range(n):
    table[i][0] = y_val[i]

# Create backward difference table
table = create_table(x_val, y_val, n, table)

for row in table:
    print(row)

inp = float(input("Enter the value to be interpolated at: "))

y = table[n-1][0]
u = (x - x_val[n - 1]) / (x_val[1] - x_val[0])

for i in range(1, n):
    y = y + (((calculate(u, i)) * (table[n-1][i])) / fact(i))

pol = simplify(y)
res = pol.subs(x, inp)
print("Result: ", res)

x_pts = np.linspace(min(x_val), max(x_val), 100)
y_pts = []

for pt in x_pts:
    y_pts.append(pol.subs(x, pt))

x_val.append(inp)
y_val.append(res)

fig, ax = plt.subplots()
plt.plot(x_pts, y_pts)
plt.scatter(x_val, y_val, color='red')

first_diff = diff(pol)
res1 = first_diff.subs(x, inp)
print("Res_1: ", res1)

second_diff = diff(first_diff)
res2 = second_diff.subs(x, inp)
print("Res_2: ", res2)


#1891 1901 1911 1921 1931
#46 66 81 93 101
#1925
#96.836

#---------------------------------------------------------------------------

#Trapezoidal rule

from sympy import*

def f(fun,val):
  return fun.subs(x,val)

def trapezoidal(x0,xn,n,fun):
  h=(xn-x0)/n

  lis = [x0 + i * h for i in range(n + 1)]
  n_lis=len(lis)
  #print(lis)

  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  for i in range(1,n_lis-1):
    sum2=sum2+f(fun,lis[i])

  sum2=2*sum2

  integ=(h/2)*(sum1+sum2)
  return integ


x=Symbol('x')
inp=input("Enter the function to be integrated: ")
eq=Eq(eval(inp),0)
fun=eq.lhs
lower_limit=float(input("Enter the lower limit: "))
upper_limit=float(input("Enter the upper limit: "))
sub_interval=int(input("Enter the number of sub-intervals:"))

result=trapezoidal(lower_limit,upper_limit,sub_interval,fun)
print(result)

#1/(1 + x**2)
#0
#1
#output:0.784240766617816

#----------------------------------------------------------------------

#Simpsons_1_3 rule

from sympy import*

def f(fun,val):
  return fun.subs(x,val)

def simpson_3(x0,xn,n,fun):
  h=(xn-x0)/n
  lis = [round(x0 + i * h,2) for i in range(n + 1)]
  n_lis=len(lis)
  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  sum3=0
  for i in range(1,n_lis-1):
    if(i%2==0): #even ordinates as indexing starts from 0
      sum2=sum2+f(fun,lis[i])
    else: #odd ordinates as indexing starts from 0
      sum3=sum3+f(fun,lis[i])

  sum2=2*sum2
  sum3=4*sum3
  integ=(h/3)*(sum1+sum2+sum3)
  return integ

x=Symbol('x')
inp=input("Enter the function to be integrated: ")
eq=Eq(eval(inp),0)
fun=eq.lhs
lower_limit=float(input("Enter the lower limit: "))
upper_limit=float(input("Enter the upper limit: "))
sub_interval=int(input("Enter the number of sub-intervals:"))

result=simpson_3(lower_limit,upper_limit,sub_interval,fun)
print(result)

#exp(-x)*sqrt(x)
#0.5
#0.7
#4
#output:0.0848271185743080

#-----------------------------------------------------------------------

#Simpsons_3_8 rule

from sympy import*

def f(fun,val):
  return fun.subs(x,val)


def simpson_3_8(x0,xn,n,fun):
  h=(xn-x0)/n
  lis = [round(x0 + i * h,2) for i in range(n + 1)]
  #print(lis)
  n_lis=len(lis)
  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  sum3=0
  for i in range(1,n_lis-1):
    if(i==3):
      continue
    #print(f(fun,lis[i]))
    sum2=sum2+f(fun,lis[i])

  for i in range(3,n-3+1,3):
    sum3=sum3+f(fun,lis[i])

  sum2=3*sum2
  sum3=2*sum3

  integ=((3*h)/8)*(sum1+sum2+sum3)
  return integ

x=Symbol('x')
inp=input("Enter the function to be integrated: ")
eq=Eq(eval(inp),0)
fun=eq.lhs
lower_limit=float(input("Enter the lower limit: "))
upper_limit=float(input("Enter the upper limit: "))
sub_interval=int(input("Enter the number of sub-intervals:"))

result=simpson_3_8(lower_limit,upper_limit,sub_interval,fun)
print(result)

#1/(1+x**2)
#0
#6
#6
#output: 1.35708083649260


#-----------------------------------------------------------------------

"""When a train is moving at 30 metres per second stream is shut off and
brakes are applied. The speed of the train (V ) in metres per second
after t seconds is given by
t 0 5 10 15 20 25 30 35 40
V 30 24 19.5 16 13.6 11.7 10.0 8.5 7.0
Using Simpson’s rule and Trapezoidal rule determine the distance moved
by the train in 40 secs."""

#Question 3) Problem sheet:

import matplotlib.pyplot as plt
from sympy import *
import numpy as np


def f(fun,val):
  return fun.subs(x,val)

def trapezoidal(x0,xn,n,fun):
  h=(xn-x0)/n

  lis = [x0 + i * h for i in range(n + 1)]
  n_lis=len(lis)
  #print(lis)

  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  for i in range(1,n_lis-1):
    sum2=sum2+f(fun,lis[i])

  sum2=2*sum2

  integ=(h/2)*(sum1+sum2)
  return integ

def simpson_3(x0,xn,n,fun):
  h=(xn-x0)/n
  lis = [round(x0 + i * h,2) for i in range(n + 1)]
  n_lis=len(lis)
  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  sum3=0
  for i in range(1,n_lis-1):
    if(i%2==0): #even ordinates as indexing starts from 0
      sum2=sum2+f(fun,lis[i])
    else: #odd ordinates as indexing starts from 0
      sum3=sum3+f(fun,lis[i])

  sum2=2*sum2
  sum3=4*sum3
  integ=(h/3)*(sum1+sum2+sum3)
  return integ

def simpson_3_8(x0,xn,n,fun):
  h=(xn-x0)/n
  lis = [round(x0 + i * h,2) for i in range(n + 1)]
  #print(lis)
  n_lis=len(lis)
  sum1=f(fun,x0)+f(fun,xn)
  sum2=0
  sum3=0
  for i in range(1,n_lis-1):
    if(i==3):
      continue
    #print(f(fun,lis[i]))
    sum2=sum2+f(fun,lis[i])

  for i in range(3,n-3+1,3):
    sum3=sum3+f(fun,lis[i])

  sum2=3*sum2
  sum3=2*sum3

  integ=((3*h)/8)*(sum1+sum2+sum3)
  return integ


def create_table(x,y,n,table):
  for col in range(1,n):
    for row in range(n-col):
      table[row][col]=round(table[row+1][col-1]-table[row][col-1],2)

  return table

def calculate(u,n):
  temp = u
  for i in range(1, n):
      temp = temp * (u - i)
  return temp

def fact(n):
    f = 1;
    for i in range(2, n + 1):
        f *= i
    return f

x_val=list(map(float,input("Enter x values:").split()))
y_val=list(map(float,input("Enter y values:").split()))
x,y=symbols('x y')

n=len(x_val)

table=[[0.0 for i in range(n)]for i in range(n)]
for i in range(n):
    table[i][0]=y_val[i]

#creating forward difference
table=create_table(x_val,y_val,n,table)

#for row in table:
  #print(row)


y=table[0][0]
u=(x-x_val[0])/(x_val[1]-x_val[0])

for i in range(1,n):
  y = y + (((calculate(u, i)) * (table[0][i])) / fact(i));

pol=simplify(y)

lower_limit=min(x_val)
upper_limit=max(x_val)
sub_interval=len(x_val)-1
result1=trapezoidal(lower_limit,upper_limit,sub_interval,pol)
print("Trapezoidal : ", result1)
result2=simpson_3(lower_limit,upper_limit,sub_interval,pol)
print("Simpson 1/3 : ", result2)
result3=simpson_3_8(lower_limit,upper_limit,sub_interval,pol)
print("Simpson 3/8 : ", result3)


#0 5 10 15 20 25 30 35 40
#30 24 19.5 16 13.6 11.7 10.0 8.5 7.0


#-------------------------------------------------------------------------
#----------------------------------------------------------------------

#WS7

#EULER
import math


def differential_equation(x, y):
    return (x**2 - y**2 )/(x**2 + y**2)


x0 = float(input("Enter the initial x-value : "))
y0 = float(input("Enter the initial y-value : "))
x_target = float(input("Enter the target x-value : "))


h = 0.1


x = x0
y = y0


while x < x_target:
    y_next = y + h * differential_equation(x, y)
    x += h
    y = y_next


print(f'y({x_target}) ≈ {y:.6f}')




MODIFIED EULER


def differential_equation(x, y):
    return x**2 - y + 2*x + 4*y


x0 = float(input("Enter the initial x-value (e.g., 0.3): "))
y0 = float(input("Enter the initial y-value (e.g., -0.18): "))
x_target = float(input("Enter the target x-value (e.g., 2.5): "))


h = 0.1


x = x0
y = y0


while x < x_target:

    y_pred = y + h * differential_equation(x, y)
    

    y_next = y + 0.5 * h * (differential_equation(x, y) + differential_equation(x + h, y_pred))
    
    x += h
    y = y_next


print(f'y({x_target}) ≈ {y:.6f}')



#------------------------------------------------------------------

#RUNGE KUTTA THIRD ORDER


def differential_equation(x, y):
    return (x * y) / (x**2 + y**2 - x * y)


x0 = 0
y0 = 10


x_target = 3


h = 0.1


x = x0
y = y0


while x < x_target:
    k1 = h * differential_equation(x, y)
    k2 = h * differential_equation(x + 0.5 * h, y + 0.5 * k1)
    k3 = h * differential_equation(x + h, y - k1 + 2 * k2)
    
    y_next = y + (1/6) * (k1 + 4 * k2 + k3)
    x += h
    y = y_next

print(f'y({x_target}) ≈ {y:.6f}')

#-------------------------------------------------------------------------------

#RUNGE KUTTA FOURTH ORDER

# Define the differential equation
def differential_equation(x, y):
    return x * y + x - y**3

# Initial conditions
x0 = 1
y0 = -4

# Target x-value
x_target = 2.5

# Step size
h = 0.1

# Initialize variables
x = x0
y = y0

# Runge-Kutta method of fourth order to approximate the solution
while x < x_target:
    k1 = h * differential_equation(x, y)
    k2 = h * differential_equation(x + 0.5 * h, y + 0.5 * k1)
    k3 = h * differential_equation(x + 0.5 * h, y + 0.5 * k2)
    k4 = h * differential_equation(x + h, y + k3)
    
    y_next = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
    x += h
    y = y_next

# Print the result
print(f'y({x_target}) ≈ {y:.6f}')


#--------------------------------------------------------------------
#--------------------------------------------------------------------


#WS8


# 1. Milne’s predictor corrector method and Adam Bashforth method

import numpy as np
import math

def f(x,y):
    return ((1 + x**2)*y**2)/2

h = 0.1

x0 = 0
y0 = 1

x1 = 0.1
y1 = 1.06

x2 = x0 + h
y2 = 1.12

x3 = x1 + h
y3 = 1.21

x4 = x2 + h

#Adams Predictor formula
yp4 = y3 + h*(55*f(x3,y3) - 59*f(x2,y2) + 37*f(x1,y1) - 9*f(x0,y0))/24
x4 = x3 + h
fp4 = f(x4,yp4)
print('Adams Predictor formula')
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)

#Adams Corrector formula
yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*fp4 )/24
f4 = f(x4,yc4)
print( '\nAdams Corrector formula')
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )

yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*f4 )/24
print ('\nRefined y = ',yc4)

print("-----------------------------")



#Milne Predictor formula
yp4 = y0 + 4*h*(2*f(x1,y1) - f(x2,y2) + 2*f(x3,y3))/3
x4 = x3 + h
fp4 = f(x4,yp4)
print("Milne’s predictor corrector method")
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)

#Simpson Corrector formula
yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + fp4)/3
f4 = f(x4,yc4)
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )


yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + f4)/3
print ('\n Refined y = ',yc4)

print("-----------------------------")


#----------------------------------------------------------------------

#2. using runge kutta third order

def runge_kutta_third_order(f, x0, y0, h, num_steps):
    x = x0
    y = y0

    for _ in range(num_steps):
        k1 = h * f(x, y)
        k2 = h * f(x + h / 2, y + k1 / 2)
        k3 = h * f(x + h, y - k1 + 2 * k2)

        y = y + (k1 + 4 * k2 + k3) / 6
        x = x + h

    return x, y

def f(x, y):
    return x*y + x**2

x0 = 0
y0 = 1
h = 0.1

xVals = [0.1, 0.2, 0.3]
yVals = []
for i in xVals:
    num_steps = int((i - x0) / h)
    xf, yf = runge_kutta_third_order(f, x0, y0, h, num_steps)
    yVals.append(yf)


x1, x2, x3 = 0.1, 0.2, 0.3
y1, y2, y3 = yVals[0], yVals[1], yVals[2]

x4 = x2 + h

#Adams Predictor formula
yp4 = y3 + h*(55*f(x3,y3) - 59*f(x2,y2) + 37*f(x1,y1) - 9*f(x0,y0))/24
x4 = x3 + h
fp4 = f(x4,yp4)
print('Adams Predictor formula')
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)
#Adams Corrector formula
yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*fp4 )/24
f4 = f(x4,yc4)
print( '\nAdams Corrector formula')
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )

yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*f4 )/24
print ('\nRefined y = ',yc4)

print("-----------------------------")
#Milne Predictor formula
yp4 = y0 + 4*h*(2*f(x1,y1) - f(x2,y2) + 2*f(x3,y3))/3
x4 = x3 + h
fp4 = f(x4,yp4)
print("Milne’s predictor corrector method")
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)
#Simpson Corrector formula
yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + fp4)/3
f4 = f(x4,yc4)
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )


yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + f4)/3
print ('\n Refined y = ',yc4)

print("-----------------------------")


#-------------------------------------------------------------------


# 3.using kutta fourth order
def runge_kutta_fourth_order(f, x0, y0, h, num_steps):
    x = x0
    y = y0

    for _ in range(num_steps):
        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        x = x + h

    return x, y

def f(x, y):
    return x*y + x**2

x0 = 0
y0 = 1
h = 0.1

xVals = [0.1, 0.2, 0.3]
yVals = []
for i in xVals:
    num_steps = int((i - x0) / h)
    xf, yf = runge_kutta_third_order(f, x0, y0, h, num_steps)
    yVals.append(yf)


x1, x2, x3 = 0.1, 0.2, 0.3
y1, y2, y3 = yVals[0], yVals[1], yVals[2]

x4 = x2 + h

#Adams Predictor formula
yp4 = y3 + h*(55*f(x3,y3) - 59*f(x2,y2) + 37*f(x1,y1) - 9*f(x0,y0))/24
x4 = x3 + h
fp4 = f(x4,yp4)
print('Adams Predictor formula')
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)
#Adams Corrector formula
yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*fp4 )/24
f4 = f(x4,yc4)
print( '\nAdams Corrector formula')
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )

yc4 = y3 + h*( f(x1,y1) - 5*f(x2,y2) + 19*f(x3,y3) + 9*f4 )/24
print ('\nRefined y = ',yc4)

print("-----------------------------")
#Milne Predictor formula
yp4 = y0 + 4*h*(2*f(x1,y1) - f(x2,y2) + 2*f(x3,y3))/3
x4 = x3 + h
fp4 = f(x4,yp4)
print("Milne’s predictor corrector method")
print( 'y predicted Value = ',yp4 )
print ('f(y predicted)  = ',fp4)
#Simpson Corrector formula
yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + fp4)/3
f4 = f(x4,yc4)
print( 'y Corrected Value = ',yc4)
print( 'f(y Corrected) = ',f4 )


yc4 = y2 + h*( f(x2,y2) + 4*f(x3,y3) + f4)/3
print ('\n Refined y = ',yc4)

print("-----------------------------")


#---------------------------------------------------------------



#4. euler
def f(x, y):
    return (2 - y**2) / (5 * x)

def euler(x0, y, h, x):
    temp = 0
    while x0 < x:
        temp = y
        y = y + h * f(x0, y)
        x0 = x0 + h
    return y0


def milne_predictor_corrector(x0, y0, h, xn):
    y1 = modified_euler(x0, y0, h, x0 + h)
    y2 = modified_euler(x0, y0, h, x0 + 2*h)
    y3 = modified_euler(x0, y0, h, x0 + 3*h)
    while x0 < xn:
        y4_predict = y0 + (4*h/3) * (2*f(x0, y0) - f(x0 + h, y1) + 2*f(x0 + 2*h, y2))
        y4_correct = y2 + (h/3) * (f(x0 + 2*h, y2) + 4*f(x0 + h, y1) + f(x0, y0))
        y0, y1, y2, y3 = y1, y2, y3, y4_correct
        x0 = x0 + h
    return y4_correct

def adam_bashforth(x0, y0, h, xn):
    y1 = modified_euler(x0, y0, h, x0 + h)
    y2 = modified_euler(x0, y0, h, x0 + 2*h)
    y3 = modified_euler(x0, y0, h, x0 + 3*h)
    while x0 < xn:
        y4_predict = y3 + (h/24) * (55*f(x0 + 3*h, y3) - 59*f(x0 + 2*h, y2) + 37*f(x0 + h, y1) - 9*f(x0, y0))
        y4_correct = y3 + (h/24) * (9*f(x0 + h, y4_predict) + 19*f(x0 + 3*h, y3) - 5*f(x0 + 2*h, y2) + f(x0, y0))
        y0, y1, y2, y3 = y1, y2, y3, y4_correct
        x0 = x0 + h
    return y4_correct

y0 = 1
x0 = 4
h = 0.1
xn1 = 4.1
xn2 = 4.2
xn3 = 4.3
xn4 = 4.4

y1 = euler(x0, y0, h, xn1)
y2 = euler(x0, y0, h, xn2)
y3 = euler(x0, y0, h, xn3)
y4_milne = milne_predictor_corrector(x0, y0, h, xn4)
y4_adam = adam_bashforth(x0, y0, h, xn4)

print(f"y({xn1}) = {y1}")
print(f"y({xn2}) = {y2}")
print(f"y({xn3}) = {y3}")
print(f"y({xn4}) by Milne's predictor corrector method = {y4_milne}")
print(f"y({xn4}) by Adam Bashforth method = {y4_adam}")


#----------------------------------------------------------------------------


#5. modified euler

def f(x, y):
    return (2 - y**2) / (5 * x)

def modified_euler(x0, y0, h, xn):
    while x0 < xn:
        y0 = y0 + h * f(x0 + h/2, y0 + (h/2) * f(x0, y0))
        x0 = x0 + h
    return y0

def milne_predictor_corrector(x0, y0, h, xn):
    y1 = modified_euler(x0, y0, h, x0 + h)
    y2 = modified_euler(x0, y0, h, x0 + 2*h)
    y3 = modified_euler(x0, y0, h, x0 + 3*h)
    while x0 < xn:
        y4_predict = y0 + (4*h/3) * (2*f(x0, y0) - f(x0 + h, y1) + 2*f(x0 + 2*h, y2))
        y4_correct = y2 + (h/3) * (f(x0 + 2*h, y2) + 4*f(x0 + h, y1) + f(x0, y0))
        y0, y1, y2, y3 = y1, y2, y3, y4_correct
        x0 = x0 + h
    return y4_correct

def adam_bashforth(x0, y0, h, xn):
    y1 = modified_euler(x0, y0, h, x0 + h)
    y2 = modified_euler(x0, y0, h, x0 + 2*h)
    y3 = modified_euler(x0, y0, h, x0 + 3*h)
    while x0 < xn:
        y4_predict = y3 + (h/24) * (55*f(x0 + 3*h, y3) - 59*f(x0 + 2*h, y2) + 37*f(x0 + h, y1) - 9*f(x0, y0))
        y4_correct = y3 + (h/24) * (9*f(x0 + h, y4_predict) + 19*f(x0 + 3*h, y3) - 5*f(x0 + 2*h, y2) + f(x0, y0))
        y0, y1, y2, y3 = y1, y2, y3, y4_correct
        x0 = x0 + h
    return y4_correct

y0 = 1
x0 = 4
h = 0.1
xn1 = 4.1
xn2 = 4.2
xn3 = 4.3
xn4 = 4.4

y1 = modified_euler(x0, y0, h, xn1)
y2 = modified_euler(x0, y0, h, xn2)
y3 = modified_euler(x0, y0, h, xn3)
y4_milne = milne_predictor_corrector(x0, y0, h, xn4)
y4_adam = adam_bashforth(x0, y0, h, xn4)

print(f"y({xn1}) = {y1}")
print(f"y({xn2}) = {y2}")
print(f"y({xn3}) = {y3}")
print(f"y({xn4}) by Milne's predictor corrector method = {y4_milne}")
print(f"y({xn4}) by Adam Bashforth method = {y4_adam}")

#------------------------------------------------------------------
#FORMULA


#EULER yn = yn-1 + hf(x0 + n-1 h,yn-1)
#MODIFIED EULER  y2 = y1 + (h/2)[f(x0 + h,y1)+f(x0 + 2h,y2)]
#RUNGE KUTTA THIRD ORDER

# k1 = hf(x0,y0)
# k2 = hf(x0 + 1/2 h,y0 + 1/2k1)
# k' = hf(x0+h,y0+k1)
# k3 = hf(x0 + h,y0 + k')
# k = 1/6(k1 + 4k2 + k3)

#RUNGE KUTTA FOURTH ORDER

# k1 = hf(x0,y0)
# k2 = hf(x0 + 1/2 h,y0 + 1/2k1)
# k3 = hf(x0 + 1/2h,y0+1/2k2)
# k4 = hf(x0+h,y0+k3)
# k = 1/6(k1+2k2+2k3+k4)
