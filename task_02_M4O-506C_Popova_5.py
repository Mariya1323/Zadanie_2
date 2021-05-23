import requests
import math
import scipy.special
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
import requests as rqst
import json
import csv

# Рассчитать ЭПР
# Построить график
# Сохранить результаты в файл

# Radar cross section
def RCS(lam, r):
    summ = 0
    kr = 2 * math.pi * r / lam
    # Задаем значения функций Бесселя для n = 0 для первой итерации
    J_prev = spherical_jn (0, kr)
    Y_prev = spherical_yn (0, kr)
    H_prev = J_prev + 1j * Y_prev
    for n in range(1, 50):
        # Вычисляем значения функций Бесселя для текущей n
        J_now = spherical_jn (n, kr)
        #J_prev = spherical_jn (n - 1, kr)
        Y_now = spherical_yn (n, kr)
        #Y_prev = spherical_yn (n - 1, kr)
        H_now = J_now + 1j * Y_now
        #H_prev = J_prev + 1j * Y_prev
        # Считаем коэффициенты a и b
        a = J_now / H_now
        b = (kr * J_prev - n * J_now) / (kr * H_prev - n * H_now)
        summ += ((-1) ** n) * (n + 0.5) * (b - a)
        # Переносим значения функций Бесселя на следующий шаг
        J_prev = J_now
        Y_prev = Y_now
        H_prev = H_now
    return lam * lam * np.abs(summ) * np.abs(summ) / math.pi

def graf(lam, p):
  plt.plot(lam, p)
  plt.ylabel('RCS')
  plt.xlabel('lambda')
  plt.grid()
  plt.show()

def graf_freq(f, p):
  plt.plot(f, p)
  plt.ylabel('RCS, [м2]')
  plt.xlabel('f, [Гц]')
  plt.grid()
  plt.show()

# Скачать файл с вариантом задания 

def download(url):
  r=rqst.get(url)
  return r.text

# Разобрать прочитанные данные и найти нужные данные

# Возврат данных по номеру варианта
def var(text,nomervar):
  t=text.splitlines()
  return t[nomervar]

# Запись в txt-файл
def txtfile(ls, s):
    with open('task_02_M40-506C_05.txt', 'w') as txt:
        txt.write('f, [Гц]{0:4}RSC, [м2]\n'.format(' '))
        for i in range(len(ls)):
            txt.write(str('{0:.0f}{1:4}{2}\n'.format(ls[i], ' ', s[i])))
    
  
if __name__ == '__main__':
    txt = download('https://jenyay.net/uploads/Student/Modelling/task_02_01.txt')
    nomervar = 5
    line = var(txt, nomervar - 1)
    print(line)
    L = line.split(';')
    for i in range(len(L)):
        L[i] = float(L[i].split('=')[1])
    D, fmin, fmax = L[:]
    s = []
    ls = np.linspace(fmin,fmax,500)
    for f in ls:
        p=RCS(3e8 / f, D/2)
        s.append(p)
    
    graf_freq(ls,s)
    txtfile(ls, s)
