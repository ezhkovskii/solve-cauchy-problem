import math

from scipy.misc import derivative
import sympy


def power_series(x0, x, n, init_cond, der_one_func):
    """
    Составление степенного ряда 
    
    Args:
        x0 - в какой точке находим производную
        x - точка x, для нахождение y(x) - значения степенного ряда
        n - размер степенного ряда (до какой производной строим степ. ряд)
        init_cond - начальные условия задачи Коши
        der_one_func - функция первой производной

    Return:
        Значение y(x) - решение степенного ряда
    """

    # Значения перед иксами в степенном ряде
    values = []

    # y(x0)
    values.append(init_cond(x0))
    # y'(x0)
    values.append(der_one_func(x0))

    # вычисления значений следующих производных
    for der in range(2, n+1):
        order = der + 1 if der % 2 == 0 else der + 2
        value = derivative(func=der_one_func, x0=x0, n=der, order=order)
        values.append(value)
    
    # Построение степенного ряда
    res = 0
    for i in range(len(values)):
        res += values[i] * x**i / math.factorial(i)

    return res


def successive_approximation_method(x0, x, n, init_cond):
    """
    Метод последовательных приближений 
    
    Args:
        x0 - в какой точке находим производную
        x - точка x, для нахождение y(x)
        n - размер приближения
        init_cond - функция начальные условия задачи Коши

    Return:
        Значение y(x) - решение
    """

    X= sympy.Symbol('X')
    y0 = yN = init_cond()
    for i in range(1, n+1):        
        yN = y0 + sympy.integrate(yN, (X, x0, X))

    return yN.subs(X, x)


def euler_method(der_one_func, X0, Y0, step, left=None, right=None, length_segment=None):
    """
    Метод Эйлера 
    
    Args:
        der_one_func - функция первой производной
        X0, Y0 - начальные условия
        step - шаг
        left - начало отрезка
        right - конец отрезка
        length_segment - длина отрезка

    Return:
        список y для построения графика
    """
    y_list = []
    x_list = []
    y = Y0
    x = X0

    if length_segment is None and left is not None and right is not None:
        length_segment = int((right - left) / step)

    for i in range(length_segment+1):
        y += step * der_one_func(x)
        x += step
        y_list.append(y)
        x_list.append(x)
    
    return y_list, x_list


def runge_kutte_method(der_one_func, X0, Y0, step, left=None, right=None, length_segment=None):
    """
    Метод Рунге - Кутты
    
    Args:
        der_one_func - функция первой производной
        X0, Y0 - начальные условия
        step - шаг
        left - начало отрезка
        right - конец отрезка
        length_segment - длина отрезка

    Return:
        список y для построения графика
    """
    y_list = []
    x_list = []
    y = Y0
    x = X0

    if length_segment is None and left is not None and right is not None:
        length_segment = int((right - left) / step)

    for i in range(length_segment+1):
        k1 = der_one_func(x)
        k2 = der_one_func(x + step / 2)
        k3 = der_one_func(x + step / 2)
        k4 = der_one_func(x + step)
        y_delta = step / 6 * (k1 + 2*k2 + 2*k3 + k4)
        y += y_delta
        x += step
        y_list.append(y)
        x_list.append(x)
    
    return y_list, x_list
    
def runge_rule(y, y_2, p):
    """
    Оценка точности численного решения ОДУ правилом Рунге

    Args:
        y - список игриков полученных с шагом h
        y_2 - список игриков полученных с шагом h/2
        p - порядок точности использованного численного метода

    Return:
        порядок точности
    """
    numerator = max([abs(i1 - i2) for i1, i2 in zip(y, y_2)])
    denominator = 2**p - 1

    return numerator / denominator


    
