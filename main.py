import math

import matplotlib.pyplot as plt

from methods import power_series, successive_approximation_method, euler_method, runge_kutte_method, runge_rule


STEP = 0.1
LEFT = 0
RIGHT = 1

X0 = 0
Y0 = 1

# Делать ли оценку погрешности по правилу Рунге для каждого метода
ESTIMATION_ERROR = True 


def der_one_func(x):
    """
    Функция y первой производной
    """
    return 2 + math.exp(x) * math.sin(x)

def func_init_cond(x0=X0):
    """
    Начальное условие
    """
    return Y0


if __name__=="__main__":
    # Массивы Y для вывода графиков
    y_power_series_2 = []
    y_power_series_4 = []
    y_sam_1 = []
    y_sam_2 = []

    # Вычисление массива X  
    length_segment = int((RIGHT - LEFT) / STEP)
    array_x = [ i * STEP  for i in range(length_segment+1) ] # todo: np.arange

    # Решение задачи Коши несколькими методами
    for step in array_x:
        y_power_series_2.append(power_series(X0, step, 2, func_init_cond, der_one_func))
        y_power_series_4.append(power_series(X0, step, 4, func_init_cond, der_one_func))
        y_sam_1.append(successive_approximation_method(X0, step, 1, func_init_cond))
        y_sam_2.append(successive_approximation_method(X0, step, 2, func_init_cond))
    y_euler_1, x_list_1 = euler_method(der_one_func, X0, Y0, STEP, length_segment=length_segment)
    y_euler_2, x_list_2 = euler_method(der_one_func, X0, Y0, STEP/2, left=LEFT, right=RIGHT)
    y_r_k_1, x_list_1 = runge_kutte_method(der_one_func, X0, Y0, STEP, length_segment=length_segment)
    y_r_k_2, x_list_2 = runge_kutte_method(der_one_func, X0, Y0, STEP/2, left=LEFT, right=RIGHT)
    
    # Рисование графиков
    plt.plot(array_x, y_power_series_2, label='Степенной ряд, n = 2')
    plt.plot(array_x, y_power_series_4, label='Степенной ряд, n = 4')
    plt.plot(array_x, y_sam_1, label='Метод последовательных приближений, n = 1')
    plt.plot(array_x, y_sam_2, label='Метод последовательных приближений, n = 2')
    plt.plot(x_list_1, y_euler_1, label='Метод Эйлера, с шагом h = 0.1')
    plt.plot(x_list_2, y_euler_2, label='Метод Эйлера, с шагом h = 0.05')
    plt.plot(x_list_1, y_r_k_1, label='Метод Рунге-Кутты, с шагом h = 0.1')
    plt.plot(x_list_2, y_r_k_2, label='Метод Рунге-Кутты, с шагом h = 0.05')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Задание 4')
    plt.legend()
    plt.show()

    if ESTIMATION_ERROR:
        STEP = STEP / 2

        length_segment = int((RIGHT - LEFT) / STEP)
        array_x = [ i * STEP for i in range(length_segment+1) ]

        y_power_series_2_runge = []
        y_power_series_4_runge = []
        y_sam_1_runge = []
        y_sam_2_runge = []

        for step in array_x:
            y_power_series_2_runge.append(power_series(X0, step, 2, func_init_cond, der_one_func))
            y_power_series_4_runge.append(power_series(X0, step, 4, func_init_cond, der_one_func))
            y_sam_1_runge.append(successive_approximation_method(X0, step, 1, func_init_cond))
            y_sam_2_runge.append(successive_approximation_method(X0, step, 2, func_init_cond))
        y_euler_1_runge, x_list_1 = euler_method(der_one_func, X0, Y0, STEP, length_segment=length_segment)
        y_euler_2_runge, x_list_2 = euler_method(der_one_func, X0, Y0, STEP/2, left=LEFT, right=RIGHT)
        y_r_k_1_runge, x_list_1 = runge_kutte_method(der_one_func, X0, Y0, STEP, length_segment=length_segment)
        y_r_k_2_runge, x_list_2 = runge_kutte_method(der_one_func, X0, Y0, STEP/2, left=LEFT, right=RIGHT)

        print('Оценка погрешности по правилу Рунге')
        print('\nСтепенной ряд, n = 2')
        print(runge_rule(y_power_series_2, y_power_series_2_runge, 2))
        print('\nСтепенной ряд, n = 4')
        print(runge_rule(y_power_series_4, y_power_series_4_runge, 4))
        print('\nМетод последовательных приближений, n = 1')
        print(runge_rule(y_sam_1, y_sam_1_runge, 1))
        print('\nМетод последовательных приближений, n = 2')
        print(runge_rule(y_sam_2, y_sam_2_runge, 2))
        print('\nМетод Эйлера, с шагом h = 0.1')
        print(runge_rule(y_euler_1, y_euler_1_runge, 1))
        print('\nМетод Эйлера, с шагом h = 0.05')
        print(runge_rule(y_euler_2, y_euler_2_runge, 1))
        print('\nМетод Рунге-Кутты, с шагом h = 0.1')
        print(runge_rule(y_r_k_1, y_r_k_1_runge, 4))
        print('\nМетод Рунге-Кутты, с шагом h = 0.05')
        print(runge_rule(y_r_k_2, y_r_k_2_runge, 4))