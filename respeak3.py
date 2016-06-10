import numpy as np
import matplotlib.pyplot as mpl
import os
from scipy.optimize import curve_fit
import math
import sys
import csv


def func(f, R0, f0, dW, m, b):
    # noinspection PyUnresolvedReferences
    return (R0 + np.power((2 * (f - f0)) / dW, 2)) / (1 + np.power((2 * (f - f0)) / dW, 2)) + m * f + b


def back_line(f, m, b):
    return m * f + b


def main(argv):

    try:
        res_peak_folder = sys.argv[1]
        print('res_peak_folder = ', res_peak_folder, '\n')
    except IndexError:
        print('no folder command line argument\n')
        res_peak_folder = input('path containing TRMC data: ')

    file_name = os.path.join(res_peak_folder, 'respeak.txt')

    res_peak = np.genfromtxt(file_name, delimiter=',', skip_header=12)

    # now finds min frequency automatically
    guess = [.1, 0, 1e7, 1e-8, 0]

    x_data = res_peak[:, 0]
    # noinspection PyUnresolvedReferences
    y_data = np.power(res_peak[:, 1], 2)/ 50
    # noinspection PyUnresolvedReferences
    y_data = np.divide(y_data[:], np.amax(y_data[:]))  # normalize resonance peak so maximum at 1

    guess[1] = x_data[np.argmin(y_data)]

    background_x = np.array(x_data[0:50])  # + res_peak[-1:-9, 0]])
    background_y = np.array(y_data[0:50])  # + res_peak[-1:-9, 1]])
    background_x = np.append(background_x, x_data[-50:])
    background_y = np.append(background_y, y_data[-50:])

    line_guess = [-2e9, 0]

    try:
        p_opt, p_cov = curve_fit(back_line, np.array(background_x), np.array(background_y), p0=line_guess)
        print(p_opt)
        # print(p_cov)
        linear_fit = back_line(x_data, p_opt[0], p_opt[1])
        np.savetxt(os.path.join(res_peak_folder, 'linear_background.txt'), linear_fit)
        mpl.plot(x_data, linear_fit)

    except RuntimeError:
        print("Error - curve_fit failed")
        return 1

    mpl.plot(x_data, y_data)

    # mpl.show()
    try:
        p_opt, p_cov = curve_fit(func, x_data.T, y_data.T, p0=guess)
        print(p_opt)
        # print(p_cov)
        fit = func(x_data, p_opt[0], p_opt[1], p_opt[2], p_opt[3], p_opt[4])
        np.savetxt(os.path.join(res_peak_folder, 'res_peakfit.txt'), fit)
        mpl.plot(x_data, fit)

        print(linear_fit[np.argmin(y_data)])
        print('uncorrected R0 = ' + str(np.amin(y_data)))
        R0 = np.amin(y_data) / linear_fit[np.argmin(y_data)]
        Q = p_opt[1] / p_opt[2]
        responseTime = (Q / math.pi / p_opt[1])
        if Q < 0:
            print('\n\nWARNING: Q<0. DO NOT USE THESE PARAMETERS FOR FITTING. '
                  'TRY AGAIN. resonance_params.txt not written.\n\n')
        else:
            with open(os.path.join(res_peak_folder, 'resonance_params.csv'), 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
                writer.writerow(['Q'] + ['R0'] + ['Response Time'] + ['f0'])
                writer.writerow([Q] + [R0] + [responseTime] + [p_opt[1]])
                print('Q = ' + str(Q) + '\nR0 = ' + str(R0) + '\nresponseTime = ' + str(responseTime) + '\nf0 = ' + str(p_opt[1]))
                # np.savetxt(res_peak_folder + '\\res_peakfit.txt',fit, delimiter=',')
                # np.savetxt(res_peak_folder + '\\y_data.txt',y_data, delimiter=',')

    except RuntimeError:
        print("Error - curve_fit failed")
        return 1

    mpl.savefig(os.path.join(res_peak_folder, 'res_peak.png'))
    mpl.show()

if __name__ == '__main__':
    main(sys.argv[1:])


