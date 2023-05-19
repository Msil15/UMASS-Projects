import numpy as np

x_i = [4.0, 4.2, 4.5, 4.7, 5.1, 5.5, 5.9, 6.3, 6.8, 7.1]
y_i = [102.56, 113.18, 130.11, 142.05, 167.53, 195.14, 224.87, 256.73, 299.50, 326.72]

x_2 = [-2, 0, 4, 2]
y_2 = [4, 2, 10, 8]

def least_squares(n, x_i, y_i):
    m = len(x_i)
    A_T = [[0 for m in range(m)]for n in range(n+1)]

    for j in range(m):
        for i in range(n + 1):
            A_T[i][j] = x_i[j]**(i)


    A_T = np.asmatrix(A_T)
    A = np.matrix.transpose(A_T)
    y_i = np.asmatrix(y_i).T

    left = np.matmul(A_T, A)
    left_inv = np.linalg.inv(left)

    right = np.matmul(A_T, y_i)

    c = left_inv*right

    Error_List = []

    for x in x_i:
        term_list = []

        for j in range(n + 1):
            term = c[j]*(x**j)
            term_list.append(term)

        error_term = (sum(term_list) - y_i[x_i.index(x)])**2

        Error_List.append(error_term)

    Error_complete = sum(Error_List)


    print(c)
    print(Error_complete)

    return

least_squares(1, x_i, y_i)
print('')
least_squares(2, x_i, y_i)
print('')
least_squares(3, x_i, y_i)
print('')

def least_squares_exp(n, x_i, y_i):
    m = len(x_i)
    A_T = [[0 for m in range(m)]for n in range(n+1)]

    for j in range(m):
        for i in range(n + 1):
            A_T[i][j] = x_i[j]**(i)


    A_T = np.asmatrix(A_T)
    A = np.matrix.transpose(A_T)
    y_i = np.log(y_i)
    y_i = np.asmatrix(y_i).T

    left = np.matmul(A_T, A)
    left_inv = np.linalg.inv(left)

    right = np.matmul(A_T, y_i)

    c = left_inv*right

    Error_List = []

    for x in x_i:
        term_list = []

        for j in range(n + 1):
            term = c[j]*(x**j)
            term_list.append(term)

        error_term = (sum(term_list) - y_i[x_i.index(x)])**2

        Error_List.append(error_term)

    Error_complete = sum(Error_List)


    print(c)
    print(Error_complete)

    return

least_squares_exp(1, x_i, y_i)
print('')

def least_squares_power(n, x_i, y_i):
    m = len(x_i)

    x_i = np.log(x_i).tolist()

    A_T = [[0 for m in range(m)]for n in range(n+1)]

    for j in range(m):
        for i in range(n + 1):
            A_T[i][j] = x_i[j]**(i)


    A_T = np.asmatrix(A_T)
    A = np.matrix.transpose(A_T)
    y_i = np.log(y_i)
    y_i = np.asmatrix(y_i).T

    left = np.matmul(A_T, A)
    left_inv = np.linalg.inv(left)

    right = np.matmul(A_T, y_i)

    c = left_inv*right

    Error_List = []

    for x in x_i:
        term_list = []

        for j in range(n + 1):
            term = c[j]*(x**j)
            term_list.append(term)

        error_term = (sum(term_list) - y_i[x_i.index(x)])**2

        Error_List.append(error_term)

    Error_complete = sum(Error_List)


    print(c)
    print(Error_complete)

    return

least_squares_power(1, x_i, y_i)