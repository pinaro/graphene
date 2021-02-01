import numpy as np
from math import log, ceil, sqrt

# Constants
PYBLT_PATH = '/home/pinar/IBLT-optimization/'
NUM = 5000
START = 0.001
MAX = 10e18

class search_params:
    # the main method here is solve_a.
    # that method calls these first two methods

    # This computes the number of bytes in a bloom Filter
    # this code snippet is ripped or based on the pybloom class,
    # from the point where it allocates the bit array given a specified FPR
    # if we changed to another BF, we'd have to change this function.
    # this method is where using ceiling functions happens.
    def __init__(self):
        self.params = []
        with open(PYBLT_PATH+'param.export.0.995833.2018-07-17.csv') as f:
            # Skip first line
            line = f.readline()
            line = f.readline()
            while line:
                items, hedge, numhashes, size, p = line.split(',')
                self.params.append((int(items), float(hedge), int(numhashes), int(size), float(p)))
                line = f.readline()

    def bf_num_bytes(self, error_rate, capacity):
        assert error_rate != 1
        if np.isnan(error_rate):
            return 10e19
        assert error_rate > 0 and error_rate <= 1
        num_slices = int(ceil(log(1.0/error_rate, 2)))
        #print('num', num_slices)
        bits_per_slice = int(ceil((capacity*abs(log(error_rate)))/(num_slices*(log(2)**2))))
        return num_slices * bits_per_slice / 8.0

    def total(self, a, fpr, n, y):
        """ Total cost of graphene for given parameters"""
        if ceil(a)+y-1 >= len(self.params): # can't compute any longer
            return MAX, MAX
        assert a > 0, "a==%f <= 0" % a
        s = self.bf_num_bytes(fpr, n)
        #print('a', a)
        #print('size', self.params[a-1][3])
        i = self.params[ceil(a)+y-1][3] * 12
        total = s + i
        return total, ceil(self.params[ceil(a)+y-1][3])

    def solve_a(self, m, n, x, y=0):
        assert x+y == n and x <= m
        denom = float(m) - x
        if denom == 0:
            denom = float(1)
        assert denom > 0
        min_a = START
        min_total, min_iblt_rows = self.total(min_a, min_a/denom, n, y)
        search_space_for_a = np.linspace(start=START, stop=denom, endpoint=False, num=NUM)
        #print('Params for search:')
        #print('m', m)
        #print('n', n)
        #print('x', x)
        #print('y', y)
        #print(search_space_for_a)
        for c in search_space_for_a:
            #print('c', c)
            t, iblt_rows = self.total(c, c/denom, n, y)
            if t == MAX:
                break
            if t < min_total:
                min_total = t
                min_a = ceil(c)
                min_fpr = c/denom
                min_iblt_rows = iblt_rows
        assert min_fpr <= 1 and min_fpr > 0, "a out of bounds"
        return min_a, min_fpr, min_iblt_rows

    def CB_bound(self, a, fpr, B=1e-3):
        s = (-1 * log(B)) / a
        temp =  sqrt((s*(s+8)))
        delta_1 = 0.5 * (s + temp)
        delta_2 = 0.5 * (s - temp)
        assert delta_1 >= 0
        #print('delta_2', delta_2)
        assert delta_2 <= 0
        return (1+delta_1) * a

    def CB_solve_a(self, m, n, x, y, bound):
        assert x+y == n and x <= m
        denom = float(m) - x
        if denom == 0:
            denom = float(1)
        assert denom > 0
        min_a = self.CB_bound(START, START/denom, bound)
        min_total, min_iblt_rows = self.total(min_a, START/denom, n, 0)
        search_space_for_a = np.linspace(start=START, stop=denom, endpoint=False, num=NUM)
        #print('Params for search:')
        #print('m', m)
        #print('n', n)
        #print('x', x)
        #print('y', y)
        #print(search_space_for_a)
        for c in search_space_for_a:
            #print('c', c)
            b = self.CB_bound(c, c/denom, bound)
            t, iblt_rows = self.total(b, c/denom, n, 0)
            if t == MAX:
                break
            if t < min_total:
                min_total = t
                min_a = ceil(b)
                min_fpr = c/denom
                min_iblt_rows = iblt_rows
        assert min_fpr <= 1 and min_fpr > 0, "a out of bounds"
        return min_a, min_fpr, min_iblt_rows
