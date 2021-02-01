import numpy as np
from scipy import special
from math import log, ceil, sqrt, exp

# Constants
PYBLT_PATH = '/home/pinar/IBLT-optimization/'
NUM = 5000
START = 0.001

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
        #if np.isnan(error_rate):
        #    return 10e19
        assert error_rate > 0 and error_rate < 1
        num_slices = int(ceil(log(1.0/error_rate, 2)))
        #print('num', num_slices)
        bits_per_slice = int(ceil((capacity*abs(log(error_rate)))/(num_slices*(log(2)**2))))
        return num_slices * bits_per_slice / 8.0

    def total(self, a, fpr, n, y):
        """ Total cost of graphene for given parameters"""
        #print('a', a)
        #print('size', self.params[a-1][3])
        assert a > 0, "a==%f <= 0" % a
        s = self.bf_num_bytes(fpr, n)
        if ceil(a)+y-1 >= len(self.params): # difference too much
            tmp = (ceil(a)+y) * 1.362549
            rows = ceil(tmp)
            i = rows *  12
        else:
            rows = self.params[ceil(a)+y-1][3]
            i = rows * 12
        total = s + i
        return total, rows

    def solve_a(self, m, n, x, y):
        assert x <= m
        denom = float(m) - x
        if denom == 0:
            denom = float(1)
        assert denom > 0
        min_a = ceil(START)
        min_fpr = START/denom
        min_total, min_iblt_rows = self.total(START, START/denom, n, y)
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
            if t < min_total:
                min_total = t
                min_a = ceil(c)
                min_fpr = c/denom
                min_iblt_rows = iblt_rows
        assert min_fpr <= 1 and min_fpr > 0, "a out of bounds"
        return min_a, min_fpr, min_iblt_rows

    def CB_bound(self, a, fpr, bound):
        #a = (m-x) * fpr
        s = (-1 * log(bound)) / a
        temp =  sqrt((s*(s+8)))
        delta_1 = 0.5 * (s + temp)
        delta_2 = 0.5 * (s - temp)
        assert delta_1 >= 0
        #print('delta_2', delta_2)
        assert delta_2 <= 0
        return (1+delta_1) * a

    def CB_solve_a(self, m, n, x, y, bound):
        #assert x+y == n and x <= m
        assert x <= m
        denom = float(m) - x
        if denom == 0:
            denom = float(1)
        assert denom > 0
        min_a = self.CB_bound(START, START/denom, bound)
        min_fpr = START/denom
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
            if t < min_total:
                min_total = t
                min_a = ceil(b)
                min_fpr = c/denom
                min_iblt_rows = iblt_rows
        assert min_fpr <= 1 and min_fpr > 0, "a out of bounds"
        return min_a, min_fpr, min_iblt_rows

    def compute_delta(self, z, x, m, fpr):
        temp = (m-x) * fpr
        temp = (z-x) / temp
        return temp - 1

    def compute_RHS(self, delta, m, x, fpr):
        #num = exp(delta)
        #print('delta', delta)
        #denom = (1+delta)**(1+delta)
        #print('denom', denom)
        #base = num / denom
        exponent = (m-x) * fpr
        base = exp(delta-(1+delta) * log(1+delta))
        #cb=base**(negatives*f_s)
        return base ** exponent

    def compute_binomial_prob(self, m, x, z, fpr):
        trials = m-x
        success = z-x
        term_one = fpr ** success
        term_two = (1-fpr) ** (trials-success)
        term_three = special.binom(trials, success)  # binom coefficient
        return (term_one * term_two * term_three)

    def search_x_star(self, z, mempool_size, fpr, bound, blk_size):
        total = 0
        x_star = 0
        val = min(z, blk_size)
        for x in range(0, val+1):
            delta = self.compute_delta(z, x, mempool_size, fpr)
            result = self.compute_RHS(delta, mempool_size, x, fpr)
            #result = self.compute_binomial_prob(mempool_size, x, z, fpr)
            if (total + result > bound):
                if x == 0:
                    x_star = x
                else:
                    x_star = x - 1
                break
            else:
                total += result
        # print('x_star', x_star)
        return x_star

#with open('/Users/Pinar/Desktop/empirical-dist/ping-pong.csv', 'r') as f:
#    with open('/Users/Pinar/Desktop/empirical-dist/real-data.csv', 'w+') as f:
#        ind = 1
#        for line in f.readlines():
#            lst = line.split('\t')
#            if len(lst) != 24:
#                print('len', len(lst))
#                print('lst', lst)
#                print('line_num', ind)
#                break
#            ind += 1
#
# with open('/Users/Pinar/Desktop/empirical-dist/ping-pong-10000.csv', 'r') as f:
#    with open('/Users/Pinar/Desktop/empirical-dist/short-half.csv', 'a') as f2:
#        with open('/Users/Pinar/Desktop/empirical-dist/long-half.csv', 'a') as f3:
#            for line in f.readlines():
#                lst = line.split('\t')
#                if len(lst) == 24:
#                    f2.write(line)
#                else:
#                    f3.write(line)
#
# with open('/Users/Pinar/Desktop/empirical-dist/ping-pong-temp.csv', 'r') as f:
#     with open('/Users/Pinar/Desktop/empirical-dist/short-half.csv', 'w+') as f2:
#         with open('/Users/Pinar/Desktop/empirical-dist/long-half.csv', 'w+') as f3:
#             cnt_thousand = 0
#             cnt_hundred = 0
#             total = 0
#             for line in f.readlines():
#                 lst = line.split('\t')
#                 assert lst[0] == '200' or lst[0] == '2000' or lst[0] == '10000'
#                 if lst[0] != '10000':
#                     total += 1
#                     if lst[0] == '2000':
#                         cnt_thousand += 1
#                     else:
#                         cnt_hundred += 1
#                     if len(lst) == 24:
#                         f2.write(line)
#                     else:
#                         f3.write(line)
#             print('200 cnt:', cnt_hundred)
#             print('2000 cnt:', cnt_thousand)
#             print('total:', total)
