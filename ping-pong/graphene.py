import sys
import random
import numpy as np
from math import ceil, exp
from pybloom_live import BloomFilter
sys.path.append('.')
from search_params import search_params

PYBLT_PATH = '/home/pinar/IBLT-optimization/'
sys.path.append(PYBLT_PATH)
from pyblt import PYBLT
PYBLT.set_parameter_filename(PYBLT_PATH+'param.export.0.995833.2018-07-17.csv')

# Constants
PATH = '.'
TXN_SHORT_BYTES_CB = 6
TXN_SHORT_BYTES = 8
BITSIZE = TXN_SHORT_BYTES * 8 # number of bytes per transaction, converted to bits per transaction
TAU = 12 # bytes per IBLT cell

BLK_SIZE = [200, 2000, 10000] # num of txns in block
B = [1/240]
FRACTION = np.arange(1, 0, -0.05)
NUM_TRIAL = 1000

def create_mempools(mempool_size, fraction, blk_size, num_doesnt_have):
    # create a random set of transactions in a mempool
    blk = [random.getrandbits(BITSIZE) for x in range(blk_size)]
    num_has = int(blk_size * fraction) # fraction of txns receiver has
    #print('num has', num_has)
    in_blk = random.sample(blk, num_has)
    # num_doesnt_have = mempool_size - num_has # rest of txns in mempool
    in_mempool = [random.getrandbits(BITSIZE) for x in range(num_doesnt_have)]
    receiver_mempool = in_blk + in_mempool
    return blk, receiver_mempool

# Check whether bk was reconstructed properly
def decode_blk(result, passed, blk):
    not_in_blk = set()
    in_blk = set()
    for key in result:
        if result[key][1] == 1:
            not_in_blk.add(key)
        elif result[key][1] == -1:
            in_blk.add(key)
    possibly_in_blk = set(passed)
    possibly_in_blk.difference_update(not_in_blk)
    reconstructed_blk = list(in_blk.union(possibly_in_blk))
    flag = set(reconstructed_blk) == set(blk)
    return flag, in_blk

def try_ping_pong(first_IBLT, second_IBLT, in_blk, not_in_blk):
    flag = False
    first_boolean, result = second_IBLT.list_entries()
    if len(result) == 0:
        second_boolean, result = first_IBLT.list_entries()
        flag = True
    if len(result) == 0: # both results are zero
        if first_boolean == True and second_boolean == True:
            return True, in_blk, not_in_blk
        else:
            return False, in_blk, not_in_blk
    for key in result:
        if result[key][1] == 1:
            not_in_blk.add(key)
            second_IBLT.erase(key, 0x0)
            first_IBLT.erase(key, 0x0)
        elif result[key][1] == -1:
            in_blk.add(key)
            second_IBLT.insert(key, 0x0)
            first_IBLT.insert(key, 0x0)
    if flag == True:
        return try_ping_pong(first_IBLT, second_IBLT, in_blk, not_in_blk)
    else:
        return try_ping_pong(second_IBLT, first_IBLT, in_blk, not_in_blk)

def trial(fd):
    params = search_params()
    for blk_size in BLK_SIZE:
        true_false_positives = blk_size
        for bound in B:
            for fraction in FRACTION:

                # True_positives is the number of txns in the blk the receiver has
                true_positives = int(blk_size * fraction)
                mempool_size = true_false_positives + true_positives

                print('Running %d trials for parameter combination: extra txns in mempool %d blk size %d CB bound %f fraction %f' % (NUM_TRIAL, true_false_positives, blk_size, bound, fraction))

                # Size of Compact block (inv + getdata)
                getdata = (1-fraction) * blk_size * TXN_SHORT_BYTES_CB
                inv = blk_size * TXN_SHORT_BYTES_CB
                compact = inv + getdata

                for i in range(NUM_TRIAL):
                    blk, receiver_mempool = create_mempools(mempool_size, fraction, blk_size, true_false_positives)

                    # Sender creates BF of blk
                    a, fpr_sender, iblt_rows_first = params.CB_solve_a(mempool_size, blk_size, blk_size, 0, bound)
                    bloom_sender = BloomFilter(blk_size, fpr_sender)
                    tmp = blk_size + 0.5
                    exponent = (-bloom_sender.num_slices*tmp) / (bloom_sender.num_bits-1)
                    real_fpr_sender = (1-exp(exponent)) ** bloom_sender.num_slices
                    #exponent = (-bloom_sender.num_slices*blk_size) / bloom_sender.num_bits
                    #tmp = (1-exp(exponent)) ** bloom_sender.num_slices
                    #real_fpr_sender = max(tmp, fpr_sender)
                    #assert real_fpr_sender >= fpr_sender

                    # Sender creates IBLT of blk
                    iblt_sender_first = PYBLT(a, TXN_SHORT_BYTES)

                    # Add to BF and IBLT
                    for txn in blk:
                        bloom_sender.add(txn)
                        iblt_sender_first.insert(txn, 0x0)

                    # Receiver computes how many items pass through BF of sender and creates IBLT
                    iblt_receiver_first = PYBLT(a, TXN_SHORT_BYTES)
                    Z = []
                    for txn in receiver_mempool:
                        if txn in bloom_sender:
                            Z.append(txn)
                            iblt_receiver_first.insert(txn, 0x0) #(id and content)
                    z = len(Z)
                    observed_false_positives = z - true_positives

                    # Eppstein subtraction
                    T = iblt_receiver_first.subtract(iblt_sender_first)
                    boolean, result = T.list_entries()
                    #assert boolean == False

                    # Check whether decoding successful
                    if boolean == True:
                        flag, in_blk = decode_blk(result, Z, blk)

                        # Each component of graphene blk size
                        first_IBLT = (iblt_rows_first * TAU)
                        first_BF = (bloom_sender.num_bits / 8.0)
                        extra = (len(in_blk) * TXN_SHORT_BYTES)
                        # Compute size of Graphene block
                        graphene = first_IBLT + first_BF + extra

                        fd.write(str(true_false_positives)+'\t'+str(blk_size)+'\t'+str(bound)+'\t'+str(fraction)+'\t'+str(mempool_size)+'\t'+str(fpr_sender)+'\t'+str(real_fpr_sender)+'\t'+str(0)+'\t'+str(a)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(z)+'\t'+str(0)+'\t'+str(observed_false_positives)+'\t'+str(boolean and flag)+'\t'+str(False)+'\t'+str(graphene)+'\t'+str(first_IBLT)+'\t'+str(first_BF)+'\t'+str(0)+'\t'+str(0)+'\t'+str(extra)+'\t'+str(iblt_rows_first)+'\t'+str(0)+'\t'+str(compact)+'\n')
                    else:
                        # Receiver creates BF of txns that passed through sender's BF
                        # print('z', z)
                        # print('bound', bound)
                        x_star = params.search_x_star(z, mempool_size, real_fpr_sender, bound, blk_size)
                        temp = (mempool_size - x_star) * real_fpr_sender
                        y_star = params.CB_bound(temp, real_fpr_sender, bound)
                        #print('y_star', y_star)
                        y_star = ceil(y_star)
                        b, fpr_receiver, iblt_rows_second = params.solve_a(blk_size, z, x_star, y_star)

                        bloom_receiver = BloomFilter(z, fpr_receiver)
                        for txn in Z:
                            bloom_receiver.add(txn)

                        # Receiver determines IBLT size
                        iblt_sender_second = PYBLT(b+y_star, TXN_SHORT_BYTES)
                        # Sender creates IBLT of blk again and sends txns that do not pass through BF of receiver
                        count = 0
                        for txn in blk:
                            iblt_sender_second.insert(txn, 0x0)
                            if txn not in bloom_receiver:
                                T.insert(txn, 0x0) # add txns just received to subtracted IBLT
                                Z = Z + [txn] # sends the txn to the receiver
                                count = count+1

                        iblt_receiver_second = PYBLT(b+y_star, TXN_SHORT_BYTES)
                        for txn in Z:
                            iblt_receiver_second.insert(txn, 0x0)

                        # Eppstein subtraction
                        T_second = iblt_receiver_second.subtract(iblt_sender_second)
                        boolean, result = T_second.list_entries()
                        #print(boolean)
                        #print('Z', z)

                        # Check whether blk was reconstructed properly
                        flag, in_blk = decode_blk(result, Z, blk)

                        final = False
                        if boolean == False or flag == False:
                            final, in_blk, not_in_blk = try_ping_pong(T, T_second, set(), set())
                            #print('Ping pong result', final)
                            if final == True:
                                possibly_in_blk = set(Z)
                                possibly_in_blk.difference_update(not_in_blk)
                                reconstructed_blk = list(in_blk.union(possibly_in_blk))
                                assert set(reconstructed_blk) == set(blk)

                        # Each component of graphene blk size
                        first_IBLT = (iblt_rows_first * TAU)
                        first_BF = (bloom_sender.num_bits / 8.0)
                        second_IBLT = (iblt_rows_second * TAU)
                        second_BF = (bloom_receiver.num_bits / 8.0)
                        extra = (len(in_blk) * TXN_SHORT_BYTES)
                        # Compute size of Graphene block
                        graphene = first_IBLT + first_BF + second_IBLT + second_BF + extra

                        fd.write(str(true_false_positives)+'\t'+str(blk_size)+'\t'+str(bound)+'\t'+str(fraction)+'\t'+str(mempool_size)+'\t'+str(fpr_sender)+'\t'+str(real_fpr_sender)+'\t'+str(fpr_receiver)+'\t'+str(a)+'\t'+str(b)+'\t'+str(x_star)+'\t'+str(y_star)+'\t'+str(z)+'\t'+str(count)+'\t'+str(observed_false_positives)+'\t'+str(boolean and flag)+'\t'+str(final)+'\t'+str(graphene)+'\t'+str(first_IBLT)+'\t'+str(first_BF)+'\t'+str(second_IBLT)+'\t'+str(second_BF)+'\t'+str(extra)+'\t'+str(iblt_rows_first)+'\t'+str(iblt_rows_second)+'\t'+str(compact)+'\n')

                    fd.flush()

def main():
    fd =open(PATH+'/results/ping-pong-%d.csv' % (random.getrandbits(25)),'w')
    trial(fd)
    fd.close()

if __name__=="__main__":
    main()
