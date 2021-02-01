import sys
import random
import numpy as np
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
B = [1-(239/240)]
MULTIPLE_BLK = np.arange(0, 5.5, 0.5)
NUM_TRIAL = 1000

def create_mempools(mempool_size, blk_size):
    # create a random set of transactions in a mempool
    blk = [random.getrandbits(BITSIZE) for x in range(blk_size)]
    in_blk = blk.copy()
    num_doesnt_have = mempool_size - blk_size # rest of txns in mempool
    in_mempool = [random.getrandbits(BITSIZE) for x in range(num_doesnt_have)]
    receiver_mempool = in_blk + in_mempool
    return blk, receiver_mempool

# Check whether blk was reconstructed properly
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

def trial(fd):
    params = search_params()
    for multiple in MULTIPLE_BLK:
        for blk_size in BLK_SIZE:
            for bound in B:
                true_false_positives = int(multiple * blk_size)
                mempool_size = true_false_positives + blk_size
                print('Running %d trials for parameter combination: multiple of blk %d blk size %d' % (NUM_TRIAL, multiple, blk_size))
                for x in range(NUM_TRIAL):
                    blk, receiver_mempool = create_mempools(mempool_size, blk_size)

                    # Sender creates BF of blk
                    a, fpr_sender, iblt_rows_first = params.solve_a(mempool_size, blk_size, blk_size, 0)
                    bloom_sender = BloomFilter(blk_size, fpr_sender)

                    # Sender creates IBLT of blk
                    iblt_sender_first = PYBLT(a, TXN_SHORT_BYTES)

                    # Add to BF and IBLT
                    for txn in blk:
                        bloom_sender.add(txn)
                        iblt_sender_first.insert(txn, 0x0)

                    # Receiver computes how many items pass through BF of sender and creates IBLT
                    iblt_receiver_first = PYBLT(a, TXN_SHORT_BYTES)
                    passed = []
                    for txn in receiver_mempool:
                        if txn in bloom_sender:
                            passed.append(txn)
                            iblt_receiver_first.insert(txn, 0x0) #(id and content)
                    observed_false_positives = len(passed) - blk_size

                    # Eppstein subtraction
                    T = iblt_receiver_first.subtract(iblt_sender_first)
                    boolean, result = T.list_entries()

                    # Check whether decoding successful
                    flag, in_blk = decode_blk(result, passed, blk)
                    # Each component of graphene blk size
                    first_IBLT = (iblt_rows_first * TAU)
                    first_BF = (bloom_sender.num_bits / 8.0)
                    # Compute size of Graphene block
                    graphene = first_IBLT + first_BF

                    # Size of Compact block (inv + getdata)
                    # getdata = (1-fraction) * BLK_SIZE * TXN_SHORT_BYTES
                    getdata = 0
                    inv = blk_size * TXN_SHORT_BYTES_CB
                    compact = inv + getdata
                    #print('getdata', getdata)
                    #print('Pinar', (len(in_blk) * TXN_SHORT_BYTES))
                    #print((boolean and flag))
                    #assert getdata == (len(in_blk) * TXN_SHORT_BYTES)

                    fd.write(str(mempool_size)+'\t'+str(blk_size)+'\t'+str(bound)+'\t'+str(fpr_sender)+'\t'+str(a)+'\t'+str(iblt_rows_first)+'\t'+str(true_false_positives)+'\t'+str(observed_false_positives)+'\t'+str(compact)+'\t'+str(boolean and flag)+'\t'+str(graphene)+'\t'+str(first_IBLT)+'\t'+str(first_BF)+'\t'+str(multiple)+'\n')
                    fd.flush()

def main():
    fd =open(PATH+'/results/experiment1-noCB-%d.csv' % (random.getrandbits(25)),'w')
    trial(fd)
    fd.close()

if __name__=="__main__":
    main()
