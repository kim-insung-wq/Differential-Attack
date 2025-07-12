from ctypes import *
import os
from sys import platform
from sys import version_info
import struct
import math
import multiprocessing as mt
import numpy as np
import pickle
import mpmath

__SYSTEM_BIT    = struct.calcsize("P") * 8
__THIS_FILE_DIR = os.path.dirname(os.path.realpath(__file__))
__PYTHON_VER    = version_info.major

if __SYSTEM_BIT == 64:
    __LIB_PATH = os.path.join(__THIS_FILE_DIR,"../","x64","Release")
else:#default
    __LIB_PATH = os.path.join(__THIS_FILE_DIR,"../","x64","Release")

if platform == "linux" or platform == "linux2":
    # linux
    __LIB_PATH = os.path.join(__LIB_PATH,"libFLY.so")
    libastbb_srbn = CDLL(__LIB_PATH)
elif platform == "win32":
    # Windows
    __LIB_PATH = os.path.join(__LIB_PATH,"FLY.dll")
    libastbb_srbn = CDLL(__LIB_PATH)


class c_uint32_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_uint32 * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_int32_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_int32 * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_int_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_int * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_uint8_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_uint8 * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_uint64_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_uint64 * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_int16_array:
    def __init__(self, length=1, values = None):
        self.length = length
        self.array  = (c_int16 * length)()
        if values != None:
            list_values = list(values)
            for idx in range(len(list_values)):
                self.array[idx] = (list_values[idx])

class c_ubyte_array:
    def __init__(self, length=1):
        self.length = length
        self.array  = (c_ubyte * length)()

class c_double_array:
    def __init__(self, length=1):
        self.length = length
        self.array  = (c_double * length)()


NEGA     = -1
ZERO     = 0
POSI     = 1
class CORR_t(Structure):
    _fields_ = [
        ("sign",        c_int),
        ("magnitude",   c_double),
        ]



NUM_SBOX_IN_A_STATE_THRESHOLD = 8
SBOX_BIT_SIZE_THRESHOLD       = 8
SBOX_CARDINALITY_THRESHOLD    = 2**(SBOX_BIT_SIZE_THRESHOLD)


class SRBPN_INFO_t(Structure):
    _fields_ = [
        ("ALG_NAME",            c_char_p),
        ("SBOX_BIT_SIZE",       c_int),
        ("SBOX_CARDINALITY",    c_int),
        ("SBOX",                POINTER(c_uint8)),

        ("NUM_SBOX_IN_A_STATE", c_int),
        ("OFFSET",              POINTER(c_int)),
        ]


class BitPattern64(Structure):
    _fields_ = [
        ("must_be_1_mask",    c_uint8 * SBOX_BIT_SIZE_THRESHOLD),
        ("must_be_0_mask",    c_uint8 * SBOX_BIT_SIZE_THRESHOLD),
        ("dont_care_mask",    c_uint8 * SBOX_BIT_SIZE_THRESHOLD),
        ]



Prep_Dif_Trail_Searching          = libastbb_srbn.Prep_Dif_Trail_Searching
Best_Trail_Prob_Only              = libastbb_srbn.Best_Trail_Prob_Only
Best_Trail_Prob_All               = libastbb_srbn.Best_Trail_Prob_All
Best_Trail_Prob_IO                = libastbb_srbn.Best_Trail_Prob_IO
Computation_LS_of_DDT             = libastbb_srbn.Computation_LS_of_DDT
Computation_LS_of_INV_DDT         = libastbb_srbn.Computation_LS_of_INV_DDT



def int_to_bytes(int_data, data_len):
    
    return int_data.to_bytes(data_len, byteorder='big', signed = False)

def bytes_to_int(bytes_data):
        
    return int.from_bytes(bytes_data, byteorder = 'big', signed=False)

def int_to_cbytes(int_data, data_len):
    cbytes = (c_ubyte * data_len)()
    tmp = int_data
    for idx in range(0, data_len):
        cbytes[data_len - idx - 1] = (tmp & 0xff)
        tmp = tmp >> 8
    return cbytes

def cbytes_to_int(cbytes_data):
    tmp = bytes(cbytes_data)
    return bytes_to_int(tmp)

def cbytes_buf(data_len):   
    return (c_ubyte*data_len)()


def int_to_cshorts(int_data, data_len):
    cshorts = (c_int16 * data_len)()
    tmp = int_data
    for idx in range(0, data_len):
        cshorts[idx] = (tmp & 0xffff)
        tmp = tmp >> 16
    return cshorts

def cshorts_buf(data_len):   
    return (c_int16 * data_len)()



def int_to_clonglongs(int_data, data_len):
    cshorts = (c_int64 * data_len)()
    tmp = int_data
    for idx in range(0, data_len):
        cshorts[idx] = (tmp & 0xffffffffffffffff)
        tmp = tmp >> 64
    return cshorts

def clonglongs_buf(data_len):   
    return (c_int64 * data_len)()




def int_to_binstr(int_data, bitsize):
    return format(int_data, 'b').zfill(bitsize)

def int_to_hexstr(int_data, bitsize):
    return format(int_data, 'X').zfill(int(math.ceil(bitsize/4)))


def print_out_opt(ct_arry, ct_len, output_opt = 'int'):
    ct_int   = cbytes_to_int(ct_arry)
    ct_bytes = bytes(ct_arry)
    opt = output_opt.lower()
    if opt == 'int':
        return ct_int
    elif opt== 'bytes':
        return ct_bytes
    elif opt == 'cbytes':
        return ct_arry
    elif opt == 'hexlist':
        return list(ct_bytes)
    elif opt == 'binstr':
        return int_to_binstr(ct_int, ct_len*8)
    elif opt == 'hexstr':
        return int_to_hexstr(ct_int, ct_len*8)
    else:
        print(
        '''\
# output_opt must be one of the following options
# 'int'    : 136792598789324718765670228683992083246
# 'bytes'  : b'f\\xe9K\\xd4\\xef\\x8a,;\\x88L\\xfaY\\xca4+.'
# 'cbytes' : <CipherTools.c_ubyte_Array_16 object at 0x000001994C70B3C8>
# 'hexlist': [102, 233, 75, 212, 239, 138, 44, 59, 136, 76, 250, 89, 202, 52, 43, 46]
# 'binstr' : '01100110111010010100101111010100111011111000101000101100001110111000100001001100111110100101100111001010001101000010101100101110'
# 'hexstr' : '66e94bd4ef8a2c3b884cfa59ca342b2e' '''
        )
        return False


class SRBPN_cipher:
    def __init__(self, sbox, offset, num_sbox_in_a_state, algname):
        #size check
        self.SBOX                   = list(sbox)
        self.OFFSET                 = list(offset)
        self.SBOX_CARDINALITY       = len(self.SBOX)
        self.NUM_OFFSET             = len(self.OFFSET)
        self.NUM_SBOX_IN_A_STATE    = num_sbox_in_a_state
        
        tmp_bin_str = bin(self.SBOX_CARDINALITY)
        if(tmp_bin_str.count("1") != 1):
            raise TypeError("S-box cardinality must be 2^{*}")
        self.SBOX_BIT_SIZE = int(math.log2(self.SBOX_CARDINALITY))
        if(self.SBOX_BIT_SIZE != self.NUM_OFFSET):
            raise TypeError("S-box bit size must be eqaul to the number of offsets")
        if(max(self.OFFSET) >= self.NUM_SBOX_IN_A_STATE):
            raise TypeError("Invalid offset!-offsets must be smaller than the number of sboxes in a state")


        self.ct_algname = create_string_buffer(bytes(algname, "ascii"), 256)

        self.ct_sbox    = c_uint8_array(self.SBOX_CARDINALITY, self.SBOX).array

        self.ct_offset  = c_int_array(self.NUM_OFFSET, self.OFFSET).array

        CIP_INFO = SRBPN_INFO_t(
            cast(self.ct_algname, c_char_p),
            c_int(self.SBOX_BIT_SIZE),
            c_int(self.SBOX_CARDINALITY),
            cast(self.ct_sbox, POINTER(c_uint8)),
            c_int(self.NUM_SBOX_IN_A_STATE),
            cast(self.ct_offset, POINTER(c_int))
        )

        Prep_Dif_Trail_Searching(pointer(CIP_INFO))
        Computation_LS_of_DDT()
        Computation_LS_of_INV_DDT()

    def get_BDP(self, target_round, prev_rsts):

        if len(prev_rsts) < target_round:
            raise TypeError("Previous Best Differential Probabilities are Required")
        
        ct_prv_prob_rst = c_double_array(target_round).array
        ct_prob_rst     = c_double(0)

        ##prev_corr_setting
        ct_prv_prob_rst[0] = c_double(0)
        for i in range(1, target_round):
            ct_prv_prob_rst[i] =  c_double(prev_rsts[i])

        Best_Trail_Prob_Only(
            pointer(ct_prob_rst),
            c_int(target_round),
            cast(ct_prv_prob_rst, POINTER(c_double)),
        )

        self.BDP_Results = []
        for i in range(target_round):
            self.BDP_Results.append(ct_prv_prob_rst[i])
        self.BDP_Results.append(ct_prob_rst.value)

        return self.BDP_Results
        

    def get_all_BDT(self, target_round, given_prob):

        if len(given_prob) != target_round + 1:
            raise TypeError("The number of the given differential probabilities must be equal to the target round")

        ct_given_prob_rst = c_double_array(target_round + 1).array

        ct_given_prob_rst[0] = c_double(0)
        for i in range(1, target_round + 1):
            ct_given_prob_rst[i] =  c_double(given_prob[i])
        num_bdts = Best_Trail_Prob_All(target_round, ct_given_prob_rst)
        return num_bdts


    def get_DT_IO(self, gap, target_round, given_prob, instate, oustate):
        if (type(instate) == int):
            instate = self.int_to_list(instate)
        if (type(oustate) == int):
            oustate = self.int_to_list(oustate)
            
        if len(instate) != self.NUM_SBOX_IN_A_STATE:
            raise TypeError("Num of words in the input state must be equal to the number of sboxes in a state")
            
        if len(oustate) != self.NUM_SBOX_IN_A_STATE:
            raise TypeError("Num of words in the output state must be equal to the number of sboxes in a state")
        
        ct_instate = c_uint8_array(NUM_SBOX_IN_A_STATE_THRESHOLD).array
        for i in range(self.NUM_SBOX_IN_A_STATE):
            ct_instate[i] = c_uint8(instate[i])
        ct_oustate = c_uint8_array(NUM_SBOX_IN_A_STATE_THRESHOLD).array
        for i in range(self.NUM_SBOX_IN_A_STATE):
            ct_oustate[i] = c_uint8(oustate[i])

        ct_given_prob_rst = c_double_array(target_round + 1).array

        ct_given_prob_rst[0] = c_double(0)
        for i in range(1, target_round + 1):
            ct_given_prob_rst[i] =  c_double(given_prob[i])

        num_trails = Best_Trail_Prob_IO(c_double(gap), ct_instate, ct_oustate, target_round, ct_given_prob_rst)

        return num_trails

    def int_to_list(self, x, outlen = 8):

        out = [0 for i in range(outlen)] 
        tmp = x
        for idx in range(outlen):
            out[outlen - 1 - idx] = x & 0xff
            x = x >> 8
        return out

    def endian_change(self, hex_val, word_size = 8):
        new_hex_val = 0
        for idx in range(word_size):
            new_hex_val = new_hex_val ^ (hex_val & 0xff)
            if idx < (word_size - 1):
                new_hex_val = new_hex_val << 8
                hex_val = hex_val >> 8
        return new_hex_val

    def list_to_int(self, x, outlen = 8):
        hex_val = 0
        for idx in range(outlen):
            hex_val = hex_val ^ x[idx]
            if idx < outlen - 1:
                hex_val = hex_val << 8
        return hex_val
    
    def INVPERM(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        perm_out_arry = cbytes_buf(perm_out_len)

        libastbb_srbn.INV_PERM(perm_out_arry, perm_in_arry)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))

        return list(perm_out_arry)
      
    def PERM(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        perm_out_arry = cbytes_buf(perm_out_len)

        libastbb_srbn.PERM(perm_out_arry, perm_in_arry)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))


        return list(perm_out_arry)
    
    def cal_key_recov_comp_11R(self, diff_in, diff_out, prob, num_of_right_pair, c_e):
        if (type(diff_in) == list):
            diff_in = self.list_to_int(diff_in)
        if (type(diff_out) == list):
            diff_out = self.list_to_int(diff_out)
        diff_in_len = 8
        diff_out_len = 8
        diff_in_arry = int_to_cbytes(diff_in, diff_in_len)
        diff_out_arry = int_to_cbytes(diff_out, diff_out_len)

        c_prob = c_double(prob)
        c_num_of_right_pair = c_int(num_of_right_pair)
        c_c_e = c_double(c_e)

        comp = libastbb_srbn.Cal_Key_Recovery_Comp_11R(diff_in_arry, diff_out_arry, c_prob, c_num_of_right_pair, c_c_e)

        #print(f"SUB_IN  : {list(diff_in_arry)}")
        #print(f"SUB_OUT : {list(diff_out_arry)}")
        #print(f"PROB    : {c_prob}")
        #print(f"Complexity : {comp}")
        return comp
    
    def cal_key_recov_comp_12R(self, diff_in, diff_out, prob, num_of_right_pair, c_e):
        if (type(diff_in) == list):
            diff_in = self.list_to_int(diff_in)
        if (type(diff_out) == list):
            diff_out = self.list_to_int(diff_out)
        diff_in_len = 8
        diff_out_len = 8
        diff_in_arry = int_to_cbytes(diff_in, diff_in_len)
        diff_out_arry = int_to_cbytes(diff_out, diff_out_len)

        c_prob = c_double(prob)
        c_num_of_right_pair = c_int(num_of_right_pair)
        c_c_e = c_double(c_e)

        comp = libastbb_srbn.Cal_Key_Recovery_Comp_12R(diff_in_arry, diff_out_arry, c_prob, c_num_of_right_pair, c_c_e)

        #print(f"SUB_IN  : {list(diff_in_arry)}")
        #print(f"SUB_OUT : {list(diff_out_arry)}")
        #print(f"PROB    : {c_prob}")
        #print(f"Complexity : {comp}")
        return comp
        
    def cal_involved_bits(self, diff_in, diff_out):
        if (type(diff_in) == list):
            diff_in = self.list_to_int(diff_in)
        if (type(diff_out) == list):
            diff_out = self.list_to_int(diff_out)
        diff_in_len = 8
        diff_out_len = 8
        diff_in_arry = int_to_cbytes(diff_in, diff_in_len)
        diff_out_arry = int_to_cbytes(diff_out, diff_out_len)

        involved_bits = libastbb_srbn.Cal_Num_of_Involved_Bits(diff_in_arry, diff_out_arry)

        return involved_bits


class Latex_tab(SRBPN_cipher):
    def __init__(self, sbox, offset, num_sbox_in_a_state, algname):
        SRBPN_cipher.__init__(self, sbox, offset, num_sbox_in_a_state, algname)
    
    def bitmask_to_string(self, m0, m1, dc):
        bitstring = ""
        for byte0, byte1, bytedc in zip(m0, m1, dc):
            for i in range(8):  
                bit = 1 << (7-i)
                if byte1 & bit:
                    bitstring += '1'
                elif byte0 & bit:
                    bitstring += '0'
                elif bytedc & bit:
                    bitstring += '*'
                else:
                    bitstring += '?'  # undefined, for debug
        return bitstring
        
    def Cal_1r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)

        libastbb_srbn.PY_Cal_1r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_1r_BP.restype = None
        libastbb_srbn.PY_Cal_1r_BP(perm_in_arry, bp)
        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))


        return out_str

    def Cal_BP_1r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)

        libastbb_srbn.PY_Cal_BP_1r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_BP_1r_BP.restype = None
        libastbb_srbn.PY_Cal_BP_1r_BP(perm_in_arry, bp)
        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))


        return out_str


    def Cal_2r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)
        bp_ptr = pointer(bp)
        libastbb_srbn.PY_Cal_2r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_2r_BP.restype = None
        libastbb_srbn.PY_Cal_2r_BP(perm_in_arry, bp)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))

        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)


        return out_str


    def Cal_BP_2r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)
        bp_ptr = pointer(bp)
        libastbb_srbn.PY_Cal_BP_2r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_BP_2r_BP.restype = None
        libastbb_srbn.PY_Cal_BP_2r_BP(perm_in_arry, bp)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))

        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)


        return out_str

    def Cal_INV_1r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)

        libastbb_srbn.PY_Cal_INV_1r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_INV_1r_BP.restype = None
        libastbb_srbn.PY_Cal_INV_1r_BP(perm_in_arry, bp)
        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))


        return out_str

    def Cal_INV_BP_1r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)

        libastbb_srbn.PY_Cal_INV_BP_1r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_INV_BP_1r_BP.restype = None
        libastbb_srbn.PY_Cal_INV_BP_1r_BP(perm_in_arry, bp)
        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))


        return out_str


    def Cal_INV_2r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)
        bp_ptr = pointer(bp)
        libastbb_srbn.PY_Cal_INV_2r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_INV_2r_BP.restype = None
        libastbb_srbn.PY_Cal_INV_2r_BP(perm_in_arry, bp)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))

        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)


        return out_str


    def Cal_INV_BP_2r_BP(self, perm_in, output_opt = 'int'):
        if (type(perm_in) == int):
            perm_in = self.int_to_list(perm_in)
        perm_in_len = 8
        perm_out_len = 8
        perm_in_tmp = 0
        for idx in range(8):
            perm_in_tmp = perm_in_tmp ^ perm_in[idx]
            if idx !=7:
                perm_in_tmp = perm_in_tmp << 8


        perm_in_arry = int_to_cbytes(perm_in_tmp, perm_in_len)
        #perm_out_arry = cbytes_buf(perm_out_len)
        bp_must_1_mask = cbytes_buf(perm_out_len)
        bp_must_0_mask = cbytes_buf(perm_out_len)
        bp_dont_c_mask = cbytes_buf(perm_out_len)
        bp = BitPattern64(bp_must_1_mask, bp_must_0_mask, bp_dont_c_mask)
        bp_ptr = pointer(bp)
        libastbb_srbn.PY_Cal_INV_BP_2r_BP.argtypes = [POINTER(c_uint8), POINTER(BitPattern64)]
        libastbb_srbn.PY_Cal_INV_BP_2r_BP.restype = None
        libastbb_srbn.PY_Cal_INV_BP_2r_BP(perm_in_arry, bp)
        #print(list(perm_in_arry))
        #print(list(perm_out_arry))

        out_str = self.bitmask_to_string(bp.must_be_0_mask, bp.must_be_1_mask, bp.dont_care_mask)


        return out_str

def int_to_list(x, outlen = 8):

    out = [0 for i in range(outlen)] 
    tmp = x
    for idx in range(outlen):
        out[outlen - 1 - idx] = x & 0xff
        x = x >> 8
    return out

def endian_change(hex_val, word_size = 8):
    new_hex_val = 0
    for idx in range(word_size):
        new_hex_val = new_hex_val ^ (hex_val & 0xff)
        if idx < (word_size - 1):
            new_hex_val = new_hex_val << 8
            hex_val = hex_val >> 8
    return new_hex_val

def list_to_int(x, outlen = 8):
    hex_val = 0
    for idx in range(outlen):
        hex_val = hex_val ^ x[idx]
        if idx < outlen - 1:
            hex_val = hex_val << 8
    return hex_val
