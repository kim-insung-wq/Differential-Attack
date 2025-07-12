#include "astbb_srbn.h"
#include "dc_prob.h"

//(INNER) Cipher Information
char      ALG_NAME[256];
GEN_CNT_t SBOX_BIT_SIZE;
GEN_CNT_t SBOX_CARDINALITY;
SUB_WRD_t SBOX[SBOX_CARDINALITY_THRESHOLD] = { 0, };

GEN_CNT_t NUM_SBOX_IN_A_STATE;
GEN_CNT_t OFFSET[SBOX_BIT_SIZE_THRESHOLD] = { -1, };

SUB_WRD_t	SUB_WRD_MASK;
CIP_WRD_t	CIP_WRD_MASK;

PROB_t			MAX_PROB; //non-one probability
PROB_t			MIN_PROB; //non-zero probability
DP_I_O_PROB_t	DP_IOP[SBOX_CARDINALITY_THRESHOLD * SBOX_CARDINALITY_THRESHOLD];
DP_I_PROB_t		DP_IP_FO[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
DP_O_PROB_t		DP_OP_FI[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		DP_NUM_IOP_NONZERO;
SUB_CNT_t		DP_NUM_IOP_MAX;
SUB_CNT_t		DP_NUM_IP_FO_NONZERO[SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		DP_NUM_IP_FO_MAX[SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		DP_NUM_OP_FI_NONZERO[SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		DP_NUM_OP_FI_MAX[SBOX_CARDINALITY_THRESHOLD];
PROB_t			PROB_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		INV_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
SUB_CNT_t		Inv_LS_DDT[SBOX_CARDINALITY_THRESHOLD];		// Invariant linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 0 for all x.
SUB_CNT_t		Com_LS_DDT[SBOX_CARDINALITY_THRESHOLD];		// Complementary linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 1 for all x.
SUB_CNT_t		LS_DDT[SBOX_CARDINALITY_THRESHOLD];			// Invariant + Complementary e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = constant for all x.
SUB_CNT_t		Inv_LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD]; // Invariant linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 0 for all x.
SUB_CNT_t		Com_LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD]; // Complementary linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 1 for all x.
SUB_CNT_t		LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD];     // Invariant + Complementary e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = constant for all x.
SUB_CNT_t		Standard_Vec[SBOX_BIT_SIZE_THRESHOLD];

void CS_FROM_AS(CIP_STATE_t out, ANA_STATE_t in)
{
	//L_Shift     7          6          5          4          3          2          1          0 
	//  c\a     in[0]  ||  in[1]  ||  in[2]  ||  in[3]  ||  in[4]  ||  in[5]  ||  in[6]  ||  in[7]
	//out[0] =    0          0          0          0          0          0          0          0      LSB     
	//out[1] =    1          1          1          1          1          1          1          1     
	//out[.] =   ...        ...        ...        ...        ...        ...        ...        ...    
	//out[7] =    7          7          7          7          7          7          7          7      MSB
	int S_idx, X_idx;

	ANA_WRD_t S_pivot = (ANA_WRD_t)1;
	for (X_idx = 0; X_idx < SBOX_BIT_SIZE; X_idx++)
	{
		CIP_WRD_t X_pivot = ((CIP_WRD_t)1 << ((NUM_SBOX_IN_A_STATE - 1)));

		out[X_idx] = 0;
		for (S_idx = 0; S_idx < NUM_SBOX_IN_A_STATE; S_idx++)
		{
			if ((in[S_idx] & S_pivot) != 0)
			{
				out[X_idx] += X_pivot;
			}
			X_pivot = X_pivot >> 1;
		}
		S_pivot = S_pivot << 1;
	}
}

void AS_FROM_CS(ANA_STATE_t out, CIP_STATE_t in)
{
	//L_Shift     7           6            5          4           3           2           1           0 
	//  c\a     out[0]  ||  out[1]  ||  out[2]  ||  out[3]  ||  out[4]  ||  out[5]  ||  out[6]  ||  out[7]
	//in[0]  =    0           0           0           0           0           0           0           0      LSB     
	//in[1]  =    1           1           1           1           1           1           1           1     
	//in[.]  =   ...         ...         ...         ...         ...         ...         ...         ...    
	//in[7]  =    7           7           7           7           7           7           7           7      MSB
	int X_idx, S_idx;

	CIP_WRD_t X_pivot = ((CIP_WRD_t)1 << ((NUM_SBOX_IN_A_STATE - 1)));

	for (S_idx = 0; S_idx < NUM_SBOX_IN_A_STATE; S_idx++)
	{
		ANA_WRD_t S_pivot = (ANA_WRD_t)1;

		out[S_idx] = 0;
		for (X_idx = 0; X_idx < SBOX_BIT_SIZE; X_idx++)
		{
			if ((in[X_idx] & X_pivot) != 0)
			{
				out[S_idx] += S_pivot;
			}
			S_pivot = S_pivot << 1;
		}

		X_pivot = X_pivot >> 1;
	}
}


#define CIP_ROL(cip_wrd, offset) (((cip_wrd) << ((offset)%(NUM_SBOX_IN_A_STATE)))&(CIP_WRD_MASK)  | ((cip_wrd) >> ( (NUM_SBOX_IN_A_STATE) - ((offset)%(NUM_SBOX_IN_A_STATE))))&(CIP_WRD_MASK))
#define CIP_ROR(cip_wrd, offset) (((cip_wrd) >> ((offset)%(NUM_SBOX_IN_A_STATE)))&(CIP_WRD_MASK)  | ((cip_wrd) << ( (NUM_SBOX_IN_A_STATE) - ((offset)%(NUM_SBOX_IN_A_STATE))))&(CIP_WRD_MASK))

int pattern_integrity(BitPattern64* a) {
	for (int i = 0; i < SBOX_BIT_SIZE_THRESHOLD; ++i)
	{
		if ((a->must_be_1_mask[i] ^ a->must_be_0_mask[i] ^ a->dont_care_mask[i]) != 0xFF)
			return 1;
	}
	return 0;
}

void BP_PERM(BitPattern64* a)
{
	BitPattern64 b = { 0 };
	BitPattern64* b_ptr = &b;
	if (pattern_integrity(a))
	{
		printf("impossible pattern is detected.\n");
	}

	else
	{
		PERM(b_ptr->must_be_0_mask, a->must_be_0_mask);
		PERM(b_ptr->must_be_1_mask, a->must_be_1_mask);
		PERM(b_ptr->dont_care_mask, a->dont_care_mask);

	}
	*a = *b_ptr;
}
void BP_INV_PERM(BitPattern64* a)
{
	BitPattern64 b = { 0 };
	BitPattern64* b_ptr = &b;
	if (pattern_integrity(a))
	{
		printf("impossible pattern is detected.\n");
	}

	else
	{
		INV_PERM(b_ptr->must_be_0_mask, a->must_be_0_mask);
		INV_PERM(b_ptr->must_be_1_mask, a->must_be_1_mask);
		INV_PERM(b_ptr->dont_care_mask, a->dont_care_mask);

	}
	*a = *b_ptr;
}


void PERM(ANA_STATE_t out, ANA_STATE_t in)
{
	CIP_STATE_t X;
	int idx;
	CS_FROM_AS(X, in);

	for (idx = 0; idx < SBOX_BIT_SIZE; idx++)
	{
		X[idx] = CIP_ROL(X[idx], OFFSET[idx]);
	}

	AS_FROM_CS(out, X);
}

void INV_PERM(ANA_STATE_t out, ANA_STATE_t in)
{
	CIP_STATE_t X;
	int idx;
	CS_FROM_AS(X, in);

	for (idx = 0; idx < SBOX_BIT_SIZE; idx++)
	{
		X[idx] = CIP_ROR(X[idx], OFFSET[idx]);
	}

	AS_FROM_CS(out, X);
}

void PY_Cal_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_OUT = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_out : %02x\n", DIFF_OUT[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_DDT[DIFF_OUT[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_DDT[DIFF_OUT[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	
}

void PY_Cal_BP_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_OUT = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_out : %02x\n", DIFF_OUT[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_DDT[DIFF_OUT[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_DDT[DIFF_OUT[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_PERM(_nr_BP);
}



void PY_Cal_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_OUT = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_out : %02x\n", DIFF_OUT[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_DDT[DIFF_OUT[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_DDT[DIFF_OUT[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_PERM(_nr_BP); //Consider the BitPattern after 1 round.
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{
		if (_nr_BP->must_be_0_mask[idx] != 0b11111111)
		{
			_nr_BP->must_be_0_mask[idx] = 0;
			_nr_BP->must_be_1_mask[idx] = 0;
			_nr_BP->dont_care_mask[idx] = 0b11111111; // Consider the BitPattern after 2 rounds.
			num_of_inactive_sboxes = num_of_inactive_sboxes - 1;

		}
	}
}

void PY_Cal_BP_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_OUT = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_out : %02x\n", DIFF_OUT[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_DDT[DIFF_OUT[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_DDT[DIFF_OUT[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_PERM(_nr_BP); //Consider the BitPattern after 1 round.
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{
		if (_nr_BP->must_be_0_mask[idx] != 0b11111111)
		{
			_nr_BP->must_be_0_mask[idx] = 0;
			_nr_BP->must_be_1_mask[idx] = 0;
			_nr_BP->dont_care_mask[idx] = 0b11111111; // Consider the BitPattern after 2 rounds.
			num_of_inactive_sboxes = num_of_inactive_sboxes - 1;

		}
	}
	BP_PERM(_nr_BP);
}

void PY_Cal_INV_1r_BP(ANA_STATE_t _nr_SUB_IN, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_IN = { 0, };
	INV_PERM(DIFF_IN, _nr_SUB_IN);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_IN : %02x\n", DIFF_IN[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_INV_DDT[DIFF_IN[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_INV_DDT[DIFF_IN[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}

}

void PY_Cal_INV_BP_1r_BP(ANA_STATE_t _nr_SUB_IN, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_IN = { 0, };
	INV_PERM(DIFF_IN, _nr_SUB_IN);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_IN : %02x\n", DIFF_IN[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_INV_DDT[DIFF_IN[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_INV_DDT[DIFF_IN[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_INV_PERM(_nr_BP);
}



void PY_Cal_INV_2r_BP(ANA_STATE_t _nr_SUB_IN, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_IN = { 0, };
	INV_PERM(DIFF_IN, _nr_SUB_IN);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_IN : %02x\n", DIFF_IN[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_INV_DDT[DIFF_IN[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_INV_DDT[DIFF_IN[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_INV_PERM(_nr_BP); //Consider the BitPattern after 1 round.
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{
		if (_nr_BP->must_be_0_mask[idx] != 0b11111111)
		{
			_nr_BP->must_be_0_mask[idx] = 0;
			_nr_BP->must_be_1_mask[idx] = 0;
			_nr_BP->dont_care_mask[idx] = 0b11111111; // Consider the BitPattern after 2 rounds.
			num_of_inactive_sboxes = num_of_inactive_sboxes - 1;

		}
	}
}

void PY_Cal_INV_BP_2r_BP(ANA_STATE_t _nr_SUB_IN, BitPattern64* _nr_BP)
{
	int idx;
	int num_of_inactive_sboxes = NUM_SBOX_IN_A_STATE_THRESHOLD;
	ANA_STATE_t DIFF_IN = { 0, };
	INV_PERM(DIFF_IN, _nr_SUB_IN);
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{

		//printf("@@@@@\n");
		//printf("diff_IN : %02x\n", DIFF_IN[idx]);
		_nr_BP->must_be_0_mask[idx] = Inv_LS_INV_DDT[DIFF_IN[idx]];
		//printf("must_be_0_mask : %02x\n", _nr_BP->must_be_0_mask[idx]);
		_nr_BP->must_be_1_mask[idx] = Com_LS_INV_DDT[DIFF_IN[idx]];
		_nr_BP->dont_care_mask[idx] = 0b11111111 ^ _nr_BP->must_be_0_mask[idx] ^ _nr_BP->must_be_1_mask[idx];
	}
	BP_INV_PERM(_nr_BP); //Consider the BitPattern after 1 round.
	for (idx = 0; idx < NUM_SBOX_IN_A_STATE_THRESHOLD; idx++)
	{
		if (_nr_BP->must_be_0_mask[idx] != 0b11111111)
		{
			_nr_BP->must_be_0_mask[idx] = 0;
			_nr_BP->must_be_1_mask[idx] = 0;
			_nr_BP->dont_care_mask[idx] = 0b11111111; // Consider the BitPattern after 2 rounds.
			num_of_inactive_sboxes = num_of_inactive_sboxes - 1;

		}
	}
	BP_INV_PERM(_nr_BP);
}


///////////////////////////////////






void Dif_trail_fprintf(char filename[512], char opt_str[8], DC_1ROUND_CHAR_t * out_trail, int round, PROB_t tot_prob)
{
	int __FF__;
	int i;
	char ana_word_template[256] = { 0, };
	char cip_word_template[256] = { 0, };
	CIP_STATE_t cipher_wise_state;
	FILE * ofp;
	ofp = fopen(filename, opt_str);
	fprintf(ofp, "[%2d] Round Best : 2^{%0.4lf} \n", round, tot_prob);
	fprintf(ofp, "========================\n");

	i = 0;
	while (SBOX_BIT_SIZE > (i * 4))
	{
		i++;
	}
	sprintf(ana_word_template, "%%0%dX ", i);
	i = 0;
	//while (NUM_SBOX_IN_A_STATE > (i * 4))
	//{
	//	i++;
	//}
	//sprintf(cip_word_template, "%%0%dX ", i);


	for (__FF__ = 1; __FF__ <= round; __FF__++)
	{
		//SBOX IN
		CS_FROM_AS(cipher_wise_state, out_trail[__FF__].sub_i);
		//fprintf(ofp, "R %2d: SBOXIN              | ", __FF__);
		//for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
		//{
		//	fprintf(ofp, cip_word_template, cipher_wise_state[i]);
		//}
		//fprintf(ofp, "\t");

		fprintf(ofp, "R %2d: SBOXIN(Sbox-wise)   | ", __FF__);
		for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
		{
			fprintf(ofp, ana_word_template, out_trail[__FF__].sub_i[i]);
		}
		fprintf(ofp, "\n");

		//SBOX OUT
		CS_FROM_AS(cipher_wise_state, out_trail[__FF__].sub_o);

		//fprintf(ofp, "R %2d: SBOXOUT             | ", __FF__);

		//for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
		//{
		//	fprintf(ofp, cip_word_template, cipher_wise_state[i]);
		//}
		//fprintf(ofp, "\t");

		fprintf(ofp, "R %2d: SBOXOUT(Sbox-wise)  | ", __FF__);
		for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
		{
			fprintf(ofp, ana_word_template, out_trail[__FF__].sub_o[i]);
		}
		fprintf(ofp, "\n");


		//DIFF OUT
		CS_FROM_AS(cipher_wise_state, out_trail[__FF__].dif_o);
		//fprintf(ofp, "R %2d: DIFFOUT             | ", __FF__);
		//for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
		//{
		//	fprintf(ofp, cip_word_template, cipher_wise_state[i]);
		//}
		//fprintf(ofp, "\t");

		fprintf(ofp, "R %2d: DIFFOUT(Sbox-wise)  | ", __FF__);
		for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
		{
			fprintf(ofp, ana_word_template, out_trail[__FF__].dif_o[i]);
		}
		fprintf(ofp, "\n");


		fprintf(ofp, "R %2d: PROB     | 2^{%0.4lf}\n", __FF__, out_trail[__FF__].p);

		fprintf(ofp, "========================\n");
	}
	fclose(ofp);
}

void Dif_IO_fprintf(char filename[512], char opt_str[8], DC_1ROUND_CHAR_t * out_trail, int round, PROB_t tot_prob)
{
	int i;
	char ana_word_template[256] = { 0, };
	char cip_word_template[256] = { 0, };
	CIP_STATE_t cipher_wise_state;
	FILE * ofp;
	ofp = fopen(filename, opt_str);

	i = 0;
	while (SBOX_BIT_SIZE > (i * 4))
	{
		i++;
	}
	sprintf(ana_word_template, "%%0%dX", i);
	i = 0;
	while (NUM_SBOX_IN_A_STATE > (i * 4))
	{
		i++;
	}
	sprintf(cip_word_template, "%%0%dX", i);

	//1ROUND SBOX IN
	CS_FROM_AS(cipher_wise_state, out_trail[1].sub_i);
	//fprintf(ofp, "R %2d: SBOXIN              | ", 1);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		//fprintf(ofp, cip_word_template, cipher_wise_state[i]);
	}
	//fprintf(ofp, "\t");

	//fprintf(ofp, "R %2d: SBOXIN(Sbox-wise)   | ", 1);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		fprintf(ofp, ana_word_template, out_trail[1].sub_i[i]);
	}
	fprintf(ofp, "\n");

	//TARGET ROUND SBOX OUT
	CS_FROM_AS(cipher_wise_state, out_trail[round].sub_o);

	//fprintf(ofp, "R %2d: SBOXOUT             | ", round);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		//fprintf(ofp, cip_word_template, cipher_wise_state[i]);
	}
	//fprintf(ofp, "\t");

	//fprintf(ofp, "R %2d: SBOXOUT(Sbox-wise)  | ", round);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		fprintf(ofp, ana_word_template, out_trail[round].sub_o[i]);
		//fprintf(ofp, ana_word_template, out_trail[round].dif_o[i]);
	}
	fprintf(ofp, " | ");
	fprintf(ofp, "%f", tot_prob);
	fprintf(ofp, "\n");
	
	fclose(ofp);
}


void SDif_IO_fprintf(char filename[512], char opt_str[8], DC_1ROUND_CHAR_t* out_trail, int round, PROB_t tot_prob)
{
	int i;
	char ana_word_template[256] = { 0, };
	char cip_word_template[256] = { 0, };
	CIP_STATE_t cipher_wise_state;
	FILE* ofp;
	ofp = fopen(filename, opt_str);

	i = 0;
	while (SBOX_BIT_SIZE > (i * 4))
	{
		i++;
	}
	sprintf(ana_word_template, "%%0%dX", i);
	i = 0;
	while (NUM_SBOX_IN_A_STATE > (i * 4))
	{
		i++;
	}
	sprintf(cip_word_template, "%%0%dX", i);

	//1ROUND SBOX IN
	CS_FROM_AS(cipher_wise_state, out_trail[1].sub_i);
	//fprintf(ofp, "R %2d: SBOXIN              | ", 1);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		//fprintf(ofp, cip_word_template, cipher_wise_state[i]);
	}
	//fprintf(ofp, "\t");

	//fprintf(ofp, "R %2d: SBOXIN(Sbox-wise)   | ", 1);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		fprintf(ofp, ana_word_template, out_trail[1].sub_i[i]);
	}
	fprintf(ofp, "\n");

	//TARGET ROUND SBOX OUT
	CS_FROM_AS(cipher_wise_state, out_trail[round].sub_o);

	//fprintf(ofp, "R %2d: SBOXOUT             | ", round);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		//fprintf(ofp, cip_word_template, cipher_wise_state[i]);
	}
	//fprintf(ofp, "\t");

	//fprintf(ofp, "R %2d: SBOXOUT(Sbox-wise)  | ", round);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		fprintf(ofp, ana_word_template, out_trail[round].sub_o[i]);
		//fprintf(ofp, ana_word_template, out_trail[round].dif_o[i]);
	}
	fprintf(ofp, "\n");




	fclose(ofp);
}



void IO_printf(char file_in[512], char file_ou[512],  ANA_STATE_t _1r_sub_in, ANA_STATE_t _nr_sub_out, int round)
{
	int i;
	char ana_word_template[256] = { 0, };
	char ana_word_template2[256] = { 0, };
	char cip_word_template[256] = { 0, };
	CIP_STATE_t cipher_wise_state;
	ANA_STATE_t   _nr_dif_out;
	PERM(_nr_dif_out, _nr_sub_out);

	i = 0;
	while (SBOX_BIT_SIZE > (i * 4))
	{
		i++;
	}
	sprintf(ana_word_template, "%%0%dX ", i);
	sprintf(ana_word_template2, "%%0%dX", i);
	i = 0;
	while (NUM_SBOX_IN_A_STATE > (i * 4))
	{
		i++;
	}
	sprintf(cip_word_template, "%%0%dX ", i);

	//1ROUND SBOX IN
	CS_FROM_AS(cipher_wise_state, _1r_sub_in);
	
	/*
	printf("R %2d: SBOXIN              | ", 1);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		printf(cip_word_template, cipher_wise_state[i]);
	}
	printf("\t");
	*/
	
	//printf("R %2d: SBOXIN(Sbox-wise)   | ", 1);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		char tmpstr[512] = { 0, };
		//printf(ana_word_template, _1r_sub_in[i]);
		sprintf(tmpstr, ana_word_template2, _1r_sub_in[i]);
		strcat(file_in, tmpstr);
	}
	//printf("\n");

	//TARGET ROUND SBOX OUT
	CS_FROM_AS(cipher_wise_state, _nr_sub_out);

	/*
	printf("R %2d: SBOXOUT             | ", round);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		printf(cip_word_template, cipher_wise_state[i]);
	}
	printf("\t");
	*/
	//printf("R %2d: SBOXOUT(Sbox-wise)  | ", round);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		char tmpstr[512] = { 0, };
		//printf(ana_word_template, _nr_sub_out[i]);
		sprintf(tmpstr, ana_word_template2, _nr_sub_out[i]);
		strcat(file_ou, tmpstr);
	}
	//printf("\n");

	//TARGET ROUND DIFF OUT
	CS_FROM_AS(cipher_wise_state, _nr_dif_out);

	/*
	printf("R %2d: DIFFOUT             | ", round);
	for (i = (SBOX_BIT_SIZE - 1); i >= 0; i--)
	{
		printf(cip_word_template, cipher_wise_state[i]);
	}
	printf("\t");

	printf("R %2d: DIFFOUT(Sbox-wise)  | ", round);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		printf(ana_word_template, _nr_dif_out[i]);
	}
	printf("\n");
	*/
}