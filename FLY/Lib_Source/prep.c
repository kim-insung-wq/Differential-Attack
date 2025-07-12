#include "astbb_srbn.h"
#include "dc_prob.h"



////extern(decleared in global.c)

//Write here
extern char      ALG_NAME[256];
extern GEN_CNT_t SBOX_BIT_SIZE;
extern GEN_CNT_t SBOX_CARDINALITY;
extern SUB_WRD_t SBOX[SBOX_CARDINALITY_THRESHOLD];

extern GEN_CNT_t NUM_SBOX_IN_A_STATE;
extern GEN_CNT_t OFFSET[SBOX_BIT_SIZE_THRESHOLD];

extern SUB_WRD_t SUB_WRD_MASK;
extern CIP_WRD_t CIP_WRD_MASK;

extern PROB_t			MAX_PROB; //non-one probability
extern PROB_t			MIN_PROB; //non-zero probability
extern DP_I_O_PROB_t	DP_IOP[SBOX_CARDINALITY_THRESHOLD * SBOX_CARDINALITY_THRESHOLD];
extern DP_I_PROB_t		DP_IP_FO[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern DP_O_PROB_t		DP_OP_FI[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		DP_NUM_IOP_NONZERO;
extern SUB_CNT_t		DP_NUM_IOP_MAX;
extern SUB_CNT_t		DP_NUM_IP_FO_NONZERO[SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		DP_NUM_IP_FO_MAX[SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		DP_NUM_OP_FI_NONZERO[SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		DP_NUM_OP_FI_MAX[SBOX_CARDINALITY_THRESHOLD];
extern PROB_t			PROB_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t		INV_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
//Here
//SUB_CNT_t DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];

void set_global_cip_info(SRBPN_INFO_t * cipher_info)
{
	int i;
	memcpy(ALG_NAME, cipher_info->ALG_NAME, sizeof(char) * 256);
	SBOX_BIT_SIZE	 = cipher_info->SBOX_BIT_SIZE;
	SBOX_CARDINALITY = cipher_info->SBOX_CARDINALITY;
	memcpy(SBOX, cipher_info->SBOX, sizeof(SUB_WRD_t)*SBOX_CARDINALITY);

	NUM_SBOX_IN_A_STATE = cipher_info->NUM_SBOX_IN_A_STATE;
	memcpy(OFFSET, cipher_info->OFFSET, sizeof(GEN_CNT_t)*SBOX_BIT_SIZE);
	
	//INIT THE MASK
	SUB_WRD_MASK = (SUB_WRD_t)0;
	for (i = 0; i < SBOX_BIT_SIZE; i++)
	{
		SUB_WRD_MASK = SUB_WRD_MASK ^ ((SUB_WRD_t)1 << i);
	}

	//INIT THE MASK
	CIP_WRD_MASK = (CIP_WRD_t)0;
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		CIP_WRD_MASK = CIP_WRD_MASK ^ ((CIP_WRD_t)1 << i);
	}
}




//Compare Function for descending order.
int descending_P(const void *first, const void *second)
{
	PROB_t first_p = *((PROB_t*)first);
	PROB_t second_p = *((PROB_t*)second);
	return COMP_PROB(second_p, first_p);
}

//Compare Function for descending order.
int descending_IOP(const void *first, const void *second)
{
	const DP_I_O_PROB_t * t_first = first;
	const DP_I_O_PROB_t * t_second = second;
	PROB_t first_p = t_first->p;
	PROB_t second_p = t_second->p;
	return COMP_PROB(second_p, first_p);
}

//Compare Function for descending order.
int descending_OP(const void *first, const void *second)
{
	const DP_O_PROB_t * t_first = first;
	const DP_O_PROB_t * t_second = second;
	PROB_t first_p = t_first->p;
	PROB_t second_p = t_second->p;
	return COMP_PROB(second_p, first_p);
}

//Compare Function for descending order.
int descending_IP(const void *first, const void *second)
{
	const DP_I_PROB_t * t_first = first;
	const DP_I_PROB_t * t_second = second;
	PROB_t first_p = t_first->p;
	PROB_t second_p = t_second->p;
	return COMP_PROB(second_p, first_p);
}


int Compute_DDT(void)
{
	SUB_CNT_t i;
	SUB_CNT_t delta_i;
	SUB_CNT_t delta_o;

	//init
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			DDT[delta_i][delta_o] = 0;
		}
	}

	//compute DDT
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++) // difference
	{
		for (i = 0; i < SBOX_CARDINALITY; i++) // plaintext
		{
			SUB_WRD_t in = (SUB_WRD_t)i;
			SUB_WRD_t delta_o = SBOX[in] ^ SBOX[in^delta_i];
			DDT[delta_i][delta_o]++;
		}
	}

	MAX_PROB = ZERO_PROB;
	MIN_PROB = ONE_PROB;
	//compute PROB_DDT
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			SUB_CNT_t this_DDT_ent = DDT[delta_i][delta_o];

			if (this_DDT_ent == 0)
			{
				PROB_DDT[delta_i][delta_o] = (PROB_t)ZERO_PROB;
			}
			else //non-zero
			{
				if (this_DDT_ent == SBOX_CARDINALITY)
				{
					PROB_DDT[delta_i][delta_o] = (PROB_t)ONE_PROB;
				}
				else
				{
					PROB_DDT[delta_i][delta_o] = log2((PROB_t)this_DDT_ent / (PROB_t)SBOX_CARDINALITY);
				}

				if (COMP_PROB(MIN_PROB, PROB_DDT[delta_i][delta_o]) == LEFT)
				{
					MIN_PROB = PROB_DDT[delta_i][delta_o];
				}
			}
			if ((delta_i != 0) || (delta_o != 0))
			{
				if (COMP_PROB(MAX_PROB, PROB_DDT[delta_i][delta_o]) == RIGHT)
				{
					MAX_PROB = PROB_DDT[delta_i][delta_o];
				}
			}
			
		}
	}
	return 0;
}

int Compute_INV_DDT(void)
{
	SUB_WRD_t INV_SBOX[SBOX_CARDINALITY_THRESHOLD];

	SUB_CNT_t i;
	SUB_CNT_t delta_i;
	SUB_CNT_t sbox_i;
	SUB_CNT_t delta_o;
	//
	// GEN_INV_Sbox
	for (sbox_i = 0; sbox_i < SBOX_CARDINALITY; sbox_i++)
	{
		INV_SBOX[SBOX[sbox_i]] = sbox_i;
	}
	

	//init
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			INV_DDT[delta_i][delta_o] = 0;
		}
	}

	//compute DDT
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++) // difference
	{
		for (i = 0; i < SBOX_CARDINALITY; i++) // plaintext
		{
			SUB_WRD_t in = (SUB_WRD_t)i;
			SUB_WRD_t delta_o = INV_SBOX[in] ^ INV_SBOX[in ^ delta_i];
			INV_DDT[delta_i][delta_o]++;
		}
	}
	return 0;
}



int Prep_DDT_Info(void)
{
	SUB_CNT_t delta_i;
	SUB_CNT_t delta_o;
	int idx;
	PROB_t the_prob;

	// DDT information rearrange
	// all I_O_PROB
	// O_PROB with fixed input
	// I_PROB with fixed output
	//rearrange
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			//all I_O_PROB
			the_prob = PROB_DDT[delta_i][delta_o];
			DP_IOP[delta_i * SBOX_CARDINALITY + delta_o].i = (SUB_WRD_t)delta_i;
			DP_IOP[delta_i * SBOX_CARDINALITY + delta_o].o = (SUB_WRD_t)delta_o;
			DP_IOP[delta_i * SBOX_CARDINALITY + delta_o].p = (PROB_t)the_prob;
			//O_PROB with fixed input
			DP_OP_FI[delta_i][delta_o].o = (SUB_WRD_t)delta_o;
			DP_OP_FI[delta_i][delta_o].p = (PROB_t)the_prob;
			//I_PROB with fixed output
			DP_IP_FO[delta_o][delta_i].i = (SUB_WRD_t)delta_i;
			DP_IP_FO[delta_o][delta_i].p = (PROB_t)the_prob;
		}
	}

	//sorting all I_O_PROB
	qsort((void*)DP_IOP, SBOX_CARDINALITY * SBOX_CARDINALITY, sizeof(DP_I_O_PROB_t), descending_IOP);

	//sorting O_PROB with fixed input
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		qsort((void*)DP_OP_FI[delta_i], SBOX_CARDINALITY, sizeof(DP_O_PROB_t), descending_OP);
	}
	//sorting I_PROB with fixed output
	for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
	{
		qsort((void*)DP_IP_FO[delta_o], SBOX_CARDINALITY, sizeof(DP_I_PROB_t), descending_IP);
	}

	//Counting the numbers of each information
	/*
	*  1) The number of (Input,Output) with NON-zero probability.
	*  2) The number of (Input,Output) with the best probability.
	*  3) The number of (Input)/(Output) with NON-zero probability and each fixed (Output)/(Input).
	*  4) The number of (Input)/(Output) with the best probability and each fixed (Output)/(Input).
	*/

	////init
	DP_NUM_IOP_NONZERO = 0;
	DP_NUM_IOP_MAX = 0;
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		DP_NUM_OP_FI_NONZERO[delta_i] = 0;
		DP_NUM_OP_FI_MAX[delta_i] = 0;
	}
	for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
	{
		DP_NUM_IP_FO_NONZERO[delta_o] = 0;
		DP_NUM_IP_FO_MAX[delta_o] = 0;
	}

	////Counting IO_PROB
	for (idx = 0; idx < SBOX_CARDINALITY * SBOX_CARDINALITY; idx++)
	{
		PROB_t IO_prob = DP_IOP[idx].p;
		PROB_t max_IO_prob = DP_IOP[0].p;
		if (COMP_PROB(IO_prob, ZERO_PROB) != EQUAL)
		{
			DP_NUM_IOP_NONZERO++;
		}
		if (COMP_PROB(max_IO_prob, IO_prob) == EQUAL)
		{
			DP_NUM_IOP_MAX++;
		}
	}

	////Counting O_PROB
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (idx = 0; idx < SBOX_CARDINALITY; idx++)
		{
			PROB_t O_prob = DP_OP_FI[delta_i][idx].p;
			PROB_t max_O_prob = DP_OP_FI[delta_i][0].p;
			if (COMP_PROB(O_prob, ZERO_PROB) != EQUAL)
			{
				DP_NUM_OP_FI_NONZERO[delta_i]++;
			}
			if (COMP_PROB(max_O_prob, O_prob) == EQUAL)
			{
				DP_NUM_OP_FI_MAX[delta_i]++;
			}
		}
	}

	////Counting I_PROB
	for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
	{
		for (idx = 0; idx < SBOX_CARDINALITY; idx++)
		{
			PROB_t I_prob = DP_IP_FO[delta_o][idx].p;
			PROB_t max_I_prob = DP_IP_FO[delta_o][0].p;
			if (COMP_PROB(I_prob, ZERO_PROB) != EQUAL)
			{
				DP_NUM_IP_FO_NONZERO[delta_o]++;
			}
			if (COMP_PROB(max_I_prob, I_prob) == EQUAL)
			{
				DP_NUM_IP_FO_MAX[delta_o]++;
			}
		}
	}

	return 0;
}


void Prep_Dif_Trail_Searching(SRBPN_INFO_t * cipher_info)
{
	set_global_cip_info(cipher_info);
	Compute_DDT();
	//Compute_INV_DDT();
	Prep_DDT_Info();
}

