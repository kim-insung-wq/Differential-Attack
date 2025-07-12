#include "astbb_srbn.h"
#include "dc_prob.h"


void PERM(ANA_STATE_t out, ANA_STATE_t in);

//Cipher Information : read-only here(declared in global.c)
extern char      ALG_NAME[256];
extern GEN_CNT_t SBOX_BIT_SIZE;
extern GEN_CNT_t SBOX_CARDINALITY;
extern GEN_CNT_t NUM_SBOX_IN_A_STATE;
//DDT analysis results : read-only here(declared in global.c)
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

#include "active_map.h"

void Dif_IO_fprintf(char filename[512], char opt_str[8], DC_1ROUND_CHAR_t* out_trail, int round, PROB_t tot_prob);

//used here only
static PROB_t* DP_BOUNDS = NULL;
static DC_1ROUND_CHAR_t* DIFF_TRAIL_IN_PROG = NULL;
static DC_1ROUND_CHAR_t* DIFF_TRAIL_FOR_OUT = NULL;
static PROB_t				DP_BOUND_IN_PROG;
static FLAG_t				TOUCH_THE_LEAF;
static char					OUT_FILE_NAME[512];
static NUM_TRAIL_t			TOTAL_NUM_TRAILS = 0;
static WEIGHT_t COMP_THRESHOLD = 100;

static DEV_INLINE PROB_t put_probs(PROB_t* p)
{
	int i;
	PROB_t P;
	P = p[0];
	for (i = 1; i < NUM_SBOX_IN_A_STATE; i++)
	{
		MUL_PROB(&P, P, p[i]);
	}
	return P;
}

static DEV_INLINE void init_probs(PROB_t* p)
{
	int i;
	for (i = 0; i < NUM_SBOX_IN_A_STATE_THRESHOLD; i++)
	{
		p[i] = ONE_PROB;
	}
}


static DEV_INLINE int bound_checker_prob(PROB_t check)
{
	if (COMP_PROB(check, DP_BOUND_IN_PROG) == RIGHT)
		return EXCEED_BOUND;
	else
		return UNDER_BOUND;
}

static DEV_INLINE int expect_state_prob_at_1_j(PROB_t* expected_prob, TRUNC_STATE_t AT, int j)
{
	int first_active = THE_LAST_WORD;
	PROB_t tmp_cum_prob = ONE_PROB;
	int wrd_idx;
	if (j == ROUND_START)
	{
		for (wrd_idx = 0; wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
		{
			if ((AT & ((TRUNC_STATE_t)0x1 << wrd_idx)) != 0)
			{
				MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
				if (first_active == THE_LAST_WORD)
				{
					first_active = wrd_idx;
				}
			}
		}
	}
	else if (j < (NUM_SBOX_IN_A_STATE - 1))
	{
		for (wrd_idx = (j + 1); wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
		{
			if ((AT & ((TRUNC_STATE_t)0x1 << wrd_idx)) != 0)
			{
				MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
				if (first_active == THE_LAST_WORD)
				{
					first_active = wrd_idx;
				}
			}
		}
	}
	else // if (j == (NUM_SBOX_IN_A_STATE - 1))
	{
		//return initial out
	}
	*expected_prob = tmp_cum_prob;
	return first_active;
}



static DEV_INLINE int expect_state_prob_at_j(PROB_t* expected_prob, ANA_STATE_t X, int j)
{
	int first_active = THE_LAST_WORD;
	PROB_t tmp_cum_prob = ONE_PROB;
	int wrd_idx;
	if (j == ROUND_START)
	{
		for (wrd_idx = 0; wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
		{
			if (X[wrd_idx] != 0)
			{
				MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
				if (first_active == THE_LAST_WORD)
				{
					first_active = wrd_idx;
				}
			}
		}
	}
	else if (j < (NUM_SBOX_IN_A_STATE - 1))
	{
		for (wrd_idx = (j + 1); wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
		{
			if (X[wrd_idx] != 0)
			{
				MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
				if (first_active == THE_LAST_WORD)
				{
					first_active = wrd_idx;
				}
			}
		}
	}
	else // if (j == (NUM_SBOX_IN_A_STATE - 1))
	{
		//return initial out
	}
	*expected_prob = tmp_cum_prob;
	return first_active;
}

static DEV_INLINE void expect_next_state_prob_at_j(PROB_t* expected_next_prob, ANA_STATE_t Y, int j)
{
	ANA_STATE_t copied = { 0, };
	ANA_STATE_t dif_part;
	PROB_t tmp_cum_prob = ONE_PROB;
	int wrd_idx;
	memcpy(copied, Y, sizeof(SUB_WRD_t) * (j + 1));

	PERM(dif_part, copied);

	for (wrd_idx = 0; wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
	{
		if (dif_part[wrd_idx] != 0)
		{
			MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
		}
	}

	*expected_next_prob = tmp_cum_prob;
}


static void recursive_construct_state(ANA_STATE_t in_var_state, ANA_STATE_t ou_var_state, ANA_STATE_t in_pivot_state, ANA_STATE_t ou_pivot_state, int wrd_idx, int target_round)
{
	if (wrd_idx == NUM_SBOX_IN_A_STATE * 2)
	{
		//done and print
		memcpy(DIFF_TRAIL_FOR_OUT[1].sub_i, in_var_state, sizeof(ANA_STATE_t));
		memcpy(DIFF_TRAIL_FOR_OUT[target_round].sub_o, ou_var_state, sizeof(ANA_STATE_t));
		PERM(DIFF_TRAIL_FOR_OUT[target_round].dif_o, DIFF_TRAIL_FOR_OUT[target_round].sub_o);
		TOTAL_NUM_TRAILS++;
		PROB_t current_trail_prob = DIFF_TRAIL_FOR_OUT[target_round].p;
		Dif_IO_fprintf(OUT_FILE_NAME, "a", DIFF_TRAIL_FOR_OUT, target_round, current_trail_prob);
	}
	else
	{
		if (wrd_idx < NUM_SBOX_IN_A_STATE)
		{
			//determine in_var_state[wrd_idx] 
			int in_wrd = wrd_idx;
			ANA_WRD_t fix_sbox_out = in_pivot_state[in_wrd];
			if (fix_sbox_out == 0)
			{
				in_var_state[in_wrd] = 0;
				recursive_construct_state(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round);
			}
			else
			{
				int max_idx;
				for (max_idx = 0; max_idx < DP_NUM_IP_FO_MAX[fix_sbox_out]; max_idx++)
				{
					in_var_state[in_wrd] = DP_IP_FO[fix_sbox_out][max_idx].i;
					recursive_construct_state(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round);
				}
			}
		}
		else
		{
			//determine ou_var_state[wrd_idx - NUM_SBOX_IN_A_STATE]
			int ou_wrd = wrd_idx - NUM_SBOX_IN_A_STATE;
			ANA_WRD_t fix_sbox_in = ou_pivot_state[ou_wrd];
			if (fix_sbox_in == 0)
			{
				ou_var_state[ou_wrd] = 0;
				recursive_construct_state(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round);
			}
			else
			{
				int max_idx;
				for (max_idx = 0; max_idx < DP_NUM_OP_FI_MAX[fix_sbox_in]; max_idx++)
				{
					ou_var_state[ou_wrd] = DP_OP_FI[fix_sbox_in][max_idx].o;
					recursive_construct_state(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round);
				}
			}
		}



	}
}


static void recursive_construct_state_with_prob(ANA_STATE_t in_var_state, ANA_STATE_t ou_var_state, ANA_STATE_t in_pivot_state, ANA_STATE_t ou_pivot_state, int wrd_idx, int target_round, PROB_t tot_prob)
{
	if (wrd_idx == NUM_SBOX_IN_A_STATE * 2)
	{
		//done and print
		memcpy(DIFF_TRAIL_FOR_OUT[1].sub_i, in_var_state, sizeof(ANA_STATE_t));
		memcpy(DIFF_TRAIL_FOR_OUT[target_round].sub_o, ou_var_state, sizeof(ANA_STATE_t));
		PERM(DIFF_TRAIL_FOR_OUT[target_round].dif_o, DIFF_TRAIL_FOR_OUT[target_round].sub_o);
		TOTAL_NUM_TRAILS++;
		Dif_IO_fprintf(OUT_FILE_NAME, "a", DIFF_TRAIL_FOR_OUT, target_round, tot_prob);
	}
	else
	{
		if (wrd_idx < NUM_SBOX_IN_A_STATE)
		{
			//determine in_var_state[wrd_idx] 
			int in_wrd = wrd_idx;
			ANA_WRD_t fix_sbox_out = in_pivot_state[in_wrd];
			if (fix_sbox_out == 0)
			{
				in_var_state[in_wrd] = 0;
				recursive_construct_state_with_prob(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round, tot_prob);
			}
			else
			{
				int max_idx;
				for (max_idx = 0; max_idx < DP_NUM_IP_FO_MAX[fix_sbox_out]; max_idx++)
				{
					in_var_state[in_wrd] = DP_IP_FO[fix_sbox_out][max_idx].i;
					recursive_construct_state_with_prob(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round, tot_prob);
				}
			}
		}
		else
		{
			//determine ou_var_state[wrd_idx - NUM_SBOX_IN_A_STATE]
			int ou_wrd = wrd_idx - NUM_SBOX_IN_A_STATE;
			ANA_WRD_t fix_sbox_in = ou_pivot_state[ou_wrd];
			if (fix_sbox_in == 0)
			{
				ou_var_state[ou_wrd] = 0;
				recursive_construct_state_with_prob(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round, tot_prob);
			}
			else
			{
				int max_idx;
				for (max_idx = 0; max_idx < DP_NUM_OP_FI_MAX[fix_sbox_in]; max_idx++)
				{
					ou_var_state[ou_wrd] = DP_OP_FI[fix_sbox_in][max_idx].o;
					recursive_construct_state_with_prob(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, wrd_idx + 1, target_round, tot_prob);
				}
			}
		}



	}
}





static void expansion_to_available_trails(int target_round)
{
	ANA_STATE_t in_pivot_state, ou_pivot_state;
	ANA_STATE_t in_var_state, ou_var_state;
	//input pivot state : 1-round SUB OUT
	memcpy(in_pivot_state, DIFF_TRAIL_FOR_OUT[1].sub_o, sizeof(ANA_STATE_t));

	//output pivot state : target_round-round SUB_IN
	memcpy(ou_pivot_state, DIFF_TRAIL_FOR_OUT[target_round].sub_i, sizeof(ANA_STATE_t));

	//start constructing
	recursive_construct_state(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, 0, target_round);
	//printf("Num of Found Best Trails : %llu\r", TOTAL_NUM_TRAILS);
}

static void expansion_to_available_trails_with_prob(int target_round, PROB_t tot_prob)
{
	ANA_STATE_t in_pivot_state, ou_pivot_state;
	ANA_STATE_t in_var_state, ou_var_state;
	//input pivot state : 1-round SUB OUT
	memcpy(in_pivot_state, DIFF_TRAIL_FOR_OUT[1].sub_o, sizeof(ANA_STATE_t));

	//output pivot state : target_round-round SUB_IN
	memcpy(ou_pivot_state, DIFF_TRAIL_FOR_OUT[target_round].sub_i, sizeof(ANA_STATE_t));

	//start constructing
	recursive_construct_state_with_prob(in_var_state, ou_var_state, in_pivot_state, ou_pivot_state, 0, target_round, tot_prob);
	//printf("Num of Found Best Trails : %llu\r", TOTAL_NUM_TRAILS);
}





NUM_TRAIL_t Best_Trail_Prob_All(GEN_CNT_t target_round, PROB_t* given_round_bdp);

void SRBPN_Round1_Prob_All(int target_round);
void SRBPN_Round1_j_Prob_All(ANA_STATE_t X1, ANA_STATE_t Y1, PROB_t P1[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, TRUNC_STATE_t AT1, int j);

void SRBPN_Roundi_Prob_All(PROB_t BEF_DET_PROB, int target_round, int i);
void SRBPN_Roundi_j_Prob_All(ANA_STATE_t Xi, ANA_STATE_t Yi, PROB_t Pi[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, int i, int j);

void SRBPN_Roundn_Prob_All(PROB_t BEF_DET_PROB, int n);
void SRBPN_Roundn_j_Prob_All(ANA_STATE_t Xn, ANA_STATE_t Yn, PROB_t Pn[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n, int j);


NUM_TRAIL_t Best_Trail_Prob_All(GEN_CNT_t target_round, PROB_t* given_round_bdp)
{
	int round_idx;
	FILE* ofp;

	//Step 0 - Prep
	DP_BOUNDS = (PROB_t*)malloc(sizeof(PROB_t) * (target_round + 1));
	DIFF_TRAIL_IN_PROG = (DC_1ROUND_CHAR_t*)malloc(sizeof(DC_1ROUND_CHAR_t) * (target_round + 1));
	DIFF_TRAIL_FOR_OUT = (DC_1ROUND_CHAR_t*)malloc(sizeof(DC_1ROUND_CHAR_t) * (target_round + 1));
	sprintf(OUT_FILE_NAME, "%s_All_BDTs_R%d.txt", ALG_NAME, target_round);
	ofp = fopen(OUT_FILE_NAME, "w"); //reset
	fclose(ofp);
	TOTAL_NUM_TRAILS = 0;

	//Step 1 - Set Best Differential Probabilities of Round BDTs
	if (target_round == 1)
	{
		goto done;//we don't have to search the trail(trivial case)
	}
	else// target_round >=2
	{
		for (round_idx = 1; round_idx <= target_round; round_idx++)
		{
			DP_BOUNDS[round_idx] = given_round_bdp[round_idx];
		}
	}
	//Start searching (Looooook! Only target round is applied)
	DP_BOUND_IN_PROG = DP_BOUNDS[target_round];


	//Step 2(main) - Start searching
	/*Finding Trail*/
	SRBPN_Round1_Prob_All(target_round);
	/***************/

done:
	free(DP_BOUNDS);			DP_BOUNDS = NULL;
	free(DIFF_TRAIL_IN_PROG);	DIFF_TRAIL_IN_PROG = NULL;
	free(DIFF_TRAIL_FOR_OUT);	DIFF_TRAIL_FOR_OUT = NULL;

	return TOTAL_NUM_TRAILS;
}



void SRBPN_Round1_Prob_All(int target_round)
{
	ANA_STATE_t		X1;
	ANA_STATE_t		Y1;
	PROB_t			P1[NUM_SBOX_IN_A_STATE_THRESHOLD];
	GEN_CNT_t		NUM_AT1;
	TRUNC_STATE_t	AT1;
	PROB_t			_1_round_expected;
	PROB_t			_1_bound;
	int				next_word_idx;
	PROB_t			DET_PROB = ONE_PROB;

	for (NUM_AT1 = 1; NUM_AT1 <= NUM_SBOX_IN_A_STATE; NUM_AT1++)
	{
		FLAG_t is_last;
		//Active Map Setting
		Init_AT1(NUM_AT1);
		do
		{
			is_last = Next_AT1(&AT1, TRUE);
			next_word_idx = expect_state_prob_at_1_j(&_1_round_expected, AT1, ROUND_START); //1R
			MUL_PROB(&_1_bound, _1_round_expected, DP_BOUNDS[target_round - 1]);			//[2~n]R

			//pruning
			if (bound_checker_prob(_1_bound) == EXCEED_BOUND)
				return; //we can get any better bound due to num_at1 is sorted by ascending order

			//init
			memset(X1, 0, sizeof(ANA_STATE_t));
			memset(Y1, 0, sizeof(ANA_STATE_t));
			init_probs(P1);

			SRBPN_Round1_j_Prob_All(X1, Y1, P1, DET_PROB, target_round, AT1, next_word_idx);

		} while (is_last == THIS_IS_NOT_THE_LAST);
	}
}


void SRBPN_Round1_j_Prob_All(ANA_STATE_t X1, ANA_STATE_t Y1, PROB_t P1[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, TRUNC_STATE_t AT1, int j)
{
	SUB_CNT_t	_1_j;
	SUB_WRD_t	X1_j;
	SUB_WRD_t	Y1_j;
	PROB_t		P1_j;
	int			next_word_idx;
	PROB_t		_1_round_expected;
	PROB_t		_2_round_expected;
	PROB_t		_1_bound;
	PROB_t      _2_bound;
	PROB_t		DET_PROB;

	// this word is always active
	for (_1_j = 0x1; _1_j < SBOX_CARDINALITY; _1_j++)
	{
		X1_j = DP_IP_FO[_1_j][0].i;
		Y1_j = (SUB_WRD_t)_1_j;
		P1_j = DP_IP_FO[_1_j][0].p;

		MUL_PROB(&DET_PROB, BEF_DET_PROB, P1_j);

		next_word_idx = expect_state_prob_at_1_j(&_1_round_expected, AT1, j);
		MUL_PROB(&_1_bound, _1_round_expected, DET_PROB);			//1R
		MUL_PROB(&_1_bound, _1_bound, DP_BOUNDS[target_round - 1]); //[2~n]R

		//pruning
		if (bound_checker_prob(_1_bound) == EXCEED_BOUND)
			continue; //need to check next value on j word

		//store current word
		X1[j] = X1_j;
		Y1[j] = Y1_j;
		P1[j] = P1_j;


		MUL_PROB(&_2_bound, _1_round_expected, DET_PROB);				//1R
		expect_next_state_prob_at_j(&_2_round_expected, Y1, j);
		MUL_PROB(&_2_bound, _2_round_expected, _2_bound);				//2R

		if (target_round > 2)
		{
			MUL_PROB(&_2_bound, _2_bound, DP_BOUNDS[target_round - 2]); //[3~n]R
		}

		//pruning
		if (bound_checker_prob(_2_bound) == EXCEED_BOUND)
			continue; //need to check next value on j word

		if (next_word_idx == THE_LAST_WORD) //when this word is the last word
		{
			memcpy(DIFF_TRAIL_IN_PROG[1].sub_i, X1, sizeof(ANA_STATE_t));
			memcpy(DIFF_TRAIL_IN_PROG[1].sub_o, Y1, sizeof(ANA_STATE_t));
			PERM(DIFF_TRAIL_IN_PROG[1].dif_o, DIFF_TRAIL_IN_PROG[1].sub_o);
			DIFF_TRAIL_IN_PROG[1].p = put_probs(P1);

			if (target_round == 2)
			{
				//Move to the last round
				SRBPN_Roundn_Prob_All(DET_PROB, target_round);
			}
			else
			{
				//Move to next(2) round
				SRBPN_Roundi_Prob_All(DET_PROB, target_round, 2);
			}
		}
		else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
		{
			//Move to next word
			SRBPN_Round1_j_Prob_All(X1, Y1, P1, DET_PROB, target_round, AT1, next_word_idx);
		}
	}
}

void SRBPN_Roundi_Prob_All(PROB_t BEF_DET_PROB, int target_round, int i)
{
	ANA_STATE_t Xi;
	ANA_STATE_t Yi;
	PROB_t		Pi[NUM_SBOX_IN_A_STATE_THRESHOLD];
	PROB_t		_i_round_expected;
	PROB_t		_i_bound;
	int			next_word_idx;

	memcpy(Xi, DIFF_TRAIL_IN_PROG[i - 1].dif_o, sizeof(ANA_STATE_t));

	_i_bound = BEF_DET_PROB;														//[1~(i-1)]R
	next_word_idx = expect_state_prob_at_j(&_i_round_expected, Xi, ROUND_START);	//iR
	MUL_PROB(&_i_bound, _i_bound, _i_round_expected);
	MUL_PROB(&_i_bound, _i_bound, DP_BOUNDS[target_round - i]);						//[(i+1)~n]R

	if (bound_checker_prob(_i_bound) != EXCEED_BOUND)
	{
		//init the sub out state and probabilities
		memset(Yi, 0, sizeof(ANA_STATE_t));
		init_probs(Pi);

		SRBPN_Roundi_j_Prob_All(Xi, Yi, Pi, BEF_DET_PROB, target_round, i, next_word_idx);
	}
}


void SRBPN_Roundi_j_Prob_All(ANA_STATE_t Xi, ANA_STATE_t Yi, PROB_t Pi[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, int i, int j)
{
	SUB_CNT_t	_i_j;
	SUB_WRD_t	Xi_j = Xi[j]; //It is fixed from before
	SUB_WRD_t	Yi_j;
	PROB_t		Pi_j;
	int			next_word_idx;
	PROB_t		_i_round_expected;
	PROB_t		_i1_round_expected;
	PROB_t		_i_bound;
	PROB_t		_i1_bound;
	PROB_t		DET_PROB;

	for (_i_j = 0; _i_j < DP_NUM_OP_FI_NONZERO[Xi_j]; _i_j++)
	{
		Yi_j = DP_OP_FI[Xi_j][_i_j].o;
		Pi_j = DP_OP_FI[Xi_j][_i_j].p;

		MUL_PROB(&DET_PROB, BEF_DET_PROB, Pi_j);


		next_word_idx = expect_state_prob_at_j(&_i_round_expected, Xi, j);
		MUL_PROB(&_i_bound, _i_round_expected, DET_PROB);						//[1~i]R
		MUL_PROB(&_i_bound, _i_bound, DP_BOUNDS[target_round - i]);				//[i+1~n]R

		//pruning by using descending order
		if (bound_checker_prob(_i_bound) == EXCEED_BOUND)
			break;

		//store current word
		Yi[j] = Yi_j;
		Pi[j] = Pi_j;


		MUL_PROB(&_i1_bound, _i_round_expected, DET_PROB);							//[1~i]R
		expect_next_state_prob_at_j(&_i1_round_expected, Yi, j);
		MUL_PROB(&_i1_bound, _i1_bound, _i1_round_expected);						//i+1R

		if (i < (target_round - 1))
		{
			MUL_PROB(&_i1_bound, _i1_bound, DP_BOUNDS[target_round - 1 - i]);		//[i+2~n]R
		}


		//pruning
		if (bound_checker_prob(_i1_bound) == EXCEED_BOUND)
			continue; //next value

		if (next_word_idx == THE_LAST_WORD)  //when this word is the last word
		{
			memcpy(DIFF_TRAIL_IN_PROG[i].sub_i, Xi, sizeof(ANA_STATE_t));
			memcpy(DIFF_TRAIL_IN_PROG[i].sub_o, Yi, sizeof(ANA_STATE_t));
			PERM(DIFF_TRAIL_IN_PROG[i].dif_o, DIFF_TRAIL_IN_PROG[i].sub_o);
			DIFF_TRAIL_IN_PROG[i].p = put_probs(Pi);
			if (i == (target_round - 1))
			{
				SRBPN_Roundn_Prob_All(DET_PROB, target_round);
			}
			else if (i < (target_round - 1))
			{
				SRBPN_Roundi_Prob_All(DET_PROB, target_round, i + 1);
			}
		}
		else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
		{
			//Move to next word
			SRBPN_Roundi_j_Prob_All(Xi, Yi, Pi, DET_PROB, target_round, i, next_word_idx);
		}
	}
}


void SRBPN_Roundn_Prob_All(PROB_t BEF_DET_PROB, int n)
{
	ANA_STATE_t Xn;
	ANA_STATE_t Yn;
	PROB_t		Pn[NUM_SBOX_IN_A_STATE_THRESHOLD];
	PROB_t		_n_round_expected;
	PROB_t		_n_bound;
	int			next_word_idx;

	memcpy(Xn, DIFF_TRAIL_IN_PROG[n - 1].dif_o, sizeof(ANA_STATE_t));

	_n_bound = BEF_DET_PROB;														//[1~(n-1)]R
	next_word_idx = expect_state_prob_at_j(&_n_round_expected, Xn, ROUND_START);	//nR
	MUL_PROB(&_n_bound, _n_bound, _n_round_expected);

	if (bound_checker_prob(_n_bound) != EXCEED_BOUND)
	{
		//init the sub out state and probabilities
		memset(Yn, 0, sizeof(ANA_STATE_t));
		init_probs(Pn);

		SRBPN_Roundn_j_Prob_All(Xn, Yn, Pn, BEF_DET_PROB, n, next_word_idx);
	}
}


void SRBPN_Roundn_j_Prob_All(ANA_STATE_t Xn, ANA_STATE_t Yn, PROB_t Pn[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n, int j)
{
	SUB_WRD_t	Xn_j = Xn[j]; //It is fixed from before
	SUB_WRD_t	Yn_j = DP_OP_FI[Xn_j][0].o;
	PROB_t		Pn_j = DP_OP_FI[Xn_j][0].p; //max prob with fixed in
	int			next_word_idx;
	PROB_t		_n_round_expected;
	PROB_t		_n_bound;
	PROB_t		DET_PROB;

	MUL_PROB(&DET_PROB, BEF_DET_PROB, Pn_j);
	next_word_idx = expect_state_prob_at_j(&_n_round_expected, Xn, j);
	MUL_PROB(&_n_bound, _n_round_expected, DET_PROB);						//[1~n]R

	//pruning
	if (bound_checker_prob(_n_bound) == EXCEED_BOUND)
		return;

	Yn[j] = Yn_j;
	Pn[j] = Pn_j;


	if (next_word_idx == THE_LAST_WORD)  //when this word is the last word
	{
		int L;
		WEIGHT_t w_key_recov_comp;
		
		//Copy 1~n-1-Round states to Out trail 
		for (L = 1; L <= n - 1; L++)
		{
			memcpy(DIFF_TRAIL_FOR_OUT[L].sub_i, DIFF_TRAIL_IN_PROG[L].sub_i, sizeof(ANA_STATE_t));
			memcpy(DIFF_TRAIL_FOR_OUT[L].sub_o, DIFF_TRAIL_IN_PROG[L].sub_o, sizeof(ANA_STATE_t));
			memcpy(DIFF_TRAIL_FOR_OUT[L].dif_o, DIFF_TRAIL_IN_PROG[L].dif_o, sizeof(ANA_STATE_t));
			DIFF_TRAIL_FOR_OUT[L].p = DIFF_TRAIL_IN_PROG[L].p;
		}
		//Copy n-Round state to Out trail
		memcpy(DIFF_TRAIL_FOR_OUT[n].sub_i, Xn, sizeof(ANA_STATE_t));
		memcpy(DIFF_TRAIL_FOR_OUT[n].sub_o, Yn, sizeof(ANA_STATE_t));
		PERM(DIFF_TRAIL_FOR_OUT[n].dif_o, Yn);
		DIFF_TRAIL_FOR_OUT[n].p = put_probs(Pn);

			//Expansion the current trail to all available(constructable) trails
			expansion_to_available_trails_with_prob(n, DET_PROB);
			

	}
	else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
	{
		//Move to next word
		SRBPN_Roundn_j_Prob_All(Xn, Yn, Pn, DET_PROB, n, next_word_idx);
	}
}











