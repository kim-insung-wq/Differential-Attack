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

extern PROB_t			PROB_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];

#include "active_map.h"


void IO_printf(char file_in[512], char file_ou[512], ANA_STATE_t _1r_sub_in, ANA_STATE_t _nr_sub_out, int round);
void Dif_trail_fprintf(char filename[512], char opt_str[8], DC_1ROUND_CHAR_t* out_trail, int round, PROB_t tot_prob);

//used here only
static PROB_t* DP_BOUNDS = NULL;
static DC_1ROUND_CHAR_t* DIFF_TRAIL_IN_PROG = NULL;
static DC_1ROUND_CHAR_t* DIFF_TRAIL_FOR_OUT = NULL;
static PROB_t				DP_BOUND_IN_PROG;
static FLAG_t				TOUCH_THE_LEAF;
static char					OUT_FILE_NAME[512];
static uint64_t				TOTAL_NUM_TRAILS = 0;
static ANA_STATE_t			_1R_SUB_IN;
static ANA_STATE_t			_nR_SUB_OUT;


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

static DEV_INLINE FLAG_t expect_next_state_prob_at_n_1_j(PROB_t* expected_next_prob, ANA_STATE_t Y, int j)
{
	ANA_STATE_t copied = { 0, };
	ANA_STATE_t dif_part;
	PROB_t tmp_cum_prob = ONE_PROB;
	int wrd_idx;
	memcpy(copied, Y, sizeof(SUB_WRD_t) * (j + 1));

	PERM(dif_part, copied);

	//Check Last Round Active Truncation As well
	for (wrd_idx = 0; wrd_idx < NUM_SBOX_IN_A_STATE; wrd_idx++)
	{
		if (dif_part[wrd_idx] != 0)
		{
			if (_nR_SUB_OUT[wrd_idx] == 0)
			{
				return FALSE;
			}
			MUL_PROB(&tmp_cum_prob, tmp_cum_prob, MAX_PROB);
		}
	}

	*expected_next_prob = tmp_cum_prob;
	return TRUE;
}


NUM_TRAIL_t Best_Trail_Prob_IO(PROB_t bound_gap, ANA_STATE_t _1r_sub_in, ANA_STATE_t _nr_sub_out, GEN_CNT_t target_round, PROB_t* given_round_bdp);

void SRBPN_Roundi_Prob_IO(PROB_t BEF_DET_PROB, int target_round, int i);
void SRBPN_Roundi_j_Prob_IO(ANA_STATE_t Xi, ANA_STATE_t Yi, PROB_t Pi[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, int i, int j);

void SRBPN_Roundn_1_Prob_IO(PROB_t BEF_DET_PROB, int n_1);
void SRBPN_Roundn_1_j_Prob_IO(ANA_STATE_t Xn_1, ANA_STATE_t Yn_1, PROB_t Pn_1[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n_1, int j);

void SRBPN_Roundn_Prob_IO(PROB_t BEF_DET_PROB, int n);
void SRBPN_Roundn_j_Prob_IO(ANA_STATE_t Xn, ANA_STATE_t Yn, PROB_t Pn[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n, int j);


NUM_TRAIL_t Best_Trail_Prob_IO(PROB_t bound_gap, ANA_STATE_t _1r_sub_in, ANA_STATE_t _nr_sub_out, GEN_CNT_t target_round, PROB_t* given_round_bdp)
{
	int round_idx;
	char file_in[512] = { 0 }, file_ou[512] = { 0 };
	FILE* ofp = NULL;

	//Step 0 - Prep
	DP_BOUNDS = (PROB_t*)malloc(sizeof(PROB_t) * (target_round + 1));
	DIFF_TRAIL_IN_PROG = (DC_1ROUND_CHAR_t*)malloc(sizeof(DC_1ROUND_CHAR_t) * (target_round + 1));
	DIFF_TRAIL_FOR_OUT = (DC_1ROUND_CHAR_t*)malloc(sizeof(DC_1ROUND_CHAR_t) * (target_round + 1));
	memcpy(_1R_SUB_IN, _1r_sub_in, sizeof(ANA_STATE_t));
	memcpy(_nR_SUB_OUT, _nr_sub_out, sizeof(ANA_STATE_t));
	IO_printf(file_in, file_ou, _1R_SUB_IN, _nR_SUB_OUT, target_round);
	DP_BOUND_IN_PROG = given_round_bdp[target_round] - bound_gap;

	sprintf(OUT_FILE_NAME, "%s_DT_IO_R%d_[%s-%s][%0.1lf-%0.1lf].txt",
		ALG_NAME,
		target_round,
		file_in, file_ou,
		-given_round_bdp[target_round], -DP_BOUND_IN_PROG
	);
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

	memcpy(DIFF_TRAIL_IN_PROG[0].dif_o, _1R_SUB_IN, sizeof(ANA_STATE_t));
	memcpy(DIFF_TRAIL_IN_PROG[target_round].sub_o, _nR_SUB_OUT, sizeof(ANA_STATE_t));
	PERM(DIFF_TRAIL_IN_PROG[target_round].dif_o, DIFF_TRAIL_IN_PROG[target_round].sub_o);

	//Step 2(main) - Start searching
	/*Finding Trail*/
	if (target_round == 2)
	{
		SRBPN_Roundn_1_Prob_IO(ONE_PROB, 1);
	}
	else
	{
		SRBPN_Roundi_Prob_IO(ONE_PROB, target_round, 1);
	}

done:
	free(DP_BOUNDS);			DP_BOUNDS = NULL;
	free(DIFF_TRAIL_IN_PROG);	DIFF_TRAIL_IN_PROG = NULL;
	free(DIFF_TRAIL_FOR_OUT);	DIFF_TRAIL_FOR_OUT = NULL;

	return TOTAL_NUM_TRAILS;
}




void SRBPN_Roundi_Prob_IO(PROB_t BEF_DET_PROB, int target_round, int i)
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

		SRBPN_Roundi_j_Prob_IO(Xi, Yi, Pi, BEF_DET_PROB, target_round, i, next_word_idx);
	}
}


void SRBPN_Roundi_j_Prob_IO(ANA_STATE_t Xi, ANA_STATE_t Yi, PROB_t Pi[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int target_round, int i, int j)
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

			if (i == (target_round - 2))
			{
				SRBPN_Roundn_1_Prob_IO(DET_PROB, i + 1);
			}
			else if (i < (target_round - 2))
			{
				SRBPN_Roundi_Prob_IO(DET_PROB, target_round, i + 1);
			}
		}
		else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
		{
			//Move to next word
			SRBPN_Roundi_j_Prob_IO(Xi, Yi, Pi, DET_PROB, target_round, i, next_word_idx);
		}
	}
}


void SRBPN_Roundn_1_Prob_IO(PROB_t BEF_DET_PROB, int n_1)
{
	ANA_STATE_t Xn_1;
	ANA_STATE_t Yn_1;
	PROB_t		Pn_1[NUM_SBOX_IN_A_STATE_THRESHOLD];
	PROB_t		_n_1_round_expected;
	PROB_t		_n_1_bound;
	int			next_word_idx;

	memcpy(Xn_1, DIFF_TRAIL_IN_PROG[n_1 - 1].dif_o, sizeof(ANA_STATE_t));

	_n_1_bound = BEF_DET_PROB;														    //[1~(n-2)]R
	next_word_idx = expect_state_prob_at_j(&_n_1_round_expected, Xn_1, ROUND_START);
	MUL_PROB(&_n_1_bound, _n_1_bound, _n_1_round_expected);								//(n-1)R
	MUL_PROB(&_n_1_bound, _n_1_bound, DP_BOUNDS[1]);									//(n)R

	if (bound_checker_prob(_n_1_bound) != EXCEED_BOUND)
	{
		//init the sub out state and probabilities
		memset(Yn_1, 0, sizeof(ANA_STATE_t));
		init_probs(Pn_1);

		SRBPN_Roundn_1_j_Prob_IO(Xn_1, Yn_1, Pn_1, BEF_DET_PROB, n_1, next_word_idx);
	}
}


void SRBPN_Roundn_1_j_Prob_IO(ANA_STATE_t Xn_1, ANA_STATE_t Yn_1, PROB_t Pn_1[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n_1, int j)
{
	SUB_CNT_t	_n_1_j;
	SUB_WRD_t	Xn_1_j = Xn_1[j]; //It is fixed from before
	SUB_WRD_t	Yn_1_j;
	PROB_t		Pn_1_j;
	int			next_word_idx;
	PROB_t		_n_1_round_expected;
	PROB_t		_n_11_round_expected;
	PROB_t		_n_1_bound;
	PROB_t		_n_11_bound;
	PROB_t		DET_PROB;

	for (_n_1_j = 0; _n_1_j < DP_NUM_OP_FI_NONZERO[Xn_1_j]; _n_1_j++)
	{
		Yn_1_j = DP_OP_FI[Xn_1_j][_n_1_j].o;
		Pn_1_j = DP_OP_FI[Xn_1_j][_n_1_j].p;

		MUL_PROB(&DET_PROB, BEF_DET_PROB, Pn_1_j);


		next_word_idx = expect_state_prob_at_j(&_n_1_round_expected, Xn_1, j);
		MUL_PROB(&_n_1_bound, _n_1_round_expected, DET_PROB);			//[1~(n-1)]R
		MUL_PROB(&_n_1_bound, _n_1_bound, DP_BOUNDS[1]);				//(n)R

		//pruning by using descending order
		if (bound_checker_prob(_n_1_bound) == EXCEED_BOUND)
			break;

		//store current word
		Yn_1[j] = Yn_1_j;
		Pn_1[j] = Pn_1_j;


		MUL_PROB(&_n_11_bound, _n_1_round_expected, DET_PROB);				//[1~(n-1)]R

		if (expect_next_state_prob_at_n_1_j(&_n_11_round_expected, Yn_1, j) == FALSE)
			continue; //Active Truncation Check

		MUL_PROB(&_n_11_bound, _n_11_bound, _n_11_round_expected);			//(n)R

		//pruning
		if (bound_checker_prob(_n_11_bound) == EXCEED_BOUND)
			continue; //next value

		if (next_word_idx == THE_LAST_WORD)  //when this word is the last word
		{
			memcpy(DIFF_TRAIL_IN_PROG[n_1].sub_i, Xn_1, sizeof(ANA_STATE_t));
			memcpy(DIFF_TRAIL_IN_PROG[n_1].sub_o, Yn_1, sizeof(ANA_STATE_t));
			PERM(DIFF_TRAIL_IN_PROG[n_1].dif_o, DIFF_TRAIL_IN_PROG[n_1].sub_o);
			DIFF_TRAIL_IN_PROG[n_1].p = put_probs(Pn_1);

			SRBPN_Roundn_Prob_IO(DET_PROB, n_1 + 1);
		}
		else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
		{
			//Move to next word
			SRBPN_Roundn_1_j_Prob_IO(Xn_1, Yn_1, Pn_1, DET_PROB, n_1, next_word_idx);
		}
	}
}



void SRBPN_Roundn_Prob_IO(PROB_t BEF_DET_PROB, int n)
{
	ANA_STATE_t Xn;
	ANA_STATE_t Yn;
	PROB_t		Pn[NUM_SBOX_IN_A_STATE_THRESHOLD];
	PROB_t		_n_round_expected;
	PROB_t		_n_bound;
	int			next_word_idx;
	int			word_idx;

	memcpy(Xn, DIFF_TRAIL_IN_PROG[n - 1].dif_o, sizeof(ANA_STATE_t));

	//Active Truncation Check
	for (word_idx = 0; word_idx < NUM_SBOX_IN_A_STATE; word_idx++)
	{
		uint8_t ck_flag1 = Xn[word_idx] == 0 ? (uint8_t)1 : (uint8_t)0;
		uint8_t ck_flag2 = _nR_SUB_OUT[word_idx] == 0 ? (uint8_t)1 : (uint8_t)0;

		//The Active Truncation does not match
		if ((ck_flag1 ^ ck_flag2) == 1)
		{
			return;
		}
	}

	_n_bound = BEF_DET_PROB;														//[1~(n-1)]R
	next_word_idx = expect_state_prob_at_j(&_n_round_expected, Xn, ROUND_START);	//nR
	MUL_PROB(&_n_bound, _n_bound, _n_round_expected);

	if (bound_checker_prob(_n_bound) != EXCEED_BOUND)
	{
		//init the sub out state and probabilities
		memset(Yn, 0, sizeof(ANA_STATE_t));
		init_probs(Pn);

		SRBPN_Roundn_j_Prob_IO(Xn, Yn, Pn, BEF_DET_PROB, n, next_word_idx);
	}
}




void SRBPN_Roundn_j_Prob_IO(ANA_STATE_t Xn, ANA_STATE_t Yn, PROB_t Pn[NUM_SBOX_IN_A_STATE_THRESHOLD], PROB_t BEF_DET_PROB, int n, int j)
{
	SUB_WRD_t	Xn_j = Xn[j]; //It is fixed from before
	SUB_WRD_t	Yn_j = _nR_SUB_OUT[j]; //It is fixed by the out
	PROB_t		Pn_j = PROB_DDT[Xn_j][Yn_j];
	int			next_word_idx;
	PROB_t		_n_round_expected;
	PROB_t		_n_bound;
	PROB_t		DET_PROB;

	if (Pn_j == ZERO_PROB)
		return;

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

		TOTAL_NUM_TRAILS++;
		//printf("Num of Found Best Trails : %llu\r", TOTAL_NUM_TRAILS);
		Dif_trail_fprintf(OUT_FILE_NAME, "a", DIFF_TRAIL_FOR_OUT, n, DET_PROB);
	}
	else if (next_word_idx < NUM_SBOX_IN_A_STATE) //when this word is not the last word
	{
		//Move to next word
		SRBPN_Roundn_j_Prob_IO(Xn, Yn, Pn, DET_PROB, n, next_word_idx);
	}
}











