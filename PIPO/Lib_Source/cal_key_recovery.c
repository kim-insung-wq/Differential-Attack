#include "astbb_srbn.h"
#include "dc_prob.h"

#ifdef _MSC_VER
#include <intrin.h>
#define __builtin_popcount __popcnt
#endif

//extern(decleared in global.c)
//Write here
extern char      ALG_NAME[256];
extern GEN_CNT_t SBOX_BIT_SIZE;
extern GEN_CNT_t SBOX_CARDINALITY;
extern SUB_WRD_t SBOX[SBOX_CARDINALITY_THRESHOLD];
extern GEN_CNT_t NUM_SBOX_IN_A_STATE;
extern GEN_CNT_t OFFSET[SBOX_BIT_SIZE_THRESHOLD];
extern SUB_WRD_t SUB_WRD_MASK;
extern CIP_WRD_t CIP_WRD_MASK;
extern SUB_CNT_t INV_DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t DDT[SBOX_CARDINALITY_THRESHOLD][SBOX_CARDINALITY_THRESHOLD];
extern SUB_CNT_t Inv_LS_DDT[SBOX_CARDINALITY_THRESHOLD]; // Invariant linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 0 for all x.
extern SUB_CNT_t Com_LS_DDT[SBOX_CARDINALITY_THRESHOLD]; // Complementary linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 1 for all x.
extern SUB_CNT_t LS_DDT[SBOX_CARDINALITY_THRESHOLD]; // Invariant + Complementary e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = constant for all x.
extern SUB_CNT_t Inv_LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD]; // Invariant linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 0 for all x.
extern SUB_CNT_t Com_LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD]; // Complementary linear structure e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = 1 for all x.
extern SUB_CNT_t LS_INV_DDT[SBOX_CARDINALITY_THRESHOLD]; // Invariant + Complementary e.g. f(x) + f(x + (0,0,0,0,0,0,0,1)) = constant for all x.
extern SUB_CNT_t Standard_Vec[SBOX_BIT_SIZE_THRESHOLD];



//static PROB_t* DP_BOUNDS = NULL;
//static PROB_t	 DP_BOUND_IN_PROG;
//static char	 OUT_FILE_NAME[512];

// Assume that the cipher consists of eight 8-bit S-boxes.


typedef struct {
	WEIGHT_t value;
	unsigned int index;
} Element;


static int compare(const void* a, const void* b)
{
	WEIGHT_t diff = ((Element*)a)->value - ((Element*)b)->value;
	if (diff > 0) return 1;
	else if (diff < 0) return -1;
	else return 0;
}

static WEIGHT_t sumArray(WEIGHT_t* arr, int size) {
	WEIGHT_t sum = 0;
	int i;
	for (i = 0; i < size; i++) {
		sum += arr[i];
	}
	return sum;
}

void Computation_LS_of_DDT()
{
	SUB_CNT_t i;
	SUB_CNT_t delta_i;
	SUB_CNT_t delta_o;
	SUB_CNT_t _0, _1;
	//init
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			Inv_LS_DDT[delta_i] = 0; //Initialization to zeros...
			Com_LS_DDT[delta_i] = 0; //Initialization to zeros...
			LS_DDT[delta_i] = 0; //Initialization to zeros...
		}
	}
	Inv_LS_DDT[0] = 0b11111111;
	LS_DDT[0] = 0b11111111;
	for (i = 0; i < SBOX_BIT_SIZE; i++)
	{
		Standard_Vec[i] = (1 << i);
	}

	for (delta_i = 1; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		_0 = 0b00000000;
		_1 = 0b11111111;
		for (delta_o = 1; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			if (DDT[delta_i][delta_o] != 0)
			{
				_0 = _0 | delta_o;
				_1 = _1 & delta_o;
			}
		}
		for (i = 0; i < SBOX_BIT_SIZE; i++)
		{
			if ((Standard_Vec[i] | _0) != _0)
			{
				Inv_LS_DDT[delta_i] = Inv_LS_DDT[delta_i] ^ Standard_Vec[i]; //  0.
				LS_DDT[delta_i] = LS_DDT[delta_i] ^ Standard_Vec[i];
				
#ifdef DEBUG
				printf("\Inv_LS_DDT : (%02x, %02x)\n", delta_i, Standard_Vec[i]);
#endif // DEBUG


			}

			if ((Standard_Vec[i] & _1) != 0)
			{
				Com_LS_DDT[delta_i] = Com_LS_DDT[delta_i] ^ Standard_Vec[i]; //  1.
				LS_DDT[delta_i] = LS_DDT[delta_i] ^ Standard_Vec[i];
#ifdef DEBUG
				printf("\Com_LS_DDT : (%02x, %02x)\n", delta_i, Standard_Vec[i]);
#endif // DEBUG
			}
		}

	}

#if defined DEBUG
	int idx;
	printf("\n********Inv_LS_DDT********\n");
	for (idx = 0; idx < SBOX_CARDINALITY; idx++)
	{
		printf("%02x ", Inv_LS_DDT[idx]);
		if ((idx + 1) % 8 == 0)
			printf("\n");

	}
	printf("\n********Com_LS_DDT********\n");
	for (idx = 0; idx < SBOX_CARDINALITY; idx++)
	{
		printf("%02x ", Com_LS_DDT[idx]);
		if ((idx + 1) % 8 == 0)
			printf("\n");

	}


#endif


}

void Computation_LS_of_INV_DDT()
{
	SUB_CNT_t i;
	SUB_CNT_t delta_i;
	SUB_CNT_t delta_o;
	SUB_CNT_t _0, _1;
	SUB_CNT_t INV_Standard_Vec[SBOX_BIT_SIZE_THRESHOLD];
	//init
	for (delta_i = 0; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		for (delta_o = 0; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			Inv_LS_INV_DDT[delta_i] = 0; //Initialization to zeros...
			Com_LS_INV_DDT[delta_i] = 0; //Initialization to zeros...
			LS_INV_DDT[delta_i] = 0; //Initialization to zeros...
		}
	}
	Inv_LS_INV_DDT[0] = 0b11111111;
	LS_INV_DDT[0] = 0b11111111;
	for (i = 0; i < SBOX_BIT_SIZE; i++)
	{
		INV_Standard_Vec[i] = (1 << i);
	}

	for (delta_i = 1; delta_i < SBOX_CARDINALITY; delta_i++)
	{
		_0 = 0b00000000;
		_1 = 0b11111111;
		for (delta_o = 1; delta_o < SBOX_CARDINALITY; delta_o++)
		{
			if (INV_DDT[delta_i][delta_o] != 0)
			{
				_0 = _0 | delta_o;
				_1 = _1 & delta_o;
			}
		}
		for (i = 0; i < SBOX_BIT_SIZE; i++)
		{
			if ((INV_Standard_Vec[i] | _0) != _0)
			{
				Inv_LS_INV_DDT[delta_i] = Inv_LS_INV_DDT[delta_i] ^ INV_Standard_Vec[i]; //  0.
				LS_INV_DDT[delta_i] = LS_INV_DDT[delta_i] ^ INV_Standard_Vec[i];


			}

			if ((INV_Standard_Vec[i] & _1) != 0)
			{
				Com_LS_INV_DDT[delta_i] = Com_LS_INV_DDT[delta_i] ^ INV_Standard_Vec[i]; //  1.
				LS_INV_DDT[delta_i] = LS_INV_DDT[delta_i] ^ INV_Standard_Vec[i];

			}
		}

	}



}



static uint8_t Count_Trivial_Sbox(ANA_STATE_t DIFFERENCE)
{
	int i;
	uint8_t cnt = 0;
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		if (DIFFERENCE[i] == 0)
		{
			cnt++;
		}
	}
	return cnt;
}

static int* Obt_index_NonTrvial_Sbox(ANA_STATE_t arr, GEN_CNT_t size) 
{

	int nonTrivialCount = 0;
	int i;
	for (i = 0; i < size; i++)
	{
		if (arr[i] != 0)
		{
			nonTrivialCount++;
		}
	}

	int* nonTrivial_Indices = (int*)malloc(nonTrivialCount * sizeof(int));
	if (nonTrivial_Indices == NULL)
	{
		printf("Memory allocation failed.\n");
		return -1;
	}

	int j = 0;
	for (i = 0; i < size; i++)
	{
		if (arr[i] != 0)
		{
			if (j < nonTrivialCount)
			{
				nonTrivial_Indices[j] = i;
				j++;
			}
			else
			{
				printf("Unexpected overrun error.\n");
				free(nonTrivial_Indices);
				return -1;
			}
		}
	}

	return nonTrivial_Indices;
}

static int* Obt_index_NonTrvial_Sbox_from_BP(BitPattern64* _nr_BP, GEN_CNT_t size)
{
	int nonTrivialCount = 0;
	int i;
	for (i = 0; i < size; i++)
	{
		if (_nr_BP->must_be_0_mask[i] != 0b11111111)
		{
			nonTrivialCount++;
		}
	}

	int* nonTrivial_Indices_from_BP = (int*)malloc(nonTrivialCount * sizeof(int));
	if (nonTrivialCount > 0)
	{
		if (nonTrivial_Indices_from_BP == NULL)
		{
			printf("Memory allocation failed.\n");
			return -1;
		}
		
		int j = 0;
		for (i = 0; i < size; i++)
		{
			if (_nr_BP->must_be_0_mask[i] != 0b11111111)
			{
				if (j < nonTrivialCount)
				{
					nonTrivial_Indices_from_BP[j] = i;
					j++;
				}
				else
				{
					printf("Unexpected overrun error.\n");
					free(nonTrivial_Indices_from_BP);
					return -1;
				}
			}
		}
	}

	return nonTrivial_Indices_from_BP;

}




static void Cal_BP_1r(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
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



int Cal_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP)
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
	//printf("num_of_inactive_sboxes : %d\n", num_of_inactive_sboxes);
	return num_of_inactive_sboxes;
}


int is_in(int arr[], int size, int target) {
	if (arr == NULL)
	{
		return 0;
	}
	for (int i = 0; i < size; i++) {
		if (arr[i] == target) {
			return 1;
		}
	}
	return 0;
}

int Cal_Num_of_Involved_Bits(ANA_STATE_t _1r_SUB_IN, ANA_STATE_t _nr_SUB_OUT)
{
	ANA_STATE_t DIFF_IN = { 0, };
	ANA_STATE_t DIFF_OUT = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	INV_PERM(DIFF_IN, _1r_SUB_IN);
	int num_of_fwd_guessed_sboxes_1r = NUM_SBOX_IN_A_STATE - Count_Trivial_Sbox(DIFF_IN);
	int num_of_bwd_guessed_sboxes_1r = NUM_SBOX_IN_A_STATE - Count_Trivial_Sbox(DIFF_OUT);
	int num_of_bwd_guessed_sboxes_2r = NUM_SBOX_IN_A_STATE;
	int num_of_involved_bits = SBOX_BIT_SIZE * num_of_fwd_guessed_sboxes_1r + SBOX_BIT_SIZE * num_of_bwd_guessed_sboxes_2r + (SBOX_BIT_SIZE - num_of_fwd_guessed_sboxes_1r) * num_of_bwd_guessed_sboxes_1r;
#if defined DEBUG

	printf("\nCount_Trivial_Sbox(DIFF_IN) : %d\n", Count_Trivial_Sbox(DIFF_IN));
	printf("\NUM_SBOX_IN_A_STATE : %d\n", NUM_SBOX_IN_A_STATE);
	printf("\nnum_of_fwd_guessed_sboxes_1r : %d\n", num_of_fwd_guessed_sboxes_1r);
	printf("\nnum_of_bwd_guessed_sboxes_1r : %d\n", num_of_bwd_guessed_sboxes_1r);
	printf("\nnum_of_bwd_guessed_sboxes_2r : %d\n", num_of_bwd_guessed_sboxes_2r);
	printf("\nnum_of_involved_bits : %d\n", num_of_involved_bits);
#endif	
	return num_of_involved_bits;
}

WEIGHT_t Cal_Key_Recovery_Comp(ANA_STATE_t _1r_SUB_IN, ANA_STATE_t _nr_SUB_OUT, WEIGHT_t GIVEN_PROB, int num_of_right_pair, double c_e)
{
	int i;
	BitPattern64 _2r_BP = { 0 };
	BitPattern64 _1r_BP = { 0 };
	ANA_STATE_t DIFF_OUT = { 0, };
	ANA_STATE_t DIFF_IN = { 0, };
	ANA_STATE_t Determined_bit_2r = { 0, };
	ANA_STATE_t Determined_bit_1r = { 0, };
	PERM(DIFF_OUT, _nr_SUB_OUT);
	INV_PERM(DIFF_IN, _1r_SUB_IN);
	int num_of_fwd_guessed_sboxes_1r = NUM_SBOX_IN_A_STATE - Count_Trivial_Sbox(DIFF_IN);
	WEIGHT_t d_in = SBOX_BIT_SIZE * num_of_fwd_guessed_sboxes_1r;
	int num_of_bwd_guessed_sboxes_1r = NUM_SBOX_IN_A_STATE - Count_Trivial_Sbox(DIFF_OUT);
	int num_of_bwd_guessed_sboxes_2r = NUM_SBOX_IN_A_STATE;
	int num_of_bwd_nontrivial_sboxes_2r = NUM_SBOX_IN_A_STATE - Cal_2r_BP(_nr_SUB_OUT, &_2r_BP);
	double ciphertext_filtering_ratio = SBOX_BIT_SIZE * (NUM_SBOX_IN_A_STATE - num_of_bwd_nontrivial_sboxes_2r);
	WEIGHT_t tot_num_of_guessed_sboxes = num_of_fwd_guessed_sboxes_1r + num_of_bwd_guessed_sboxes_1r + num_of_bwd_guessed_sboxes_2r;
	WEIGHT_t guessing_1r_fwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t guessing_1r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t guessing_2r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t filtering_1r_fwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t filtering_1r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t filtering_2r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	WEIGHT_t w_tot_time_comp;
	double time_comp_inner;
	WEIGHT_t w_time_comp_inner;
	Cal_BP_1r(_nr_SUB_OUT, &_1r_BP);
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		Determined_bit_1r[i] = LS_DDT[DIFF_OUT[i]];
	}
	PERM(Determined_bit_2r, Determined_bit_1r);

#if defined DEBUG
	printf("\nDIFF_IN : ");
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		printf("%02X ", DIFF_IN[i]);
	}
	printf("\nDIFF_OUT : ");
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		printf("%02X ", DIFF_OUT[i]);
	}
	printf("\nDetermined_bit_1r : ");
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		printf("%02X ", Determined_bit_1r[i]);
	}

	printf("\nDetermined_bit_2r : ");
	for (i = 0; i < NUM_SBOX_IN_A_STATE; i++)
	{
		printf("%02X ", Determined_bit_2r[i]);
	}
#endif

	int* nonTrivial_Indices_1r_bwd = Obt_index_NonTrvial_Sbox(DIFF_OUT, NUM_SBOX_IN_A_STATE);

	


	for (i = 0; i < num_of_fwd_guessed_sboxes_1r; i++)
	{
		guessing_1r_fwd[i] = SBOX_BIT_SIZE;
		filtering_1r_fwd[i] = SBOX_BIT_SIZE;
	}


	for (i = 0; i < num_of_bwd_guessed_sboxes_1r; i++)
	{
		guessing_1r_bwd[i] = SBOX_BIT_SIZE - num_of_fwd_guessed_sboxes_1r;
		filtering_1r_bwd[i] = SBOX_BIT_SIZE - __builtin_popcount(Determined_bit_1r[nonTrivial_Indices_1r_bwd[i]]);
	}
	free(nonTrivial_Indices_1r_bwd);

	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		guessing_2r_bwd[i] = SBOX_BIT_SIZE;
	}

	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		if (_1r_BP.must_be_0_mask[i] == 0b11111111)
		{
			filtering_2r_bwd[i] = 0;
		}
		
		else 
		{
			filtering_2r_bwd[i] = __builtin_popcount(Determined_bit_2r[i]);
		}
		
	}

#if defined DEBUG
	printf("\guessing_1r_fwd : ");
	for (i = 0; i < num_of_fwd_guessed_sboxes_1r; i++)
	{
		printf("%f, ", guessing_1r_fwd[i]);
	}
	printf("\n");

	printf("\filtering_1r_fwd : ");
	for (i = 0; i < num_of_fwd_guessed_sboxes_1r; i++)
	{
		printf("%f, ", filtering_1r_fwd[i]);
	}
	printf("\n");

	printf("\guessing_1r_bwd : ");
	for (i = 0; i < num_of_bwd_guessed_sboxes_1r; i++)
	{
		printf("%f, ", guessing_1r_bwd[i]);
	}
	printf("\n");

	printf("\filtering_1r_bwd : ");
	for (i = 0; i < num_of_bwd_guessed_sboxes_1r; i++)
	{
		printf("%f, ", filtering_1r_bwd[i]);
	}
	printf("\n");

	printf("\guessing_2r_bwd : ");
	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		printf("%f, ", guessing_2r_bwd[i]);
	}
	printf("\n");

	printf("\filtering_2r_bwd : ");
	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		printf("%f, ", filtering_2r_bwd[i]);
	}
	printf("\n");

	printf("\n# of Determined_bit_2r : ");
	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		printf("%d, ", __builtin_popcount(Determined_bit_2r[i]));
	}
	printf("\n");
#endif



	Element comp_arry_1r_fwd[NUM_SBOX_IN_A_STATE_THRESHOLD];
	Element comp_arry_1r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD];
	Element comp_arry_2r_bwd[NUM_SBOX_IN_A_STATE_THRESHOLD];

	for (i = 0; i < NUM_SBOX_IN_A_STATE_THRESHOLD; i++)
	{
		if (i < num_of_fwd_guessed_sboxes_1r)
		{
			comp_arry_1r_fwd[i].value = guessing_1r_fwd[i] - filtering_1r_fwd[i];
		}
		else
		{
			comp_arry_1r_fwd[i].value = 100000;
		}
		comp_arry_1r_fwd[i].index = i;
	}

	for (i = 0; i < NUM_SBOX_IN_A_STATE_THRESHOLD; i++)
	{
		if (i < num_of_bwd_guessed_sboxes_1r)
		{
			comp_arry_1r_bwd[i].value = guessing_1r_bwd[i] - filtering_1r_bwd[i];
		}
		else
		{
			comp_arry_1r_bwd[i].value = 100000;
		}
		comp_arry_1r_bwd[i].index = i;
	}
	for (i = 0; i < NUM_SBOX_IN_A_STATE_THRESHOLD; i++)
	{
		if (i < num_of_bwd_guessed_sboxes_2r)
		{
			comp_arry_2r_bwd[i].value = guessing_2r_bwd[i] - filtering_2r_bwd[i];
		}
		else
		{
			comp_arry_2r_bwd[i].value = 100000;
		}
		comp_arry_2r_bwd[i].index = i;
	}
	qsort(comp_arry_1r_fwd, NUM_SBOX_IN_A_STATE_THRESHOLD, sizeof(Element), compare);
	qsort(comp_arry_1r_bwd, NUM_SBOX_IN_A_STATE_THRESHOLD, sizeof(Element), compare);
	qsort(comp_arry_2r_bwd, NUM_SBOX_IN_A_STATE_THRESHOLD, sizeof(Element), compare);

	WEIGHT_t sequence[3 * NUM_SBOX_IN_A_STATE_THRESHOLD] = { 0, };
	
	for (i = 0; i < num_of_fwd_guessed_sboxes_1r; i++)
	{
		sequence[i] = comp_arry_1r_fwd[i].value;
	}
	
	for (i = 0; i < num_of_bwd_guessed_sboxes_2r; i++)
	{
		sequence[num_of_fwd_guessed_sboxes_1r + i] = comp_arry_2r_bwd[i].value;
	}

	for (i = 0; i < num_of_bwd_guessed_sboxes_1r; i++)
	{
		sequence[num_of_fwd_guessed_sboxes_1r + num_of_bwd_guessed_sboxes_2r + i] = comp_arry_1r_bwd[i].value;
	}

#if defined DEBUG
	printf("\nsequence : ");
	for (i = 0; i < tot_num_of_guessed_sboxes; i++)
	{
		printf("%f ", sequence[i]);
	}
	printf("\n");
#endif



	time_comp_inner = 0;
	
	for (i = 0; i < tot_num_of_guessed_sboxes; i++)
	{
		time_comp_inner = time_comp_inner + pow(2, sumArray(sequence, i+1));

	}

	time_comp_inner = 1 + pow(2.0, -c_e) * time_comp_inner;

	w_time_comp_inner = log2(time_comp_inner);



	w_tot_time_comp = log2(num_of_right_pair) + 1 - GIVEN_PROB - ciphertext_filtering_ratio + d_in + w_time_comp_inner;
#if defined DEBUG
	printf("\nw_time_comp_inner : %.20f", w_time_comp_inner);
	printf("\nw_time_comp_inner : %.20f\n", w_time_comp_inner);
	printf("log2(RIGHT)+1              : %.20f\n", log2(num_of_right_pair) + 1);
	printf("GIVEN_PROB                 : %.20f\n", GIVEN_PROB);
	printf("ciphertext_filtering_ratio :   %.20f\n", ciphertext_filtering_ratio);
	printf("structure_size             : %.20f\n", d_in);
	printf("w_time_comp_inner          : %.20f\n", w_time_comp_inner);
	printf("\nw_tot_time_comp : %f\n", w_tot_time_comp);
#endif

	return w_tot_time_comp;
}
