#include "astbb_srbn.h"
#include "active_map.h"



//read-only
extern GEN_CNT_t NUM_SBOX_IN_A_STATE;

//here only
static SUB_WRD_t		INNER_1ROUND_ACTIVE_MAP[NUM_SBOX_IN_A_STATE_THRESHOLD];
static TRUNC_STATE_t	MASK = 0;

void Init_AT1(GEN_CNT_t NUM_AT1)
{
	GEN_CNT_t sbox_idx;
	memset(INNER_1ROUND_ACTIVE_MAP, 0, sizeof(SUB_WRD_t)*(size_t)NUM_SBOX_IN_A_STATE_THRESHOLD);
	for (sbox_idx = 0; sbox_idx < NUM_AT1; sbox_idx++)
	{
		INNER_1ROUND_ACTIVE_MAP[sbox_idx] = (SUB_WRD_t)1; //The First(smallest) Representive AT1 Map with the NUM_AT1 
	}

	//INIT THE MASK
	MASK = (TRUNC_STATE_t)0;
	for (sbox_idx = 0; sbox_idx < NUM_SBOX_IN_A_STATE; sbox_idx++)
	{
		MASK = MASK ^ ((TRUNC_STATE_t)1 << sbox_idx);
	}
}

DEV_INLINE TRUNC_STATE_t ACTIVE_MAP_TO_AT(void)
{
	TRUNC_STATE_t rst = 0x0;
	GEN_CNT_t word_idx;
	//compute i-round active maprst
	for (word_idx = 0; word_idx < NUM_SBOX_IN_A_STATE; word_idx++)
	{
		if (INNER_1ROUND_ACTIVE_MAP[word_idx] != 0x0)
		{
			rst |= ((TRUNC_STATE_t)0x1 << word_idx);
		}
	}
	return rst;
}

DEV_INLINE void Swap(SUB_WRD_t * x, SUB_WRD_t * y)
{
	SUB_WRD_t temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

DEV_INLINE void Sorting_Bin_Array(GEN_CNT_t from, GEN_CNT_t size)
{
	SUB_WRD_t * from_array = INNER_1ROUND_ACTIVE_MAP + from;
	GEN_CNT_t num_ones = 0;
	GEN_CNT_t i;
	for (i = 0; i < size; i++)
	{
		if (from_array[i] == 1)
		{
			num_ones++;
		}
		from_array[i] = 0;
	}
	for (i = 0; i < num_ones; i++)
	{
		from_array[i] = 1;
	}
}

DEV_INLINE FLAG_t Get_Next(void)
{
	GEN_CNT_t i;
	// Find the rightmost active word
	// which is left to non-active word(e.g, 111'1'0,) 
	// Let us call it 'first active word'
	// Look!! Here NUM_SBOX_IN_A_STATE must be > 2
	for (i = NUM_SBOX_IN_A_STATE - 2; i >= 0; --i)
		if (INNER_1ROUND_ACTIVE_MAP[i] > INNER_1ROUND_ACTIVE_MAP[i + 1])
			break;

	//when 00000011111
	if (i == -1)
	{
		return THIS_IS_THE_LAST;
	}
	else
	{
		Swap(&INNER_1ROUND_ACTIVE_MAP[i], &INNER_1ROUND_ACTIVE_MAP[i + 1]);
		Sorting_Bin_Array(i + 1, NUM_SBOX_IN_A_STATE - i - 1);
		return THIS_IS_NOT_THE_LAST;
	}
}


DEV_INLINE void AT_ROL(TRUNC_STATE_t * roted_at, int rot_offset, TRUNC_STATE_t cur_at)
{
	TRUNC_STATE_t out_at = (TRUNC_STATE_t)0;
	out_at = (cur_at << rot_offset) & MASK;
	out_at = out_at | ((cur_at >> (NUM_SBOX_IN_A_STATE - rot_offset)) & MASK);
	*roted_at = out_at;
}

DEV_INLINE FLAG_t check_repre_on_rot_sym(void)
{
	TRUNC_STATE_t cur_at = ACTIVE_MAP_TO_AT();
	int rot_offset;
	TRUNC_STATE_t roted_at;
	FLAG_t repre_rst = TRUE;
	for (rot_offset = 1; rot_offset < NUM_SBOX_IN_A_STATE; rot_offset++)
	{
		AT_ROL(&roted_at, rot_offset, cur_at);

		if (roted_at < cur_at)
		{
			repre_rst = FALSE;
			break;
		}
	}
	return repre_rst;
}


FLAG_t Next_AT1(TRUNC_STATE_t * NEXT_AT1, FLAG_t rot_sym_on)
{
	*NEXT_AT1 = ACTIVE_MAP_TO_AT();

	if (rot_sym_on == TRUE)
	{
		FLAG_t check_last	= THIS_IS_NOT_THE_LAST;
		FLAG_t check_repre	= FALSE;
		while ((check_last == THIS_IS_NOT_THE_LAST) && (check_repre == FALSE))
		{
			check_last	= Get_Next();
			check_repre = check_repre_on_rot_sym();
		}
		return check_last;
	}
	else // if (rot_sym_on == FALSE)
	{
		return Get_Next();
	}
}










