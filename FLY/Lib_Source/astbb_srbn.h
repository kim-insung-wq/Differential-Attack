#ifndef __ASTBB_SRBN_H__
#define __ASTBB_SRBN_H__
#define _CRT_SECURE_NO_WARNINGS


#ifdef __cplusplus
extern "C" {
#endif

#if defined _MSC_VER
	//Visual Studio
#ifdef _DEVELOPMENT
#define DEV_DEFINE __declspec(dllexport)
#else
#define DEV_DEFINE __declspec(dllimport)
#endif
#elif defined __GNUC__
	//GCC
#ifdef _DEVELOPMENT
#define DEV_DEFINE __attribute__ ((visibility("default")))
#else
	//nothing to define
#define DEV_DEFINE 
#endif
#endif

#if defined	__NO_INLINE__
#define DEV_INLINE //nothing
#else
#define DEV_INLINE inline
#endif


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>



/*in utils*/
enum time_flag { START_TIME = 1, END_TIME, ONLY_PRINT };
void Realtime_Print(int check_flag);

#define UNDER_BOUND			(1)
#define EXCEED_BOUND		(-1)

#define ROUND_START			(-3)
#define LAST_ROUND_START    (-2)
#define THE_LAST_WORD		(-1)

#define TRUE				1
#define FALSE				0


//0<= 2^{PROB_t} <=1 always i.e, always negative
typedef double PROB_t;
typedef double WEIGHT_t;

//#define DEBUG

//#define RESULT_ONLY
#define NUM_SBOX_IN_A_STATE_THRESHOLD	8
#define SBOX_BIT_SIZE_THRESHOLD			8
#define SBOX_CARDINALITY_THRESHOLD		(1<<SBOX_BIT_SIZE_THRESHOLD)


typedef	uint8_t				ANA_WRD_t;  //col bitsize(S-box bitsize) maximum		--> smaller or eqaul to 10(<= sizeof(ANA_WRD_t)*8)
typedef uint64_t			CIP_WRD_t;  //row bitsize(num sbox in a state) maximum	--> smaller or equal to 64(<= sizeof(CIP_WRD_t)*8)
typedef uint64_t			TRUNC_STATE_t;
typedef ANA_WRD_t			SUB_WRD_t;
typedef SUB_WRD_t			ANA_STATE_t[NUM_SBOX_IN_A_STATE_THRESHOLD];
typedef CIP_WRD_t			CIP_STATE_t[SBOX_BIT_SIZE_THRESHOLD];
typedef int					GEN_CNT_t;
typedef int					SUB_CNT_t;
typedef int					DIF_CNT_t;
typedef int					FLAG_t;
typedef uint64_t			NUM_TRAIL_t;



typedef struct 
{
	char      * ALG_NAME;

	GEN_CNT_t SBOX_BIT_SIZE;
	GEN_CNT_t SBOX_CARDINALITY;
	SUB_WRD_t * SBOX;   //SBOX_CARDINALITY

	GEN_CNT_t NUM_SBOX_IN_A_STATE;
	GEN_CNT_t * OFFSET; //NUM_SBOX_IN_A_STAT


}SRBPN_INFO_t;

typedef struct {
	uint8_t must_be_1_mask[SBOX_BIT_SIZE_THRESHOLD];
	uint8_t must_be_0_mask[SBOX_BIT_SIZE_THRESHOLD];
	uint8_t dont_care_mask[SBOX_BIT_SIZE_THRESHOLD]; // dont_care_mask ^ must_be_0_mask ^ must_be_1_mask == 0b11111111
} BitPattern64;



DEV_DEFINE void Computation_LS_of_INV_DDT();
DEV_DEFINE void Computation_LS_of_DDT();
DEV_DEFINE WEIGHT_t Cal_Key_Recovery_Comp_11R(ANA_STATE_t _1r_SUB_IN, ANA_STATE_t _nr_SUB_OUT, WEIGHT_t GIVEN_PROB, int num_of_right_pair, double c_e);
DEV_DEFINE WEIGHT_t Cal_Key_Recovery_Comp_12R(ANA_STATE_t _1r_SUB_IN, ANA_STATE_t _nr_SUB_OUT, WEIGHT_t GIVEN_PROB, int num_of_right_pair, double c_e);
DEV_DEFINE int Cal_Num_of_Involved_Bits(ANA_STATE_t _1r_SUB_IN, ANA_STATE_t _nr_SUB_OUT);
DEV_DEFINE void Prep_Dif_Trail_Searching(SRBPN_INFO_t* cipher_info);
DEV_DEFINE void Best_Trail_Prob_Only(PROB_t* target_bdp_rst, GEN_CNT_t target_round, PROB_t* prev_round_bdp);
DEV_DEFINE NUM_TRAIL_t Best_Trail_Prob_All(GEN_CNT_t target_round, PROB_t* given_round_bdp);
DEV_DEFINE NUM_TRAIL_t Best_Trail_Prob_IO(PROB_t bound_gap, ANA_STATE_t _1r_sub_in, ANA_STATE_t _nr_sub_out, GEN_CNT_t target_round, PROB_t* given_round_bdp);

DEV_DEFINE void PERM(ANA_STATE_t out, ANA_STATE_t in);
DEV_DEFINE void INV_PERM(ANA_STATE_t out, ANA_STATE_t in);


/* For Generateing Latex Tables*/
DEV_DEFINE void BP_PERM(BitPattern64* a);
DEV_DEFINE void BP_INV_PERM(BitPattern64* a);
DEV_DEFINE void PY_Cal_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_BP_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_BP_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_INV_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_INV_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_INV_BP_1r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
DEV_DEFINE void PY_Cal_INV_BP_2r_BP(ANA_STATE_t _nr_SUB_OUT, BitPattern64* _nr_BP);
#ifdef __cplusplus
}
#endif /*extern "C"*/
#endif /*__ASTBB_SRBN_H__*/



