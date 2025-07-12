#ifndef  _ACTIVE_MAP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "astbb_srbn.h"


#define THIS_IS_NOT_THE_LAST  ((FLAG_t)(0))
#define THIS_IS_THE_LAST	  ((FLAG_t)(1))
void	Init_AT1(GEN_CNT_t NUM_AT1);
FLAG_t	Next_AT1(TRUNC_STATE_t * NEXT_AT1, FLAG_t rot_sym_on);

#ifdef __cplusplus
}
#endif /*extern "C"*/

#endif /*_ACTIVE_MAP_H_*/
