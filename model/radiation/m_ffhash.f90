module m_ffhash
  implicit none
#define FFH_KEY_TYPE character(len=40)
#define FFH_KEY_IS_STRING
#define FFH_VAL_TYPE real*8
#include "./ffhash/ffhash_inc.f90"
end module m_ffhash
