#define YES 1
#define NO 0

#define LATTICEEASY 1
#define HLATTICE1 2
#define HLATTICE2 3

#define MINKOWSKI_BACKGROUND 1
#define FRW_BACKGROUND 2
#define FRW_PERTURB_FIXED 3
#define FRW_PERTURB_ADAPTIVE 4

#define SYMPLECTIC_2ND 1
#define SYMPLECTIC_4TH 2
#define SYMPLECTIC_6TH 3

#define UNCHANGED 1
#define OUTPUT_TT 2
#define OUTPUT_ADAPTIVE 3


#define LOOP  do k=1,n;do j=1,n;do i=1,n
#define ENDLOOP enddo; enddo; enddo
#define LOOPFLD LOOP; do fld=1,ns
#define ENDLOOPFLD ENDLOOP; enddo 
#define FLDLOOP do fld=1,ns; LOOP
#define ENDFLDLOOP enddo; ENDLOOP 
#define DEFINE_IND integer(IB):: i, j, k, fld
#define DEFINE_ALLIND integer(IB):: i,j,k,i1,j1,k1,i2,j2,k2, fld

!!since V2.0 this feature is off
#define MATCH_CONFIGURATION_SPACE_FLUC NO



!! Slice cuts
#define YZ 1
#define XY 2
#define XZ 3