#if WANTGW || METRIC_OPTION == FRW_PERTURB_FIXED || METRIC_OPTION == FRW_PERTURB_ADAPTIVE
#define HIJ_DEFINED YES
#else
#define HIJ_DEFINED NO
#endif

#if METRIC_OPTION == FRW_PERTURB_FIXED || METRIC_OPTION == FRW_PERTURB_ADAPTIVE
#define METRIC_PERTURB YES
#else
#define METRIC_PERTURB NO
#endif

!!defining discrete derivatives on the grids
#define GRID(f,i,j,k) f(sind(i),sind(j),sind(k))
#define GRID_FLD(f,fld,i,j,k) f(fld,sind(i),sind(j),sind(k))

#if DIS_SCHEME == LATTICEEASY

#define LAPLACIAN(i,j,k) (GRID_FLD(fields_f,:,i+1,j,k)+GRID_FLD(fields_f,:,i-1,j,k)+GRID_FLD(fields_f,:,i,j+1,k)+GRID_FLD(fields_f,:,i,j-1,k)+GRID_FLD(fields_f,:,i,j,k+1)+GRID_FLD(fields_f,:,i,j,k-1) - 6.*GRID_FLD(fields_f,:,i,j,k))

#elif DIS_SCHEME == HLATTICE1

#define TWO_DF_X(f,i,j,k) (GRID(f,i+1,j,k)-GRID(f,i-1,j,k))
#define TWO_DF_Y(f,i,j,k) (GRID(f,i,j+1,k)-GRID(f,i,j-1,k))
#define TWO_DF_Z(f,i,j,k) (GRID(f,i,j,k+1)-GRID(f,i,j,k-1))

#define FOUR_D2F_XX(f,i,j,k) (GRID(f,i+2,j,k)+GRID(f,i-2,j,k)-2.*GRID(f,i,j,k))
#define FOUR_D2F_YY(f,i,j,k) (GRID(f,i,j+2,k)+GRID(f,i,j-2,k)-2.*GRID(f,i,j,k))
#define FOUR_D2F_ZZ(f,i,j,k) (GRID(f,i,j,k+2)+GRID(f,i,j,k-2)-2.*GRID(f,i,j,k))
#define FOUR_D2F_XY(f,i,j,k) (GRID(f,i+1,j+1,k)+GRID(f,i-1,j-1,k)-GRID(f,i+1,j-1,k)-GRID(f,i-1,j+1,k))
#define FOUR_D2F_XZ(f,i,j,k) (GRID(f,i+1,j,k+1)+GRID(f,i-1,j,k-1)-GRID(f,i+1,j,k-1)-GRID(f,i-1,j,k+1))
#define FOUR_D2F_YZ(f,i,j,k) (GRID(f,i,j+1,k+1)+GRID(f,i,j-1,k-1)-GRID(f,i,j+1,k-1)-GRID(f,i,j-1,k+1))

#define FOUR_NABLA2(f,i,j,k) (GRID(f,i+2,j,k)+GRID(f,i-2,j,k)+GRID(f,i,j+2,k)+GRID(f,i,j-2,k)+GRID(f,i,j,k+2)+GRID(f,i,j,k-2)-6.*GRID(f,i,j,k))

#define TWO_DF_X_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+1,j,k)-GRID_FLD(f,fld,i-1,j,k))
#define TWO_DF_Y_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j+1,k)-GRID_FLD(f,fld,i,j-1,k))
#define TWO_DF_Z_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j,k+1)-GRID_FLD(f,fld,i,j,k-1))

#define FOUR_D2F_XX_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+2,j,k)+GRID_FLD(f,fld,i-2,j,k)-2.*GRID_FLD(f,fld,i,j,k))
#define FOUR_D2F_YY_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j+2,k)+GRID_FLD(f,fld,i,j-2,k)-2.*GRID_FLD(f,fld,i,j,k))
#define FOUR_D2F_ZZ_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j,k+2)+GRID_FLD(f,fld,i,j,k-2)-2.*GRID_FLD(f,fld,i,j,k))
#define FOUR_D2F_XY_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+1,j+1,k)+GRID_FLD(f,fld,i-1,j-1,k)-GRID_FLD(f,fld,i+1,j-1,k)-GRID_FLD(f,fld,i-1,j+1,k))
#define FOUR_D2F_XZ_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+1,j,k+1)+GRID_FLD(f,fld,i-1,j,k-1)-GRID_FLD(f,fld,i+1,j,k-1)-GRID_FLD(f,fld,i-1,j,k+1))
#define FOUR_D2F_YZ_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j+1,k+1)+GRID_FLD(f,fld,i,j-1,k-1)-GRID_FLD(f,fld,i,j+1,k-1)-GRID_FLD(f,fld,i,j-1,k+1))

#define FOUR_NABLA2_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+2,j,k)+GRID_FLD(f,fld,i-2,j,k)+GRID_FLD(f,fld,i,j+2,k)+GRID_FLD(f,fld,i,j-2,k)+GRID_FLD(f,fld,i,j,k+2)+GRID_FLD(f,fld,i,j,k-2)-6.*GRID_FLD(f,fld,i,j,k))

#define TWO_DFLD_X(i,j,k) TWO_DF_X_FLD(fields_f,:,i,j,k)
#define TWO_DFLD_Y(i,j,k) TWO_DF_Y_FLD(fields_f,:,i,j,k)
#define TWO_DFLD_Z(i,j,k) TWO_DF_Z_FLD(fields_f,:,i,j,k)

#define FOUR_LAP(i,j,k)  FOUR_NABLA2_FLD(fields_f,:,i,j,k)

#define FOUR_D2FLD_XY(i,j,k) FOUR_D2F_XY_FLD(fields_f,:,i,j,k)
#define FOUR_D2FLD_XZ(i,j,k) FOUR_D2F_XZ_FLD(fields_f,:,i,j,k)
#define FOUR_D2FLD_YZ(i,j,k) FOUR_D2F_YZ_FLD(fields_f,:,i,j,k)

#define TWO_DU_X(fld,i,j,k) TWO_DF_X_FLD(metric_h,fld,i,j,k)
#define TWO_DU_Y(fld,i,j,k) TWO_DF_Y_FLD(metric_h,fld,i,j,k)
#define TWO_DU_Z(fld,i,j,k) TWO_DF_Z_FLD(metric_h,fld,i,j,k)

#define TWO_DV_X(fld,i,j,k) TWO_DF_X_FLD(metric_h,fld+3,i,j,k)
#define TWO_DV_Y(fld,i,j,k) TWO_DF_Y_FLD(metric_h,fld+3,i,j,k)
#define TWO_DV_Z(fld,i,j,k) TWO_DF_Z_FLD(metric_h,fld+3,i,j,k)

#define FOUR_D2U_XX(fld,i,j,k) FOUR_D2F_XX_FLD(metric_h,fld,i,j,k)
#define FOUR_D2U_YY(fld,i,j,k) FOUR_D2F_YY_FLD(metric_h,fld,i,j,k)
#define FOUR_D2U_ZZ(fld,i,j,k) FOUR_D2F_ZZ_FLD(metric_h,fld,i,j,k)
#define FOUR_D2U_XY(fld,i,j,k) FOUR_D2F_XY_FLD(metric_h,fld,i,j,k)
#define FOUR_D2U_XZ(fld,i,j,k) FOUR_D2F_XZ_FLD(metric_h,fld,i,j,k)
#define FOUR_D2U_YZ(fld,i,j,k) FOUR_D2F_YZ_FLD(metric_h,fld,i,j,k)

#define FOUR_D2V_XX(fld,i,j,k) FOUR_D2F_XX_FLD(metric_h,fld+3,i,j,k)
#define FOUR_D2V_YY(fld,i,j,k) FOUR_D2F_YY_FLD(metric_h,fld+3,i,j,k)
#define FOUR_D2V_ZZ(fld,i,j,k) FOUR_D2F_ZZ_FLD(metric_h,fld+3,i,j,k)
#define FOUR_D2V_XY(fld,i,j,k) FOUR_D2F_XY_FLD(metric_h,fld+3,i,j,k)
#define FOUR_D2V_XZ(fld,i,j,k) FOUR_D2F_XZ_FLD(metric_h,fld+3,i,j,k)
#define FOUR_D2V_YZ(fld,i,j,k) FOUR_D2F_YZ_FLD(metric_h,fld+3,i,j,k)

#elif DIS_SCHEME == HLATTICE2

#define TWELVE_DF_X(f,i,j,k) (8.*(GRID(f,i+1,j,k)-GRID(f,i-1,j,k))-(GRID(f,i+2,j,k)-GRID(f,i-2,j,k)))
#define TWELVE_DF_Y(f,i,j,k) (8.*(GRID(f,i,j+1,k)-GRID(f,i,j-1,k))-(GRID(f,i,j+2,k)-GRID(f,i,j-2,k)))
#define TWELVE_DF_Z(f,i,j,k) (8.*(GRID(f,i,j,k+1)-GRID(f,i,j,k-1))-(GRID(f,i,j,k+2)-GRID(f,i,j,k-2)))

#define I44_D2F_XX(f,i,j,k) (GRID(f,i+4,j,k)+GRID(f,i-4,j,k)+16.*(GRID(f,i+1,j,k)+GRID(f,i-1,j,k)-GRID(f,i+3,j,k)-GRID(f,i-3,j,k)+4.*(GRID(f,i+2,j,k)+GRID(f,i-2,j,k)))-130.*GRID(f,i,j,k))
#define I44_D2F_YY(f,i,j,k) (GRID(f,i,j+4,k)+GRID(f,i,j-4,k)+16.*(GRID(f,i,j+1,k)+GRID(f,i,j-1,k)-GRID(f,i,j+3,k)-GRID(f,i,j-3,k)+4.*(GRID(f,i,j+2,k)+GRID(f,i,j-2,k)))-130.*GRID(f,i,j,k))
#define I44_D2F_ZZ(f,i,j,k) (GRID(f,i,j,k+4)+GRID(f,i,j,k-4)+16.*(GRID(f,i,j,k+1)+GRID(f,i,j,k-1)-GRID(f,i,j,k+3)-GRID(f,i,j,k-3)+4.*(GRID(f,i,j,k+2)+GRID(f,i,j,k-2)))-130.*GRID(f,i,j,k))
#define I44_D2F_XY(f,i,j,k) ((GRID(f,i+2,j+2,k)+GRID(f,i-2,j-2,k)-GRID(f,i+2,j-2,k)-GRID(f,i-2,j+2,k))+8.*(8.*(GRID(f,i+1,j+1,k)+GRID(f,i-1,j-1,k)-GRID(f,i+1,j-1,k)-GRID(f,i-1,j+1,k))-(GRID(f,i+2,j+1,k)+GRID(f,i-2,j-1,k)-GRID(f,i+2,j-1,k)-GRID(f,i-2,j+1,k))-(GRID(f,i+1,j+2,k)+GRID(f,i-1,j-2,k)-GRID(f,i-1,j+2,k)-GRID(f,i+1,j-2,k))))
#define I44_D2F_XZ(f,i,j,k) ((GRID(f,i+2,j,k+2)+GRID(f,i-2,j,k-2)-GRID(f,i+2,j,k-2)-GRID(f,i-2,j,k+2))+8.*(8.*(GRID(f,i+1,j,k+1)+GRID(f,i-1,j,k-1)-GRID(f,i+1,j,k-1)-GRID(f,i-1,j,k+1))-(GRID(f,i+2,j,k+1)+GRID(f,i-2,j,k-1)-GRID(f,i+2,j,k-1)-GRID(f,i-2,j,k+1))-(GRID(f,i+1,j,k+2)+GRID(f,i-1,j,k-2)-GRID(f,i+1,j,k-2)-GRID(f,i-1,j,k+2))))
#define I44_D2F_YZ(f,i,j,k) ((GRID(f,i,j+2,k+2)+GRID(f,i,j-2,k-2)-GRID(f,i,j+2,k-2)-GRID(f,i,j-2,k+2))+8.*(8.*(GRID(f,i,j+1,k+1)+GRID(f,i,j-1,k-1)-GRID(f,i,j+1,k-1)-GRID(f,i,j-1,k+1))-(GRID(f,i,j+2,k+1)+GRID(f,i,j-2,k-1)-GRID(f,i,j+2,k-1)-GRID(f,i,j-2,k+1))-(GRID(f,i,j+1,k+2)+GRID(f,i,j-1,k-2)-GRID(f,i,j+1,k-2)-GRID(f,i,j-1,k+2))))

#define I44_NABLA2(f,i,j,k) (GRID(f,i+4,j,k)+GRID(f,i-4,j,k)+GRID(f,i,j+4,k)+GRID(f,i,j-4,k)+GRID(f,i,j,k+4)+GRID(f,i,j,k-4)+16.*(GRID(f,i+1,j,k)+GRID(f,i-1,j,k)-GRID(f,i+3,j,k)-GRID(f,i-3,j,k)+GRID(f,i,j+1,k)+GRID(f,i,j-1,k)-GRID(f,i,j+3,k)-GRID(f,i,j-3,k)+GRID(f,i,j,k+1)+GRID(f,i,j,k-1)-GRID(f,i,j,k+3)-GRID(f,i,j,k-3)+4.*(GRID(f,i+2,j,k)+GRID(f,i-2,j,k)+GRID(f,i,j+2,k)+GRID(f,i,j-2,k)+GRID(f,i,j,k+2)+GRID(f,i,j,k-2)))-390.*GRID(f,i,j,k))

#define TWELVE_DF_X_FLD(f,fld,i,j,k) (8.*(GRID_FLD(f,fld,i+1,j,k)-GRID_FLD(f,fld,i-1,j,k))-(GRID_FLD(f,fld,i+2,j,k)-GRID_FLD(f,fld,i-2,j,k)))
#define TWELVE_DF_Y_FLD(f,fld,i,j,k) (8.*(GRID_FLD(f,fld,i,j+1,k)-GRID_FLD(f,fld,i,j-1,k))-(GRID_FLD(f,fld,i,j+2,k)-GRID_FLD(f,fld,i,j-2,k)))
#define TWELVE_DF_Z_FLD(f,fld,i,j,k) (8.*(GRID_FLD(f,fld,i,j,k+1)-GRID_FLD(f,fld,i,j,k-1))-(GRID_FLD(f,fld,i,j,k+2)-GRID_FLD(f,fld,i,j,k-2)))

#define I44_D2F_XX_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+4,j,k)+GRID_FLD(f,fld,i-4,j,k)+16.*(GRID_FLD(f,fld,i+1,j,k)+GRID_FLD(f,fld,i-1,j,k)-GRID_FLD(f,fld,i+3,j,k)-GRID_FLD(f,fld,i-3,j,k)+4.*(GRID_FLD(f,fld,i+2,j,k)+GRID_FLD(f,fld,i-2,j,k)))-130.*GRID_FLD(f,fld,i,j,k))
#define I44_D2F_YY_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j+4,k)+GRID_FLD(f,fld,i,j-4,k)+16.*(GRID_FLD(f,fld,i,j+1,k)+GRID_FLD(f,fld,i,j-1,k)-GRID_FLD(f,fld,i,j+3,k)-GRID_FLD(f,fld,i,j-3,k)+4.*(GRID_FLD(f,fld,i,j+2,k)+GRID_FLD(f,fld,i,j-2,k)))-130.*GRID_FLD(f,fld,i,j,k))
#define I44_D2F_ZZ_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i,j,k+4)+GRID_FLD(f,fld,i,j,k-4)+16.*(GRID_FLD(f,fld,i,j,k+1)+GRID_FLD(f,fld,i,j,k-1)-GRID_FLD(f,fld,i,j,k+3)-GRID_FLD(f,fld,i,j,k-3)+4.*(GRID_FLD(f,fld,i,j,k+2)+GRID_FLD(f,fld,i,j,k-2)))-130.*GRID_FLD(f,fld,i,j,k))
#define I44_D2F_XY_FLD(f,fld,i,j,k) ((GRID_FLD(f,fld,i+2,j+2,k)+GRID_FLD(f,fld,i-2,j-2,k)-GRID_FLD(f,fld,i+2,j-2,k)-GRID_FLD(f,fld,i-2,j+2,k))+8.*(8.*(GRID_FLD(f,fld,i+1,j+1,k)+GRID_FLD(f,fld,i-1,j-1,k)-GRID_FLD(f,fld,i+1,j-1,k)-GRID_FLD(f,fld,i-1,j+1,k))-(GRID_FLD(f,fld,i+2,j+1,k)+GRID_FLD(f,fld,i-2,j-1,k)-GRID_FLD(f,fld,i+2,j-1,k)-GRID_FLD(f,fld,i-2,j+1,k))-(GRID_FLD(f,fld,i+1,j+2,k)+GRID_FLD(f,fld,i-1,j-2,k)-GRID_FLD(f,fld,i-1,j+2,k)-GRID_FLD(f,fld,i+1,j-2,k))))
#define I44_D2F_XZ_FLD(f,fld,i,j,k) ((GRID_FLD(f,fld,i+2,j,k+2)+GRID_FLD(f,fld,i-2,j,k-2)-GRID_FLD(f,fld,i+2,j,k-2)-GRID_FLD(f,fld,i-2,j,k+2))+8.*(8.*(GRID_FLD(f,fld,i+1,j,k+1)+GRID_FLD(f,fld,i-1,j,k-1)-GRID_FLD(f,fld,i+1,j,k-1)-GRID_FLD(f,fld,i-1,j,k+1))-(GRID_FLD(f,fld,i+2,j,k+1)+GRID_FLD(f,fld,i-2,j,k-1)-GRID_FLD(f,fld,i+2,j,k-1)-GRID_FLD(f,fld,i-2,j,k+1))-(GRID_FLD(f,fld,i+1,j,k+2)+GRID_FLD(f,fld,i-1,j,k-2)-GRID_FLD(f,fld,i+1,j,k-2)-GRID_FLD(f,fld,i-1,j,k+2))))
#define I44_D2F_YZ_FLD(f,fld,i,j,k) ((GRID_FLD(f,fld,i,j+2,k+2)+GRID_FLD(f,fld,i,j-2,k-2)-GRID_FLD(f,fld,i,j+2,k-2)-GRID_FLD(f,fld,i,j-2,k+2))+8.*(8.*(GRID_FLD(f,fld,i,j+1,k+1)+GRID_FLD(f,fld,i,j-1,k-1)-GRID_FLD(f,fld,i,j+1,k-1)-GRID_FLD(f,fld,i,j-1,k+1))-(GRID_FLD(f,fld,i,j+2,k+1)+GRID_FLD(f,fld,i,j-2,k-1)-GRID_FLD(f,fld,i,j+2,k-1)-GRID_FLD(f,fld,i,j-2,k+1))-(GRID_FLD(f,fld,i,j+1,k+2)+GRID_FLD(f,fld,i,j-1,k-2)-GRID_FLD(f,fld,i,j+1,k-2)-GRID_FLD(f,fld,i,j-1,k+2))))

#define I44_NABLA2_FLD(f,fld,i,j,k) (GRID_FLD(f,fld,i+4,j,k)+GRID_FLD(f,fld,i-4,j,k)+GRID_FLD(f,fld,i,j+4,k)+GRID_FLD(f,fld,i,j-4,k)+GRID_FLD(f,fld,i,j,k+4)+GRID_FLD(f,fld,i,j,k-4)+16.*(GRID_FLD(f,fld,i+1,j,k)+GRID_FLD(f,fld,i-1,j,k)-GRID_FLD(f,fld,i+3,j,k)-GRID_FLD(f,fld,i-3,j,k)+GRID_FLD(f,fld,i,j+1,k)+GRID_FLD(f,fld,i,j-1,k)-GRID_FLD(f,fld,i,j+3,k)-GRID_FLD(f,fld,i,j-3,k)+GRID_FLD(f,fld,i,j,k+1)+GRID_FLD(f,fld,i,j,k-1)-GRID_FLD(f,fld,i,j,k+3)-GRID_FLD(f,fld,i,j,k-3)+4.*(GRID_FLD(f,fld,i+2,j,k)+GRID_FLD(f,fld,i-2,j,k)+GRID_FLD(f,fld,i,j+2,k)+GRID_FLD(f,fld,i,j-2,k)+GRID_FLD(f,fld,i,j,k+2)+GRID_FLD(f,fld,i,j,k-2)))-390.*GRID_FLD(f,fld,i,j,k))

#define TWELVE_DFLD_X(i,j,k) TWELVE_DF_X_FLD(fields_f,:,i,j,k)
#define TWELVE_DFLD_Y(i,j,k) TWELVE_DF_Y_FLD(fields_f,:,i,j,k)
#define TWELVE_DFLD_Z(i,j,k) TWELVE_DF_Z_FLD(fields_f,:,i,j,k)

#define I44_LAP(i,j,k)  I44_NABLA2_FLD(fields_f,:,i,j,k)

#define I44_D2FLD_XY(i,j,k) I44_D2F_XY_FLD(fields_f,:,i,j,k)
#define I44_D2FLD_XZ(i,j,k) I44_D2F_XZ_FLD(fields_f,:,i,j,k)
#define I44_D2FLD_YZ(i,j,k) I44_D2F_YZ_FLD(fields_f,:,i,j,k)

#define TWELVE_DU_X(fld,i,j,k) TWELVE_DF_X_FLD(metric_h,fld,i,j,k)
#define TWELVE_DU_Y(fld,i,j,k) TWELVE_DF_Y_FLD(metric_h,fld,i,j,k)
#define TWELVE_DU_Z(fld,i,j,k) TWELVE_DF_Z_FLD(metric_h,fld,i,j,k)

#define TWELVE_DV_X(fld,i,j,k) TWELVE_DF_X_FLD(metric_h,fld+3,i,j,k)
#define TWELVE_DV_Y(fld,i,j,k) TWELVE_DF_Y_FLD(metric_h,fld+3,i,j,k)
#define TWELVE_DV_Z(fld,i,j,k) TWELVE_DF_Z_FLD(metric_h,fld+3,i,j,k)

#define I44_D2U_XX(fld,i,j,k) I44_D2F_XX_FLD(metric_h,fld,i,j,k)
#define I44_D2U_YY(fld,i,j,k) I44_D2F_YY_FLD(metric_h,fld,i,j,k)
#define I44_D2U_ZZ(fld,i,j,k) I44_D2F_ZZ_FLD(metric_h,fld,i,j,k)
#define I44_D2U_XY(fld,i,j,k) I44_D2F_XY_FLD(metric_h,fld,i,j,k)
#define I44_D2U_XZ(fld,i,j,k) I44_D2F_XZ_FLD(metric_h,fld,i,j,k)
#define I44_D2U_YZ(fld,i,j,k) I44_D2F_YZ_FLD(metric_h,fld,i,j,k)

#define I44_D2V_XX(fld,i,j,k) I44_D2F_XX_FLD(metric_h,fld+3,i,j,k)
#define I44_D2V_YY(fld,i,j,k) I44_D2F_YY_FLD(metric_h,fld+3,i,j,k)
#define I44_D2V_ZZ(fld,i,j,k) I44_D2F_ZZ_FLD(metric_h,fld+3,i,j,k)
#define I44_D2V_XY(fld,i,j,k) I44_D2F_XY_FLD(metric_h,fld+3,i,j,k)
#define I44_D2V_XZ(fld,i,j,k) I44_D2F_XZ_FLD(metric_h,fld+3,i,j,k)
#define I44_D2V_YZ(fld,i,j,k) I44_D2F_YZ_FLD(metric_h,fld+3,i,j,k)

#endif

!!defining the metric
#if METRIC_PERTURB
#define H11(i,j,k) metric_h(1,sind(i),sind(j),sind(k))
#define H22(i,j,k) metric_h(2,sind(i),sind(j),sind(k))
#define H33(i,j,k) metric_h(3,sind(i),sind(j),sind(k))
#define H23(i,j,k) metric_h(4,sind(i),sind(j),sind(k))
#define H13(i,j,k) metric_h(5,sind(i),sind(j),sind(k))
#define H12(i,j,k) metric_h(6,sind(i),sind(j),sind(k))
#define H32(i,j,k) H23(i,j,k)
#define H31(i,j,k) H13(i,j,k)
#define H21(i,j,k) H12(i,j,k)

#define M_U(fld,i,j,k) metric_h(fld,sind(i),sind(j),sind(k))
#define M_V(fld,i,j,k) metric_h(fld+3,sind(i),sind(j),sind(k))
#define MEU(i,j,k) metric_h(1:3,sind(i),sind(j),sind(k))
#define MEV(i,j,k) metric_h(4:6,sind(i),sind(j),sind(k))

#define MEPIU(i,j,k) metric_p(1:3,sind(i),sind(j),sind(k))
#define MEPIV(i,j,k) metric_p(4:6,sind(i),sind(j),sind(k))
#define TRH(i,j,k) (sum(metric_h(1:3,sind(i),sind(j),sind(k))))
!!this is actually sqrt(g)
#define DETG(i,j,k) fast_exp(TRH(i,j,k)/2._dl)
!!local 1/a^2 approximation

#define GUP11(i,j,k) (fast_exp(-H11(i,j,k))+0.5_dl*(H12(i,j,k)**2*fast_exp(-(H22(i,j,k)+H11(i,j,k)*2.)/3.)+H13(i,j,k)**2*fast_exp(-(H33(i,j,k)+H11(i,j,k)*2.)/3.)))
#define GUP22(i,j,k) (fast_exp(-H22(i,j,k))+0.5_dl*(H21(i,j,k)**2*fast_exp(-(H11(i,j,k)+H22(i,j,k)*2.)/3.)+H23(i,j,k)**2*fast_exp(-(H33(i,j,k)+H22(i,j,k)*2.)/3.)))
#define GUP33(i,j,k) (fast_exp(-H33(i,j,k))+0.5_dl*(H31(i,j,k)**2*fast_exp(-(H11(i,j,k)+H33(i,j,k)*2.)/3.)+H32(i,j,k)**2*fast_exp(-(H22(i,j,k)+H33(i,j,k)*2.)/3.)))
#define GUP23(i,j,k) (-H23(i,j,k)*fast_exp((-H22(i,j,k)-H33(i,j,k))/2._dl)+0.5_dl*H21(i,j,k)*H13(i,j,k)*fast_exp(-TRH(i,j,k)/3._dl))
#define GUP13(i,j,k) (-H13(i,j,k)*fast_exp((-H11(i,j,k)-H33(i,j,k))/2._dl)+0.5_dl*H12(i,j,k)*H23(i,j,k)*fast_exp(-TRH(i,j,k)/3._dl))
#define GUP12(i,j,k) (-H12(i,j,k)*fast_exp((-H11(i,j,k)-H22(i,j,k))/2._dl)+0.5_dl*H13(i,j,k)*H32(i,j,k)*fast_exp(-TRH(i,j,k)/3._dl))
#define DETGUP11(i,j,k) (DETG(i,j,k)*GUP11(i,j,k))
#define DETGUP22(i,j,k) (DETG(i,j,k)*GUP22(i,j,k))
#define DETGUP33(i,j,k) (DETG(i,j,k)*GUP33(i,j,k))
#define DETGUP23(i,j,k) (DETG(i,j,k)*GUP23(i,j,k))
#define DETGUP13(i,j,k) (DETG(i,j,k)*GUP13(i,j,k))
#define DETGUP12(i,j,k) (DETG(i,j,k)*GUP12(i,j,k))
#define CACH_DG(i,j,x) sipar%cach_dg(sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP11(i,j,x) sipar%cach_dgup(1,sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP22(i,j,x) sipar%cach_dgup(2,sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP33(i,j,x) sipar%cach_dgup(3,sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP23(i,j,x) sipar%cach_dgup(4,sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP13(i,j,x) sipar%cach_dgup(5,sind(i),sind(j),sipar%cach_ind(x))
#define CACH_DGUP12(i,j,x) sipar%cach_dgup(6,sind(i),sind(j),sipar%cach_ind(x))
#else
#define H11(i,j,k) 0._dl
#define H22(i,j,k) 0._dl
#define H33(i,j,k) 0._dl
#define H23(i,j,k) 0._dl
#define H13(i,j,k) 0._dl
#define H12(i,j,k) 0._dl
#define H32(i,j,k) H23(i,j,k)
#define H31(i,j,k) H13(i,j,k)
#define H21(i,j,k) H12(i,j,k)
#define DETG(i,j,k) 1._dl
#define GUP11(i,j,k) 1._dl
#define GUP22(i,j,k) 1._dl
#define GUP33(i,j,k) 1._dl
#define GUP23(i,j,k) 0._dl
#define GUP13(i,j,k) 0._dl
#define GUP12(i,j,k) 0._dl
#define DETGUP11(i,j,k) 1._dl
#define DETGUP22(i,j,k) 1._dl
#define DETGUP33(i,j,k) 1._dl
#define DETGUP23(i,j,k) 0._dl
#define DETGUP13(i,j,k) 0._dl
#define DETGUP12(i,j,k) 0._dl
#define CACH_DG(i,j,x) 1._dl
#define CACH_DGUP11(i,j,x) 1._dl
#define CACH_DGUP22(i,j,x) 1._dl
#define CACH_DGUP33(i,j,x) 1._dl
#define CACH_DGUP23(i,j,x) 0._dl
#define CACH_DGUP13(i,j,x) 0._dl
#define CACH_DGUP12(i,j,x) 0._dl
#endif
