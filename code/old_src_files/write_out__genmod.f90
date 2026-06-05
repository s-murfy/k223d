        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul  6 18:38:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_OUT__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT(AMESH,OUT_TYPE,OUT_FILE,PDF,EN_VAR,    &
     &OUTPUT_LEVEL)
              USE TYPEDEF
              TYPE (MESH) :: AMESH
              CHARACTER(LEN=3) :: OUT_TYPE
              CHARACTER(LEN=30) :: OUT_FILE
              REAL(KIND=4) :: PDF(AMESH%QUAKEELEMNO)
              CHARACTER(LEN=20) :: EN_VAR
              INTEGER(KIND=4) :: OUTPUT_LEVEL
            END SUBROUTINE WRITE_OUT
          END INTERFACE 
        END MODULE WRITE_OUT__genmod
