        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul  6 18:38:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_MODEL_PARAMETERS__genmod
          INTERFACE 
            SUBROUTINE WRITE_MODEL_PARAMETERS(MODEL,EN_VAR)
              USE LATERATION
              TYPE (MODEL_PARAM), INTENT(IN) :: MODEL
              CHARACTER(LEN=20) :: EN_VAR
            END SUBROUTINE WRITE_MODEL_PARAMETERS
          END INTERFACE 
        END MODULE WRITE_MODEL_PARAMETERS__genmod
