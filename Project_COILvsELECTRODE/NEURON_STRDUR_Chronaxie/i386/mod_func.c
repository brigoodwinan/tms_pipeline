#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," kvz_nature.mod");
    fprintf(stderr," naz_nature.mod");
    fprintf(stderr," vsource.mod");
    fprintf(stderr," xtra.mod");
    fprintf(stderr, "\n");
  }
  _kvz_nature_reg();
  _naz_nature_reg();
  _vsource_reg();
  _xtra_reg();
}
