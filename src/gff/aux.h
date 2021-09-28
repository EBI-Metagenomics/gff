#ifndef GFF_AUX_H
#define GFF_AUX_H

#include "gff/target.h"

struct gff_aux
{
    char *begin;
    char *pos;
    char *end;
    char id[GFF_ID_MAX];
};

#endif
