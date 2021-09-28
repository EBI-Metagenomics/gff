#ifndef AUX_H
#define AUX_H

#include "gff/aux.h"

static inline void aux_init(struct gff_aux *aux)
{
    aux->begin = NULL;
    aux->pos = NULL;
    aux->end = NULL;
    aux->id[0] = '\0';
}

#endif
