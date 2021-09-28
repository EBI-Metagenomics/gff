#ifndef GFF_TARGET_H
#define GFF_TARGET_H

#include "gff/export.h"
#include "gff/rc.h"
#include <stdio.h>

#define GFF_ID_MAX 64
#define GFF_DESC_MAX 128
#define GFF_SEQ_MAX (1024 * 1024)

struct gff;

struct gff_target
{
    char const *id;
    char const *desc;
    char const *seq;
};

static inline struct gff_target gff_target(char const *id, char const *desc,
                                           char const *seq)
{
    return (struct gff_target){id, desc, seq};
}

struct __gff_target
{
    char id[GFF_ID_MAX];
    char desc[GFF_DESC_MAX];
    char seq[GFF_SEQ_MAX];
};

#endif
