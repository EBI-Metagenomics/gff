#ifndef GFF_ELEM_H
#define GFF_ELEM_H

#include "gff/feature.h"
#include "gff/region.h"

enum gff_elem_type
{
    GFF_ELEM_UNKNOWN,
    GFF_ELEM_VERSION,
    GFF_ELEM_REGION,
    GFF_ELEM_FEATURE,
};

#define GFF_VERSION_SIZE 16

struct gff_elem
{
    enum gff_elem_type type;
    union
    {
        char version[GFF_VERSION_SIZE];
        struct gff_region region;
        struct gff_feature feature;
    };
};

static inline void gff_elem_init(struct gff_elem *elem)
{
    elem->type = GFF_ELEM_UNKNOWN;
    elem->version[0] = '\0';
    gff_region_init(&elem->region);
    gff_feature_init(&elem->feature);
}

#endif
