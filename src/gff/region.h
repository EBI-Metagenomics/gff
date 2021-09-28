#ifndef GFF_REGION_H
#define GFF_REGION_H

#include <assert.h>
#include <stddef.h>

#define GFF_REGION_NAME_SIZE 32
#define GFF_REGION_START_SIZE 16
#define GFF_REGION_END_SIZE 16

#define GFF_REGION_SIZE                                                        \
    (GFF_REGION_NAME_SIZE + GFF_REGION_START_SIZE + GFF_REGION_END_SIZE)

struct gff_region
{
    struct
    {
        char const *name;
        char const *start;
        char const *end;
    };
    char buffer[GFF_REGION_SIZE];
};

static inline void gff_region_init(struct gff_region *region)
{
    region->name = region->buffer;
    region->start = NULL;
    region->end = NULL;
    region->buffer[0] = '\0';
}

#endif
