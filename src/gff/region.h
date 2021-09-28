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
    union
    {
        struct
        {
            char name[GFF_REGION_NAME_SIZE];
            char start[GFF_REGION_START_SIZE];
            char end[GFF_REGION_END_SIZE];
        };
        char buffer[GFF_REGION_SIZE];
    };
};

static_assert(offsetof(struct gff_region, name) == 0, "disallow padding");
static_assert(offsetof(struct gff_region, start) == GFF_REGION_NAME_SIZE,
              "disallow padding");
static_assert(offsetof(struct gff_region, end) ==
                  GFF_REGION_NAME_SIZE + GFF_REGION_START_SIZE,
              "disallow padding");
static_assert(offsetof(struct gff_region, buffer) == 0, "disallow padding");

static inline void gff_region_init(struct gff_region *region)
{
    region->name[0] = '\0';
    region->start[0] = '\0';
    region->end[0] = '\0';
}

#endif
