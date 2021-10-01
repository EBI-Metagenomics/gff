#ifndef REGION_H
#define REGION_H

#include "gff/region.h"

static inline bool region_set_name(struct gff_region *region, char const *name)
{
    size_t n = gff_strlcpy(region->name, name, GFF_REGION_NAME_SIZE);
    region->start = region->name + n + 1;
    return 0 < n && n < GFF_REGION_NAME_SIZE;
}

static inline bool region_set_start(struct gff_region *region,
                                    char const *start)
{
    size_t n = gff_strlcpy(region->start, start, GFF_REGION_START_SIZE);
    region->end = region->start + n + 1;
    return 0 < n && n < GFF_REGION_START_SIZE;
}

static inline bool region_set_end(struct gff_region *region, char const *end)
{
    size_t n = gff_strlcpy(region->end, end, GFF_REGION_END_SIZE);
    return 0 < n && n < GFF_REGION_END_SIZE;
}

#endif
