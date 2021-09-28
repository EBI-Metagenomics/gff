#ifndef GFF_REGION_H
#define GFF_REGION_H

#define GFF_REGION_NAME_SIZE 32
#define GFF_REGION_START_SIZE 16
#define GFF_REGION_END_SIZE 16

struct gff_region
{
    char name[GFF_REGION_NAME_SIZE];
    char start[GFF_REGION_START_SIZE];
    char end[GFF_REGION_END_SIZE];
};

static inline void gff_region_init(struct gff_region *region)
{
    region->name[0] = '\0';
    region->start[0] = '\0';
    region->end[0] = '\0';
}

#endif
