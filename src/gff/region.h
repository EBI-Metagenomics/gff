#ifndef GFF_REGION_H
#define GFF_REGION_H

#include "gff/strlcpy.h"
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
        char *name;
        char *start;
        char *end;
    };
    char buffer[GFF_REGION_SIZE];
};

#define GFF_REGION_INIT(self, n, s, e)                                         \
    (struct gff_region)                                                        \
    {                                                                          \
        {self.buffer, self.buffer + sizeof(n),                                 \
         self.buffer + sizeof(n) + sizeof(s)},                                 \
            n "\0" s "\0" e                                                    \
    }

static inline void gff_region_init(struct gff_region *region)
{
    region->name = region->buffer;
    region->start = NULL;
    region->end = NULL;
    region->buffer[0] = '\0';
}

#endif
