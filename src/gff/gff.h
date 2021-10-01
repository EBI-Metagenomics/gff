#ifndef GFF_GFF_H
#define GFF_GFF_H

#include "gff/elem.h"
#include "gff/error.h"
#include "gff/export.h"
#include "gff/rc.h"
#include "gff/tok.h"
#include <stdio.h>

enum gff_mode
{
    GFF_READ,
    GFF_WRITE,
};

struct gff
{
    FILE *restrict fd;
    enum gff_mode mode;
    struct gff_elem elem;
    unsigned state;
    struct gff_tok tok;
    char *pos;
    bool version_written;
    char error[GFF_ERROR_SIZE];
};

GFF_API void gff_init(struct gff *gff, FILE *restrict fd, enum gff_mode mode);

GFF_API enum gff_rc gff_read(struct gff *gff);

GFF_API void gff_clearerr(struct gff *gff);

GFF_API enum gff_rc gff_write(struct gff *gff);

GFF_API bool gff_set_version(struct gff *gff, char const *val);

GFF_API bool gff_set_region(struct gff *gff, char const *name,
                            char const *start, char const *end);

static inline struct gff_feature *gff_set_feature(struct gff *gff)
{
    gff->elem.type = GFF_ELEM_FEATURE;
    gff_feature_init(&gff->elem.feature);
    return &gff->elem.feature;
}

#endif
