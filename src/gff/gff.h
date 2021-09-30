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

static inline bool gff_set_version(struct gff *gff, char const *val)
{
    if (val == NULL)
    {
        gff->elem.type = GFF_ELEM_VERSION;
        gff->elem.version[0] = '3';
        gff->elem.version[1] = '\0';
        return true;
    }
    size_t n = gff_strlcpy(gff->elem.version, val, GFF_VERSION_SIZE);
    bool ok = n > 0 && n < GFF_VERSION_SIZE;
    if (ok) gff->elem.type = GFF_ELEM_VERSION;
    return ok;
}

static inline bool gff_set_region(struct gff *gff, char const *name,
                                  char const *start, char const *end)
{
    gff_region_init(&gff->elem.region);
    if (!gff_rset_name(&gff->elem.region, name)) return false;
    if (!gff_rset_start(&gff->elem.region, start)) return false;
    bool ok = gff_rset_end(&gff->elem.region, end);
    if (ok) gff->elem.type = GFF_ELEM_REGION;
    return ok;
}

static inline struct gff_feature *gff_set_feature(struct gff *gff)
{
    gff->elem.type = GFF_ELEM_FEATURE;
    gff_feature_init(&gff->elem.feature);
    return &gff->elem.feature;
}

#endif
