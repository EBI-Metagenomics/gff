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
    char error[GFF_ERROR_SIZE];
};

GFF_API void gff_init(struct gff *fa, FILE *restrict fd, enum gff_mode mode);

GFF_API enum gff_rc gff_read(struct gff *fa);

GFF_API void gff_clearerr(struct gff *fa);

GFF_API enum gff_rc gff_write(struct gff *fa, struct gff_elem const *elem);

#endif
