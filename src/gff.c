#include "gff/gff.h"
#include "error.h"
#include "fsm.h"
#include "region.h"
#include "tok.h"
#include "unused.h"
#include <errno.h>
#include <string.h>

void gff_init(struct gff *gff, FILE *restrict fd, enum gff_mode mode)
{
    gff->fd = fd;
    gff->mode = mode;
    gff_elem_init(&gff->elem);
    fsm_init(&gff->state);
    gff->pos = NULL;
    gff->version_written = false;
    gff->error[0] = '\0';
    tok_init(&gff->tok, gff->error);
}

enum gff_rc gff_read(struct gff *gff)
{
    if (gff->state == STATE_END) return GFF_ENDFILE;

    if (gff->state != STATE_BEGIN && gff->state != STATE_PAUSE)
        return error(GFF_RUNTIMEERROR, gff->error, "unexpected gff_read call");

    gff_elem_init(&gff->elem);
    enum state initial_state = gff->state;
    do
    {
        enum gff_rc rc = GFF_SUCCESS;
        if ((rc = tok_next(&gff->tok, gff->fd))) return rc;

        if ((gff->state = fsm_next(gff->state, &gff->tok, &gff->elem,
                                   &gff->pos)) == STATE_ERROR)
            return GFF_PARSEERROR;

    } while (gff->state != STATE_PAUSE && gff->state != STATE_END);

    if (gff->state == STATE_END)
    {
        assert(initial_state == STATE_BEGIN || initial_state == STATE_PAUSE);
        unused(initial_state);
        return GFF_ENDFILE;
    }

    return GFF_SUCCESS;
}

void gff_clearerr(struct gff *gff) { gff->error[0] = '\0'; }

bool gff_set_version(struct gff *gff, char const *val)
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

bool gff_set_region(struct gff *gff, char const *name, char const *start,
                    char const *end)
{
    gff_region_init(&gff->elem.region);
    if (!region_set_name(&gff->elem.region, name)) return false;
    if (!region_set_start(&gff->elem.region, start)) return false;
    bool ok = region_set_end(&gff->elem.region, end);
    if (ok) gff->elem.type = GFF_ELEM_REGION;
    return ok;
}

enum gff_rc write_version(struct gff *gff);
enum gff_rc write_region(struct gff *gff, struct gff_region const *region);
enum gff_rc write_feature(struct gff *gff, struct gff_feature const *feat);

enum gff_rc gff_write(struct gff *gff)
{
    if (!gff->version_written && gff->elem.type != GFF_ELEM_VERSION)
        return error(GFF_ILLEGALARG, gff->error, "write version first");

    if (gff->elem.type == GFF_ELEM_REGION)
        return write_region(gff, &gff->elem.region);
    else if (gff->elem.type == GFF_ELEM_FEATURE)
        return write_feature(gff, &gff->elem.feature);
    else if (gff->elem.type == GFF_ELEM_VERSION)
        return write_version(gff);

    return error(GFF_ILLEGALARG, gff->error, "GFF_ELEM_UNKNOWN element type");
}

enum gff_rc write_version(struct gff *gff)
{
    if (fprintf(gff->fd, "##gff-version %s\n", gff->elem.version) < 0)
        return error_io(gff->error, errno);
    gff->version_written = true;
    return GFF_SUCCESS;
}

#define reg_write(field) fprintf(gff->fd, "%s", reg->field)

enum gff_rc write_region(struct gff *gff, struct gff_region const *reg)
{
    if (fprintf(gff->fd, "##sequence-region") < 0)
        return error_io(gff->error, errno);
    if (fputc(' ', gff->fd) == EOF) return error_io(gff->error, errno);
    if (reg_write(name) < 0) return error_io(gff->error, errno);
    if (fputc(' ', gff->fd) == EOF) return error_io(gff->error, errno);
    if (reg_write(start) < 0) return error_io(gff->error, errno);
    if (fputc(' ', gff->fd) == EOF) return error_io(gff->error, errno);
    if (reg_write(end) < 0) return error_io(gff->error, errno);
    if (fputc('\n', gff->fd) == EOF) return error_io(gff->error, errno);
    return GFF_SUCCESS;
}

#define feat_write(field) fprintf(gff->fd, "%s", feat->field)
#define check_empty(gff, feat, field)                                          \
    if (feat->field[0] == '\0')                                                \
        return error(GFF_ILLEGALARG, gff->error, "empty " #field);

enum gff_rc write_feature(struct gff *gff, struct gff_feature const *feat)
{
    check_empty(gff, feat, seqid);
    check_empty(gff, feat, source);
    check_empty(gff, feat, type);
    check_empty(gff, feat, start);
    check_empty(gff, feat, end);
    check_empty(gff, feat, score);
    check_empty(gff, feat, strand);
    check_empty(gff, feat, phase);
    check_empty(gff, feat, attrs);

    if (feat_write(seqid) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(source) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(type) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(start) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(end) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(score) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(strand) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(phase) < 0) return error_io(gff->error, errno);
    if (fputc('\t', gff->fd) == EOF) return error_io(gff->error, errno);

    if (feat_write(attrs) < 0) return error_io(gff->error, errno);

    if (fputc('\n', gff->fd) == EOF) return error_io(gff->error, errno);

    return GFF_SUCCESS;
}
