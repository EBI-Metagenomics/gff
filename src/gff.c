#include "gff/gff.h"
#include "error.h"
#include "fsm.h"
#include "tok.h"
#include <errno.h>
#include <string.h>

void gff_init(struct gff *gff, FILE *restrict fd, enum gff_mode mode)
{
    gff->fd = fd;
    gff->mode = mode;
    gff_elem_init(&gff->elem);
    fsm_init(&gff->state);
    gff->error[0] = '\0';
    tok_init(&gff->tok, gff->error);
}

enum gff_rc gff_read(struct gff *gff)
{
    if (gff->state == STATE_END) return GFF_ENDFILE;

    if (gff->state != STATE_BEGIN && gff->state != STATE_PAUSE)
        return error_runtime(gff->error, "unexpected %s call", __func__);

    gff_elem_init(&gff->elem);
    enum state initial_state = gff->state;
    do
    {
        enum gff_rc rc = GFF_SUCCESS;
        if ((rc = tok_next(&gff->tok, gff->fd))) return rc;

        if ((gff->state = fsm_next(gff->state, &gff->tok, &gff->aux,
                                   &gff->elem)) == STATE_ERROR)
            return GFF_PARSEERROR;

    } while (gff->state != STATE_PAUSE && gff->state != STATE_END);

    if (gff->state == STATE_END)
    {
        assert(initial_state == STATE_BEGIN || initial_state == STATE_PAUSE);
        return GFF_ENDFILE;
    }

    return GFF_SUCCESS;
}

void gff_clearerr(struct gff *gff) { gff->error[0] = '\0'; }

#if 0
enum gff_rc gff_write(struct gff *gff, struct gff_elem *elem)
{
    if (fprintf(gff->fd, ">%s", tgt.id) < 0)
        return error_io("failed to write", errno);

    if (tgt.desc[0] && fprintf(gff->fd, " %s", tgt.desc) < 0)
        return error_io("failed to write", errno);

    for (char const *c = tgt.seq; *c; ++c)
    {
        if (((c - tgt.seq) % ncols) == 0)
        {
            if (fputc('\n', gff->fd) == EOF)
                return error_io("failed to write", errno);
        }
        if (fputc(*c, gff->fd) == EOF)
            return error_io("failed to write", errno);
    }
    if (fputc('\n', gff->fd) == EOF) return error_io("failed to write", errno);
    return GFF_SUCCESS;
}
#endif

enum gff_rc write_region(struct gff *gff, struct gff_region const *region);
enum gff_rc write_feature(struct gff *gff, struct gff_feature const *feat);

#define nl_wfail() error_io("failed to write newline", errno)
#define tab_wfail() error_io("failed to write tab", errno)
#define feat_wfail(field) error_io("failed to write feature " field, errno)
#define feat_write(field) fprintf(gff->fd, "%s", feat->field)

enum gff_rc write_region(struct gff *gff, struct gff_region const *reg)
{
    if (fprintf(gff->fd, "%s", reg->name) < 0)
        return error_io("failed to write region name", errno);

    if (fprintf(gff->fd, "%s", reg->start) < 0)
        return error_io("failed to write region start", errno);

    if (fprintf(gff->fd, "%s", reg->end) < 0)
        return error_io("failed to write region end", errno);

    if (fputc('\n', gff->fd) == EOF) return nl_wfail();
    return GFF_SUCCESS;
}

enum gff_rc write_feature(struct gff *gff, struct gff_feature const *feat)
{
    if (feat_write(seqid) < 0) return feat_wfail("seqid");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(source) < 0) return feat_wfail("source");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(type) < 0) return feat_wfail("type");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(start) < 0) return feat_wfail("start");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(end) < 0) return feat_wfail("end");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(score) < 0) return feat_wfail("score");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(strand) < 0) return feat_wfail("strand");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(phase) < 0) return feat_wfail("phase");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    if (feat_write(attrs) < 0) return feat_wfail("attrs");
    if (fputc('\t', gff->fd) == EOF) return tab_wfail();

    return GFF_SUCCESS;
}
