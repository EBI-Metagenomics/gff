#include "gff/gff.h"
#include "error.h"
#include "fsm.h"
#include "tok.h"
#include <errno.h>
#include <string.h>

static void buffer_init(struct gff *gff)
{
    gff->buffer.id[0] = '\0';
    gff->buffer.desc[0] = '\0';
    gff->buffer.seq[0] = '\0';
}

void gff_init(struct gff *gff, FILE *restrict fd, enum gff_mode mode)
{
    gff->fd = fd;
    gff->mode = mode;
    gff->target.id = gff->buffer.id;
    gff->target.desc = gff->buffer.desc;
    gff->target.seq = gff->buffer.seq;
    buffer_init(gff);
    fsm_init(&gff->state);
    gff->error[0] = '\0';
    tok_init(&gff->tok, gff->error);
}

enum gff_rc gff_read(struct gff *gff)
{
    if (gff->state == STATE_END) return GFF_ENDFILE;

    if (gff->state != STATE_BEGIN && gff->state != STATE_PAUSE)
        return error_runtime(gff->error, "unexpected %s call", __func__);

    buffer_init(gff);
    if (gff->state == STATE_PAUSE) strcpy(gff->buffer.id, gff->aux.id);

    enum state initial_state = gff->state;
    do
    {
        enum gff_rc rc = GFF_SUCCESS;
        if ((rc = tok_next(&gff->tok, gff->fd))) return rc;

        if ((gff->state = fsm_next(gff->state, &gff->tok, &gff->aux,
                                   &gff->buffer)) == STATE_ERROR)
            return GFF_PARSEERROR;

    } while (gff->state != STATE_PAUSE && gff->state != STATE_END);

    if (gff->state == STATE_END && initial_state == STATE_BEGIN)
        return GFF_ENDFILE;

    return GFF_SUCCESS;
}

void gff_clearerr(struct gff *gff) { gff->error[0] = '\0'; }

enum gff_rc gff_write(struct gff *gff, struct gff_target tgt, unsigned ncols)
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
