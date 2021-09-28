#include "fsm.h"
#include "aux.h"
#include "error.h"
#include "gff/aux.h"
#include "gff/elem.h"
#include "gff/tok.h"
#include "tok.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GFF_MATCH_EXCESS_SIZE 5

struct args
{
    struct gff_tok *tok;
    enum state state;
    struct gff_aux *aux;
    struct gff_elem *elem;
};

struct trans
{
    enum state const next;
    enum gff_rc (*action)(struct args *a);
};

static enum gff_rc nop(struct args *a) { return GFF_SUCCESS; }

static enum gff_rc unexpect_eof(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected end-of-file");
}

static enum gff_rc unexpect_tok(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number, "unexpected token");
}

static enum gff_rc unexpect_pragma(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number, "unexpected pragma");
}

static enum gff_rc unexpect_id(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number, "unexpected id");
}

static enum gff_rc unexpect_nl(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected newline");
}

static enum gff_rc desc_begin(struct args *a);
static enum gff_rc desc_cont(struct args *a);
static enum gff_rc read_pragma(struct args *a);
static enum gff_rc read_pragma_id(struct args *a);
static enum gff_rc read_version(struct args *a);
static enum gff_rc set_version_type(struct args *a);
static enum gff_rc set_feature_type(struct args *a);
static enum gff_rc read_feat_seqid(struct args *a);
static enum gff_rc read_feat_source(struct args *a);
static enum gff_rc read_feat_type(struct args *a);
static enum gff_rc read_feat_start(struct args *a);
static enum gff_rc read_feat_end(struct args *a);
static enum gff_rc read_feat_score(struct args *a);
static enum gff_rc read_feat_strand(struct args *a);
static enum gff_rc read_feat_phase(struct args *a);
static enum gff_rc read_feat_attrs(struct args *a);

static struct trans const transition[][6] = {
    [STATE_BEGIN] = {[TOK_NL] = {STATE_BEGIN, &unexpect_nl},
                     [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_VERSION] = {STATE_VERSION, &nop},
                     [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                     [TOK_EOF] = {STATE_END, &nop}},
    [STATE_VERSION] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                       [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_WORD] = {STATE_VERSION_NL, &read_version},
                       [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_VERSION_NL] = {[TOK_NL] = {STATE_PAUSE, &set_version_type},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_NL] = {[TOK_NL] = {STATE_NL, &nop},
                  [TOK_PRAGMA] = {STATE_PAUSE, &store_id},
                  [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                  [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                  [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                  [TOK_EOF] = {STATE_END, &nop}},
    [STATE_PAUSE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                     [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_WORD] = {STATE_FEAT_SOURCE, &read_feat_seqid},
                     [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_SOURCE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_WORD] = {STATE_FEAT_TYPE, &read_feat_source},
                           [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_TYPE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                         [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                         [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                         [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                         [TOK_WORD] = {STATE_FEAT_START, &read_feat_type},
                         [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_START] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_WORD] = {STATE_FEAT_END, &read_feat_start},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                        [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                        [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                        [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                        [TOK_WORD] = {STATE_FEAT_SCORE, &read_feat_end},
                        [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_SCORE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_WORD] = {STATE_FEAT_STRAND, &read_feat_score},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_STRAND] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_WORD] = {STATE_FEAT_PHASE, &read_feat_strand},
                           [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_PHASE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_WORD] = {STATE_FEAT_ATTRS, &read_feat_phase},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_ATTRS] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_WORD] = {STATE_FEAT_NL, &read_feat_attrs},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_NL] = {[TOK_NL] = {STATE_PAUSE, &set_feature_type},
                       [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                       [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                   [TOK_PRAGMA] = {STATE_ERROR, &unexpect_id},
                   [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                   [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                   [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                   [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_ERROR] = {[TOK_NL] = {STATE_ERROR, &nop},
                     [TOK_PRAGMA] = {STATE_ERROR, &nop},
                     [TOK_VERSION] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_REGION] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_WORD] = {STATE_ERROR, &nop},
                     [TOK_EOF] = {STATE_ERROR, &nop}},
};

static char state_name[][14] = {
    [STATE_BEGIN] = "BEGIN",
    [STATE_VERSION] = "VERSION",
    [STATE_VERSION_NL] = "VERSION_NL",
    [STATE_REGION_NAME] = "NAME",
    [STATE_REGION_START] = "START",
    [STATE_REGION_END] = "END",
    [STATE_NL] = "NL",
    [STATE_PAUSE] = "PAUSE",
    [STATE_END] = "END",
    [STATE_ERROR] = "ERROR",
};

enum state fsm_next(enum state state, struct gff_tok *tok, struct gff_aux *aux,
                    struct gff_elem *elem)
{
    unsigned row = (unsigned)state;
    unsigned col = (unsigned)tok->id;
    struct trans const *const t = &transition[row][col];
    struct args args = {tok, state, aux, elem};
    if (t->action(&args)) return STATE_ERROR;
    return t->next;
}

char const *fsm_name(enum state state) { return state_name[state]; }

static enum gff_rc read_desc(struct args *a)
{
    a->aux->pos = memccpy(a->aux->pos - 1, a->tok->value, '\0',
                          (unsigned long)(a->aux->end - a->aux->pos));
    if (!a->aux->pos)
        return error_parse(a->tok->error, a->tok->line.number,
                           "too long description");
    return GFF_SUCCESS;
}

static enum gff_rc desc_begin(struct args *a)
{
    assert(a->tok->id == TOK_WORD || a->tok->id == TOK_PRAGMA);
    /* a->aux->begin = a->tgt->desc; */
    a->aux->pos = a->aux->begin + 1;
    a->aux->end = a->aux->begin + GFF_DESC_MAX;
    return read_desc(a);
}

static enum gff_rc desc_cont(struct args *a)
{
    assert(a->tok->id == TOK_WORD || a->tok->id == TOK_PRAGMA);
    *(a->aux->pos - 1) = ' ';
    a->aux->pos++;
    return read_desc(a);
}

static enum gff_rc read_pragma(struct args *a)
{
    assert(a->tok->id == TOK_PRAGMA);
#if 0
    char const *ptr = memccpy(a->tgt->id, a->tok->value + 1, '\0', GFF_ID_MAX);
    if (!ptr)
        return error_parse(a->tok->error, a->tok->line.number, "too long id");
    if (ptr - a->tgt->id == 1)
        return error_parse(a->tok->error, a->tok->line.number, "empty id");
#endif
    return GFF_SUCCESS;
}

static enum gff_rc read_pragma_id(struct args *a)
{
    assert(a->tok->id == TOK_PRAGMA);
    a->elem->type = GFF_VERSION;
    char const *ptr =
        memccpy(a->elem->version, a->tok->value + 1, '\0', GFF_VERSION_SIZE);
    if (!ptr)
        return error_parse(a->tok->error, a->tok->line.number, "too long id");
    if (ptr - a->elem->version == 1)
        return error_parse(a->tok->error, a->tok->line.number, "empty id");
    return GFF_SUCCESS;
}

static enum gff_rc read_version(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    char const *ptr =
        memccpy(a->elem->version, a->tok->value, '\0', GFF_VERSION_SIZE);
    if (!ptr)
        return error_parse(a->tok->error, a->tok->line.number,
                           "too long version");
    if (ptr - a->elem->version == 1)
        return error_parse(a->tok->error, a->tok->line.number, "empty version");
    return GFF_SUCCESS;
}

static enum gff_rc set_version_type(struct args *a)
{
    assert(a->tok->id == TOK_NL);
    a->elem->type = GFF_VERSION;
    return GFF_SUCCESS;
}

static enum gff_rc set_feature_type(struct args *a)
{
    assert(a->tok->id == TOK_NL);
    a->elem->type = GFF_FEATURE;
    return GFF_SUCCESS;
}

static enum gff_rc read_feat_seqid(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    char const *ptr =
        memccpy(a->elem->version, a->tok->value, '\0', GFF_FEATURE_SEQID_SIZE);
    if (!ptr)
        return error_parse(a->tok->error, a->tok->line.number, "too long seid");
    if (ptr - a->elem->version == 1)
        return error_parse(a->tok->error, a->tok->line.number, "empty seqid");
    return GFF_SUCCESS;
}

static enum gff_rc read_feat_source(struct args *a) {}

static enum gff_rc read_feat_type(struct args *a) {}

static enum gff_rc read_feat_start(struct args *a) {}

static enum gff_rc read_feat_end(struct args *a) {}

static enum gff_rc read_feat_score(struct args *a) {}

static enum gff_rc read_feat_strand(struct args *a) {}

static enum gff_rc read_feat_phase(struct args *a) {}

static enum gff_rc read_feat_attrs(struct args *a) {}

static enum gff_rc seq_begin(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    /* a->aux->begin = a->tgt->seq; */
    a->aux->pos = a->aux->begin + 1;
    a->aux->end = a->aux->begin + GFF_SEQ_MAX;
    return seq_cont(a);
}

static enum gff_rc seq_cont(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    a->aux->pos = memccpy(a->aux->pos - 1, a->tok->value, '\0',
                          (unsigned long)(a->aux->end - a->aux->pos));
    if (!a->aux->pos)
        return error_parse(a->tok->error, a->tok->line.number,
                           "too long sequence");
    return GFF_SUCCESS;
}

static enum gff_rc store_id(struct args *a)
{
    assert(a->tok->id == TOK_PRAGMA);
    char const *ptr = memccpy(a->aux->id, a->tok->value + 1, '\0', GFF_ID_MAX);
    if (!ptr)
        return error_parse(a->tok->error, a->tok->line.number, "too long id");
    if (ptr - a->aux->id == 1)
        return error_parse(a->tok->error, a->tok->line.number, "empty id");
    return GFF_SUCCESS;
}
