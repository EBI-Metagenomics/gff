#include "fsm.h"
#include "error.h"
#include "gff/aux.h"
#include "gff/elem.h"
#include "gff/tok.h"
#include "tok.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

static enum gff_rc unexpect_comment(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected comment");
}

static enum gff_rc unexpect_pragma(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected directive");
}

static enum gff_rc unexpect_version(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected version directive");
}

static enum gff_rc unexpect_region(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected region directive");
}

static enum gff_rc unexpect_fasta(struct args *a)
{
    return error_parse(a->tok->error, a->tok->line.number,
                       "unexpected fasta directive");
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

static enum gff_rc read_version(struct args *a);
static enum gff_rc read_region_name(struct args *a);
static enum gff_rc read_region_start(struct args *a);
static enum gff_rc read_region_end(struct args *a);
static enum gff_rc read_seqid(struct args *a);
static enum gff_rc read_source(struct args *a);
static enum gff_rc read_type(struct args *a);
static enum gff_rc read_start(struct args *a);
static enum gff_rc read_end(struct args *a);
static enum gff_rc read_score(struct args *a);
static enum gff_rc read_strand(struct args *a);
static enum gff_rc read_phase(struct args *a);
static enum gff_rc read_attrs_init(struct args *a);
static enum gff_rc read_attrs_cont(struct args *a);
static enum gff_rc set_version_type(struct args *a);
static enum gff_rc set_region_type(struct args *a);
static enum gff_rc set_feature_type(struct args *a);

static struct trans const transition[][8] = {
    [STATE_BEGIN] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                     [TOK_COMMENT] = {STATE_COMMENT, &nop},
                     [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_VERSION] = {STATE_VERSION, &nop},
                     [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                     [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                     [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                     [TOK_EOF] = {STATE_END, &nop}},
    [STATE_VERSION] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                       [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                       [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                       [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                       [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                       [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                       [TOK_WORD] = {STATE_VERSION_NL, &read_version},
                       [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_VERSION_NL] = {[TOK_NL] = {STATE_PAUSE, &set_version_type},
                          [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                          [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                          [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_COMMENT] = {[TOK_NL] = {STATE_COMMENT_END, &nop},
                       [TOK_COMMENT] = {STATE_COMMENT, &nop},
                       [TOK_PRAGMA] = {STATE_COMMENT, &nop},
                       [TOK_VERSION] = {STATE_COMMENT, &nop},
                       [TOK_REGION] = {STATE_COMMENT, &nop},
                       [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                       [TOK_WORD] = {STATE_COMMENT, &nop},
                       [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_COMMENT_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_COMMENT] = {STATE_COMMENT, &nop},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                           [TOK_REGION] = {STATE_REGION_NAME, &nop},
                           [TOK_FASTA] = {STATE_END, &nop},
                           [TOK_WORD] = {STATE_FEAT_SOURCE, &read_seqid},
                           [TOK_EOF] = {STATE_END, &nop}},
    [STATE_REGION_NAME] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                           [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                           [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                           [TOK_WORD] = {STATE_REGION_START, &read_region_name},
                           [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_REGION_START] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                            [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                            [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                            [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                            [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                            [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                            [TOK_WORD] = {STATE_REGION_END, &read_region_start},
                            [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_REGION_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                          [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                          [TOK_WORD] = {STATE_REGION_NL, &read_region_end},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_REGION_NL] = {[TOK_NL] = {STATE_PAUSE, &set_region_type},
                         [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                         [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                         [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                         [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                         [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                         [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                         [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_SOURCE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                           [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                           [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                           [TOK_WORD] = {STATE_FEAT_TYPE, &read_source},
                           [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_TYPE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                         [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                         [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                         [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                         [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                         [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                         [TOK_WORD] = {STATE_FEAT_START, &read_type},
                         [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_START] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                          [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                          [TOK_WORD] = {STATE_FEAT_END, &read_start},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                        [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                        [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                        [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                        [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                        [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                        [TOK_WORD] = {STATE_FEAT_SCORE, &read_end},
                        [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_SCORE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                          [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                          [TOK_WORD] = {STATE_FEAT_STRAND, &read_score},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_STRAND] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                           [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                           [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                           [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                           [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                           [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                           [TOK_WORD] = {STATE_FEAT_PHASE, &read_strand},
                           [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_PHASE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                          [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                          [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                          [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                          [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                          [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                          [TOK_WORD] = {STATE_FEAT_ATTRS_INIT, &read_phase},
                          [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_ATTRS_INIT] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                               [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                               [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                               [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                               [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                               [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                               [TOK_WORD] = {STATE_FEAT_ATTRS_CONT,
                                             &read_attrs_init},
                               [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_FEAT_ATTRS_CONT] = {[TOK_NL] = {STATE_PAUSE, &set_feature_type},
                               [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                               [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                               [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                               [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                               [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                               [TOK_WORD] = {STATE_FEAT_ATTRS_CONT,
                                             &read_attrs_cont},
                               [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_PAUSE] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                     [TOK_COMMENT] = {STATE_COMMENT, &nop},
                     [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                     [TOK_REGION] = {STATE_REGION_NAME, &nop},
                     [TOK_FASTA] = {STATE_END, &nop},
                     [TOK_WORD] = {STATE_FEAT_SOURCE, &read_seqid},
                     [TOK_EOF] = {STATE_END, &nop}},
    [STATE_END] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                   [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                   [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                   [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                   [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                   [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                   [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                   [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
    [STATE_ERROR] = {[TOK_NL] = {STATE_ERROR, &unexpect_nl},
                     [TOK_COMMENT] = {STATE_ERROR, &unexpect_comment},
                     [TOK_PRAGMA] = {STATE_ERROR, &unexpect_pragma},
                     [TOK_VERSION] = {STATE_ERROR, &unexpect_version},
                     [TOK_REGION] = {STATE_ERROR, &unexpect_region},
                     [TOK_FASTA] = {STATE_ERROR, &unexpect_fasta},
                     [TOK_WORD] = {STATE_ERROR, &unexpect_tok},
                     [TOK_EOF] = {STATE_ERROR, &unexpect_eof}},
};

static char state_name[][16] = {[STATE_BEGIN] = "BEGIN",
                                [STATE_VERSION] = "VERSION",
                                [STATE_VERSION_NL] = "VERSION_NL",
                                [STATE_COMMENT] = "COMMENT",
                                [STATE_COMMENT_END] = "COMMENT_END",
                                [STATE_REGION_NAME] = "REGION_NAME",
                                [STATE_REGION_START] = "REGION_START",
                                [STATE_REGION_END] = "REGION_END",
                                [STATE_REGION_NL] = "REGION_NL",
                                [STATE_FEAT_SOURCE] = "FEAT_SOURCE",
                                [STATE_FEAT_TYPE] = "FEAT_TYPE",
                                [STATE_FEAT_START] = "FEAT_START",
                                [STATE_FEAT_END] = "FEAT_END",
                                [STATE_FEAT_SCORE] = "FEAT_SCORE",
                                [STATE_FEAT_STRAND] = "FEAT_STRAND",
                                [STATE_FEAT_PHASE] = "FEAT_PHASE",
                                [STATE_FEAT_ATTRS_INIT] = "FEAT_ATTRS_INIT",
                                [STATE_FEAT_ATTRS_CONT] = "FEAT_ATTRS_CONT",
                                [STATE_PAUSE] = "PAUSE",
                                [STATE_END] = "END",
                                [STATE_ERROR] = "ERROR"};

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

static enum gff_rc tokcpy0(char *dst, struct gff_tok *tok, size_t count,
                           char const *name, char **ptr);
static enum gff_rc tokcpy(char *dst, struct gff_tok *tok, size_t count,
                          char const *name);

static enum gff_rc read_version(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    return tokcpy(a->elem->version, a->tok, GFF_VERSION_SIZE, "version");
}

static enum gff_rc read_region_name(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_region *r = &a->elem->region;
    gff_region_init(r);
    enum gff_rc rc = tokcpy(r->buffer, a->tok, GFF_REGION_NAME_SIZE, "name");
    if (rc) return rc;
    r->start = r->name + strlen(r->name) + 1;
    return rc;
}

static enum gff_rc read_region_start(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_region *r = &a->elem->region;
    size_t n = (size_t)(r->start - r->name);
    enum gff_rc rc =
        tokcpy(r->buffer + n, a->tok, GFF_REGION_START_SIZE, "start");
    if (rc) return rc;
    r->end = r->start + strlen(r->start) + 1;
    return rc;
}

static enum gff_rc read_region_end(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_region *r = &a->elem->region;
    size_t n = (size_t)(r->end - r->name);
    return tokcpy(r->buffer + n, a->tok, GFF_REGION_END_SIZE, "end");
}

static enum gff_rc set_version_type(struct args *a)
{
    assert(a->tok->id == TOK_NL);
    a->elem->type = GFF_VERSION;
    return GFF_SUCCESS;
}

static enum gff_rc set_region_type(struct args *a)
{
    assert(a->tok->id == TOK_NL);
    a->elem->type = GFF_REGION;
    return GFF_SUCCESS;
}

static enum gff_rc set_feature_type(struct args *a)
{
    assert(a->tok->id == TOK_NL);
    a->elem->type = GFF_FEATURE;
    return GFF_SUCCESS;
}

static enum gff_rc read_seqid(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->seqid, a->tok, GFF_FEATURE_SEQID_SIZE, "seqid");
}

static enum gff_rc read_source(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->source, a->tok, GFF_FEATURE_SOURCE_SIZE, "source");
}

static enum gff_rc read_type(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->type, a->tok, GFF_FEATURE_TYPE_SIZE, "type");
}

static enum gff_rc read_start(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->start, a->tok, GFF_FEATURE_START_SIZE, "start");
}

static enum gff_rc read_end(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->end, a->tok, GFF_FEATURE_END_SIZE, "end");
}

static enum gff_rc read_score(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->score, a->tok, GFF_FEATURE_SCORE_SIZE, "score");
}

static enum gff_rc read_strand(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->strand, a->tok, GFF_FEATURE_STRAND_SIZE, "strand");
}

static enum gff_rc read_phase(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    return tokcpy(f->phase, a->tok, GFF_FEATURE_PHASE_SIZE, "phase");
}

static enum gff_rc read_attrs_init(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    a->aux->pos = a->elem->feature.attrs;
    /* struct gff_feature *f = &a->elem->feature; */
    return tokcpy0(a->aux->pos, a->tok, GFF_FEATURE_ATTRS_SIZE, "attributes",
                   &a->aux->pos);
}

static enum gff_rc read_attrs_cont(struct args *a)
{
    assert(a->tok->id == TOK_WORD);
    struct gff_feature *f = &a->elem->feature;
    *(a->aux->pos - 1) = ' ';
    if (GFF_FEATURE_ATTRS_SIZE + f->attrs < a->aux->pos)
        return error_parse(a->tok->error, a->tok->line.number,
                           "too long attributes");
    size_t n = (size_t)(GFF_FEATURE_ATTRS_SIZE - (a->aux->pos - f->attrs));
    return tokcpy0(a->aux->pos, a->tok, n, "attributes", &a->aux->pos);
}

static enum gff_rc tokcpy0(char *dst, struct gff_tok *tok, size_t count,
                           char const *name, char **ptr)
{
    *ptr = memccpy(dst, tok->value, '\0', count);
    if (!ptr)
        return error_parse(tok->error, tok->line.number, "too long %s", name);
    if (*ptr - dst == 1)
        return error_parse(tok->error, tok->line.number, "empty %s", name);
    return GFF_SUCCESS;
}

static enum gff_rc tokcpy(char *dst, struct gff_tok *tok, size_t count,
                          char const *name)
{
    char *ptr = NULL;
    return tokcpy0(dst, tok, count, name, &ptr);
}
