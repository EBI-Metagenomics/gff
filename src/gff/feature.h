#ifndef GFF_FEATURE_H
#define GFF_FEATURE_H

#include "gff/strlcpy.h"

#define GFF_FEATURE_SEQID_SIZE 64
#define GFF_FEATURE_SOURCE_SIZE 32
#define GFF_FEATURE_TYPE_SIZE 16
#define GFF_FEATURE_START_SIZE 16
#define GFF_FEATURE_END_SIZE 16
#define GFF_FEATURE_SCORE_SIZE 32
#define GFF_FEATURE_STRAND_SIZE 2
#define GFF_FEATURE_PHASE_SIZE 8
#define GFF_FEATURE_ATTRS_SIZE 256

struct gff_feature
{
    char seqid[GFF_FEATURE_SEQID_SIZE];
    char source[GFF_FEATURE_SOURCE_SIZE];
    char type[GFF_FEATURE_TYPE_SIZE];
    char start[GFF_FEATURE_START_SIZE];
    char end[GFF_FEATURE_END_SIZE];
    char score[GFF_FEATURE_SCORE_SIZE];
    char strand[GFF_FEATURE_STRAND_SIZE];
    char phase[GFF_FEATURE_PHASE_SIZE];
    char attrs[GFF_FEATURE_ATTRS_SIZE];
};

static inline void gff_feature_init(struct gff_feature *feature)
{
    feature->seqid[0] = '\0';
    feature->source[0] = '\0';
    feature->type[0] = '\0';
    feature->start[0] = '\0';
    feature->end[0] = '\0';
    feature->score[0] = '\0';
    feature->strand[0] = '\0';
    feature->phase[0] = '\0';
    feature->attrs[0] = '\0';
}

#define GFF_FEATURE_SET(feat, field, FIELD, val)                               \
    gff_strlcpy((feat)->field, val, GFF_FEATURE_##FIELD##_SIZE)

#define GFF_FSET_SEQID(feat, val) GFF_FEATURE_SET(feat, seqid, SEQID, val)
#define GFF_FSET_SOURCE(feat, val) GFF_FEATURE_SET(feat, source, SOURCE, val)
#define GFF_FSET_TYPE(feat, val) GFF_FEATURE_SET(feat, type, TYPE, val)
#define GFF_FSET_START(feat, val) GFF_FEATURE_SET(feat, start, START, val)
#define GFF_FSET_END(feat, val) GFF_FEATURE_SET(feat, end, END, val)
#define GFF_FSET_SCORE(feat, val) GFF_FEATURE_SET(feat, score, SCORE, val)
#define GFF_FSET_STRAND(feat, val) GFF_FEATURE_SET(feat, strand, STRAND, val)
#define GFF_FSET_PHASE(feat, val) GFF_FEATURE_SET(feat, phase, PHASE, val)
#define GFF_FSET_ATTRS(feat, val) GFF_FEATURE_SET(feat, attrs, ATTRS, val)

#endif
