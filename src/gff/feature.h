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

#define GFF_FEATURE_DEF(f, F)                                                  \
    static inline bool gff_feature_set_##f(struct gff_feature *feat,           \
                                           char const *val)                    \
    {                                                                          \
        size_t n = gff_strlcpy(feat->f, val, GFF_FEATURE_##F##_SIZE);          \
        return 0 < n && n < GFF_FEATURE_##F##_SIZE;                            \
    }

GFF_FEATURE_DEF(seqid, SEQID);
GFF_FEATURE_DEF(source, SOURCE);
GFF_FEATURE_DEF(type, TYPE);
GFF_FEATURE_DEF(start, START);
GFF_FEATURE_DEF(end, END);
GFF_FEATURE_DEF(score, SCORE);
GFF_FEATURE_DEF(strand, STRAND);
GFF_FEATURE_DEF(phase, PHASE);
GFF_FEATURE_DEF(attrs, ATTRS);

#endif
