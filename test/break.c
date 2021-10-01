#include "gff/gff.h"
#include "hope/hope.h"

#define REP0(X)
#define REP1(X) X
#define REP2(X) REP1(X) X
#define REP3(X) REP2(X) X
#define REP4(X) REP3(X) X
#define REP5(X) REP4(X) X
#define REP6(X) REP5(X) X
#define REP7(X) REP6(X) X
#define REP8(X) REP7(X) X
#define REP9(X) REP8(X) X
#define REP10(X) REP9(X) X

#define REP(HUNDREDS, TENS, ONES, X)                                           \
    REP##HUNDREDS(REP10(REP10(X))) REP##TENS(REP10(X)) REP##ONES(X)

void break_version(void);
void break_region(void);
void break_feature(void);

int main(void)
{
    break_version();
    break_region();
    break_feature();
    return hope_status();
}

void break_version(void)
{
    FILE *fd = fopen(TMPDIR "/version.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    COND(!gff_set_version(&gff, ""));
    COND(!gff_set_version(&gff, REP(0, 1, 6, "1")));
    COND(gff_set_version(&gff, REP(0, 1, 5, "1")));
    EQ(gff_write(&gff), GFF_SUCCESS);

    fclose(fd);
}

void break_region(void)
{
    FILE *fd = fopen(TMPDIR "/region.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    COND(gff_set_version(&gff, "3.1.26"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(!gff_set_region(&gff, "", "1", "1497228"));
    COND(!gff_set_region(&gff, REP(0, 3, 2, "X"), "1", "1497228"));
    COND(gff_set_region(&gff, REP(0, 3, 1, "X"), "1", "1497228"));

    COND(!gff_set_region(&gff, "X", "", "1497228"));
    COND(!gff_set_region(&gff, "X", REP(0, 1, 6, "1"), "1497228"));
    COND(gff_set_region(&gff, "X", REP(0, 1, 5, "1"), "1497228"));

    COND(!gff_set_region(&gff, "X", "1", ""));
    COND(!gff_set_region(&gff, "X", "1", REP(0, 1, 6, "2")));
    COND(gff_set_region(&gff, "X", "1", REP(0, 1, 5, "2")));

    EQ(gff_write(&gff), GFF_SUCCESS);

    fclose(fd);
}

void break_feature(void)
{
    FILE *fd = fopen(TMPDIR "/feature.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    COND(gff_set_version(&gff, "3.1.26"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    struct gff_feature *feat = gff_set_feature(&gff);

    COND(gff_feature_set_seqid(feat, "AE014075.1:190-252|dna"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_source(feat, "iseq"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_type(feat, "."));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_start(feat, "1"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_end(feat, "63"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_score(feat, "0.0"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_strand(feat, "+"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_phase(feat, "."));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(gff_feature_set_attrs(feat, "ID=item1;Target_alph=dna"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(!gff_feature_set_seqid(feat, ""));
    COND(!gff_feature_set_seqid(feat, REP(0, 6, 4, "x")));
    COND(gff_feature_set_seqid(feat, REP(0, 6, 3, "x")));

    COND(!gff_feature_set_source(feat, ""));
    COND(!gff_feature_set_source(feat, REP(0, 3, 2, "x")));
    COND(gff_feature_set_source(feat, REP(0, 3, 1, "x")));

    COND(!gff_feature_set_type(feat, ""));
    COND(!gff_feature_set_type(feat, REP(0, 1, 6, "x")));
    COND(gff_feature_set_type(feat, REP(0, 1, 5, "x")));

    COND(!gff_feature_set_start(feat, ""));
    COND(!gff_feature_set_start(feat, REP(0, 1, 6, "x")));
    COND(gff_feature_set_start(feat, REP(0, 1, 5, "x")));

    COND(!gff_feature_set_end(feat, ""));
    COND(!gff_feature_set_end(feat, REP(0, 1, 6, "x")));
    COND(gff_feature_set_end(feat, REP(0, 1, 5, "x")));

    COND(!gff_feature_set_score(feat, ""));
    COND(!gff_feature_set_score(feat, REP(0, 3, 2, "x")));
    COND(gff_feature_set_score(feat, REP(0, 3, 1, "x")));

    COND(!gff_feature_set_strand(feat, ""));
    COND(!gff_feature_set_strand(feat, REP(0, 0, 2, "x")));
    COND(gff_feature_set_strand(feat, REP(0, 0, 1, "x")));

    COND(!gff_feature_set_phase(feat, ""));
    COND(!gff_feature_set_phase(feat, REP(0, 0, 8, "x")));
    COND(gff_feature_set_phase(feat, REP(0, 0, 7, "x")));

    COND(!gff_feature_set_attrs(feat, ""));
    COND(!gff_feature_set_attrs(feat, REP(2, 5, 6, "x")));
    COND(gff_feature_set_attrs(feat, REP(2, 5, 5, "x")));

    fclose(fd);
}
