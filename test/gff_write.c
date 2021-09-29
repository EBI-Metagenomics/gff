#include "gff/gff.h"
#include "hope/hope.h"

int main(void)
{
    FILE *fd = fopen(TMPDIR "/example1.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    gff_set_version(&gff);
    enum gff_rc rc = gff_write(&gff);
    EQ(rc, GFF_SUCCESS);

    struct gff_feature *feat = gff_set_feature(&gff);

    GFF_FEATURE_SET(feat, seqid, SEQID, "AE014075.1:190-252|dna");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, source, SOURCE, "iseq");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, type, TYPE, ".");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, start, START, "1");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, end, END, "63");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, score, SCORE, "0.0");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, strand, STRAND, "+");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, phase, PHASE, ".");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, attrs, ATTRS,
                    "ID=item1;Target_alph=dna;Profile_name=Leader_Thr;Profile_"
                    "alph=dna;Profile_acc=PF08254.12;Window=0;Bias=17.5;E-"
                    "value=2.9e-14;Epsilon=0.01;Score=38.8");
    EQ(gff_write(&gff), GFF_SUCCESS);

    feat = gff_set_feature(&gff);

    GFF_FEATURE_SET(feat, seqid, SEQID, "AE014075.1:534-908|dna");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, source, SOURCE, "iseq");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, type, TYPE, ".");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, start, START, "1");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, end, END, "306");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, score, SCORE, "0.0");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, strand, STRAND, "+");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, phase, PHASE, ".");
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    GFF_FEATURE_SET(feat, attrs, ATTRS,
                    "ID=item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph="
                    "dna;Profile_acc=PF01797.17;Window=0;Bias=0.0;E-value=1.7e-"
                    "29;Epsilon=0.01;Score=88.6");
    EQ(gff_write(&gff), GFF_SUCCESS);

    fclose(fd);

    FILE *actual = fopen(TMPDIR "/example1.gff", "r");
    FILE *desired = fopen(ASSETS "/example1.gff", "r");
    NOTNULL(actual);
    NOTNULL(desired);
    EQ(actual, desired);
    fclose(desired);
    fclose(actual);

    return hope_status();
}
