#include "gff/gff.h"
#include "hope/hope.h"

void test_example1(void);
void test_example2(void);

int main(void)
{
    test_example1();
    /* test_example2(); */
    return hope_status();
}

void test_example1(void)
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

    GFF_FSET_SEQID(feat, "AE014075.1:534-908|dna");
    GFF_FSET_SOURCE(feat, "iseq");
    GFF_FSET_TYPE(feat, ".");
    GFF_FSET_START(feat, "1");
    GFF_FSET_END(feat, "306");
    GFF_FSET_SCORE(feat, "0.0");
    GFF_FSET_STRAND(feat, "+");
    GFF_FSET_PHASE(feat, ".");
    GFF_FSET_ATTRS(feat,
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
}

void test_example2(void)
{
    FILE *fd = fopen(TMPDIR "/example2.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    gff_set_version(&gff);
    gff.elem.type = GFF_VERSION;
    /* if (gff_strlcpy(args->output_file, arg, PATH_MAX) >= PATH_MAX) */
    /* gff.elem.version[0] = {'3'}; */
    gff.elem.version[0] = '3';
    gff.elem.version[1] = '\0';

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

    GFF_FSET_SEQID(feat, "AE014075.1:534-908|dna");
    GFF_FSET_SOURCE(feat, "iseq");
    GFF_FSET_TYPE(feat, ".");
    GFF_FSET_START(feat, "1");
    GFF_FSET_END(feat, "306");
    GFF_FSET_SCORE(feat, "0.0");
    GFF_FSET_STRAND(feat, "+");
    GFF_FSET_PHASE(feat, ".");
    GFF_FSET_ATTRS(feat,
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
}
