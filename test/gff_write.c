#include "gff/gff.h"
#include "hope/hope.h"

void test_example1(void);
void test_example4(void);
void test_wrong_usage(void);

int main(void)
{
    test_example1();
    test_example4();
    test_wrong_usage();
    return hope_status();
}

void test_example1(void)
{
    FILE *fd = fopen(TMPDIR "/example1.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    COND(gff_set_version(&gff, NULL));
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

    COND(gff_feature_set_attrs(
        feat, "ID=item1;Target_alph=dna;Profile_name=Leader_Thr;Profile_"
              "alph=dna;Profile_acc=PF08254.12;Window=0;Bias=17.5;E-"
              "value=2.9e-14;Epsilon=0.01;Score=38.8"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    feat = gff_set_feature(&gff);

    COND(gff_feature_set_seqid(feat, "AE014075.1:534-908|dna"));
    COND(gff_feature_set_source(feat, "iseq"));
    COND(gff_feature_set_type(feat, "."));
    COND(gff_feature_set_start(feat, "1"));
    COND(gff_feature_set_end(feat, "306"));
    COND(gff_feature_set_score(feat, "0.0"));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(
        feat, "ID=item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph="
              "dna;Profile_acc=PF01797.17;Window=0;Bias=0.0;E-value=1.7e-"
              "29;Epsilon=0.01;Score=88.6"));
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

void test_example4(void)
{
    FILE *fd = fopen(TMPDIR "/example4.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    COND(gff_set_version(&gff, "3.1.26"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_set_region(&gff, "ctg123", "1", "1497228"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    struct gff_feature *feat = gff_set_feature(&gff);
    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "gene"));
    COND(gff_feature_set_start(feat, "1000"));
    COND(gff_feature_set_end(feat, "9000"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat, "ID=gene00001;Name=EDEN"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "TF_binding_site"));
    COND(gff_feature_set_start(feat, "1000"));
    COND(gff_feature_set_end(feat, "1012"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat, "ID=tfbs00001;Parent=gene00001"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "mRNA"));
    COND(gff_feature_set_start(feat, "1050"));
    COND(gff_feature_set_end(feat, "9000"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat,
                               "ID=mRNA00001;Parent=gene00001;Name=EDEN.1"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "mRNA"));
    COND(gff_feature_set_start(feat, "1050"));
    COND(gff_feature_set_end(feat, "9000"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat,
                               "ID=mRNA00002;Parent=gene00001;Name=EDEN.2"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "mRNA"));
    COND(gff_feature_set_start(feat, "1300"));
    COND(gff_feature_set_end(feat, "9000"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat,
                               "ID=mRNA00003;Parent=gene00001;Name=EDEN.3"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "exon"));
    COND(gff_feature_set_start(feat, "1300"));
    COND(gff_feature_set_end(feat, "1500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat, "ID=exon00001;Parent=mRNA00003"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "exon"));
    COND(gff_feature_set_start(feat, "1050"));
    COND(gff_feature_set_end(feat, "1500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(
        gff_feature_set_attrs(feat, "ID=exon00002;Parent=mRNA00001,mRNA00002"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "exon"));
    COND(gff_feature_set_start(feat, "3000"));
    COND(gff_feature_set_end(feat, "3902"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(
        gff_feature_set_attrs(feat, "ID=exon00003;Parent=mRNA00001,mRNA00003"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "exon"));
    COND(gff_feature_set_start(feat, "5000"));
    COND(gff_feature_set_end(feat, "5500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(
        feat, "ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "exon"));
    COND(gff_feature_set_start(feat, "7000"));
    COND(gff_feature_set_end(feat, "9000"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(
        feat, "ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "1201"));
    COND(gff_feature_set_end(feat, "1500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "3000"));
    COND(gff_feature_set_end(feat, "3902"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "5000"));
    COND(gff_feature_set_end(feat, "5500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "7000"));
    COND(gff_feature_set_end(feat, "7600"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "1201"));
    COND(gff_feature_set_end(feat, "1500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "5000"));
    COND(gff_feature_set_end(feat, "5500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "7000"));
    COND(gff_feature_set_end(feat, "7600"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "3301"));
    COND(gff_feature_set_end(feat, "3902"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "5000"));
    COND(gff_feature_set_end(feat, "5500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "1"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "7000"));
    COND(gff_feature_set_end(feat, "7600"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "1"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "3391"));
    COND(gff_feature_set_end(feat, "3902"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "0"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "5000"));
    COND(gff_feature_set_end(feat, "5500"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "1"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    COND(gff_feature_set_seqid(feat, "ctg123"));
    COND(gff_feature_set_source(feat, "."));
    COND(gff_feature_set_type(feat, "CDS"));
    COND(gff_feature_set_start(feat, "7000"));
    COND(gff_feature_set_end(feat, "7600"));
    COND(gff_feature_set_score(feat, "."));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "1"));
    COND(gff_feature_set_attrs(
        feat, "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    fclose(fd);

    FILE *actual = fopen(TMPDIR "/example4.gff", "r");
    FILE *desired = fopen(ASSETS "/example4.gff", "r");
    NOTNULL(actual);
    NOTNULL(desired);
    EQ(actual, desired);
    fclose(desired);
    fclose(actual);
}

void test_wrong_usage(void)
{
    FILE *fd = fopen(TMPDIR "/wrong_usage1.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    struct gff_feature *feat = gff_set_feature(&gff);
    COND(gff_feature_set_seqid(feat, "AE014075.1:190-252|dna"));
    COND(gff_feature_set_source(feat, "iseq"));
    COND(gff_feature_set_type(feat, "."));
    COND(gff_feature_set_start(feat, "1"));
    COND(gff_feature_set_end(feat, "63"));
    COND(gff_feature_set_score(feat, "0.0"));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(feat, "value=2.9e-14"));
    EQ(gff_write(&gff), GFF_ILLEGALARG);

    COND(!gff_set_version(&gff, ""));
    EQ(gff_write(&gff), GFF_ILLEGALARG);
    EQ(gff.error, "Runtime error: write version first");
    gff_clearerr(&gff);
    EQ(gff.error, "");

    COND(gff_set_version(&gff, NULL));
    EQ(gff_write(&gff), GFF_SUCCESS);

    feat = gff_set_feature(&gff);
    COND(gff_feature_set_seqid(feat, "AE014075.1:534-908|dna"));
    COND(gff_feature_set_source(feat, "iseq"));
    COND(gff_feature_set_type(feat, "."));
    COND(gff_feature_set_start(feat, "1"));
    COND(gff_feature_set_end(feat, "306"));
    COND(gff_feature_set_score(feat, "0.0"));
    COND(!gff_feature_set_strand(feat, "++++"));
    COND(gff_feature_set_strand(feat, "+"));
    COND(gff_feature_set_phase(feat, "."));
    COND(gff_feature_set_attrs(
        feat, "ID=item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_alph="
              "dna;Profile_acc=PF01797.17;Window=0;Bias=0.0;E-value=1.7e-"
              "29;Epsilon=0.01;Score=88.6"));
    EQ(gff_write(&gff), GFF_SUCCESS);

    fclose(fd);
}
