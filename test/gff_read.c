#include "gff/gff.h"
#include "hope/hope.h"

void test_read_empty(void);
void test_read_example1(void);
void test_read_example3(void);
void test_read_example4(void);
void test_read_example5(void);
void test_read_damaged1(void);
void test_read_damaged2(void);
void test_read_damaged3(void);

int main(void)
{
    test_read_empty();
    test_read_example1();
    test_read_example3();
    test_read_example4();
    test_read_example5();
    test_read_damaged1();
    test_read_damaged2();
    test_read_damaged3();
    return hope_status();
}

void test_read_empty(void)
{
    FILE *fd = fopen(ASSETS "/empty.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        i++;
    }
    EQ(i, 0);
    EQ(rc, GFF_ENDFILE);

    fclose(fd);
}

static enum gff_elem_type ex1_type[] = {GFF_VERSION, GFF_FEATURE, GFF_FEATURE,
                                        GFF_FEATURE};

static struct gff_feature ex1_feature[2] = {
    (struct gff_feature){"AE014075.1:190-252|dna", "iseq", ".", "1", "63",
                         "0.0", "+", ".",
                         "ID=item1;Target_alph=dna;Profile_name=Leader_Thr;"
                         "Profile_alph=dna;Profile_"
                         "acc=PF08254.12;Window=0;Bias=17.5;E-value=2.9e-14;"
                         "Epsilon=0.01;Score=38."
                         "8"},
    (struct gff_feature){"AE014075.1:534-908|dna", "iseq", ".", "1", "306",
                         "0.0", "+", ".",
                         "ID=item2;Target_alph=dna;Profile_name=Y1_Tnp;Profile_"
                         "alph=dna;Profile_acc=PF01797.17;Window=0;Bias=0.0;E-"
                         "value=1.7e-29;Epsilon=0.01;Score=88.6"}};

static enum gff_elem_type ex4_type[25] = {
    GFF_VERSION, GFF_REGION,  GFF_FEATURE, GFF_FEATURE, GFF_FEATURE,
    GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE,
    GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE,
    GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE,
    GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE, GFF_FEATURE,
};

static struct gff_region ex4_region[1] = {
    GFF_REGION_INIT(ex4_region[0], "ctg123", "1", "1497228")};

static struct gff_feature ex4_feature[23] = {
    (struct gff_feature){"ctg123", ".", "gene", "1000", "9000", ".", "+", ".",
                         "ID=gene00001;Name=EDEN"},
    (struct gff_feature){"ctg123", ".", "TF_binding_site", "1000", "1012", ".",
                         "+", ".", "ID=tfbs00001;Parent=gene00001"},
    (struct gff_feature){"ctg123", ".", "mRNA", "1050", "9000", ".", "+", ".",
                         "ID=mRNA00001;Parent=gene00001;Name=EDEN.1"},
    (struct gff_feature){"ctg123", ".", "mRNA", "1050", "9000", ".", "+", ".",
                         "ID=mRNA00002;Parent=gene00001;Name=EDEN.2"},
    (struct gff_feature){"ctg123", ".", "mRNA", "1300", "9000", ".", "+", ".",
                         "ID=mRNA00003;Parent=gene00001;Name=EDEN.3"},
    (struct gff_feature){"ctg123", ".", "exon", "1300", "1500", ".", "+", ".",
                         "ID=exon00001;Parent=mRNA00003"},
    (struct gff_feature){"ctg123", ".", "exon", "1050", "1500", ".", "+", ".",
                         "ID=exon00002;Parent=mRNA00001,mRNA00002"},
    (struct gff_feature){"ctg123", ".", "exon", "3000", "3902", ".", "+", ".",
                         "ID=exon00003;Parent=mRNA00001,mRNA00003"},
    (struct gff_feature){"ctg123", ".", "exon", "5000", "5500", ".", "+", ".",
                         "ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003"},
    (struct gff_feature){"ctg123", ".", "exon", "7000", "9000", ".", "+", ".",
                         "ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003"},
    (struct gff_feature){"ctg123", ".", "CDS", "1201", "1500", ".", "+", "0",
                         "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"},
    (struct gff_feature){"ctg123", ".", "CDS", "3000", "3902", ".", "+", "0",
                         "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"},
    (struct gff_feature){"ctg123", ".", "CDS", "5000", "5500", ".", "+", "0",
                         "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"},
    (struct gff_feature){"ctg123", ".", "CDS", "7000", "7600", ".", "+", "0",
                         "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1"},
    (struct gff_feature){"ctg123", ".", "CDS", "1201", "1500", ".", "+", "0",
                         "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"},
    (struct gff_feature){"ctg123", ".", "CDS", "5000", "5500", ".", "+", "0",
                         "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"},
    (struct gff_feature){"ctg123", ".", "CDS", "7000", "7600", ".", "+", "0",
                         "ID=cds00002;Parent=mRNA00002;Name=edenprotein.2"},
    (struct gff_feature){"ctg123", ".", "CDS", "3301", "3902", ".", "+", "0",
                         "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"},
    (struct gff_feature){"ctg123", ".", "CDS", "5000", "5500", ".", "+", "1",
                         "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"},
    (struct gff_feature){"ctg123", ".", "CDS", "7000", "7600", ".", "+", "1",
                         "ID=cds00003;Parent=mRNA00003;Name=edenprotein.3"},
    (struct gff_feature){"ctg123", ".", "CDS", "3391", "3902", ".", "+", "0",
                         "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"},
    (struct gff_feature){"ctg123", ".", "CDS", "5000", "5500", ".", "+", "1",
                         "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"},
    (struct gff_feature){"ctg123", ".", "CDS", "7000", "7600", ".", "+", "1",
                         "ID=cds00004;Parent=mRNA00003;Name=edenprotein.4"}};

static enum gff_elem_type ex5_type[3] = {
    GFF_VERSION,
    GFF_FEATURE,
    GFF_FEATURE,
};

static struct gff_feature ex5_feature[2] = {
    (struct gff_feature){"J02448", "GenBank", "region", "1", "6407", ".", "+",
                         ".", "ID=J02448;Name=J02448;Is_circular=true;"},
    (struct gff_feature){"J02448", "GenBank", "CDS", "6006", "7238", ".", "+",
                         "0", "ID=geneII;Name=II;Note=protein II;"}};

static void eq_feat(struct gff_feature const *actual,
                    struct gff_feature const *desired)
{
    EQ(actual->seqid, desired->seqid);
    EQ(actual->source, desired->source);
    EQ(actual->type, desired->type);
    EQ(actual->start, desired->start);
    EQ(actual->end, desired->end);
    EQ(actual->score, desired->score);
    EQ(actual->strand, desired->strand);
    EQ(actual->phase, desired->phase);
    EQ(actual->attrs, desired->attrs);
}

static void eq_region(struct gff_region const *actual,
                      struct gff_region const *desired)
{
    EQ(actual->name, desired->name);
    EQ(actual->start, desired->start);
    EQ(actual->end, desired->end);
}

void test_read_example1(void)
{
    FILE *fd = fopen(ASSETS "/example1.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    unsigned i_feat = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        EQ(gff.elem.type, ex1_type[i]);

        if (ex1_type[i] == GFF_VERSION)
            EQ(gff.elem.version, "3");
        else if (ex1_type[i] == GFF_REGION)
        {
        }
        else if (ex1_type[i] == GFF_FEATURE)
            eq_feat(&gff.elem.feature, &ex1_feature[i_feat++]);
        else
            COND(false);

        i++;
    }
    EQ(i, 3);
    EQ(rc, GFF_ENDFILE);

    fclose(fd);
}

void test_read_example3(void)
{
    FILE *fd = fopen(ASSETS "/example3.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    enum gff_rc rc = gff_read(&gff);
    EQ(gff.elem.type, GFF_VERSION);
    EQ(gff.elem.version, "3");
    rc = gff_read(&gff);
    EQ(rc, GFF_ENDFILE);

    fclose(fd);
}

void test_read_example4(void)
{
    FILE *fd = fopen(ASSETS "/example4.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    unsigned i_region = 0;
    unsigned i_feat = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        EQ(gff.elem.type, ex4_type[i]);

        if (gff.elem.type == GFF_VERSION)
            EQ(gff.elem.version, "3.1.26");
        else if (gff.elem.type == GFF_REGION)
            eq_region(&gff.elem.region, &ex4_region[i_region++]);
        else if (gff.elem.type == GFF_FEATURE)
            eq_feat(&gff.elem.feature, &ex4_feature[i_feat++]);
        else
            COND(false);

        i++;
    }
    EQ(i, 25);
    EQ(rc, GFF_ENDFILE);

    fclose(fd);
}

void test_read_example5(void)
{
    FILE *fd = fopen(ASSETS "/example5.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    unsigned i_feat = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        EQ(gff.elem.type, ex5_type[i]);

        if (gff.elem.type == GFF_VERSION)
            EQ(gff.elem.version, "3.1.26");
        else if (gff.elem.type == GFF_FEATURE)
            eq_feat(&gff.elem.feature, &ex5_feature[i_feat++]);
        else
            COND(false);

        i++;
    }
    EQ(i, 3);
    EQ(rc, GFF_ENDFILE);

    fclose(fd);
}

void test_read_damaged1(void)
{
    FILE *fd = fopen(ASSETS "/damaged1.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        i++;
    }
    EQ(i, 0);
    EQ(rc, GFF_PARSEERROR);
    EQ(gff.error, "Parse error: unexpected token: line 1");
    gff_clearerr(&gff);
    EQ(gff.error, "");

    fclose(fd);
}

void test_read_damaged2(void)
{
    FILE *fd = fopen(ASSETS "/damaged2.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        i++;
    }
    EQ(i, 0);
    EQ(rc, GFF_PARSEERROR);
    EQ(gff.error, "Parse error: unexpected id: line 2");
    gff_clearerr(&gff);
    EQ(gff.error, "");

    fclose(fd);
}

void test_read_damaged3(void)
{
    FILE *fd = fopen(ASSETS "/damaged3.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        i++;
    }
    EQ(i, 0);
    EQ(rc, GFF_PARSEERROR);
    EQ(gff.error, "Parse error: unexpected token: line 4");
    gff_clearerr(&gff);
    EQ(gff.error, "");

    fclose(fd);
}
