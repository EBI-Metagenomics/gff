#include "gff/gff.h"
#include "hope/hope.h"

void test_read_empty(void);
void test_read_example1(void);
void test_read_damaged1(void);
void test_read_damaged2(void);
void test_read_damaged3(void);

int main(void)
{
    test_read_empty();
    test_read_example1();
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

static char *mix_id[] = {"LCBO", "MCHU", "gi|5524211|gb|AAD44166.1|",
                         "gi|5524211|gb|AAD44166.1|"};

static char *mix_desc[] = {
    "- Prolactin precursor - Bovine",
    "- Calmodulin - Human, rabbit, bovine, rat, and chicken",
    "cytochrome b [Elephas maximus maximus]",
    "cytochrome b [Elephas maximus maximus]"};

static char *mix_seq[] = {
    "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYI"
    "HDLSSEMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLL"
    "RSWNDPLYHLVTEVRGMKGAPDAILSRAIEIEEENKRLLEGMEMIFGQVIPGAKETEPY"
    "PVWSGLPSLQTKDEDARYSAFYNLLHCLRRDSSKIDTYLKLLNCRIIYNNNC",
    "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDAD"
    "GNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLT"
    "DEEVDEMIREADIDGDGQVNYEEFVQMMTAK*",
    "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV"
    "EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG"
    "LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL"
    "GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX"
    "IENY",
    "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGT"};

void test_read_example1(void)
{
    FILE *fd = fopen(ASSETS "/example1.gff", "r");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_READ);

    unsigned i = 0;
    enum gff_rc rc = GFF_SUCCESS;
    while (!(rc = gff_read(&gff)))
    {
        EQ(gff.target.id, mix_id[i]);
        EQ(gff.target.desc, mix_desc[i]);
        EQ(gff.target.seq, mix_seq[i]);
        i++;
    }
    EQ(i, 4);
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
