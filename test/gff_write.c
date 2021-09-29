#include "gff/gff.h"
#include "hope/hope.h"

#if 0
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

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))
#endif

int main(void)
{
    FILE *fd = fopen(TMPDIR "/example1.gff", "w");
    NOTNULL(fd);

    struct gff gff;
    gff_init(&gff, fd, GFF_WRITE);

    gff_set_version(&gff);
    enum gff_rc rc = gff_write(&gff);
    EQ(rc, GFF_SUCCESS);

#if 0
    for (unsigned i = 0; i < ARRAY_SIZE(mix_id); ++i)
    {
        gff_write(&gff, gff_target(mix_id[i], mix_desc[i], mix_seq[i]), 60);
    }
#endif

    fclose(fd);

#if 0
    FILE *actual = fopen(TMPDIR "/mix.gff", "r");
    FILE *desired = fopen(ASSETS "/desired_mix.gff", "r");
    NOTNULL(actual);
    NOTNULL(desired);
    EQ(actual, desired);
    fclose(desired);
    fclose(actual);
#endif

    return hope_status();
}
