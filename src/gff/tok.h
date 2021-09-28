#ifndef GFF_TOK_H
#define GFF_TOK_H

#include <stdbool.h>

#define GFF_TOK_LINE_MAX 512

struct gff_tok
{
    unsigned id;
    char const *value;
    struct
    {
        char data[GFF_TOK_LINE_MAX];
        unsigned number;
        bool consumed;
        char *ctx;
    } line;
    char *error;
};

#endif
