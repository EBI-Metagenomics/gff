#ifndef TOK_H
#define TOK_H

#include "gff/rc.h"
#include <stdio.h>

struct gff_tok;

enum tok_id
{
    TOK_NL,
    TOK_COMMENT,
    TOK_PRAGMA,
    TOK_VERSION,
    TOK_REGION,
    TOK_FASTA,
    TOK_WORD,
    TOK_EOF,
};

void tok_init(struct gff_tok *tok, char *error);
enum gff_rc tok_next(struct gff_tok *tok, FILE *restrict fd);

#endif
