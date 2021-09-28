#include "gff/tok.h"
#include "error.h"
#include "gff/error.h"
#include "tok.h"
#include <string.h>

#define DELIM " \t\r"

static void add_space_before_newline(char line[GFF_TOK_LINE_MAX]);
static enum gff_rc next_line(FILE *restrict fd, char error[GFF_ERROR_SIZE],
                             char line[GFF_TOK_LINE_MAX]);

void tok_init(struct gff_tok *tok, char *error)
{
    tok->id = TOK_NL;
    tok->value = tok->line.data;
    memset(tok->line.data, '\0', GFF_TOK_LINE_MAX);
    tok->line.number = 0;
    tok->line.consumed = true;
    tok->line.ctx = NULL;
    tok->error = error;
}

enum gff_rc tok_next(struct gff_tok *tok, FILE *restrict fd)
{
    enum gff_rc rc = GFF_SUCCESS;

    if (tok->line.consumed)
    {
        if ((rc = next_line(fd, tok->error, tok->line.data)))
        {
            if (rc == GFF_ENDFILE)
            {
                tok->value = NULL;
                tok->id = TOK_EOF;
                tok->line.data[0] = '\0';
                return GFF_SUCCESS;
            }
            return rc;
        }
        tok->value = strtok_r(tok->line.data, DELIM, &tok->line.ctx);
        tok->line.number++;
    }
    else
        tok->value = strtok_r(NULL, DELIM, &tok->line.ctx);

    if (!tok->value) return GFF_PARSEERROR;

    if (!strcmp(tok->value, "\n"))
        tok->id = TOK_NL;
    else if (tok->value[0] == '#' && tok->value[1] == '#')
    {
        if (!strcmp(tok->value, "##gff-version"))
            tok->id = TOK_VERSION;
        else if (!strcmp(tok->value, "##sequence-region"))
            tok->id = TOK_REGION;
        else
            tok->id = TOK_PRAGMA;
    }
    else
        tok->id = TOK_WORD;

    tok->line.consumed = tok->id == TOK_NL;

    return GFF_SUCCESS;
}

static enum gff_rc next_line(FILE *restrict fd, char error[GFF_ERROR_SIZE],
                             char line[GFF_TOK_LINE_MAX])
{
    if (!fgets(line, GFF_TOK_LINE_MAX - 1, fd))
    {
        if (feof(fd)) return GFF_ENDFILE;

        return error_io(error, ferror(fd));
    }

    add_space_before_newline(line);
    return GFF_SUCCESS;
}

static void add_space_before_newline(char line[GFF_TOK_LINE_MAX])
{
    unsigned n = (unsigned)strlen(line);
    if (n > 0)
    {
        if (line[n - 1] == '\n')
        {
            line[n - 1] = ' ';
            line[n] = '\n';
            line[n + 1] = '\0';
        }
        else
        {
            line[n - 1] = '\n';
            line[n] = '\0';
        }
    }
}
