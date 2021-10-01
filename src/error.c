#include "error.h"
#include "gff/error.h"
#include "unused.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

static char prefix[][20] = {
    "", "", "IO error:", "Runtime error:", "Parse error:", "Illegal argument:"};

enum gff_rc error(enum gff_rc rc, char *dst, char const *msg)
{
    int n = snprintf(dst, GFF_ERROR_SIZE, "%s %s", prefix[rc], msg);
    assert(0 < n && n < GFF_ERROR_SIZE);
    unused(n);
    return rc;
}

enum gff_rc error_io(char *dst, int errnum)
{
    char errstr[32] = {0};
    int rc = strerror_r(errnum, errstr, GFF_ERROR_SIZE);
    assert(!rc);
    unused(rc);
    int n = snprintf(dst, GFF_ERROR_SIZE, "%s %s", prefix[GFF_IOERROR], errstr);
    assert(0 < n && n < GFF_ERROR_SIZE);
    unused(n);
    return GFF_IOERROR;
}

enum gff_rc error_parse(char *dst, unsigned line, char const *msg)
{
    int n = snprintf(dst, GFF_ERROR_SIZE, "%s %s: line %d",
                     prefix[GFF_PARSEERROR], msg, line);
    assert(0 < n && n < GFF_ERROR_SIZE);
    unused(n);
    return GFF_PARSEERROR;
}
