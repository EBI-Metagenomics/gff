#ifndef ERROR_H
#define ERROR_H

#include "gff/rc.h"

enum gff_rc error(enum gff_rc rc, char *dst, char const *msg);
enum gff_rc error_io(char *dst, int errnum);
enum gff_rc error_parse(char *dst, unsigned line, char const *msg);

#endif
