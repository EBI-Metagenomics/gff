#ifndef ERROR_H
#define ERROR_H

#include "gff/rc.h"

enum gff_rc error_io(char *dst, int errnum);
enum gff_rc error_runtime(char *dst, char const *fmt, ...);
enum gff_rc error_parse(char *dst, unsigned line, char const *fmt, ...);

#endif
