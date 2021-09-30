#ifndef GFF_STRLCPY_H
#define GFF_STRLCPY_H

#include "gff/export.h"
#include <stdbool.h>
#include <stddef.h>

GFF_API size_t gff_strlcpy(char *dst, const char *src, size_t len);

#endif
