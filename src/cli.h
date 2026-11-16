#ifndef CLI_H
#define CLI_H

#include <stddef.h>

typedef enum { OPT_INT, OPT_DOUBLE, OPT_STRING } opt_type_t;

typedef struct {
	const char	*key;
	opt_type_t	type;
	size_t	offset;
	double	min, max;
	const char	*help;
} opt_t;



#define OPT(S, field, lo, hi, help) \
	{ #field, \
	  _Generic(((S *)0)->field, int: OPT_INT, double: OPT_DOUBLE, const char *: OPT_STRING), \
	  offsetof(S, field), (lo), (hi), (help) }

int cli_parse(int argc, char **argv, void *cfg, const opt_t *opts, size_t nopts);
#endif
