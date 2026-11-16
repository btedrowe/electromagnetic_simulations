#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

static void print_usage(const char *prog, const void *cfg, const opt_t *opts, size_t n) {
	fprintf(stderr, "Usage: %s [key=value ...]\n", prog);
	for (size_t i = 0; i < n; i++) {
		const void *p = (const char *)cfg + opts[i].offset;
		fprintf(stderr, " %-8s %s (default: ", opts[i].key, opts[i].help);
		
		switch (opts[i].type) {
		case OPT_INT:
			fprintf(stderr, "%d", *(const int *)p);
			break;
		case OPT_DOUBLE:
			fprintf(stderr, "%g", *(const double *)p);
			break;
		case OPT_STRING:
			fprintf(stderr, "%s", *(const char *const *)p ? *(const char *const *)p : "(none)");
			break;
		}
		
		fprintf(stderr, ")\n");
	}
}

int cli_parse(int argc, char **argv, void *cfg, const opt_t *opts, size_t nopts) {
	for (int i = 1; i < argc; i++) {
		const char *arg = argv[i];
		
		if (!strcmp(arg, "-h") || !strcmp(arg, "--help")) {
			print_usage(argv[0], cfg, opts, nopts);
			return 2;
		}
		
		const char *eq = strchr(arg, '=');
		if (!eq || eq == arg) {
			fprintf(stderr, "bad argument '%s' (expected key=value)\n", arg);
			return 1;
		}
		size_t klen = (size_t)(eq - arg);
		const char *val = eq + 1;
		
		const opt_t *o = NULL;
		for (size_t j = 0; j < nopts; j++) {
			if (strlen(opts[j].key) == klen && 
			    !strncmp(opts[j].key, arg, klen)) {
				o = &opts[j];
				break;
			}
		}
		if (!o) {
			fprintf(stderr, "unknown key '%.*s'\n", (int)klen, arg);
			return 1;
		}
		
		void *dst = (char *)cfg + o->offset;
		char *end;
		errno = 0;
		
		switch (o->type) {
		case OPT_INT: {
			long v = strtol(val, &end, 10);
			if (end == val || *end || errno ||
			    v < o->min || v > o->max) {
				fprintf(
					stderr, 
					"%s=%s: expected an integer in [%g, %g]\n", 
					o->key,
					val,
					o->min,
					o->max
				);
				return 1;
			}
			*(int *)dst = (int)v;
			break;
		}
		case OPT_DOUBLE: {
			double v = strtod(val, &end);
			if (end == val || *end || errno || 
			    !isfinite(v) || v < o->min || v > o->max) {
				fprintf(
					stderr, 
					"%s=%s: expected a number in [%g, %g]\n",
					o->key, val, o->min, o->max
				);
				return 1;
			}
			*(double *)dst = v;
			break;
			
		}
		case OPT_STRING:
			if (*val == '\0') {
				fprintf(stderr, "%s=%s: expected non-empty value", o->key, val);
				return 1;
			}
			*(const char **)dst = val;
			break;
		}
	}
	
	return 0;
}
