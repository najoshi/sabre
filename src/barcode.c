#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

int strncmp_with_mismatch (const char *s1, const char *s2, register size_t n, register size_t mismatch) {

	register unsigned char u1, u2;
	int cnt=0;

	while (n-- > 0) {
		u1 = (unsigned char) *s1++;
		u2 = (unsigned char) *s2++;

		if (u1 != u2) {
			cnt++;
			if (cnt > mismatch) return u1 - u2;
		}

		if (u1 == '\0') return 0;
	}

	return 0;
}
