#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "khashl.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 134217728U)

KHASHL_MAP_INIT(static,
                rankcount_t,
                rankcount,
                int,
                uint64_t,
                kh_hash_uint32,
                kh_eq_generic)

/**
 * Parses the lcarank from the end of the line in-place.
 * NOTE: Modifies the input string.
 * * @param line Pointer to the mutable line buffer (e.g., from fgets/getline)
 * @return Pointer to the start of the rank inside the buffer, or NULL if malformed.
 */
char* getrank(char *line) {
  if (!line) return NULL;
  int end = strlen(line) - 1;
  // 1. Strip trailing newlines or carriage returns
  while (end >= 0 && (line[end] == '\n' || line[end] == '\r')) {
    end--;
  }
  // 2. The last character must be the closing quote of the rank
  if (end < 0 || line[end] != '"') return NULL;
  // Replace closing quote with null terminator to cut the string here
  line[end] = '\0'; 
  // 3. Walk backward to find the opening quote of the rank
  int open = end - 1;
  while (open >= 0 && line[open] != '"') {
    open--;
  }
  // 4. Validation: must find the quote, and it must be preceded by a colon ':'
  if (open <= 0 || line[open - 1] != ':') {
    return NULL; 
  }
  // Return pointer to the character right after the opening quote
  return &line[open + 1];
}

int main(int argc, char *argv[])
{
	int ret = -1;
	if (argc < 2) goto error;
  rankcount_t *rankmap = rankcount_init();
  gzFile fp = gzopen(argv[1], "r"); 
	kstream_t *ks1 = ks_init(fp);
  kstring_t kstr1 = {0};
	char *tok;
	while ( (ks_getuntil(ks1, '\n', &kstr1, 0)) >= 0 ) {
    if (kstr1.l == 0) break;
    tok = strtok(kstr1.s, "\t");
    tok = strtok(NULL, "\t");
		fprintf(stdout, "%s\n", getrank(tok) );
	  fflush(stdout);
	}
	goto exit;

 

	exit:
    ret = 0;
		if (rankmap) rankcount_destroy(rankmap);
		if (ks1) ks_destroy(ks1);
		if (fp) gzclose(fp);
		return ret;
	error:
	  return ret;
}