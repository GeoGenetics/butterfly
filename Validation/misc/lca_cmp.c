#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khashl.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 134217728U)

KHASHL_SET_INIT(static,
                readset_t,
                readset,
                uint64_t,
                kh_hash_uint64,
                kh_eq_generic)

static void usage(void)
{
  fprintf(stderr, "Usage: cut -f1,6 <file1> | lca_cmp /dev/stdin <taxid> <file2>\n");
  fprintf(stderr, "Compares two lca files and outputs file2 reads whose assigned taxid in file1 is <taxid>\n");
  fprintf(stderr, "Arguments:\n");
  fprintf(stderr, "  <file1> : Path to the first input file (gzipped).\n");
  fprintf(stderr, "  <taxid> : Taxonomic ID to filter reads from file1 Use 0 to include all.\n");
  fprintf(stderr, "  <file2> : Path to the second input file (gzipped).\n");
}

static uint64_t hashfun_n(const char *s, size_t len)
{
  const unsigned char *p = (const unsigned char *)s;
  uint64_t h = 1469598103934665603ULL;
  while (len--) {
    h ^= (uint64_t)*p++;
    h *= 1099511628211ULL;
  }
  return h;
}

/**
 * Parses the lcarank from the end of the line in-place.
 * NOTE: Modifies the input string.
 * * @param line Pointer to the mutable line buffer (e.g., from fgets/getline)
 * @return Pointer to the start of the rank inside the buffer, or NULL if malformed.
 */
char* getrank(char *line)
{
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

/**
 * Extracts the taxid from the line and returns it as an unsigned integer.
 * * @param line The input line string.
 * @return The taxid as an unsigned int, or 0 if the line is malformed.
 */
unsigned int gettaxid(const char *line)
{
    if (!line) return 0;
    // 1. Find the tab character that separates Column 1 and Column 2
    const char *tab = strchr(line, '\t');
    if (!tab) return 0; // Malformed line: no tab found
    // 2. The taxid starts exactly one character after the tab
    const char *taxid_start = tab + 1;
    // 3. Quick validation: ensure the first character is actually a digit
    if (*taxid_start < '0' || *taxid_start > '9') {
        return 0;
    }
    // 4. Convert to unsigned long and cast to unsigned int.
    // strtoul will gracefully stop reading when it hits the ':'
    return (unsigned int)strtoul(taxid_start, NULL, 10);
}

int main(int argc, char *argv[])
{
	int ret = -1, VERBOSE = 0;
	gzFile fp = NULL, fp2 = NULL;
	kstream_t *ks1 = NULL, *ks2 = NULL;
  kstring_t kstr1 = {0}, kstr2 = {0};
	if (argc < 4) goto error;
  if (argc > 4) VERBOSE = 1;
  unsigned int taxid = (unsigned int)strtoul(argv[2], NULL, 10);
	readset_t *readset = readset_init();
	fp  = gzopen(argv[1], "r");
  fp2 = gzopen(argv[3], "r");

  if (!fp || !fp2) {
    fprintf(stderr, "Failed to open %s\n", fp ? argv[2] : argv[1]);
    goto exit;
  }
	ks1 = ks_init(fp);
  ks2 = ks_init(fp2);
  if (!ks1 || !ks2) {
    fprintf(stderr, "Failed to initialize stream for %s\n", ks1 ? argv[2] : argv[1]);
    goto exit;
  }

	unsigned int count = 0;
  while ( (ks_getuntil(ks1, '\n', &kstr1, 0)) >= 0 ) {
    char *tab;
    if (kstr1.l == 0) break;
		count++;
    if ( (taxid == gettaxid(kstr1.s) ) || (taxid == 0) ) {
      tab = strchr(kstr1.s, '\t');
      if (!tab) continue;
			int absent;
			readset_put(readset, hashfun_n(kstr1.s, (size_t)(tab - kstr1.s)), &absent);
			if (absent < 0) {
				fprintf(stderr, "Failed to insert hashed read id into hash set\n");
				goto exit;
			}
		  if (!(count % 1000000)  && VERBOSE) {
				fprintf(stderr, "%u reads processed so far %u assigned\n", count, kh_size(readset));
			  fflush(stderr);
			}
		}
	}
  count = 0;
  fprintf(stderr, "Now checking %s for reads assigned in %s\n", argv[3], argv[1]);
  fflush(stderr);
  while ( (ks_getuntil(ks2, '\n', &kstr2, 0)) >= 0 ) {
    char *tab = strchr(kstr2.s, '\t');
    khint_t k;
    if (kstr2.l == 0) break;
    if (!tab) continue;
    k = readset_get(readset, hashfun_n(kstr2.s, (size_t)(tab - kstr2.s)));
    if (k != kh_end(readset)) {
      count++;
      fputs(kstr2.s, stdout);
      fputc('\n', stdout);
    }
  }
  fprintf(stderr, "%u reads in %s\n", kh_size(readset), argv[1]);
  fprintf(stderr, "\tof which %u (%.3f) were found in %s\n", count, (double)count / kh_size(readset) * 100, argv[3]);
  ret = 0;
	exit:
		if (readset) readset_destroy(readset);
		if (ks1) ks_destroy(ks1);
		if (ks2) ks_destroy(ks2);
		if (fp) gzclose(fp);
		if (fp2) gzclose(fp2);
    free(kstr1.s);
    free(kstr2.s);
		return ret;
	error:
	  usage();
	  return ret;
}
