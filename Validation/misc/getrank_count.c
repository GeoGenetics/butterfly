#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khashl.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 134217728U)

KHASHL_MAP_INIT(static,
                rankcount_t,
                rankcount,
                char *,
                uint64_t,
                kh_hash_str,
                kh_eq_str)

static void usage(void)
{
  fprintf(stderr, "Usage: cut -f1,6 <file1> | getrank_count /dev/stdin\n");
  fprintf(stderr, "Counts the occurrences of each rank in the input file.\n");
  fprintf(stderr, "Arguments:\n");
  fprintf(stderr, "  <file1> : Path to the first input file (gzipped).\n");
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

int main(int argc, char *argv[])
{
  int ret = -1;
  if (argc < 2) goto error;
  rankcount_t *rankmap = rankcount_init();
  gzFile fp = gzopen(argv[1], "r");
  if (!fp) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
    goto exit;
  }
  kstream_t *ks1 = ks_init(fp);
  if (!ks1) {
    fprintf(stderr, "Failed to initialize stream for %s\n", argv[1]);
    goto exit;
  }
  kstring_t kstr1 = {0};
  char *tok;
  ks_getuntil(ks1, '\n', &kstr1, 0);
  ks_getuntil(ks1, '\n', &kstr1, 0);
  while ( (ks_getuntil(ks1, '\n', &kstr1, 0)) >= 0 ) {
    if (kstr1.l == 0) break;
    tok = strtok(kstr1.s, "\t");
    tok = strtok(NULL, "\t");
    char *rank =  getrank(tok);
    khint_t k;
    int absent;
    k = rankcount_put(rankmap, rank, &absent);
    if (absent) {
      kh_key(rankmap, k) = strdup(rank);
      if (!kh_key(rankmap, k)) {
        fprintf(stderr, "Out of memory while copying rank string\n");
        goto exit;
      }
      kh_val(rankmap, k) = 0;
    }
    kh_val(rankmap, k)++;
    fflush(stdout);
  }
  goto exit;
  exit:
    ret = 0;
    if (rankmap) {
      khint_t k;
      kh_foreach(rankmap, k) {
        fprintf(stderr, "%s\t%lu\n", kh_key(rankmap, k), kh_val(rankmap, k));
        free((char *)kh_key(rankmap, k));
      }
      rankcount_destroy(rankmap);
    }
    if (ks1) ks_destroy(ks1);
    if (fp) gzclose(fp);
    if (ret) usage();
    return ret;
  error:
    usage();
    return ret;
}
