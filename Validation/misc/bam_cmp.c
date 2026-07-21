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
  fprintf(stderr, "bam_cmp <(samtools view <file1> | cut -f1) <(samtools view <file2> | cut -f1)\n");
  fprintf(stderr, "Compares two bam files and outputs intersection, file1 exclusive, and file2 exclusive read ids\n");
  fprintf(stderr, "Arguments:\n");
  fprintf(stderr, "  <file1> : Stream to file1 read ids.\n");
  fprintf(stderr, "  <file2> : Stream to file2 read ids.\n");
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

int main(int argc, char *argv[])
{
	int ret = -1, VERBOSE = 0;
	gzFile fp1 = NULL, fp2 = NULL;
	kstream_t *ks1 = NULL, *ks2 = NULL;
  kstring_t kstr1 = {0}, kstr2 = {0};
	if (argc < 3) goto error;
  if (argc > 3) VERBOSE = 1;
	readset_t *readset = readset_init();
	fp1 = gzopen(argv[1], "r");
  fp2 = gzopen(argv[2], "r");
  if (!fp1 || !fp2) {
    fprintf(stderr, "Failed to open %s\n", fp1 ? argv[2] : argv[1]);
    goto exit;
  }
	ks1 = ks_init(fp1);
  ks2 = ks_init(fp2);
  if (!ks1 || !ks2) {
    fprintf(stderr, "Failed to initialize stream for %s\n", ks1 ? argv[2] : argv[1]);
    goto exit;
  }
  //stream file 1
	unsigned int count = 0;
  while ( (ks_getuntil(ks1, '\n', &kstr1, 0)) >= 0 ) {
    if (kstr1.l == 0) break;
		count++;
	  fprintf(stderr, "%s\n", kstr1.s);
  }
  fprintf(stderr, "%u reads in file1\n", kh_size(readset));
  ret = 0;
	exit:
		if (readset) readset_destroy(readset);
		if (ks1) ks_destroy(ks1);
		if (ks2) ks_destroy(ks2);
		if (fp1) gzclose(fp1);
		if (fp2) gzclose(fp2);
    free(kstr1.s);
    free(kstr2.s);
		return ret;
	error:
	  usage();
	  return ret;
}
