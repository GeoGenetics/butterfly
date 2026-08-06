#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khashl.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 134217728U)

KHASHL_MAP_INIT(static,
                readset_t,
                readset,
                uint64_t,
                uint8_t,
                kh_hash_uint64,
                kh_eq_generic)

static void usage(void)
{
  fprintf(stderr, "bam_cmp <(samtools view <file1> | cut -f1) <(samtools view <file2> | cut -f1)\n");
  fprintf(stderr, "Compares two bam files and outputs file1 exclusive, intersection and file2 exclusive counts\n");
  fprintf(stderr, "\tAdditionally, prints readids of intersection.\n");
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
  int ret = -1;
  gzFile fp1 = NULL, fp2 = NULL, fp_intersection = NULL, fp_exclusive1 = NULL, fp_exclusive2 = NULL;
  kstream_t *ks1 = NULL, *ks2 = NULL;
  kstring_t kstr1 = {0}, kstr2 = {0};
  uint8_t _printexcl2 = 0;
  if (argc < 3) goto error;
  if (argc > 3) _printexcl2 = 1;
  readset_t *readset1 = readset_init();
  readset_t *readset2 = readset_init();
  readset_t *readseti = readset_init();
  fp1 = gzopen(argv[1], "r");
  fp2 = gzopen(argv[2], "r");
  fp_intersection = gzopen("intersection.txt.gz", "w9");
  //fp_exclusive1   = gzopen("exclusive1.txt.gz", "w9");
  fp_exclusive2   = gzopen("exclusive2.txt.gz", "w9");
  if (!fp1 || !fp2 || !fp_intersection || !fp_exclusive2) {
    fprintf(stderr, "Failed to open files\n");
    goto exit;
  }
  ks1 = ks_init(fp1);
  ks2 = ks_init(fp2);
  if (!ks1 || !ks2) {
    fprintf(stderr, "Failed to initialize stream for %s\n", ks1 ? argv[2] : argv[1]);
    goto exit;
  }
  //stream file 1 and create hash set
  while ( (ks_getuntil(ks1, '\n', &kstr1, 0)) >= 0 ) {
    if (kstr1.l == 0) break;
    int absent;
    khint_t k = readset_put(readset1, hashfun_n(kstr1.s, strlen(kstr1.s)), &absent);
    if (absent) kh_val(readset1, k) = 0;
  }
  uint64_t n_intersection = 0, n_file2 = 0;
  //Stream file 2 and check for intersection and exclusive reads
  while ( (ks_getuntil(ks2, '\n', &kstr2, 0)) >= 0 ) {
    if (kstr2.l == 0) break;
    int absent;
    khint_t k = readset_put(readset2, hashfun_n(kstr2.s, strlen(kstr2.s)), &absent);
    if (absent) kh_val(readset2, k) = 0;
    k = readset_get(readset1, hashfun_n(kstr2.s, strlen(kstr2.s)) );
    if (k != kh_end(readset1)) { //Intersection
      khint_t ki = readset_put(readseti, hashfun_n(kstr2.s, strlen(kstr2.s)), &absent);
      if (absent) {
        gzprintf(fp_intersection, "%s\n", kstr2.s);
        kh_val(readseti, ki) = 0;
        n_intersection++;
      }
      kh_val(readset1, k) = 1;
    }
    else { //File2 exclusive
      k = readset_get(readset2, hashfun_n(kstr2.s, strlen(kstr2.s)) );
      if (!kh_val(readset2, k)) {
        n_file2++;
        if ( _printexcl2 ) gzprintf(fp_exclusive2, "%s\n", kstr2.s);
      }
      kh_val(readset2, k) = 1;
    }
  }
  khint_t k;
  uint64_t n_file1 = 0;
  kh_foreach(readset1, k) { //File1 exclusive count
    if (kh_val(readset1, k) == 0) n_file1++;
  }
  fprintf(stdout, "%u\t%lu\t%lu\t%u\t%lu\n", kh_size(readset1), n_file1, n_intersection, kh_size(readset2), n_file2);
  ret = 0;
  exit:
    if (readset1) readset_destroy(readset1);
    if (readset2) readset_destroy(readset2);
    if (readseti) readset_destroy(readseti);
    if (ks1) ks_destroy(ks1);
    if (ks2) ks_destroy(ks2);
    if (fp1) gzclose(fp1);
    if (fp2) gzclose(fp2);
    if (fp_intersection) gzclose(fp_intersection);
    if (fp_exclusive1) gzclose(fp_exclusive1);
    if (fp_exclusive2) gzclose(fp_exclusive2);
    free(kstr1.s);
    free(kstr2.s);
    return ret;
  error:
    usage();
    return ret;
}
