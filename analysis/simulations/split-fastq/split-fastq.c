
#include <stdio.h>
#include <string.h>

#include "common.h"
#include "parse.h"

int main(int argc, char* argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: split-fastq READ_LENGTH reads.fastq\n");
        return 1;
    }

    unsigned long read_length = strtoul(argv[1], NULL, 10);

    FILE* f = fopen_or_die(argv[2], "r");

    char* c = strstr(argv[2], ".fastq");
    if (c) *c = '\0';

    char *fn1, *fn2;
    asprintf(&fn1, "%s_1.fastq", argv[2]);
    asprintf(&fn2, "%s_2.fastq", argv[2]);

    FILE* fout1 = fopen_or_die(fn1, "w");
    FILE* fout2 = fopen_or_die(fn2, "w");

    fprintf(stderr, "writing to %s and %s.\n", fn1, fn2);

    fastq_t* fqf = fastq_create(f);
    seq_t* seq = seq_create();

    while (fastq_read(fqf, seq)) {
        if (seq->seq.n != read_length) continue;

        if (seq->id1.s[seq->id1.n - 1] == '2') {
            seq->id1.s[seq->id1.n - 4] = '\0'; // HACk
            fastq_print(fout2, seq);
        }
        else {
            seq->id1.s[seq->id1.n - 4] = '\0'; // HACK
            fastq_print(fout1, seq);
        }
    }

    seq_free(seq);
    fastq_free(fqf);
    fclose(f);
    fclose(fout1);
    fclose(fout2);

    return 0;
}
