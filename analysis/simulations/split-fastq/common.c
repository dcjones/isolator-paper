
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>

#include "common.h"

#ifndef O_BINARY
#define O_BINARY 0
#endif


void or_die(int b, const char* msg)
{
    if (b == 0) {
        fputs(msg, stderr);
        exit(1);
    }
}


void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(1);
    }
    return p;
}


void* realloc_or_die(void* ptr, size_t n)
{
    void* p = realloc(ptr, n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(1);
    }
    return p;
}


FILE* fopen_or_die(const char* path, const char* mode)
{
    FILE* f = fopen(path, mode);
    if (f == NULL) {
        fprintf(stderr, "Can not open file %s with mode %s.\n", path, mode);
        exit(1);
    }
    return f;
}


/* Open a file for writing, creating it if it doesn't exist, and complaining if
 * it does. */
FILE* open_without_clobber(const char* filename)
{
    int fd = open(filename, O_WRONLY | O_CREAT | O_BINARY | O_EXCL,
                  S_IRUSR | S_IWUSR);

    if (fd == -1) {
        if (errno == EEXIST) {
            fprintf(stderr, "Refusing to overwrite %s.\n", filename);
            exit(EXIT_FAILURE);
        }
        else {
            fprintf(stderr, "Cannot open %s for writing.\n", filename);
            exit(EXIT_FAILURE);
        }
    }

    FILE* f = fdopen(fd, "wb");
    if (f == NULL) {
        fprintf(stderr, "Cannot open %s for writing.\n", filename);
        exit(EXIT_FAILURE);
    }

    return f;
}


