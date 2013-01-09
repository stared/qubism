// Routines for parsing 
#ifndef PARSING_HEADER
#define PARSING_HEADER

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

char* get_until(const char*, const char[]);
char* get_token(const char*, const char[], int);
bool is_here(const char *, const char *);
bool is_there(const char *, const char *);
int get_int(const char *, const char [], int);
double get_double(const char *, const char [], int);
int find_first(const char *, char );
char* part(const char *, int , int );
int get_order(char*& command, char*& arg, char*& program);
char *compose(char*,char*);
void right_add(char*&, char*);
void left_add(char*&, char*);

void add_char(char **, size_t *, int , char );
ssize_t get_line(char **, size_t *, FILE *);
int get_next_line(char **, FILE *);

#endif
