/******************************************************************************
* Header for utils.c
******************************************************************************/
#ifndef _UTILS_H
#define _UTILS_H

/* Flags */
#define TRUE  1
#define FALSE 0

#define FAILURE     0   /* function/command failed   */
#define SUCCESS     1   /* function/command succeded */

#define MKCOPY   1
#define NOCOPY   0

#define DEFAULT       1
#define SKIP_BLOCKS   0

/* Aliases */
//#define AND  &&
//#define OR   ||
//#define EQ   ==
//#define NEQ  !=

/* Macros  */
#define is_white_space(c)  (((c)==' ') || ((c)=='\t'))
#define new_array(n, type) ((type *) mem_alloc((n)*sizeof(type),(type *)NULL))
#define mod_array(ptr, n, type) ((type *) mem_alloc((n)*sizeof(type),ptr))
//#ifndef max
//#define max(x,y) ((x) > (y) ? (x) : (y))
//#define min(x,y) ((x) < (y) ? (x) : (y))
//#endif

/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// function prototypes
///////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
void print_error(char *msg,...);

void *mem_alloc(size_t nbytes, void *ptr);

char *read_line(FILE *infile);

char *get_line(char *prompt);

char *get_line_default(char *prompt, char *default_val);

FILE  *open_file(char *fname);

void  close_file(FILE *fid);

char *create_string(int len);

char *char_to_string(char ch);

char *copy_string(char *s);

char *concat_2(char *s1, char *s2);

char *concat_3(char *s1, char *s2, char *s3);

char *concat_strings(char *str1, ...);

char *concat_int_to_string(char *str, int val, int prec);

char *dot_delim_string(char *s1);

char *sub_string(char *s, int p1, int p2);

int string_equal(char *s1, char *s2);

int string_compare(char *s1, char *s2);

int string_len(char *s);

int find_char(char ch, char *text, int start);

int find_string(char *text, char *str, int start);

int find_string_skip_blocks(char *s, char *str, int p1);

int is_open_block( char c);

int is_close_block( char c);

char *string_replace(char *str, char *fstr, char *rstr,int *num_rep);

char *to_lower(char *s);

char *to_upper(char *s);

char *num_to_string(double d);

double string_to_num(char *s);

int  string_to_int(char *s);

double *string_to_dbl_array(char *str, int *n);

double *string_to_dbl_array_fmt(char *str, int *nrow, int *ncol );

int *string_to_int_array_fmt(char *str, int *nrow, int *ncol );

int num_words_in_string(char *str);

char *get_word_from_string(char *str, int which);

char *strip_white_space(char *string);

char *clean_string(char *string, int cp_flag);

char **create_str_arr(int num); 

char **init_str_arr(char *str1,...);

char **add_str_to_arr(char **strarr, char *str);

char **remove_string_from_arr(char **strarr, int pos);

char *get_str_from_arr(char **strarr, int j);

int num_str_in_arr(char **strarr);

char **split_string(char *string, char *delim, int *num, int block_flag); 

int clear_str_arr(char **strarr);

int check_str_arr_for_dups(char **strarr);

int find_string_in_arr(char **strarr, char *str);

int *sort_str_arr(char **strs);

int init_print_buff(int use_std_out);

int buff_print(char *f, ...);

int disp_print_buff(int clear);

char **get_print_buff();

int clear_print_buff();

#endif
/*****************************************************************************/
