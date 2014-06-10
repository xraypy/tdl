/******************************************************************************
* Simple Utility Functions
*
* Authors/modifications
* ---------------------
* Tom Trainor, 7-24-02
* Note several of the functions here are based on material
* from the textbook "The Art and Science of C" by E. Roberts.
*
*
* Notes
* -----
* Note many of the functions below (e.g. string functions and string
* array functions) allocate new memory.  It is the callers
* responsibility to free the allocated memory.
*
* Todo
* ----
* Add cmd processing and data list fncs from data shell (as standalone utils)
******************************************************************************/

//#include <sys/types.h>
//#include <sys/stat.h>
//#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "utils.h"

/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// Utility Functions (Errors, Mem allocation, I/O)
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
* print_error()
*
* Parameters
* ---------
* - error msg (NULL is OK).  you can also pass a formatted string, e.g.
*   print_error("Reading file %s failed\n", fname);
*
* Returns
* -------
* - None
******************************************************************************/
void print_error(char *msg,...)
{
    va_list args;
    int max_err_msg = 500;
    char errbuf[501];
    char *errmsg;
    int errlen;

    va_start(args, msg);          // put va_ptr to the first arg after msg
    vsprintf(errbuf, msg, args);  // put contents of msg and args into errbuf
    va_end(args);                 // set the va_ptr to NULL
    errlen = string_len(errbuf);

    if (errlen > max_err_msg) {
        fprintf(stderr, "Error: Error Message too long\n");
    }
    errmsg = malloc(errlen + 1);
    if (errmsg == NULL) {
        errmsg = "No memory available";
    } else {
        strcpy(errmsg, errbuf);
    }
    fprintf(stderr, "%s\n", errmsg);
    perror("Error: ");
}

/******************************************************************************
* mem_alloc()
*
* Parameters
* ---------
* - number of bytes to allocate
* - ptr to current memory location if expanding
*
* Returns
* -------
* - pointer to new/expanded memory
*
* Notes
* -----
* This function calls malloc or realloc depending on ptr (if
* ptr is NULL then malloc is used)
* Use this function through the new_array or mod_array macros
* defined in the header file
******************************************************************************/
void *mem_alloc(size_t nbytes, void *ptr)
{
    void *result;

    //buff_print("%d \n",nbytes);
    if(nbytes == 0){
        return(NULL);
    }
    if (ptr == NULL){
        result = malloc(nbytes);
        if (result == NULL) print_error("mem_alloc: No memory available malloc");
    } else {
        result = realloc(ptr, nbytes);
        if (result == NULL) print_error("mem_alloc: No memory available relloc");
    }
    //buff_print("mem_alloc allocated %d bytes\n", nbytes);
    return (result);
}

/******************************************************************************
* read_line(), get_line()
*
* Parameters
* ---------
* - pointer to a file
*
* Returns
* -------
* - a string
*
* Notes
* -----
* read_line operates by reading characters from the file
* into a dynamically allocated buffer.  If the buffer becomes
* full before the end of the line is reached, a new buffer
* twice the size of the previous one is allocated.
* This function returns NULL for EOF
*
* get_line is a simple wrapper using stdin for the file;
* all the work is done by read_line.
*
* INITIAL_BUFFER_SIZE -- Initial buffer size for read_line
******************************************************************************/
#define INITIAL_BUFFER_SIZE 100
char *read_line(FILE *infile)
{
    char *line, *nline;
    int n, ch, size;

    n = 0;
    size = INITIAL_BUFFER_SIZE;
    line = create_string(size);
    while ((ch = getc(infile)) != '\n' && ch != EOF ) {
        if (n == size) {
            size *= 2;
            nline = create_string(size);
            strncpy(nline, line, n);
            free(line);
            line = nline;
        }
        /* Add this to read windows and unix files */
        /* ie always ignore '\r' */
        if (ch != '\r') line[n++] = ch;
    }
    if (n == 0 && ch == EOF) {
        free(line);
        return (NULL);
    }
    line[n] = '\0';
    nline = create_string(n);
    strcpy(nline, line);
    free(line);
    return (nline);
}

char *get_line(char *prompt)
{
    buff_print("%s", prompt);
    return (read_line(stdin));
}

char *get_line_default(char *prompt, char *default_val)
{
    char *temp1=NULL, *temp2=NULL, *ret;
    int   len;

    buff_print("%s (%s) >", prompt, default_val);
    temp1 = read_line(stdin);
    temp2 = strip_white_space(temp1);
    free(temp1);

    len = string_len(temp2);
    if (len == 0) {
        ret = copy_string( default_val );
        free(temp2);
    }else{
        ret = temp2 ;
    }
    return(ret);
}

/******************************************************************************
* open_file(), close_file()
* Open and close files in c
*
* Parameters
* ---------
* - file name
*
* Returns
* -------
* - pointer to file or NULL if not found
*
* Notes
* -----
* Simple wrapper around fopen
*
/*****************************************************************************/
FILE  *open_file(char *fname)
{
    FILE  *fid;
    // open the file if the name isnt NULL
    if (fname != NULL){
        if ((fid = fopen(fname,"r")) == NULL){
            printf("File not found %s\n",fname);
        }
        return(fid);
    } else {
        return(NULL);
    }
}

void  close_file(FILE *fid)
{
    if (fid != NULL){
        fclose(fid);
        fid = NULL;
    }
}

/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// String Functions
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
* create_string()
* Create a new string
*
* Parameters
* ---------
* - length of string
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* This function dynamically allocates space for a string of
* len characters, leaving room for the null character at the
* end.
*******************************************************************************/
char *create_string(int len)
{
    //return ( (char *) mem_alloc( (len + 1), NULL ) );
    return ( (char *) mem_alloc( (len + 1)*sizeof(char), NULL ) );
}

/******************************************************************************
* char_to_string()
* Convert a char to a string
*
* Parameters
* ---------
* - single character
*
* Returns
* -------
* - a new string
*
******************************************************************************/
char *char_to_string(char ch)
{
    char *result;
    result = create_string(1);
    result[0] = ch;
    result[1] = '\0';
    return (result);
}

/******************************************************************************
* copy_string()
* Copy a string
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* copy_string copies the string s into dynamically allocated
* storage and returns the new string.
******************************************************************************/
char *copy_string(char *s)
{
    char *newstr;
    //if (s == NULL) print_error("NULL string passed to copy_string");
    if (s == NULL) return (NULL);
    newstr = create_string(string_len(s));
    strcpy(newstr, s);
    return (newstr);
}

/******************************************************************************
* concat_2()
* Concatenate 2 strings
*
* Parameters
* ---------
* - two strings
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* This function concatenates two strings by joining them end
* to end.  For example, concat_2("ABC", "DE") returns the string
* "ABCDE".
******************************************************************************/
char *concat_2(char *s1, char *s2)
{
    char  *s;
    int len1, len2, len;
    int s1_ok = TRUE;
    int s2_ok = TRUE;

    if(s1 == NULL){
        len1 = 0;
        s1_ok = FALSE;
    }else {
        len1 = string_len(s1);
    }
    if(s2 == NULL){
        len2 = 0;
        s1_ok = FALSE;
    }else {
        len2 = string_len(s2);
    }
    len = len1 + len2;
    if (len == 0) return( (char *)NULL );
    s = create_string(len);
    if (s1_ok) strcpy(s, s1);
    if (s2_ok) strcpy(s + len1, s2);
    return (s);
}

/******************************************************************************
* concat_3()
* Concatenate 3 strings
*
* Parameters
* ---------
* - 3 strings
*
* Returns
* -------
* - a new string
******************************************************************************/
char *concat_3(char *s1, char *s2, char *s3)
{
    char *s;
    int len1, len2, len3, len;
    int s1_ok = TRUE;
    int s2_ok = TRUE;
    int s3_ok = TRUE;

    if(s1 == NULL){
        len1 = 0;
        s1_ok = FALSE;
    }else {
        len1 = string_len(s1);
    }
    if(s2 == NULL){
        len2 = 0;
        s2_ok = FALSE;
    }else {
        len2 = string_len(s2);
    }
    if(s3 == NULL){
        len3 = 0;
        s3_ok = FALSE;
    }else {
        len3 = string_len(s3);
    }
    len = len1 + len2 + len3;
    if (len == 0) return( (char *)NULL );
    s = create_string(len);
    if (s1_ok) strcpy(s, s1);
    if (s2_ok) strcpy(s + len1, s2);
    if (s3_ok) strcpy(s + len1 + len2, s3);
    return (s);
}

/******************************************************************************
* concat_strings()
* Concatenate an arbitrary number of strings
*
* Parameters
* ---------
* - variable list of strings
* - Pass NULL as the terminating argument?
*
* Returns
* -------
* - a new string
******************************************************************************/
char *concat_strings(char *str1, ...)
{
    int      count = 0;
    va_list  marker;
    char     *new_str=NULL, *tmp1=NULL, *tmp2=NULL;

    if (str1 == NULL) return (NULL);
    new_str = copy_string(str1);

    va_start( marker, str1 );     /* Initialize variable arguments. */
    while(1){
        tmp1 = va_arg( marker, char* );
        if (tmp1 == NULL) break;
        tmp2  = concat_2(new_str, tmp1);
        free(new_str);
        new_str = tmp2;
        count++;
    }
    va_end( marker );    /* Reset variable arguments. */
    return(new_str);
}

/******************************************************************************
* concat_int_to_string()
* Concatenate an integer and string
*
* Parameters
* ---------
* - string
* - integer
* - precision of the integer
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* Example:  tmp = concat_int_to_string('Filename',9,3);
* tmp => 'Filname009'
******************************************************************************/
char *concat_int_to_string(char *str, int val, int prec)
{
    char buffer[30], *ret_str;
    int  num;
    if (str == NULL)    return(NULL);
    num = sprintf(buffer,"%.*d",prec,val);
    ret_str = concat_2(str,buffer);
    return(ret_str);
}

/******************************************************************************
* dot_delim_string()
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* Takes a whitespace delimited string and make it dot delimited
* there is NO leading or trailing dot
******************************************************************************/
char *dot_delim_string(char *s1)
{

    int num_words, j;
    char *ret_str;

    if (s1 == NULL ) print_error("NULL string passed to dot_delim_string");
    num_words = num_words_in_string(s1);
    ret_str = get_word_from_string(s1,1);
    if (num_words > 1){
        for (j=2;j<=num_words;j++){
            // note memorey leak, should use a temp string
            // so can free old ret_str bf reassigning.....
            ret_str = concat_2(char_to_string('.'),get_word_from_string(s1,j) );
        }
    }
    return(ret_str);
}

/******************************************************************************
* sub_string
* Create a sub string
*
* Parameters
* ---------
* - string
* - integer indicies of start (p1) and end (p2)
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* sub_string returns a copy of the substring of s consisting
* of the characters between index positions p1 and p2,
* inclusive.  The following special cases apply:
*
* 1. If p1 is less than 0, it is assumed to be 0.
* 2. If p2 is greater than the index of the last string
*    position, which is string_len(s) - 1, then p2 is
*    set equal to string_len(s) - 1.
* 3. If p2 < p1, sub_string returns the empty string.
******************************************************************************/
char *sub_string(char *s, int p1, int p2)
{
    int len;
    char *result;

    //if (s == NULL) print_error("NULL string passed to sub_string");
    if (s == NULL) return(NULL);
    len = string_len(s);
    if (p1 < 0) p1 = 0;
    if (p2 >= len) p2 = len - 1;
    len = p2 - p1 + 1;
    //if (len < 0) len = 0;
    if (len < 1) {
        return(NULL);
    }
    result = create_string(len);
    strncpy(result, s + p1, len);
    result[len] = '\0';
    return (result);
}

/******************************************************************************
* string_equal
* Determines if two strings are equal
*
* Parameters
* ---------
* - two strings
*
* Returns
* -------
* This function returns TRUE if the strings s1 and s2 are
* equal.  For the strings to be considered equal, every
* character in one string must precisely match the
* corresponding character in the other.  Uppercase and
* lowercase characters are considered to be different.
******************************************************************************/
int string_equal(char *s1, char *s2)
{
    if (s1 == NULL || s2 == NULL) {
        //print_error("NULL string passed to string_equal");
        return(FALSE);
    }
    if (strcmp(s1, s2) == 0){
        return (TRUE);
    } else {
        return (FALSE);
    }
}

/******************************************************************************
* string_compare()
* Null safe version of strcmp
*
* Parameters
* ---------
* - two strings
*
* Returns
* -------
*   0 if strings are equal
*   > 0 if s1 > s2
*   < 0 of s1 < s2
*
*******************************************************************************/
int string_compare(char *s1, char *s2){
    if (s1 == NULL || s2 == NULL) {
        return(0);
    }
    return(strcmp(s1, s2));
}

/******************************************************************************
* string_len()
* Safe strlen function
*
* Parameters
* ---------
* - string
*
* Returns
* -------
* - length of string
******************************************************************************/
int string_len(char *s)
{
    int len = 0;
    if (s == NULL ) {
        return(len);
    }
    len = (int) strlen(s);
    return (len);
}

/******************************************************************************
* find_char
* Find the index of a character in a string
*
* Parameters
* ---------
* - character
* - string
* - position to start search
*
* Returns
* -------
* Beginning at position start in the string text, this
* function searches for the character ch and returns the
* first index at which it appears or -1 if no match is
* found.
******************************************************************************/
int find_char(char ch, char *text, int start)
{
    char *cptr;

    //if (text == NULL) print_error("NULL string passed to find_char");
    if (text == NULL) return(-1);
    if (start < 0) start = 0;
    if (start > string_len(text)) return (-1);
    cptr = strchr(text + start, ch);
    if (cptr == NULL) return (-1);
    return ((int) (cptr - text));
}

/******************************************************************************
* find_string()
* Find a string in a string
*
* Parameters
* ---------
* - string to search in
* - string to search for
* - start position for search
*
* Returns
* -------
* Beginning at position start in the string text, this
* function searches for the string str and returns the
* first index at which it appears or -1 if no match is
* found.
******************************************************************************/
int find_string(char *text, char *str, int start)
{
    char *cptr;

    if (str == NULL) return(-1); //print_error("NULL   pattern string in find_string");
    if (text == NULL) return(-1); //print_error("NULL text string in   find_string");
    if (start < 0) start = 0;
    if (start > string_len(text)) return (-1);
    cptr = strstr(text + start, str);
    if (cptr == NULL) return (-1);
    return ((int) (cptr - text));
}

/******************************************************************************
* find_string_skip_blocks
*
* Parameters
* ---------
* - string to search in
* - string to search for
* - start position for search
*
* Returns
* -------
* Same as find_string except it ignores occurences of str
* within the blocks (),{},[], ""
*
* This needs some testing, at least with "" blocks
*
******************************************************************************/
int find_string_skip_blocks(char *s, char *str, int p1){

    int j, b, btype;
    int p2, ok, s_len, s2_len;

    if (s == NULL) return(-1); //print_error("NULL string in find_string_skip_blocks");
    if (str == NULL) return(-1); //print_error("NULL string in find_string_skip_blocks");

    s_len = string_len(s);
    s2_len = string_len(str);

    if (p1 < 0) p1 = 0;
    if (p1 > s_len - 1) return(-1);

    b = 0; p2 = -1;
    j = p1;
    ok = 1;

    while ( (ok) && (s[j] !=    '\0') ) {
        if (  strncmp( s+j, str, s2_len) == 0   ) {
            p2 =    j;
            ok =    0;
        } else  if ( btype=is_open_block(s[j]) ) {
            //buff_print("%d\n",btype);
            b = 1; j++;
            while( (b != 0 ) && ( s[j] != '\0') ){
                // look for closed first
                // this should allow for quotes to be used as
                // open and closed blocks ?
                if ( is_close_block(s[j]) == btype ){
                    b--;
                } else if ( is_open_block(s[j]) == btype ){
                    b++;
                }
                j++;
            }
        } else {
            j++;
        }
    }
    return( p2 );
}

int is_open_block( char c){
    if(c == '(' ) return(1);
    if(c == '[' ) return(2);
    if(c == '{' ) return(3);
    if(c == '"' ) return(4);
    return(0);
}

int is_close_block( char c){
    if(c == ')' ) return(1);
    if(c == ']' ) return(2);
    if(c == '}' ) return(3);
    if(c == '"' ) return(4);
    return(0);
}

/******************************************************************************
* string_replace()
* Replace a string
*
* Parameters
* ---------
* - a text string
* - fstr: string to search for
* - rstr: replacement string
* - numrep: count on how many times the string was replaced
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* This replaces all instances of fstr with rstr in str.
* Note this Frees str!!!!
*
******************************************************************************/
char *string_replace(char *str, char *fstr,
                     char *rstr,int *num_rep)
{
    char *tmp1, *tmp2, *tmp3;
    int  p1=0, len, ok = TRUE;
    if ((str == NULL)||(fstr == NULL) ||(rstr == NULL)){
        return(NULL);
    }
    if (strcmp(fstr,rstr)==0) return(str);
    len = string_len(fstr);
    *num_rep = 0;
    while (ok){
        p1 = find_string(str,fstr,0);
        if (p1 < 0) {
            ok = FALSE;
        } else if (p1 > 0){
            tmp1 = sub_string(str,0,p1-1);
            tmp2 = str + p1 + len;
        } else {
            tmp1 = NULL;
            tmp2 = str + len;
        }
        if(ok){
            tmp3 = concat_3(tmp1,rstr,tmp2);
            free(tmp1);
            free(str);
            str = tmp3;
            *num_rep = *num_rep +1;
        }
    }
    return(str);
}

/******************************************************************************
* to_lower()
* Convert a string to lower
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* This function returns a new string with all
* alphabetic characters converted to lower case.
******************************************************************************/
char *to_lower(char *s)
{
    char *result;
    int i;

    if (s == NULL) {
        print_error("NULL string passed to to_lower");
    }
    result = create_string(string_len(s));
    for (i = 0; s[i] != '\0'; i++) result[i] = tolower(s[i]);
    result[i] = '\0';
    return (result);
}

/******************************************************************************
* to_upper
* Convert a string to upper
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* This function returns a new string with all
* alphabetic characters converted to upper case.
******************************************************************************/
char *to_upper(char *s)
{
    char *result;
    int i;

    if (s == NULL) {
        print_error("NULL string passed to to_upper");
    }
    result = create_string(string_len(s));
    for (i = 0; s[i] != '\0'; i++) result[i] = toupper(s[i]);
    result[i] = '\0';
    return (result);
}

/******************************************************************************
* num_to_string(), string_to_num()
* String to numbers a visa versa
*
* Parameters
* ---------
* - string or number to be converted
*
* Returns
* -------
* - string or number
*
* Notes
* -----
* Note: MaxDigits must be larger than the maximum
* number of digits that can appear in a number.
******************************************************************************/
#define MaxDigits 30

char *num_to_string(double d)
{
    char buffer[MaxDigits];
    sprintf(buffer, "%G", d);
    return (copy_string(buffer));
}

double string_to_num(char *s)
{
    double result;
    char dummy;

    if (s == NULL) print_error("NULL string passed to StringToReal");
    if (sscanf(s, " %lg %c", &result, &dummy) != 1) {
        print_error("string_to_num called on illegal number %s", s);
        result = 0.0;
    }
    return (result);
}

int string_to_int(char *s)
{
    int result;
    char dummy;

    if (s == NULL) print_error("NULL string passed to StringToReal");
    if (sscanf(s, " %d %c", &result, &dummy) != 1) {
        print_error("string_to_int called on illegal number %s", s);
        result = 0;
    }
    return (result);
}

/******************************************************************************
* string_to_dbl_array()
* Convert a string to an array of doubles
*
* Parameters
* ---------
* - a string
* - n is the count on the length of the array
*
* Returns
* -------
* - a new array of doubles
*
* Notes
* -----
* Take a string of whitespace delim doubles and make into an array
* Note any words in the string which cant be interpretted as a number
* result in a zero value (this is a property of atof)!
* Could use the string_to_num function to include error checking??
*
******************************************************************************/
double *string_to_dbl_array(char *str, int *n )
{
    int     num_words, j;
    double *arr;
    char   *temp;
    //char *dummy;

    num_words = num_words_in_string(str);
    *n = num_words;
    if (num_words == 0) return( (double *)NULL );

    arr = new_array(num_words, double);
    for(j=0;j<num_words;j++){
        temp = get_word_from_string(str,j+1);
        arr[j] = atof(temp);
        //arr[j] = strtod(temp,&dummy);
        //buff_print("%s\n",temp);
        free(temp);
    }
    return(&arr[0]);
}

/******************************************************************************
* string_to_dbl_array_fmt()
* Convert a formatted string to an array of doubles
*
* Parameters
* ---------
* - a formatted string
* - nrow is number of rows in converted arrays
* - ncols is the number of columns
*
* Returns
* -------
* - new array
*
* Notes
* -----
* This will take a string with formatting and convert it to an arrays
* Semi-colons are used to deliniate rows of an array.
******************************************************************************/
double *string_to_dbl_array_fmt(char *str, int *nrow, int *ncol )
{
    int     n,j,k;
    double *arr;
    char   *temp;
    char  **strarr;
    //char *dummy;

    if (str == NULL){
        *nrow = 0;
        *ncol = 0;
        return (NULL);
    }

    strarr    = split_string(str, ";",   nrow, DEFAULT);
    strarr[0] = clean_string(strarr[0], NOCOPY);
    *ncol     = num_words_in_string(strarr[0]);

    for (j=1;j<*nrow;j++){
        strarr[j] = clean_string(strarr[j], NOCOPY);
        if( *ncol != num_words_in_string(strarr[j])){
            buff_print("Error: rows and cols dont match\n");
            *nrow = 1; *ncol = 1;
            clear_str_arr(strarr);
            return (NULL);
        }
    }

    n = (*nrow) * (*ncol);
    if (n == 0) {
        clear_str_arr(strarr);
        return( (double *)NULL );
    }
    arr = new_array(n, double);

    n = 0;
    for(j=0;j<*nrow;j++){
        for(k=0;k<*ncol;k++){
            temp = get_word_from_string(strarr[j],k+1);
            arr[n] = atof(temp);
            free(temp);
            n++;
        }
    }
    clear_str_arr(strarr);
    return(&arr[0]);
}

/******************************************************************************
* string_to_int_array_fmt()
* Convert a formated string to an array of ints
*
* Parameters
* ---------
* - a formatted string
* - nrow is number of rows in converted arrays
* - ncols is the number of columns
*
* Returns
* -------
* - new array
*
* Notes
* -----
* This will take a string with formatting and convert it to an arrays
* Semi-colons are used to deliniate rows of an array.
******************************************************************************/
int *string_to_int_array_fmt(char *str, int *nrow, int *ncol )
{
    int     n,j,k;
    int    *arr;
    char   *temp;
    char  **strarr;

    if (str == NULL){
        *nrow = 0;
        *ncol = 0;
        return (NULL);
    }

    strarr    = split_string(str, ";",   nrow, DEFAULT);
    strarr[0] = clean_string(strarr[0], NOCOPY);
    *ncol = num_words_in_string(strarr[0]);
    for (j=1;j<*nrow;j++){
        strarr[j] = clean_string(strarr[j], NOCOPY);
        if( *ncol != num_words_in_string(strarr[j])){
            buff_print("Error: rows and cols dont match\n");
            *nrow = 1; *ncol = 1;
            return (NULL);
        }
    }

    n = (*nrow) * (*ncol);
    if (n == 0) {
        clear_str_arr(strarr);
        return( (int *)NULL );
    }
    arr = new_array(n, int);
    n = 0;

    for(j=0;j<*nrow;j++){
        for(k=0;k<*ncol;k++){
            temp = get_word_from_string(strarr[j],k+1);
            arr[n] = atoi(temp);
            free(temp);
            n++;
        }
    }
    clear_str_arr(strarr);
    return(&arr[0]);
}

/*****************************************************************************
* num_words_in_string
* Return the number of words in a string
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* - number of words
*
* Notes
* -----
* Get the number of whitespace delimited words in a string.
* this assumes that all the cntrl character (eg "\n", "\r" etc..)
* have been stripped from str already
******************************************************************************/
int num_words_in_string(char *str)
{
    int num=0;
    int j=0;
    int ok=TRUE;
    int len;

    if (str == NULL) return (0);

    len = string_len(str);

    /* skip through leading whitespace */
    /* note is_white_space('\0') should eval to false */
    while( is_white_space(str[j]) ) j++;
    if (str[j] == '\0') return(0);

    /* apparently we are at a word */
    while (ok) {
        while( !is_white_space(str[j]) && (j<len) ) j++;
        num++;

        while ( is_white_space(str[j]) && (j<len) ) j++;

        if( (str[j] == '\0') || ( j>= len) ) ok=FALSE;
    }
    return num;
}

/******************************************************************************
* get_word_from_string
* Get a specific word from a srting
*
* Parameters
* ---------
* - a string
* - which word
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* Get the "which" word from a string, where the words are space
* delimited.  note this starts with 1, ie which = 1 is the first word
* in the string, the returned string should have NO whitespace
* Note this rets a new string.
******************************************************************************/
char *get_word_from_string(char *str, int which)
{
    int num_words = 0;
    int num=0;
    int j=0;
    int idx_l, idx_h;
    int cont;
    int len;

    if (str == NULL) return(NULL);
    len = string_len(str);
    num_words = num_words_in_string(str);
    /* if the last word is only a single char then it
    misses that in the count ??*/

    if (which > num_words) return NULL;

    /* skip through leading whitespace */
    while( is_white_space(str[j]) ) j++;
    if (str[j] == '\0') return NULL;

    /* apparently we are at a word */
    cont = TRUE;
    while (cont) {
        idx_l = j;
        while( !is_white_space(str[j]) && (j<len) ) j++;
        num++;
        idx_h = j-1;
        if (num == which){
            cont = FALSE;
        } else {
            while (is_white_space(str[j]) && (j < len)) j++;
            if( (str[j] == '\n') || (str[j] == '\0')  || ((j+1)>len) ) cont=FALSE;
            // this last line seems wrong but it works?????
        }
    }
    return (sub_string(str,idx_l,idx_h));
}

/******************************************************************************
* strip_white_space()
* Strip white space from a string
*
* Parameters
* ---------
* - a string
*
* Returns
* -------
* - a new string
*
* Notes
* -----
* This function returns a ptr to the first non-whitespace char of a string
* and insert a '\0' after the last non-whitespace char of string
* This function creates a new string.
*
******************************************************************************/
char *strip_white_space (char *in_string)
{
    char *s, *t, *string, *ret_str;
    if(in_string == NULL) return(NULL);
    string = copy_string(in_string);
    for (s = string; is_white_space(*s); s++);
    if (*s == 0) return (NULL);
    t = s + string_len(s) - 1;
    while (t > s && is_white_space (*t))  t--;
    *++t = '\0';
    ret_str = copy_string(s);
    free(string);
    return (ret_str);
}

/******************************************************************************
* clean_string()
* Clean up a string by removing white space
*
* Parameters
* ---------
* - a string
* - flag to indicate creating a new string
*
* Returns
* -------
* - a string
*
* Notes
* -----
* Replace all ',' '(' ')' '\n', '\r' etc..  with white space
* Note cp_flag determines if a new copy of the string is made
******************************************************************************/
char *clean_string (char *string, int cp_flag)
{
    char *s ;
    int j,len;

    if (cp_flag == MKCOPY){
        s = copy_string(string);
    }else{
        s = string;
    }
    len = string_len(s);

    /* replace all extras with white space */
    for  (j=0;j<len;j++) {
        switch( s[j] ) {
        case '(':
            s[j] = ' ';
            break;
        case ')':
            s[j] = ' ';
            break;
        case '[':
            s[j] = ' ';
            break;
        case ']':
            s[j] = ' ';
            break;
        case '{':
            s[j] = ' ';
            break;
        case '}':
            s[j] = ' ';
            break;
        case ';':
            s[j] = ' ';
            break;
        case ',':
            s[j] = ' ';
            break;
        }
    }
    return s;
}

/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// String Array Functions
// Note: The last element in a string array should be a NULL string
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
* create_str_arr
* Create a string array
*
* Parameters
* ---------
* - number of strings
*
* Returns
* -------
* This returns an array of num+1 elements all initilized to NULL
* Note the smallest possible array has two null entries.
******************************************************************************/
char **create_str_arr(int num)
{
    int j;
    char **strarr = NULL;

    if(num<=0) num = 1;
    //strarr = (char **) malloc((num+1)*sizeof(char *));
    strarr = new_array(num+1, char *);
    for (j=0; j < num+1; j++){
        strarr[j] = (char *)NULL;
    }
    return( strarr );
}

/******************************************************************************
* init_str_arr()
* Initialize a string array
*
* Parameters
* ---------
* - variable list of strings
*
* Returns
* -------
* This return a str array given variable number of str inputs.
* use NULL to terminate the array eg.
*    strarr =   init_str_arr(str1, str2, str3, NULL);
*
* Notes
* -----
* The strings in the string array are new copies of the passed in strings
******************************************************************************/
char **init_str_arr(char *str1,...)
{
    int        count = 0;
    va_list    marker;
    char      *tmp;
    char     **strarr;

    if (str1 == NULL) return (NULL);

    strarr    = create_str_arr(1);
    strarr[0] = copy_string(str1);
    va_start( marker, str1 );    /* Initialize variable arguments. */
    while(1){
        tmp = va_arg( marker, char* );
        if (tmp == NULL) break;
        strarr  = add_str_to_arr(strarr, tmp);
        count++;
    }
    va_end( marker );             /* Reset variable arguments.  */
    return(strarr);
}

/******************************************************************************
* add_str_to_arr()
* Add a string to a string array
*
* Parameters
* ---------
* - string array
* - new string
*
* Returns
* -------
* - string array
*
* Notes
* -----
* This function will add the new string to the end of the strarr
* and returns the pointer to the new strarr.
* Note that the old pointer is invalid since we use mod_array
* Pass a NULL val for strarr to create a new one
******************************************************************************/
char **add_str_to_arr(char **strarr, char *str)
{
    int last,strarr_len;

    //
    if (str == NULL) {
        return (strarr);
    }

    // pass a NULL strarr to create a new one
    if(strarr == NULL){
        strarr = create_str_arr(1);
    }

    // last is the idx of the NULL string in the original strarr
    last = 0;
    while( strarr[last] != NULL )   last++;

    // strarr_len is the num of strings we want in the new strarr
    // note the strarr size will be one larger than this since
    // it always includes a NULL string at the end.
    strarr_len = last+1;
    strarr = mod_array(strarr, strarr_len + 1, char *);
    strarr[last]   = copy_string(str);
    strarr[last+1] = NULL;

    return(strarr);
}

/******************************************************************************
* remove_string_from_arr()
* Remove a string from a sring array
*
* Parameters
* ---------
* - string array
* - which string to remove
*
* Returns
* -------
* - new string array
*
* Notes
* -----
* Given a string array and an index, remove the string at the
* given position and return the smaller array
*
******************************************************************************/
char **remove_string_from_arr(char **strarr, int pos)
{
    int num, j;
    char **new_strarr;

    if (pos < 0) return(strarr);
    num = num_str_in_arr(strarr);
    if (pos >= num) return(strarr);

    new_strarr = create_str_arr(num-1);
    for (j=0;j<num;j++){
        if (j != pos){
            new_strarr = add_str_to_arr(new_strarr, strarr[j]);
        }
    }
    clear_str_arr(strarr);
    return(new_strarr);
}

/******************************************************************************
* get_str_from_arr()
* Get the given str from the strarr.
*
* Parameters
* ---------
* - string array
* - which string
*
* Returns
* -------
* - new string
*
* Notes
* -----
* This returns a copy of the given string
******************************************************************************/
char *get_str_from_arr(char **strarr, int j)
{
    int   num;
    //char *ret_str;
    if (strarr == NULL) return (NULL);
    num = num_str_in_arr(strarr);
    if (num == 0) return(NULL);
    if ( (j < 0) || (j > num - 1) ) return (NULL);
    return( copy_string(strarr[j]) );
}

/******************************************************************************
* num_str_in_arr
* Number of strings in a string array
*
* Parameters
* ---------
* - string array
*
* Returns
* -------
* - number of strings
******************************************************************************/
int num_str_in_arr(char **strarr)
{
    int j;
    if (strarr == NULL) return (0);
    // last is the idx of the NULL string terminating the Arr
    j = 0;
    while( strarr[j] != NULL ) j++;
    return(j);
}

/******************************************************************************
* split_string()
* Split a string and return a string array
*
* Parameters
* ---------
* - string
* - delim
* - num is the number of strings converted
* - block flag indicates if 'blocks' should be ignored
*
* Returns
* -------
* - a new string array
*
* Notes
* -----
* This splits the given string between the delim strings.
* If the first segment of str contains the delim it doesnt get added
* similiar if its at the end.  If there is no delim at the start or end
* these still get added.  eg:
*    str = "this.should.get.delim", delim = "."
* returns a strarr with four elements:
*     this, should , get,  delim
*
* Note the last element of the strarr is a null string
*
* Use block_flag = SKIP_BLOCKS to ignore instances of delim
* inside blocks, use block_flag = DEFAULT to look inside blocks
* (defined in header)
*
* This is pretty good but could use a bit of work...
*
* Note modified this so that every string has whitespace
* stripped from start and end.....
*
* Note improve option: if delim == NULL then split based on whitespace
* ie return blocks as one element of the strarr if block_flag = TRUE
*
******************************************************************************/
char **split_string (char *string, char *delim,
                       int *num, int block_flag)
{
    int    len, delim_len,j;
    int    ok, p1, p2, num_splits;
    char **strarr;
    char  *temp1, *temp2;

    *num = 0;
    if(string == NULL) return (NULL);

    // if delim == NULL split based on whitespace
    // note right now this doesnt use the block thing....
    if(delim == NULL){
        *num = num_words_in_string(string);
        if (*num == 0) return (NULL);
        strarr = create_str_arr(*num);
        for (j=0;j<*num; j++){
            strarr[j] = get_word_from_string(string,j+1);
        }
        return(strarr);
    }

    // otherwise split based on delim
    len = string_len(string);
    delim_len = string_len(delim);

    if( (len == 0) || (delim_len == 0) ){
        *num = 0;
        //return ( &strarr[0]   );
        return(NULL);
    }

    strarr = create_str_arr(1);

    ok = TRUE;
    p1 = p2 = num_splits = 0;
    while(ok){
        if (block_flag == DEFAULT){
            p2 =   find_string(string,delim, p1);
        }else if(block_flag == SKIP_BLOCKS){
            p2 =   find_string_skip_blocks(string,delim, p1);
        }
        // if p2 is greater than zero than it found delim in the str
        // btwn p1 and p2
        if(p2>0){
            temp1 = sub_string(string,p1,p2-1);
            temp2 = strip_white_space(temp1);
            if (temp2!=NULL){
                strarr = add_str_to_arr(strarr, temp2);
                num_splits ++;
            }
            p1 = p2 + delim_len;
            free(temp1);
            free(temp2);

        // if p2 = 0 than the delim is the fist char so skip ahead
        } else if(p2 == 0) {
            p1 = p2 + delim_len;

        // otherwise p2 = -1 which means it cant find it
        } else {
            // add the last segment if there is anything left
            if (p1 < len ){
                temp1 = sub_string(string,p1,len-1);
                temp2 = strip_white_space(temp1);
                if (temp2!=NULL){
                    strarr = add_str_to_arr(strarr,temp2);
                    num_splits ++;
                }
                free(temp1);
                free(temp2);
            }
            ok = FALSE;
        }
    }
    *num = num_splits;
    return( &strarr[0] );
}

/******************************************************************************
* clear_str_arr
* Free the memory allocated for a string array
*
* Parameters
* ---------
* - string array
*
* Returns
* -------
* - number of strings freed
******************************************************************************/
int clear_str_arr(char **strarr)
{
    int j;
    if (strarr == NULL) return(0);
    // last is the idx of the NULL string terminating the Arr
    j = 0;
    while( strarr[j] != NULL ) {
        free(strarr[j]);
        j++;
    }
    free(strarr);
    return(j);
}

/*****************************************************************************
* check_str_arr_for_dups
* See if a strarr has any duplicates
*
* Parameters
* ---------
* - string array
*
* Returns
* -------
* - TRUE/FALSE
******************************************************************************/
int check_str_arr_for_dups(char **strarr)
{
    int j, n, k, ret;
    n = num_str_in_arr(strarr);
    if (n == 0) return(FALSE);

    for (j=0;j<n-1;j++){
        for (k=j+1; k<n; k++){
            ret = string_equal(strarr[j], strarr[k]);
            if (ret == TRUE) return(TRUE);
        }
    }
    return(FALSE);
}

/******************************************************************************
* find_string_in_arr()
* Return idx of a string in a strarr, -1 if not found
*
* Parameters
* ---------
* - string array
* - search string
*
* Returns
* -------
* - index of search string in string array
*   or -1 if not found
******************************************************************************/
int find_string_in_arr(char **strarr, char *str)
{
    int j, n, ret;
    n = num_str_in_arr(strarr);
    if (n == 0) return(-1);
    for (j=0;j<n-1;j++){
        ret = string_equal(strarr[j], str);
        if (ret == TRUE) return(TRUE);
    }
    return(-1);
}

/******************************************************************************
* sort_str_arr()
* Sort a string array
*
* Parameters
* ---------
* - string array
*
* Returns
* -------
* This returns a NEW integer array thats has the indicies of the
* string array listed in alphabetical order
*
******************************************************************************/
int *sort_str_arr(char **strs)
{
    int j,m,cnt,t, *idx;
    char *s1, *s2;

    m = num_str_in_arr(strs);
    if (m < 1) return(NULL);
    idx = new_array(m, int);
    for (j=0;j<m;j++) idx[j] = j;

    cnt = 1;
    while (cnt != 0){
        cnt = 0;
        j = 0;
        while (j < m-1){
            s1 = strs[ idx[j] ];
            s2 = strs[ idx[j+1] ];
            if (string_compare(s1, s2) > 0 ){
                t = idx[j];
                idx[j]=idx[j+1];
                idx[j+1]= t;
                cnt = cnt + 1;
            }
            j = j+1;
        }
    }
    return(idx);
}
/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//Functions  for printing to a global buffer
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
* init_print_buff
* Initialize the print buffer
*
* Parameters
* ---------
* - flag for use of stdout
*
* Returns
* -------
* - SUCCESS/FAILURE
*
* Notes
* -----
* Note PRINTBUFF need to be initialized by the
* user of the library by calling this fcn.
******************************************************************************/
typedef struct {
   int    use_std_out;
   char **list;
} PrintBuffer_t;
PrintBuffer_t PRINTBUFF;
int init_print_buff(int use_std_out)
{
    PRINTBUFF.list = NULL;
    PRINTBUFF.use_std_out = use_std_out;
    return(SUCCESS);
}

/******************************************************************************
* buff_print()
* Print to a global buffer
*
* Parameters
* ---------
* - same as printf()
*
* Notes
* -----
* This will print to the global buffer:
*   char  PRINTBUFF.list;
*
* The flag: PRINT_BUFF.use_std_out
*   TRUE = print to standard out, behaves just like printf
*   FALSE = print to the global buffer
*
* Need to call init_print_buff() function above before calling this....
*
*******************************************************************************/
#define MAX_P_CHAR  8192
int buff_print(char *f, ...)
{
    int ret;
    va_list arg_ptr;
    char temp[MAX_P_CHAR];

    // set arg_ptr to the head of the list
    va_start(arg_ptr,f);

    // if USE_STDOUT == TRUE just print to stdout
    if (PRINTBUFF.use_std_out == TRUE){
        ret = vprintf(f, arg_ptr );
        return(ret);
    }

    // print to the temp variable
    ret = vsprintf(temp, f, arg_ptr );

    // if read more chars than MAX_P_CHAR, got some serious problems!
    if (ret > MAX_P_CHAR){
        printf("Error in buff_print, buffer overflow!\n");
        return(-1);
    } else if (ret < 1) {
        return( ret);
    }
    // add the string in temp to the global PRINTBUFF.list
    PRINTBUFF.list = add_str_to_arr(PRINTBUFF.list, temp);
    va_end( arg_ptr );
    return(ret);
}

/******************************************************************************
* disp_print_buff()
* Print values in global buffer to sdtout
*
* Parameters
* ---------
* - clear: Flag to indicate if buffer should be cleared
*
******************************************************************************/
int disp_print_buff(int clear)
{
    int n,j;
    n = num_str_in_arr(PRINTBUFF.list);
    for (j=0;j<n;j++){
        printf("%s",PRINTBUFF.list[j]);
    }
    if (clear == TRUE) {
        clear_str_arr(PRINTBUFF.list);
        PRINTBUFF.list = NULL;
    }
    return(SUCCESS);
}

/******************************************************************************
* get_print_buff()
* Get the print buffer
*
* Parameters
* ---------
* - none
*
* Returns
* -------
* - string array with current print buffer
*
******************************************************************************/
char **get_print_buff()
{
    //disp_print_buff(FALSE);
    return(PRINTBUFF.list);
}

/******************************************************************************
* clear_print_buff()
* Clear the print buffer
*
* Parameters
* ---------
* - None
*
* Notes
* -----
* - clears the string array associated with the print buffer
******************************************************************************/
int clear_print_buff()
{
    clear_str_arr(PRINTBUFF.list);
    PRINTBUFF.list = NULL;
    return(SUCCESS);
}

/*****************************************************************************/

