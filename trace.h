/***************************************************************************************
 * Filename: trace.h
 *
 * Description: Header file containing macros for tracing information to the log file
 *              trace.log.
 *
 *              This allows the user to compile in "debug mode".  In order to do this,
 *              use the compile option -DDEBUG_BUILD. If this compile option is not set 
 *              then all macros in this file (TRACE_XXX) are compiled out.
 ***************************************************************************************/

#ifndef __TRACE_H_
#define __TRACE_H_

/***************************************************************************************
 * Debug macros
 ***************************************************************************************/
#ifdef DEBUG_BUILD

/* Trace file */
#define TRACE_FILE "trace.log"

/* Trace at start of main function */
#define TRACE_BEGIN(PROG_NAME) \
	FILE *file_ptr; \
	file_ptr = fopen(TRACE_FILE, "w"); \
	fprintf(file_ptr, "Trace file for %s\n", PROG_NAME); \
	fprintf(file_ptr, "main {\n"); \
	fclose(file_ptr);

/* Trace at start of function */
#define TRACE_START(FN_NAME) \
	FILE *file_ptr; \
	file_ptr = fopen(TRACE_FILE, "a"); \
	fprintf(file_ptr, "%s {\n", FN_NAME); \
	fclose(file_ptr);

/* Trace at end of function */
#define TRACE_END \
	file_ptr = fopen(TRACE_FILE, "a"); \
	fprintf(file_ptr, "}\n"); \
	fclose(file_ptr);

/* Trace text in middle of function */
#define TRACE_TEXT(TEXT) \
	file_ptr = fopen(TRACE_FILE, "a"); \
	fprintf(file_ptr, "    "); \
	fprintf(file_ptr, TEXT); \
	fprintf(file_ptr, "\n"); \
	fclose(file_ptr);

/* Trace text and value in middle of function */
#define TRACE_TEXT_VAL(TEXT, VAL) \
	file_ptr = fopen(TRACE_FILE, "a"); \
	fprintf(file_ptr, "    "); \
	fprintf(file_ptr, TEXT, (VAL)); \
	fprintf(file_ptr, "\n"); \
	fclose(file_ptr);

#else /* DEBUG_BUILD */

/***************************************************************************************
 * Release macros - all set to do nothing
 ***************************************************************************************/
#define TRACE_BEGIN(FN_NAME)
#define TRACE_START(FN_NAME)
#define TRACE_END
#define TRACE_TEXT(TEXT)
#define TRACE_TEXT_VAL(TEXT, VAL)

#endif /* DEBUG_BUILD */

#endif /* __TRACE_H_ */
