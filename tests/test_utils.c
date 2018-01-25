#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include "CUnit/Basic.h"
#include "utils.h"

static FILE* tmpFile = NULL;

int init_utils(void) {
	/*if ((tmpFile = mkstemp("test_utils_")) == -1) {
		perror("init_utils");
      return -1;
   }*/
	return 0;
}

int clean_utils(void) {
   /*if (!fclose(tmpFile)) {
		perror("clean_utils");
      return -1;
   }
   tmpFile = NULL;*/
   return 0;
}

void testLOWERCASE(void) {
	const char template[] = "abcd efghi\t 01234\njklmno  !@#$\% pqrst";
	char *buffer;
	int buflen;
	int i, j, pos, result;
	
	buflen = strlen(template);
	buffer = (char*) malloc((buflen+1) * sizeof(char));
	srand((int)getpid());
	for (i = 0; i < 10; i++) {
		strcpy(buffer, template);
		// Change randomly 5 characters
		for (j = 0; j < 5; j++) {
			pos = (float)rand()/RAND_MAX*buflen;
			buffer[pos] = toupper(buffer[pos]);
		}
		result = make_lowercase(buffer);
		CU_ASSERT(result >= 0);
		CU_ASSERT(result <= 5);
		CU_ASSERT(strcmp(template, buffer) == 0);
	}
}

void testSTRIPLINE(void) {
	const char *templates[10] = {
		"\t\tabcd", "\tabcd   ", "\nabcd\t", "\nabcd", "  abcd",
		" abcd\n", "  ab\nd  ", "ab\tc\n", "\tab\td\t", "a  b" };
	char buffer[20];
	int buflen;
	int i, j, pos, result;
	
	for (i = 0; i < 10; i++) {
		strcpy(buffer, templates[i]);
		buflen = stripline(buffer);
		CU_ASSERT(buflen == 4);
	}
}

int main() {
   CU_pSuite pSuite = NULL;

	/* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry()) return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite("Utils", init_utils, clean_utils);
   if (pSuite == NULL) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   if (CU_add_test(pSuite, "test of make_lowercase()", testLOWERCASE) == NULL) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   if (CU_add_test(pSuite, "test of stripline()", testSTRIPLINE) == NULL) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}
