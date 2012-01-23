/***************************************************************************

    moltools-python

    Python module for manipulation of atomic coordinates
    Copyright (C) 2012, Borys Szefczyk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************/




#include "moltools.h"


/*****************************************************************************
 *  This is a helper function that reads a line from the file desctiptor fd  *
 *  and returns an allocated string, including the new-line character.       *
 *  If the read operation return no characters an empty string is returned.  *
 *  In any case, the string is terminated with the '\0' character.           *
 *  Make sure, you free the allocated memory after use!                      *
 *****************************************************************************/

#define CHUNK 10
char *readline(FILE *fd) {
	int length, last;
	long int offset;
	char buffer[CHUNK];
	char *line;

	/* Initial offset in the file */
	if ( (offset = ftell(fd)) == -1) return NULL;

	/* Get the first chunk of data; if nothing read, make it empty */
	if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
	last = strlen(buffer) - 1;
	length = last + 1;
	//printf("First chunk:%s:\n", buffer);

	/* Continue reading because this is not the end of the line yet */
	while( buffer[last] != '\n' && last != -1) {

		/* Get the next chunk; if empty - finish here */
		if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
		last = strlen(buffer) - 1;
		length += last + 1;
		//printf("Next  chunk:%s:\n", buffer);
	}
	//printf("Last = %d, length = %d\n", (int)last, length);

	/* Rewind the file */
	if ( fseek(fd, offset, SEEK_SET) != 0 ) return NULL;

	/* Allocate the memory */
	line = (char*) malloc( length + 1 );

	/* Realliy read the data */
	if ( fgets(line, length + 1, fd) == NULL ) buffer[0] = '\0';
	//printf("Line:%s:\n", line);

	return line;
}




int make_lowercase(char *line) {
	int length, i, c, count = 0;

	length = strlen(line);
	for(i = 0; i < length; i++ ) {
		c = tolower(line[i]);
		if ( c != line[i] ) count++;
		line[i] = c;
	}

	return count;
}



int stripline(char *line) {
	int length, i, start, end;
	char c;

	length = strlen(line);
	start = 0;
	while( start < length ) {
		c = line[start++];
		if ( c != ' ' && c != '\t' && c != '\n' && c != '\r' ) break;
	}
	start -= 1;
	end = length-1;
	while( end > 0 ) {
		c = line[end--];
		if ( c != ' ' && c != '\t' && c != '\n' && c != '\r' ) break;
	}
	end += 1;
	//printf("%d %d %d\n", length, start, end);
	length = end-start+1;
	for(i = 0; i < length; i++) {
		line[i] = line[i+start];
	}
	line[i] = '\0';

	return length;
}
