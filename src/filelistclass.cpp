#include "tomompi.h"

//_____________________________________________________________________________________

FileListClass::FileListClass (char *file_name)
{
	strcpy(this->file_name, file_name);

	next_in_list = NULL;
}

//_____________________________________________________________________________________

void FileListClass::Print (void)
{
	printf ("Found file %s\n", file_name);
}

//_____________________________________________________________________________________

FileListClass::~FileListClass ()
{
}

//_____________________________________________________________________________________
