#ifndef ERRORLOGCLASS_H_
#define ERRORLOGCLASS_H_

#include <string.h>
#include <iostream>
#include <fstream>
using namespace std;

//__________________________________________________________________________

class errorlogclass
{
private:
char        error_file_name[256],
            error_file_path[256],
            error_file_location[512];
ofstream    *error_file;

public:
	errorlogclass(): error_file(NULL){ strcpy (error_file_name, "none");};

    void setErrorFileLocation (const char *file_path, const char *file_name);
    void addError (const char *error_str, const char *error_location);
    void addAutoResolution (const char *resolution_str);
    void close (void);

	virtual ~errorlogclass(){
         if (error_file != NULL){
             error_file->close ();
             delete (error_file);
         }
    };
};

//__________________________________________________________________________

#endif /*ERRORLOGCLASS_H_*/
