//---------------------------------------------------------------------------
#ifndef NexusExceptionClassH
#define NexusExceptionClassH
//---------------------------------------------------------------------------

#include "string.h"

class NexusExceptionClass
{
public:
    NexusExceptionClass (NexusExceptionClass &old_class){
        strcpy (error_string, old_class.getErrorString());
	    strcpy (error_type, old_class.getErrorType());
    };

	NexusExceptionClass (char *error){
        strcpy (error_string, error);
    };

	NexusExceptionClass (char *error, char *type){
        strcpy (error_string, error);
	    strcpy (error_type, type);
    };

    char *getErrorString (void){
        return (error_string);
    };

    char *getErrorType (void){
        return (error_type);
    };

private:
    char	error_string[256],
    		error_type[256];
};

#endif

