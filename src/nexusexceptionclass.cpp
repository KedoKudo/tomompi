//---------------------------------------------------------------------------
#pragma hdrstop

#include "string.h"
#include "nexusexceptionclass.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
//---------------------------------------------------------------------------

char *NexusExceptionClass::getErrorString ()
{
	return (error_string);
}

//---------------------------------------------------------------------------

char *NexusExceptionClass::getErrorType ()
{
	return (error_type);
}

//---------------------------------------------------------------------------
