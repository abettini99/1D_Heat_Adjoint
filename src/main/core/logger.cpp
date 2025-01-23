#include "logger.hpp"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>


void logOutput(logLevel level, const char* message, ...){
    const char* levelStrings[6] = { "[FATAL]: ", 
                                    "[ERROR]: ", 
                                    "[WARN] : ", 
                                    "[INFO] : ", 
                                    "[DEBUG]: ", 
                                    "[TRACE]: "
                                  };
    b8 isError = level < LOG_LEVEL_WARN;

    // avoid dynamic memory allocation
    // 16k char limit, but its OK
    // TODO: clean this up
    const i32 msgLength = 16000;
    char outMessage[msgLength];
    memset(outMessage, 0, sizeof(outMessage));

    // format original message
    // NOTE: MS headers override the GCC/Clang va_list type with "typedef char* va_list" sometimes.
    // Results in a srange error. Workaround is to use __builtin_va_list, tye type GCC/Clang va_start
    // expects.

    // take arg list, create pointer
    __builtin_va_list argPtr; // instead of va_list

    // take all arguments after message -> append all to message
    va_start(argPtr, message);              
    vsnprintf(outMessage, msgLength, message, argPtr);
    va_end(argPtr);

    // prepend some stuff
    char outMessage2[msgLength];
    sprintf(outMessage2, "%s%s\n", levelStrings[level], outMessage);
    
    // TODO: Platform-specific output
    printf("%s", outMessage2);
    
};

