#pragma once

#include "definesStandard.hpp"

/** Enable logging of warning statements */
#define LOG_WARN_ENABLED  1
/** Enable logging of info statements */
#define LOG_INFO_ENABLED  1
/** Enable logging of debug statements */
#define LOG_DEBUG_ENABLED 1
/** Enable logging of trace statements */
#define LOG_TRACE_ENABLED 1

#if RELEASE == 1
/** DISABLED --> RELEASE flag has been enabled */ // TODO: Add Doxygen ref?
#define LOG_DEBUG_ENABLED 0
/** DISABLED --> RELEASE flag has been enabled */ // TODO: Add Doxygen ref?
#define LOG_TRACE_ENABLED 0
#endif

/* list of logging flags (with an equivalant numeric value) */
typedef enum logLevel{
    LOG_LEVEL_FATAL = 0, /**< used for events that cause fatal crash */
    LOG_LEVEL_ERROR = 1, /**< used for events that are recoverable */
    LOG_LEVEL_WARN  = 2, /**< used for events that may cause issues */
    LOG_LEVEL_INFO  = 3, /**< used for normal events */
    LOG_LEVEL_DEBUG = 4, /**< used for debugging events */
    LOG_LEVEL_TRACE = 5, /**< used for tracing events */
} logLevel;



/************************************************************************************************************************
*  @brief   Writes to the console a log message of \p level severity.
* 
*  @details A message function purely intended for logging purposes. 
* 
*           Example:
* 
*           @code{.cpp}
*           logOutput(LOG_LEVEL_FATAL, "I took a bullet to the %s and i am %s.", "head", "hurt")
*           @endcode
* 
*           would print out:
* 
*           [FATAL] I took a bullet to the head and i am hurt.
*
*  @param level   a level from the @ref logLevel struct, used purely for visual purposes
*  @param message a message with format specifics included (%s, %f, etc... )
*  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
* 
*  @return None
************************************************************************************************************************/
void logOutput(logLevel level, const char* message, ... );



/************************************************************************************************************************ 
 *  @brief Throws a fatal-error log message to the prompt. Does not crash the code. Preferrably used in a function that
 *  invokes a crash.
 * 
 *  @param message a message to throw before the code crashs. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define FATAL_MSG(message, ...) logOutput(LOG_LEVEL_FATAL, message, ##__VA_ARGS__);



#ifndef ERROR_MSG
/************************************************************************************************************************ 
 *  @brief Throws an error log message to the prompt. 
 * 
 *  @param message a message to throw when the code errors. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define ERROR_MSG(message, ...) logOutput(LOG_LEVEL_ERROR, message, ##__VA_ARGS__);
#endif



#if LOG_WARN_ENABLED == 1
/************************************************************************************************************************ 
 *  @brief Throws a warning log message to the prompt. 
 * 
 *  @param message a message to throw when the code issues a warning. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define WARN_MSG(message, ...)  logOutput(LOG_LEVEL_WARN,  message, ##__VA_ARGS__);
#else
/** DISABLED --> @ref LOG_WARN_ENABLED is set to 0 */
#define WARN_MSG(message, ...)
#endif



#if LOG_INFO_ENABLED == 1
/************************************************************************************************************************ 
 *  @brief Throws an info log message to the prompt. 
 * 
 *  @param message a message to throw when there is information to be displayed. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define INFO_MSG(message, ...)  logOutput(LOG_LEVEL_INFO,  message, ##__VA_ARGS__);
#else
/** DISABLED --> @ref LOG_INFO_ENABLED is set to 0 */
#define INFO_MSG(message, ...)
#endif



#if LOG_DEBUG_ENABLED == 1
/************************************************************************************************************************ 
 *  @brief Throws a debug log message to the prompt. 
 * 
 *  @param message a message to throw when there is debug information to be displayed. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define DEBUG_MSG(message, ...) logOutput(LOG_LEVEL_DEBUG, message, ##__VA_ARGS__);
#else
/** DISABLED --> @ref LOG_DEBUG_ENABLED is set to 0 */
#define DEBUG_MSG(message, ...)
#endif



#if LOG_TRACE_ENABLED == 1
/************************************************************************************************************************ 
 *  @brief Throws a trace log message to the prompt. 
 * 
 *  @param message a message to throw when there is trace information to be displayed. Message uses format specifiers.
 *  @param ...     a variadic list that gets appended to the message in place of the format specifiers.
 * 
 *  @return None
 ************************************************************************************************************************/ 
#define TRACE_MSG(message, ...) logOutput(LOG_LEVEL_TRACE, message, ##__VA_ARGS__);
#else
/** DISABLED --> @ref LOG_TRACE_ENABLED is set to 0 */
#define TRACE_MSG(message, ...)
#endif
