// Miscellaneous macros and constexpr functions

#ifndef UTILS_H
#define UTILS_H

#include <string>

// Macro to evaluate member function using a pointer to it
// - One should generally be highly allergic to macros, but isocpp.org states that 
//   this is one of the very few exceptions to the rule
//
#define CALL_MEMBER_FXN(object,pMethod) ((object).*(pMethod))

// Calculate n! recursively
constexpr int factorial(int n) {
	return n > 0 ? n*factorial(n-1) : 1;
}

// Convert 'x' to a string using arcane preprocessor tricks
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

// Prints the file name and line number where it's expanded
#define LOCATION_IN_SOURCE __FILE__ ":" TO_STRING(__LINE__)

// Prints the file name and line number where it's expanded
#define LOCATION_IN_SOURCE_STRING std::string(LOCATION_IN_SOURCE)

// Prints the "pretty" name of the function along in which the macro is expanded,
// along with the associated file name and line number
#define FANCY_FUNCTION "\"" + std::string(__PRETTY_FUNCTION__) + "\" " \
		"(" + LOCATION_IN_SOURCE_STRING + ")"

#endif /* UTILS_H */
