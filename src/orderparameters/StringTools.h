/* StringTools.h
 *
 * ABOUT: Object with convenient routines for dealing with strings
 * DEVELOPMENT:
 */

#pragma once

#ifndef STRING_TOOLS_H
#define STRING_TOOLS_H

// Standard headers
#include <algorithm> // for std::transform
#include <array>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class StringTools
{
 public:
	StringTools();
	~StringTools();

	// Converts "yes/no", "true/false", etc. to the appropriate bool
	bool stringToBool(const std::string& str) const;

	// Convert a string to all lowercase
	std::string toLowercase(const std::string& str) const;

	// Returns the file extension of the given file name (if it can be found),
	// else it returns an empty string
	// - File extension is everything from the last period (.) to the end
	// - e.g. if file = "myFile.dat", returns "dat"
	//        if file = "myFile", returns ""
	std::string getFileExtension(const std::string& file) const;

	// Removes leading and trailing whitespace
	// - Solution of StackOverflow user "GManNickG" (Nov 25 '09 at 16:29)
	// - See: https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
	std::string trimWhitespace(
		const std::string& str,
		const std::string& whitespace = " \t" // defines "whitespace"
	) const
	{
		const auto strBegin = str.find_first_not_of(whitespace);
		if (strBegin == std::string::npos) { return ""; }

		const auto strEnd   = str.find_last_not_of(whitespace);
		const auto strRange = strEnd - strBegin + 1;

		return str.substr(strBegin, strRange);
	}

	// Split whitespace-delimited string into tokens
	void split(
		const std::string& str, 
		std::vector<std::string>& tokens
	) const;

	// Trim the comment off the end of a string, where the beginning of the comment is
	// *any* of the characters in "comment_chars" (e.g. the character #)
	void removeTrailingCommentFromLine(
		const std::string& str,
		const std::string& comment_chars,
		std::string& trimmed_str,
		std::string& comment
	) const;

	/* Reads a block of text from an input stream of the form (e.g.):
	 *
	 *     BlockToken some_tokens {
	 *        text and numbers   # I am a comment
   *        text# another comment
	 *        # Yet another comment
   *     }
	 *
	 * Notes on proper usage
	 * - If "opening_token" is found, everything from the input location of
	 *   the stream until "opening_token" is discarded
	 * - Nested blocks are not allowed
	 *   - If "opening_token" is found more than once, an exception is thrown
	 * - Any text after the "closing_token" on the same line is discarded!
	 */
	std::string readAndUnwrapInputBlock(
		std::istream& input_stream,
		const std::string& comment_chars = "#",
		const std::string& opening_token = "{", 
		const std::string& closing_token = "}"
	) const;

	//----- Conversions -----//

	// Convert a std::string to a value of the given type
	template<typename T> 
	void stringToValue(const std::string& str, T& value) const {
		std::stringstream ss( str );
		ss >> value;

		if ( ss.fail() ) {
			std::stringstream err_ss;
			err_ss << "unable to convert " << str << " to a number\n";
			throw std::runtime_error( err_ss.str() );
		}
	}

	template<typename T, std::size_t dim>
	void stringsToArray( const std::vector<std::string>& strings, std::array<T,dim>& arr ) const {
		if ( strings.size() != dim ) {
			throw std::runtime_error("size of vector of strings does not match output array size");
		}

		for ( unsigned i=0; i<dim; ++i ) {
			stringToValue(strings[i], arr[i]);
		}
	}

	template<typename T>
	void stringsToVector( const std::vector<std::string>& strings, std::vector<T>& vec ) const {
		// Ensure the output vector has the correct size
		int len = strings.size();
		vec.resize(len);

		for ( int i=0; i<len; ++i ) {
			stringToValue(strings[i], vec[i]);
		}
	}
};


#endif /* STRING_TOOLS_H */
