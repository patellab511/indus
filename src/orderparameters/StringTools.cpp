/* StringTools.cpp
 *
 */

// Standard headers

// Project headers
#include "StringTools.h"

// Constructor
StringTools::StringTools() { }


// Destructor
StringTools::~StringTools() { }


// Converts a string containing "yes/no", "true/false", etc. to the appropriate bool
bool StringTools::stringToBool(const std::string& str) const
{
	std::string lowercase_str = this->toLowercase(str);

	if ( lowercase_str == "yes" || lowercase_str == "true" || 
		   lowercase_str == "1"   || lowercase_str[0] == 'y' ) { 
		return true; 
	}
	else if ( lowercase_str == "no" || lowercase_str == "false" ||
	          lowercase_str == "0"  || lowercase_str[0] == 'n' ) { 
		return false;
	}
	else {
		std::string errorMsg = "Error in StringTools::stringToBool - Could not interpret = \"" 
							             + str + "\" as true/yes/1 or false/no/0. Check your input.\n";
		throw std::runtime_error(errorMsg);
	}
}

// Convert a string to all lowercase
std::string StringTools::toLowercase(const std::string& str) const
{
	std::string lowercase_str = str;
	std::transform(lowercase_str.begin(), lowercase_str.end(), lowercase_str.begin(), ::tolower);
	return lowercase_str;
}


std::string StringTools::getFileExtension(const std::string& file) const
{
	size_t last_period = file.find_last_of(".");
	if ( last_period != std::string::npos ) {
		return file.substr(last_period + 1);
	}
	else {
		return std::string("");
	}
}


void StringTools::split(const std::string& str, std::vector<std::string>& tokens) const
{
	tokens.clear();
	std::stringstream ss(str);
	std::string token;
	while ( ss >> token ) {
		tokens.push_back(token);
	}
}


void StringTools::removeTrailingCommentFromLine(
	const std::string& str, const std::string& comment_chars,
	std::string& trimmed_str, std::string& comment
) const
{
	size_t pos_comment = str.find_first_of(comment_chars);

	if ( pos_comment != std::string::npos ) {
		if ( pos_comment > 0 ) {
			size_t trimmed_str_length = pos_comment - 1;
			trimmed_str = str.substr(0, trimmed_str_length);
		}
		else {
			// Comment is the entire line
			trimmed_str = "";
		}
		comment = str.substr(pos_comment + comment_chars.size());
	}
	else {
		// No comment
		trimmed_str = str;
		comment     = "";
	}
}

std::string StringTools::readAndUnwrapInputBlock(
	std::istream& input_stream, const std::string& comment_chars,
	const std::string& opening_token, const std::string& closing_token) const
{
	std::string unwrapped_block;

	bool found_opening_token = false, found_closing_token = false;

	std::string line, token, trimmed_line, comment;
	
	while ( std::getline(input_stream, line) ) {
		// Clean up the string
		removeTrailingCommentFromLine(line, comment_chars,
		                              trimmed_line, comment);
		trimmed_line = trimWhitespace(trimmed_line);
		if ( trimmed_line.empty() ) { continue; }

		// Read the line, token-by-token
		std::stringstream ss(trimmed_line);
		while ( ss >> token ) {
			if ( token == closing_token ) {
				found_closing_token = true;
				break;
			}
			else if ( token == opening_token ) {
				if ( found_opening_token ) {
					const std::string err_msg = 
						"Found \"" + opening_token + "\" inside a block. "
						+ "Nested blocks are not allowed.\n";
					throw std::runtime_error(err_msg);
				}
				// Discard "spare" tokens from before the block
				found_opening_token = true;
				unwrapped_block = "";
			}
			else {
				unwrapped_block += " " + token;
			}
		}

		if ( found_closing_token ) { break; }
	}

	// Ensure that the closing token was found
	if ( not found_closing_token ) {
		const std::string err_msg = "Expected \"" + closing_token 
		                             + "\" at the end of the block.\n";
		throw std::runtime_error(err_msg);
	}

	return unwrapped_block;
}


// Template specialization for std::string: direct copy is possible
template<>
void StringTools::stringToValue<std::string>(const std::string& str, std::string& value) const {
	value = str;
}


