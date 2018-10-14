#include "InputParser.h"

//-------------------------//
//----- ParameterPack -----//
//-------------------------//

ParameterPack::ParameterPack() {}

ParameterPack::ParameterPack(const std::string& pack_name) 
 : name(pack_name) {}

std::ostream& ParameterPack::print(std::ostream& os, const std::string& indent) const
{
	int num_levels = 3;
	std::vector<std::string> indents(num_levels, indent);
	for ( int i=0; i<3; ++i ) {
		for ( int j=0; j<i; ++j ) {
			indents[i] += "  ";
		}		
	}

	os << indents[0] << this->name << " = {\n";

	// Values
	os << indents[1] << "values = {\n";
	for ( auto& entry : this->values ) {
		const std::string& key   = entry.first;
		const std::string& value = entry.second;
		
		// One value per line
		os << indents[2] << key << " = " << value << "\n";
	}
	os << indents[1] << "}\n";

	// Vectors
	os << indents[1] << "vectors = {\n";
	for ( auto& entry : this->vectors ) {
		const std::string& key = entry.first;
		const std::vector<std::string>& this_vector = entry.second;

		// One vector per line
		os << indents[2] << key << " = [";
		int num_values = this_vector.size();
		for ( int i=0; i<num_values; ++i ) {
			os << " " << this_vector[i];
		}
		os << " ]\n";
	}
	os << indents[1] << "}\n";

	// ParameterPacks
	os << indents[1] << "parameter_packs = {\n";
	for ( auto& entry : this->parameter_packs ) {
		const ParameterPack& parameter_pack = entry.second;
		parameter_pack.print(os, indents[2]);
	}
	os << indents[1] << "}\n";

	// End the outer ParameterPack scope
	os << indents[0] << "}\n";

	return os;
}


const std::string* ParameterPack::findValue(
		const std::string& key, const KeyType& key_type
) const
{
	const std::vector<const std::string*> value_ptrs = findValues(key, key_type);
	int num_values = value_ptrs.size();
	if ( num_values == 1 ) {
		return value_ptrs[0];
	}
	else if ( num_values == 0 ) {
		return nullptr;
	}
	else {
		throw DuplicateKeyException(key, num_values, this->name);
	}
}


const std::vector<std::string>* ParameterPack::findVector(
		const std::string& key, const KeyType& key_type
) const
{
	const std::vector<const std::vector<std::string>*> vector_ptrs = findVectors(key, key_type);
	int num_vectors = vector_ptrs.size();
	if ( num_vectors == 1 ) {
		return vector_ptrs[0];
	}
	else if ( num_vectors == 0 ) {
		return nullptr;
	}
	else {
		throw DuplicateKeyException(key, num_vectors, this->name);
	}
}


const ParameterPack* ParameterPack::findParameterPack(
		const std::string& key, const KeyType& key_type
) const
{
	const std::vector<const ParameterPack*> parameter_pack_ptrs = findParameterPacks(key, key_type);
	int num_parameter_packs = parameter_pack_ptrs.size();
	if ( num_parameter_packs == 1 ) {
		return parameter_pack_ptrs[0];
	}
	else if ( num_parameter_packs == 0 ) {
		return nullptr;
	}
	else {
		throw DuplicateKeyException(key, num_parameter_packs, this->name);
	}
}


std::vector<const std::string*> ParameterPack::findValues( 
		const std::string& key, const KeyType& key_type ) const
{
	// Get iterators bounding the matches
	auto it_pair = this->values.equal_range(key);

	// Check the number of matches
	int num_matches = std::distance(it_pair.first, it_pair.second);
	if ( num_matches == 0 and key_type == KeyType::Required ) {
		throw KeyNotFoundException(key, this->name);
	}

	// For user convenience, return a vector of ptrs to the underlying objects
	// rather than a pair of map iterators
	std::vector<const std::string*> value_ptrs;
	for ( auto it = it_pair.first; it != it_pair.second; ++it ) {
		value_ptrs.push_back( &(it->second) );
	}

	return value_ptrs;
}


std::vector<const std::vector<std::string>*> ParameterPack::findVectors( 
		const std::string& key, const KeyType& key_type ) const
{
	// Get iterators bounding the matches
	auto it_pair = this->vectors.equal_range(key);

	// Check the number of matches
	int num_matches = std::distance(it_pair.first, it_pair.second);
	if ( num_matches == 0 and key_type == KeyType::Required ) {
		throw KeyNotFoundException(key, this->name);
	}

	// For user convenience, return a vector of ptrs to the underlying objects
	// rather than a pair of map iterators
	std::vector<const std::vector<std::string>*> vector_ptrs;
	for ( auto it = it_pair.first; it != it_pair.second; ++it ) {
		vector_ptrs.push_back( &(it->second) );
	}

	return vector_ptrs;
}


std::vector<const ParameterPack*> ParameterPack::findParameterPacks(
		const std::string& key, const KeyType& key_type ) const
{
	// Get iterators bounding the matches
	auto it_pair = this->parameter_packs.equal_range(key);

	// Check the number of matches
	int num_matches = std::distance(it_pair.first, it_pair.second);
	if ( num_matches == 0 and key_type == KeyType::Required ) {
		throw KeyNotFoundException(key, this->name);
	}

	// For user convenience, return a vector of ptrs to the underlying objects
	// rather than a pair of map iterators
	std::vector<const ParameterPack*> parameter_pack_ptrs;
	for ( auto it = it_pair.first; it != it_pair.second; ++it ) {
		parameter_pack_ptrs.push_back( &(it->second) );
	}

	return parameter_pack_ptrs;
} 


//-----------------------//
//----- InputParser -----//
//-----------------------//


InputParser::InputParser() {}


// Extract the next non-comment token from the TokenStream
// TODO parsing strings which consist of multiple tokens
//  - e.g. "atom_type OW" should be treated as one token, not two
TokenStream::Status TokenStream::getNextToken( std::string& token ) 
{
	token.clear();
	std::string line;

	if ( line_stream_ >> token ) {
		// Successfully extracted a token
		if ( token.find_first_of(comment_chars_) == std::string::npos ) {
			// Valid token
			return TokenStream::Status::Success;
		}
		else {
			// Comment encountered: discard this token and the rest of the line
			line_stream_.str("");
			line_stream_.clear();

			// Try again
			return this->getNextToken(token);
		}
	}
	else if ( getline(io_stream_, line) ) {
		// Move on to next line
		line_stream_.clear();   // Reset the stringstream state
		line_stream_.str(line); // Update content

		// Try to get the next token
		return this->getNextToken(token);
	}
	else if ( io_stream_.eof() ) {
		// Done with the stream
		return TokenStream::Status::EndOfStream;
	}
	else {
		throw std::runtime_error("bad state encountered while trying to get next TokenStream token");
	}
}


void InputParser::parseFile( 
	const std::string file, ParameterPack& parameter_pack
) const
{
	// Set up the ParameterPack
	parameter_pack.clear();
	parameter_pack.name = file;

	// Open a TokenStream for the file
	std::fstream ifs(file);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("unable to open input file \"" + file + "\"");
	}
	TokenStream token_stream(ifs);

	using Status = TokenStream::Status;
	Status status;

	// Read the file
	std::string token;
	while (1) {
		status = parseNextEntry(token_stream, parameter_pack);

		if ( status == Status::EndOfStream ) {
			// Done
			break;
		}
	}

	/*
	// DEBUG Print raw tokenization
	while ( token_stream.getNextToken(token) != Status::EndOfStream ) {
		std::cout << "\"" << token << "\"\n";
	}
	*/

	ifs.close();

	return;
}


TokenStream::Status InputParser::parseNextEntry(
	TokenStream& token_stream, ParameterPack& parameter_pack
) const
{
	// Working variables
	std::string key, delimiter, token;

	// Stream parsing status
	using Status = TokenStream::Status;
	Status status;

	// Read the key
	status = token_stream.getNextToken(key);
	if ( status == Status::EndOfStream ) {
		// Done reading the file
		return status;
	}
	else if ( token == "]" ) {
		// End of the current vector
		return Status::ClosingBracket;
	}
	else if ( key == "}" ) {
		// End of the current ParameterPack
		return Status::ClosingBrace;
	}

	// Next token is the "=" delimiter
	status = token_stream.getNextToken(delimiter);
	if ( status != Status::Success ) {
		throw IncompleteEntryException(key, token_stream);
	}
	else if ( delimiter != "=" ) {
		throw MissingDelimiterException(key, token_stream);
	}

	// Next token is the value (or the start of a new sub-object)
	status = token_stream.getNextToken(token);
	if ( status != Status::Success ) {
		throw IncompleteEntryException(key, token_stream);
	}
	else if ( token == "[" ) {
		// Beginning of a new vector
		auto new_vec_pair_it = parameter_pack.vectors.insert( 
				ParameterPack::VectorPair( key, std::vector<std::string>() ) 
		);
		std::vector<std::string>& new_vec = new_vec_pair_it->second;

		// Parse the vector
		status = parseVector(token_stream, new_vec);
		if ( status != Status::ClosingBracket ) {
			throw MissingBracketException(key, token_stream);
		}
	}
	else if ( token == "{" ) {
		// Beginning of a new ParameterPack
		auto new_pack_pair_it = parameter_pack.parameter_packs.insert( 
				ParameterPack::ParameterPackPair( key, ParameterPack() ) 
		);
		ParameterPack& new_pack = new_pack_pair_it->second;
		new_pack.name = key;

		// Parse the ParameterPack
		status = parseParameterPack(token_stream, new_pack);
		if ( status != Status::ClosingBrace ) {
			throw MissingBraceException(key, token_stream);
		}
	}
	else {
		// Simple value
		parameter_pack.values.insert( ParameterPack::ValuePair(key, token) );
	}

	return Status::Success;
}


TokenStream::Status InputParser::parseVector(
	TokenStream& token_stream, std::vector<std::string>& vec
) const
{
	vec.clear();

	using Status = TokenStream::Status;
	Status status;

	// Parse entries in this vector
	std::string token;
	while (1) {
		status = token_stream.getNextToken(token);

		if ( token == "]" ) {
			// Done reading this vector
			return Status::ClosingBracket;
		}
		else if ( status == Status::EndOfStream ) {
			// Unexpected end of stream
			return status;
		}
		else {
			// TODO check for tokens which shouldn't be in a vector
			vec.push_back( token );
		}
	}

	// Program should never get here
	return status;
}


TokenStream::Status InputParser::parseParameterPack(
	TokenStream& token_stream, ParameterPack& parameter_pack
) const
{
	parameter_pack.clear();

	using Status = TokenStream::Status;
	Status status;

	// Parse entries in this ParameterPack
	while (1) {
		status = parseNextEntry(token_stream, parameter_pack);

		if ( status == Status::ClosingBrace or status == Status::EndOfStream ) {
			// Done
			break;
		}
	}

	return status;
}
