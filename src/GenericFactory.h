

#ifndef GENERIC_FACTORY_H
#define GENERIC_FACTORY_H

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <map>
#include <sstream>
#include <string>

// Sources:
//   Jim Hyslop and Herb Sutter 
//     http://www.drdobbs.com/conversations-abstract-factory-template/184403786
//
//   StackOverflow
//     https://stackoverflow.com/questions/24119808/derived-class-selection-based-on-string-c

template <class    BaseType,
          typename DerivedTypeKey = const std::string,
          typename DerivedTypeInput = const std::string>
class GenericFactory
{
	// Typedefs
	using DerivedTypeConstructor = std::function< BaseType*(DerivedTypeInput&) >;
	using Registry               = std::map< DerivedTypeKey, DerivedTypeConstructor >;

 public:
	// Returns a handle to the static GenericFactory for this type
	// - Uses the "construct on first use" idiom
	static GenericFactory& factory() { 
		static GenericFactory factory_instance;
		return factory_instance;
	};

	// Constructor (and mandatory operators TODO)
	GenericFactory() {};
	//GenericFactory(const GenericFactory&); // Not yet implemented
	//GenericFactory &operator=(const GenericFactory&); // Not yet implemented

	// Registers methods to construct derived types
	void registerCreateFunction(DerivedTypeKey& key, DerivedTypeConstructor constructor_ptr) {
		// See whether the mapping already exists
		typename Registry::const_iterator registry_entry = registry_.find(key);
		if ( registry_entry == registry_.end() ) {
			// New mapping
			registry_[key] = constructor_ptr;
		}
		else {
			// Changing an existing mapping!
			std::stringstream err_msg;
			err_msg << "Error in GenericFactory: key \"" << key << "\" already exists!\n";
			//std::cout << err_msg.str();
			throw std::runtime_error( err_msg.str() ); //TODO better way to handle when an entry is registered multiple times?
		}
	}

	// Creates an instance of DerivedType by passing "input" to its constructor,
	// then returns a ptr (to BaseType!) to that instance
	BaseType* create(DerivedTypeKey& key, DerivedTypeInput& input) const {
    typename Registry::const_iterator registry_entry = registry_.find(key);
    if ( registry_entry != registry_.end() ) {
        return registry_entry->second(input);
    }
		else {
			throw std::runtime_error("ERROR: GenericFactory key not found!");
		}
	}

 private:
	Registry registry_;
};


// Template-class method for registering <Base,Derived,Key,Input> tuples
template <class    BaseType, 
          class    DerivedType, 
          typename DerivedTypeKey = const std::string,   // Key which maps to desired DerivedType
          typename DerivedTypeInput = const std::string> // Type of input to DerivedType
class RegisterInFactory
{
 public:
	// Constructor registers the map from "key" to a lambda that produces an 
	// instance of Derived : Base using "const DerivedTypeInput& input"

	RegisterInFactory(DerivedTypeKey& key) {
		GenericFactory<BaseType,DerivedTypeKey,DerivedTypeInput>::factory().registerCreateFunction(
			// Keyword
			key,
			// Lambda to produce the derived class instance
			[](DerivedTypeInput& input){ return ( new DerivedType(input) ); }
		);
	}
};


/*
// Clunky macro to force static registration
// - Doesn't seem to help!
#define REGISTERME_STATIC_CLASS_NAME(str) str##_##RegisterMe##_##__FILE__##_##__LINE__

// Macro to register a derived Shape with the ShapeFactory
// - Registers the shape in the constructor of a static class, which is auto-initialized
//   by the time main() is called
#define REGISTER_CLASS(BaseClass,DerivedClass,keyword) \
	static class REGISTERME_STATIC_CLASS_NAME(DerivedClass) \
	{ \
	 public: \
		REGISTERME_STATIC_CLASS_NAME(DerivedClass)() { \
			RegisterInFactory<BaseClass,DerivedClass,std::string> registerMe(keyword); \
		} \
		~REGISTERME_STATIC_CLASS_NAME(DerivedClass)() {} \
	} REGISTERME_STATIC_CLASS_NAME(DerivedClass);
*/

#endif /* GENERIC_FACTORY_H */
