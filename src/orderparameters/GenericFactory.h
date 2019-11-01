// GenericFactory
// - Implements a factory for registering and creating objects with a
//   shared base class based on a key (typically a string)

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
          typename DerivedTypeKey,
          typename... DerivedTypeArgs>
class GenericFactory
{
	// Typedefs
	using DerivedTypeConstructor = std::function< BaseType*(DerivedTypeArgs&&...) >;
	using Registry               = std::map< DerivedTypeKey, DerivedTypeConstructor >;

 public:
	// Returns a handle to the static GenericFactory for this type
	// - Uses the "construct on first use" idiom
	static GenericFactory& factory() { 
		static GenericFactory factory_instance;
		return factory_instance;
	};

	// Constructor and key operators (TODO?)
	//GenericFactory() {};
	//GenericFactory(const GenericFactory&);
	//GenericFactory& operator=(const GenericFactory&);

	// Registers methods to construct derived types
	void registerCreateFunction(const DerivedTypeKey& key, DerivedTypeConstructor constructor_ptr) {
		// See whether the mapping already exists
		typename Registry::const_iterator registry_entry = registry_.find(key);
		if ( registry_entry == registry_.end() ) {
			// New mapping
			registry_.insert( std::make_pair(key, constructor_ptr) );

			//std::cout << "GenericFactory: registered key \"" << key << "\"" << std::endl;  // FIXME DEBUG
		}
		else {
			// Changing an existing mapping is not allowed (likely a programmer mistake)
			// TODO better way to handle when an entry is registered multiple times? Compile-time checks?
			std::stringstream err_msg;
			err_msg << "Error in GenericFactory: key \"" << key << "\" already exists!\n";
			throw std::runtime_error( err_msg.str() ); 
		}
	}

	// Creates an instance of DerivedType by passing "input" to its constructor,
	// then returns a ptr (to BaseType!) to that instance
	BaseType* create(const DerivedTypeKey& key, DerivedTypeArgs&&... input) const {
    typename Registry::const_iterator registry_entry = registry_.find(key);
    if ( registry_entry != registry_.end() ) {
        return registry_entry->second( std::forward<DerivedTypeArgs>(input)... );
    }
		else {
			throw std::runtime_error("ERROR: GenericFactory key " + key + " not found");
		}
	}

	const Registry& get_registry() const {
		return registry_;
	}

 private:
	Registry registry_;
};


// Template-class method for registering <Base,Derived,Key,Input> tuples
template <class    BaseType, 
          class    DerivedType, 
          typename DerivedTypeKey,      // Key which maps to desired DerivedType
          typename... DerivedTypeArgs>  // Types of inputs to DerivedType
class RegisterInFactory
{
 public:
	// Constructor registers the map from "key" to a lambda that produces an 
	// instance of Derived : Base using "const DerivedTypeArgs& input"
	RegisterInFactory(const DerivedTypeKey& key) {
		GenericFactory<BaseType,DerivedTypeKey,DerivedTypeArgs...>::factory().registerCreateFunction(
			// Keyword
			key,
			// Lambda to produce the derived class instance
			[](DerivedTypeArgs&&... input){ return ( new DerivedType( std::forward<DerivedTypeArgs>(input)... ) ); }
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
