#ifndef  RCConfiguration_H
#define  RCConfiguration_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/JetDefinition.hh"

enum ConfigTypes{
	LARGE_R,
	LARGE_R_GROOMED,
	RECLUSTERED_FIXED,
	RECLUSTERED_FLOATING,
	ITER_RECLUSTERED_FIXED,
	ITER_RECLUSTERED_FLOATING,
	LARGE_R_JVF,
	NO_TYPE
};


class RCConfiguration{



private:

	// list of possible configuration types


	// storing what type of configuration we represent
	ConfigTypes 			myType;

	// possible configuration parameters
	float 					pt_frac_float;
	float 					pt_fixed_cut;
	float 					large_r;
	float 					small_r;
	float 					iter_small_r; // This is the *smaller* of the two
	float 					jvf_cut; 
	bool					massless;
	bool					extraSub;
	fastjet::JetAlgorithm	large_algorithm;
	fastjet::JetAlgorithm 	small_algorithm;
	fastjet::JetAlgorithm 	iter_small_algorithm; // This is the *smaller* of the two

	fastjet::Filter*  		filterer;
	fastjet::JetAlgorithm   filtering_algorithm;
	fastjet::Selector		filtering_selector;

public:


	~RCConfiguration(){if(filterer) delete filterer;}

	RCConfiguration(){ filterer = 0; Reset(); }

	// function to clear everything
	void Reset(){
		pt_frac_float 			= -1;
		pt_fixed_cut  			= -1;
		large_r       			= -1;
		small_r       			= -1;
		iter_small_r  			= -1;
		jvf_cut       			= -1;
		massless 				= false;
		large_algorithm 		= fastjet::undefined_jet_algorithm;
		small_algorithm 		= fastjet::undefined_jet_algorithm;
		iter_small_algorithm 	= fastjet::undefined_jet_algorithm;
		filtering_algorithm     = fastjet::undefined_jet_algorithm;
		myType					= NO_TYPE;
		extraSub 				= false;
	}

	// all the possible configuration settings.
	// this seems to be the most expressive/simple way of declaring the configs.
	// note that JVF is an optional argument in all the small_r configs
	// if not provided, will default to ignoring; if provided, will use corr_jvf

	void SetLargeR(const fastjet::JetAlgorithm _large_algorithm, const float _large_r){
		Reset();

		large_algorithm = _large_algorithm;
		large_r 		= _large_r;

		myType 			= LARGE_R;
	}

	//TODO: this needs a better setup for filtering, can't just pass a trimmer because it
	//needs rho, which we really should do internally
	void SetLargeRGroomed(const fastjet::JetAlgorithm _large_algorithm, const float _large_r, const fastjet::JetAlgorithm _filtering_algorithm, fastjet::Selector _filtering_selector, bool doSubtraction = false){
		//fastjet::Filter* _filterer){
		Reset();

		large_algorithm 	= _large_algorithm;
		large_r 			= _large_r;
		//filterer 		= _filterer;

		filtering_algorithm = _filtering_algorithm;
		filtering_selector  = _filtering_selector;
		myType  		= LARGE_R_GROOMED;

		extraSub		= doSubtraction;
	}

	void SetLargeRJVF(const fastjet::JetAlgorithm _large_algorithm, const float _large_r, 
			const fastjet::JetAlgorithm _small_algorithm, const float _small_r, const float _pt_frac_float, const float _jvf_cut){
		Reset();

		large_algorithm 	 = _large_algorithm;
		large_r 	    	 = _large_r;
		small_r 			 = _small_r;
		small_algorithm  	 = _small_algorithm;
		pt_frac_float 		 = _pt_frac_float;
		jvf_cut 			 = _jvf_cut;

		myType 			     = LARGE_R_JVF;

	}

	void SetReclusteredFixed(const fastjet::JetAlgorithm _large_algorithm, const float _large_r,
	  	const fastjet::JetAlgorithm _small_algorithm, const float _small_r, const float _pt_fixed_cut,
	  	bool _massless = false, float _jvf_cut = -999){
		Reset();

		large_algorithm = _large_algorithm;
		large_r 		= _large_r;

		small_algorithm = _small_algorithm;
		small_r 		= _small_r;

		pt_fixed_cut    = _pt_fixed_cut;
		jvf_cut 		= _jvf_cut;

		massless 		= _massless;

		myType 			= RECLUSTERED_FIXED;
	}

	void SetReclusteredFloating(const fastjet::JetAlgorithm _large_algorithm, const float _large_r,
	  	const fastjet::JetAlgorithm _small_algorithm, const float _small_r, const float _pt_frac_float,
	  	bool _massless = false, float _jvf_cut = -999){
		Reset();

		large_algorithm = _large_algorithm;
		large_r 		= _large_r;

		small_algorithm = _small_algorithm;
		small_r 		= _small_r;

		pt_frac_float   = _pt_frac_float;
		jvf_cut 		= _jvf_cut;

		massless 		= _massless;

		myType 			= RECLUSTERED_FLOATING;
	}

	void SetIterativeReclusteredFixed(const fastjet::JetAlgorithm _large_algorithm, const float _large_r,
	  	const fastjet::JetAlgorithm _small_algorithm, const float _small_r, const float _pt_fixed_cut,
	  	const fastjet::JetAlgorithm _iter_small_algorithm, const float _iter_small_r, float _jvf_cut = -999){
		Reset();

		large_algorithm = _large_algorithm;
		large_r 		= _large_r;

		small_algorithm = _small_algorithm;
		small_r 		= _small_r;

		pt_fixed_cut    = _pt_fixed_cut;
		jvf_cut 		= _jvf_cut;

		iter_small_algorithm = _iter_small_algorithm;
		iter_small_r         = _iter_small_r;

		myType 			= ITER_RECLUSTERED_FIXED;
	}


	void SetIterativeReclusteredFloating(const fastjet::JetAlgorithm _large_algorithm, const float _large_r,
	  	const fastjet::JetAlgorithm _small_algorithm, const float _small_r, const float _pt_fixed_cut,
	  	const fastjet::JetAlgorithm _iter_small_algorithm, const float _iter_small_r, float _jvf_cut = -999){
		Reset();

		large_algorithm = _large_algorithm;
		large_r 		= _large_r;

		small_algorithm = _small_algorithm;
		small_r 		= _small_r;

		pt_fixed_cut    = _pt_fixed_cut;
		jvf_cut 		= _jvf_cut;

		iter_small_algorithm = _iter_small_algorithm;
		iter_small_r         = _iter_small_r;

		myType 			= ITER_RECLUSTERED_FIXED;
	}

	// time for getters. boring.

	const ConfigTypes GetConfigType(){
		return myType;
	}

	const float GetLargeR(){
		return large_r;
	}

	const fastjet::JetAlgorithm GetLargeAlgorithm(){
		return large_algorithm;
	}

	const float GetSmallR(){
		if(myType != LARGE_R && myType != LARGE_R_GROOMED && myType != NO_TYPE)
			return small_r;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE SMALL_R FOR TYPE " << myType << std::endl; 

		return -1;
	}

	const fastjet::JetAlgorithm GetSmallAlgorithm(){
		if(myType != LARGE_R && myType != LARGE_R_GROOMED && myType != NO_TYPE)
			return small_algorithm;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE SMALL_ALGORITHM FOR TYPE " << myType << std::endl; 
		return fastjet::undefined_jet_algorithm;
	}

	const fastjet::JetAlgorithm GetIterAlgorithm(){
		if(myType == ITER_RECLUSTERED_FLOATING || myType == ITER_RECLUSTERED_FIXED)
			return iter_small_algorithm;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE ITER_SMALL_ALGORITHM FOR TYPE " << myType << std::endl; 
		return fastjet::undefined_jet_algorithm;
	}

	const float GetFixedCut(){
		if(myType == RECLUSTERED_FIXED || myType == ITER_RECLUSTERED_FIXED)
			return pt_fixed_cut;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE PT_FIXED_CUT FOR TYPE " << myType << std::endl;
		return -1;
	}

	const float GetFloatingCut(){
		if(myType == RECLUSTERED_FLOATING || myType == ITER_RECLUSTERED_FLOATING || myType == LARGE_R_JVF)
			return pt_frac_float;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE PT_FIXED_CUT FOR TYPE " << myType << std::endl;
		return -1;
	}

	const fastjet::JetAlgorithm GetFilteringAlgorithm(){
		if(myType == LARGE_R_GROOMED)
			return filtering_algorithm;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE FILTERING_ALGORITHM FOR TYPE " << myType << std::endl;
		return fastjet::undefined_jet_algorithm;
	}

	fastjet::Selector GetFilteringSelector(){
		if(myType == LARGE_R_GROOMED)
			return filtering_selector;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE filtering_selector FOR TYPE " << myType << std::endl;
		return fastjet::Selector();
	}

	const float GetJVFCut(){
		if(myType == RECLUSTERED_FIXED || myType == RECLUSTERED_FLOATING || myType == LARGE_R_JVF)
			return jvf_cut;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE jvf_cut FOR TYPE " << myType << std::endl;
		return -999;
	}

	const bool GetMassless(){
		if(myType == RECLUSTERED_FIXED || myType == RECLUSTERED_FLOATING)
			return massless;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE massless FOR TYPE " << myType << std::endl;
		return false;
	}

	const bool GetExtraSub(){
		if(myType == LARGE_R_GROOMED)
			return extraSub;
		std::cout << " ERROR: ATTEMPTING TO RETRIEVE extraSub FOR TYPE " << myType << std::endl;
		return false;
	}

	// abandoned efforts to do configurations with overloaded constructors. not very expressive.

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), myType(LARGE_R), filterer(0){}

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r,
	//  const &fastjet::Filter _filterer):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), filterer(_filterer), myType(LARGE_R_TRIMMED){}

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r,
	//  const fastjet::algorithm _small_algorithm, const float _small_r, const float _pt_fixed_cut,
	//  const float _jvf_cut = -10):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), small_algorithm(_small_algorithm), small_r(_small_r),
	// 	pt_fixed_cut(_pt_fixed_cut), jvf_cut(_jvf_cut), myType(RECLUSTERED_FIXED) {}

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r,
	//  const fastjet::algorithm _small_algorithm, const float _small_r, const float _pt_frac_float, bool _float,
	//  const float _jvf_cut = -10):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), small_algorithm(_small_algorithm), small_r(_small_r),
	// 	pt_frac_float(_pt_frac_float), jvf_cut(_jvf_cut), myType(RECLUSTERED_FLOATING) {}

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r,
	//  const fastjet::algorithm _small_algorithm, const float _small_r, const float _pt_fixed_cut,
	//  const fastjet::algorithm _iter_small_algorithm, const float _iter_small_r, const float _jvf_cut = -10):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), small_algorithm(_small_algorithm), small_r(_small_r),
	// 	pt_fixed_cut(_pt_fixed_cut), iter_small_algorithm(_iter_small_algorithm), jvf_cut(_jvf_cut), 
	// 	myType(ITER_RECLUSTERED_FIXED) {}

	// RCConfiguration( const fastjet::algorithm _large_algorithm, const float _large_r,
	//  const fastjet::algorithm _small_algorithm, const float _small_r, const float _pt_frac_float,
	//  const fastjet::algorithm _iter_small_algorithm, const float _iter_small_r, bool _float, const float _jvf_cut = -10):
	// 	large_algorithm(_large_algorithm), large_r(_large_r), small_algorithm(_small_algorithm), small_r(_small_r),
	// 	pt_frac_float(_pt_frac_float), iter_small_algorithm(_iter_small_algorithm), jvf_cut(_jvf_cut), 
	// 	myType(ITER_RECLUSTERED_FLOATING) {}

	// fastjet::algorithm getLargeAlgorithm(){
	// 	return large_algorithm;
	// }




};


#endif
