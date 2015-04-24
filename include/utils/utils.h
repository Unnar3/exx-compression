namespace EXX{
	namespace utils{

		template<typename T>
		T l2_norm( T a ){
			T power = 2.0;
			return std::sqrt( std::pow( a, power ) + std::pow( a, power ) );
		}

	}
}