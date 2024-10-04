#pragma once
#include "config.hpp"
#include <vector>

namespace YADRO_TEST {
	enum class vector_map_direction {
		lengthwise,
		crosswise
	};

	class matrix_vector_mapper {
	public:
		typedef complex_t& (matrix_vector_mapper::* getter)(complex_t*, size_t, size_t);
		typedef const complex_t& (matrix_vector_mapper::* const_getter)(const complex_t*, size_t, size_t) const;
	private:
		mutable  uint32_t _count;
		mutable uint32_t _length;
		mutable vector_map_direction _dir;

		static getter getters[2];
		static const_getter ñ_getters[2];

		complex_t& getLine(complex_t* ptr, size_t sequence_index, size_t element_index);
		complex_t& getColumn(complex_t* ptr, size_t sequence_index, size_t element_index);

		const complex_t& getLine(const complex_t* ptr, size_t sequence_index, size_t element_index) const;
		const complex_t& getColumn(const complex_t* ptr, size_t sequence_index, size_t element_index) const;

	public:
		matrix_vector_mapper(const size_t vecSize, const uint32_t sequence_length, vector_map_direction dir);
		
		matrix_vector_mapper(const size_t vecSize, const uint32_t sequence_count, const uint32_t sequence_length,  vector_map_direction dir);
		
		
		complex_t& get(complex_t* ptr, size_t sequence_index, size_t element_index);
		const complex_t& get(const complex_t* ptr, size_t sequence_index, size_t element_index) const;

		void change_direction() const;
		void swap_dimensions() const;

		size_t count() const;
		size_t length() const;
	};
}

