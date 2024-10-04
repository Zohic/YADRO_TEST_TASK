#include "matrix_vecor_mapper.hpp"

using namespace YADRO_TEST;

YADRO_TEST::matrix_vector_mapper::getter YADRO_TEST::matrix_vector_mapper::getters[2] = {
	&YADRO_TEST::matrix_vector_mapper::getLine,
	&YADRO_TEST::matrix_vector_mapper::getColumn
};

YADRO_TEST::matrix_vector_mapper::const_getter YADRO_TEST::matrix_vector_mapper::ñ_getters[2] = {
	&YADRO_TEST::matrix_vector_mapper::getLine,
	&YADRO_TEST::matrix_vector_mapper::getColumn
};

matrix_vector_mapper::matrix_vector_mapper(const size_t vecSize, const uint32_t sequence_length, vector_map_direction dir) {
	if (vecSize == 0 || vecSize % sequence_length != 0)
		throw std::exception("wrong matrix mapping to a vector, sizes don't match, index would have been out of range");
	_length = sequence_length;
	_count = vecSize / sequence_length;
	_dir = dir;
}
matrix_vector_mapper::matrix_vector_mapper(const size_t vecSize, const uint32_t sequence_count, const uint32_t sequence_length,  vector_map_direction dir) {
	if (sequence_length * sequence_count != vecSize)
		throw std::exception("wrong matrix mapping to a vector, sizes don't match, index would have been out of range");
	_length = sequence_length;
	_count = sequence_count;
	_dir = dir;
}

complex_t& matrix_vector_mapper::getLine(complex_t* ptr, size_t sequence_index, size_t element_index) {
	if (sequence_index >= _count || element_index >= _length)
		throw std::exception("out of range");
	return ptr[_length * sequence_index + element_index];
}
complex_t& matrix_vector_mapper::getColumn(complex_t* ptr, size_t sequence_index, size_t element_index) {
	if (sequence_index >= _count || element_index >= _length)
		throw std::exception("out of range");
	return ptr[_count * element_index + sequence_index];
}
complex_t& matrix_vector_mapper::get(complex_t* ptr, size_t sequence_index, size_t element_index) {
	return (this->*getters[(size_t)_dir])(ptr, sequence_index, element_index);
}

const complex_t& matrix_vector_mapper::getLine(const complex_t* ptr, size_t sequence_index, size_t element_index) const {
	if (sequence_index >= _count || element_index >= _length)
		throw std::exception("out of range");
	return ptr[_length * sequence_index + element_index];
}
const complex_t& matrix_vector_mapper::getColumn(const complex_t* ptr, size_t sequence_index, size_t element_index) const {
	if (sequence_index >= _count || element_index >= _length)
		throw std::exception("out of range");
	return ptr[_count * element_index + sequence_index];
}
const complex_t& matrix_vector_mapper::get(const complex_t* ptr, size_t sequence_index, size_t element_index) const {
	return (this->*ñ_getters[(size_t)_dir])(ptr, sequence_index, element_index);
}

void matrix_vector_mapper::change_direction() const {
	_dir = (vector_map_direction)(!(uint8_t)_dir);
}
void matrix_vector_mapper::swap_dimensions() const {
	std::swap(_count, _length);
}

size_t matrix_vector_mapper::count() const{
	return _count;
}
size_t matrix_vector_mapper::length() const{
	return _length;
}
