#pragma once

#ifndef UNTITLED_GF_H
#define UNTITLED_GF_H

#include <vector>
#include <string>

//! element of field
class Element {
private:
	int zetch;
	int primitive_degree;
	std::vector<int> value;

	friend class GF;

	int into_number(int characteristic);

public:
	/*! returns element presented in polynomial form
	 * \return string
	*/
	std::string print();

	/*! returns primitive degree of element(if element is zero returns zero)
	 * \return int
	*/
	[[nodiscard]] int get_primitive_degree() const;


	/*! creates element from vector of values, values must be non - negative
	 * \param v is vector of coefficients of element with highter degree right(at he tail of vector), coefficients must be non - negative, and smaller than characteristic of field
	 * \param degree is degree of field
	*/
	Element(std::vector<int> const& v, int degree);

	//! copy constructor
	Element(Element const& e);

	//! sets the element same as rvalue
	Element& operator=(Element const& e);

};

//! Galoi field 
class GF {
private:
	std::vector<unsigned int> primitive_polynom;
	unsigned int characteristic;
	unsigned int degree;
	int number_of_elements;
	std::vector<Element> field;
	friend class Polynom;


public:
	/*! creates galoi field from parameters, building galoi field bigger than 2^20 is unoptimal with this library
	 * \param characteristic is characteristic of field
	 * \param degree is degree of field
	 * \param polynom is vector of coefficients primitive with highter degree right(at the tail of vector),the polynom must be normalized and coefficients must be non - negative, if polynom is nonprimitive field does not work correctly
	*/
	GF(unsigned int characteristic, unsigned int degree, std::vector<int> const& polynom);

	/*! returns element of field that is result of division of first element on second
	 * \param e1 divisible element of field
	 * \param e2 is divider element of field
	 * \return result of division as element of field
	*/
	Element division(Element& e1, Element& e2);


	/*! returns element of field that is result of multiplying of first element on second
	 * \param e1  is first multiplier element of field
	 * \param e2 is second multiplier element of field
	 * \return result of multiplying as element of field
	*/
	Element multiply(Element& e1, Element e2);

	/*! returns element of field that is result of substruction of first element with second
	 * \param e1 is minuend element of field
	 * \param e2 is subtrahend element of field
	 * \return result of substruction as element of field
	*/
	Element sub(Element& e1, Element& e2);

	/*! returns element of field that is result of sum of first element with second
	 * \param e1 is first summand element of field
	 * \param e2 is second summand element of field
	 * \return result of sum as element of field
	*/
	Element add(Element& e1, Element& e2);

	/*! returns element of field that is indexed as first parametr
	 * \param i is index of element in field, index must be not smaller than zero and smaller than characteristic in power of degree without one
	 * \return element of this index
	*/
	Element get_element(int i);

	/*! returns field by elements presented in polynomial form, every element is on new line
	 * \return string
	*/
	std::string print();

	/*! returns number of elements in this field
	 * \return int
	*/
	int get_number_of_elements();

	/*! returns true if element is zero int his field, else returns false
	 * \param e is element of field, element must be from this field
	 * \return int
	*/
	bool is_zero(Element& e);
};

//!Polynom with coefficients that are elements of galoi field
class Polynom {
private:
	std::vector<Element> value;
	GF* field;

	void recalc_degree();

public:
	/*! creates polynom from field and coefficients of the polynom with coefficient of higher degree right(at the tail of vector)
	 * \param field is galoi field
	 * \param primitive_degrees is vector of indexes of coefficients of the polynom with index coefficient of higher degree right(at the tail of vector), coefficients must be non - negative,indexes must be not smaller than zero and smaller than characteristic in power of degree without one, last index must not be zero
	*/
	Polynom(GF& field, std::vector<int>& primitive_degrees);

	//! copy constructor
	Polynom(const Polynom& p);

	/*! returns pointer on polynom that is result of sum of polynoms
	 * \param p summand polynom, polynoms must be from same field
	 * \return polynom pointer
	*/
	Polynom* operator+(Polynom& p);

	/*! returns pointer on polynom that is result of substruction of polynoms
	 * \param p subtrahend polynom, polynoms must be from same field
	 * \return polynom pointer
	*/
	Polynom* operator-(Polynom& p);

	/*! returns pointer on polynom that is result of division of polynoms
	 * \param p divider polynom, polynoms must be from same field
	 * \return polynom pointer
	*/
	Polynom* operator/(Polynom& p);

	/*! returns pointer on polynom that is result of multiplying of polynoms
	 * \param p multiplier polynom, polynoms must be from same field
	 * \return polynom pointer
	*/
	Polynom* operator*(Polynom& p);


	/*! returns polynom with coefficients as elements presented in polynomial form
	 * \return string
	*/
	std::string print();

	/*! returns polynom with coefficients as elements presented as primitive degreece in field
	* \return string
	*/
	std::string print_as_degree();

	/*! returns element of calculating this polynom by this element
	 * \param e element of field, element must be from the same field as polynom
	 * \return element
	*/
	Element calculate(Element& e);

	~Polynom();
};

	/*! return true if polynom is primitive
	 * \param characteristic is characteristic of field
	 * \param degree is degree of field
	 * \param polynom is vector of coefficients primitive with highter degree right(at the tail of vector)
	*/
bool is_polynom_primitive(unsigned int characteristic, unsigned int degree, std::vector<int> const& polynom);
#endif //UNTITLED_GF_H