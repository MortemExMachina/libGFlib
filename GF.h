#pragma once

#ifndef UNTITLED_GF_H
#define UNTITLED_GF_H

#include <vector>
#include <string>

class Element {
private:
	int zetch;
	int primitive_degree;
	std::vector<int> value;

	friend class GF;

	int into_number(int characteristic);

public:
	std::string print();

	[[nodiscard]] int get_primitive_degree() const;

	Element(std::vector<int> const& v, int degree);

	Element(Element const& e);

	Element& operator=(Element const& e);

};

class GF {
private:
	std::vector<unsigned int> primitive_polynom;
	unsigned int characteristic;
	unsigned int degree;
	int number_of_elements;
	std::vector<Element> field;
	friend class Polynom;


public:
	GF(unsigned int characteristic, unsigned int degree, std::vector<int> const& polynom);

	Element division(Element& e1, Element& e2);

	Element multiply(Element& e1, Element e2);

	Element sub(Element& e1, Element& e2);

	Element add(Element& e1, Element& e2);

	Element get_element(int i);

	std::string print();

	int get_number_of_elements();

	bool is_zero(Element& e);
};

class Polynom {
private:
	std::vector<Element> value;
	GF* field;

	void recalc_degree();

public:
	Polynom(GF& field, std::vector<int>& primitive_degrees);

	Polynom(const Polynom& p);

	Polynom* operator+(Polynom& p);

	Polynom* operator-(Polynom& p);

	Polynom* operator/(Polynom& p);

	Polynom* operator*(Polynom& p);

	std::string print();
	std::string print_as_degree();
	Element calculate(Element& e);

	~Polynom();
};

#endif //UNTITLED_GF_H