#include <cmath>
#include <string>
#include "GFCaller.h"
#include  <set>


[[nodiscard]] int Element::get_primitive_degree() const {
	return primitive_degree;
}

Element::Element(std::vector<int> const& v, int degree) {
	value = v;
	primitive_degree = degree;
	zetch = 0;
}

Element::Element(Element const& e) {
	value = e.value;
	primitive_degree = e.primitive_degree;
	zetch = e.zetch;
}

std::string Element::print() {
	std::string s = std::string();
	for (int i = value.size() - 1; i > 1; --i) {
		if (value[i] != 0) {
			s += std::to_string(value[i]) + "a^" + std::to_string(i) + '+';
		}
	}
	if (value.size() > 1 && value[1] != 0) {
		s += std::to_string(value[1]) + "a+";
	}
	if (value[0] != 0 || s.size() == 0) {
		s += std::to_string(value[0]);
	}
	if (s[s.size() - 1] == '+') {
		s.pop_back();
	}
	return s;
}

int Element::into_number(int characteristic) {
	size_t res = 0;
	res += value[0];
	for (int i = 1; i < value.size(); ++i) {
		res += value[i] * characteristic;
		characteristic *= characteristic;
	}
	return res;
}

Element& Element::operator=(Element const& e) = default;


GF::GF(unsigned int characteristic, unsigned int degree, std::vector<int> const& polynom) {//highest degree right
	this->characteristic = characteristic;
	this->degree = degree;
	std::vector<int> primitive = polynom;
	for (int i = 0; i < degree; i++) {
		primitive[i] = (-primitive[i] + characteristic) % characteristic;
	}
	number_of_elements = pow(characteristic, degree);
	field = std::vector<Element>();
	std::vector<int> v = std::vector<int>(degree, 0);
	Element e = Element(v, 0);
	field.push_back(e);
	v.clear();
	v = std::vector<int>(degree, 0);
	v[0] = 1;
	e = Element(v, 0);
	field.push_back(e);
	std::vector<int> prev = std::move(v);
	for (size_t i = 0; i < number_of_elements - 2; ++i) {
		v = std::vector<int>(degree, 0);
		for (unsigned int k = degree - 1; k > 0; --k) {
			v[k] = (prev[k - 1] + primitive[k] * prev[degree - 1]) % characteristic;
		}
		v[0] = (primitive[0] * prev[degree - 1]) % characteristic;
		prev.clear();
		prev = v;
		e = Element(v, i + 1);
		field.push_back(e);
		v.clear();
	}
	std::vector<size_t> values = std::vector<size_t>(number_of_elements);
	std::vector<size_t> values_z = std::vector<size_t>(number_of_elements);
	values[0] = 0;
	values[1] = 1;
	for (size_t i = 1; i < number_of_elements; ++i) {
		values[i] = field[i].into_number(characteristic);
		if (values[i] % characteristic == characteristic - 1) {
			values_z[i] = values[i] - characteristic + 1;
		}
		else {
			values_z[i] = values[i] + 1;
		}
	}
	field[0].zetch = 1;
	size_t calc_size = (number_of_elements - 1) / 2 + 1;
	if (characteristic != 2) {
		field[calc_size].zetch = 0;
	}
	else {
		field[1].zetch = 0;
		++calc_size;
	}
	for (size_t i = 1; i < calc_size; ++i) {
		for (size_t j = 1; j < number_of_elements; ++j) {
			if (values_z[i] == values[j]) {
				field[i].zetch = j;
				break;
			}
		}
	}
	if (characteristic != 2) {
		++calc_size;
	}
	for (size_t i = calc_size; i < number_of_elements; ++i) {
		field[i].zetch = (field[number_of_elements - i + 1].zetch + i - 1) % (number_of_elements - 1);
	}
}

Element GF::division(Element& e1, Element& e2) {
	if (is_zero(e2) || is_zero(e1)) {
		return field[0];
	}
	if (e1.primitive_degree - e2.primitive_degree >= 0) {
		return field[(e1.primitive_degree - e2.primitive_degree + 1)];
	}
	else {
		return field[(e1.primitive_degree - e2.primitive_degree + number_of_elements)];
	}
}

Element GF::multiply(Element& e1, Element e2) {
	if (is_zero(e1) || is_zero(e2)) {
		return field[0];
	}
	if ((e1.primitive_degree + e2.primitive_degree + 1) < number_of_elements) {
		return field[(e1.primitive_degree + e2.primitive_degree + 1)];
	}
	else {
		return field[(e1.primitive_degree + e2.primitive_degree - number_of_elements + 2)];
	}
}

Element GF::add(Element& e1, Element& e2) {
	if (is_zero(e1)) {
		return e2;
	}
	if (is_zero(e2)) {
		return e1;
	}
	if (e1.primitive_degree > e2.primitive_degree) {
		int dif = e1.primitive_degree - e2.primitive_degree;
		return field[(e2.primitive_degree + field[dif + 1].zetch) % number_of_elements];
	}
	else {
		int dif = e2.primitive_degree - e1.primitive_degree;
		return field[(e1.primitive_degree + field[dif + 1].zetch) % number_of_elements];
	}
}

Element GF::sub(Element& e1, Element& e2) {
	if (this->characteristic != 2) {
		if (is_zero(e1)) {
			if (is_zero(e2)) {
				return field[0];
			}
			return field[(e2.primitive_degree + 1 + (number_of_elements / 2) + 1) % number_of_elements];
		}
		else {
			if (e1.zetch == e2.zetch) {
				return field[0];
			}
			return add(field[(e1.primitive_degree + (number_of_elements / 2) + 1) % number_of_elements], e2);
		}
	}
	else {
		return add(e1, e2);
	}
}

std::string GF::print() {
	std::string s = std::string();
	int j = 0;
	for (Element e : field) {
		s += e.print() + '\n';
		j++;
	}
	return s;
}
int GF::get_number_of_elements()
{
	return number_of_elements;
}
bool GF::is_zero(Element& e) {
	return e.zetch == field[0].zetch;
}

Element GF::get_element(int i) {
	return field[i];
}

Polynom::Polynom(GF& field, std::vector<int>& primitive_degrees) {
	value = std::vector<Element>();
	for (size_t primitive_degree : primitive_degrees) {
		value.push_back(field.field[primitive_degree]);
	}
	this->field = &field;
}

Polynom* Polynom::operator+(Polynom& p) {
	Polynom* result = nullptr;
	if (value.size() > p.value.size()) {
		result = new Polynom(*this);
		for (int i = 0; i < p.value.size(); ++i) {
			result->value[i] = field->add(result->value[i], p.value[i]);
		}
	}
	else {
		result = new Polynom(p);
		for (int i = 0; i < value.size(); ++i) {
			result->value[i] = field->add(result->value[i], value[i]);
		}
	}
	result->recalc_degree();
	return result;
}

Polynom* Polynom::operator-(Polynom& p) {
	Polynom* result = nullptr;
	if (value.size() > p.value.size()) {
		result = new Polynom(*this);
		for (int i = 0; i < p.value.size(); ++i) {
			result->value[i] = field->sub(result->value[i], p.value[i]);
		}
	}
	else {
		result = new Polynom(p);
		for (int i = 0; i < value.size(); ++i) {
			result->value[i] = field->sub(result->value[i], value[i]);
		}
	}
	result->recalc_degree();
	return result;
}

Polynom* Polynom::operator/(Polynom& p) {
	std::vector<int> v;
	if (p.value.size() > value.size()) {
		v = { 0 };
		return new Polynom(*field, v);
	}
	v = std::vector<int>(value.size(), 0);
	Polynom* result = new Polynom(*field, v);
	size_t i = p.value.size();
	Polynom tmp = Polynom(*this);
	for (int k = tmp.value.size() - i; k >= 0; --k) {
		result->value[k] = field->division(tmp.value[k + i - 1], p.value[i - 1]);
		for (int j = k + i - 1; j >= k; j--) {
			Element sub = field->multiply(p.value[j - k], result->value[k]);
			tmp.value[j] = field->sub(tmp.value[j], sub);
		}
	}
	return result;
}

Polynom* Polynom::operator*(Polynom& p) {
	std::vector<int> v = std::vector<int>(value.size() + p.value.size() - 1, 0);
	Polynom* result = new Polynom(*field, v);
	v.clear();
	for (size_t i = 0; i < value.size(); ++i) {
		for (size_t j = 0; j < p.value.size(); ++j) {
			Element e = field->multiply(value[i], p.value[j]);
			result->value[i + j] = field->add(result->value[i + j], e);
		}
	}
	return result;
}

std::string Polynom::print() {
	std::string s = std::string();
	for (int i = value.size() - 1; i > 1; --i) {
		if (!field->is_zero(value[i])) {
			s += '(' + value[i].print() + ')' + "x^" + std::to_string(i) + '+';
		}
	}
	if (value.size() > 1 && !field->is_zero(value[1])) {
		s += '(' + value[1].print() + ')' + "x+";
	}
	if (!field->is_zero(value[0]) || s.size() == 0) {
		s += value[0].print();
	}
	if (s[s.size() - 1] == '+') {
		s.pop_back();
	}
	return s;
}

std::string Polynom::print_as_degree() {
	std::string s = std::string();
	for (int i = value.size() - 1; i > 1; --i) {
		if (!field->is_zero(value[i])) {
			s += std::to_string(value[i].get_primitive_degree()) + "x^" + std::to_string(i) + '+';
		}
	}
	if (value.size() > 1 && !field->is_zero(value[1])) {
		s += std::to_string(value[1].get_primitive_degree()) + "x+";
	}
	if (!field->is_zero(value[0]) || s.size() == 0) {
		s += std::to_string(value[0].get_primitive_degree());
	}
	if (s[s.size() - 1] == '+') {
		s.pop_back();
	}
	return s;
}

Element Polynom::calculate(Element& e) {
	Element result = Element(value[value.size() - 1]);
	for (size_t i = value.size() - 2; i > 0; --i) {
		result = field->multiply(result, e);
		result = field->add(result, value[i]);
	}
	return result;
}

Polynom::Polynom(const Polynom& p) {
	value = std::vector<Element>(p.value);
	this->field = p.field;
}

Polynom::~Polynom() {
	value.clear();
	field = nullptr;
}

void Polynom::recalc_degree() {
	for (size_t i = value.size() - 1; i > 0; --i) {
		if (field->is_zero(value[i])) {
			value.pop_back();
		}
		else {
			break;
		}
	}
	if (value.size() == 0) {
		value.push_back(field->field[0]);
	}
}

bool is_polynom_primitive(unsigned int characteristic, unsigned int degree, std::vector<int> const& polynom) {
	std::vector<int> primitive = polynom;
	for (int i = 0; i < degree; i++) {
		primitive[i] = (-primitive[i] + characteristic) % characteristic;
	}
	std::vector<int> prev = std::vector<int>(degree, 0);
	prev[0] = 1;
	int number_of_elements = pow(characteristic, degree);
	for (size_t i = 0; i < number_of_elements - 2; ++i) {
		std::vector<int>  v = std::vector<int>(degree, 0);
		for (unsigned int k = degree - 1; k > 0; --k) {
			v[k] = (prev[k - 1] + primitive[k] * prev[degree - 1]) % characteristic;
		}
		v[0] = (primitive[0] * prev[degree - 1]) % characteristic;
		if (v[0] == 1) {
			bool is_one = true;
			for (int i = 1;i < v.size();++i) {
				if (v[i] != 0) {
					is_one = false;
				}
			}
			if (is_one) {
				return false;
			}
		}
		prev = std::move(v);
		v.clear();
	}
	return true;
}