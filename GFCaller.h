#pragma once
#ifndef GFLIB_GFCALLER_H

#define GFLIB_GFCALLER_H

#include "GF.h"
#include <comdef.h>
#include<comutil.h>

extern "C" __declspec(dllexport) Element * __stdcall  CreateElement(int const* v, size_t v_size, size_t degree) {
	std::vector<int> vector = std::vector<int>(v_size);
	for (int i = 0; i < v_size; ++i) {
		vector[i] = v[i];
	}
	return new Element(vector, degree);
}

extern "C" __declspec(dllexport) Element * __stdcall  CreateElementCopy(Element const* e) {
	return new Element(*e);
}

extern "C" __declspec(dllexport) size_t __stdcall ElementGetPrimitiveDegree(Element const* e) {
	return e->get_primitive_degree();
}

extern "C" __declspec(dllexport) BSTR __stdcall ElementPrint(Element * e) {
	return _com_util::ConvertStringToBSTR(e->print().c_str());
}

extern "C" __declspec(dllexport) void __stdcall DeleteElement(Element * e) {
	delete e;
}



extern "C" __declspec(dllexport) GF * __stdcall
CreateGF(unsigned int characteristic, unsigned int degree, int const* polynom, size_t polynom_size) {
	std::vector<int> vector = std::vector<int>(polynom_size);
	for (int i = 0; i < polynom_size; ++i) {
		vector[i] = polynom[i];
	}
	return new GF(characteristic, degree, vector);
}
extern "C" __declspec(dllexport) Element * __stdcall GFDivision(GF * field, Element * e1, Element * e2) {
	return new Element(field->division(*e1, *e2));
}
extern "C" __declspec(dllexport) Element * __stdcall GFAdd(GF * field, Element * e1, Element * e2) {
	return new Element(field->add(*e1, *e2));
}
extern "C" __declspec(dllexport) Element * __stdcall GFSub(GF * field, Element * e1, Element * e2) {
	return new Element(field->sub(*e1, *e2));
}
extern "C" __declspec(dllexport) Element * __stdcall GFMultiply(GF * field, Element * e1, Element * e2) {
	return new Element(field->multiply(*e1, *e2));
}

extern "C" __declspec(dllexport) Element * __stdcall GFGetElement(GF * field, size_t i) {
	return new Element(field->get_element(i));
}
extern "C" __declspec(dllexport) BSTR __stdcall GFPrint(GF * field) {
	return _com_util::ConvertStringToBSTR(field->print().c_str());
}
extern "C" __declspec(dllexport) size_t __stdcall GFGetNumber(GF * field) {
	return field->get_number_of_elements();
}
extern "C" __declspec(dllexport) bool __stdcall GFIsElementZero(GF * field, Element * e) {
	return field->is_zero(*e);
}
extern "C" __declspec(dllexport) void __stdcall DeleteGF(GF * field) {
	delete field;
}

extern "C" __declspec(dllexport) Polynom * __stdcall CreatePolynom(GF * field, int const* primitive_degrees, int primitive_degrees_size) {
	std::vector<int> vector = std::vector<int>(primitive_degrees_size);
	for (size_t i = 0; i < primitive_degrees_size; ++i) {
		vector[i] = primitive_degrees[i];
	}
	return new Polynom(*field, vector);
}
extern "C" __declspec(dllexport) Polynom * __stdcall PolynomOperatorPlus(Polynom * p1, Polynom * p2) {
	return *p1 + *p2;
}
extern "C" __declspec(dllexport) Polynom * __stdcall PolynomOperatorSub(Polynom * p1, Polynom * p2) {
	return *p1 - *p2;
}
extern "C" __declspec(dllexport) Polynom * __stdcall PolynomOperatorDiv(Polynom * p1, Polynom * p2) {
	return *p1 / *p2;
}
extern "C" __declspec(dllexport) Polynom * __stdcall PolynomOperatorMul(Polynom * p1, Polynom * p2) {
	return *p1 * *p2;
}
extern "C" __declspec(dllexport) Element * __stdcall PolynomCalculate(Polynom * p, Element * e) {
	return new Element(p->calculate(*e));
}

extern "C" __declspec(dllexport) BSTR __stdcall PolynomPrint(Polynom * p) {
	return _com_util::ConvertStringToBSTR(p->print().c_str());
}

extern "C" __declspec(dllexport) BSTR __stdcall PolynomPrintAsDegree(Polynom * p) {
	return _com_util::ConvertStringToBSTR(p->print_as_degree().c_str());
}
extern "C" __declspec(dllexport) void __stdcall DeletePolynom(Polynom * p) {
	delete p;
}
extern "C" __declspec(dllexport) bool __stdcall IsPolynomPrimitive(unsigned int characteristic, unsigned int degree, int const* polynom, size_t polynom_size) {
	std::vector<int> vector = std::vector<int>(polynom_size);
	for (int i = 0; i < polynom_size; ++i) {
		vector[i] = polynom[i];
	}
	return is_polynom_primitive(characteristic, degree, vector);
}

#endif //GFLIB_GFCALLER_H