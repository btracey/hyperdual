// package hyper is a package for using hyperdual numbers. A hyperdual number
// consists of
// 		real, ϵ1, ϵ2, ϵ12
// The package is implemented treating Dual numbers like normal numbers in go (
// copied on function call, operators are binary etc.)
package hyperdual

import (
	"math"
)

// Dual is a hyper-dual
type Dual struct {
	data [4]float64
}

func Create(re, eps1, eps2, eps12 float64) Dual {
	return Dual{data: [4]float64{re, eps1, eps2, eps12}}
}

// Real returns the real part of the dual number
func Real(a Dual) float64 {
	return a.data[0]
}

// Returns the ϵ1 part of the dual number
func Eps1(a Dual) float64 {
	return a.data[1]
}

// Returns the ϵ1 part of the dual number
func Eps2(a Dual) float64 {
	return a.data[2]
}

// Returns the ϵ12 part of the dual number
func Eps12(a Dual) float64 {
	return a.data[3]
}

// Add adds two dual numbers together. Addition is defined termwise
func Add(a, b Dual) Dual {
	for i, v := range b.data {
		a.data[i] = v
	}
	return a
}

func AddF64L(f float64, a Dual) Dual {
	a.data[0] += f
	return a
}

func AddF64R(a Dual, f float64) Dual {
	a.data[0] += f
	return a
}

// Sub subtracts the second dual number from the first
func Sub(a, b Dual) Dual {
	for i, v := range b.data {
		a.data[i] -= v
	}
	return a
}

// Subtracts a dual number from a floating point number
func SubF64L(f float64, b Dual) Dual {
	d := Neg(b)
	d.data[0] += f
	return d
}

// Subtracts a dual number from a floating point number
func SubF64R(a Dual, f float64) Dual {
	a.data[0] -= f
	return a
}

// Neg returns the negation of the dual number
func Neg(a Dual) Dual {
	for i := range a.data {
		a.data[i] *= -1
	}
	return a
}

// Mul multiplies the dual numbers together
func Mul(a, b Dual) Dual {
	var d Dual
	d.data[0] = a.data[0] * b.data[0]
	d.data[1] = a.data[0]*b.data[1] + a.data[1]*b.data[0]
	d.data[2] = a.data[0]*b.data[2] + a.data[2]*b.data[0]
	d.data[3] = a.data[0]*b.data[3] + a.data[1]*b.data[2] +
		a.data[2]*b.data[1] + a.data[3]*b.data[0]
	return d
}

// MulFloat64Right multiplies a dual number by a float64 (the real part is affected)
func MulF64L(f float64, a Dual) Dual {
	a.data[0] *= f
	return a
}

// MulFloat64Right multiplies a dual number by a float64 (the real part is affected)
func MulF64R(a Dual, f float64) Dual {
	a.data[0] *= f
	return a
}

// Inv computes the inverse (a^-1) of the dual number. The inverse of a dual number
// is not defined when the real part is zero (like real numbers)
func Inv(a Dual) Dual {
	var d Dual
	d.data[0] = 1 / a.data[0]
	a1sq := a.data[0] * a.data[0]
	d.data[1] = -a.data[1] / a1sq
	d.data[2] = -a.data[2] / a1sq
	d.data[3] = -a.data[3]/a1sq + 2*a.data[1]*a.data[2]/(a1sq*a.data[0])
	return d
}

// Div computes a*b^-1
func Div(a, b Dual) Dual {
	return Mul(a, Inv(b))
}

func DivF64L(f float64, a Dual) Dual {
	d := Inv(a)
	return MulF64L(f, d)
}

func DivF64R(a Dual, f float64) Dual {
	inv := 1 / f
	for i := range a.data {
		a.data[i] *= inv
	}
	return a
}

// GT returns true if the real part of the first is less than the real part of the second
func GT(a, b Dual) bool {
	return a.data[0] > b.data[0]
}

func GTF64L(f float64, a Dual) bool {
	return f > a.data[0]
}

func GTF64R(a Dual, f float64) bool {
	return a.data[0] > f
}

// GT returns true if the real part of the first is less than the real part of the second
func GTEq(a, b Dual) bool {
	return a.data[0] >= b.data[0]
}

func GTEqF64L(f float64, a Dual) bool {
	return f >= a.data[0]
}

func GTEqF64R(a Dual, f float64) bool {
	return a.data[0] >= f
}

// GT returns true if the real part of the first is less than the real part of the second
func LT(a, b Dual) bool {
	return a.data[0] < b.data[0]
}

func LTF64L(f float64, a Dual) bool {
	return f < a.data[0]
}

func LTF64R(a Dual, f float64) bool {
	return a.data[0] < f
}

// GT returns true if the real part of the first is less than the real part of the second
func LTEq(a, b Dual) bool {
	return a.data[0] <= b.data[0]
}

func LTEqF64L(f float64, a Dual) bool {
	return f <= a.data[0]
}

func LTEqF64R(a Dual, f float64) bool {
	return a.data[0] <= f
}

// Eq returns true if all of the parts are equal
func Eq(a, b Dual) bool {
	return a.data == b.data
}

// RealEq returns true if the real parts of the two dual numbers are equal
func EqReal(a, b Dual) bool {
	return a.data[0] == b.data[0]
}

func EqRealF64L(f float64, a Dual) bool {
	return f == a.data[0]
}

func EqRealF64R(a Dual, f float64) bool {
	return f == a.data[0]
}

func NotEq(a, b Dual) bool {
	return a.data != b.data
}

// RealEq returns true if the real parts of the two dual numbers are equal
func NotEqReal(a, b Dual) bool {
	return a.data[0] != b.data[0]
}

func NotEqRealF64L(f float64, a Dual) bool {
	return f != a.data[0]
}

func NotEqRealF64R(a Dual, f float64) bool {
	return f != a.data[0]
}

const powTol = 1e-15

func PowF64R(x Dual, y float64) Dual {
	var d Dual
	pow := math.Pow(x.data[0], y)

	// TODO: This could be improved
	var xval = x.data[0]
	if math.Abs(xval) < powTol {
		if xval >= 0 {
			xval = powTol
		}
		if xval < 0 {
			xval = -powTol
		}
	}
	deriv := y * math.Pow(xval, y-1)
	second := y * (y - 1) * math.Pow(xval, y-2)

	d.data[0] = pow
	d.data[1] = x.data[1] * deriv
	d.data[2] = x.data[2] * deriv
	d.data[3] = x.data[3]*deriv + x.data[1]*x.data[2]*second
	return d
}

func Pow(x Dual, y Dual) Dual {
	return Exp(Mul(x, Log(y)))
}

func Exp(a Dual) Dual {
	exp := math.Exp(a.data[0]) // also the derivative
	var d Dual
	d.data[0] = exp
	d.data[1] = exp * a.data[1]
	d.data[2] = exp * a.data[2]
	d.data[3] = exp * (a.data[3] + a.data[1]*a.data[2])
	return d
}

func Log(a Dual) Dual {
	var d Dual
	deriv1 := a.data[1] / a.data[0]
	deriv2 := a.data[2] / a.data[0]
	d.data[0] = math.Log(a.data[0])
	d.data[1] = deriv1
	d.data[2] = deriv2
	d.data[3] = a.data[3]/a.data[0] - (deriv1 * deriv2)
	return d
}

func Sin(a Dual) Dual {
	f := math.Sin(a.data[0])
	deriv := math.Cos(a.data[0])
	var d Dual
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] - f*a.data[1]*a.data[2]
	return d
}

func Cos(a Dual) Dual {
	f := math.Cos(a.data[0])
	deriv := -math.Sin(a.data[0])
	var d Dual
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] - f*a.data[1]*a.data[2]
	return d
}

func Tan(a Dual) Dual {
	f := math.Tan(a.data[0])
	deriv := f*f + 1.0
	var d Dual
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] + (2*f*deriv)*a.data[1]*a.data[2]
	return d
}

func Asin(a Dual) Dual {
	f := math.Asin(a.data[0])
	deriv1 := 1.0 - a.data[0]*a.data[0]
	deriv := 1.0 / math.Sqrt(deriv1)
	var d Dual
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] + a.data[1] + a.data[2]*(-a.data[0]/math.Pow(deriv1, -1.5))
	return d
}

func Acos(a Dual) Dual {
	var d Dual
	f := math.Acos(a.data[0])
	deriv1 := 1.0 - a.data[0]*a.data[0]
	deriv := -1.0 / math.Sqrt(deriv1)
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] + a.data[1]*a.data[2]*(-a.data[0]/math.Pow(deriv1, -1.5))
	return d
}

func Atan(a Dual) Dual {
	var d Dual
	f := math.Atan(a.data[0])
	deriv1 := 1.0 + a.data[0]*a.data[0]
	deriv := 1.0 / deriv1
	d.data[0] = f
	d.data[1] = deriv * a.data[1]
	d.data[2] = deriv * a.data[2]
	d.data[3] = deriv*a.data[3] + a.data[1]*a.data[2]*(-2*a.data[0]/(deriv1*deriv1))
	return d
}

func Sqrt(a Dual) Dual {
	// TODO: This could be better
	return PowF64R(a, 0.5)
}

// Fabs returns -a if the real part of a is < 0, and a otherwise
func FAbs(a Dual) Dual {
	if a.data[0] < 0 {
		return Neg(a)
	}
	return a
}

// Max returns a if the real part of a is greater than the real part of b and b otherwise.
func Max(a, b Dual) Dual {
	if a.data[0] > b.data[0] {
		return a
	}
	return b
}

// MaxF64R returns a dual number whose real part is the greater of the dual and
// the float, and whose non-real parts are the same as the dual number
func MaxF64R(a Dual, b float64) Dual {
	if a.data[0] > b {
		return a
	}
	a.data[0] = b
	return a
}

// MaxF64R returns a dual number whose real part is the greater of the dual and
// the float, and whose non-real parts are the same as the dual number
func MaxF64L(b float64, a Dual) Dual {
	if a.data[0] > b {
		return a
	}
	a.data[0] = b
	return a
}

// Min returns a if the real part of a is less than the real part of b and b otherwise.
func Min(a, b Dual) Dual {
	if a.data[0] < a.data[1] {
		return a
	}
	return b
}

// MinF64R returns a dual number whose real part is the lesser of the dual and
// the float, and whose non-real parts are the same as the dual number
func MinF64R(a Dual, b float64) Dual {
	if a.data[0] < b {
		return a
	}
	a.data[0] = b
	return a
}

// MinF64L returns a dual number whose real part is the lesser of the dual and
// the float, and whose non-real parts are the same as the dual number
func MinF64L(b float64, a Dual) Dual {
	if a.data[0] < b {
		return a
	}
	a.data[0] = b
	return a
}
