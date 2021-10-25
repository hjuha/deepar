#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

#define EPSILON 0.000000001
#define ZERO_DECIMALS 20
#define TYPE 2 // 0 - real numbers; 1 - complex numbers; 2 - quaternions

using namespace std;

class Matrix;
class Complex;
std::ostream& operator<< (std::ostream&, const Complex&);
std::ostream& operator<< (std::ostream&, const Matrix&);

long double time_since(chrono::steady_clock::time_point start) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(chrono::steady_clock::now() - start).count() / 1000.;
}

inline long double log_sum_exp(const long double x1, const long double x2) {
	if (isnan(x1)) return x2;
	if (isnan(x2)) return x1;
	long double x = max(x1, x2);
	return x + log(exp(x1 - x) + exp(x2 - x));
}

inline long double log_sub_exp(const long double x1, const long double x2) {
	if (isnan(x1)) return x2;
	if (isnan(x2)) return x1;
	long double x = max(x1, x2);
	return x + log(exp(x1 - x) - exp(x2 - x));
}

class AbsLog {
public:
	short sign;
	long double value;

	AbsLog(short sign, long double d) {
		this->sign = sign;
		this->value = d;
		if (d < -ZERO_DECIMALS) {
			this->sign = 0;
			this->value = 0;
		}
	}

	AbsLog(long double d) {
		if (abs(d) < EPSILON) {
			this->sign = 0;
			this->value = 0;
		} else {
			value = log(abs(d));
			if (d < 0) {
				this->sign = -1;
			} else {
				this->sign = 1;
			}
		}
	}

	AbsLog() {
		this->sign = 0;
		this->value = 0;
	}

	AbsLog& operator=(const AbsLog z){
		this->sign = z.sign;
		this->value = z.value;
		return *this;
	}

	void sqrt() {
		this->value /= 2;
	}

	bool operator<(const AbsLog& val) const {
        if (this->sign != val.sign) return this->sign < val.sign;
        if (this->sign == 0) return false;
        if (this->sign == -1) return this->value > val.value;
        return this->value < val.value;
    }
};

inline AbsLog operator+(AbsLog lhs, AbsLog rhs) {
    if (!lhs.sign) return rhs;
    if (!rhs.sign) return lhs;
    if (lhs.sign == rhs.sign) {
    	return AbsLog(lhs.sign, log_sum_exp(lhs.value, rhs.value));
    }
    long double value = log_sub_exp(max(lhs.value, rhs.value), min(lhs.value, rhs.value));
    if (lhs.value > rhs.value) return AbsLog(lhs.sign, value);
    else return AbsLog(rhs.sign, value);
}

inline AbsLog operator-(AbsLog lhs, AbsLog rhs) {
    if (!lhs.sign) return AbsLog(-rhs.sign, rhs.value);
    if (!rhs.sign) return lhs;
    if (lhs.sign != rhs.sign) {
    	return AbsLog(lhs.sign, log_sum_exp(lhs.value, rhs.value));
    }
    long double value = log_sub_exp(max(lhs.value, rhs.value), min(lhs.value, rhs.value));
    if (lhs.value > rhs.value) return AbsLog(lhs.sign, value);
    else return AbsLog(-rhs.sign, value);
}

inline AbsLog operator*(AbsLog lhs, AbsLog rhs) {
    if (!lhs.sign) return lhs;
    if (!rhs.sign) return rhs;
    if (lhs.sign == rhs.sign) return AbsLog(1, lhs.value + rhs.value);
    return AbsLog(-1, lhs.value + rhs.value);
}

inline AbsLog operator/(AbsLog lhs, AbsLog rhs) {
    if (!lhs.sign) return lhs;
    if (lhs.sign == rhs.sign) return AbsLog(1, lhs.value - rhs.value);
    return AbsLog(-1, lhs.value - rhs.value);
}

class Complex {
public:
	AbsLog real;
	AbsLog imag;

	Complex(long double real, long double imag) {
		this->real = AbsLog(real);
		this->imag = AbsLog(imag);
	}

	Complex(AbsLog real, AbsLog imag) {
		this->real = real;
		this->imag = imag;
	}

	Complex(Complex* c) {
		this->real = c->real;
		this->imag = c->imag;
	}

	Complex() {
		this->real = AbsLog(0);
		this->imag = AbsLog(0);
	}

	Complex& operator=(const Complex z){
		this->real = z.real;
		this->imag = z.imag;
		return *this;
	}

	Complex conj() {
		Complex c = Complex(this);
		c.imag = AbsLog(-1) * c.imag;
		return c;
	}

	void sqrt() {
		this->real.sqrt();
	}

	bool zero() {
		return (this->real.value < -ZERO_DECIMALS && this->imag.value < -ZERO_DECIMALS) || (this->real.sign == 0 && this->imag.sign == 0);
	}

	friend ostream& operator<<(ostream&, const Complex&);
};

inline Complex operator+(Complex lhs, Complex rhs) {
	return Complex(lhs.real + rhs.real, lhs.imag + rhs.imag);
}

inline Complex operator-(Complex lhs, Complex rhs) {
	return Complex(lhs.real - rhs.real, lhs.imag - rhs.imag);
}

inline Complex operator*(Complex lhs, Complex rhs) {
	return Complex(lhs.real * rhs.real - lhs.imag * rhs.imag, lhs.real * rhs.imag + lhs.imag * rhs.real);
}

inline Complex operator/(Complex lhs, Complex rhs) {
	return Complex((lhs.real * rhs.real + lhs.imag * rhs.imag) / (rhs.real * rhs.real + rhs.imag * rhs.imag), (lhs.imag * rhs.real - lhs.real * rhs.imag) / (rhs.real * rhs.real + rhs.imag * rhs.imag));
}

class Matrix {
public:
	int n;
	vector<vector<Complex>> v;

	Matrix(vector<vector<Complex>> v) {
		n = v.size();
		this->v = v;
	}

	friend ostream& operator<<(ostream&, const Matrix&);

	Complex determinant() {
		vector<vector<Complex>> w = v;
		int mul = 1;
		for (int i = 0; i < n; i++) {
			if (w[i][i].zero()) {
				for (int j = i + 1; j < n; j++) {
					if (!w[j][i].zero()) {
						mul *= -1;
						for (int k = 0; k < n; k++) {
							swap(w[i][k], w[j][k]);
						}
						break;
					}
				}
				if (w[i][i].zero()) return Complex(0, 0);
			}
			for (int j = i + 1; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (i != k) w[j][k] = w[j][k] - w[j][i] / w[i][i] * w[i][k];
				}
				w[j][i] = Complex(0, 0);
			}
		}
		Complex det(1, 0);
		for (int i = 0; i < n; i++) {
			det = det * w[i][i];
		}
		return Complex(mul, 0) * det;
	}
};

std::ostream& operator<< (std::ostream &out, const Complex& c) {
	out<<"("<<(c.real.sign * exp(c.real.value))<<", "<<(c.imag.sign * exp(c.imag.value))<<")";
	return out;
}

std::ostream& operator<< (std::ostream &out, const Matrix& matrix) {
	for (int i = 0; i < matrix.n; i++) {
		for (int j = 0; j < matrix.n; j++) {
			out<<(matrix.v[i][j])<<' ';
		}
		out<<endl;
	}
	return out;
}

AbsLog godsil(Matrix m) {
	for (int i = 0; i < m.n; i++) {
		for (int j = 0; j < m.n; j++) {
			m.v[i][j].sqrt();
			if (rand() & 1) m.v[i][j] = m.v[i][j] * Complex(-1, 0);
		}
	}
	Complex det = m.determinant();
	return (det * det.conj()).real;
}

AbsLog karmarkar(Matrix m) {
	for (int i = 0; i < m.n; i++) {
		for (int j = 0; j < m.n; j++) {
			m.v[i][j].sqrt();
			if (rand() & 1) {
				if (rand() & 1) m.v[i][j] = m.v[i][j] * Complex(1, 0);
				else m.v[i][j] = m.v[i][j] * Complex(-1, 0);
			} else {
				if (rand() & 1) m.v[i][j] = m.v[i][j] * Complex(0, 1);
				else m.v[i][j] = m.v[i][j] * Complex(0, -1);
			}
		}
	}
	Complex det = m.determinant();
	return (det * det.conj()).real;
}

AbsLog chien(Matrix m) {
	int n = 2 * m.n;
	vector<vector<Complex>> v(n);
	for (int i = 0; i < n; i++) {
		v[i] = vector<Complex>(n);
	}
	for (int i = 0; i < n / 2; i++) {
		for (int j = 0; j < n / 2; j++) {
			m.v[i][j].sqrt();
			Complex c((rand() & 1) ? 1 : -1, 0);
			if (rand() & 1) {
				if (rand() & 1) {
					v[2 * i][2 * j] = m.v[i][j] * c;
					v[2 * i + 1][2 * j + 1] = m.v[i][j] * c;
				} else {
					v[2 * i][2 * j] = m.v[i][j] * c * Complex(0, 1);
					v[2 * i + 1][2 * j + 1] = m.v[i][j] * c * Complex(0, -1);
				}
			} else {
				if (rand() & 1) {
					v[2 * i + 1][2 * j] = Complex(0, 0) - m.v[i][j] * c;
					v[2 * i][2 * j + 1] = m.v[i][j] * c;
				} else {
					v[2 * i + 1][2 * j] = m.v[i][j] * c * Complex(0, 1);
					v[2 * i][2 * j + 1] = m.v[i][j] * c * Complex(0, 1);
				}
			}
		}
	}
	m = Matrix(v);
	Complex det = m.determinant();
	return det.real;
}

/*
NOTE! Overflows happen rather quickly because the number of trials per experiment increases exponentially

INPUT FORMAT
============
epsilon delta time_limit
n
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
*/
int main() {
	srand(time(0));
	int n;
	long double epsilon, delta;
	long double time_limit;
	int repeats;
	int matrices;
	cin>>epsilon>>delta>>time_limit;
	cin>>n;
	
	vector<vector<Complex>> v;
	for (int i = 0; i < n; i++) {
		vector<Complex> w;
		for (int j = 0; j < n; j++) {
			long double d;
			cin>>d;
			w.push_back(Complex(d, 0));
		}
		v.push_back(w);
	}
	Matrix m(v);
	
	vector<long double> absvals;
	chrono::steady_clock::time_point start_time = chrono::steady_clock::now();			

	#if TYPE == 0
	long double critical_ratio = 3.0;
	#elif TYPE == 1
	long double critical_ratio = 2.0;
	#elif TYPE == 2
	long double critical_ratio = 1.5;
	#endif

	long long N = ceil(4.0 / epsilon / epsilon * pow(critical_ratio, n / 2.0)) + 0.05;
	long long r = ceil(21.156209733527273 * log(1 / delta)) + 0.05; // 21.15 from Chernoff bounds
	if (r % 2 == 0) r++;
	// cout<<"Repeating an experiment of "<<N<<" trials "<<r<<" times..."<<endl;
	for (int i = 0; i < r; i++) {
		AbsLog sum(0);
		long double trial = 0;
		for (int j = 0; j < N; j++) {
			trial += 1;
			
			#if TYPE == 0
			sum = sum + godsil(m);
			#elif TYPE == 1
			sum = sum + karmarkar(m);
			#elif TYPE == 2
			sum = sum + chien(m);
			#endif

			// cout<<((sum.value - log(trial)))<<endl;
			if (time_since(start_time) > time_limit) {
				break;
			}
		}
		absvals.push_back(sum.value - log(trial));
		if (time_since(start_time) > time_limit) {
			break;
		}
	}
	sort(absvals.begin(), absvals.end());
	// OUTPUT: estimate running_time
	cout<<absvals[r / 2]<<" "<<time_since(start_time)<<endl;
}