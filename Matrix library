#include <iostream>
#include <math.h>
using namespace std;

class matrix {
public:
	matrix() {
		width = 0;
		height = 0;
		values = new double[0];
	}
	friend istream& operator>> (istream& in, matrix& matrix);
	friend ostream& operator<< (ostream& in, const matrix& matrix);
	friend matrix operator+ (matrix const &A, matrix const &B);
	friend matrix operator- (matrix const &A, matrix const &B);
	friend matrix operator* (matrix const &A, matrix const &B);
	friend matrix operator* (double k, matrix const &A);
	friend matrix operator* (matrix const &A, double k);
	void set(int w, int h, double* val) {
		width = w;
		height = h;
		delete[] values;
		values = new double[w*h];
		for (int i = 0; i < w*h; i++) {
			values[i] = val[i];
		}
	}
	double get_val(int i, int j) {
		return values[(j - 1)*width + (i - 1)];
	}
	matrix get_column(int i) {
		i--;
		matrix answer;
		double* val = new double[height];

		for (int j = 0; j < height; j++) {
			val[j] = values[j*width + i];
		}

		answer.set(1, height, val);

		return answer;
	}
	matrix get_stroke(int j) {
		j--;
		matrix answer;
		double* val = new double[width];

		for (int i = 0; i < width; i++) {
			val[i] = values[j*width + i];
		}

		answer.set(width, 1, val);

		return answer;
	}
	matrix minor(int i, int j) {
		matrix answer;

		int l = 0;
		double* val = new double[(width - 1)*(height - 1)];
		for (int k = 0; k < height; k++) {
			for (int p = 0; p < width; p++) {
				if (p + 1 == i || k + 1 == j) { continue; }
				val[l] = values[k*width + p];
				l++;
			}
		}

		answer.set(width - 1, height - 1, val);

		return answer;
	}
	matrix identity() {
		int size;
		matrix answer;
		if (width <= height) {
			size = width;
		}
		else {
			size = height;
		}
		double* val = new double[size*size];

		for (int j = 0; j < size; j++) {
			for (int i = 0; i < size; i++) {
				if (i == j) {
					val[j*size + i] = 1;
				}
				else {
					val[j*size + i] = 0;
				}
			}
		}

		answer.set(size, size, val);

		return answer;
	}
	double det() {
		double answer = 0;
		if (width != height) { return NULL; }
		if (width == 1) {
			answer = values[0];
		}
		else if (width == 2) {
			answer = values[0] * values[3] - values[1] * values[2];
		}
		else {
			for (int i = 0; i < width; i++) {
				answer += pow(-1, i)*get_val(i + 1, 1)*minor(i + 1, 1).det();
			}
		}
		return answer;
	}
	matrix ks_transvection(int a, double k) {
		matrix stroke, answer;

		double* val = new double[width*height];

		stroke = k*get_stroke(a);
		a--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (j == a) {
					val[j*width + i] = stroke.values[i];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix kc_transvection(int a, double k) {
		matrix column, answer;

		double* val = new double[width*height];

		column = k*get_column(a);
		a--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (i == a) {
					val[j*width + i] = column.values[j];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix swaps_transvection(int a, int b) {
		matrix stroke1, stroke2, answer;

		double* val = new double[width*height];

		stroke1 = get_stroke(a);
		stroke2 = get_stroke(b);
		a--; b--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (j == a) {
					val[j*width + i] = stroke2.values[i];
				}
				else if (j == b) {
					val[j*width + i] = stroke1.values[i];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix swapc_transvection(int a, int b) {
		matrix column1, column2, answer;

		double* val = new double[width*height];

		column1 = get_column(a);
		column2 = get_column(b);
		a--; b--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (i == a) {
					val[j*width + i] = column2.values[j];
				}
				else if (i == b) {
					val[j*width + i] = column1.values[j];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix s_transvection(int a, int b, double k) {
		matrix stroke, answer;

		double* val = new double[width*height];

		stroke = get_stroke(a) + k*get_stroke(b);
		a--; b--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (j == a) {
					val[j*width + i] = stroke.values[i];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix c_transvection(int a, int b, double k) {
		matrix column, answer;

		double* val = new double[width*height];

		column = get_column(a) + k*get_column(b);
		a--; b--;

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (i == a) {
					val[j*width + i] = column.values[j];
				}
				else {
					val[j*width + i] = values[j*width + i];
				}
			}
		}

		answer.set(width, height, val);

		return answer;

	}
	matrix transpose() {
		matrix A;
		A.set(width, height, values);

		double* val = new double[width*height];
		int k = 0;

		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				val[k] = A.values[j*width + i];
				k++;
			}
		}

		A.set(height, width, val);

		return A;
	}
	matrix alg_comp() {
		matrix A;
		A.set(width, height, values);
		double* val = new double[width*height];

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				val[j*width + i] = pow(-1, i + j)*A.minor(i + 1, j + 1).det();
			}
		}

		A.set(width, height, val);

		return A;
	}
	matrix power(int n) {
		matrix A, answer;
		A.set(width, height, values);
		if (n > 1 && width == height) {
			answer = A;
			for (int i = 1; i < n; i++) {
				answer = answer * A;
			}
		}
		else if (n == 1) {
			answer = A;
		}
		else if (n == 0 && width == height) {
			answer = A.identity();
		}
		else if (n == -1 && width == height && A.det() != 0) {
			answer = (1. / A.det())*A.alg_comp().transpose();
		}

		return answer;
	}
private:
	int width, height;
	double* values = new double[0];
};

istream& operator>> (istream& in, matrix& matrix) {
	cout << "---------------insert matrix" << endl;
	cout << "width: ";
	cin >> matrix.width;
	cout << "height: ";
	cin >> matrix.height;
	delete[] matrix.values;
	matrix.values = new double[matrix.width*matrix.height];

	for (int l = 0; l < matrix.width*matrix.height; l++) {
		cout << (int)(l / matrix.width) + 1 << "|" << l - (((int)(l / matrix.width))*matrix.width) + 1 << " : ";
		in >> matrix.values[l];
	}
	cout << "---------------end insert matrix" << endl;
	return in;
}

ostream& operator<< (ostream& os, const matrix& matrix) {
	for (int j = 0; j < matrix.height; j++) {
		for (int i = 0; i < matrix.width; i++) {
			os << matrix.values[matrix.width*j + i] << '\t';
		}
		os << endl;
	}

	return os;
}

matrix operator+ (matrix const &A, matrix const &B) {
	matrix answer;

	double* val = new double[A.width*A.height];

	for (int i = 0; i < A.width*A.height; i++) {
		val[i] = A.values[i] + B.values[i];
	}

	answer.set(A.width, A.height, val);

	return answer;
}

matrix operator- (matrix const &A, matrix const &B) {
	matrix answer;

	double* val = new double[A.width*A.height];

	for (int i = 0; i < A.width*A.height; i++) {
		val[i] = A.values[i] - B.values[i];
	}

	answer.set(A.width, A.height, val);

	return answer;
}

matrix operator* (matrix const &A, matrix const &B) {
	matrix answer;
	int w = B.width;
	int h = A.height;

	double* val = new double[w*h];
	int m = 0;

	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			for (int k = 0; k < A.width; k++) {
				if (k == 0) { val[m] = 0; }
				//val[m] += round(A.values[j*A.width + k] * B.values[k*B.width + i]*10000) / 10000;
				val[m] += A.values[j*A.width + k] * B.values[k*B.width + i];
			}
			m++;
		}
	}

	answer.set(w, h, val);

	return answer;
}

matrix operator* (double k, matrix const &A) {
	matrix answer;
	double* val = new double[A.width*A.height];

	for (int j = 0; j < A.height; j++) {
		for (int i = 0; i < A.width; i++) {
			//val[j*A.width + i] = round(k * A.values[j*A.width + i] * 10000) / 10000;
			val[j*A.width + i] = k * A.values[j*A.width + i];
		}
	}

	answer.set(A.width, A.height, val);

	return answer;
}

matrix operator* (matrix const &A, double k) {
	matrix answer;
	double* val = new double[A.width*A.height];

	for (int j = 0; j < A.height; j++) {
		for (int i = 0; i < A.width; i++) {
			//val[j*A.width + i] = round(A.values[j*A.width + i] * k * 10000) / 10000;
			val[j*A.width + i] = A.values[j*A.width + i] * k;
		}
	}

	answer.set(A.width, A.height, val);

	return answer;
}
