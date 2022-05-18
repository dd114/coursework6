#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

using namespace std;
class Impes {
protected:
	int numberOfPointByT, numberOfPointByX, numberOfPointByZ;
	double stepByT, stepByX, stepByZ;
	double a, b, c, d;
	double m, muw, muo, Bw, Bo, bw, bo, bc;
	double (*pressureTEqual0)(double, double);
	double (*saturationTEqual0)(double, double);

	double (*fa)(double, double), (*fb)(double, double), (*fc)(double, double), (*fd)(double, double);
	double (*alfa1)(double, double, double), (*alfa2)(double, double, double), (*alfa3)(double, double, double), (*alfa4)(double, double, double), 
		(*beta1)(double, double, double), (*beta2)(double, double, double), (*beta3)(double, double, double), (*beta4)(double, double, double);

	vector<vector<vector<double>>> saturation, p, k0, kw, ko, phasePermeability, k;

public:
	Impes() { //individual options //almost all options have to variation in[0, 1]. For example x, z, sigma, pressure and etc //grid will be 1x1x1
		this->numberOfPointByT = 1e+1;
		this->numberOfPointByX = 1e+1;
		this->numberOfPointByZ = 1e+1;

		this->m = 1.;
		this->muo = 1.;
		this->muw = 1.;

		this->bw = 1.;
		this->bc = 1.;
		this->Bw = bw + bc;
		this->Bo = bo + bc;




	
		setInitialBoundaryConditions();
	
	}

	void Calculate() {
		for (int n = 0; n < numberOfPointByT - 1; n++) { //is defined without boundaries!!!
			for (int i = 1; i < numberOfPointByX - 1; i++) {
				for (int j = 1; j < numberOfPointByZ - 1; j++) {
					saturation[n + 1][i][j] = saturation[n][i][j] + omega(n, i, j) * stepByT;
				}
			}
		}

		//printArray(saturation);
		printArray(p);
	}

	double omega(int n, int i, int j) { 
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		double meanSigmaXIP = (sigmaX("w", n, i, j) + sigmaX("w", n, i + 1, j)) / 2.;
		double meanSigmaXIN = (sigmaX("w", n, i - 1, j) + sigmaX("w", n, i, j)) / 2.;

		double meanSigmaXJP = (sigmaX("w", n, i, j) + sigmaX("w", n, i, j + 1)) / 2.;
		double meanSigmaXJN = (sigmaX("w", n, i, j - 1) + sigmaX("w", n, i, j)) / 2.;
		
		return 1. / m * (
			1. / (stepByX * stepByX) * (meanSigmaXIP * p[n + 1][i + 1][j] - (meanSigmaXIP + meanSigmaXIN) * p[n + 1][i][j] + meanSigmaXIN * p[n + 1][i - 1][j]) +
			1. / (stepByZ * stepByZ) * (meanSigmaXJP * p[n + 1][i][j + 1] - (meanSigmaXJP + meanSigmaXJN) * p[n + 1][i][j] + meanSigmaXJN * p[n + 1][i][j - 1]) -
			N("w", n, i, j) - Bw * saturation[n][i][j] * (p[n + 1][i][j] - p[n][i][j]) / stepByT
						 );
	}

	double sigmaX(string phase, int n, int i, int j) {
		return sigma(phase, n, i, j);
	}

	double sigmaZ(string phase, int n, int i, int j) {
		return sigma(phase, n, i, j);
	}

	double sigma(string phase, int n, int i, int j) {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "w") {
			return kw[n][i][j] / muw;
		} else if (phase == "o") {
			return ko[n][i][j] / muo;
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	double N(string phase, int n, int i, int j) {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "w") {
			return 1.;
		} else if (phase == "o") {
			return 1.;
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	void setK(int n) {
		double cw = 0.1, co = 0.3, ak = 1., bk = 1.;

		for (int i = 0; i < numberOfPointByX; i++) {
			for (int j = 0; j < numberOfPointByZ; j++) {
				double x = stepByX * i;
				double z = stepByZ * j;

				double s = saturation[n][i][j];

				kw[n][i][j] = (s - cw) * (1. - ak * (1. - s));
				ko[n][i][j] = (1. - s - co) * (1. - bk * s);

				//k[n][i][j] = k0[n][i][j] * ������� �������������(������������) // must be filled
			}
		}
	}

	vector<double> tridiagonalSolution(const vector<vector<double>>& matrix1, const vector<double>& matrix2) {
		assert(matrix1.size() == matrix2.size() && "Sizes match");

		double y;
		int N = matrix2.size();
		int N1 = N - 1;
		vector<double> a(N), B(N), matRes(N);

		y = matrix1[0][0];
		a[0] = -matrix1[0][1] / y;
		B[0] = matrix2[0] / y;
		for (int i = 1; i < N1; i++) {
			y = matrix1[i][i] + matrix1[i][i - 1] * a[i - 1];
			a[i] = -matrix1[i][i + 1] / y;
			B[i] = (matrix2[i] - matrix1[i][i - 1] * B[i - 1]) / y;
		}

		matRes[N1] = (matrix2[N1] - matrix1[N1][N1 - 1] * B[N1 - 1]) / (matrix1[N1][N1] + matrix1[N1][N1 - 1] * a[N1 - 1]);
		for (int i = N1 - 1; i >= 0; i--) {
			matRes[i] = a[i] * matRes[i + 1] + B[i];
		}

		return matRes;

	}

	void setInitialBoundaryConditions() {
		this->stepByT = 1. / (numberOfPointByT - 1.);
		this->stepByX = 1. / (numberOfPointByX - 1.);
		this->stepByZ = 1. / (numberOfPointByZ - 1.);

		this->saturation = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->p = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));

		this->k0 = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->kw = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->ko = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->k = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		//this->phasePermeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));



		this->pressureTEqual0 = [](double x, double z) {
			return x + z;
		};
		this->saturationTEqual0 = [](double x, double z) {
			return x + z;
		};
		
		for (int i = 0; i < numberOfPointByX; i++) {
			for (int j = 0; j < numberOfPointByZ; j++) {
				double x = stepByX * i;
				double z = stepByZ * j;

				p[0][i][j] = pressureTEqual0(x, z);
				saturation[0][i][j] = saturationTEqual0(x, z);
			}
		}
		


		this->fa = [](double t, double z) { //let these functions set the boundary conditions for the pressure
			return 1.;
		};
		this->fb = [](double t, double z) {
			return 1.;
		};
		this->fc = [](double t, double x) {
			return 1.;
		};
		this->fd = [](double t, double x) {
			return 1.;
		};

		
		for (int n = 0; n < numberOfPointByT; n++) { //definition of pressure on boundaries
			double t = stepByT * n;

				for (int j = 0; j < numberOfPointByZ; j++) {
					double z = stepByZ * j;
					p[n][0][j] = fa(z, t);
					p[n][numberOfPointByX - 1][j] = fb(z, t);
				}

				for (int i = 0; i < numberOfPointByX; i++) {
					double x = stepByX * i;
					p[n][i][0] = fc(x, t);
					p[n][i][numberOfPointByZ - 1] = fd(x, t);
				}
			
		}

		this->alfa1 = [](double t, double x, double z) {
			return 1.;
		};
		this->alfa2 = [](double t, double x, double z) {
			return 1.;
		};
		this->alfa3 = [](double t, double x, double z) {
			return 1.;
		};
		this->alfa4 = [](double t, double x, double z) {
			return 1.;
		};

		this->beta1 = [](double t, double x, double z) {
			return 1.;
		};
		this->beta2 = [](double t, double x, double z) {
			return 1.;
		};
		this->beta3 = [](double t, double x, double z) {
			return 1.;
		};
		this->beta4 = [](double t, double x, double z) {
			return 1.;
		};


	}

	virtual void setNumberOfPointByT(int numberOfPointByT) {
		this->numberOfPointByT = numberOfPointByT;
	}

	virtual void setNumberOfPointByX(int numberOfPointByX) {
		this->numberOfPointByX = numberOfPointByX;
	}

	virtual void setNumberOfPointByZ(int numberOfPointByZ) {
		this->numberOfPointByZ = numberOfPointByZ;
	}

	template <typename T>
	void printArray(const vector<vector<vector<T>>>& matrix1) {
		cout << matrix1.size() << "x" << matrix1[0].size() << "x" << matrix1[0][0].size() << endl;
		for (int i = 0; i < matrix1.size(); i++) {

			cout << "LAYER = " << i << endl;
			for (int j = 0; j < matrix1[i].size(); j++) {
				for (int k = 0; k < matrix1[i][j].size(); k++) {
					cout << matrix1[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}

	}

	template <typename T>
	void printArray(const vector<vector<T>>& matrix1) {
		cout << matrix1.size() << "x" << matrix1[0].size() << endl;
		for (int i = 0; i < matrix1.size(); i++) {
			for (int j = 0; j < matrix1[i].size(); j++) {
				cout << matrix1[i][j] << " ";
			}
			cout << endl;
		}

	}

	template <typename T>
	void printArray(const vector<T>& matrix1) {
		cout << matrix1.size() << endl;
		for (int i = 0; i < matrix1.size(); i++) {
			cout << matrix1[i] << endl;
		}

	}




};

