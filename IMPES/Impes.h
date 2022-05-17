#pragma once
#include <iostream>
#include <vector>

using namespace std;
class Impes {
protected:
	int numberOfPointByT, numberOfPointByX, numberOfPointByZ;
	int stepByT, stepByX, stepByZ;
	double a, b, c, d;
	double (*pressureTEqual0)(double, double);
	double (*saturationTEqual0)(double, double);

	double (*fa)(double, double), (*fb)(double, double), (*fc)(double, double), (*fd)(double, double);
	double (*alfa1)(double, double, double), (*alfa2)(double, double, double), (*alfa3)(double, double, double), (*alfa4)(double, double, double), 
		(*beta1)(double, double, double), (*beta2)(double, double, double), (*beta3)(double, double, double), (*beta4)(double, double, double);

	vector<vector<vector<double>>> saturation, pressure, k0, kw, ko, phasePermeability,	permeability;

public:
	Impes() { //individual options //almost all options have to variation in[0, 1]. For example x, z, sigma pressure and etc //grid will be 1x1x1
		this->numberOfPointByT = 1e+1;
		this->numberOfPointByX = 1e+1;
		this->numberOfPointByZ = 1e+1;

		setInitialBoundaryConditions();


		this->stepByT = 1 / (numberOfPointByT - 1);
		this->stepByX = 1 / (numberOfPointByX - 1);
		this->stepByZ = 1 / (numberOfPointByZ - 1);

		this->saturation = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->pressure = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		
		this->k0 = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->kw = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->ko = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->phasePermeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->permeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
	
		
	
	}

	void Calculate(double t) {
		//double omega = 
	}

	double omega(int n, int i, int j) {
		//return 1 / m * (Kw)
	}

	void setInitialBoundaryConditions() {
		this->pressureTEqual0 = [](double x, double z) {
			return x + z;
		};
		this->saturationTEqual0 = [](double x, double z) {
			return x + z;
		};



		this->fa = [](double t, double z) {
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

