#pragma once
#include <iostream>
#include <vector>

using namespace std;
class Impes {
protected:
	int numberOfPointByT, numberOfPointByX, numberOfPointByZ;
	double (*pressureTEqual0)(double, double);
	double (*saturationTEqual0)(double, double);
	double (*alfa1)(double, double, double), (*alfa2)(double, double, double), (*alfa3)(double, double, double), (*alfa4)(double, double, double), 
		(*beta1)(double, double, double), (*beta2)(double, double, double), (*beta3)(double, double, double), (*beta4)(double, double, double);

	vector<vector<vector<double>>> saturation, pressure, k0, phasePermeability,	permeability;

public:
	Impes() { // individual options
		this->numberOfPointByT = 1e+1;
		this->numberOfPointByX = 1e+1;
		this->numberOfPointByZ = 1e+1;

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


		this->saturation = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->pressure = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		
		this->k0 = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->phasePermeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->permeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
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

