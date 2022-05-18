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

	vector<vector<vector<double>>> saturationW, p, k0, kw, ko, phasePermeability, k;

public:
	Impes() { //individual options //almost all options have to variation in [0, 1]. For example x, z, sigma, pressure and etc //grid will be 1x1x1
		this->numberOfPointByT = 1e+1;
		this->numberOfPointByX = 1e+1;
		this->numberOfPointByZ = 1e+1;

		this->m = 1.;
		this->muo = 1.;
		this->muw = 1.;

		this->bw = 1.;
		this->bo = 1.;
		this->bc = 1.;
		this->Bw = bw + bc;
		this->Bo = bo + bc;



	
		setInitialBoundaryConditions();
	
	}

	void Calculate() {






		for (int n = 0; n < numberOfPointByT - 1; n++) { //is defined without boundaries!!!
			for (int i = 1; i < numberOfPointByX - 1; i++) {
				for (int j = 1; j < numberOfPointByZ - 1; j++) {
					saturationW[n + 1][i][j] = saturationW[n][i][j] + omega(n, i, j) * stepByT;
				}
			}
			setKwAndKo(n + 1);
		}

		//printArray(saturationW);
		printArray(p);
	}

	double omega(int n, int i, int j) { 
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		double meanSigmaXIP = (sigmaX(n, i, j, "w") + sigmaX(n, i + 1, j, "w")) / 2.;
		double meanSigmaXIN = (sigmaX(n, i - 1, j, "w") + sigmaX(n, i, j, "w")) / 2.;

		double meanSigmaXJP = (sigmaX(n, i, j, "w") + sigmaX(n, i, j + 1, "w")) / 2.;
		double meanSigmaXJN = (sigmaX(n, i, j - 1, "w") + sigmaX(n, i, j, "w")) / 2.;

		double answer = 1. / m * (
			1. / (stepByX * stepByX) * (meanSigmaXIP * p[n + 1][i + 1][j] - (meanSigmaXIP + meanSigmaXIN) * p[n + 1][i][j] + meanSigmaXIN * p[n + 1][i - 1][j]) +
			1. / (stepByZ * stepByZ) * (meanSigmaXJP * p[n + 1][i][j + 1] - (meanSigmaXJP + meanSigmaXJN) * p[n + 1][i][j] + meanSigmaXJN * p[n + 1][i][j - 1]) -
			N(n, i, j, "w") - B(n, i, j, "w") * saturationW[n][i][j] * (p[n + 1][i][j] - p[n][i][j]) / stepByT
			);

		return answer;
	}

	double sigmaX(int n, int i, int j, string phase = "general") {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		int numberPreLastPointInArrayX = numberOfPointByX - 1 - 1;
		int numberPreLastPointInArrayZ = numberOfPointByZ - 1 - 1;

		int I = i, J = j;

		if (I > numberPreLastPointInArrayX) {
			I = numberPreLastPointInArrayX;
		} else if (I == 0) {
			I = 1;
		}

		if (J > numberPreLastPointInArrayZ) {
			J = numberPreLastPointInArrayZ;
		} else if (J == 0) {
			J = 1;
		}

		if (phase == "w") {
			return k0[n][I][J] * kw[n][I][J] / muw;
		} else if (phase == "o") {
			return k0[n][I][J] * ko[n][I][J] / muo;
		} else {
			return k0[n][I][J] * kw[n][I][J] / muw + k0[n][I][J] * ko[n][I][J] / muo;
			//cerr << "PHASE ERROR" << endl;
			//throw "PHASE ERROR";
		}

	}

	double sigmaZ(int n, int i, int j, string phase = "general") {
		return sigmaX(n, i, j, phase) / 10.;
	}

	double N(int n, int i, int j, string phase = "general") {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "w") {
			return 1.;
		} else if (phase == "o") {
			return 1.;
		} else if (phase == "general") {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	double B(int n, int i, int j, string phase = "general") {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "w") {
			return Bw;
		} else if (phase == "o") {
			return Bo;
		} else if (phase == "general") {
			return Bw * saturationW[n][i][j] + Bo * (1 - saturationW[n][i][j]);
			//cerr << "PHASE ERROR" << endl;
			//throw "PHASE ERROR";
		}

	}

	void setKwAndKo(int n) {
		double cw = 0.1, co = 0.3, ak = 1., bk = 1.;

		for (int i = 0; i < numberOfPointByX; i++) {
			for (int j = 0; j < numberOfPointByZ; j++) {
				double x = stepByX * i;
				double z = stepByZ * j;

				double s = saturationW[n][i][j];

				kw[n][i][j] = (s - cw) * (1. - ak * (1. - s));
				ko[n][i][j] = (1. - s - co) * (1. - bk * s);

				//k[n][i][j] = k0[n][i][j] * ������� �������������(������������) // must be filled
			}
		}
	}


	void setInitialBoundaryConditions() {
		this->stepByT = 1. / (numberOfPointByT - 1.);
		this->stepByX = 1. / (numberOfPointByX - 1.);
		this->stepByZ = 1. / (numberOfPointByZ - 1.);

		this->saturationW = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->p = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));

		this->k0 = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ, 1)));
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
				saturationW[0][i][j] = saturationTEqual0(x, z);
			}
		}

		setKwAndKo(0);

		


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

	vector<double> SeidelMethod(const vector<vector<double>>& matrix1, const vector<double>& matrix2, double eps) {
		// ��������� ������ �������� �������
		int size = matrix1.size();

		// ����� ������� ������� � �������, ��������� �� 
		// �������� ������������ �����
		vector <vector<double>> matrix = matrix1;

		//// ������� ����� ����� ������ (size) x (size + 1),
		//// c ������ ������� ��������� ������    
		//matrix.resize(size);
		//for (int i = 0; i < size; i++) {
		//	matrix[i].resize(size + 1);

		//	for (int j = 0; j < size + 1; j++) {
		//		cin >> matrix[i][j];
		//	}
		//}


		// ������ ������ �������� ����������� �� ���������� ��������,
		// ������ �������� ����� ����� ����� � �������, �.�. size,
		// ������ �������� ������ ���������� ��������� ��� ������
		vector <double> previousVariableValues(size, 0.0);

		// ����� ��������� ������������ ������� �� ��� ���, 
		// ���� �� ����� ���������� ����������� ��������    
		while (true) {
			// ������ ������ �������� ����������� �� ������� ����       
			vector <double> currentVariableValues(size);

			// ��������� �������� ����������� �� ������� ��������
			// � ������������ � �������������� ���������
			for (int i = 0; i < size; i++) {
				// �������������� i-�� ����������� ��������� 
				// ���������� ����� i-�� ������ �������
				currentVariableValues[i] = matrix2[i];

				// �������� ����� �� ���� �������� �� i-�� �����������
				for (int j = 0; j < size; j++) {
					// ��� j < i ����� ������������ ��� �����������
					// �� ���� �������� �������� �����������
					if (j < i) {
						currentVariableValues[i] -= matrix[i][j] * currentVariableValues[j];
					}

					// ��� j > i ���������� �������� � ������� ��������
					if (j > i) {
						currentVariableValues[i] -= matrix[i][j] * previousVariableValues[j];
					}
				}

				// ����� �� ����������� ��� i-�� �����������
				currentVariableValues[i] /= matrix[i][i];
			}

			// ��������� ������� ����������� ������������ ���������� ��������
			double error = 0.0;

			for (int i = 0; i < size; i++) {
				error += abs(currentVariableValues[i] - previousVariableValues[i]);
			}

			// ���� ����������� �������� ����������, �� ��������� �������
			if (error < eps) {
				break;
			}

			// ��������� � ��������� ��������, ��� 
			// ��� ������� �������� ����������� 
			// ���������� ���������� �� ���������� ��������
			previousVariableValues = currentVariableValues;
		}

		//// ������� ��������� �������� ����������� � 8 ������� ��������
		//for (int i = 0; i < size; i++) {
		//	printf("%.8llf ", previousVariableValues[i]);
		//}

		return previousVariableValues;
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
		cout << matrix1.size() << "x1" << endl;
		for (int i = 0; i < matrix1.size(); i++) {
			cout << matrix1[i] << endl;
		}

	}




};

