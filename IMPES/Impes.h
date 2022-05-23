#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>

using namespace std;
class Impes {

protected:
	int numberOfPointByT, numberOfPointByX, numberOfPointByZ;
	double stepByT, stepByX, stepByZ;
	double a, b, c, d, t;
	double m, muw, muo, Bw, Bo, bw, bo, bc;
	double Lx, Lz, p0, B0, sigmaX0, sigmaZ0, sigmaZW0, sigmaXW0, A, n_;
	double (*pressureTEqual0)(double, double);
	double (*saturationTEqual0)(double, double);

	double (*fa)(double, double), (*fb)(double, double), (*fc)(double, double), (*fd)(double, double);
	double (*alfa1)(double, double, double), (*alfa2)(double, double, double), (*alfa3)(double, double, double), (*alfa4)(double, double, double), 
		(*beta1)(double, double, double), (*beta2)(double, double, double), (*beta3)(double, double, double), (*beta4)(double, double, double);

	vector<vector<vector<double>>> saturationW, p, k0, kw, ko, phasePermeability, k;


	string pointSchedule = "with points pointtype 5";
	string lineSchedule = "with lines lw 3";
	string path = "\"C:/Users/student/Desktop/albert/coursework6/IMPES";

public:
	Impes() { //individual options //almost all options have to variation in [0, 1]. For example x, z, sigma, pressure and etc //grid will be 1x1x1
		this->numberOfPointByT = 100;
		this->numberOfPointByX = 10;
		this->numberOfPointByZ = 10;

		this->m = .24;
		this->muo = 40.;
		this->muw = 9.;

		this->bw = 3.7e-10;
		this->bo = 7.4e-9;
		this->bc = 4.5e-10;
		this->Bw = bw + bc;
		this->Bo = bo + bc;

		this->a = 0;
		this->c = 0;
		this->b = 100;
		this->d = 10;
		this->t = 1;

		this->Lx = 200;
		this->Lz = 10;

		


	
		setInitialBoundaryConditions();
	
	}

	void Calculate() {
		
		for (int n = 0; n < numberOfPointByT - 1; n++) { //is defined without boundaries!!!

			//pressure
			int numberOfVariablesOnX = (numberOfPointByZ - 2);
			int numberOfVariablesOnZ = (numberOfPointByX - 2);
			int numberOfVariables = numberOfVariablesOnX  * numberOfVariablesOnZ;
			vector<vector<double>> matrix1(numberOfVariables, vector<double>(numberOfVariables, 1));
			vector<double> matrix2(numberOfVariables, 0);

			for (int i = 1; i < numberOfPointByX - 1; i++) { 
				for (int j = 1; j < numberOfPointByZ - 1; j++) {
					double meanSigmaXIP = (sigmaX(n, i, j) + sigmaX(n, i + 1, j)) / 2.;
					double meanSigmaXIN = (sigmaX(n, i - 1, j) + sigmaX(n, i, j)) / 2.;

					double meanSigmaXJP = (sigmaX(n, i, j) + sigmaX(n, i, j + 1)) / 2.;
					double meanSigmaXJN = (sigmaX(n, i, j - 1) + sigmaX(n, i, j)) / 2.;


					double meanSigmaZIP = (sigmaZ(n, i, j) + sigmaZ(n, i + 1, j)) / 2.;
					double meanSigmaZIN = (sigmaZ(n, i - 1, j) + sigmaZ(n, i, j)) / 2.;

					double meanSigmaZJP = (sigmaZ(n, i, j) + sigmaZ(n, i, j + 1)) / 2.;
					double meanSigmaZJN = (sigmaZ(n, i, j - 1) + sigmaZ(n, i, j)) / 2.;

					vector<double> checkPoint(5, 0);
					double ak = D(n, i - 1, j) * stepByT / (stepByX * stepByX) * meanSigmaXIN;
					double bk = A * D(n, i, j - 1) * stepByT / (stepByZ * stepByZ) * meanSigmaZJN;

					double ck = D(n, i, j) * stepByT / (stepByX * stepByX) * (meanSigmaXIP + meanSigmaXIN) + 
								D(n, i, j) * stepByT / (stepByZ * stepByZ) * A * (meanSigmaZJP + meanSigmaZJN) + 1;

					double dk = A * D(n, i, j + 1) * stepByT / (stepByZ * stepByZ) * meanSigmaZJP;
					double ek = D(n, i + 1, j) * stepByT / (stepByX * stepByX) * meanSigmaXIP;

					double hk = n_ * N(n, i, j) * D(n, i, j) * stepByT - p[n][i][j];

					if (j == 1) {
						checkPoint[1] = 0;
						hk -= bk * p[n][i][j - 1];
					} else if (j == numberOfPointByZ - 2) {
						checkPoint[3] = 0;
						hk -= dk * p[n][i][j + 1];
					}

					if (i == 1) {
						checkPoint[0] = 0;
						hk -= ak * p[n][i - 1][j];
					} else if (i == numberOfPointByX - 2) {
						checkPoint[4] = 0;
						hk -= ek * p[n][i + 1][j];
					}


					int currentPositionX = (i - 1) * numberOfVariablesOnX;
					int currentPositionZ = (j - 1);
					int currentPosition = currentPositionX + currentPositionZ;

					matrix2[currentPosition] = hk;


					matrix1[currentPosition][currentPosition] += ck;
					if (j != 1)
						matrix1[currentPosition][currentPosition - 1] += bk;
					if (j != numberOfPointByZ - 2)
						matrix1[currentPosition][currentPosition + 1] += dk;
					if (i != 1)
						matrix1[currentPosition][currentPosition - numberOfVariablesOnX] += ak;
					if (i != numberOfPointByX - 2)
						matrix1[currentPosition][currentPosition + numberOfVariablesOnX] += ek;
				}
			}

			//printArray(matrix1);
			cout << "LINE = " << __LINE__ << endl;

			vector<double> pAnswer = SeidelMethod(matrix1, matrix2, 0.01);

			//printArray(pAnswer);
			cout << "n = " << n << endl;

			for (int i = 1; i < numberOfPointByX - 1; i++) {
				for (int j = 1; j < numberOfPointByZ - 1; j++) {
					//int currentPositionX = (i - 1) * numberOfVariablesOnX;
					//int currentPositionZ = (j - 1);
					//int currentPosition = currentPositionX + currentPositionZ;

					p[n + 1][i][j] = pAnswer[(i - 1) * numberOfVariablesOnX + (j - 1)];
				}
			}



			//saturation
			for (int i = 1; i < numberOfPointByX - 1; i++) { 
				for (int j = 1; j < numberOfPointByZ - 1; j++) {
					saturationW[n + 1][i][j] = saturationW[n][i][j] + omega(n, i, j) * stepByT;
				}
			}
			setKwAndKo(n + 1);
		}

		printArray(saturationW);
		//printArray(p);
	}

	double omega(int n, int i, int j) { 
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		double meanSigmaXIP = (sigmaX(n, i, j, "w") + sigmaX(n, i + 1, j, "w")) / 2.;
		double meanSigmaXIN = (sigmaX(n, i - 1, j, "w") + sigmaX(n, i, j, "w")) / 2.;

		double meanSigmaXJP = (sigmaX(n, i, j, "w") + sigmaX(n, i, j + 1, "w")) / 2.;
		double meanSigmaXJN = (sigmaX(n, i, j - 1, "w") + sigmaX(n, i, j, "w")) / 2.;

	   
		double meanSigmaZIP = (sigmaZ(n, i, j, "w") + sigmaZ(n, i + 1, j, "w")) / 2.;
		double meanSigmaZIN = (sigmaZ(n, i - 1, j, "w") + sigmaZ(n, i, j, "w")) / 2.;

		double meanSigmaZJP = (sigmaZ(n, i, j, "w") + sigmaZ(n, i, j + 1, "w")) / 2.;
		double meanSigmaZJN = (sigmaZ(n, i, j - 1, "w") + sigmaZ(n, i, j, "w")) / 2.;

		double E1 = (meanSigmaXIP * p[n + 1][i + 1][j] - (meanSigmaXIP + meanSigmaXIN) * p[n + 1][i][j] + meanSigmaXIN * p[n + 1][i - 1][j]);
		double E2 = (meanSigmaZJP * p[n + 1][i][j + 1] - (meanSigmaZJP + meanSigmaZJN) * p[n + 1][i][j] + meanSigmaZJN * p[n + 1][i][j - 1]);
		double diffPressure = (p[n + 1][i][j] - p[n][i][j]);

		double valueB = -m * B(n, i, j, "w") * saturationW[n][i][j] * diffPressure / stepByT;
		double valueN = N(n, i, j, "w");

		double answer = B0 * p0 / m * (
			1. / (stepByX * stepByX) * E1 +
			1. / (stepByZ * stepByZ) * A * E2 -
			n_ * N(n, i, j, "w") - m * B(n, i, j, "w") * saturationW[n][i][j] * diffPressure / stepByT
			);

		//cout << answer << endl;
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
		} else if (phase == "general") {
			return k0[n][I][J] * kw[n][I][J] / muw + k0[n][I][J] * ko[n][I][J] / muo;
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	double sigmaZ(int n, int i, int j, string phase = "general") {
		return sigmaX(n, i, j, phase) / 10.;
	}

	double N(int n, int i, int j, string phase = "general") {
		double coeff = 1.5;
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "w") {
			return -(coeff);
		}
		else if (phase == "o") {
			return -(coeff);
		}
		else if (phase == "general") {
			return -(coeff + coeff);
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	double B(int n, int i, int j, string phase = "general") {
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
			return Bw;
		} else if (phase == "o") {
			return Bo;
		} else if (phase == "general") {
			return Bw * saturationW[n][I][J] + Bo * (1 - saturationW[n][I][J]);
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
		}

	}

	double D(int n, int i, int j, string phase = "general") {
		double t = stepByT * n;
		double x = stepByX * i;
		double z = stepByZ * j;

		if (phase == "general") {
			return 1. / (m * B(n, i, j));
		} else {
			cerr << "PHASE ERROR" << endl;
			throw "PHASE ERROR";
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

				//k[n][i][j] = k0[n][i][j] * фазовая проницаемость(насыщенность) // must be filled
			}
		}
	}

	void printSaturationWPermeabilityW() {
		int layerOfZ = numberOfPointByZ - 2;
		vector<double> first(numberOfPointByT);
		vector<double> second(numberOfPointByT);

		for (int i = 0; i < numberOfPointByT; i++) {
			first[i] = i;
			second[i] = saturationW[i][numberOfPointByX - 2][layerOfZ];
		}
		makeFileForGraph(first, second, "printSaturationWChangesOverTime.txt");
		drawGraph("printSaturationWChangesOverTime.txt", "graph");
		//printArray(first);
	}

	void setInitialBoundaryConditions() {
		this->stepByT = 1 / (numberOfPointByT - 1.);
		this->stepByX = (b - a) / (numberOfPointByX - 1.);
		this->stepByZ = (d - c) / (numberOfPointByZ - 1.);

		this->saturationW = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->p = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));

		this->k0 = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ, 1e-12)));
		this->kw = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->ko = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		this->k = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));
		//this->phasePermeability = vector<vector<vector<double>>>(numberOfPointByT, vector<vector<double>>(numberOfPointByX, vector<double>(numberOfPointByZ)));

		//this->X0 = k0[0][0][0] / mu;

		this->pressureTEqual0 = [](double x, double z) { //must be added
			return 0.7 * 1.78e+7;
		};
		this->saturationTEqual0 = [](double x, double z) { //must be added
			return 0.4;
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

		

		//must be added
		this->fa = [](double t, double z) { //let these functions set the boundary conditions for the pressure
			return 0.7 * 1.78e+7;
		};
		this->fb = [](double t, double z) {
			return 0.7 * 1.78e+7;
		};
		this->fc = [](double t, double x) { 
			return 1.78e+7;
			//return 1.;
		};
		this->fd = [](double t, double x) {
			return 0.5 * 1.78e+7;
			//return 1.;
		};

		
		for (int n = 1; n < numberOfPointByT; n++) { //definition of pressure and saturation on boundaries
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

					saturationW[n][i][0] = 1.;
					//saturationW[n][i][numberOfPointByZ - 1] = 0.4;
				}
			
		}

		this->sigmaX0 = 0;
		this->sigmaXW0 = 0;

		this->sigmaZ0 = 0;
		this->sigmaZW0 = 0;

		this->p0 = 0;
		this->B0 = 0;


			for (int i = 0; i < numberOfPointByX; i++) {
				for (int j = 0; j < numberOfPointByZ; j++) {
					if (B(0, i, j) > B0)
						B0 = B(0, i, j);
					if (sigmaX(0, i, j) > sigmaX0)
						sigmaX0 = sigmaX(0, i, j);
					if (sigmaZ(0, i, j) > sigmaZ0)
						sigmaZ0 = sigmaZ(0, i, j);
					if (p[0][i][j] > p0)
						p0 = p[0][i][j];
				}
			}

			//this->A = Lx * Lx * sigmaZ0 / (Lz * Lz * sigmaX0);
			//this->n_ = Lx * Lx / (p0 * sigmaX0);
			//cout << "coeffs have been added" << endl;

			this->A = 1;
			this->n_ = 1;
	}

	vector<double> SeidelMethod(const vector<vector<double>>& matrix1, const vector<double>& matrix2, double eps) {
		// Считываем размер вводимой матрицы
		int size = matrix1.size();

		// Будем хранить матрицу в векторе, состоящем из 
		// векторов вещественных чисел
		vector <vector<double>> matrix = matrix1;

		// Введем вектор значений неизвестных на предыдущей итерации,
		// размер которого равен числу строк в матрице, т.е. size,
		// причем согласно методу изначально заполняем его нулями
		vector <double> previousVariableValues(size, 0.0);

		// Будем выполнять итерационный процесс до тех пор, 
		// пока не будет достигнута необходимая точность    
		while (true) {
			// Введем вектор значений неизвестных на текущем шаге       
			vector <double> currentVariableValues(size);

			// Посчитаем значения неизвестных на текущей итерации
			// в соответствии с теоретическими формулами
			for (int i = 0; i < size; i++) {
				// Инициализируем i-ую неизвестную значением 
				// свободного члена i-ой строки матрицы
				currentVariableValues[i] = matrix2[i];

				// Вычитаем сумму по всем отличным от i-ой неизвестным
				for (int j = 0; j < size; j++) {
					// При j < i можем использовать уже посчитанные
					// на этой итерации значения неизвестных
					if (j < i) {
						currentVariableValues[i] -= matrix[i][j] * currentVariableValues[j];
					}

					// При j > i используем значения с прошлой итерации
					if (j > i) {
						currentVariableValues[i] -= matrix[i][j] * previousVariableValues[j];
					}
				}

				// Делим на коэффициент при i-ой неизвестной
				currentVariableValues[i] /= matrix[i][i];
			}

			// Посчитаем текущую погрешность относительно предыдущей итерации
			double error = 0.0;

			for (int i = 0; i < size; i++) {
				error += abs(currentVariableValues[i] - previousVariableValues[i]);
			}

			// Если необходимая точность достигнута, то завершаем процесс
			if (error < eps) {
				break;
			}

			// Переходим к следующей итерации, так 
			// что текущие значения неизвестных 
			// становятся значениями на предыдущей итерации
			previousVariableValues = currentVariableValues;
		}

		return previousVariableValues;
	}

	void drawGraph(string name, string text) {
		string str;
		str = path + name + "\" using 1:2 title \"" + text + "\" " + pointSchedule + ";";
		std::ofstream graphic("file");
		cout << str << endl;
		//graphic << "plot " << str << "; pause mouse keypress" << "\n";
		graphic << "set border 3" << endl;
		graphic << "plot " << str << endl;
		graphic.close();
		system("gnuplot -persist file");

		graphic.close();
	}

	template <typename T, typename K>
	void makeFileForGraph(const vector<T>& first, const vector<K>& second, const string fileName) {
		//cout << "file" << endl;

		ofstream file(fileName);

		assert(first.size() == second.size() && "Size mismatch");

		for (int i = 0; i < first.size(); i++) {
			//cout << "file = " << i << endl;

			file << first[i] << " " << second[i] << endl;
		}


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

