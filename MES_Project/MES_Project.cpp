#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <array>
#include <algorithm>

using namespace std;

void Xi(double* x, double* b, double** A, int m, int n)
{
	if (m < 1)
		return;

	x[m - 1] = b[m - 1];
	int j = n - 1;
	for (int i = 0; i < n - (m); i++)
	{
		x[m - 1] -= A[m - 1][n - 1 - i] * x[j];
		j--;
	}
	x[m - 1] /= A[m - 1][m - 1];
	Xi(x, b, A, m - 1, n);
}
double* ElimGauss(double** A, double* b, int n)
{
	double a11 = A[0][0], buf, buf2;
	double* x = new double[n];

	for (int k = 0; k < n; k++)
	{
		a11 = A[k][k];
		for (int j = k + 1; j < n; j++)
		{
			buf = A[j][k];
			for (int i = k; i < n; i++)
			{
				A[j][i] -= (buf * A[k][i]) / a11;
			}
			b[j] -= (buf * b[k]) / a11;
		}
	}

	Xi(x, b, A, n, n);

	return x;
}

double PochodnaKsi(double eta, int j) {
	double pochodna;
	switch (j) {
	case 1:
		pochodna = (eta - 1) / 4; //dla N1
		return pochodna;
		break;
	case 2:
		pochodna = (1 - eta) / 4; //dla N2
		return pochodna;
		break;
	case 3:
		pochodna = (eta + 1) / 4; //dla N3
		return pochodna;
		break;
	case 4:
		pochodna = (-eta - 1) / 4; //dla N4
		return pochodna;
		break;
	}
}
double PochodnaEta(double ksi, int j) {
	double pochodna;
	switch (j) {
	case 1:
		pochodna = (ksi - 1) / 4; //dla N1
		return pochodna;
		break;
	case 2:
		pochodna = (-ksi - 1) / 4; //dla N2
		return pochodna;
		break;
	case 3:
		pochodna = (ksi + 1) / 4; //dla N3
		return pochodna;
		break;
	case 4:
		pochodna = (1 - ksi) / 4; //dla N4
		return pochodna;
		break;
	}
}
vector<double> FunKsztaltu(double ksi, double eta) {
	vector<double> N(4);
	N[0] = 0.25 * (1 + ksi) * (1 + eta);
	N[1] = 0.25 * (1 - ksi) * (1 + eta);
	N[2] = 0.25 * (1 - ksi) * (1 - eta);
	N[3] = 0.25 * (1 + ksi) * (1 - eta);
	return N;
	//
	// N2----|------N1
	// |     |      |
	// |  Ele|m x   |
	//-|-----|------|-ksi
	// |     |      |
	// |     |      |
	// N3----|------N4
	//      eta
}
vector<double> FunKsztaltuPoOBJ(double ksi, double eta) {
	vector<double> N(4);
	N[0] = 0.25 * (1 - ksi) * (1 - eta);
	N[1] = 0.25 * (1 + ksi) * (1 - eta);
	N[2] = 0.25 * (1 + ksi) * (1 + eta);
	N[3] = 0.25 * (1 - ksi) * (1 + eta);
	return N;
	//
	// N4----|------N3
	// |     |      |
	// |  Ele|m x   |
	//-|-----|------|-ksi
	// |     |      |
	// |     |      |
	// N1----|------N2
	//      eta
}

double PochodnaDoJakobianu(double** val, vector<double> coord, int pc) { //oblicza pochodna dla podanego punktu całkowania, przyjmuje wartość eta lub ksi oraz współrzędne x lub y
	//cout << endl << "[X/Y] = " << coord[0] << " " << coord[1] << " " << coord[2] << " " << coord[3] << endl;
	//cout<<endl<<"dN/dKsilubEta = "<<val[pc][0] <<" "<<val[pc][1]<<" "<< val[pc][2]<<" "<< val[pc][3]<<endl;
	//double wynik = (val[pc][0] * coord[0]) + (val[pc][1] * coord[1]) + (val[pc][2] * coord[2]) + (val[pc][3] * coord[3])
	double wynik = 0;
	for (int i = 0; i < 4; i++) {
		wynik += val[pc][i] * coord[i];
	}
	return wynik;
}
double** Oblicz_macierzHdlaPc(double** tabdNdx, double** tabdNdy, double con, double* detJ, int pc) {
	double** macierzH = new double* [4];
	for (int i = 0; i < 4; i++) {
		macierzH[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			macierzH[i][j] = con * ((tabdNdx[pc][i] * tabdNdx[pc][j]) + (tabdNdy[pc][i] * tabdNdy[pc][j])) * detJ[pc];
		}
	}
	return macierzH;
}
double** Oblicz_macierzHbc(double* ksi, double* eta, int alfa, double detJ, int pc) {
	double** macierzHbc = new double* [4];
	//double* N = new double[4];
	vector<double>N(4);
	for (int i = 0; i < 4; i++) {
		macierzHbc[i] = new double[4];
	}
	//wypelnienie macierzy 0 aby dla elementow bez warunkow brzegowych funkcja nie zwracala smieci
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			macierzHbc[i][j] = 0;
		}
	}
	for (int o = 0; o < pc; o++) {
		N = FunKsztaltu(ksi[o], eta[o]);
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzHbc[i][j] += ((N[i] * N[j]) * alfa) * detJ;
			}
		}
	}
	return macierzHbc;
}
double* Oblicz_vecP(double* ksi, double* eta, int alfa, double detJ, int pc, int tot) {
	double* vecP=new double[4];
	vector<double>N(4);
	for (int i = 0; i < 4; i++) {
		vecP[i] = 0;
	}
	for (int i = 0; i < pc; i++) {
		N = FunKsztaltu(ksi[i], eta[i]);
		for (int j = 0; j < 4; j++) {
			vecP[j] += alfa * N[j] * tot * detJ;
		}
	}
	return vecP;
}
double** Oblicz_macierzC_dlaPc(double** tabN, int spcHeat, int ro, double* detJ, int pc) {
	double** macierzC = new double* [4];
	for (int i = 0; i < 4; i++) {
		macierzC[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			macierzC[i][j] = spcHeat * ro * (tabN[pc][i] * tabN[pc][j]) * detJ[pc];
		}
	}
	return macierzC;
}
class Elem4 {
public:
	class Bok {
	public:
		vector<double> x;
		vector<double> y;
		/*vector<double> ksi;
		vector<double> eta;*/
		double* ksi;
		double* eta;
		int pc=0;
		int alfa=0;
		double** macierzHbc;
		double l=0;
		double detJ=0;
		vector<int>lokalnyWektorWezlowZBc;
		double* vecP; //Vektor Obciążeń
		int tot=0;
		Bok() {
			macierzHbc = new double* [4];
			vecP = new double [4];
			for (int i = 0; i < 4; i++) {
				macierzHbc[i] = new double[4];
				vecP[i] = 0;
			}
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					macierzHbc[i][j] = 0;
				}
			}
		}
		void CalcHBC(double* ksi, double* eta) {
			double** HBClokal = Oblicz_macierzHbc(ksi, eta, alfa, detJ, pc);
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					macierzHbc[i][j] += HBClokal[i][j];
				}
			}
		}
		void ClearHBC() {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					macierzHbc[i][j] = 0;
				}
			}
		}
		void CalcP(double* ksi, double* eta) {
			double* vecP_Lokal = Oblicz_vecP(ksi, eta, alfa, detJ, pc, tot);
			for (int i = 0; i < 4; i++) {
				vecP[i] += vecP_Lokal[i];
			}
			
			/*for (int i = 0; i < 4; i++) {
				vecP[i] += Oblicz_vecP(ksi, eta, alfa, detJ, pc, tot)[i];
			}*/
		}
		void ClearP() {
			for (int i = 0; i < 4; i++) {
				vecP[i] = 0;
			}
		}

		void UpdateArg(vector<double> x, vector<double> y, vector<int> lokalnyWektorWezlowZBc, int alfa, int pc, int tot) {
			this->x = x;
			this->y = y;
			this->lokalnyWektorWezlowZBc = lokalnyWektorWezlowZBc;
			this->alfa = alfa;
			this->pc = pc;
			this->tot = tot;
			ksi = new double[pc];
			eta = new double[pc];
			if (pc == 2) {
				double local = 1 / sqrt(3);

				ksi[0] = local;
				ksi[1] = -local;

				eta[0] = local;
				eta[1] = -local;
				/*ksi.push_back(local);
				ksi.push_back(-local);

				eta.push_back(local);
				eta.push_back(-local);*/
			}
			if (pc == 3) {
				double local = sqrt(3.0 / 5.0);
				ksi[0] = local;
				ksi[1] = 0;
				ksi[2] = -local;

				eta[0] = local;
				eta[1] = 0;
				eta[2] = -local;
				/*ksi.push_back(local);
				ksi.push_back(0);
				ksi.push_back(-local);

				eta.push_back(local);
				eta.push_back(0);
				eta.push_back(-local);*/
			}
			if (pc == 4) {
				double local[] = { (0.339981),(0.861136) };
				ksi[0] = local[0];
				ksi[1] = local[1];
				ksi[2] = -local[0];
				ksi[3] = -local[1];

				eta[0] = local[0];
				eta[1] = local[1];
				eta[2] = -local[0];
				eta[3] = -local[1];
				/*ksi.push_back(local[0]);
				ksi.push_back(local[1]);
				ksi.push_back(-local[0]);
				ksi.push_back(-local[1]);

				eta.push_back(local[0]);
				eta.push_back(local[1]);
				eta.push_back(-local[0]);
				eta.push_back(-local[1]);*/
			}
		}
		void CalcData() {
			//gora
			if (count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 0) && count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 1)) {
				l = sqrt(pow(abs(x[0] - x[1]), 2) + pow(abs(y[0] - y[1]), 2));
				detJ = l / 2;
				double* etaLoc = new double[pc];
				for (size_t i = 0; i < pc; i++)
				{
					etaLoc[i] = 0;
				}
				//eta.clear();
				for (int i = 0; i < pc; i++) {
					etaLoc[i] = 1;
					//eta.push_back(1);
				}
				CalcHBC(ksi, etaLoc);
				CalcP(ksi, etaLoc);
			}
			//lewa
			if (count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 1) && count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 2)) {
				l = sqrt(pow(abs(x[1] - x[2]), 2) + pow(abs(y[1] - y[2]), 2));
				detJ = l / 2;
				double* ksiLoc = new double[pc];
				for (size_t i = 0; i < pc; i++)
				{
					ksiLoc[i] = 0;
				}
				//ksi.clear();
				for (int i = 0; i < pc; i++) {
					ksiLoc[i] = -1;
					//ksi.push_back(-1);
				}
				CalcHBC(ksiLoc, eta);
				CalcP(ksiLoc, eta);
			}
			//dol
			if (count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 2) && count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 3)) {
				//cout <<endl<< "Jestem na dole" << endl;
				l = sqrt(pow(abs(x[2] - x[3]), 2) + pow(abs(y[2] - y[3]), 2));
				detJ = l / 2;
				double* etaLoc = new double[pc];
				for (size_t i = 0; i < pc; i++)
				{
					etaLoc[i] = 0;
				}
				//eta.clear();
				for (int i = 0; i < pc; i++) {
					etaLoc[i] = -1;
					//eta.push_back(-1);
				}
				CalcHBC(ksi, etaLoc);
				CalcP(ksi, etaLoc);
			}
			//prawo
			if (count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 0) && count(lokalnyWektorWezlowZBc.begin(), lokalnyWektorWezlowZBc.end(), 3)) {
				l = sqrt(pow(abs(x[0] - x[3]), 2) + pow(abs(y[0] - y[3]), 2));
				detJ = l / 2;
				double* ksiLoc = new double[pc];
				for (size_t i = 0; i < pc; i++)
				{
					ksiLoc[i] = 0;
				}
				//ksi.clear();
				for (int i = 0; i < pc; i++) {
					ksiLoc[i] = 1;
					//ksi.push_back(1);
				}
				CalcHBC(ksiLoc, eta);
				CalcP(ksiLoc, eta);
			}
		}
		~Bok() {
			for (int i = 0; i < 4; i++) {
				delete[] macierzHbc[i];
			}
			delete[] macierzHbc;
			delete[] vecP;
			delete[] ksi;
			delete[] eta;
		}

	};
	int n = 0;
	//[]wiersz []kolumna
	double** tab_dN_przez_dKsi;	//tablica dla pochodnych dN/dKsi
	double** tab_dN_przez_dEta;	//tablica dla pochodnych dN/dEta
	double** tabN;				//tablica z funkcjami kształtu
	double** tabdNdx;
	double** tabdNdy;

	double** macierzHPC;
	double** macierzHBC;
	double** finalMatrixH;
	double** macierzH_plus_HBC;
	double** macierzC_PC;
	double** finalMacierzC;

	double** macierzHpom;
	double* vecPpom;
	vector<double> x;
	vector<double> y;

	int k;
	double* detJ;
	double* eta;
	double* ksi;
	int nN;
	int* idElem = new int[4];
	double* vecP;
	vector <int> tabWezlowBc;
	int iloscBc_Dla_Elementu;
	int alfa;
	int tot; //temperatura otoczenia
	int spcHeat;
	int ro;
	//float tempInit;
	int symTime, symStepTime; 
	Bok bok;

	float minTemp, maxTemp;
	double* temp1;
	double* temp0;
	Elem4(vector<int> tabWezlowBc, int alfa, int tot, int spcHeat, int ro, float tempInit, int symTime, int symStepTime) {
		this->tabWezlowBc = tabWezlowBc;
		this->alfa = alfa;
		this->tot = tot;
		this->spcHeat = spcHeat;
		this->ro = ro;
		this->temp0 = temp0;
		this->symTime = symTime;
		this->symStepTime = symStepTime;
	}

	void UpdateArg(int n, int k, vector<double> x, vector<double> y, int nN, int* idElem) {
		this->n = n;
		this->x = x;
		this->y = y;
		this->k = k;
		this->nN = nN;
		this->idElem = idElem;
		detJ = new double[pow(n, 2)];
		eta = new double[pow(n, 2)];
		ksi = new double[pow(n, 2)];
		vecP = new double[4];
		vector <int> tabWezlowBcLokal;//Vektor przechowywujący lokalne id węzłów z BC dla elementu
		iloscBc_Dla_Elementu = 0;

		//wypelnienie tablic ksi i eta wartościami odpowiednio dla wybranej ilości punktów całkowania
		if (n == 2) {
			eta[0] = (-(1 / (sqrt(3))));	//PC1
			eta[1] = (-(1 / (sqrt(3))));	//PC2
			eta[2] = (1 / (sqrt(3)));		//PC3
			eta[3] = (1 / (sqrt(3)));		//PC4

			ksi[0] = (-(1 / (sqrt(3))));
			ksi[1] = (1 / (sqrt(3)));
			ksi[2] = (-(1 / (sqrt(3))));
			ksi[3] = (1 / (sqrt(3)));
		}
		else if (n == 3) {
			eta[0] = (-sqrt((3.0 / 5.0)));
			eta[1] = (-sqrt((3.0 / 5.0)));
			eta[2] = (-sqrt((3.0 / 5.0)));
			eta[3] = 0;
			eta[4] = 0;
			eta[5] = 0;
			eta[6] = (sqrt((3.0 / 5.0)));
			eta[7] = (sqrt((3.0 / 5.0)));
			eta[8] = (sqrt((3.0 / 5.0)));

			ksi[0] = (-sqrt((3.0 / 5.0)));
			ksi[1] = 0;
			ksi[2] = (sqrt((3.0 / 5.0)));
			ksi[3] = (-sqrt((3.0 / 5.0)));
			ksi[4] = 0;
			ksi[5] = (sqrt((3.0 / 5.0)));
			ksi[6] = (-sqrt((3.0 / 5.0)));
			ksi[7] = 0;
			ksi[8] = (sqrt((3.0 / 5.0)));
		}
		else if (n == 4) {
			eta[0] = (-0.861136);
			eta[1] = (-0.861136);
			eta[2] = (-0.861136);
			eta[3] = (-0.861136);
			eta[4] = (-0.339981);
			eta[5] = (-0.339981);
			eta[6] = (-0.339981);
			eta[7] = (-0.339981);
			eta[8] = (0.339981);
			eta[9] = (0.339981);
			eta[10] = (0.339981);
			eta[11] = (0.339981);
			eta[12] = (0.861136);
			eta[13] = (0.861136);
			eta[14] = (0.861136);
			eta[15] = (0.861136);

			ksi[0] = (-0.861136);
			ksi[1] = (-0.339981);
			ksi[2] = (0.339981);
			ksi[3] = (0.861136);
			ksi[4] = (-0.861136);
			ksi[5] = (-0.339981);
			ksi[6] = (0.339981);
			ksi[7] = (0.861136);
			ksi[8] = (-0.861136);
			ksi[9] = (-0.339981);
			ksi[10] = (0.339981);
			ksi[11] = (0.861136);
			ksi[12] = (-0.861136);
			ksi[13] = (-0.339981);
			ksi[14] = (0.339981);
			ksi[15] = (0.861136);
		}

		//Tworzenie miejsca w pamieci dla tablic
		tab_dN_przez_dKsi = new double* [pow(n, 2)];
		tab_dN_przez_dEta = new double* [pow(n, 2)];
		tabdNdx = new double* [pow(n, 2)];
		tabdNdy = new double* [pow(n, 2)];
		tabN = new double* [pow(n, 2)];
		macierzHPC = new double* [4];
		macierzHBC = new double* [4];
		finalMatrixH = new double* [4];
		macierzH_plus_HBC = new double* [4];
		macierzC_PC = new double* [4];
		finalMacierzC = new double* [4];

		macierzHpom=new double*[nN];
		vecPpom = new double[nN];
		temp1 = new double[nN];
		temp0 = new double[nN];
		for (int i = 0; i < nN; i++) {
			vecPpom[i] = 0;
			macierzHpom[i] = new double[nN];

		}
		for (int i = 0; i < pow(n, 2); i++) {
			tab_dN_przez_dKsi[i] = new double[4];
			tab_dN_przez_dEta[i] = new double[4];
			tabdNdx[i] = new double[4];
			tabdNdy[i] = new double[4];
			tabN[i] = new double[4];
		}
		for (int i = 0; i < 4; i++) {
			macierzHPC[i] = new double[4];
			finalMatrixH[i] = new double[4];
			macierzHBC[i] = new double[4];
			macierzH_plus_HBC[i] = new double[4];
			macierzC_PC[i] = new double[4];
			finalMacierzC[i] = new double[4];
			for (int j = 0; j < tabWezlowBc.size(); j++) {
				if (idElem[i] == tabWezlowBc[j]) {
					tabWezlowBcLokal.push_back(i);
					iloscBc_Dla_Elementu++;
				}
			}
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				macierzHpom[i][j] = 0;
			}
		}
		for (int i = 0; i < nN; i++) {
			temp1[i] = 0;
			temp0[i] = 0;
		}
		bok.UpdateArg(x, y, tabWezlowBcLokal, alfa, n, tot);
		bok.CalcData();
	}
	void WypelnijTab(double** matrixAgg, double* vecPAgg, double** matrixCAgg) {
		double dX_przez_dKsi, dY_przez_dKsi, dX_przez_dEta, dY_przez_dEta;	//zmienne potrzebne do obliczenia detJ
		for (int i = 0; i < pow(n, 2); i++) {
			for (int j = 0; j < 4; j++) {
				tab_dN_przez_dKsi[i][j] = PochodnaKsi(eta[i], j + 1);
				tab_dN_przez_dEta[i][j] = PochodnaEta(ksi[i], j + 1);
			}
		}
		for (int i = 0; i < pow(n, 2); i++) {
			dX_przez_dKsi = PochodnaDoJakobianu(tab_dN_przez_dKsi, x, i);
			dY_przez_dKsi = PochodnaDoJakobianu(tab_dN_przez_dKsi, y, i);
			dX_przez_dEta = PochodnaDoJakobianu(tab_dN_przez_dEta, x, i);
			dY_przez_dEta = PochodnaDoJakobianu(tab_dN_przez_dEta, y, i);
			detJ[i] = (dX_przez_dKsi * dY_przez_dEta) - (dX_przez_dEta * dY_przez_dKsi);
			//cout << endl << "detJ: " << detJ[i];
			for (int j = 0; j < 4; j++) {
				tabdNdx[i][j] = ((1 / detJ[i]) * dY_przez_dEta) * tab_dN_przez_dKsi[i][j] + ((1 / detJ[i]) * (-dY_przez_dKsi)) * tab_dN_przez_dEta[i][j];
				tabdNdy[i][j] = ((1 / detJ[i]) * (-dX_przez_dEta)) * tab_dN_przez_dKsi[i][j] + ((1 / detJ[i]) * dX_przez_dKsi) * tab_dN_przez_dEta[i][j];
			}
		}
		for (int i = 0; i < pow(n, 2); i++) {
			vector<double> lokN = FunKsztaltuPoOBJ(ksi[i], eta[i]);
			for (int j = 0; j < 4; j++) {
				tabN[i][j] = lokN[j];
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				finalMatrixH[i][j] = 0;
				finalMacierzC[i][j] = 0;
			}
		}
		for (int i = 0; i < pow(n, 2); i++) {
			macierzHPC = Oblicz_macierzHdlaPc(tabdNdx, tabdNdy, k, detJ, i);
			//cout << endl << "Macierz C dla punktu calkowania nr " << i + 1 << endl;
			macierzC_PC = Oblicz_macierzC_dlaPc(tabN, spcHeat, ro, detJ, i);
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					finalMatrixH[j][k] += macierzHPC[j][k];
					finalMacierzC[j][k] += macierzC_PC[j][k];
					//cout << macierzHPC[j][k] << " ";
					//cout << macierzC_PC[j][k] << " ";
				}
				//cout << endl;
			}
		}
		/*cout << endl << "Wspolrzedne x elementu: ";
		for (int i = 0; i < x.size(); i++) {
			cout << x[i] << " ";
		}
		cout << endl << "Wspolrzedne y elementu: ";
		for (int i = 0; i < y.size(); i++) {
			cout << y[i] << " ";
		}*/

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzHBC[i][j] = 0;
				macierzH_plus_HBC[i][j] = 0;
			}
		}
		macierzHBC = bok.macierzHbc;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				macierzH_plus_HBC[i][j] += macierzHBC[i][j] + finalMatrixH[i][j];
			}
		}
		for (int i = 0; i < 4; i++) {
			vecP[i] = 0;
		}
		vecP = bok.vecP;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				//matrixAgg[idElem[i] - 1][idElem[j] - 1] += finalMatrixH[i][j];
				matrixAgg[idElem[i] - 1][idElem[j] - 1] += macierzH_plus_HBC[i][j];
				matrixCAgg[idElem[i] - 1][idElem[j] - 1] += finalMacierzC[i][j];
			}
			//matrixAgg[idElem[i]]
		}
		for (int i = 0; i < 4; i++) {
			vecPAgg[idElem[i] - 1] += vecP[i];
		}
	}
	void CalcPomocniczeTab(double** matrixAgg, double* vecPAgg, double** matrixCAgg, double* tempinit) {
		double* pomVec = new double[nN];
		for (int i = 0; i < nN; i++) {
			pomVec[i] = 0;
			//cout << "Temp init: "<<tempinit[i] << endl;
		}
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				pomVec[i] += (matrixCAgg[i][j] / symStepTime)*tempinit[j];
			}
		}
		//cout << "Sym step time: " << symStepTime << endl;
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				macierzHpom[i][j] = matrixAgg[i][j] + (matrixCAgg[i][j] / symStepTime);
				vecPpom[i] = vecPAgg[i]+pomVec[i];
			}
		}
		
	
		/*cout <<endl<< "Matrix H+ C/dt" << endl;
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				cout << macierzHpom[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "wektor P+(c/dt)*t0 " << endl;
		for (int i = 0; i < nN; i++) {
			cout << vecPpom[i] << " ";
		}*/
		//temp1 = temp0;
		//maxTemp=max_element(temp1[0], temp1[nN]);
		//cout << "MAX TEMP " << maxTemp << endl;
	}
	void CalcTemp() {
		temp1 = ElimGauss(macierzHpom, vecPpom, nN);
		/*for (int i = 0; i < nN; i++) {
			temp1[i] = ElimGauss(macierzHpom, vecPpom, nN)[i];
		}*/
		/*for (int i = 0; i < nN; i++) {
			cout << "Temp: " << temp1[i] << endl;;
		}*/
		maxTemp = temp1[0];
		minTemp = temp1[0];
		for (int i = 0; i < nN; i++) {
			//cout << "Temp 1: " << temp1[i] << endl;
			if (temp1[i] > maxTemp)
			{
				maxTemp = temp1[i];
			}
			if (temp1[i] < minTemp) {
				minTemp = temp1[i];
			}
		}
		/*cout << "Max temp: " << maxTemp;
		cout << "Min temp: " << minTemp;*/
	}
	void WypiszTab() {

		/*cout << endl << "Eta | dN1/dksi | dN2/dksi | dN3/dksi | dN4/dksi" << endl;
		for (int i = 0; i < pow(n, 2); i++) {
			cout << eta[i] << " ";
			for (int j = 0; j < 4; j++) {
				cout << tab_dN_przez_dKsi[i][j] << "  ";
			}
			cout << endl;
		}

		cout << endl << "Ksi | dN1/deta | dN2/deta | dN3/deta | dN4/deta" << endl;
		for (int i = 0; i < pow(n, 2); i++) {
			cout << ksi[i] << " ";
			for (int j = 0; j < 4; j++) {
				cout << tab_dN_przez_dEta[i][j] << "  ";
			}
			cout << endl;
		}*/
		/*cout << endl << "pc | dN1/dx | dN2/dx | dN3/dx | dN4/dx" << endl;
		for (int i = 0; i < pow(n, 2); i++) {
			cout << i + 1 << " ";
			for (int j = 0; j < 4; j++) {
				cout << tabdNdx[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "pc | dN1/dy | dN2/dy | dN3/dy | dN4/dy" << endl;
		for (int i = 0; i < pow(n, 2); i++) {
			cout << i + 1 << " ";
			for (int j = 0; j < 4; j++) {
				cout << tabdNdy[i][j] << " ";
			}
			cout << endl;
		}*/

		//Wypisanie macierzy H dla pojedynczego punktu całkowania


		cout << endl << "Koncowa macierz H: " << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << finalMatrixH[i][j] << " ";
			}
			cout << endl;
		}
		/*cout << endl << "Macierz HBC: " << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << macierzHBC[i][j] << " ";
			}
			cout << endl;
		}*/
		cout << endl << "Macierz HBC + H: " << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << macierzH_plus_HBC[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Vektor P: " << endl;
		for (int i = 0; i < 4; i++) {
			cout << vecP[i] << " ";
		}
		/*cout << endl << "Fun ksztaltu po objetosci";
		cout << endl << "pc | N1 | N2 | N3 | N4" << endl;
		for (int i = 0; i < pow(n, 2); i++) {
			cout << i + 1 << " ";
			for (int j = 0; j < 4; j++) {
				cout << tabN[i][j] << " ";
			}
			cout << endl;
		}*/
		cout << endl << endl << "Macierz C" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << finalMacierzC[i][j] << " ";
			}
			cout << endl;
		}

		bok.ClearHBC();
		bok.ClearP();
	}
	void WypiszAgg(double** matrixAgg, double* vecPAgg, double** matrixCAgg) {
		cout << endl << endl << "Agregacja: " << endl;
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				cout << matrixAgg[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << endl << "Agregacja macierzy C: " << endl;
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				cout << matrixCAgg[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Wektor P:" << endl;
		for (int i = 0; i < nN; i++) {
			cout << vecPAgg[i] << " ";
		}
		cout << endl;
	}
	~Elem4() {
		for (int i = 0; i < pow(n, 2); i++) {
			delete[] tab_dN_przez_dKsi[i];
			delete[] tab_dN_przez_dEta[i];
			delete[] tabdNdx[i];
			delete[] tabdNdy[i];
			delete[] tabN[i];
		}
		for (int i = 0; i < 4; i++) {
			delete[] macierzHPC[i];
			delete[] finalMatrixH[i];
			//delete[] macierzHBC[i];
			delete[] macierzH_plus_HBC[i];
			delete[] macierzC_PC[i];
			delete[] finalMacierzC[i];
		}
		for (int i = 0; i < nN; i++) {
			delete[] macierzHpom[i];
		}
		delete[] macierzHpom;
		delete[] tab_dN_przez_dKsi;
		delete[] tab_dN_przez_dEta;
		delete[] tabdNdx;
		delete[] tabdNdy;
		delete[] tabN;
		delete[] macierzHPC;
		delete[] finalMatrixH;
		delete[] eta;
		delete[] ksi;
		delete[] macierzH_plus_HBC;
		delete[] macierzC_PC;
		delete[] finalMacierzC;
		delete[] vecPpom;
		delete[] detJ;
		delete[] temp0;
	}
};

struct Node
{
	float x, y; //współrzędne
	float temp; //temperatura
	bool BC; //warunek czy zachodzi warunek brzegowy na węźle
};
struct Element
{
	int id[4];
};
struct Grid
{
	int nN; //liczba węzłów
	int nE; //liczba elementów
};
struct GlobalData
{
	int t, dt;	//czas symulacji, przyrost czasu 
	int spcHeat;//ciepło właściwe	
	int ro;		//gęstość
	int k;		//współczynnik przwodności cieplnej
	int alfa;	//alfa do warunku brzegowego konwekcji	
	int tot;	//temperatura otoczenia
};
string LastWordFromLine(string line)
{
	string word;
	istringstream ss(line);
	while (ss >> word);
	return word;
}
bool LoadFile(string fileName)
{
	ifstream file;
	file.open(fileName.c_str());
	if (!file.good())
		return false;
	GlobalData data;
	Node node;
	Element element;
	Grid grid;
	string line;
	if (!file.fail())
	{
		getline(file, line);
		data.t = stoi(LastWordFromLine(line));			//pobranie ostatniego słowa oraz zamienienie go na liczbę stałą
		cout << "czas symulacji = " << data.t << endl;

		getline(file, line);
		data.dt = stoi(LastWordFromLine(line));
		cout << "przyrost czasu = " << data.dt << endl;

		getline(file, line);
		data.k = stoi(LastWordFromLine(line));
		cout << "wspolczynnik przewodnosci cieplnej = " << data.k << endl;

		getline(file, line);
		data.alfa = stoi(LastWordFromLine(line));
		cout << "alfa do warunku brzegowego konwekcji = " << data.alfa << endl;

		getline(file, line);
		data.tot = stoi(LastWordFromLine(line));
		cout << "temperatura otoczenia = " << data.tot << endl;

		getline(file, line);
		node.temp = stof(LastWordFromLine(line));
		cout << "temperatura poczatkowa = " << node.temp << endl;

		getline(file, line);
		data.ro = stoi(LastWordFromLine(line));
		cout << "gestosc = " << data.ro << endl;

		getline(file, line);
		data.spcHeat = stoi(LastWordFromLine(line));
		cout << "cieplo wlasciwe = " << data.spcHeat << endl;

		getline(file, line);
		grid.nN = stoi(LastWordFromLine(line));
		Node* nodes = new Node[grid.nN];
		cout << "liczba wezlow = " << grid.nN << endl;

		getline(file, line);
		grid.nE = stoi(LastWordFromLine(line));
		Element* elements = new Element[grid.nE];
		cout << "liczba elementow = " << grid.nE << endl;

		cout << endl << "Wspolrzedne wezlow: " << endl;
		for (int i = 0; i < grid.nN; i++)
		{
			getline(file, line, ',');	//pobranie całej linijki do momentu przecinka
			getline(file, line, ',');
			float x = stof(line);
			getline(file, line);
			float y = stof(LastWordFromLine(line));
			nodes[i].x = x;
			nodes[i].y = y;
			cout << "Wezel " << i + 1 << ": x = " << nodes[i].x << ", y = " << nodes[i].y << endl;
		}
		getline(file, line);
		cout << endl << "Id elementow: " << endl;
		for (int i = 0; i < grid.nE; i++)
		{
			getline(file, line, ',');
			getline(file, line, ',');
			int a = stoi(line);
			getline(file, line, ',');
			int b = stoi(line);
			getline(file, line, ',');
			int c = stoi(line);
			getline(file, line, '\n');
			int d = stoi(line);
			elements[i].id[0] = a;
			elements[i].id[1] = b;
			elements[i].id[2] = c;
			elements[i].id[3] = d;

			cout << "e" << i + 1 << " id = [" << elements[i].id[0] << " " << elements[i].id[1] << " " << elements[i].id[2] << " " << elements[i].id[3] << "]" << endl;
		}
		getline(file, line);

		while (true)
		{
			if (!file.fail())
			{
				getline(file, line, ',');
				int id = stoi(line) - 1;
				nodes[id].BC = true;
			}
			else
				break;
		}
		int rozmiar = 0;
		cout << "Warunek brzegowy zachodzi w wezlach nr: " << endl;
		for (int i = 0; i < grid.nN; i++)
		{
			if (nodes[i].BC == true)
			{
				rozmiar++;
				cout << i + 1 << ", ";
			}
		}
		vector<int> tabWezlowBc(rozmiar, 0);
		int j = 0;
		for (int i = 0; i < grid.nN; i++) {
			if (nodes[i].BC == true)
			{
				tabWezlowBc[j] = i + 1;
				j++;
			}
		}
		/*cout << endl << "Wspolrzedne pierwszego elementu: " << endl;
		cout << "x: " << nodes[elements[0].id[0] - 1].x << ", " << nodes[elements[0].id[1] - 1].x << ", " << nodes[elements[0].id[2] - 1].x << ", " << nodes[elements[0].id[3] - 1].x << endl;
		cout << "y: " << nodes[elements[0].id[0] - 1].y << ", " << nodes[elements[0].id[1] - 1].y << ", " << nodes[elements[0].id[2] - 1].y << ", " << nodes[elements[0].id[3] - 1].y << endl;*/

		double** matrixHAgg = new double* [grid.nN];
		double** matrixCAgg = new double* [grid.nN];
		double* vecPAgg = new double[grid.nN];
		for (int i = 0; i < grid.nN; i++) {
			matrixHAgg[i] = new double[grid.nN];
			matrixCAgg[i] = new double[grid.nN];
		}
		for (int i = 0; i < grid.nN; i++) {
			vecPAgg[i] = 0;
			for (int j = 0; j < grid.nN; j++) {
				matrixHAgg[i][j] = 0;
				matrixCAgg[i][j] = 0;
			}
		}
		
		int n;
		//int l;
		//for (;;) {
			cout << "wybierz schemat calkowania" << endl;
			cout << "2. - 2 punktowy, 3. - 3 punktowy, 4. - 4 punktowy, inna liczba = koniec" << endl;
			cin >> n;
			if (n == 2 || n == 3 || n == 4) {
				Elem4 elem(tabWezlowBc, data.alfa, data.tot, data.spcHeat, data.ro, node.temp, data.t, data.dt);
				for (int i = 0; i < grid.nE; i++) {
					vector<double> elemX = { nodes[elements[i].id[0] - 1].x,nodes[elements[i].id[1] - 1].x,nodes[elements[i].id[2] - 1].x,nodes[elements[i].id[3] - 1].x };		//wspolrzedne x punktu siatki
					vector<double> elemY = { nodes[elements[i].id[0] - 1].y,nodes[elements[i].id[1] - 1].y,nodes[elements[i].id[2] - 1].y,nodes[elements[i].id[3] - 1].y };		//wspolrzedne y punkty siatki
					//Elem4 elem(n, data.k, elemX, elemY, grid.nN, elements[i].id);
					cout << endl << "==========Element " << i + 1 << " =========";
					elem.UpdateArg(n, data.k, elemX, elemY, grid.nN, elements[i].id);
					elem.WypelnijTab(matrixHAgg, vecPAgg, matrixCAgg);
					elem.WypiszTab();
					//system("pause");
				}
				elem.WypiszAgg(matrixHAgg, vecPAgg, matrixCAgg);
				double* initialTemp = new double[grid.nN];
				for (int i = 0; i < grid.nN; i++) {
					initialTemp[i] = node.temp;
				}
				elem.CalcPomocniczeTab(matrixHAgg, vecPAgg, matrixCAgg, initialTemp);
				elem.CalcTemp();
				cout <<endl<< "Czas | MinTemp | MaxTemp " << endl;
				//cout << 50 << " " << elem.minTemp << " " << elem.maxTemp << endl;
				for (int i = data.dt; i <= data.t; i = i + data.dt) {
					elem.CalcTemp();
					initialTemp = elem.temp1;
					elem.CalcPomocniczeTab(matrixHAgg, vecPAgg, matrixCAgg, initialTemp);
					cout << i <<" "<<elem.minTemp<<" " <<elem.maxTemp<< endl;
				}
				delete[] initialTemp;
			}
			else {
				return 0;
			}
			
		
		//}
		for (int i = 0; i < grid.nN; i++) {
			delete[] matrixHAgg[i];
			delete[] matrixCAgg[i];
		}
		delete[]nodes;
		delete[]elements;
		delete[] matrixHAgg;
		delete[]matrixCAgg;
		delete[]vecPAgg;
	}
	else
		return -1;
	return true;
}


int main() {
	//if (!LoadFile("Test1_4_4.txt"))
	if (!LoadFile("Test2_4_4_MixGrid.txt"))
	//if (!LoadFile("Test3_31_31_kwadrat.txt"))
		cout << "Fail to load" << endl;

	return 1;
}