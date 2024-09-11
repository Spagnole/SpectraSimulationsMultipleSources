#ifndef MY_CLASS_H // include guard
#define MY_CLASS_H


const int NSources = 10; //this number can be changed if needed
int used_sources = 0;
vector<string> source_name;

bool PrintPeakWidths = false;
bool PrintEscPeakData = false;

class MyTransition{
		public:
				int Index;
				double lvlEn; //0
				double gammaEn; //1
				double finalLvl; //2
				double gammaInt; //3
				double gammaIntError; //4
				double lvlPop; //5
				double gammaBr; //6
				double ICC; //internal conversion coefficient
				//string Multipolarity; //E1, E2, M1, M2, etc..
		
		void Display(){
				cout << Index << "\t" <<lvlEn << "\t" << gammaEn << "\t" << finalLvl << "\t" <<
				gammaInt << "\t" <<	gammaIntError << "\t" << lvlPop << "\t" << gammaBr << endl;
		}
};
vector<MyTransition> gammaList[NSources];

class MyLevels{
public:
		int Index;
		double lvlEn;
		double lvlPop;
		void Display(){
				cout << Index << "\t" << lvlEn << "\t" << lvlPop << endl;
		}
};
vector<MyLevels> levelList[NSources];

class gammagamma{
public:
		double gammaOne;
		double gammaTwo;
		double coincidences;
		void Display(){
				cout << gammaOne << "\t" << gammaTwo << "\t" << coincidences << endl;
		}
};
vector<gammagamma> coincList[NSources];

class Cube {		
public:
		int IndexOne;
		int IndexTwo;
		int IndexThree;
		double gammaOne;
		double gammaTwo;
		double gammaThree;
		double triples;
		vector<MyTransition> MyPath;
		void Display(){
				cout << IndexOne << " " << gammaOne << "\t" << IndexTwo << " " << gammaTwo << "\t"
				<< IndexThree << " " << gammaThree << "\t" << triples << " initial\t";// << endl;
				if( MyPath.size() != 0) MyPath.at(0).Display();
				else cout << endl;
		}
};
vector<Cube> cubeList[NSources];



TFile *fRealData; //file contaaining real data

TF1 *fWidth; //linear function to fit peak widths as a function of energy
TF1 *fEscPeak;	//quadratic function to fit escape peak intensities as a function of energy


TH1D *hRealSpectra;	//real experimental histogram
TH1D *hSimPeaks[NSources]; //simulated peak intensities
TH1D *hEscPeaks[NSources]; //simulated escape peak intensities
TH1D *hSimSource[NSources]; //simulated full energy and single escape peaks for given source
TH1D *hBkgr; //simulated background extracted from real spectra
TH1D *hFullSim; //full simulated spectra including simulated full-energy and single-escape peaks on top of background
TH1D *hSimPeaks_Bkgr;  //simulated spectra using only simulated full-energy peaks on top of background
TH1D *hEscPeaks_Bkgr; //simulated spectra using only simulated single-escape peaks on top of background

TGraph *gEff; // = new TGraph(); //graph for efficiency curve
TGraphErrors *gSigma; // = new TGraphErrors(); //graph for peak widths, !!!Please use sigma not FWHM
TGraphErrors *gEscPeaks; // = new TGraphErrors(); //graph to determine escape peak intensity



#endif
