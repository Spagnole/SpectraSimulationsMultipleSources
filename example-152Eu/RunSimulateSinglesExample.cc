#include "../SimulateSingles.cc"

void RunFullSimulation(){

	ReadDecayScheme("152Sm.dat","152Sm");
	ReadDecayScheme("152Gd.dat","152Gd");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hgE_152Eu");
	GetSpectrumBackground(hRealSpectra);
	FillSimulation(0);
	FillSimulation(1);
	FillEscapePeaks(0);
	FillEscapePeaks(1);
	BuildSimuledSpectra();
	
}


