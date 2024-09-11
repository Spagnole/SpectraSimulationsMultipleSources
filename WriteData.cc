#include "ReadInputData.hh"

#ifndef WriteFunctions
#define WriteFunctions

void printGammaList(int SourceNum = 0){
		for(auto MyGamma : gammaList[SourceNum] ) MyGamma.Display();
}

void printLevels(int source_number = 0){
		for(auto level : levelList[source_number]) level.Display();
}

void printLevelsReverse(int source_number = 0){
		for(auto level = levelList[source_number].rbegin(); level != levelList[source_number].rend(); level++)	level->Display();
}

void printCoincidences(int source_number){
		for(auto gg : coincList[source_number]) gg.Display();
}

void printCubeList(int SourceNum = 0){
		for(auto ggg : cubeList[SourceNum] ) ggg.Display();
}

//write the gammaList to a file
void WriteGammaListToFile(const std::string& filename, int SourceNum = 0) {
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
				std::cerr << "Error opening file: " << filename << std::endl;
				return;
		}

		// Write header (optional)
		outfile << "Index\tlvlEn\tgammaEn\tfinalLvl\tgammaInt\tgammaIntError\tlvlPop\tgammaBr\n";
		
		// Iterate through gammaList[0] and write each MyTransition to the file
		for (const auto& transition : gammaList[SourceNum]) {
				outfile << transition.Index << "\t"
								<< transition.lvlEn << "\t"
								<< transition.gammaEn << "\t"
								<< transition.finalLvl << "\t"
								<< transition.gammaInt << "\t"
								<< transition.gammaIntError << "\t"
								<< transition.lvlPop << "\t"
								<< transition.gammaBr << "\n";
		}

		outfile.close();
		std::cout << "Gamma list successfully written to " << filename << std::endl;
}

//write the gammaList to a file
void WriteLevelListToFile(const std::string& filename, int SourceNum = 0) {
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
				std::cerr << "Error opening file: " << filename << std::endl;
				return;
		}

		// Write header (optional)
		outfile << "Index\tlvlEn\tlvlPop\n";
		
		// Iterate through gammaList[0] and write each MyTransition to the file
		for (const auto& Level : levelList[SourceNum]) {
				outfile << Level.Index << "\t"
								<< Level.lvlEn << "\t"
								<< Level.lvlPop << "\n";
		}

		outfile.close();
		std::cout << "Level list successfully written to " << filename << std::endl;
}



#endif
