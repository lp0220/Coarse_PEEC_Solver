#include "Mod_Read_File.h"
#include "Mod_Solver_Freq.h"
#include "Mod_Solver_Time.h"
#include <fstream>
#include <iomanip>
#include <string>

static void SaveMatrixToTxt(const std::vector<std::vector<double>>& matrix, const std::string& file_path) {
	std::ofstream fout(file_path);
	if (!fout.is_open()) {
		std::cerr << "Fail to write matrix file: " << file_path << std::endl;
		return;
	}

	fout << std::scientific << std::setprecision(12);
	for (const auto& row : matrix) {
		for (size_t j = 0; j < row.size(); ++j) {
			fout << row[j];
			if (j + 1 < row.size()) {
				fout << " ";
			}
		}
		fout << '\n';
	}
}

int main(int argc, char* argv[]) {
	if (argc > 1) {
		// Integrated mode: data directory provided via command line
		std::string data_path = argv[1];
		if (!data_path.empty() && data_path.back() != '\\' && data_path.back() != '/') {
			data_path += "\\";
		}
		MAP_PATH = data_path;
		SET_FILE = data_path + "set.txt";
		DIELECTRIC_FILE = data_path + "DIELECTRIC.txt";
		SOURCE_PATH = data_path;

		Read_Setting();
		Read_Dielectric();

		std::string mode = (argc > 2) ? argv[2] : "freq";
		if (mode == "time") {
			Time_Solver();
		}
		else {
			Freq_Solver();
			SaveMatrixToTxt(A_E, MAP_PATH + "A_E.txt");
			SaveMatrixToTxt(A_ET, MAP_PATH + "A_ET.txt");
		}
	}
	else {
		// Original standalone mode
		Read_File();
		Freq_Solver();
		SaveMatrixToTxt(A_E, MAP_PATH + "A_E.txt");
		SaveMatrixToTxt(A_ET, MAP_PATH + "A_ET.txt");
	}
	return 0;
}