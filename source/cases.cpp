#include "cases.hpp"

const std::vector<double> ref_lvl = {0.08, 0.06, 0.05, 0.04, 0.035};

void tryMesh(const std::string& filemesh, const std::string& fileres, const std::vector<int>& b_in, const functionN& fN, const double& beta)
{
	FastMarch f(filemesh, b_in, beta, fN);
	pairXY res;
	std::ofstream fp;

	f.solve(100000);
	res = f.getSurface(300);

	fp.open(fileres, std::ofstream::out);
	
	fp << "Total Volume: " << trapInt(res) << ", Initial surface: " << res[0].second << std::endl;
	fp << "#" << std::endl;	
	for(const auto&[t, A] : res){
		fp << t << " " << A << std::endl;
	}
	fp.close();
}

void tryMesh(const MeshInfo& minfo, const std::string& fileres, const std::vector<int>& b_in, const functionN& fN, const double& beta)
{
	FastMarch f(minfo, b_in, beta, fN);
	pairXY res;
	std::ofstream fp;

	f.solve(100000);
	res = f.getSurface(300);

	fp.open(fileres, std::ofstream::out);
	
	fp << "Total Volume: " << trapInt(res) << ", Initial surface: " << res[0].second << std::endl;
	fp << "#" << std::endl;	
	for(const auto&[t, A] : res){
		fp << t << " " << A << std::endl;
	}
	fp.close();
}

std::vector<int> str2vec(const std::vector<std::string>& s)
{
	std::vector<int> res;

	for(const std::string& ss : s){
		try{
			res.emplace_back(std::stoi(ss));
		}
		catch(const std::invalid_argument& is) {
			std::cout << "Error reading the file." << std::endl;
			std::cout << "Invalid characters in file: " << is.what() << std::endl;
			exit(-1);
		}
	}
	return res;
}

void writeThings(const std::string& filename, const FastMarch* const f, const pairXY& res)
{
	std::fstream fp;
	uint i = 0;

	while(true){
		if(!fs::exists(filename + "_" + std::to_string(i) + ".txt") && !fs::exists(filename + "_info" + std::to_string(i) + ".txt")) break;
		++i;
	}

	fp.open(filename + "_" + std::to_string(i) + ".txt", std::ofstream::out);

	std::cout << "Info\t: Writing to file . . ." << std::endl;

	for(const auto&[wb, As] : res){
		fp << wb << "\t" << As << std::endl;
	}

	fp.close();
	fp.open(filename + "_info" + std::to_string(i) + ".txt", std::ofstream::out);
	fp << "Mesh info:\nNumber of nodes: " << f->getNnode() << std::endl;
	fp << "Number of elements: " << f->getNelem() << std::endl;
	//fp << "Read time: " << std::chrono::duration_cast<std::chrono::seconds>(mid1 - start).count() << " seconds.\n" << std::endl;
	
	fp << "Solution info:\nTotal volume: " << std::fixed << trapInt(res) << ", Initial surface: " << std::fixed << res[0].second << std::endl;
	//fp << "Solving time: " << std::chrono::duration_cast<std::chrono::seconds>(mid2 - mid1).count() << " seconds." << std::endl;
	//fp << "Geting surface time: " << std::chrono::duration_cast<std::chrono::seconds>(end - mid2).count() << " seconds." << std::endl;
	fp.close();
	std::cout << "Info\t: Done." << std::endl;
}

void printHelp()
{
	std::cout << "Possible program usage, parenthesis indicate optional arguments:\n" << std::endl;
	std::cout << ".\\eikonal.exe FILE.NAME (-s XXX) (-f XXX) (-p) (-m) (-r XXX) (-rf XXX) (-e)\n" << std::endl;
	std::cout << "-s: number of sample surfaces" << std::endl;
	std::cout << "-f: flame spreading speed (if applicable)" << std::endl;
	std::cout << "-p: do not plot the graph after even if script is available" << std::endl;
	std::cout << "-m: force meshing even if .diff file is detected" << std::endl;
	std::cout << "-r: mesh refinement level from 0 (extra coarse) to 4 (extra fine)" << std::endl;
	std::cout << "-rf: give the mesh size factor directly instead of the pre-defined levels" << std::endl;
	//std::cout << "-e: export mesh file (only if mesh is made by the program)" << std::endl;
	//std::cout << "\nExample:" << std::endl;
	//std::cout << ".\\eikonal.exe TPS.diff -n 150 -f 200\n" << std::endl; 
}

void readFileSolve(const int& argc, char* argv[])
{
	FastMarch *f;
	MeshInfo minfo;
	std::string filename, aux;
	std::fstream fp;
	std::vector<int> b_in = {}, b_start = {};
	double aux_d = 100.0;
	std::vector<std::string> aux_s, sep;
	functionN fN = fTrue;
	pairXY res;
	int n_surf = 150, aux_i;
	bool force_mesh = false, exp = false, prt = true;
	std::chrono::time_point<std::chrono::high_resolution_clock> start, mid1, mid2, end;

	/*
	Taking care of input arguments
	*/
	if(argc < 2){
		std::cout << "Error\t: No file was given." << std::endl;
		printHelp();
		return;
	}
	for(int i = 2; i < argc; ++i){
		aux = convertCharP(argv[i]);
		// no ploting
		if(aux == "-p") prt = false;
		// force meshing, even if .diff file is present
		else if(aux == "-m") force_mesh = true;
		// export mesh
		//else if(aux == "-e") exp = true;
		// refinement factor
		else if(aux == "-r"){
			++i;
			if(i == argc){
				std::cout << "\nArgument -r requires a numeral to select the refinement level." << std::endl;
				std::cout << "-r XXX" << std::endl;
			}
			aux = convertCharP(argv[i]);
			try{
				aux_i = std::stoi(argv[i]);
			}
			catch(const std::invalid_argument& is){
				std::cout << "\nInvalid numeral given for the refinement level." << std::endl;
				return;
			}

			if(aux_i < 0 || aux_i >= (int)ref_lvl.size()){
				std::cout << "\nRefinement level must be between 0 and " << ref_lvl.size()-1 << "." << std::endl;
				return;
			}
			minfo.mSize = ref_lvl[aux_i];
		}
		else if(aux == "-rf"){
			++i;
			if(i == argc){
				std::cout << "\nArgument -rf requires a numeral to select the mesh size factor." << std::endl;
				std::cout << "-r XXX" << std::endl;
			}
			aux = convertCharP(argv[i]);
			try{
				aux_d = std::stod(argv[i]);
			}
			catch(const std::invalid_argument& is){
				std::cout << "\nInvalid numeral given for the mesh size factor." << std::endl;
				return;
			}

			if(aux_d < 0 || aux_d >= 1){
				std::cout << "\nMesh size factor must be between 0 and 1." << std::endl;
				return;
			}
			minfo.mSize = aux_d;
		}
		// flame spreading speed
		else if(aux == "-f"){
			++i;
			if(i == argc){
				std::cout << "\nArgument -f requires a numeral argument for the flame spreading speed." << std::endl;
				std::cout << "-f XXX" << std::endl;
				return;
			}
			aux = convertCharP(argv[i]);
			try {
				aux_d = std::stod(aux);
			}
			catch(const std::invalid_argument& is){
				std::cout << "\nInvalid numeral given for flame spreading speed." << std::endl;
				return;
			}
			if(aux_d < 0){
				std::cout << "\nFlame spreading speed must be greater or equal to zero." << std::endl;
				return;
			}
			//std::cout << "\nFlame spreading speed is " << aux_d << "." << std::endl; 
		}
		// number of sample surfaces in the file
		else if(aux == "-s"){
			++i;
			if(i == argc){
				std::cout << "\nArgument -n requires a numeral argument for the number of surface samples." << std::endl;
				std::cout << "-n XXX" << std::endl;
				return;
			}
			aux = convertCharP(argv[i]);
			try {
				n_surf = std::stoi(aux);
			}
			catch(const std::invalid_argument& is){
				std::cout << "\nInvalid numeral given for surface sample number." << std::endl;
				return;
			}
			if(n_surf < 2){
				std::cout << "\nNumber of surface samples must be at least two." << std::endl;
				return;
			}
			//std::cout << "\nNumber of surface samples is " << n_surf << "." << std::endl; 
		}
		// unrecognized argument
		else{
			std::cout << "\nUnrecognized argument \"" << aux << "\"" << std::endl;
			printHelp();
			return;
		}
	}





	// get filename without exention and open geometry file
	filename = remExt(convertCharP(argv[1]));
	fp.open(filename + ".geo", std::ifstream::in);
	minfo.f = filename;

	if(!fp){
		std::cout << "Error\t: Geometry file " << filename << ".geo could not be found." << std::endl;
		return;
	}

	if(!(fs::exists(fs::path(filename + ".diff")) || force_mesh)){
		std::cout << "Warning\t: No mesh file detected, will generate own mesh." << std::endl;
		force_mesh = true;
	}

	// get physical surfaces
	while(fp.good()){
		aux = readLine(fp);
		if(aux.find("Physical Surface") == 0){
			aux_s.emplace_back(aux);
		}
	}

	fp.close();

	if(!aux_s.size()){
		std::cout << "Warning\t: No physical groups detected. Entire boundary surface will be considered as initial condition." << std::endl;
	}
	else if(force_mesh){
		for(const std::string& s : aux_s){
			if(s.find("Surface(\"I\"") != std::string::npos){
				minfo.in_tag = std::stoi(s.substr(21));
				b_in = {minfo.in_tag};
			}
			else if(s.find("Surface(\"U\"") != std::string::npos){
				minfo.un_tag = std::stoi(s.substr(21));
				fN = [minfo](const Node& n){
					if(findIn(n.bound, minfo.un_tag) != -1) return true;
					return false;
				};
			}
			if(!(minfo.un_tag == -1 || minfo.in_tag == -1)) break;
		}
	}
	else{
		for(const std::string& s : aux_s){
			if(s.find("Surface(\"I\"") != std::string::npos){
				aux = s.substr(s.find('{')+1, s.find('}')-s.find('{')-1);
				b_in = str2vec(spltString(aux, ','));
			}
			else if(s.find("Surface(\"U\"") != std::string::npos){
				aux = s.substr(s.find('{')+1, s.find('}')-s.find('{')-1);
				b_start = str2vec(spltString(aux, ','));

			}
		}
	}

	start = std::chrono::high_resolution_clock::now();
	if(force_mesh) f = new FastMarch(minfo, b_in, aux_d, fN);
	else f = new FastMarch(filename + ".diff", b_in, aux_d, fN);	
	mid1 = std::chrono::high_resolution_clock::now();
	f->solve(100000);
	mid2 = std::chrono::high_resolution_clock::now();
	res = f->getSurface(n_surf);
	end = std::chrono::high_resolution_clock::now();

	//std::cout << "Info\t: Writing files. . ." << std::endl;
	writeThings(filename, f, res);
	//std::cout << "Info\t: Done.\n" << std::endl;

	std::cout << ((force_mesh) ? "Info\t: Meshing time: " : "Reading mesh time: ") << std::chrono::duration_cast<std::chrono::seconds>(mid1 - start).count() << " seconds." << std::endl;
	std::cout << "Info\t: Solving time: " << std::chrono::duration_cast<std::chrono::seconds>(mid2 - mid1).count() << " seconds." << std::endl;
	std::cout << "Info\t: Getting surface time: " << std::chrono::duration_cast<std::chrono::seconds>(end - mid2).count() << " seconds." << std::endl;
	std::cout << "Info\t: Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds." << std::endl;

	if(exp){
		f->exportMesh();
	}

	delete f;

	if(prt){
		if(fs::exists(fs::path(".\\plotBurnArea.py"))){
			aux = "python .\\plotBurnArea.py \"" + filename + ".txt\"";
			system(aux.c_str());
		}
		else{
			std::cout << "\nError\t: Plotting script could't be found." << std::endl;
		}
	}
}







void performFlameSpreading(const std::string& filename, const std::string& baselog, const std::vector<double>& beta, const std::vector<double>& xs)
{
	MeshInfo minfo;
	std::fstream fp(filename + ".geo", std::fstream::in);
	std::string aux;
	std::vector<int> b_in = {};
	functionN fN;
	pairXY res;

	if(!fp) return;

	// get physical surfaces
	while(fp.good()){
		aux = readLine(fp);
		if(aux.find("Physical Surface(\"I\"") == 0){
			minfo.in_tag = std::stoi(aux.substr(21));
			b_in = {minfo.in_tag};
			break;
		}
	}
	fp.close();

	minfo.f = filename;
	FastMarch f(minfo, b_in, 0, fTrue);

	for(const double& x : xs){
		fN = [x](const Node& n){
			return n.coord[1] > x;
		};
		for(const double& b : beta){
			std::cout << "\nInfo\t: Performing now geometry at velocity " << b << " and position " << x << "." << std::endl;
			f.reset(b_in, b, fN);
			f.solve(100000);
			res = f.getSurface(300);

			fp.open(baselog + std::to_string(x) + "_" + std::to_string(b) + ".txt", std::fstream::out);
			for(const auto&[wb, As] : res){
				fp << wb << "\t" << As << std::endl;
			}
			fp.close();
		}
	}
}