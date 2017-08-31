// Copyright: Benedikt Bergmann
// Setup gives .root-Files as output. This script includes tools to evaluate and analyze this data. Most things are made by Benedikt Bergmann, however Daniel may decide to add or restructure to enhance usability. Changes by him are marked as such by comment-lines.
//To load, open a ROOT-session and type ".L EvaluateCoincidences.cpp". Now the fruits of our humble work are at your disposal! Use them wisely...


using namespace std;

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void defineHistogram(TH1* h, std::string x, std::string y, std::string title, bool adjust_range, double titleSize = 0.05, bool rebin = false, double average_entries_per_bin = 10.)
{
	h->GetXaxis()->SetTitle(x.c_str());
	h->GetYaxis()->SetTitle(y.c_str());
	h->SetStats(kFALSE);
	h->SetLineWidth(2);
	h->SetTitle(title.c_str());
	h->GetYaxis()->SetDecimals();
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitleSize(titleSize);
	h->GetYaxis()->SetTitleSize(titleSize);
	h->GetXaxis()->SetLabelSize(titleSize);
	h->GetYaxis()->SetLabelSize(titleSize);
	h->SetTitleSize(titleSize);
	
	h->GetXaxis()->SetNdivisions(503, true);
	
	h->GetYaxis()->SetNdivisions(505, true);

	int lastBin = h->GetNbinsX();
	if(adjust_range == true) 
	{
		for(int i = 1; i <= h->GetNbinsX(); i++) 
		{
			if(h->GetBinContent(i) > 0)
				lastBin = i;
		}
		lastBin= (int) lastBin * 1.2;

		if(rebin == true) 
		{
			int bin_factor = (average_entries_per_bin * lastBin / h->GetEntries());
			if(bin_factor > 1) 
			{
				h->Rebin(bin_factor);
				lastBin = lastBin / bin_factor * 1.2;
			}
		}
		h->GetXaxis()->SetRange(0,lastBin);
	}
 }

TFile* f;
TTree* t_data;
TTree* t_info;

//Opens the ROOT file and loads the tree with the clustered data!
void LoadTree(std::string filename)
{
	f = new TFile(filename.c_str(), "OPEN");
	t_data = (TTree*) f->Get("clusteredData");
   // cout << (TTree)f->Get(f->GetListOfKeys()->At(0)->GetName()) << endl;
}

//Daniel: Prints out Info

void GiveInfo(std::string filename)
{
    f = new TFile(filename.c_str(), "OPEN");
    t_info = (TTree*) f->Get("InfoTree_layer");
    t_info->SetMakeClass(1);
    
    double serie_start_time;
    double acq_time;
    string chipboard_id;
    double start_time;
    double threshold;
    
    t_info->GetBranch("serie_start_time")->SetAddress(serie_start_time);
    t_info->GetBranch("acq_time")->SetAddress(acq_time);
    t_info->GetBranch("chipboard_id")->SetAddress(chipboard_id);
    t_info->GetBranch("start_time")->SetAddress(start_time);
    t_info->GetBranch("threshold")->SetAddress(threshold);
    
    cout << "Serie Start: " << serie_start_time << endl;
    cout << "Acquisition time: " << acq_time << endl;
    cout << "Chipboard: " << chipboard_id << endl;
    cout << "Start time (differs from Serie Start!): " << start_time << endl;
    cout << "Threshold: " << threshold << endl;
}

//Plots a set of individual events based on some attributes
//Daniel: Also TH1F to compare ToA Values of certain columns. Seems to be a bug here
//Daniel: Also added: couts to ensure it working. At the end it calculates runtime, which is really neat, but i think we need more control
//Daneil: Also added: Counting of cluster of conversion electrons; Method of counting is questionable, do not trust result; preferably use cuts in root-command and take NumberOfEvents
void PlotXYC(std::string cuts = "", int eventNo_afterCut = -1, int _maximumNumberOfEntries = -1, std::string _title = "")
{
	cout<<"begin program!"<<endl;
	time_t start_time = time(0);
	cout<<"maxNumberOfEvents: "<<_maximumNumberOfEntries<<endl;
	
	gStyle->SetPaintTextFormat("4.1f");
	cuts = (cuts).c_str();

	bool energy_cali = false;

	set_plot_style();

    int countConversion = 0;
    
	const int maxNumberOfClusters = 65536;
	int clstrSize;
	short PixX[maxNumberOfClusters];
	short PixY[maxNumberOfClusters];
	float ToT_keV[maxNumberOfClusters];
	int ToT[maxNumberOfClusters];
	double ToA[maxNumberOfClusters];
    double clstrVolume_keV;
	long int _eventNo;
	if(t_data->GetListOfBranches()->Contains("ToT_keV")){ t_data->GetBranch("ToT_keV")->SetAddress(ToT_keV); energy_cali = true;}

	double minToA;
	double maxToA;
	int triggerNo;
	long coinc_group;

	t_data->GetBranch("clstrSize")->SetAddress(&clstrSize);
	t_data->GetBranch("PixX")->SetAddress(PixX);
	t_data->GetBranch("PixY")->SetAddress(PixY);
	t_data->GetBranch("ToA")->SetAddress(ToA);
	t_data->GetBranch("ToT")->SetAddress(ToT);
	t_data->GetBranch("min_ToA")->SetAddress(&minToA);
	t_data->GetBranch("coincidence_group")->SetAddress(&coinc_group);
    t_data->GetBranch("clstrVolume_keV")->SetAddress(clstrVolume_keV);

	TH2F* pix = new TH2F("pixelmatrix", "", 256, -0.5,255.5,256,-0.5,255.5);
	TH2D* pix_ToA = new TH2D("pixel_ToA", "", 256, -0.5,255.5,256,-0.5,255.5);
	TH2F* pix_hits = new TH2F("pix_hits", "", 256,-0.5,255.5, 256,-0.5,255.5);
    TH1D* pix_ToA_1D_bad = new TH1D("pixel_ToA bad","",256,-0.5,255.5);
    TH1D* pix_Toa_1D_good = new TH1D("pixel_ToA good","",256,-0.5,255.5);
    TH2D* pix_ToA_test = new TH2D("test","",256,-0.5,255.5,256,-0.5,255.5);
    TH2D* pix_ToA_test2 = new TH2D("test2","",256,-0.5,255.5,256,-0.5,255.5);

	int numberOfEvents = 0; 
	int number_pixels = 0;

	double min_TOA;
	double max_TOA;
    
    //ofstream strtoagood;
    //strtoagood.open("toagood2.dat");
    //ofstream strtoabad;
    //strtoabad.open("toabad2.dat");
    //ifstream toacorrect;
    //toacorrect.open("toacorrect.dat");
    //double toacorrection[255];
    //for (int v=0; v<256; v++) {
    //    toacorrect>>toacorrection[v];
    //}

	t_data->Draw(">>myList", cuts.c_str(), "entryList");
    
    cout<<"Starting big loop now!" << endl;

//Daniel: To solve offset in ToA histo, 'dif' wasn't plotted. Instead ToA[j] is now in.
//Daniel: Redaction: dif is now in again, as well as several testplots, but may be omitted for interest. :P
    
    double ToA_dif_raw[256][256] = { 0 };
    
	for(unsigned int i = 0; i < myList->GetN(); i++)
	{
		if(eventNo_afterCut == -1 || i == eventNo_afterCut)
		{
			t_data->GetEntry(myList->GetEntry(i));
			cerr.precision(24);
			if(minToA < min_TOA || i == 0){
                min_TOA = minToA;
			}

			numberOfEvents++;
			for(short int j = 0; j < clstrSize; j++)
			{
				number_pixels++;
				if(energy_cali == false)
					pix->Fill(PixX[j],PixY[j],ToT[j]);
				else{ 
					if(ToT_keV[j] < 0)
						ToT_keV[j] = 0.1;
					pix->Fill(PixX[j],PixY[j],ToT_keV[j]);
				}

				double dif;
                //if(PixX[j]==193 || PixX[j]==194){ToA[j] = ToA[j] - toacorrection[PixY[j]];}
				dif = ToA[j] - minToA + 0.00001;
                
                /*
                if (PixX[j]==193 || PixX[j]==194) {
                    pix_ToA->Fill(PixX[j], PixY[j], dif-toacorrection[PixY[j]]);
                }
                */
               // else{
                    pix_ToA->Fill(PixX[j], PixY[j], dif);
                ToA_dif_raw[PixX[j]][PixY[j]] = ToA_dif_raw[PixX[j]][PixY[j]] + dif;
               // }
				
				pix_hits->Fill(PixX[j], PixY[j], 1);
                pix_ToA_test->Fill(PixX[j],PixY[j], ToA[j]);
			}
			if(i == eventNo_afterCut){break;}
		}
		if(numberOfEvents >= _maximumNumberOfEntries && _maximumNumberOfEntries > 0){break;}
        //cout << "I'm filling histos. Right now i am at " << i << "from " << myList.GetN() << endl;
        
        if (clstrSize>2 && clstrVolume_keV>750 && clstrVolume_keV < 900) {
            countConversion++;
        }
        
	}
    
    cout << "Ending big loop now!" << endl;
    
    for (int k=0; k<256; k++) {
        double bad_bin1 = pix_ToA->GetBinContent(193,k);
        double bad_bin2 = pix_ToA->GetBinContent(194,k);
        double good_bin1 = pix_ToA->GetBinContent(192,k);
        double good_bin2 = pix_ToA->GetBinContent(195,k);
        pix_ToA_1D_bad->Fill(k,bad_bin1);
        pix_Toa_1D_good->Fill(k,good_bin1);
        //strtoagood << good_bin2 << endl;
        //strtoabad << bad_bin2 << endl;
    }
    
    //To do ToAcorrection in program:
    for (int z=0; z<256 ; z++) {
        double ex3 = ((ToA_dif_raw[192][z] + ToA_dif_raw[193][z])/2) - ((ToA_dif_raw[191][z] + ToA_dif_raw[194][z])/2); // Mittelwert Differenz --> "mittlerer Offset"
        ToA_dif_raw[192][z]= ToA_dif_raw[192][z] - ex3;
        ToA_dif_raw[193][z]= ToA_dif_raw[193][z] - ex3;
        //ToA_dif_raw[192][z]= ToA_dif_raw[191][z];
        //ToA_dif_raw[193][z]= ToA_dif_raw[194][z];
    }
    
    for (int m=0; m<256; m++) {
        for (int n=0; n<256; n++) {
            pix_ToA_test2->Fill(m, n, ToA_dif_raw[m][n]);
        }
    }
    // -----
    
	TCanvas* c = new TCanvas(("c_"+cuts).c_str(), "result canvas", 900, 1300);
    TCanvas* c_toa = new TCanvas("ToA","control canvas",1300,620);
	c->Divide(1,2);
    c_toa->Divide(2,1);

	for(short ic = 1; ic <= 2; ic++) 
	{
		c->GetPad(ic)->SetLeftMargin(0.12);
		c->GetPad(ic)->SetBottomMargin(0.12);
		c->GetPad(ic)->SetRightMargin(0.1);
		gPad->Modified();
	}

	TPaveText* pt = new TPaveText(0.83, 0.91 ,0.97, 0.99, "brNDC");
	pt->SetFillColor(kWhite);
	if(energy_cali == false)
		defineHistogram(pix, "Pixel No. X", "Pixel No.Y", "ToT (counts)", false);
	else
		defineHistogram(pix, "Pixel No. X", "Pixel No.Y", "ToT (keV)", false);

	defineHistogram(pix_ToA, "Pixel No. X", "Pixel No.Y", "#DeltaToA (ns)", false);

	if(energy_cali == true)
	{
		pix->SetMinimum(0);
		gPad->SetLogz();
	}
	//pix_ToA->SetMinimum(min_TOA);
	//pix_ToA->SetMaximum(max_TOA);

	c->cd(1);
	pix->Draw("COLZ9");
	c->cd(2);
	pix_ToA_test2->Draw("COLZ9");
    
    c_toa->cd(1);
    pix_Toa_1D_good->Draw();
    c_toa->cd(2);
    pix_ToA_1D_bad->Draw();
    //strtoagood.close();
    //strtoabad.close();

	time_t endtime = time(0);
	cout<<"Found "<<numberOfEvents<<" Events! and number pixels = "<<number_pixels<<endl;
    cout<<"Found "<<countConversion<<" conversion electrons!"<<endl;
    cout<<"runtime: "<<endtime-start_time << endl;
    cout<<"Disclaimer: In ToA plot, a thin stripe (2 pixel rows) at around x=190 is standing out. The source of this problem has yet to be determined."<<endl;
}


//Select coincidence group with 2 clusters and look at their energy and time.
//Daniel: changed cuts to measure conversion electrons in coincidence with K-fluorescences. Plots are not to be trusted now!
void ScatterPlotEnergiesCoincidenceGroups(std::string filename){
    int countConversion =0;
	f = new TFile(filename.c_str(), "OPEN");
	TTree* t = (TTree*) f->Get("clusteredData");

	set_plot_style();

	float clstrVolume_keV;
	long coinc_group_no;
	int clstrSize;
	double clstrTime;

	//load the needed branches;
	t->GetBranch("clstrVolume_keV")->SetAddress(&clstrVolume_keV); //energy per cluster
	t->GetBranch("coincidence_group")->SetAddress(&coinc_group_no); 
	t->GetBranch("clstrSize")->SetAddress(&clstrSize);
	t->GetBranch("min_ToA")->SetAddress(&clstrTime);

	//Create and entry list of interesting events;
	t->Draw(">>myList", "coincidence_group_size == 2", "entryList");

	//containers for the coincidence group results
	double energies[2] = {0};
	double time[2] = {0};
	int size[2] = {0};
	
	//Define the histrograms
	TH1F* h_dt = new TH1F("h_dt", "", 100, 0, 156);
	TH1F* h_e = new TH1F("h_e", "", 150, 0, 150);
	//scatter plot of energies of event 1 and event2 
	TH2F* h_e1e2 = new TH2F("h_e1e2", "", 1000, 0, 1000, 1000, 0, 1000);

	long coinc_group_no_temp = -1;
	//int maxEvents = 10000;

	int j = 0;
	for(long i = 0; i < myList->GetN(); i++){
		//Go through the entry list and do some analysis;
		t->GetEntry(myList->GetEntry(i));

		if(coinc_group_no == coinc_group_no_temp){
			size[1] = clstrSize;
			time[1] = clstrTime;
			energies[1] = clstrVolume_keV;
            
//Daniel: HIER DIE CUTS FUER COINCIDENCES WOHOOO!!! :-< :-< :-< :-< :-< :-< :-< :-< :-< <--DAS SIND SCHEREN SCHNIP SCHNIP
            
/*
			if(size[0] <= 2 && size[1] <= 2 && energies[0] + energies[1] < 18){
				h_e1e2->Fill(energies[0], energies[1]);
				h_dt->Fill(fabs(time[1] - time[0]));
				h_e->Fill(energies[0] + energies[1]);
                countConversion++;
			}
 */
            if(size[0] <= 2 && size[1] > 2 && energies[0] <10 && energies[1] > 750 && energies[1] < 900){
                h_e1e2->Fill(energies[0], energies[1]);
                h_dt->Fill(fabs(time[1] - time[0]));
                h_e->Fill(energies[0] + energies[1]);
                countConversion++;
            }
		}
		else{
			energies[0] = clstrVolume_keV;
			size[0] = clstrSize;
			time[0] = clstrTime;
			coinc_group_no_temp = coinc_group_no;
		}
	}

	//Define Canvas and plot the results
	TCanvas** c = new TCanvas*[3];
	for(int j = 0; j < 3; j++){
		c[j] = new TCanvas();
		c[j]->SetBottomMargin(0.11);
		c[j]->SetLeftMargin(0.11);
	}

	c[0]->cd();
	defineHistogram(h_e1e2, "Energy [keV]", "Energy [keV]", "Coinc. Group Size: 2, cluster size <= 2", false);
	h_e1e2->GetXaxis()->SetTitleOffset(0.9);
	h_e1e2->GetYaxis()->SetTitleOffset(0.9);
	h_e1e2->Draw("COLZ");
	
	c[1]->cd();
	defineHistogram(h_dt, "#Deltat (ns)", "Number of events", "Coinc. Group Size: 2, cluster size <= 2", false);
	h_dt->GetXaxis()->SetTitleOffset(0.9);
	h_dt->GetYaxis()->SetTitleOffset(0.9);
	h_dt->Draw("");

	c[2]->cd();
	defineHistogram(h_e, "#SigmaE (keV)", "Number of events", "Coinc. Group Size: 2, cluster size <= 2", false);
	h_e->GetXaxis()->SetTitleOffset(0.9);
	h_e->GetYaxis()->SetTitleOffset(0.9);
	h_e->Draw("");
    cout << "Found " << countConversion << " conversion electrons in coincidence with fluorescences!" << endl;
}


//Daniel: nlines_energy_tot is for binning. OG:1000
int nlines_energy_tot = 4000;
//Daniel: max_energy_tot is maximum energy. OG:2000
int max_energy_tot = 2000;


//TH1F* PlotEnergySpectrum(std::string fn, int nbins = 1000, float xmin = 0, float xmax = 2000){
TH1F* PlotEnergySpectrum(std::string fn, int nbins = 4000, float xmin = 0, float xmax = 2000){
	TFile* f = new TFile(fn.c_str(), "OPEN");
	TTree* t = (TTree*) f->Get("clusteredData");
	TH1F* h = new TH1F("h", "", nbins, xmin, xmax);
	//h->SetDirectory(0);
	t->Draw("clstrVolume_keV >> h");
	return h;
	//f->Close();
	//delete f;
}

//Creates and saves the integrated energy spectrum
//Daniel: Also introduces function to integrate over x-values. And prints histogram content into .dat-File




TH1F* h_tot = new TH1F("h_all", "Mn-54: 10 h", nlines_energy_tot, 0, max_energy_tot);
void SaveEnergySpectrum(int run_start, int run_end, std::string out = ""){
    ofstream energy_output;
    energy_output.open("Energiespektrum_total.dat");
	TCanvas* c = new TCanvas();
	c->SetLogy();
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	
	TH1F* h_run;

	for(int i = run_start; i <= run_end; i++){
		//TH1F* h = new TH1F("h" , "", 1000, 0, 2000);
		std::stringstream filename;
		filename << "run" << std::setw(2) << std::setfill('0') << i << ".root";
		cerr << filename.str().c_str() << endl;
		h_run = PlotEnergySpectrum(filename.str().c_str());
		h_tot->Add(h_run);
		h_run->Clear();
	}
	defineHistogram(h_tot, "Cluster volume (keV)", "No. of events", "Mn-54: 10 h", false);
	

	h_tot->Draw();
	if(out.compare("") != 0){
		c->Print( (out+".pdf").c_str(), "pdf");
		TFile* f_out = new TFile( (out+".root").c_str(), "RECREATE");
		h_tot->Write("h_tot", TObject::kOverwrite);
		f_out->Close();
		delete f_out;
	}
    for (int p=0; p <= nlines_energy_tot; p++) {
        energy_output << p*(max_energy_tot/nlines_energy_tot) << "\t" << h_tot -> GetBinContent(p) << endl;
    }
    energy_output.close();
}

//Daniel: Gives possibility to integrate over x-Values in total energy spectrum

void IntegralEnergy(int xmin, int xmax)
{
    
    TAxis *axis = h_tot->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h_tot->Integral(bmin,bmax);
    cout << integral << endl;
    
}

//Daniel: RETRACTEDRETRACTEDRETRACTEDRETRACTED
//Daniel: The SaveEnergySpectrum-function can sum up all runXX.root Files and evaluate them.
// In the following section, the same tactic as used in this function will be used for PlotXYC with different, potentially interesting Cuts.

/*
void PlotSmallCoincidences(int run_start, int run_end, std::string out = "")
{
    TCanvas* c = new TCanvas();
    c->SetLogy();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    
    TH1F* h_run;
    
    for(int i = run_start; i <= run_end; i++){
        //TH1F* h = new TH1F("h" , "", 1000, 0, 2000);
        std::stringstream filename;
        filename << "run" << std::setw(2) << std::setfill('0') << i << ".root";
        cerr << filename.str().c_str() << endl;
        h_run = PlotEnergySpectrum(filename.str().c_str());
        h_tot->Add(h_run);
        h_run->Clear();
    }
    defineHistogram(h_tot, "Cluster volume (keV)", "No. of events", "Mn-54: 10 h", false);
    
    
    h_tot->Draw();
    if(out.compare("") != 0){
        c->Print( (out+".pdf").c_str(), "pdf");
        TFile* f_out = new TFile( (out+".root").c_str(), "RECREATE");
        h_tot->Write("h_tot", TObject::kOverwrite);
        f_out->Close();
        delete f_out;
    }
}
*/

//Daniel: RETRACTEDRETRACTEDRETRACTEDRETRACTED

//-------------------------------------------------------------------------------------------



