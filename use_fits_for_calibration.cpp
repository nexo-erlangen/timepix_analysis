//Copyright: Thomas Gleixner, Tobi 2014
//Dieses Programm übernimmt die Gaussfits der Kalibrationsmessungen und gibt eine Kalibrationsdatei aus, mit der ToT zuverlässig in Energie umgewandelt werden kann.
//Weitere Erklärungen folgen...


#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <limits>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include <sstream>
#include <time.h>
using namespace std;
double variableauslesendouble(char *, char *, int);
double variableauslesendouble(char * varname, char * line, int begincomment)
{
    std::string tempstr = line;
    std::string varnames = varname;
    std::stringstream tempstream;
    double number;
    tempstream.str(tempstr.substr(varnames.length(), begincomment
            - varnames.length()));
    tempstream >> number;
    return number;
}


int main(int argc, char** argv)
{
    cerr<< "warning: running more then one instance of this program at the same time will fuck things up because temporary files will be overwritten!\n";
    //~ double pixeldata[128][128][10][3]={0};
    
        double ****pixeldata;
    if(!(pixeldata = new double***[256]))	{ cerr << "- Not enough memory .. exit ! - " << endl; exit(-1);}
    for(int i = 0; i<256;i++)
    {
	    if(!(pixeldata[i] = new double**[256]))	{ cerr << "- Not enough memory .. exit ! - " << endl; exit(-1);}
	    for(int j = 0; j<256;j++)
	    {
		    if(!(pixeldata[i][j] = new double*[15]))	{ cerr << "- Not enough memory .. exit ! - " << endl; exit(-1);}
		    for(int k = 0; k<10;k++)
		    {
			    if(!(pixeldata[i][j][k] = new double[4]))		{ cerr << "- Not enough memory .. exit ! - " << endl; exit(-1);}
		    }
	    }
    } 

    string dummy;
    int numLines2=0;
    ifstream config_in;
    string configname=argv[1];
    config_in.open(configname.c_str());
    while (std::getline(config_in,dummy))
    {
        numLines2++;
    }
    config_in.close();
    config_in.open(configname.c_str());
    cout << configname << "\n";
    for(int kalib_points=0; kalib_points<numLines2; kalib_points++)
    {
        double energy;
        string path;
        config_in >> path;
        config_in >> energy;

        ifstream input;
        cout << "reading " << path << ", energy=" << energy << " keV\n";
        input.open((path+"fitparameterhist.dat").c_str());
	cout << (path+"fitparameterhist.dat") << "\n";
        int numLines=0;
        while (std::getline(input,dummy))
        {
            numLines++;
        }
        input.close();
        input.open((path+"fitparameterhist.dat").c_str());
        for(int i=0; i<9; i++)
        {
            input >> dummy;
        }
        int x,y;
        double amplitude, sigma, sigma_error, position, position_error, dummy_double;
        for(int pixel_number=0;pixel_number<numLines-6;pixel_number++)
        {
            input >>x;
            input >>y;
            input >>sigma;
            input >>sigma_error;
            input >>position;
            input >>position_error;
            input >>dummy_double;
            input >>dummy_double;
            input >> amplitude;
		pixeldata[x][y][kalib_points][0]=position;
		pixeldata[x][y][kalib_points][1]=energy;
		pixeldata[x][y][kalib_points][2]=position_error;
           }
	input.close();
	input.open((path+"fitparameterhist.dat").c_str());
        for(int i=0; i<9; i++)
        {
            input >> dummy;
        }
     }
     
     int dummy_calib = 0;
     
    ofstream calib_out;
    calib_out.open((configname+"calib").c_str());
     
     calib_out << "#x\t #y \t #fit_a \t #fit_a_err \t #fit_b \t #fit_b_err \t #fit_c\t #fit_c_err \t #fit_t\t #fit_t_err\n";
     
    for(int x=0; x<256; x++)
    {
        for(int y=0; y<256; y++)
        {
            bool canbefitted=true;
            for(int kalib_points=0; kalib_points<numLines2; kalib_points++)
            {
			if(pixeldata[x][y][kalib_points][0]==0){canbefitted=false;}
            }
            if(canbefitted==true)
            {
                cout << "fitting kalibration curve on pixel " << x << "\t" << y << "\n";
                ofstream fitfile;
                stringstream execute;
                execute << "temp_fitfile"<< x  << "_" << y << ".dat";
                fitfile.open(execute.str().c_str());
                for(int kalib_points=0; kalib_points<numLines2; kalib_points++)
                {
                    fitfile << pixeldata[x][y][kalib_points][1] << "\t" << pixeldata[x][y][kalib_points][0] << "\t" << pixeldata[x][y][kalib_points][2] << "\n";
                }
                fitfile.close();
                execute.str("");
                execute     << "echo -e -n \""
		<< "FIT_LIMIT = 1e-6\n"
		<< "FIT_MAXITER=1000\n"
                << "f(x)=a*x+b-(c/(x-t))\n"
		//~ << "f(x)=a*x+b\n"
                << "a=2\n"
                << "b=20\n"
                << "c=46\n"
                << "t=3\n"
                << "set fit errorvariables\n"
                << "fit  f(x) \\\"temp_fitfile"<< x  << "_" << y << ".dat\\\" u 1:2:3 via a,b,c,t\n"
		//~ << "fit  f(x) \\\"temp_fitfile"<< x  << "_" << y << ".dat\\\" u 1:2:3 via a,b\n"
                << "print \\\"a=\\\",a\nprint \\\"b=\\\",b\nprint \\\"c=\\\",c\nprint \\\"t=\\\",t\n"
                << "set terminal pdf\n"
                << "set output \\\"fit" << x << "_" << y << ".pdf\\\"\n"
                << "plot  f(x), \\\"temp_fitfile"<< x  << "_" << y << ".dat\\\" u 1:2:3 with errorbars\" | gnuplot &> temp" << x  << "_" << y << ".dat";

                system(execute.str().c_str());


                char line[500];
                char * test;
                int test2;
                ifstream getPeak;
                stringstream filename;
                filename <<  "temp" << x << "_" << y  << ".dat";
                getPeak.open(filename.str().c_str());
                double fit_a, fit_b, fit_c, fit_t;

                while (getPeak.getline(line, 500))
                {
                    char *varname;
                    test = strchr(line, '#');
                    test2 = test - line;
                    if (test2 != 0)
                    {
                        varname = "a=";
                        if (strstr(line, varname) != 0)
                        {
                            fit_a = variableauslesendouble(varname, line, test2);
                        }
                        varname = "b=";
                        if (strstr(line, varname) != 0)
                        {
                            fit_b = variableauslesendouble(varname, line, test2);
                        }
                        varname = "c=";
                        if (strstr(line, varname) != 0)
                        {
                            fit_c = variableauslesendouble(varname, line, test2);
                        }

                        varname = "t=";
                        if (strstr(line, varname) != 0)
                        {
                            fit_t = variableauslesendouble(varname, line, test2);
                        }
                    }
                }
                calib_out << x << "\t" << y << "\t" << fit_a << "\t" << dummy_calib << "\t" << fit_b << "\t" << dummy_calib << "\t" << fit_c << "\t" << dummy_calib << "\t" << fit_t << "\t" << dummy_calib << "\n";
		//~ calib_out << x << "\t" << y << "\t" << fit_a << "\t" << fit_b << "\t"  << fit_c << "\t" << fit_t << "\n";
            }
        }
    }
} 
