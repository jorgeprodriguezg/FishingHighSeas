//Code for filtering the high seas grid cells after computing the fishing effort
//Developer: Jorge P. Rodríguez
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
using namespace std;
#define lines 56774 //number of lines in the output from effort_allcells.cc
#define lineshs 98339 //number of lines in the file highseas.dat
int main(){
  int i,j,k;
  string str;
  ifstream datain("fishingtimeatcells.dat");
  ifstream datain2("highseas.dat");//Latitudes and longitudes of the cells in high seas
  ofstream dataout("fishingtimeatcellsnoEEZ.dat");
  double lat,lon,data;
  bool highsea[2*2*360*180],h;
  int cx,cy,cell,count;
  getline(datain2,str);
  count=0;
  for(i=0;i<2*2*360*180;i++)highsea[i]=0;
  for(i=0;i<lineshs-1;i++){
    datain2>>lat>>lon;
    cx=(int)(2.*(lon+180.));
    cy=(int)(2.*(lat+90.));
    cell=cx+2*360*cy;
    highsea[cell]=1;
    count++;}
  cout<<count<<endl;
  
  double max=-1.,min=10000000.;
  
  for(i=0;i<lines;i++){
    datain>>lat>>lon>>data;
    cx=(int)(2.*(lon+180.));
    cy=(int)(2.*(lat+90.));
    cell=cx+2*360*cy;
    if(highsea[cell]==1){
      dataout<<lat<<"\t"<<lon<<"\t"<<data<<endl;
      max=fmax(max,data);
      min=fmin(min,data);
    }}


  cout<<max<<"\t"<<min<<endl;//Shows the maximum and the minimum effort at high seas grid cells
  return 0;}

