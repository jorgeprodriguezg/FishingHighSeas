//Code for extracting the fishing effort in grid cells
//We consider that vessels are fishing if the reported speed is lower than 5 knots
//This includes high seas and EEZ, EEZ (maybe including vessels in or approaching port) will be removed with another program
//Developer: Jorge P. Rodríguez
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<iomanip>
using namespace std;
#define lines 244638314 //Number of lines in data file
#define R 6371. //Earth radius
#define pi 4*atan(1.)

double rad(double degree){/*Degree to radians*/
  return degree*pi/180.;}

double haversine (double deglat1, double deglon1, double deglat2, double deglon2){
  double radlat1,radlat2,a,dlat,dlon;
  radlat1=deglat1*pi/180.;
  radlat2=deglat2*pi/180.;
  dlat=radlat2-radlat1;
  dlon=(deglon2-deglon1)*pi/180.;
  double sdlat=sin(dlat/2.);
  double sdlon=sin(dlon/2.);
  a=sdlat*sdlat+cos(radlat1)*cos(radlat2)*sdlon*sdlon;
  a=2*asin(sqrt(a));
  return a*R;}

int main(){
  int i,j,k,l;
  string str;
  string fname = "fname.dat";//Inset here your filename
  ifstream datain(fname);
  ofstream dataout("fishingtimeatcells.dat");
  int id[2];//id of the vessel
  double lat[2],lon[2],datat[2],speed[2],la0,lo0;//latitude, longitude, time, speed
  double time[2*2*360*180],deltat;
  int cx,cy,cell[2];
  double d;
  for(i=0;i<2*2*360*180;i++)time[i]=0.;//Initialize the effort, grid cells of 0.5º lat x 0.5º lon
  getline(datain,str);
  datain>>j>>id[1]>>lat[1]>>lon[1]>>datat[1]>>speed[1];
  cx=(int)(2.*(lon[1]+180.));
  cy=(int)(2.*(lat[1]+90.));
  cell[1]=cx+2*360*cy;
  for(i=1;i<lines-1;i++){
      id[0]=id[1];
      lat[0]=lat[1];
      lon[0]=lon[1];
      datat[0]=datat[1];
      speed[0]=speed[1];
      cell[0]=cell[1];
      datain>>j>>id[1]>>lat[1]>>lon[1]>>datat[1]>>speed[1];
      cx=(int)(2.*(lon[1]+180.));
      cy=(int)(2.*(lat[1]+90.));
      cell[1]=cx+2*360*cy;
      if(id[1]==id[0]){
        if(speed[0]<5.){//Previous point was in fishing stage
          deltat=datat[1]-datat[0];
          d=haversine(lat[0],lon[0],lat[1],lon[1]);
          if((deltat<1.)&&(d<2000.)){//location data separated by less than one day and by less than 2000 km
            time[cell[0]]=time[cell[0]]+deltat/2;//half of the time difference to previous cell
            time[cell[1]]=time[cell[1]]+deltat/2;/*half of the time difference to current cell*/}}}}




  
    for(i=0;i<2*2*360*180;i++){
         if(time[i]>0.){
           la0=-90.+(double)((int)i/(360*2))/2+1./(2*2);
           lo0=-180.+(double)((int)i%(360*2))/2+1./(2*2);
           dataout<<la0<<"\t"<<lo0<<"\t"<<time[i]<<endl;}}
    
  return 0;}

