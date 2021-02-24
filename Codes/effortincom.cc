//Code for extracting the global fishing network between high seas fishing provinces and harbors
//Developer: Jorge P. Rodríguez
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
using namespace std;

#define lines 244638314 ////Number of lines in data file
#define nport 3684 //Number of ports in World Port Index
#define R 6371.
#define pi 4*atan(1.)
#define dmax 2000.
#define lined0 33529 //number of nodes after extracting the main provinces. 
//Provinces are computed applying Infomap to the output from matrixtrajectories.cc
#define ncoms 14 //number of fishing provinces


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
  int i,j,k,l,m;
  //1: read data about trajectories,hotspots and ports
  //1.0: read data about communities
  string str,fname;
  ifstream d0("idscoordseffortcoms.dat");//File specifying the id, coordinates, fishing (1 high seas, 0 port),effort (for harbours) and marine province
  
  string fname = "fname.dat";//Inset here your filename
  ifstream datain(fname);
  int id[2],cell[2],cx,cy;
  double lat[2],lon[2],time[2],sp[2];
  double la0,lo0;
  int com[2*2*360*180];
  bool port[2*2*360*180];
  bool fi;
  double effo;

  //1.1: hotspots and ports(and their coms)
  for(i=0;i<2*2*360*180;i++){
    com[i]=-1;
    port[i]=0;}
  getline(d0,str);
  for(i=1;i<lined0;i++){
    d0>>j>>la0>>lo0>>fi>>effo>>k;
    if(fi==1){
      cx=(int)(2.*(lo0+180.));
      cy=(int)(2.*(la0+90.));
      j=cx+2*360*cy;
      com[j]=k;}
    else{
      cx=(int)(2.*(lo0+180.));
      cy=(int)(2.*(la0+90.));
      j=cx+2*360*cy;
      port[j]=1;//We do not keep province for ports
    }}
  j=0;
  k=0;
  for(i=0;i<2*2*360*180;i++){
    j+=(bool)(com[i]!=-1);}
  cout<<"number of high seas grid cells "<<j<</*"\t"<<k<<*/endl;
    
  j=0;
  for(i=0;i<2*2*360*180;i++){
    j+=(bool)port[i];}
  cout<<"number of port grid cells "<<j<<endl;

  int last;
  int lastpt;
  bool avai,a;
  double deltat,d;
  bool in;
  double tharb[2*2*360*180][ncoms],acct[ncoms];
  int lastharb;
  for(i=0;i<2*2*180*360;i++){
    for(j=0;j<ncoms;j++){
      tharb[i][j]=0.;}}
  cout<<"Loop through dataset"<<endl;


    getline(datain,str);
    datain>>j>>id[1]>>lat[1]>>lon[1]>>time[1]>>sp[1];
    cx=(int)(2.*(lon[1]+180.));
    cy=(int)(2.*(lat[1]+90.));
    cell[1]=cx+2*360*cy;
    j=cell[1];
    lastharb=-1;
    for(m=0;m<ncoms;m++)acct[m]=0.;
    if(port[j]==1)lastharb=j;
    for(i=1;i<lines-1;i++){
        id[0]=id[1];
        lat[0]=lat[1];
        lon[0]=lon[1];
        time[0]=time[1];
        sp[0]=sp[1];
        cell[0]=cell[1];
        datain>>j>>id[1]>>lat[1]>>lon[1]>>time[1]>>sp[1];
        deltat=time[1]-time[0];
        cx=(int)(2*(lon[1]+180.0));
        cy=(int)(2*(lat[1]+90.0));
        cell[1]=cx+2*360*cy;
        j=cell[1];
        d=haversine(lat[1],lon[1],lat[0],lon[0]);
        in=(id[1]==id[0])&&(deltat<1.);
        in=(in)&&((bool)((int)(dmax/d)));
        if(in==1){//Not broken trajectory
        // Add accumulated time in previous cell(if it was a hotspot)
            k=com[cell[0]];
            if((k!=-1)&&(sp[0]<5.)){//Previous cell was a hotspot:then it accumulates fishing effort
              acct[k]+=deltat/2;}
	
            if(port[j]==1){//Arrives to a harbour: then have to update last harbour and current
                for(m=0;m<ncoms;m++){
                    tharb[j][m]+=acct[m];//Current harbour update
                    if(lastharb!=-1)tharb[lastharb][m]+=acct[m];//If there is a last harb, update it
                    acct[m]=0.;}
                lastharb=j;}
            if((com[j]!=-1)&&(sp[1]<5.)){//Current cell is a hotspot, so it accumulates fishing effort
              acct[com[j]]+=deltat/2;}}
        else{//Broken trajectories: have to distrib the accumulated for previous trajectory and clear it for the new trajectory
            for(m=0;m<ncoms;m++){
              if(lastharb!=-1){//In the previous trajectory, there is some accumulated  effort after last harbour
                //cout<<"WARNING "<<lastharb<<endl;
                tharb[lastharb][m]+=acct[m];}
                acct[m]=0.;}
            lastharb=-1;
            if(port[j]==1)lastharb=j;}}
  //Loop through all the dataset finished

  cout<<"Starting to write output"<<endl;
  ofstream dataout("globalfishingnetwork.dat");//This will be a matrix, rows are harbors, columns are a)coordinates, b)fishing effort in each province
  double sumeffort;
  for(i=0;i<2*360;i++){
    for(k=0;k<2*180;k++){
      lo0=(double)i/2-180.+0.25;
      la0=(double)k/2-90.+0.25;
      l=i+k*2*360;
      if(port[l]==1){
	sumeffort=0.;
	dataout<<la0<<"\t"<<lo0<<"\t";
	for(j=0;j<ncoms;j++){
	  dataout<<tharb[l][j]<<"\t";
	  sumeffort+=tharb[l][j];}
	dataout<<sumeffort<<endl;}}}
				 
    
  return 0;}

  

