//Code for extracting the connections in a matrix composed of high seas fishing grid cells and harbors
//We also extract the harbour associated total fishing effort
//Developer: Jorge P. Rodríguez
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
using namespace std;
#define lines 244638314 //Number of lines in data file
#define nclus 32372 //Number of high seas grid cells obtained after running filtereez.cc
#define nport 3684 //Number of ports in World Port Index
#define R 6371. //Earth radius
#define pi 4*atan(1.)
#define dmax 2000. //Maximum allowed distance between two points for not splitting the trajectory


double haversine (double deglat1, double deglon1, double deglat2, double deglon2){//Formula to calculate distances in the surface of a sphere
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
  //1: read data about trajectories,high seas grid cells (with nonzero fishing effort) and ports
  //1.1: trajectories
  cout<<"traj"<<endl;
  int id[2],cell[2],cx,cy;
  double lat[2],lon[2],time[2],sp[2];
  double la0,lo0;
  string str,fname;
  int net0[lines/100][2],E;

    //1.1: high seas grid cells
  cout<<"high seas"<<endl;
  ifstream datain2("fishingtimeatcellsnoEEZ.dat");
  int node[2*2*360*180],cl;

  double avlat[2*2*360*180],avlon[2*2*360*180];

  for(i=0;i<2*2*360*180;i++)node[i]=-1;
  for(i=0;i<nclus;i++){
    datain2>>la0>>lo0>>sp[0];
    cx=(int)(2.*(lo0+180.));
    cy=(int)(2.*(la0+90.));
    j=cx+2*360*cy;
    node[j]=i;
    avlat[i]=avlat[i]+la0;
    avlon[i]=avlon[i]+lo0;}
  cl=nclus;

  //1.2: ports
  cout<<"ports"<<endl;
  ifstream datain3("wpicoords.dat");//World Port Index
  cl=nclus;
  getline(datain3,str);
  for(i=0;i<nport;i++){
    //cout<<i<<endl;
    datain3>>j>>la0>>lo0;
    cx=(int)(2.*(lo0+180.));
    cy=(int)(2.*(la0+90.));
    j=cx+2*360*cy;
    if(node[j]==-1){//Keep only grid cells once (there are some with several ports in the same cell)
      node[j]=cl;
      avlat[cl]=la0;
      avlon[cl]=lo0;
      cl++;}}
  cout<<"Number of nodes in network "<<cl<<endl;
  int last;
  int lastpt;
  bool avai,a;
  int deg[cl];
  bool link[cl];
  double deltat,d;
  ofstream dout("linktrajectories.dat");
  bool in;
  double tharb[cl],acct;
  int lastharb;
  E=0;
  for(i=0;i<cl;i++){
    deg[i]=0;
    link[i]=0;
    tharb[i]=0.;}
  
  //Loop through the main dataset
  string fname = "fname.dat";//Inset here your filename
  ifstream datain(fname);
  int cacc;
  cout<<"Starting loop"<<endl;
  getline(datain,str);
  datain>>j>>id[1]>>lat[1]>>lon[1]>>time[1]>>sp[1];
  cx=(int)(2.*(lon[1]+180.));
  cy=(int)(2.*(lat[1]+90.));
  cell[1]=cx+2*360*cy;
  last=0;
  avai=0;
  j=cell[1];
  lastharb=-1;
  acct=0.;
  if(node[j]!=-1){
    avai=1;
    last=node[j];
    lastpt=0;
    if(node[j]>=nclus)lastharb=node[j];}
  for(i=1;i<lines-1;i++){
    id[0]=id[1];
    lat[0]=lat[1];
    lon[0]=lon[1];
    time[0]=time[1];
    sp[0]=sp[1];
    cell[0]=cell[1];
    datain>>j>>id[1]>>lat[1]>>lon[1]>>time[1]>>sp[1];
    cx=(int)(2.*(lon[1]+180.));
    cy=(int)(2.*(lat[1]+90.));
    cell[1]=cx+2*360*cy;

    //2: follow trajectories, considering links port-port, port-clus, clus-clus, only with links between consecutive events in the sequence list; directed network
    
    //Avai:{0 if there is not anything before; 1 if there was something}
    deltat=time[1]-time[0];
    j=cell[1];
    d=haversine(lat[1],lon[1],lat[0],lon[0]);
    in=(id[1]==id[0])&&(deltat<1.);
    in=(in)&&((bool)((int)(dmax/d)));//in=0 means either deltat >= 1 day or d>=2000 km, so we split the trajectory
    if(in==1){//Not broken trajectory
      // Add accumulated time in previous cell
      if((node[cell[0]]!=-1)&&(node[cell[0]]<nclus)){//Previous cell was a hotspot
        if(sp[0]<5.){//Fishing in previous cell
          acct=acct+deltat/2;}}
      
	if(node[j]!=-1){//If the visited cell is a node of our network
        	net0[E][0]=last;
        	net0[E][1]=node[j];
        	a=(avai)&&(node[j]!=last);//We have origin, and different from destination, ie no self loops
        	link[last]=link[last]||a;//Nodes in at least one link
        	link[node[j]]=link[node[j]]||a;//Nodes in at least one link
        	deg[last]=deg[last]+a;
        	E=E+a;
        	avai=1;
        	if(a==1){//a=1 means a trajectory between last and current cells
          		//If there is a trajectory arriving to a harbour, update the times
          		if(node[j]>=nclus){
           			tharb[node[j]]=tharb[node[j]]+acct;//Current harbour
            			if(lastharb!=-1)tharb[lastharb]=tharb[lastharb]+acct;//Previous harbour
            			lastharb=node[j];
            			acct=0.;//and set to 0 the fishing time
            	}
        	if((node[j]<nclus)&&(sp[1]<5.)){//This cell is a hotspot, the vessel is fishing, so it accumulates fishing effort
          		acct=acct+deltat/2;}
        	if(node[j]>=nclus)lastharb=node[j];//Do it again for the cases when avai=0, to start updating it
        	last=node[j];}
	}
    else{ //If we have new trajectories starting in cells not considered (node[cell]=-1), there will not be an available previous point
      if(lastharb!=-1)tharb[lastharb]=tharb[lastharb]+acct;
      acct=0.;
      lastharb=-1;
      if(node[j]!=-1){
        last=node[j];
        avai=1;
          if(node[j]>=nclus)lastharb=node[j];}
       else{
        avai=0;}}}
  //Loop through the main dataset finished

  cout<<"Number of links in the matrix "<<E<<endl;
  //Write the network
  int index[cl+1];
  int net1[lines/100];
  // for(i=0;i<cl;i++)cout<<deg[i]<<endl;
  cout<<"building matrix"<<endl;

  index[0]=0;
  for(i=0;i<cl;i++){
    index[i+1]=index[i]+deg[i];
    deg[i]=index[i];}

  for(i=0;i<E;i++){
    j=net0[i][0];
    k=net0[i][1];
    net1[deg[j]]=k;
    deg[j]++;}
  
  cout<<"building matrix finished"<<endl;
  int w[cl];
  ofstream dataout("matrixtrajconnections.dat");
  dataout<<"Source\tTarget\tWeight"<<endl;
  for(i=0;i<cl;i++){
    for(j=0;j<cl;j++)w[j]=0;
    for(j=index[i];j<index[i+1];j++){
      k=net1[j];
      w[k]++;}
    for(j=0;j<cl;j++){
      if(w[j]!=0)dataout<<i<<"\t"<<j<<"\t"<<w[j]<<endl;}}
  
  ofstream dataout2("matrixtraj_datanodes.dat");
  dataout2<<"Id\tLat\tLon\tFishing\tEffort"<<endl;//Fishing is 1 for high seas grid cells, 0 for harbours. Effort is the harbour effort 
  for(i=0;i<cl;i++){
    if(link[i]==1){
      if(i<nclus)dataout2<<i<<"\t"<<avlat[i]<<"\t"<<avlon[i]<<"\t1\t"<<tharb[i]<<endl;
      else dataout2<<i<<"\t"<<avlat[i]<<"\t"<<avlon[i]<<"\t0\t"<<tharb[i]<<endl;}}
    
  return 0;}

  

