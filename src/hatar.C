#include "hatar.h"


int main(){
  Point<2> * vp=new Point<2>[20000];
  int nb=0;
  char c,s[500];
  double x[500],y[500], xymin=1.e30, xymax=-1.e30;
  ifstream be ("data/proba.fig" , std::ifstream::in );
  if(be.is_open()){
    for(int i=0;i<10;i++) be.getline(s,500);
    while(!be.eof()){
      be>>x[nb];be>>y[nb];//cout<<x[nb]<<"  "<<y[nb]<<endl;
      vp[nb]=*new Point<2>(x[nb], y[nb]);
      if(x[nb]<xymin) xymin=x[nb];
      if(y[nb]<xymin) xymin=y[nb];
      if(x[nb]>xymax) xymax=x[nb];
      if(y[nb]>xymax) xymax=y[nb];
      nb++;
    }
    nb--;//egy (0,0)-t is beletesz, azt ki kell venni
    be.close();
  }else{cout<<"Nincs ilyen fájl!\n";return 0;}

  nb--;vector< Point<2> > bp(nb);
  for(int i=0;i<nb;i++) bp[i]=vp[i];
  int n=nb;
  Point<2> pont;
  Doub dy=(xymax-xymin)/10*sqrt(3.)/2.*.735;
  Doub dx=(xymax-xymin)/10*.735;

  for(int i=0;i<(xymax-xymin)/dy;i++){
    Doub py=xymin+i*dy;
    for(int j=0;j<(xymax-xymin)/dx+2;j++){
      Doub px=xymin+(j-(i % 2)*.5)*dx;
      pont=*new Point<2>(px, py);
      int wind=polywind(bp,pont);
      if(wind==1||wind==-1){
	vp[n++]=pont;
      }
    }
  }

  vector< Point<2> > racspontok(n); 
  for(int i=0;i<n;i++) racspontok[i]=vp[i];

  myDelaunay racs(racspontok,1);

  

 FILE * graphics=popen("gnuplot -persist\n","w");
  fprintf(graphics,"set term qt size 768,768 title \"Ideális folyadék áramlása két dimenzióban: példa végeselem-módszerre\"\n");
  fprintf(graphics,"unset key\n");//fprintf(graphics,"set mouse\n");
  fprintf(graphics,"set multiplot\n");
  fprintf(graphics,"set lmargin 0\n");
  fprintf(graphics,"set rmargin 0\n");
  fprintf(graphics,"set tmargin 0\n");
  fprintf(graphics,"set bmargin 0\n");
  fprintf(graphics,"set xrange[%21.6f:%21.6f]\n",xymin,xymax*1.1);
  fprintf(graphics,"set yrange[%21.6f:%21.6f]\n",xymin,xymax*1.1);



for(int i=0;i<nb;i++){
    int j=(i+1)%nb;
    int k=0,m=0,hsz=-1;
    bool felt=false;
    if(racs.linehash->get( racs.hashfn.int64(j)-racs.hashfn.int64(i),k)!=0){
      if(racs.trihash->get( racs.hashfn.int64(j) ^ racs.hashfn.int64(i) ^ racs.hashfn.int64(k),m)!=0){
	if(racs.thelist[m].stat==1){
	  racs.erasetriangle(j,i,k,-1,-1,-1);
	}else{felt=true;}
      }
    }else{felt=true;}

    if(racs.linehash->get( racs.hashfn.int64(i)-racs.hashfn.int64(j),k)!=0 ){
      if(racs.trihash->get( racs.hashfn.int64(j) ^ racs.hashfn.int64(i) ^ racs.hashfn.int64(k),m)!=0){if(racs.thelist[m].stat==1){felt=false;}}} 

    

    if(felt){
      int vi[20],ii=0;
      if(racs.pointhash->getinit(racs.hashfn.int64(i))==1){
	while(racs.pointhash->getnext(m)==1){
	  if(racs.thelist[m].stat==1){
	    int a=-1,b,ll;
	    for(ll=0;ll<3;ll++){if(racs.thelist[m].p[ll]==i){break;}}
	    a=racs.thelist[m].p[(ll+1)%3];b=racs.thelist[m].p[(ll+2)%3];
	    if(a>=nb){
	      int wind=polywind(bp,vp[a]);
	      if(wind==1||wind==-1) vi[ii++]=a;
	    }
	    if(b>=nb){
	      int wind=polywind(bp,vp[b]);
	      if(wind==1||wind==-1) vi[ii++]=b;
	    }
	    if(racs.intersect(vp[i],vp[j],vp[b],vp[a])){racs.erasetriangle(i,a,b,-1,-1,-1);}
	  }
	}
      }

      if(racs.pointhash->getinit(racs.hashfn.int64(j))==1){
	while(racs.pointhash->getnext(m)==1){
	  if(racs.thelist[m].stat==1){
	    int a=-1,b,ll;
	    for(ll=0;ll<3;ll++){if(racs.thelist[m].p[ll]==j){break;}}
	    a=racs.thelist[m].p[(ll+1)%3];b=racs.thelist[m].p[(ll+2)%3];
	    if(racs.intersect(vp[i],vp[j],vp[b],vp[a])){racs.erasetriangle(j,a,b,-1,-1,-1);}
	    for(int s=0;s<ii;s++){if((vi[s]==a||vi[s]==b)&&felt){
		racs.storetriangle(i,j,vi[s]);felt=false;break;
	      }}
	    
	  }
	}
      }
      
    
    }
    

}

  
  

  fprintf(graphics,"plot '-' w p pt 7 lc 1\n");  
  for(int i=0 ;i<n;i++){ 
    fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]);  
  }
  fprintf(graphics,"e\n");

for(int i=0 ;i<nb;i++){ 
  fprintf(graphics,"plot '-' w p pt 7 lc 3\n"); 
  fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]); 
  fprintf(graphics,"e\n");
  fflush(graphics);
  //int z;cin>>z;
  fprintf(graphics,"plot '-' w p pt 7 lc 2\n"); 
  fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]); 
  fprintf(graphics,"e\n");  
  }


  fflush(graphics);

  //goto vege;

  /*Point<2> tesztpont=*new Point<2>(10000.,5000.);

fprintf(graphics,"plot '-' w l lc 2\n");
  for(int j=0;j<nb;j++){
    Bool u=true;
    for(int k=0;k<nb;k++){ 
      if(racs.intersect(tesztpont,bp[j],bp[k],bp[(k+1)%nb])){u=false;break;}  
    }
    if(u){
      fprintf(graphics,"%21.6f %21.6f\n",tesztpont.x[0],xymax-tesztpont.x[1]); 
      fprintf(graphics,"%21.6f %21.6f\n",bp[j].x[0],xymax-bp[j].x[1]); 
      fprintf(graphics,"\n");
    }
  }
  fprintf(graphics,"e\n");*/

  //az eredeti program tesztje

  //fprintf(graphics,"plot '-' w l lc 3\n");
fprintf(graphics,"plot '-' w l lw 3 lc 2\n");
  
  for(int i=0 ;i<=nb;i++){ 
    fprintf(graphics,"%21.6f %21.6f\n",vp[i%nb].x[0],xymax-vp[i%nb].x[1]);  
  }
  fprintf(graphics,"e\n");

  for(int j=0;j<racs.thelist.size();j++){ 
      if(racs.thelist[j].stat==1){

	fprintf(graphics,"plot '-' w l lc 1\n"); 
	 for(int k=0;k<=3;k++){fprintf(graphics,"%21.6f %21.6f\n",racspontok[racs.thelist[j].p[k%3]].x[0],xymax-racspontok[racs.thelist[j].p[k%3]].x[1]);}
	 fprintf(graphics,"\n");fprintf(graphics,"e\n");
	 fflush(graphics);
	 //int z;cin>>z;
	 fprintf(graphics,"plot '-' w l lc 3\n"); 
	 for(int k=0;k<=3;k++){fprintf(graphics,"%21.6f %21.6f\n",racspontok[racs.thelist[j].p[k%3]].x[0],xymax-racspontok[racs.thelist[j].p[k%3]].x[1]);}
	 fprintf(graphics,"\n");fprintf(graphics,"e\n");
    }
    
  }
  


  /*for(int i=0;i<nb;i++){
    int j=(i+1)%nb;
    int k=0,m=0;
    if(racs.linehash->get( racs.hashfn.int64(j)-racs.hashfn.int64(i),k)!=0){
      if(racs.trihash->get( racs.hashfn.int64(j) ^ racs.hashfn.int64(i) ^ racs.hashfn.int64(k),m)!=0){
	//cout<<racs.thelist[m].d[0]<<"  "<<racs.thelist[m].d[1]<<"  "<<racs.thelist[m].d[2]<<"  "<<racs.thelist[m].stat<<endl;
	if(racs.thelist[m].stat==1){
	  if(k>=nb){
	    int wind=polywind(bp,vp[k]);
	    if(wind==1||wind==-1) continue;
	  }
	  fprintf(graphics,"plot '-' w l lc 5\n");
	  fprintf(graphics,"%21.6f %21.6f\n",vp[j].x[0],xymax-vp[j].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[k].x[0],xymax-vp[k].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[j].x[0],xymax-vp[j].x[1]);
	  fprintf(graphics,"\n");fprintf(graphics,"e\n");
	  fflush(graphics);int z;cin>>z;
	  fprintf(graphics,"plot '-' w l lc 1\n");
	  fprintf(graphics,"%21.6f %21.6f\n",vp[j].x[0],xymax-vp[j].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[k].x[0],xymax-vp[k].x[1]);
	  fprintf(graphics,"%21.6f %21.6f\n",vp[j].x[0],xymax-vp[j].x[1]);
	  fprintf(graphics,"\n");fprintf(graphics,"e\n");

	  racs.erasetriangle(j,i,k,-1,-1,-1);
	}
      }
    }
    if(racs.linehash->get( racs.hashfn.int64(j)-racs.hashfn.int64(i),k)==0 ){

    }
    }
  */

  goto vege;
   
  fprintf(graphics,"plot '-' w l lc 3\n");
  for(int j=0;j<racs.thelist.size();j++){ 
      if(racs.thelist[j].stat==1){
	 Bool u=true;
	 for(int k=0;k<=3&&u;k++){ 
	   for(int m=0;m<nb;m++){ 
	     if(racs.intersect(racspontok[racs.thelist[j].p[k%3]],racspontok[racs.thelist[j].p[(k+1)%3]],bp[m],bp[(m+1)%nb])){u=false;break;}  
	   }
	   if(!u) break;
	   Point<2> tesztpont=*new Point<2>(.5*racspontok[racs.thelist[j].p[k%3]].x[0]+.5*racspontok[racs.thelist[j].p[(k+1)%3]].x[0],.5*racspontok[racs.thelist[j].p[k%3]].x[1]+.5*racspontok[racs.thelist[j].p[(k+1)%3]].x[1]);
	   
	   if(polywind(bp,tesztpont)==0){ 
	     for(int m=0;m<nb;m++){ u=false;
	       Point<2> tesztpont2=*new Point<2>(.5*bp[m].x[0]+.5*bp[(m+1)%nb].x[0],.5*bp[m].x[1]+.5*bp[(m+1)%nb].x[1]);
	       if(fabs(tesztpont.x[0]-tesztpont2.x[0])+fabs(tesztpont.x[1]-tesztpont2.x[1])<1.e-5){u=true;break;}  
	     }
	   }
	 }
	 
	 
	 for(int k=0;k<=3&&u;k++){fprintf(graphics,"%21.6f %21.6f\n",racspontok[racs.thelist[j].p[k%3]].x[0],xymax-racspontok[racs.thelist[j].p[k%3]].x[1]);}
	 if(u) fprintf(graphics,"\n");
      
    }
    
  }
  fprintf(graphics,"e\n");



 


 vege:

  fprintf(graphics,"plot '-' w p pt 7 lc 1\n");  
  for(int i=0 ;i<n;i++){ 
    fprintf(graphics,"%21.6f %21.6f\n",vp[i].x[0],xymax-vp[i].x[1]);  
  }
  fprintf(graphics,"e\n");

  fprintf(graphics,"q\n");
  fflush(graphics);
  pclose(graphics);
}

