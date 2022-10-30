#include<delaunay.h>
struct myDelaunay: Delaunay{

  myDelaunay(vector<Point<2> > &pvec, Int options = 0);
  Bool intersect(const Point<2> &A, const Point<2> &B, const Point<2> &C, const Point<2> &D);
  Mhash<Ullong,Int,Nullhash> *pointhash;
  Int storetriangle(Int a, Int b, Int c);
  void erasetriangle(Int a, Int b, Int c, Int d0, Int d1, Int d2);

};

Bool myDelaunay::intersect(const Point<2> &A, const Point<2> &B, const Point<2> &C, const Point<2> &D){
  double vecprod1,vecprod2,prod1,prod2;
  vecprod1=(A.x[0]-B.x[0])*(C.x[1]-B.x[1])-(A.x[1]-B.x[1])*(C.x[0]-B.x[0]);
  vecprod2=(A.x[0]-B.x[0])*(D.x[1]-B.x[1])-(A.x[1]-B.x[1])*(D.x[0]-B.x[0]);
  prod1=vecprod1*vecprod2;
  vecprod1=(C.x[0]-D.x[0])*(B.x[1]-D.x[1])-(C.x[1]-D.x[1])*(B.x[0]-D.x[0]);
  vecprod2=(C.x[0]-D.x[0])*(A.x[1]-D.x[1])-(C.x[1]-D.x[1])*(A.x[0]-D.x[0]);
  prod2=vecprod1*vecprod2;
  if(prod1<0 && prod2 <0) return true;
  return false;
}


  Int myDelaunay::storetriangle(Int a, Int b, Int c){
  Int triang=Delaunay::storetriangle(a, b, c);
  pointhash->store(hashfn.int64(a),triang);pointhash->store(hashfn.int64(b),triang);pointhash->store(hashfn.int64(c),triang);
  return triang;
};

void myDelaunay::erasetriangle(Int a, Int b, Int c, Int d0, Int d1, Int d2){
  Ullong key;
  Int j;
  key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
  if (trihash->get(key,j) == 0) throw("nonexistent triangle");
  pointhash->erase(hashfn.int64(a),j);
  pointhash->erase(hashfn.int64(b),j);
  pointhash->erase(hashfn.int64(c),j);
  Delaunay::erasetriangle(a, b, c, d0, d1, d2);};

myDelaunay::myDelaunay(vector<Point<2> > &pvec, Int options):Delaunay(pvec, options){
  pointhash = new Mhash<Ullong,Int,Nullhash>(6*npts,6*npts);
  for(int j=0;j<thelist.size();j++){
    if(thelist[j].stat==1){
      int a=thelist[j].p[0],b=thelist[j].p[1],c=thelist[j].p[2];
      pointhash->store(hashfn.int64(a),j);pointhash->store(hashfn.int64(b),j);pointhash->store(hashfn.int64(c),j);
    }
  }
}
