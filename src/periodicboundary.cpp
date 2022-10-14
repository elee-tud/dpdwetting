#include "periodicboundary.hpp"
using namespace dpd;
PeriodicBoundary::PeriodicBoundary(){}
PeriodicBoundary::PeriodicBoundary(Real3D _box){
	box=_box;
	half_box=box/2.;
	neg_halfbox=-half_box;
}
PeriodicBoundary::~PeriodicBoundary(){}

Real3D PeriodicBoundary::getBoxSize() const{ return box; }

real PeriodicBoundary::getMinimumImageVectorComponentManyfold(const real& value, int i) const{
	real out=value;
	if(out>=half_box[i]){
		while(out>=half_box[i])
			out-=box[i];
	}
	else if(out<neg_halfbox[i]){
		while(out<neg_halfbox[i]){
			out+=box[i];
		}
	}
	return out;
}

Real3D PeriodicBoundary::getMinimumImageVectorManyfold(const Real3D& dist) const{
	return Real3D(getMinimumImageVectorComponentManyfold(dist[0], 0), 
	              getMinimumImageVectorComponentManyfold(dist[1], 1), 
	              getMinimumImageVectorComponentManyfold(dist[2], 2));
}

Real3D PeriodicBoundary::getMinimumImageVectorManyfold(const Real3D& pos1, const Real3D& pos2) const{
	return getMinimumImageVectorManyfold(pos1-pos2);
}

real PeriodicBoundary::getMinimumDistanceManyfold(const Real3D& dist){
	return getMinimumImageVectorManyfold(dist).abs();
}

real PeriodicBoundary::getMinimumDistanceManyfold(const Real3D& pos1, const Real3D& pos2){
	return getMinimumImageVectorManyfold(pos1-pos2).abs();
}

real PeriodicBoundary::getMinimumImageVectorComponent(const real& value, int i) const{
	real out=value;
	if(out>=half_box[i])
		out-=box[i];
	else if(out<neg_halfbox[i])
		out+=box[i];
    if(out==0.0)
        out=0.00001;
	return out;
}


Real3D PeriodicBoundary::getMinimumImageVector(const Real3D& dist) const {
	return Real3D(getMinimumImageVectorComponent(dist[0], 0),
	              getMinimumImageVectorComponent(dist[1], 1),
	              getMinimumImageVectorComponent(dist[2], 2));
}




Real3D PeriodicBoundary::getMinimumImageVector(const Real3D& pos1, const Real3D& pos2) const{
	return getMinimumImageVector(pos1-pos2);
}

real PeriodicBoundary::getMinimumDistance(const Real3D& dist) const {
	return getMinimumImageVector(dist).abs();
}
real PeriodicBoundary::getMinimumDistance(const Real3D& pos1, const Real3D& pos2) const{
	return getMinimumImageVector(pos1, pos2).abs();
}

Ivec PeriodicBoundary::isCrossingBox(const Real3D& pos1, const Real3D& pos2){
	Ivec cross(3, 0);
	Real3D dist=pos1-pos2;
	if(dist[0]>=half_box[0])
		cross[0]=-1;
	else if(dist[0]<neg_halfbox[0])
		cross[0]=1;
	if(dist[1]>=half_box[1])
		cross[1]=-1;
	else if(dist[1]<neg_halfbox[1])
		cross[1]=1;
	if(dist[2]>=half_box[2])
		cross[2]=-1;
	else if(dist[2]<neg_halfbox[2])
		cross[2]=1;
	return cross;
}

real PeriodicBoundary::getVectorComponentIntoBox(const real& value, int i) const{
    real out=value;
    if(out>=box[i]){
        while(out>=box[i]){
            out-=box[i];
        }
    }
    else if(out<0.0){
        while(out<0.0){
            out+=box[i];
        }
    }
    /*Here, the particle is brought into the box once more.*/
    if(out>=box[i])
        out-=box[i];
    else if(out<0.0)
        out+=box[i];
    return out;
}

Real3D PeriodicBoundary::getVectorIntoBox(const Real3D& vec) const{
    return Real3D(getVectorComponentIntoBox(vec[0], 0), getVectorComponentIntoBox(vec[1], 1), getVectorComponentIntoBox(vec[2], 2));
}





void PeriodicBoundary::resetBoxSize(Real3D newbox){ 
    box=newbox; 
	half_box=box/2.;
	neg_halfbox=-half_box;
    return;
}
