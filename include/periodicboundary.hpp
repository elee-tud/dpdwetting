#ifndef __PERIODICBOUNDARY__HPP
#define __PERIODICBOUNDARY__HPP
#include "real3d.hpp"
#include "types.hpp"

namespace dpd{

class PeriodicBoundary{
private:
	Real3D box;
	Real3D half_box;
	Real3D neg_halfbox;
public:
	PeriodicBoundary();
	PeriodicBoundary(Real3D _box);
	~PeriodicBoundary();

	Real3D getBoxSize() const;
    void resetBoxSize(Real3D newbox);

    real getMinimumImageVectorComponentManyfold(const real& value, int i) const;
	Real3D getMinimumImageVectorManyfold(const Real3D& dist) const;
	Real3D getMinimumImageVectorManyfold(const Real3D& pos1, const Real3D& pos2) const;
	real getMinimumDistanceManyfold(const Real3D& dist);
	real getMinimumDistanceManyfold(const Real3D& pos1, const Real3D& pos2);

	real getMinimumImageVectorComponent(const real& value, int i) const;
	Real3D getMinimumImageVector(const Real3D& dist) const ;
	Real3D getMinimumImageVector(const Real3D& pos1, const Real3D& pos2) const;
	real getMinimumDistance(const Real3D& dist) const;
	real getMinimumDistance(const Real3D& pos1, const Real3D& pos2) const;
	Ivec isCrossingBox(const Real3D& pos1, const Real3D& pos2);
    real getVectorComponentIntoBox(const real& value, int i) const;
    Real3D getVectorIntoBox(const Real3D& vec) const;
    Real3D getBox(){ return box; }
    Real3D getHalfBox(){ return half_box; }
};
};
#endif
