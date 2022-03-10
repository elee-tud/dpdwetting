#include "stress.hpp"
using namespace dpd;
void dpd::calculateStress(Particle* ptcl1, Particle* ptcl2, Real3D rij, Real3D fij){
    Rvec pairstress(9,0.0);
    pairstress[0]=fij[0]*rij[0];
    pairstress[1]=fij[1]*rij[1];
    pairstress[2]=fij[2]*rij[2];
    pairstress[3]=fij[0]*rij[1];
    pairstress[4]=fij[0]*rij[2];
    pairstress[5]=fij[1]*rij[0];
    pairstress[6]=fij[1]*rij[2];
    pairstress[7]=fij[2]*rij[0];
    pairstress[8]=fij[2]*rij[1];

    ptcl1->stress[0]+=pairstress[0];
    ptcl1->stress[1]+=pairstress[1];
    ptcl1->stress[2]+=pairstress[2];
    ptcl1->stress[3]+=pairstress[3];
    ptcl1->stress[4]+=pairstress[4];
    ptcl1->stress[5]+=pairstress[5];
    ptcl1->stress[6]+=pairstress[6];
    ptcl1->stress[7]+=pairstress[7];
    ptcl1->stress[8]+=pairstress[8];

    ptcl2->stress[0]+=pairstress[0];
    ptcl2->stress[1]+=pairstress[1];
    ptcl2->stress[2]+=pairstress[2];
    ptcl2->stress[3]+=pairstress[3];
    ptcl2->stress[4]+=pairstress[4];
    ptcl2->stress[5]+=pairstress[5];
    ptcl2->stress[6]+=pairstress[6];
    ptcl2->stress[7]+=pairstress[7];
    ptcl2->stress[8]+=pairstress[8];
    return;
}

void dpd::calculateSphericalStress(Particle* ptcl1, Particle* ptcl2, Real3D rij, Real3D fij, Real3D com, real** press, real dr){
    Real3D rad=(ptcl1->coord+ptcl2->coord)/2.-com;
    real r=rad.abs();
    real theta=acos(rad[2]/r);
    real phi=atan2(rad[1], rad[0]);
    real st=sin(theta);
    real ct=cos(theta);
    real sp=sin(phi);
    real cp=cos(phi);
    real stsp=st*sp;
    real stcp=st*cp;
    real ctsp=ct*sp;
    real ctcp=ct*cp;

    real r1=(ptcl1->coord-com).abs();
    real r2=(ptcl2->coord-com).abs();


    real T11=stcp;
    real T12=ctcp;
    real T13=-sp;
    real T21=stsp;
    real T22=ctsp;
    real T23=cp;
    real T31=ct;
    real T32=-st;

    real S11=fij[0]*rij[0];
    real S12=fij[0]*rij[1];
    real S13=fij[0]*rij[2];
    real S21=fij[1]*rij[0];
    real S22=fij[1]*rij[1];
    real S23=fij[1]*rij[2];
    real S31=fij[2]*rij[0];
    real S32=fij[2]*rij[1];
    real S33=fij[2]*rij[2];

    real srr=T11*(S11*T11+S12*T21+S13*T31)+T21*(S21*T11+S22*T21+S23*T31)+T31*(S31*T11+S32*T21+S33*T31);
    real stt=T12*(S11*T12+S12*T22+S13*T32)+T22*(S21*T12+S22*T22+S23*T32)+T32*(S31*T12+S32*T22+S33*T32);
    real spp=T13*(S11*T13+S12*T23)+T23*(S21*T13+S22*T23);


    int idx1=static_cast<int>(r1/dr);
    int idx2=static_cast<int>(r2/dr);
    


    if(idx2>idx1){
        real r12=r2-r1;
        for(int i=idx1;i<idx2;i++){
            press[0][i]+=srr/r12;
            press[1][i]+=stt/r12;
            press[2][i]+=spp/r12;
        }
    }
    else{
        real r12=r1-r2;
        for(int i=idx2;i<idx1;i++){
            press[0][i]+=srr/r12;
            press[1][i]+=stt/r12;
            press[2][i]+=spp/r12;
        }
    }






    return;
}

