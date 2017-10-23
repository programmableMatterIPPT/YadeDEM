/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@imag.fr>         *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/*************************************************************************
*  Modified version below of the original                                *
*  CohesiveFrictionalContactLaw.cpp:                                     *
*                                                                        *
*  Copyright (C) 2017 by Jakub Lengiewicz <jleng@ippt.pan.pl>            *
*  Copyright (C) 2017 by Pawel Chodkiewicz <pawel@chodkiewicz.com.pl>    * 
*  Copyright (C) 2017 by Pawel Holobut <pholob@ippt.pan.pl>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "pmCohesiveFrictionalContactLaw.hpp"
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>
#include<core/Scene.hpp>

YADE_PLUGIN((PmCohesiveFrictionalContactLaw)(PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment)(PmCohFrictMat)(PmCohFrictPhys)(PmIp2_CohFrictMat_CohFrictMat_CohFrictPhys));
CREATE_LOGGER(PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment);

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::getPlasticDissipation() {return (Real) plasticDissipation;}
void PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::initPlasticDissipation(Real initVal) {plasticDissipation.reset(); plasticDissipation+=initVal;}

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::normElastEnergy()
{
	Real normEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*>(I->phys.get());
		if (phys) {
			normEnergy += 0.5*(phys->normalForce.squaredNorm()/phys->kn);
		}
	}
	return normEnergy;
}

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::shearElastEnergy()
{
	Real shearEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*>(I->phys.get());
		if (phys) {
			shearEnergy += 0.5*(phys->shearForce.squaredNorm()/phys->ks);
		}
	}
	return shearEnergy;
}

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::bendingElastEnergy()
{
	Real bendingEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*>(I->phys.get());
		if (phys) {
			bendingEnergy += 0.5*(phys->moment_bending.squaredNorm()/phys->kr);
		}
	}
	return bendingEnergy;
}

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::twistElastEnergy()
{
	Real twistEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*>(I->phys.get());
		if (phys) {
			twistEnergy += 0.5*(phys->moment_twist.squaredNorm()/phys->ktw);
		}
	}
	return twistEnergy;
}

Real PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::totalElastEnergy()
{
	Real totalEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*>(I->phys.get());
		if (phys) {
			totalEnergy += 0.5*(phys->normalForce.squaredNorm()/phys->kn);
			totalEnergy += 0.5*(phys->shearForce.squaredNorm()/phys->ks);
			totalEnergy += 0.5*(phys->moment_bending.squaredNorm()/phys->kr);
			totalEnergy += 0.5*(phys->moment_twist.squaredNorm()/phys->ktw);
		}
	}
	return totalEnergy;
}

void PmCohesiveFrictionalContactLaw::action()
{
	if(!functor) functor=shared_ptr<PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment>(new PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment);
	functor->always_use_moment_law = always_use_moment_law;
	functor->shear_creep=shear_creep;
	functor->twist_creep=twist_creep;
	functor->creep_viscosity=creep_viscosity;
	functor->scene=scene;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		functor->go(I->geom, I->phys, I.get());
	}
}

bool PmLaw2_ScGeom6D_CohFrictPhys_CohesionMoment::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real& dt = scene->dt;
	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGeom6D* geom  = YADE_CAST<ScGeom6D*> (ig.get());
	PmCohFrictPhys* phys = YADE_CAST<PmCohFrictPhys*> (ip.get());
	Vector3r& shearForce    = phys->shearForce;
	//added
	State* de1 = Body::byId(id1,scene)->state.get();
	State* de2 = Body::byId(id2,scene)->state.get();	

	const Vector3r c1x = (geom->contactPoint - de1->pos);
	const Vector3r c2x = (geom->contactPoint - de2->pos);
	
	const Vector3r relativeVelocity = (de2->vel+(de2->angVel).cross(c2x)) - (de1->vel+(de1->angVel).cross(c1x));
	const Real normalVelocity	= (geom->normal).dot(relativeVelocity);
	const Vector3r shearVelocity	= relativeVelocity - normalVelocity*geom->normal;
	//added

	if (contact->isFresh(scene)) shearForce   = Vector3r::Zero();
	Real un     = geom->penetrationDepth;
	//modified
	//Real Fn    = phys->kn*(un-phys->unp);
	Real Fn    = phys->kn*(un-phys->unp)-phys->cn*normalVelocity;
	//modified

	if (phys->fragile && (-Fn)> phys->normalAdhesion) {
		// BREAK due to tension
		return false;
	} else {
		if ((-Fn)> phys->normalAdhesion) {//normal plasticity
			Fn=-phys->normalAdhesion;
			phys->unp = un+phys->normalAdhesion/phys->kn;
			if (phys->unpMax && phys->unp<phys->unpMax)
				return false;
		}
		phys->normalForce = Fn*geom->normal;
		//removed
		//State* de1 = Body::byId(id1,scene)->state.get();
		//State* de2 = Body::byId(id2,scene)->state.get();
		//removed
		///////////////////////// CREEP START ///////////
		if (shear_creep) shearForce -= phys->ks*(shearForce*dt/creep_viscosity);
		///////////////////////// CREEP END ////////////

		Vector3r& shearForce = geom->rotate(phys->shearForce);
		const Vector3r& dus = geom->shearIncrement();

		//Linear elasticity giving "trial" shear force
		shearForce -= phys->ks*dus;

		//Real Fs = phys->shearForce.norm();//removed
		Vector3r totalShearForce = phys->shearForce-phys->cs*shearVelocity;//added
		Real Fs = totalShearForce.norm();//added
		Real maxFs = phys->shearAdhesion;
		//if (!phys->cohesionDisablesFriction || maxFs==0)//removed
		if (!phys->cohesionDisablesFriction && Fn>0)//added
			maxFs += Fn*phys->tangensOfFrictionAngle;
		maxFs = std::max((Real) 0, maxFs);
		if (Fs  > maxFs) {//Plasticity condition on shear force
			if (phys->fragile && !phys->cohesionBroken) {
				if (shearBreaksCohesion) {//added
					phys->SetBreakingState();
					//maxFs = max((Real) 0, Fn*phys->tangensOfFrictionAngle);//removed
				}//added
				else shearForce = Vector3r::Zero();//added
			}
			maxFs = maxFs / Fs;
			//Vector3r trialForce = shearForce;//removed
			Vector3r trialForce = totalShearForce;//added
			totalShearForce *= maxFs;//added
			shearForce *= maxFs;
			if (scene->trackEnergy){
				Real sheardissip=((1/phys->ks)*(trialForce-shearForce))/*plastic disp*/ .dot(shearForce)/*active force*/;
				if(sheardissip>0) scene->energy->add(sheardissip,"shearDissip",shearDissipIx,/*reset*/false);}
			//if (Fn<0)  phys->normalForce = Vector3r::Zero();//removed
		}
		//Apply the force
		//applyForceAtContactPoint(-phys->normalForce-shearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist): Vector3r::Zero()));//removed
		applyForceAtContactPoint(-phys->normalForce-totalShearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist): Vector3r::Zero()));//added
//removed
// 		Vector3r force = -phys->normalForce-shearForce;
// 		scene->forces.addForce(id1,force);
// 		scene->forces.addForce(id2,-force);
// 		scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force));
// 		scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force));
//removed
		/// Moment law  ///
		if (phys->momentRotationLaw && (!phys->cohesionBroken || always_use_moment_law)) {
			if (!useIncrementalForm){
				if (twist_creep) {
					Real viscosity_twist = creep_viscosity * std::pow((2 * std::min(geom->radius1,geom->radius2)),2) / 16.0;
					Real angle_twist_creeped = geom->getTwist() * (1 - dt/viscosity_twist);
					Quaternionr q_twist(AngleAxisr(geom->getTwist(),geom->normal));
					Quaternionr q_twist_creeped(AngleAxisr(angle_twist_creeped,geom->normal));
					Quaternionr q_twist_delta(q_twist_creeped * q_twist.conjugate());
					geom->twistCreep = geom->twistCreep * q_twist_delta;
				}
				//added
				Vector3r relAngVel = de2->angVel-de1->angVel;
				Vector3r relAngVelTwist = geom->normal.dot(relAngVel)*geom->normal;
				Vector3r relAngVelBend = relAngVel - relAngVelTwist;
				//added
				//modified
				//phys->moment_twist = (geom->getTwist()*phys->ktw)*geom->normal;
				//phys->moment_bending = geom->getBending() * phys->kr;
				phys->moment_twist = (geom->getTwist()*phys->ktw)*geom->normal - relAngVelTwist*phys->ctw;
				phys->moment_bending = geom->getBending()*phys->kr - relAngVelBend*phys->cr;
				//modified
			}	
			else{ // Use incremental formulation to compute moment_twis and moment_bending (no twist_creep is applied)
				if (twist_creep) throw std::invalid_argument("Law2_ScGeom6D_CohFrictPhys_CohesionMoment: no twis creep is included if the incremental form for the rotations is used.");
				Vector3r relAngVel = geom->getRelAngVel(de1,de2,dt);
				// *** Bending ***//
				Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel)*geom->normal; // keep only the bending part
				Vector3r relRotBend = relAngVelBend*dt; // relative rotation due to rolling behaviour	
				// incremental formulation for the bending moment (as for the shear part)
				Vector3r& momentBend = phys->moment_bending;
				momentBend = geom->rotate(momentBend); // rotate moment vector (updated)
				momentBend = momentBend-phys->kr*relRotBend;
				// ----------------------------------------------------------------------------------------
				// *** Torsion ***//
				Vector3r relAngVelTwist = geom->normal.dot(relAngVel)*geom->normal;
				Vector3r relRotTwist = relAngVelTwist*dt; // component of relative rotation along n  FIXME: sign?
				// incremental formulation for the torsional moment
				Vector3r& momentTwist = phys->moment_twist;
				momentTwist = geom->rotate(momentTwist); // rotate moment vector (updated)
				momentTwist = momentTwist-phys->ktw*relRotTwist; // FIXME: sign?
			}
			/// Plasticity ///
			// limit rolling moment to the plastic value, if required
			Real RollMax = phys->maxRollPl*phys->normalForce.norm();
			if (RollMax>0.){ // do we want to apply plasticity?
				if (!useIncrementalForm) LOG_WARN("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not be applied correctly (the total formulation would not reproduce irreversibility).");
				Real scalarRoll = phys->moment_bending.norm();
				if (scalarRoll>RollMax){ // fix maximum rolling moment
					Real ratio = RollMax/scalarRoll;
					phys->moment_bending *= ratio;
				}	
			}
			// limit twisting moment to the plastic value, if required
			Real TwistMax = phys->maxTwistMoment.norm();
			if (TwistMax>0.){ // do we want to apply plasticity?
				if (!useIncrementalForm) LOG_WARN("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not be applied correctly (the total formulation would not reproduce irreversibility).");
				Real scalarTwist= phys->moment_twist.norm();
				if (scalarTwist>TwistMax){ // fix maximum rolling moment
					Real ratio = TwistMax/scalarTwist;
					phys->moment_twist *= ratio;
				}	
			}
			// Apply moments now
			Vector3r moment = phys->moment_twist + phys->moment_bending;
			scene->forces.addTorque(id1,-moment);
			scene->forces.addTorque(id2, moment);			
		}
		/// Moment law END       ///
	}
	return true;
}

void PmIp2_CohFrictMat_CohFrictMat_CohFrictPhys::go(const shared_ptr<Material>& b1    // PmCohFrictMat
                                        , const shared_ptr<Material>& b2 // PmCohFrictMat
                                        , const shared_ptr<Interaction>& interaction)
{
	PmCohFrictMat* sdec1 = static_cast<PmCohFrictMat*>(b1.get());
	PmCohFrictMat* sdec2 = static_cast<PmCohFrictMat*>(b2.get());
	ScGeom6D* geom = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	//modified
	const shared_ptr<Body>& bo1 = Body::byId(interaction->getId1(),scene);
	const shared_ptr<Body>& bo2 = Body::byId(interaction->getId2(),scene);

	PmState* st1 = static_cast<PmState*>(bo1->state.get());
	PmState* st2 = static_cast<PmState*>(bo2->state.get());

	int minIndex1 = st1->minIndex;
	int minIndex2 = st2->minIndex;
	int index1 = st1->index;
	int index2 = st2->index;
	int indexn1 = st1->indexn;
	int indexn2 = st2->indexn;
	bool bodiesCohere = ((minIndex1 == 0 || minIndex1 == index2 || minIndex1 == indexn2) && (minIndex2 == 0 || minIndex2 == index1 || minIndex2 == indexn1)) ? true : false;
	//modified

	//Create cohesive interractions only once
	if (setCohesionNow && cohesionDefinitionIteration==-1) cohesionDefinitionIteration=scene->iter;
	if (setCohesionNow && cohesionDefinitionIteration!=-1 && cohesionDefinitionIteration!=scene->iter) {
		cohesionDefinitionIteration = -1;
		setCohesionNow = 0;}

	if (geom) {
		if (!interaction->phys) {
			interaction->phys = shared_ptr<PmCohFrictPhys>(new PmCohFrictPhys());
			PmCohFrictPhys* contactPhysics = YADE_CAST<PmCohFrictPhys*>(interaction->phys.get());
			Real Ea 	= sdec1->young;
			Real Eb 	= sdec2->young;
			Real Va 	= sdec1->poisson;
			Real Vb 	= sdec2->poisson;
			Real Da 	= geom->radius1;
			Real Db 	= geom->radius2;
			Real fa 	= sdec1->frictionAngle;
			Real fb 	= sdec2->frictionAngle;
			Real Kn = 2.0*Ea*Da*Eb*Db/(Ea*Da+Eb*Db);//harmonic average of two stiffnesses
			//added
			Real Cna = sdec1->cn;
			Real Cnb = sdec2->cn;
			Real Csa = sdec1->cs;
			Real Csb = sdec2->cs;

			if (Cna && Cnb){ 
				Real CnaDa = Cna*Da; Real CnbDb = Cnb*Db;
				contactPhysics->cn = 2.0*CnaDa*CnbDb/(CnaDa+CnbDb);
				contactPhysics->cr = 2.0*CnaDa*Da*CnbDb*Db/(CnaDa*Da+CnbDb*Db);
			}
			else{ contactPhysics->cn = 0; contactPhysics->cr = 0; }

			if (Csa && Csb){
				Real CsaDa = Csa*Da; Real CsbDb = Csb*Db;
				contactPhysics->cs = 2.0*CsaDa*CsbDb/(CsaDa+CsbDb);
				contactPhysics->ctw = 2.0*CsaDa*Da*CsbDb*Db/(CsaDa*Da+CsbDb*Db);
			}
			else{ contactPhysics->cs = 0; contactPhysics->ctw = 0; }
			//added

			// harmonic average of alphas parameters
			Real AlphaKr = 2.0*sdec1->alphaKr*sdec2->alphaKr/(sdec1->alphaKr+sdec2->alphaKr);
			Real AlphaKtw;
			if (sdec1->alphaKtw && sdec2->alphaKtw) AlphaKtw = 2.0*sdec1->alphaKtw*sdec2->alphaKtw/(sdec1->alphaKtw+sdec2->alphaKtw);
			else AlphaKtw=0;

			Real Ks;
			if (Va && Vb) Ks = 2.0*Ea*Da*Va*Eb*Db*Vb/(Ea*Da*Va+Eb*Db*Vb);//harmonic average of two stiffnesses with ks=V*kn for each sphere
			else Ks=0;

			// Jean-Patrick Plassiard, Noura Belhaine, Frederic
			// Victor Donze, "Calibration procedure for spherical
			// discrete elements using a local moemnt law".
			contactPhysics->kr = Da*Db*Ks*AlphaKr;
			contactPhysics->ktw = Da*Db*Ks*AlphaKtw;
			contactPhysics->tangensOfFrictionAngle		= std::tan(std::min(fa,fb));

			if ((setCohesionOnNewContacts || setCohesionNow) && sdec1->isCohesive && sdec2->isCohesive)
			{
				contactPhysics->cohesionBroken = false;
				//modified
				//contactPhysics->normalAdhesion = std::min(sdec1->normalCohesion,sdec2->normalCohesion)*pow(std::min(Db, Da),2);
				//contactPhysics->shearAdhesion = std::min(sdec1->shearCohesion,sdec2->shearCohesion)*pow(std::min(Db, Da),2);
				Real normalCohesion1 [5] = {sdec1->normalCohesion,sdec1->normalCohesion1,sdec1->normalCohesion2,sdec1->normalCohesion3,sdec1->normalCohesion4};
				Real shearCohesion1 [5] = {sdec1->shearCohesion,sdec1->shearCohesion1,sdec1->shearCohesion2,sdec1->shearCohesion3,sdec1->shearCohesion4};
				Real normalCohesion2 [5] = {sdec2->normalCohesion,sdec2->normalCohesion1,sdec2->normalCohesion2,sdec2->normalCohesion3,sdec2->normalCohesion4};
				Real shearCohesion2 [5] = {sdec2->shearCohesion,sdec2->shearCohesion1,sdec2->shearCohesion2,sdec2->shearCohesion3,sdec2->shearCohesion4};
				contactPhysics->normalAdhesion = std::min(normalCohesion1[sdec2->num],normalCohesion2[sdec1->num])*pow(std::min(Db, Da),2);
				contactPhysics->shearAdhesion = std::min(shearCohesion1[sdec2->num],shearCohesion2[sdec1->num])*pow(std::min(Db, Da),2);
				//modified
				geom->initRotations(*(Body::byId(interaction->getId1(),scene)->state),*(Body::byId(interaction->getId2(),scene)->state));
			}
			contactPhysics->kn = Kn;
			contactPhysics->ks = Ks;

			contactPhysics->maxRollPl = min(sdec1->etaRoll*Da,sdec2->etaRoll*Db);
			contactPhysics->momentRotationLaw=(sdec1->momentRotationLaw && sdec2->momentRotationLaw);
			//contactPhysics->elasticRollingLimit = elasticRollingLimit;
		}
		else {// !isNew, but if setCohesionNow, all contacts are initialized like if they were newly created
			PmCohFrictPhys* contactPhysics = YADE_CAST<PmCohFrictPhys*>(interaction->phys.get());
			if ((setCohesionNow && sdec1->isCohesive && sdec2->isCohesive) || contactPhysics->initCohesion)
			{
				contactPhysics->cohesionBroken = false;
				//modified
				//contactPhysics->normalAdhesion = std::min(sdec1->normalCohesion,sdec2->normalCohesion)*pow(std::min(geom->radius2, geom->radius1),2);
				//contactPhysics->shearAdhesion = std::min(sdec1->shearCohesion,sdec2->shearCohesion)*pow(std::min(geom->radius2, geom->radius1),2);
				Real normalCohesion1 [5] = {sdec1->normalCohesion,sdec1->normalCohesion1,sdec1->normalCohesion2,sdec1->normalCohesion3,sdec1->normalCohesion4};
				Real shearCohesion1 [5] = {sdec1->shearCohesion,sdec1->shearCohesion1,sdec1->shearCohesion2,sdec1->shearCohesion3,sdec1->shearCohesion4};
				Real normalCohesion2 [5] = {sdec2->normalCohesion,sdec2->normalCohesion1,sdec2->normalCohesion2,sdec2->normalCohesion3,sdec2->normalCohesion4};
				Real shearCohesion2 [5] = {sdec2->shearCohesion,sdec2->shearCohesion1,sdec2->shearCohesion2,sdec2->shearCohesion3,sdec2->shearCohesion4};
				contactPhysics->normalAdhesion = std::min(normalCohesion1[sdec2->num],normalCohesion2[sdec1->num])*pow(std::min(geom->radius2, geom->radius1),2);
				contactPhysics->shearAdhesion = std::min(shearCohesion1[sdec2->num],shearCohesion2[sdec1->num])*pow(std::min(geom->radius2, geom->radius1),2);
				//modified
				geom->initRotations(*(Body::byId(interaction->getId1(),scene)->state),*(Body::byId(interaction->getId2(),scene)->state));
				contactPhysics->initCohesion=false;
			}
		}
		//added
		if(!bodiesCohere){
			PmCohFrictPhys* contactPhysics = YADE_CAST<PmCohFrictPhys*>(interaction->phys.get());
			contactPhysics->normalAdhesion = 0;
			contactPhysics->shearAdhesion = 0;
		}
		//added
	}
};

