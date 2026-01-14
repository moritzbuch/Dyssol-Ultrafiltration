/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#define DLL_EXPORT
#include "Ultrafiltration.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUltrafiltration();
}

//////////////////////////////////////////////////////////////////////////
/// Unit

void CUltrafiltration::CreateBasicInfo()
{
	/// Basic unit's info ///
	SetUnitName  ("Ultrafiltration");
	SetAuthorName("Buchholz");
	SetUniqueID  ("{85B875F4-0544-461C-AE52-197A1FDB4E7E}");
}

void CUltrafiltration::CreateStructure()
{
	/// Add ports ///
	port_Solvent = AddPort("Solvent", EUnitPort::INPUT);
	port_Permeat = AddPort("Permeat", EUnitPort::OUTPUT);

	/// Add unit parameters ///
	up_m_DF = AddConstRealParameter("m_DF", 1.0, "kg", "Mass during Diafiltration.", 0.0);
	up_mflow_P = AddConstRealParameter("mflow_P", 1.0, "kg/s", "Permeate mass flow.", 0.0);
	up_mflow_F = AddConstRealParameter("mflow_F", 1.0, "kg/s", "Recirculation mass flow.", 0.0);
	up_cmp_solv = AddCompoundParameter("cmp_solv", "Solvent compound.");
	up_cmp_j1 = AddCompoundParameter("cmp_R1", "Name of compound 1.");
	up_cmp_j2 = AddCompoundParameter("cmp_R2", "Name of compound 2.");
	up_cmp_j3 = AddCompoundParameter("cmp_R3", "Name of compound 3.");
	up_mode = AddComboParameter("Mode", EMode::RETFAC, { EMode::RETFAC, EMode::SOLUB }, { "Retention factor", "Solubility" }, "Mode of operation.");
	up_Rj1 = AddConstRealParameter("R1", 0.5, "kg", "Retention of compound 1.", 0.0, 1.0);
	up_Rj2 = AddConstRealParameter("R1", 0.5, "kg", "Retention of compound 2.", 0.0, 1.0);
	up_Rj3 = AddConstRealParameter("R1", 0.5, "kg", "Retention of compound 3.", 0.0, 1.0);
	up_c_solub_j1 = AddConstRealParameter("csolub_1", 0.1, "kg/kg_solv", "Solubility of compound 1.", 0.0);
	up_c_solub_j2 = AddConstRealParameter("csolub_2", 0.1, "kg/kg_solv", "Solubility of compound 2.", 0.0);
	up_c_solub_j3 = AddConstRealParameter("csolub_3", 0.1, "kg/kg_solv", "Solubility of compound 3.", 0.0);
	up_K_ls = AddConstRealParameter("K_ls", 1.0, "1/s", "Control parameter K for phase change rate between liquid and solid.", 0.0);
	up_K_sl = AddConstRealParameter("K_sl", 1.0, "1/s", "Control parameter K for phase change rate between liquid and solid.", 0.0);
	up_K_in = AddConstRealParameter("K_in", 1.0, "s", "Control parameter K for transition time to DF phase.", 0.0);


	AddParametersToGroup(up_mode, EMode::RETFAC, { up_Rj1, up_Rj2, up_Rj3 });
	AddParametersToGroup(up_mode, EMode::SOLUB, { up_c_solub_j1, up_c_solub_j2, up_c_solub_j3 });

	AddConstRealParameter("Rtol", 1-3, "-", "Relative tolerance.", 0.0, 1.0);
	AddConstRealParameter("Atol", 1-5, "-", "Relative tolerance.", 0.0, 1.0);


	/// Add holdups
	model.holdupTank = AddHoldup("Tank");
	model.streamRetentate = AddStream("Retentate");

	/// Set this unit as user data of model ///
	model.SetUserData(this);
	model.unit = static_cast<CUltrafiltration*>(model.GetUserData());
}

void CUltrafiltration::Initialize(double _time)
{
	/// Get flowsheet info
	model.M = GetCompoundsNumber();
	model.keys_cmp = GetAllCompounds();

	/// Get values of unit parameters
	model.par_m_K_DF = up_m_DF->GetValue();
	model.par_mflow_P = up_mflow_P->GetValue();
	model.par_mflow_F = up_mflow_F->GetValue();
	model.par_K_ls = up_K_ls->GetValue();
	model.par_K_sl = up_K_sl->GetValue();
	model.par_K_in = up_K_in->GetValue();

	model.mode = static_cast<EMode>(up_mode->GetValue());

	/// Get pointers to streams
	model.inletSolvent = port_Solvent->GetStream();
	model.outletPermeate = port_Permeat->GetStream();

	/// Compound information
	model.key_solv = up_cmp_solv->GetCompound();
	model.j_solv = GetCompoundIndex(model.key_solv);
	model.keys_R.push_back(up_cmp_j1->GetCompound());
	model.keys_R.push_back(up_cmp_j2->GetCompound());
	model.keys_R.push_back(up_cmp_j3->GetCompound());

	model.j_R.push_back(GetCompoundIndex(up_cmp_j1->GetCompound()));
	model.j_R.push_back(GetCompoundIndex(up_cmp_j2->GetCompound()));
	model.j_R.push_back(GetCompoundIndex(up_cmp_j3->GetCompound()));

	/// Retetion of components
	model.R_K_l.assign(model.M, 0);
	model.R_K_l[model.j_R[0]] = up_Rj1->GetValue();
	model.R_K_l[model.j_R[1]] = up_Rj2->GetValue();
	model.R_K_l[model.j_R[2]] = up_Rj3->GetValue();

	/// Solubility of components
	model.c_solub_K_l.assign(model.M, 1e+6);
	model.c_solub_K_l[model.j_R[0]] = up_c_solub_j1->GetValue();
	model.c_solub_K_l[model.j_R[1]] = up_c_solub_j2->GetValue();
	model.c_solub_K_l[model.j_R[2]] = up_c_solub_j3->GetValue();
	model.c_solub_K_l[model.j_solv] = 1; // Solvent is not limited by solubility.


	/// Additional information for plotting
	model.t_0_DF = -1;
	model.m_K_0_DF = 0;
	model.m_K_l_j_0.assign(model.M, 0);
	model.m_K_s_j_0.assign(model.M, 0);

	/// Get values of holdups
	// Masses of phases
	const double m_K_s_0 = model.SolMass(_time, model.holdupTank);      // Mass of solid in the holdup [kg].
	const double m_K_l_0 = model.LiqMass(_time, model.holdupTank);      // Mass of liquid in the holdup [kg].
	
	const double w_K_s_0 = m_K_s_0 / (m_K_s_0 + m_K_l_0); 		// Overall mass fraction of solid in the holdup [-].
	
	// compounds fractions
	const vector w_K_s_j_0 = model.SolFractions(_time, model.holdupTank);        // Mass fractions of each solid compound in holdup [-].
	const vector w_K_l_j_0 = model.LiqFractions(_time, model.holdupTank);        // Mass fractions of each liquid compound in holdup [-].


	if (model.mode == EMode::SOLUB)
	{
		// Initial compound masses in holdup before solubility check 
		const vector temp_m_K_l_j_0 = MultVector(w_K_l_j_0, m_K_l_0);
		const vector temp_m_K_s_j_0 = MultVector(w_K_s_j_0, m_K_s_0);

		const vector temp_m_K_j_0 = AddVectors(temp_m_K_l_j_0, temp_m_K_s_j_0);
		
		// Mass concentrations
		const vector temp_c_K_l_j_0 = MultVector(temp_m_K_j_0, 1. / temp_m_K_j_0[model.j_solv]); // Concentration of each liquid compound in holdup [kg/kg_solv].

		vector<double> c_K_l_j_0(model.M, 0);
		for (size_t j = 0; j < model.M; ++j)
		{
			if (j == model.j_solv) continue;
			
			c_K_l_j_0[j] = min(temp_c_K_l_j_0[j], model.c_solub_K_l[j]);
		}

		// Resulting compound masses in holdup
		model.m_K_l_j_0 = MultVector(c_K_l_j_0, temp_m_K_j_0[model.j_solv]);
		for (size_t j = 0; j < model.M; ++j)
		{
			model.m_K_s_j_0[j] = temp_c_K_l_j_0[j] * temp_m_K_j_0[model.j_solv] - model.m_K_l_j_0[j];
			if(model.m_K_s_j_0[j] < 0)
				RaiseError("Initial mass of solid component " + GetCompoundName(model.keys_cmp[j]) + " is negative. Check solubility of the compound.");
		}
	}
	else
	{
		// Compound masses in holdup
		model.m_K_l_j_0 = MultVector(w_K_l_j_0, m_K_l_0);
		model.m_K_s_j_0 = MultVector(w_K_s_j_0, m_K_s_0);
	}

	/// Initial values
	vector<double> mflow_P_s_j_0(model.M, 0); // Initial mass flow of solid components in permeate P [kg/s].
	vector<double> mflow_P_l_j_0(model.M, 0); // Initial mass flow of liquid components in permeate P [kg/s].
	
	/// Clear all state variables in model ///
	model.ClearVariables();

	/// Add state variables to the model ///
	model.i_m_K_s_j = model.AddDAEVariables(true, MultVector(w_K_s_j_0, m_K_s_0), 0, 0);
	model.i_m_K_l_j = model.AddDAEVariables(true, MultVector(w_K_l_j_0, m_K_l_0), 0, 0);
	model.i_mflow_K_ls_j = model.AddDAEVariables(false, vector<double>(model.M, 0), 0, 0);
	model.i_mflow_K_sl_j = model.AddDAEVariables(false, vector<double>(model.M, 0), 0, 0);
	model.i_mflow_P_s_j = model.AddDAEVariables(false, mflow_P_s_j_0, 0, 0);
	model.i_mflow_P_l_j = model.AddDAEVariables(false, mflow_P_l_j_0, 0, 0);
	model.i_mflow_in_l = model.AddDAEVariable(false, 0, 0, 0);
	//model.i_mflow_in_l = model.AddDAEVariable(true, 0, 0, 0);

	/// Set tolerances to the model ///
	//model.SetTolerance(1e-3, 1e-5);
	model.SetTolerance(GetConstRealParameterValue("Rtol"), GetConstRealParameterValue("Atol"));

	/// Set model to the solver ///
	if (!solver.SetModel(&model))
		RaiseError(solver.GetError());

	solver.SetMaxStep(10);

	/// Add Plots ///
	AddPlot("DF Factor", "Time", "DF Factor");
	AddCurveOnPlot("DF Factor", "DF Factor");

	AddPlot("Concentration", "Time", "Conc");
	for (size_t j = 0; j < model.M; ++j) AddCurveOnPlot("Concentration", "c_K_" + GetCompoundName(model.keys_cmp[j]));
	for (size_t j = 0; j < model.M; ++j) AddCurveOnPlot("Concentration", "c_solub_" + GetCompoundName(model.keys_cmp[j]));

	AddPlot("Total mass fraction", "Time", "Mfrac");
	for (size_t j = 0; j < model.M; ++j) AddCurveOnPlot("Total mass fraction", "w_" + GetCompoundName(model.keys_cmp[j]) + "_l");
	for (size_t j = 0; j < model.M; ++j) AddCurveOnPlot("Total mass fraction", "w_" + GetCompoundName(model.keys_cmp[j]) + "_s");

	AddPlot("Yield", "Time", "Yield");
	for (size_t j = 0; j < model.M; ++j)
	{
		if (j == model.j_solv) continue;
		AddCurveOnPlot("Yield", GetCompoundName(model.keys_cmp[j]));
	}
		

#if defined(_DEBUG) || defined(UNIT_DEBUG)
	model.iterations = 0;
#endif

}

void CUltrafiltration::Simulate(double _timeBeg, double _timeEnd)
{
	/// Run solver ///
	if (!solver.Calculate(_timeBeg, _timeEnd))
		RaiseError(solver.GetError());
}

void CUltrafiltration::SaveState()
{
	/// Save solver's state ///
	solver.SaveState();
}

void CUltrafiltration::LoadState()
{
	/// Load solver's state ///
	solver.LoadState();
}

void CUltrafiltration::Finalize()
{

}

//////////////////////////////////////////////////////////////////////////
/// Solver

void CUltrafiltrationDAEModel::Calculate(double _time, double* _vars)
{
	if (unit->HasError()) return;

	/// Clear global variables

	for (auto* v : all_to_zero) *v = 0.0;
	for (auto* v : all_to_M_vectors) v->assign(M, 0.0);

	const vector var_m_K_s_j = Slice(_vars, i_m_K_s_j);				// Currently calculated mass of solid components in tank K[kg].
	const vector var_m_K_l_j = Slice(_vars, i_m_K_l_j);				// Currently calculated mass of liquid components in tank K[kg].
	const vector var_mflow_K_ls_j = Slice(_vars, i_mflow_K_ls_j);	// Currently calculated mass flow of compenent going from liquid to solid state in tank K [kg/s].
	const vector var_mflow_K_sl_j = Slice(_vars, i_mflow_K_sl_j);	// Currently calculated mass flow of compenent going from solid to liquid state in tank K [kg/s].
	const vector var_mflow_P_s_j = Slice(_vars, i_mflow_P_s_j);		// Currently calculated mass flow of solid components in permeate P [kg/s].
	const vector var_mflow_P_l_j = Slice(_vars, i_mflow_P_l_j);		// Currently calculated mass flow of liquid components in permeate P [kg/s].
	double var_mflow_in_l = _vars[i_mflow_in_l];				// Currently calculated inlet mass flow of liquid to tank K [kg/s].

	// Incoming mass flow
	w_in_l_j = LiqFractions(_time, inletSolvent);

	m_K = VectorSum(var_m_K_s_j) + VectorSum(var_m_K_l_j);			// Total mass of components in tank K [kg].
	w_K_s = VectorSum(var_m_K_s_j) / m_K;							// Overall mass fraction of solid components in tank K [-].
	w_K_l = 1 - w_K_s;												// Overall mass fraction of liquid components in tank K [-].

	// Calculate mass fractions of solid and liquid components in tank K
	w_K_s_j = Normalized(var_m_K_s_j);
	w_K_l_j = Normalized(var_m_K_l_j);

	// Calculate mass concentrations of components in tank K
	c_K_l_j = MultVector(var_m_K_l_j, 1. / var_m_K_l_j[j_solv]);	// Concentration of each liquid compound in tank K [kg/kg_solv].

	c_K_l_j_test = MultVector(AddVectors(var_m_K_l_j, var_m_K_s_j), 1. / var_m_K_l_j[j_solv]);
	
	// Calculate mass flows of solid and liquid components in permeate P
	mflow_P_s = VectorSum(var_mflow_P_s_j);							// Total mass flow of solid components in permeate P [kg/s].
	mflow_P_l = VectorSum(var_mflow_P_l_j);							// Total mass flow of liquid components in permeate P [kg/s].
	mflow_P = mflow_P_s + mflow_P_l;										// Total mass flow in permeate P [kg/s].

	if (m_K <= par_m_K_DF)
	{
		if (t_0_DF == -1)
		{
			t_0_DF = _time;
			m_K_0_DF = m_K;
		}
		upd_mflow_in_l = (_time - t_0_DF) / par_K_in < 1 ? par_mflow_P * (_time - t_0_DF) / par_K_in : par_mflow_P;
		//upd_mflow_in_l = par_K_in * (par_m_K_DF - m_K);
	}
	else
		upd_mflow_in_l = 0;

	/// Flow to membrane
	vector<double> mflow_F_s_j(M, 0);
	vector<double> mflow_F_l_j(M, 0);
	for (size_t j = 0; j < M; ++j) mflow_F_s_j[j] = par_mflow_F * w_K_s * w_K_s_j[j];
	for (size_t j = 0; j < M; ++j) mflow_F_l_j[j] = par_mflow_F * w_K_l * w_K_l_j[j];

	/// Calculate current retentate flows
	for (size_t j = 0; j < M; ++j) mflow_R_s_j[j] = mflow_F_s_j[j] - var_mflow_P_s_j[j];
	for (size_t j = 0; j < M; ++j) mflow_R_l_j[j] = mflow_F_l_j[j] - var_mflow_P_l_j[j];

	/// Update DAE variables for solver
	// Mass balances
	for (size_t j = 0; j < M; ++j) dmdt_K_s_j[j] = -var_mflow_P_s_j[j] + var_mflow_K_ls_j[j] - var_mflow_K_sl_j[j];								// Mass balance for each solid component in tank K [kg/s].
	for (size_t j = 0; j < M; ++j) dmdt_K_l_j[j] = -var_mflow_P_l_j[j] - var_mflow_K_ls_j[j] + var_mflow_K_sl_j[j] + var_mflow_in_l * w_in_l_j[j];	// Mass balance for each liquid component in tank K [kg/s].

	
	/// Calculate permeation for solubility approach
	if (mode == EMode::SOLUB)
	{
		for (size_t j = 0; j < M; ++j) c_K_l_j[j] > c_solub_K_l[j] ? upd_mflow_K_ls_j[j] = par_K_ls * (c_K_l_j[j] - c_solub_K_l[j]) : 0;
		//for (size_t j = 0; j < M; ++j) (c_K_l_j[j] < (0.99 * c_solub_K_l[j])) && (var_m_K_s_j[j] > 1e-6) ? upd_mflow_K_sl_j[j] = par_K_sl * var_m_K_s_j[j]: 0;
		for (size_t j = 0; j < M; ++j) (c_K_l_j[j] < (0.9999999 * c_solub_K_l[j])) && (var_m_K_s_j[j] > 1e-9) ? upd_mflow_K_sl_j[j] = par_K_sl * (c_solub_K_l[j] - c_K_l_j[j]) * Clamp(pow(var_m_K_s_j[j] / 1e-5, 0.75), 0.0, 1.0) : 0;
		//for (size_t j = 0; j < M; ++j) (c_K_l_j[j] < (0.999 * c_solub_K_l[j])) && (var_m_K_s_j[j] < 1e-6 && var_m_K_s_j[j] > 0) ? upd_mflow_K_sl_j[j] = par_K_sl * (c_solub_K_l[j] - c_K_l_j[j]) * (1e-6 - var_m_K_s_j[j]) / 1e-6 * (1e-6 - var_m_K_s_j[j]) / 1e-6 : 0;
	}
	
	// Outgoing mass flows in permeate P
	switch (mode)
	{
	case EMode::RETFAC:
		for (size_t j = 0; j < M; ++j) upd_mflow_P_s_j[j] = 0;									// There is no permeation of solid material.
		for (size_t j = 0; j < M; ++j) upd_mflow_P_l_j[j] = mflow_F_l_j[j] * (1 - R_K_l[j]);	// R_K_l[j] is 0 for all j, i.e. all liquid components are permeating, except for defined compounds.
		break;
	case EMode::SOLUB:
		upd_mflow_P_s_j = vector<double>(M, 0);
		upd_mflow_P_l_j = MultVector(w_K_l_j, par_mflow_P);
		break;
	default:
		for (size_t j = 0; j < M; ++j) upd_mflow_P_s_j[j] = 0;
		for (size_t j = 0; j < M; ++j) upd_mflow_P_l_j[j] = 0;
		break;
	}

	unsigned debug = 0;
}

void CUltrafiltrationDAEModel::CalculateResiduals(double _time, double* _vars, double* _derivs, double* _res, void* _unit)
{
	Calculate(_time, _vars);

	/// Calculate residual ///
	for (size_t j = 0; j < M; ++j) _res[i_m_K_s_j[j]] = _derivs[i_m_K_s_j[j]] - dmdt_K_s_j[j];
	for (size_t j = 0; j < M; ++j) _res[i_m_K_l_j[j]] = _derivs[i_m_K_l_j[j]] - dmdt_K_l_j[j];
	//for (size_t j = 0; j < M; ++j) _res[i_mflow_K_ls_j[j]] = _derivs[i_mflow_K_ls_j[j]] - ddt_mflow_K_ls_j[j];
	for (size_t j = 0; j < M; ++j) _res[i_mflow_K_ls_j[j]] = _vars[i_mflow_K_ls_j[j]] - upd_mflow_K_ls_j[j];
	for (size_t j = 0; j < M; ++j) _res[i_mflow_K_sl_j[j]] = _vars[i_mflow_K_sl_j[j]] - upd_mflow_K_sl_j[j];
	for (size_t j = 0; j < M; ++j) _res[i_mflow_P_s_j[j]] = _vars[i_mflow_P_s_j[j]] - upd_mflow_P_s_j[j];
	for (size_t j = 0; j < M; ++j) _res[i_mflow_P_l_j[j]] = _vars[i_mflow_P_l_j[j]] - upd_mflow_P_l_j[j];
	//for (size_t j = 0; j < M; ++j) _res[i_mflow_in_l] = _derivs[i_mflow_in_l] - upd_mflow_in_l;
	for (size_t j = 0; j < M; ++j) _res[i_mflow_in_l] = _vars[i_mflow_in_l] - upd_mflow_in_l;

	#if defined(_DEBUG) || defined(UNIT_DEBUG)
		for (size_t i = 0; i < GetVariablesNumber(); ++i)
		{
			if (isnan(_res[i]))
			{
				unit->RaiseError("_res[" + to_string(i) + "] == nan");
			}
			if (isinf(_res[i]))
			{
				unit->RaiseError("_res[" + to_string(i) + "] == inf");
			}
		}
		iterations++;
	#endif
}

void CUltrafiltrationDAEModel::ResultsHandler(double _time, double* _vars, double* _derivs, void* _unit)
{
	Calculate(_time, _vars);

	unit->ShowInfo(to_string(_time) + "s...");
#if defined(_DEBUG) || defined(UNIT_DEBUG)
	unit->ShowInfo("Iterations: " + to_string(iterations));
	iterations = 0;
#endif

	const vector var_m_K_s_j = Slice(_vars, i_m_K_s_j);  // Currently calculated mass of solid components in tank K[kg].
	const vector var_m_K_l_j = Slice(_vars, i_m_K_l_j);  // Currently calculated mass of liquid components in tank K[kg].
	const vector var_mflow_P_s_j = Slice(_vars, i_mflow_P_s_j);  // Currently calculated mass flow of solid components in permeate P [kg/s].
	const vector var_mflow_P_l_j = Slice(_vars, i_mflow_P_l_j);  // Currently calculated mass flow of liquid components in permeate P [kg/s].

	holdupTank->SetCompoundsFractions(_time, EPhase::SOLID, Normalized(var_m_K_s_j));
	holdupTank->SetCompoundsFractions(_time, EPhase::LIQUID, Normalized(var_m_K_l_j));
	holdupTank->SetPhaseMass(_time, EPhase::SOLID, VectorSum(var_m_K_s_j));
	holdupTank->SetPhaseMass(_time, EPhase::LIQUID, VectorSum(var_m_K_l_j));
	holdupTank->SetPhaseMass(_time, EPhase::GAS, 0.0);
	holdupTank->SetMass(_time, VectorSum(var_m_K_s_j) + VectorSum(var_m_K_l_j));

	double mflow_P = VectorSum(var_mflow_P_s_j) + VectorSum(var_mflow_P_l_j);
	outletPermeate->CopyFromHoldup(_time, holdupTank, mflow_P);
	outletPermeate->SetPhaseMassFlow(_time, EPhase::SOLID, VectorSum(var_mflow_P_s_j));
	outletPermeate->SetPhaseMassFlow(_time, EPhase::LIQUID, VectorSum(var_mflow_P_l_j));
	outletPermeate->SetCompoundsFractions(_time, EPhase::SOLID, Normalized(var_mflow_P_s_j));
	outletPermeate->SetCompoundsFractions(_time, EPhase::LIQUID, Normalized(var_mflow_P_l_j));

	streamRetentate->CopyFromHoldup(_time, holdupTank, VectorSum(mflow_R_s_j) + VectorSum(mflow_R_s_j));
	streamRetentate->SetPhaseMassFlow(_time, EPhase::SOLID, VectorSum(mflow_R_s_j));
	streamRetentate->SetPhaseMassFlow(_time, EPhase::LIQUID, VectorSum(mflow_R_l_j));
	streamRetentate->SetCompoundsFractions(_time, EPhase::SOLID, Normalized(mflow_R_s_j));
	streamRetentate->SetCompoundsFractions(_time, EPhase::LIQUID, Normalized(mflow_R_l_j));

	/// Plotting ///
	if (_time > t_0_DF && t_0_DF > 0)
		unit->AddPointOnCurve("DF Factor", "DF Factor", _time, mflow_P * (_time - t_0_DF) / m_K_0_DF);

	for (size_t j = 0; j < M; ++j)
	{
		unit->AddPointOnCurve("Concentration", "c_K_" + unit->GetCompoundName(keys_cmp[j]), _time, c_K_l_j[j]);
		unit->AddPointOnCurve("Concentration", "c_solub_" + unit->GetCompoundName(keys_cmp[j]), _time, c_solub_K_l[j]);
	}

	for (size_t j = 0; j < M; ++j)
	{
		unit->AddPointOnCurve("Total mass fraction", "w_" + unit->GetCompoundName(keys_cmp[j]) + "_l", _time, w_K_l_j[j] * w_K_l);
		unit->AddPointOnCurve("Total mass fraction", "w_" + unit->GetCompoundName(keys_cmp[j]) + "_s", _time, w_K_s_j[j] * w_K_s);
	}

	for (size_t j = 0; j < M; ++j)
	{
		if (j == j_solv) continue;
		double yield = (m_K_l_j_0[j] + m_K_s_j_0[j]) > 0 ? ((m_K_l_j_0[j] + m_K_s_j_0[j]) - (var_m_K_l_j[j] + var_m_K_s_j[j])) / (m_K_l_j_0[j] + m_K_s_j_0[j]) : 0;
		unit->AddPointOnCurve("Yield", unit->GetCompoundName(keys_cmp[j]), _time, yield);
	}
}

std::vector<double> CUltrafiltrationDAEModel::SolFractions(double _time, const CBaseStream* _stream) const
{
	return _stream->GetCompoundsFractions(_time, EPhase::SOLID);
}

std::vector<double> CUltrafiltrationDAEModel::LiqFractions(double _time, const CBaseStream* _stream) const
{
	return _stream->GetCompoundsFractions(_time, EPhase::LIQUID);
}

double CUltrafiltrationDAEModel::SolMass(double _time, const CBaseStream* _stream) const
{
	return _stream->GetPhaseMass(_time, EPhase::SOLID);
}

double CUltrafiltrationDAEModel::LiqMass(double _time, const CBaseStream* _stream) const
{
	return _stream->GetPhaseMass(_time, EPhase::LIQUID);
}
