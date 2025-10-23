/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include "UnitDevelopmentDefines.h"

using namespace std;

class CUltrafiltration;


enum class EMode : size_t
{
	RETFAC, SOLUB
};

class CUltrafiltrationDAEModel : public CDAEModel
{
public:

	CUltrafiltration* unit{};            // Pointer to the dryer unit.


	vector<size_t> i_m_K_s_j{}; // Indices for masses of solid components in tank K [kg].
	vector<size_t> i_m_K_l_j{}; // Indices for masses of liquid components in tank K [kg].
	vector<size_t> i_mflow_K_ls_j{}; // Indices for mass flows of exchanges between solid and liquid phase [kg/s].
	vector<size_t> i_mflow_K_sl_j{}; // Indices for mass flows of exchanges between liquid and solid phase [kg/s].
	vector<size_t> i_mflow_P_s_j{};	// Indices for mass flows of solid components in permeate P [kg/s].
	vector<size_t> i_mflow_P_l_j{}; // Indices for mass flows of liquid components in permeate P [kg/s].
	size_t i_mflow_in_l{}; // Index for inlet mass flows of liquid to tank K [kg/s].

	vector<double> dmdt_K_s_j{}; // Derivative of Masses of solid components in tank K [kg].
	vector<double> dmdt_K_l_j{}; // Derivative of Masses of liquid components in tank K [kg].
	//vector<double> ddt_mflow_K_ls_j{}; // Derivative of phase change mass flows of between solid and liquid phase [kg/s].
	vector<double> upd_mflow_K_ls_j{}; // Derivative of phase change mass flows of between solid and liquid phase [kg/s].
	vector<double> upd_mflow_K_sl_j{}; // Derivative of phase change mass flows of between liquid and solid phase [kg/s].
	vector<double> upd_mflow_P_s_j{}; // Updated mass flows of solid components in permeate P [kg/s].
	vector<double> upd_mflow_P_l_j{}; // Updated mass flows of liquid components in permeate P [kg/s].
	double upd_mflow_in_l{}; // Updated inlet mass flows of liquid to tank K [kg/s].


	/// Inlet
	vector<double> w_in_l_j{}; // Mass fractions of liquid components in inlet during diafiltration [kg/kg].

	/// Tank
	double m_K{}; // Overall mass in tank K [kg].
	double w_K_s{}; // Overall mass fraction of solid components in tank K [-].
	double w_K_l{}; // Overall mass fraction of liquid components in tank K [-].
	vector<double> w_K_s_j{}; // Mass fraction of solid components in tank K [kg/kg].
	vector<double> w_K_l_j{}; // Mass fraction of liquid components in tank K [kg/kg].

	vector<double> c_K_l_j{}; // Concentration of liquid components in tank K [kg/kg_solv].
	vector<double> c_K_l_j_test{}; // Concentration of liquid components in tank K [kg/kg_solv].


	/// Retention stream
	vector<double> mflow_R_s_j{}; // Massflow of solid components in retentate line R [kg/s].
	vector<double> mflow_R_l_j{}; // Massflow of liquid components in retentate line R [kg/s].

	/// Permeate stream
	double mflow_P{}; // Overall massflow of permeate P [kg/s].
	double mflow_P_s{}; // Overall massflow of solid components in permeate P [kg/s].
	double mflow_P_l{}; // Overall massflow of liquid components in permeate P [kg/s].

	/// Compound specific
	string key_solv{}; // Key of solvent.
	size_t j_solv{}; // Index of solvent.
	
	vector<string> keys_R{}; // Keys of compounds for which retention/solubility is considered.
	vector<size_t> j_R{}; // Indices of compounds for which retention/solubility is considered.

	/// Additional information for plotting

	double t_0_DF{}; // Time point when diafiltration starts [s].
	double m_K_0_DF{}; // Total mass in tank K at the beginning of diafiltration [kg].
	vector<double> m_K_l_j_0{}; // Initial mass of solid components in tank K [kg].
	vector<double> m_K_s_j_0{}; // Initial mass of liquid components in tank K [kg].

	/// Separation specific
	// Mode Retention factor
	vector<double> R_K_l{}; // Retention of liquid components in tank K [-].
	
	// Mode Solubility
	vector<double> c_solub_K_l{}; // Solubility of liquid components in tank K [kg/kg_solv].

	/// User parameters
	double par_m_K_DF{}; // Mass during Diafiltration.
	double par_mflow_P{}; // Permeate mass flow.
	double par_mflow_F{}; // Recirculation mass flow.
	double par_K_ls{}; // Control parameter K for phase change rate between liquid and solid.
	double par_K_sl{}; // Control parameter K for phase change rate between solid and liquid.
	double par_K_in{}; // Control parameter K for transition time to DF phase.
	EMode mode{}; // Mode of operation.

	/// Flowsheet info
	size_t M{};                  // Number of compounds [#].
	vector<string> keys_cmp;    // Keys of all defined compounds.
	vector<size_t> j_cmp;		// Indices of all vapor compounds.

	/// Pointers to streams and holdups
	CHoldup* holdupTank{};          // Tank K holdup.
	CStream* inletSolvent{};		// Solvent inlet.
	CStream* streamRetentate{};		// Retentate stream.
	CStream* outletPermeate{};		// Permeate stream.

	// Vector storage for setting to zero
	// Pointers to all variables that should be set to zero.
	vector<double*> all_to_zero{  }; 

	// These vectors will be resized to the number of components and the elements set to 0.
	vector<vector<double>*> all_to_M_vectors{
		&w_in_l_j,& w_K_s_j,& w_K_l_j,& mflow_R_s_j,& mflow_R_l_j,& dmdt_K_s_j,& dmdt_K_l_j,& upd_mflow_K_ls_j, & upd_mflow_K_sl_j ,& upd_mflow_P_s_j,& upd_mflow_P_l_j
	}; 


#if defined(_DEBUG) || defined(UNIT_DEBUG)
	size_t iterations{};                   // Number of internal iterations needed for convergence of a current time point [#].
#endif

public:
	/**
	 * \brief Do actual calculations for values and the derivatives of the state variables.
	 * \param _time Time point.
	 * \param _vars State variables - guesses from the solver.
	 */
	void Calculate(double _time, double* _vars);

	/**
	 * \brief Calculates values and the derivatives of the state variables.
	 * \param _time Time point.
	 * \param _vars State variables - guesses from the solver.
	 * \param _ders Derivatives of the state variables - guesses from the solver.
	 * \param _res Residuals of the state variables - calculated by the model.
	 * \param _unit Pointer to the unit.
	 */
	void CalculateResiduals(double _time, double* _vars, double* _derivs, double* _res, void* _unit) override;
	void ResultsHandler(double _time, double* _vars, double* _derivs, void* _unit) override;

	/**
	 * Returns mass fractions of solid compounds from the stream.
	 */
	vector<double> SolFractions(double _time, const CBaseStream* _stream) const;
	/**
	 * Returns mass fractions of liquid compounds from the stream.
	 */
	vector<double> LiqFractions(double _time, const CBaseStream* _stream) const;
	/**
	* Returns mass fractions of gas compounds from the stream without vapor compounds.
	*/
	inline double SolMass(double _time, const CBaseStream* _stream) const;
	/**
	 * Returns mass/mass flow of the liquid compounds from the stream.
	 */
	inline double LiqMass(double _time, const CBaseStream* _stream) const;
	/**
	* Returns mass/mass flow of the gas compounds without vapor compounds from the stream.
	*/
};

class CUltrafiltration : public CDynamicUnit
{
private:
	CUltrafiltrationDAEModel model{};		// Model of DAE
	CDAESolver solver{};		// Solver of DAE

public:

	CConstRealUnitParameter* up_m_DF{};			// Mass during Diafiltration.
	CConstRealUnitParameter* up_mflow_P{};		// Permeate Mass flow.
	CConstRealUnitParameter* up_mflow_F{};		// Recirculation mass flow.
	CComboUnitParameter* up_mode{};				// Mode of operation.
	CCompoundUnitParameter* up_cmp_solv{};		// Solvent compound.
	CCompoundUnitParameter* up_cmp_j1{};		// Compound 1 which is reduced during filtration.
	CCompoundUnitParameter* up_cmp_j2{};		// Compound 2 which is reduced during filtration.
	CCompoundUnitParameter* up_cmp_j3{};		// Compound 3 which is reduced during filtration.

	CConstRealUnitParameter* up_Rj1{};			// Retention of compound 1.
	CConstRealUnitParameter* up_Rj2{};			// Retention of compound 2.
	CConstRealUnitParameter* up_Rj3{};			// Retention of compound 3.

	CConstRealUnitParameter* up_c_solub_j1{};	// Solubility of compound 1.
	CConstRealUnitParameter* up_c_solub_j2{};	// Solubility of compound 2.
	CConstRealUnitParameter* up_c_solub_j3{};	// Solubility of compound 3.
	CConstRealUnitParameter* up_K_ls{}; 		// Control parameter K for phase change rate between liquid and solid.
	CConstRealUnitParameter* up_K_sl{}; 		// Control parameter K for phase change rate between solid and liquid.
	CConstRealUnitParameter* up_K_in{}; 		// Control parameter K for phase change for ramping up liquid inlet flow.



	CUnitPort* port_Solvent{};
	CUnitPort* port_Permeat{};

public:
	void CreateBasicInfo() override;
	void CreateStructure() override;
	void Initialize(double _time) override;
	void Simulate(double _timeBeg, double _timeEnd) override;
	void SaveState() override;
	void LoadState() override;
	void Finalize() override;

	
};
