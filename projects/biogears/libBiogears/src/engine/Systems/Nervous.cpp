/**************************************************************************************
Copyright 2015 Applied Research Associates, Inc.
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the License
at:
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
**************************************************************************************/
#include <biogears/engine/Systems/Nervous.h>

#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/patient/actions/SEPupillaryResponse.h>
#include <biogears/cdm/properties/SEScalar0To1.h>
#include <biogears/cdm/properties/SEScalarAmountPerVolume.h>
#include <biogears/cdm/properties/SEScalarFlowCompliance.h>
#include <biogears/cdm/properties/SEScalarFlowElastance.h>
#include <biogears/cdm/properties/SEScalarFlowResistance.h>
#include <biogears/cdm/properties/SEScalarFraction.h>
#include <biogears/cdm/properties/SEScalarFrequency.h>
#include <biogears/cdm/properties/SEScalarLength.h>
#include <biogears/cdm/properties/SEScalarMassPerVolume.h>
#include <biogears/cdm/properties/SEScalarNeg1To1.h>
#include <biogears/cdm/properties/SEScalarPressure.h>
#include <biogears/cdm/properties/SEScalarPressurePerVolume.h>
#include <biogears/cdm/properties/SEScalarTime.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SEDrugSystem.h>

#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Controller/BioGears.h>
namespace BGE = mil::tatrc::physiology::biogears;

#pragma warning(disable : 4786)
#pragma warning(disable : 4275)

// #define VERBOSE
namespace biogears {
Nervous::Nervous(BioGears& bg)
  : SENervousSystem(bg.GetLogger())
  , m_data(bg)
{
  Clear();
}

Nervous::~Nervous()
{
  Clear();
}

void Nervous::Clear()
{
  SENervousSystem::Clear();

  m_Patient = nullptr;
  m_Succinylcholine = nullptr;
  m_Sarin = nullptr;
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Initializes system properties to valid homeostatic values.
//--------------------------------------------------------------------------------------------------
void Nervous::Initialize()
{
  BioGearsSystem::Initialize();
  m_FeedbackActive = false;
  m_TestBaroreceptors = false;
  m_blockActive = false;

  m_BaroreceptorFatigueScale = 0.0;
  m_BaroreceptorFrequencyBaseline_Hz = 97.0; //From Ottesen
  m_CentralVentilationDelta_L_Per_min = 0.0;
  m_ChemoreceptorFiringRate_Hz = 3.65;
  m_ChemoreceptorFiringRateSetPoint_Hz = m_ChemoreceptorFiringRate_Hz;
  m_PeripheralBloodGasInteractionBaseline_Hz = 0.0;
  m_PeripheralVentilationDelta_L_Per_min = 0.0;
  m_PreviousMeanArterialPressure_mmHg = m_data.GetCardiovascular().GetMeanArterialPressure(PressureUnit::mmHg);
  //Test Values
  m_FilteredPressure = m_PreviousMeanArterialPressure_mmHg;
  m_HeartRateModSympathetic = 0.0;
  m_HeartRateModVagal = 0.0;
  m_ElastanceMod = 0.0;
  m_ComplianceMod = 0.0;
  m_ResistanceMod = 0.0;

  double tuneTime = 10.0;

  SetBaroreceptorFrequencyComponents(std::vector<double>(3), FrequencyUnit::Hz);
  GetBaroreceptorHeartRateScale().SetValue(1.0);
  GetBaroreceptorHeartElastanceScale().SetValue(1.0);
  GetBaroreceptorResistanceScale().SetValue(1.0);
  GetBaroreceptorComplianceScale().SetValue(1.0);
  GetLeftEyePupillaryResponse().GetSizeModifier().SetValue(0);
  GetLeftEyePupillaryResponse().GetReactivityModifier().SetValue(0);
  GetRightEyePupillaryResponse().GetSizeModifier().SetValue(0);
  GetRightEyePupillaryResponse().GetReactivityModifier().SetValue(0);
  GetPainVisualAnalogueScale().SetValue(0.0);
}

bool Nervous::Load(const CDM::BioGearsNervousSystemData& in)
{
  if (!SENervousSystem::Load(in))
    return false;
  BioGearsSystem::LoadState();
  // We assume state have to be after all stabilization
  m_FeedbackActive = true;
  m_ArterialOxygenSetPoint_mmHg = in.ArterialOxygenSetPoint_mmHg();
  m_ArterialCarbonDioxideSetPoint_mmHg = in.ArterialCarbonDioxideSetPoint_mmHg();
  m_BaroreceptorFatigueScale = in.BaroreceptorFatigueScale();
  m_BaroreceptorFrequencyBaseline_Hz = in.BaroreceptorFrequencyBaseline_Hz();
  m_CentralVentilationDelta_L_Per_min = in.ChemoreceptorCentralVentilationDelta_L_Per_min();
  m_ChemoreceptorFiringRate_Hz = in.ChemoreceptorFiringRate_Hz();
  m_ChemoreceptorFiringRateSetPoint_Hz = in.ChemoreceptorFiringRateSetPoint_Hz();
  m_PeripheralBloodGasInteractionBaseline_Hz = in.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz();
  m_PeripheralVentilationDelta_L_Per_min = in.ChemoreceptorPeripheralVentilationDelta_L_Per_min();
  m_PreviousMeanArterialPressure_mmHg = in.PreviousMeanArterialPressure_mmHg();

  return true;
}
CDM::BioGearsNervousSystemData* Nervous::Unload() const
{
  CDM::BioGearsNervousSystemData* data = new CDM::BioGearsNervousSystemData();
  Unload(*data);
  return data;
}
void Nervous::Unload(CDM::BioGearsNervousSystemData& data) const
{
  SENervousSystem::Unload(data);
  data.ArterialOxygenSetPoint_mmHg(m_ArterialOxygenSetPoint_mmHg);
  data.ArterialCarbonDioxideSetPoint_mmHg(m_ArterialCarbonDioxideSetPoint_mmHg);
  data.BaroreceptorFatigueScale(m_BaroreceptorFatigueScale);
  data.BaroreceptorFrequencyBaseline_Hz(m_BaroreceptorFrequencyBaseline_Hz);
  data.ChemoreceptorCentralVentilationDelta_L_Per_min(m_CentralVentilationDelta_L_Per_min);
  data.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz(m_PeripheralBloodGasInteractionBaseline_Hz);
  data.ChemoreceptorFiringRate_Hz(m_ChemoreceptorFiringRate_Hz);
  data.ChemoreceptorFiringRateSetPoint_Hz(m_ChemoreceptorFiringRateSetPoint_Hz);
  data.ChemoreceptorPeripheralVentilationDelta_L_Per_min(m_PeripheralVentilationDelta_L_Per_min);
  data.PreviousMeanArterialPressure_mmHg(m_PreviousMeanArterialPressure_mmHg);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Initializes the nervous specific quantities
///
/// \details
/// Initializes the nervous system.
//--------------------------------------------------------------------------------------------------
void Nervous::SetUp()
{
  m_dt_s = m_data.GetTimeStep().GetValue(TimeUnit::s);
  m_normalizedGammaHeartRate = m_data.GetConfiguration().GetNormalizedHeartRateIntercept();
  m_normalizedGammaElastance = m_data.GetConfiguration().GetNormalizedHeartElastanceIntercept();
  m_normalizedGammaCompliance = m_data.GetConfiguration().GetNormalizedComplianceIntercept();
  m_normalizedGammaResistance = m_data.GetConfiguration().GetNormalizedResistanceIntercept();
  m_normalizedAlphaHeartRate = m_data.GetConfiguration().GetNormalizedHeartRateSympatheticSlope();
  m_normalizedAlphaElastance = m_data.GetConfiguration().GetNormalizedHeartElastanceSympatheticSlope();
  m_normalizedAlphaCompliance = m_data.GetConfiguration().GetNormalizedComplianceParasympatheticSlope();
  m_normalizedAlphaResistance = m_data.GetConfiguration().GetNormalizedResistanceSympatheticSlope();
  m_normalizedBetaHeartRate = m_data.GetConfiguration().GetNormalizedHeartRateParasympatheticSlope();
  m_Succinylcholine = m_data.GetSubstances().GetSubstance("Succinylcholine");
  m_Sarin = m_data.GetSubstances().GetSubstance("Sarin");
  m_Patient = &m_data.GetPatient();

  m_painStimulusDuration_s = 0.0;
  m_painVASDuration_s = 0.0;
  m_painVAS = 0.0;

  // Set when feedback is turned on
  m_ArterialOxygenSetPoint_mmHg = 0;
  //Reset after stabilization
  m_ArterialCarbonDioxideSetPoint_mmHg = 40.0;

  //Testing combined signal--move to serializeable if successful
  m_AfferentBaroreceptor_Hz = 97.0;
  m_AfferentPulmonary_Hz = 11.0;
  m_FilteredPulmonaryVenousPressure_Hz = 6.0; //May need to tune depending on BG baseline for LeftAtrium
  m_AfferentAtrial_Hz = 9.0; //Same as above
  m_AfferentThermal_Hz = 5.0; //Constant for now
  m_SympathethicNode_Hz = 4.5;
  m_SympatheticPeripheral_Hz = 4.5;
  m_Vagal_Hz = 4.5;
  m_AchSympatheticNode = 4.5;
  m_AchSympatheticPeripheral = 4.5 ;
  m_AchVagal = 4.5;
  m_AdrenalMedulla = 5.0;
  m_HeartRateEffectors = std::vector<double>{0.4, -0.2, -0.05};
  m_HeartElastanceEffectors = std::vector<double>{-0.02, 0.03, 0.001, 0.16};
  m_ResistanceEffectors = std::vector<double>{2.1, 1.4};
}

void Nervous::AtSteadyState()
{
  if (m_data.GetState() == EngineState::AtInitialStableState) {
    m_FeedbackActive = true;
  }
  if (m_data.GetState() == EngineState::AtSecondaryStableState) {
    for (auto component : GetBaroreceptorFrequencyComponents()) {
      m_BaroreceptorFrequencyBaseline_Hz += component->GetValue(FrequencyUnit::Hz);
    }
    SetBaroreceptorFrequencyComponents(std::vector<double>(3), FrequencyUnit::Hz);
  }

  // The set-points (Baselines) get reset at the end of each stabilization period.
  m_ArterialOxygenSetPoint_mmHg = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg);
  m_ArterialCarbonDioxideSetPoint_mmHg = m_data.GetBloodChemistry().GetArterialCarbonDioxidePressure(PressureUnit::mmHg);

  //Central and peripheral ventilation changes are set to 0 because patient baseline ventilation is updated to include
  //their contributions at steady state.
  m_CentralVentilationDelta_L_Per_min = 0.0;
  m_PeripheralVentilationDelta_L_Per_min = 0.0;
  //The chemoreceptor firing rate and its setpoint are reset so that central and peripheral derivatives will evaluate to 0
  //the first time step after stabilization (and will stay that way, assuming no other perturbations to blood gas levels)
  m_ChemoreceptorFiringRateSetPoint_Hz = m_ChemoreceptorFiringRate_Hz;
  m_ChemoreceptorFiringRate_Hz = m_PeripheralBloodGasInteractionBaseline_Hz;

  // The baroreceptor scales need to be reset any time the baselines are reset.
  GetBaroreceptorHeartRateScale().SetValue(1.0);
  GetBaroreceptorHeartElastanceScale().SetValue(1.0);
  GetBaroreceptorResistanceScale().SetValue(1.0);
  GetBaroreceptorComplianceScale().SetValue(1.0);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Preprocess methods for the nervous system
///
/// \details
/// Computes nervous system regulation of the body.
/// Baroreceptor and chemoreceptor feedback is computed and modifiers set in preparation for systems processing.
//--------------------------------------------------------------------------------------------------
void Nervous::PreProcess()
{
  CheckPainStimulus();
  AfferentResponse();
  SympatheticSignalProcess();
  ParasympatheticSignalProcess();
  EfferentResponse();

  //BaroreceptorFeedback();
  //UrsinoBaroreceptor();
  //ChemoreceptorFeedback();
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Nervous Process Step
///
/// \details
/// The only current Process-specific function checks the brain status to set events.
//--------------------------------------------------------------------------------------------------
void Nervous::Process()
{
  CheckNervousStatus();
  SetPupilEffects();
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Nervous PostProcess Step
///
/// \details
/// Currently no nervous postprocess methods.
//--------------------------------------------------------------------------------------------------
void Nervous::PostProcess()
{
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Calculates the baroreceptor feedback and sets the scaling parameters in the CDM
///
/// \details
/// The baroreceptor feedback function uses the current mean arterial pressure relative to the mean arterial
/// pressure set-point in order to calculate the sympathetic and parasympathetic response fractions.
/// These fractions are used to update the scaling parameters of heart rate, heart elastance, resistance and compliance
/// for each time step.
//--------------------------------------------------------------------------------------------------
/// \todo Add decompensation. Perhaps a reduction in the effect that is a function of blood volume below a threshold... and maybe time.
void Nervous::BaroreceptorFeedback()
{
  //Determine dP/dt and update last arterial pressure
  double alpha = 0.05;
  if (m_data.GetState() < EngineState::AtSecondaryStableState) {
    alpha = 0.01;
  }
  double arterialPressure_mmHg = m_data.GetCardiovascular().GetArterialPressure(PressureUnit::mmHg);
  double dMeanPressure_mmHg_Per_s = alpha * (arterialPressure_mmHg - m_PreviousMeanArterialPressure_mmHg);
  m_PreviousMeanArterialPressure_mmHg += dMeanPressure_mmHg_Per_s * m_dt_s;

  if (!m_FeedbackActive) {
    m_data.GetDataTrack().Probe("BaroreceptorComponents", GetBaroreceptorFrequencyComponents(FrequencyUnit::Hz));
    m_data.GetDataTrack().Probe("Baroreceptor-Combined", m_BaroreceptorFrequencyBaseline_Hz);
    m_data.GetDataTrack().Probe("Baroreceptor-Baseline", m_BaroreceptorFrequencyBaseline_Hz);
    m_data.GetDataTrack().Probe("MeanPressureInput", m_PreviousMeanArterialPressure_mmHg);
    return;
  }

  //New model constants
  std::vector<double> tau_s{ 0.82, 5.58, 262.0 };
  std::vector<double> gain_Hz_Per_mmHg{ 3.50, 0.73, 1.41 };
  double maxFrequency_Hz = 120.0;

  //Get previous component frequencies and total frequency (sum of components and baseline)
  std::vector<double> baroreceptorComponentFrequencies = GetBaroreceptorFrequencyComponents(FrequencyUnit::Hz);
  double combinedFrequencySignal_Hz = m_BaroreceptorFrequencyBaseline_Hz + GeneralMath::VectorSum(baroreceptorComponentFrequencies);
  m_AfferentBaroreceptor_Hz = combinedFrequencySignal_Hz; //- maxFrequency_Hz / 2.0;
  //Calcualte new frequencies
  double dComponentFrequency_Hz_Per_s = 0.0;
  double nextFrequency_Hz = 0.0;
  for (size_t pos = 0; pos < baroreceptorComponentFrequencies.size(); ++pos) {
    dComponentFrequency_Hz_Per_s = gain_Hz_Per_mmHg[pos] * dMeanPressure_mmHg_Per_s * combinedFrequencySignal_Hz * (maxFrequency_Hz - combinedFrequencySignal_Hz) / (std::pow(maxFrequency_Hz / 2.0, 2.0)) - (1.0 / tau_s[pos]) * baroreceptorComponentFrequencies[pos];
    nextFrequency_Hz = baroreceptorComponentFrequencies[pos] + dComponentFrequency_Hz_Per_s * m_dt_s;
    baroreceptorComponentFrequencies[pos] = nextFrequency_Hz;
  }
  //Push new values to CDM
  if (!SetBaroreceptorFrequencyComponents(baroreceptorComponentFrequencies, FrequencyUnit::Hz)) {
    Error("Nervous::BaroreceptorFeedback: Vector length mismatch");
  }

  m_data.GetDataTrack().Probe("BaroreceptorComponents", GetBaroreceptorFrequencyComponents(FrequencyUnit::Hz));
  m_data.GetDataTrack().Probe("Baroreceptor-Combined", combinedFrequencySignal_Hz);
  m_data.GetDataTrack().Probe("Baroreceptor-Baseline", m_BaroreceptorFrequencyBaseline_Hz);
  m_data.GetDataTrack().Probe("MeanPressureInput", m_PreviousMeanArterialPressure_mmHg);

  double meanArterialPressure_mmHg = m_data.GetCardiovascular().GetMeanArterialPressure(PressureUnit::mmHg);
  //First calculate the sympathetic and parasympathetic firing rates:
  double nu = m_data.GetConfiguration().GetResponseSlope();

  //Adjusting the mean arterial pressure set-point to account for cardiovascular drug effects
  double meanArterialPressureSetPoint_mmHg = m_data.GetPatient().GetMeanArterialPressureBaseline(PressureUnit::mmHg) //m_MeanArterialPressureNoFeedbackBaseline_mmHg
    + m_data.GetEnergy().GetExerciseMeanArterialPressureDelta(PressureUnit::mmHg);

  //Adjust the MAP set-point for baroreceptors for anesthetics, opioids, and sedatives.  Other drugs should leave set-point as is.
  for (SESubstance* sub : m_data.GetCompartments().GetLiquidCompartmentSubstances()) {
    if (!sub->HasPD())
      continue;
    if ((sub->GetClassification() == CDM::enumSubstanceClass::Anesthetic) || (sub->GetClassification() == CDM::enumSubstanceClass::Sedative) || (sub->GetClassification() == CDM::enumSubstanceClass::Opioid)) {
      meanArterialPressureSetPoint_mmHg += m_data.GetDrugs().GetMeanBloodPressureChange(PressureUnit::mmHg);
      break;
      //Only want to apply the blood pressure change ONCE (In case there are multiple sedative/opioids/etc)
      ///\TODO:  Look into a better way to implement drug classification search
    }
  }

  //Neurological effects of pain action
  if (m_data.GetActions().GetPatientActions().HasPainStimulus()) {
    double painVAS = GetPainVisualAnalogueScale().GetValue();
    painVAS *= 0.1;
    meanArterialPressureSetPoint_mmHg *= (1 + 0.65 * painVAS);
  }

  double sympatheticFraction = 1.0 / (1.0 + std::pow(meanArterialPressure_mmHg / meanArterialPressureSetPoint_mmHg, nu));
  double parasympatheticFraction = 1.0 / (1.0 + std::pow(meanArterialPressure_mmHg / meanArterialPressureSetPoint_mmHg, -nu));
  if (m_TestBaroreceptors) {
    sympatheticFraction = 1.0 / (1.0 + std::pow(combinedFrequencySignal_Hz / m_BaroreceptorFrequencyBaseline_Hz, nu));
    parasympatheticFraction = 1.0 / (1.0 + std::pow(combinedFrequencySignal_Hz / m_BaroreceptorFrequencyBaseline_Hz, -nu));
  }
  double hrParasympatheticFraction = 1.0 / (1.0 + std::pow(combinedFrequencySignal_Hz / 49.0, nu));
  double hrSympatheticFraction = 1.0 / (1.0 + std::pow(combinedFrequencySignal_Hz / 49.0, -nu));

  //Currently baroreceptor fatigue has only been tested for severe infections leading to sepsis.  We only accumulate fatigue if the sympathetic
  // outflow is above a certain threshold.  If we drop below threshold, we allow fatigue parameter to return towards 0. Note that even if infetion
  //is eliminated, the inflammation source will still be found (which we want so that inflammatory model has time to return to baseline).  
  const double fatigueThreshold = 0.65;
  const double fatigueTimeConstant_hr = 2.0;
  double fatigueInput = sympatheticFraction - fatigueThreshold;
  double dFatigueScale = 0.0;
  if (m_data.GetBloodChemistry().GetInflammatoryResponse().HasInflammationSource(CDM::enumInflammationSource::Infection)) {
    if (fatigueInput > 0.0) {
      dFatigueScale = (1.0 / fatigueTimeConstant_hr) * (fatigueInput) * (1.2 - m_BaroreceptorFatigueScale);
    } else if (m_BaroreceptorFatigueScale > ZERO_APPROX) {
      dFatigueScale = (-2.0 * m_BaroreceptorFatigueScale / fatigueTimeConstant_hr);
    }
    m_BaroreceptorFatigueScale += (dFatigueScale * m_data.GetTimeStep().GetValue(TimeUnit::hr));
  }
  
  //Calculate the normalized change in heart rate
  double normalizedHeartRate = GetBaroreceptorHeartRateScale().GetValue();
  double tauHeartRate_s = m_data.GetConfiguration().GetHeartRateDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedHeartRate = (1.0 / tauHeartRate_s) * (-normalizedHeartRate + m_normalizedAlphaHeartRate * sympatheticFraction - m_normalizedBetaHeartRate * parasympatheticFraction + m_normalizedGammaHeartRate) * m_dt_s;
  //if (m_TestBaroreceptors) {
  //  deltaNormalizedHeartRate = (1.0 / tauHeartRate_s) * (-normalizedHeartRate + m_normalizedAlphaHeartRate * hrSympatheticFraction - m_normalizedBetaHeartRate * hrParasympatheticFraction + m_normalizedGammaHeartRate) * m_dt_s;
  //}
  normalizedHeartRate += deltaNormalizedHeartRate;
  //GetBaroreceptorHeartRateScale().SetValue(normalizedHeartRate);

  //Calculate the normalized change in heart elastance
  double normalizedHeartElastance = GetBaroreceptorHeartElastanceScale().GetValue();
  double tauElastance_s = m_data.GetConfiguration().GetHeartElastanceDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedHeartElastance = (1.0 / tauElastance_s) * (-normalizedHeartElastance + m_normalizedAlphaElastance * sympatheticFraction + m_normalizedGammaElastance) * m_dt_s;
  normalizedHeartElastance += deltaNormalizedHeartElastance;
  //GetBaroreceptorHeartElastanceScale().SetValue(normalizedHeartElastance);

  //Calculate the normalized change in flow resistance for any cardiovascular resistor
  double normalizedResistance = GetBaroreceptorResistanceScale().GetValue();
  double tauResistance_s = m_data.GetConfiguration().GetSystemicResistanceDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedResistance = (1.0 / tauResistance_s) * (-normalizedResistance + (m_normalizedAlphaResistance - m_BaroreceptorFatigueScale) * sympatheticFraction + m_normalizedGammaResistance) * m_dt_s;
  normalizedResistance += (deltaNormalizedResistance);
  //GetBaroreceptorResistanceScale().SetValue(normalizedResistance);

  //Calculate the normalized change in flow compliance for any cardiovascular compliance
  double normalizedCompliance = GetBaroreceptorComplianceScale().GetValue();
  double tauCompliance_s = m_data.GetConfiguration().GetVenousComplianceDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedCompliance = (1.0 / tauCompliance_s) * (-normalizedCompliance + m_normalizedAlphaCompliance * parasympatheticFraction + m_normalizedGammaCompliance) * m_dt_s;
  normalizedCompliance += deltaNormalizedCompliance;
  GetBaroreceptorComplianceScale().SetValue(normalizedCompliance);
#ifdef VERBOSE
  m_data.GetDataTrack().Probe("normalizedHeartRate", normalizedHeartRate);
  m_data.GetDataTrack().Probe("normalizedHeartElastance", normalizedHeartElastance);
  m_data.GetDataTrack().Probe("normalizedResistance", normalizedResistance);
  m_data.GetDataTrack().Probe("normalizedCompliance", normalizedCompliance);
  m_data.GetDataTrack().Probe("meanArterialPressureSetPoint_mmHg", meanArterialPressureSetPoint_mmHg);
#endif
}

void Nervous::UrsinoBaroreceptor()
{
  //Continue to use time-weighted mean pressure
  double alpha = 0.1;
  double arterialPressure_mmHg = m_data.GetCardiovascular().GetArterialPressure(PressureUnit::mmHg);
  double dMeanPressure_mmHg_Per_s = alpha * (arterialPressure_mmHg - m_PreviousMeanArterialPressure_mmHg);
  m_PreviousMeanArterialPressure_mmHg += dMeanPressure_mmHg_Per_s * m_dt_s;
  //Afferent parameters
  double tauP = 2.077;
  double tauZ = 6.37;
  double fMin = 2.52;
  double fMax = 47.78;
  double Pn = m_data.GetPatient().GetMeanArterialPressureBaseline(PressureUnit::mmHg); //92.0;
  double kA = 11.758;
  //Efferent sympthetic
  double fES_inf = 2.10;
  double fES0 = 16.11;
  double kES = 0.0675;
  double fES_min = 2.66;
  //Efferent vagal
  double fEV0 = 3.2;
  double fEV_inf = 6.3;
  double kEV = 7.06;
  double fCS0 = 25.0;
  //Effectors
  double elastanceGain = 0.475 / 2.392;
  double resistanceGain = 0.53 / 0.78;
  double heartRateSympatheticGain = -0.13 / 0.58;
  double heartRateVagalGain = 0.09 / 0.58;
  double complianceGain = -132.5 / 1537.0;
  double tauElastance = 8.0;
  double tauResistance = 6.0;
  double tauHeartRateSympathetic = 2.0;
  double tauHeartRateVagal = 1.5;
  double tauCompliance = 20.0;
  //Afferent calculation
  double dFilterPressure = (1.0 / tauP) * (m_PreviousMeanArterialPressure_mmHg + tauZ * dMeanPressure_mmHg_Per_s - m_FilteredPressure);
  m_FilteredPressure += dFilterPressure * m_dt_s;
  double fAfferent = fMin + fMax * std::exp((m_FilteredPressure - Pn) / kA) / (1.0 + std::exp((m_FilteredPressure - Pn) / kA));
  //Efferent Sympathetic calculation
  double fEfferentSym = fES_inf + (fES0 - fES_inf) * std::exp(-kES * fAfferent);
  double fESBase = fES_inf + (fES0 - fES_inf) * std::exp(-kES * fCS0);
  //Efferent Vagal calculation
  double fEfferentVagal = fEV0 + fEV_inf * std::exp((fAfferent - fCS0) / kEV) / (1.0 + std::exp((fAfferent - fCS0) / kEV));
  double fEVBase = (fEV0 + fEV_inf) / 2.0;
  //Input to effectors
  double sympatheticInput = 0.0;
  if (fEfferentSym > fES_min) {
    sympatheticInput = std::log(fEfferentSym - fES_min + 1.0);
  }
  double inputBase = std::log(fESBase - fES_min + 1.0);
  //Elastance
  double dElastanceScale = (1.0 / tauElastance) * (-m_ElastanceMod + elastanceGain * sympatheticInput);
  double elastanceBase = elastanceGain * inputBase;
  m_ElastanceMod += dElastanceScale;
  // GetBaroreceptorHeartElastanceScale().SetValue(1.0 + m_ElastanceMod);
  //Resistance
  double dResistanceScale = (1.0 / tauResistance) * (-m_ResistanceMod + resistanceGain * sympatheticInput);
  double resistanceBase = resistanceGain * inputBase;
  m_ResistanceMod += dResistanceScale;
  // GetBaroreceptorResistanceScale().SetValue(1.0 + m_ResistanceMod);
  //Compliance
  double dComplianceScale = (1.0 / tauCompliance) * (-m_ComplianceMod + complianceGain * sympatheticInput);
  double complianceBase = complianceGain * inputBase;
  m_ComplianceMod += dComplianceScale;
  //  GetBaroreceptorComplianceScale().SetValue(1.0 + m_ComplianceMod);
  //Heart Rate
  double dHeartScaleSym = (1.0 / tauHeartRateSympathetic) * (-m_HeartRateModSympathetic + heartRateSympatheticGain * sympatheticInput);
  m_HeartRateModSympathetic += dHeartScaleSym;
  double dHeartScaleVagal = (1.0 / tauHeartRateVagal) * (-m_HeartRateModVagal + heartRateVagalGain * fEfferentVagal);
  m_HeartRateModVagal += dHeartScaleVagal;
  double heartRateBase = heartRateSympatheticGain * inputBase + heartRateVagalGain * fEVBase;
  // GetBaroreceptorHeartRateScale().SetValue(1.0 + m_HeartRateModSympathetic + m_HeartRateModVagal);

  m_data.GetDataTrack().Probe("Ursino-AfferentRate", fAfferent);
  m_data.GetDataTrack().Probe("Ursino-EfferentSym", fEfferentSym);
  m_data.GetDataTrack().Probe("Ursino-EfferentVagal", fEfferentVagal);
  m_data.GetDataTrack().Probe("Ursion-FilteredPressure", m_FilteredPressure);
  m_data.GetDataTrack().Probe("Ursino-Elastance", 1.0 - elastanceBase + m_ElastanceMod);
  m_data.GetDataTrack().Probe("Ursino-Resistance", 1.0 - resistanceBase + m_ResistanceMod);
  m_data.GetDataTrack().Probe("Ursino-Compliance", 1.0 - complianceBase + m_ComplianceMod);
  m_data.GetDataTrack().Probe("Ursino-HeartRate", 1.0 - heartRateBase + m_HeartRateModSympathetic + m_HeartRateModVagal);
  m_data.GetDataTrack().Probe("MeanPressureInput", m_PreviousMeanArterialPressure_mmHg);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Calculates the patient pain response due to stimulus, susceptibility and drugs
///
/// \details
/// A patient reacts to a noxious stimulus in a certain way. Generally this is reported as a VAS
/// scale value. This value is generally reported by the patient after the nervous system has already parsed
/// the stimulus. For a robotic manikin trainer we need to determine the nervous system and systemic responses
/// related to that stimulus
//--------------------------------------------------------------------------------------------------
void Nervous::CheckPainStimulus()
{
  //Screen for both external pain stimulus and presence of inflammation
  if (!m_data.GetActions().GetPatientActions().HasPainStimulus() && !m_data.GetBloodChemistry().GetInflammatoryResponse().HasInflammationSources()) {
    GetPainVisualAnalogueScale().SetValue(0.0);
    return;
  }

  //initialize:
  SEPainStimulus* p;
  const std::map<std::string, SEPainStimulus*>& pains = m_data.GetActions().GetPatientActions().GetPainStimuli();
  double patientSusceptability = m_Patient->GetPainSusceptibility().GetValue();
  double susceptabilityMapping = GeneralMath::LinearInterpolator(-1.0, 1.0, 0.0, 2.0, patientSusceptability); //mapping [-1,1] -> [0, 2] for scaling the pain stimulus
  double severity = 0.0;
  double painVASMapping = 0.0; //for each location map the [0,1] severity to the [0,10] VAS scale
  double tempPainVAS = 0.0; //sum, scale and store the patient score
  double CNSPainBuffer = 1.0;

  m_painVAS = GetPainVisualAnalogueScale().GetValue();

  //reset duration if VAS falls below approx zero
  if (m_painVAS == 0.0)
    m_painStimulusDuration_s = 0.0;

  //grab drug effects if there are in the body
  if (m_data.GetDrugs().HasCentralNervousResponse()) {
    double NervousScalar = 10.0;
    double CNSModifier = m_data.GetDrugs().GetCentralNervousResponse().GetValue();
    CNSPainBuffer = exp(-CNSModifier * NervousScalar);
  }

  //determine pain response from inflammation caused by burn trauma
  if (m_data.GetActions().GetPatientActions().HasBurnWound()) {
    double traumaPain = m_data.GetActions().GetPatientActions().GetBurnWound()->GetTotalBodySurfaceArea().GetValue();
    traumaPain *= 20.0; //25% TBSA burn will give pain scale = 5, 40% TBSA will give pain scale = 8.0
    tempPainVAS += (traumaPain * susceptabilityMapping * CNSPainBuffer) / (1 + exp(-m_painStimulusDuration_s + 4.0));
  }

  //iterate over all locations to get a cumulative stimulus and buffer them
  for (auto pain : pains) {
    p = pain.second;
    severity = p->GetSeverity().GetValue();
    painVASMapping = 10.0 * severity;

    tempPainVAS += (painVASMapping * susceptabilityMapping * CNSPainBuffer) / (1 + exp(-m_painStimulusDuration_s + 4.0)); //temp time will increase so long as a stimulus is present
  }

  //advance time over the duration of the stimulus

  if (severity < ZERO_APPROX)
    m_painVASDuration_s += m_dt_s;

  m_painStimulusDuration_s += m_dt_s;

  //set the VAS data:
  if (tempPainVAS > 10)
    tempPainVAS = 10.0;

  GetPainVisualAnalogueScale().SetValue(tempPainVAS);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Checks metrics in the nervous system to determine events to be thrown.  Currently includes brain status
/// and presence of fasciculation.
///
/// \details
/// Intracranial pressure is checked to determine if the patient has Intracranial Hyper/hypotension
/// Fasciculation can occur as a result of calcium/magnesium deficiency
/// (or other electrolyte imbalances),succinylcholine, nerve agents, ALS
/// Currently, only fasciculations due to the nerve agent Sarin are active.  Other causes are a subject of model improvement
//------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
void Nervous::CheckNervousStatus()
{
  //-----Check Brain Status-----------------
  double icp_mmHg = m_data.GetCardiovascular().GetIntracranialPressure().GetValue(PressureUnit::mmHg);

  //Intracranial Hypertension
  if (icp_mmHg > 25.0) // \cite steiner2006monitoring
  {
    /// \event Patient: Intracranial Hypertension. The intracranial pressure has risen above 25 mmHg.
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::IntracranialHypertension, true, m_data.GetSimulationTime());
  } else if (m_data.GetPatient().IsEventActive(CDM::enumPatientEvent::IntracranialHypertension) && icp_mmHg < 24.0) {
    /// \event Patient: End Intracranial Hypertension. The intracranial pressure has fallen below 24 mmHg.
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::IntracranialHypertension, false, m_data.GetSimulationTime());
  }

  //Intracranial Hypotension
  if (icp_mmHg < 7.0) // \cite steiner2006monitoring
  {
    /// \event Patient: Intracranial Hypotension. The intracranial pressure has fallen below 7 mmHg.
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::IntracranialHypotension, true, m_data.GetSimulationTime());
  } else if (m_data.GetPatient().IsEventActive(CDM::enumPatientEvent::IntracranialHypotension) && icp_mmHg > 7.5) {
    /// \event Patient: End Intracranial Hypotension. The intracranial pressure has risen above 7.5 mmHg.
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::IntracranialHypertension, false, m_data.GetSimulationTime());
  }

  //------Fasciculations:-------------------------------------------

  //----Fasciculations due to calcium deficiency (inactive)----------------------------------
  /*if (m_Muscleintracellular.GetSubstanceQuantity(*m_Calcium)->GetConcentration(MassPerVolumeUnit::g_Per_L) < 1.0)
    {
    /// \event Patient: Patient is fasciculating due to calcium deficiency
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, true, m_data.GetSimulationTime());
    }
    else if (m_Muscleintracellular.GetSubstanceQuantity(*m_Calcium)->GetConcentration(MassPerVolumeUnit::g_Per_L) > 3.0)
    {
    m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, false, m_data.GetSimulationTime());
    }*/

  //-----Fasciculations due to Sarin--------------------------------------------------
  //Occurs due to inhibition of acetylcholinesterase, the enzyme which breaks down the neurotransmitter acetylcholine
  double RbcAche_mol_Per_L = m_data.GetBloodChemistry().GetRedBloodCellAcetylcholinesterase(AmountPerVolumeUnit::mol_Per_L);
  double RbcFractionInhibited = 1.0 - RbcAche_mol_Per_L / (8e-9); //8 nM is the baseline activity of Rbc-Ache
  if (m_data.GetSubstances().IsActive(*m_Sarin)) {
    ///\cite nambda1971cholinesterase
    //The above study found that individuals exposed to the organophosphate parathion did not exhibit fasciculation until at least
    //80% of Rbc-Ache was inhibited.  This was relaxed to 70% because BioGears is calibrated to throw an irreversible state at
    //100% inhibition when, in actuality, a patient with 100% rbc-ache inhibition will likely survive (rbc-ache thought to act as a buffer
    //for neuromuscular ache)
    if (RbcFractionInhibited > 0.7)
      m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, true, m_data.GetSimulationTime());
    else if ((m_data.GetSubstances().IsActive(*m_Sarin)) && (RbcFractionInhibited < 0.68)) {
      //Oscillations around 70% rbc-ache inhibition are highly unlikely but give some leeway for reversal just in case
      m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, false, m_data.GetSimulationTime());
    }
  }
  //----Fasciculations due to Succinylcholine administration.---------------------------------------------------
  //No evidence exists for a correlation between the plasma concentration of succinylcholine
  //and the degree or presence of fasciculation.  Rather, it has been observed that transient fasciculation tends to occur in most patients after initial dosing
  //(particularly at a dose of 1.5 mg/kg, see refs below), subsiding once depolarization at neuromuscular synapses is accomplished.  Therefore, we model this
  //effect by initiating fasciculation when succinylcholine enters the body and removing it when the neuromuscular block level (calculated in Drugs.cpp) reaches
  //90% of maximum.  To prevent fasciculation from being re-flagged as succinylcholine leaves the body and the block dissipates, we use a sentinel (m_blockActive,
  //initialized to FALSE) so that the event cannot be triggered more than once.
  /// \cite @appiah2004pharmacology, @cite mcloughlin1994influence
  double neuromuscularBlockLevel = m_data.GetDrugs().GetNeuromuscularBlockLevel().GetValue();
  if (m_data.GetSubstances().IsActive(*m_Succinylcholine) && (neuromuscularBlockLevel > 0.0)) {
    if ((neuromuscularBlockLevel < 0.9) && (!m_blockActive))
      m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, true, m_data.GetSimulationTime());
    else {
      m_data.GetPatient().SetEvent(CDM::enumPatientEvent::Fasciculation, false, m_data.GetSimulationTime());
      m_blockActive = true;
    }
  }
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Calculates the chemoreceptor feedback and sets the scaling parameters in the CDM
///
/// \details
/// The chemoreceptor feedback function uses the current arterial partial pressure of oxygen and carbon dioxide
/// relative to the partial pressure set-points in order to calculate response signal.
/// The affected systems identify the signal and adjust accordingly. Note that chemoreception
/// is currently built into the respiratory driver; therefore, the chemoreceptor feedback only sets CV modifiers.
//--------------------------------------------------------------------------------------------------
void Nervous::ChemoreceptorFeedback()
{
  //-------Respiratory Feedback:  This is active throughtout the simulation (including stabilization)------------------------------
  //Model Parameters from
  ///\@cite Magosso2001Mathematical
  const double centralTimeConstant_s = 180.0;
  double centralGainConstant_L_Per_min_mmHg = 1.8;
  const double peripheralTimeConstant_s = 13.0;
  const double peripheralGainConstant_L_Per_min_Hz = 3.36;
  const double firingRateMin_Hz = 0.835;
  const double firingRateMax_Hz = 12.3;
  const double oxygenHalfMax_mmHg = 45.0;
  const double oxygenScale_mmHg = 29.27;
  const double gasInteractionBase = 3.0;
  const double firingRateTimeConstant_s = 2.0;
  const double tuningFactor = 1.5;

  //Note that this method uses instantaneous values of blood gas levels, not running averages
  double arterialO2Pressure_mmHg = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg);
  double arterialCO2Pressure_mmHg = m_data.GetBloodChemistry().GetArterialCarbonDioxidePressure(PressureUnit::mmHg);

  //Magosso and Ursino cite findings that central chemoreceptors are less sensitive at sub-normal levels of CO2 than to super-normal levels
  if (arterialCO2Pressure_mmHg < m_ArterialCarbonDioxideSetPoint_mmHg) {
    centralGainConstant_L_Per_min_mmHg = 0.12;
  }

  //The psi parameter captures the combined interactive effect of O2 and CO2 on the peripheral chemoreceptors.  The degree
  //of interaction varies as hypoxia deepens, with CO2 having less impact as O2 levels decrease
  double psiNum = firingRateMax_Hz + firingRateMin_Hz * exp((arterialO2Pressure_mmHg - oxygenHalfMax_mmHg) / oxygenScale_mmHg);
  double psiDen = 1.0 + exp((arterialO2Pressure_mmHg - oxygenHalfMax_mmHg) / oxygenScale_mmHg);
  double gasInteraction;
  if (arterialO2Pressure_mmHg >= 80.0) {
    gasInteraction = gasInteractionBase;
  } else if (arterialO2Pressure_mmHg >= 40.0) {
    gasInteraction = gasInteractionBase - 1.2 * (80.0 - arterialO2Pressure_mmHg) / 30.0;
  } else {
    gasInteraction = 1.4;
  }
  double psi = (psiNum / psiDen) * (gasInteraction * std::log(arterialCO2Pressure_mmHg / m_ArterialCarbonDioxideSetPoint_mmHg) + tuningFactor);

  if (m_data.GetState() < EngineState::AtSecondaryStableState) {
    //This value is continuously updated during stabilization.  When system reaches steady state, it is used to reset the value of m_ChemoreceptorFiringRate_Hz
    //in the differential equation for dFiringRate_Hz so that all derivatives reset to 0 when a stable state is achieved.
    m_PeripheralBloodGasInteractionBaseline_Hz = psi;
  }

  //Apply effects of opioids that depress central nervous activity
  double drugCNSModifier = m_data.GetDrugs().GetCentralNervousResponse().GetValue();

  double centralInput = arterialCO2Pressure_mmHg - m_ArterialCarbonDioxideSetPoint_mmHg;
  double peripheralInput = m_ChemoreceptorFiringRate_Hz - m_ChemoreceptorFiringRateSetPoint_Hz;

  //Evaluate model derivatives pertaining to change in chemoreceptor firing rate, and changes in central and peripheral contributions to ventilation
  double dFiringRate_Hz = (1.0 / firingRateTimeConstant_s) * (-m_ChemoreceptorFiringRate_Hz + psi) * m_dt_s;
  double dCentralVentilation_L_Per_min = (1.0 / centralTimeConstant_s) * (-m_CentralVentilationDelta_L_Per_min + centralGainConstant_L_Per_min_mmHg * centralInput) * m_dt_s;
  double dPeripheralVentilation_L_Per_min = (1.0 / peripheralTimeConstant_s) * (-m_PeripheralVentilationDelta_L_Per_min + peripheralGainConstant_L_Per_min_Hz * peripheralInput) * m_dt_s;

  //Calculate change in ventilation assuming no metabolic effects--The CNS modifier is applied such that at high values the chemoreceptors cannot force a change from baseline
  double nextTargetVentilation_L_Per_min = m_data.GetPatient().GetTotalVentilationBaseline(VolumePerTimeUnit::L_Per_min) + std::exp(-5.0 * drugCNSModifier) *(m_CentralVentilationDelta_L_Per_min + m_PeripheralVentilationDelta_L_Per_min);


  //Apply metabolic effects. The modifier is tuned to achieve the correct respiratory response for near maximal exercise.
  //A linear relationship is assumed for the respiratory effects due to increased metabolic exertion
  double TMR_W = m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::W);
  double BMR_W = m_data.GetPatient().GetBasalMetabolicRate(PowerUnit::W);
  double energyDeficit_W = m_data.GetEnergy().GetEnergyDeficit(PowerUnit::W);
  double metabolicFraction = (TMR_W + energyDeficit_W) / BMR_W;
  double tunedVolumeMetabolicSlope = 0.2; //Tuned fractional increase of the tidal volume due to increased metabolic rate
  double metabolicModifier = 1.0 + tunedVolumeMetabolicSlope * (metabolicFraction - 1.0);
  nextTargetVentilation_L_Per_min *= metabolicModifier;

  // Confirm that the target does not exceed the maximum ventilation (set in configuration).
  // Flag event if max is exceeded and (if event active) check to see if it has been deactivated
  double maximumPulmonaryVentilationRate = m_data.GetConfiguration().GetPulmonaryVentilationRateMaximum(VolumePerTimeUnit::L_Per_min);

  if (nextTargetVentilation_L_Per_min > maximumPulmonaryVentilationRate) {
    nextTargetVentilation_L_Per_min = maximumPulmonaryVentilationRate;
    m_Patient->SetEvent(CDM::enumPatientEvent::MaximumPulmonaryVentilationRate, true, m_data.GetSimulationTime());
  }

  if (nextTargetVentilation_L_Per_min < maximumPulmonaryVentilationRate && m_Patient->IsEventActive(CDM::enumPatientEvent::MaximumPulmonaryVentilationRate)) {
    m_Patient->SetEvent(CDM::enumPatientEvent::MaximumPulmonaryVentilationRate, false, m_data.GetSimulationTime());
  }

  //Final target ventilation
  m_data.GetRespiratory().GetTargetPulmonaryVentilation().SetValue(nextTargetVentilation_L_Per_min, VolumePerTimeUnit::L_Per_min);

  //Update values for next time step
  m_ChemoreceptorFiringRate_Hz += dFiringRate_Hz;
  m_ChemoreceptorFiringRate_Hz = std::max(0.0, m_ChemoreceptorFiringRate_Hz);
  m_CentralVentilationDelta_L_Per_min += dCentralVentilation_L_Per_min;
  m_PeripheralVentilationDelta_L_Per_min += dPeripheralVentilation_L_Per_min;

  //-----Cardiovascular Feedback:  This functionality is currently only active after stabilization.
  if (!m_FeedbackActive)
    return;

  //Heart Rate modifications
  double normalized_pO2 = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg) / m_ArterialOxygenSetPoint_mmHg;
  double normalized_pCO2 = m_data.GetBloodChemistry().GetArterialCarbonDioxidePressure(PressureUnit::mmHg) / m_ArterialCarbonDioxideSetPoint_mmHg;

  // The chemoreceptor heart rate modification function shape parameters.
  // See NervousMethodology documentation for details.
  double amax = -0.1;
  double a50 = 0.5;
  double aeta = -12.;
  double bmax = 1.;
  double b50 = 1.7;
  double beta = 18;
  double cmax = 1.;
  double c50 = 0.65;
  double ceta = -20;
  double dmax = -0.1;
  double d50 = b50;
  double deta = -aeta;

  //Calculate the normalized change in heart rate
  //double HRBaseline_per_min = m_HeartRateNoFeedbackBaseline_per_min;
  // Maximum HR delta is 1.23 times baseline. The derivation of this maximum is described in the NervousMethodology documentation
  //double maxHeartRateDelta = 1.23 * m_HeartRateNoFeedbackBaseline_per_min;
  double maxHeartRateDelta = 1.23 * m_data.GetPatient().GetHeartRateBaseline(FrequencyUnit::Per_min);
  double modifier = GeneralMath::LogisticFunction(amax, a50, aeta, normalized_pCO2);
  modifier += GeneralMath::LogisticFunction(bmax, b50, beta, normalized_pCO2);
  modifier += GeneralMath::LogisticFunction(cmax, c50, ceta, normalized_pO2);
  modifier += GeneralMath::LogisticFunction(dmax, d50, deta, normalized_pO2);

  //Apply central nervous depressant effects (currently only applies to morphine)

  if (drugCNSModifier >= 0.25)
    modifier = 0.0;

  //set to zero if below .1 percent
  if (modifier < 0.001)
    modifier = 0.0;

  GetChemoreceptorHeartRateScale().SetValue(maxHeartRateDelta * modifier);

  // Calculate the normalized change in heart elastance
  double normalizedHeartElastance = 1.0;
  /// \todo Compute and apply chemoreceptor-mediated contractility changes
  GetChemoreceptorHeartElastanceScale().SetValue(normalizedHeartElastance);
}

//--------------------------------------------------------------------------------------------------
/// \brief
/// Sets pupil size and reactivity modifiers based on drug and TBI effects
///
/// \details
/// Modifiers are on a scale between -1 (for reduction in size/reactivity) and 1 (for increase)
/// TBI effects are applied to the eye on the same side of the injury if localized or both if diffuse
/// Drug and TBI pupil effects are simply summed together
//--------------------------------------------------------------------------------------------------
void Nervous::SetPupilEffects()
{
  // Get modifiers from Drugs
  double leftPupilSizeResponseLevel = m_data.GetDrugs().GetPupillaryResponse().GetSizeModifier().GetValue();
  double leftPupilReactivityResponseLevel = m_data.GetDrugs().GetPupillaryResponse().GetReactivityModifier().GetValue();
  double rightPupilSizeResponseLevel = leftPupilSizeResponseLevel;
  double rightPupilReactivityResponseLevel = leftPupilReactivityResponseLevel;

  // Calculate the TBI response
  if (m_data.GetActions().GetPatientActions().HasBrainInjury()) {
    SEBrainInjury* b = m_data.GetActions().GetPatientActions().GetBrainInjury();

    if (b->GetSeverity().GetValue() > 0) {
      double icp_mmHg = m_data.GetCardiovascular().GetIntracranialPressure().GetValue(PressureUnit::mmHg);

      if (b->GetType() == CDM::enumBrainInjuryType::Diffuse) {
        //https://www.wolframalpha.com/input/?i=y%3D(1+%2F+(1+%2B+exp(-2.0*(x+-+24))))+from+18%3Cx%3C28
        leftPupilSizeResponseLevel += (1 / (1 + exp(-2.0 * (icp_mmHg - 20))));
        //https://www.wolframalpha.com/input/?i=y%3D-.001*pow(10,+.27*(x+-+15))+from+18%3Cx%3C28+and+-1%3Cy%3C0
        leftPupilReactivityResponseLevel += -.001 * std::pow(10, .27 * (icp_mmHg - 13));
        rightPupilSizeResponseLevel = leftPupilSizeResponseLevel;
        rightPupilReactivityResponseLevel = leftPupilReactivityResponseLevel;
      } else if (b->GetType() == CDM::enumBrainInjuryType::LeftFocal) {
        leftPupilSizeResponseLevel += (1 / (1 + exp(-2.0 * (icp_mmHg - 20))));
        leftPupilReactivityResponseLevel += -.001 * std::pow(10, .27 * (icp_mmHg - 13));
      } else if (b->GetType() == CDM::enumBrainInjuryType::RightFocal) {
        rightPupilSizeResponseLevel += (1 / (1 + exp(-2.0 * (icp_mmHg - 20))));
        rightPupilReactivityResponseLevel += -.001 * std::pow(10, .27 * (icp_mmHg - 13));
      }
    }
  }

  BLIM(leftPupilSizeResponseLevel, -1, 1);
  BLIM(leftPupilReactivityResponseLevel, -1, 1);
  BLIM(rightPupilSizeResponseLevel, -1, 1);
  BLIM(rightPupilReactivityResponseLevel, -1, 1);
  GetLeftEyePupillaryResponse().GetSizeModifier().SetValue(leftPupilSizeResponseLevel);
  GetLeftEyePupillaryResponse().GetReactivityModifier().SetValue(leftPupilReactivityResponseLevel);
  GetRightEyePupillaryResponse().GetSizeModifier().SetValue(rightPupilSizeResponseLevel);
  GetRightEyePupillaryResponse().GetReactivityModifier().SetValue(rightPupilReactivityResponseLevel);
}

void Nervous::AfferentResponse()
{
  //Generate afferent barorceptor signal (*ab = afferent baroreceptor)
  BaroreceptorFeedback();
  m_data.GetDataTrack().Probe("Afferent_Baroreceptor", m_AfferentBaroreceptor_Hz);
  //Generate afferent chemoreceptor signal (*ac = afferent chemoreceptor)
  ChemoreceptorFeedback();
  double afferentChemoreceptor_Hz = m_ChemoreceptorFiringRate_Hz;
  m_data.GetDataTrack().Probe("Afferent_Chemoreceptor", afferentChemoreceptor_Hz);


  //Generate afferent lung stretch receptor signal (*ap = afferent pulmonary)
  double tauAP_s = 2.0;
  double gainAP_Hz_Per_L = 23.29;
  double tidalVolume_L = m_data.GetRespiratory().GetTidalVolume(VolumeUnit::L);
  double dFrequencyAP_Hz = (1.0 / tauAP_s) * (-m_AfferentPulmonary_Hz + tidalVolume_L * gainAP_Hz_Per_L);
  m_AfferentPulmonary_Hz += dFrequencyAP_Hz * m_dt_s;
  m_data.GetDataTrack().Probe("Afferent_PulmonaryStretch", m_AfferentPulmonary_Hz);

  //Afferent atrial stretch receptors (*aa = afferent atrial)
  double tauAA_s = 6.37;
  double maxAA_Hz = 18.0;
  double kAA_mmHg = 3.429;
  double pulmonarySetpoint_mmHg = 6.0;
  double pulmonaryVenousPressure_mmHg = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetNode(BGE::CardiovascularNode::LeftAtrium1)->GetPressure(PressureUnit::mmHg);
  double dFilteredPulmonaryVenousPressure = (1.0 / tauAA_s) * (pulmonaryVenousPressure_mmHg - m_FilteredPulmonaryVenousPressure_Hz);
  m_FilteredPulmonaryVenousPressure_Hz += dFilteredPulmonaryVenousPressure * m_dt_s;
  double exponentAA = (m_FilteredPulmonaryVenousPressure_Hz - pulmonarySetpoint_mmHg) / kAA_mmHg;
  m_AfferentAtrial_Hz = (maxAA_Hz * std::exp(exponentAA)) / (1.0 + std::exp(exponentAA));

  m_data.GetDataTrack().Probe("AtrialVenousPressure", pulmonaryVenousPressure_mmHg);
  m_data.GetDataTrack().Probe("AtrialVenousPressure_Filtered", m_FilteredPulmonaryVenousPressure_Hz);
  m_data.GetDataTrack().Probe("Afferent_Atrial", m_AfferentAtrial_Hz);

  //Afferent thermal receptors (*AT)--may do these later, but for now leave at baseline so as to not disrupt model
  double AfferentThermal_Hz = 5.0;
  m_AfferentThermal_Hz = AfferentThermal_Hz;
}

void Nervous::SympatheticSignalProcess()
{
  //Basal sympathetic and vascular motor tone (*ES)
  double minES_Hz = 16.11; //Threshold value
  double highES_Hz = 8.35; //Halfway between value for M/F (pg. 133),
  double wES_C = 1.9; //Weight of chemoreceptor effect on tone, applies to both ESS and ESP
  double wES_Low = 3.05; //Weight of low-frequency oscillations on tone
  double rES_Low = 0.05;
  double kES_Low = 3.0; // 2.0 //With rES_Low, determine rate of low frequency oscillations in sympathetic tone
  //Sympathetic tone calc
  double sympatheticTone = highES_Hz + wES_C * m_ChemoreceptorFiringRate_Hz;
  LLIM(sympatheticTone, minES_Hz); //If tone not above threshold, use threshold value
  
  double vascularTone = 16.11;

  m_data.GetDataTrack().Probe("SympatheticToneBasal", sympatheticTone);

  //Sympathetic outflow (*ES) -- these three parameters are common to SA/AV node feedback and peripheral vascular feedback
  double infES_Hz = 2.1;
  double maxES_Hz = 60.0;
  double kES = 0.0675;
  //Sympathetic signal to SA/AV nodes
  double wSN_AB = 1.0, wSN_AC = 1.0, wSN_AP = 0.0, wSN_AA = 0.0, wSN_AT = 1.0; //Weights from arterial baroreceptors, chemoreceptors, pulmonary stretch, atrial stretch, thermal feedback
  double thetaSN = -73.43; //  Derivation:  3.57 (Ursino) - 5 (thermal offset) + 72 (offset in baroreceptor operating points)
  double exponentSN = kES * (-wSN_AB * m_AfferentBaroreceptor_Hz + wSN_AC * m_ChemoreceptorFiringRate_Hz - wSN_AP * m_AfferentPulmonary_Hz - wSN_AA * m_AfferentAtrial_Hz - wSN_AT * m_AfferentThermal_Hz - thetaSN);
  m_SympathethicNode_Hz = infES_Hz + (sympatheticTone - infES_Hz) * std::exp(exponentSN);
  ULIM(m_SympathethicNode_Hz, maxES_Hz);
  //Sympathetic signal to peripheral vascular beds
  double wSV_AB = 0.3, wSV_AC = 5.0, wSV_AP = 0.34, wSV_AA = 1.0, wSV_AT = 1.0;
  double thetaSV = 7.34;
  double exponentSV = kES * (-wSV_AB * m_AfferentBaroreceptor_Hz + wSV_AC * m_ChemoreceptorFiringRate_Hz - wSV_AP * m_AfferentPulmonary_Hz - wSV_AA * m_AfferentAtrial_Hz + wSV_AT * m_AfferentThermal_Hz - thetaSV);
  m_SympatheticPeripheral_Hz = infES_Hz + (vascularTone - infES_Hz) * std::exp(exponentSV);
  ULIM(m_SympatheticPeripheral_Hz, maxES_Hz);

  m_data.GetDataTrack().Probe("SympatheticNodeSignal", m_SympathethicNode_Hz);
  m_data.GetDataTrack().Probe("SympatheticPeripheralSignal", m_SympatheticPeripheral_Hz);

  //Tranlate sympathetic signals to non-dimensionalized neurotransmitter concentrations
  double tauSN_s = 0.72;
  double tauSV_s = 0.72;
  double fES_Res_Hz = 1.0;
  double dAchSN = (1.0 / tauSN_s) * (-m_AchSympatheticNode + m_SympathethicNode_Hz);
  double dAchSV = (1.0 / tauSV_s) * (-m_AchSympatheticPeripheral + std::max(m_SympatheticPeripheral_Hz, fES_Res_Hz));
  m_AchSympatheticNode += dAchSN * m_dt_s;
  m_AchSympatheticPeripheral += dAchSV * m_dt_s;

  m_data.GetDataTrack().Probe("SympatheticNodeACH", m_AchSympatheticNode);
  m_data.GetDataTrack().Probe("SympatheticPeripheralACH", m_AchSympatheticPeripheral);
}

void Nervous::ParasympatheticSignalProcess()
{
  //Parasympathetic / vagal signal path (*EV)
  double maxEV_Hz = 3.2;
  double patientAge = m_data.GetPatient().GetAge().GetValue(TimeUnit::yr);
  double baseEV_Hz = maxEV_Hz - 1.66e-3 * patientAge;
  double infEV_Hz = 6.3;
  double kEV_Hz = 7.06;
  double wEV_AC = 0.2, wEV_AP = 0.105, wEV_AT = 0.105;
  double thetaEV_Hz = 0.2125;  //Derivation: -0.68 (Ursino) + .105 * 5 (thermal addition) + .105 * 3.5 (pulmonary offset for BG)
  double baroreceptorBaseline_Hz = 37.0;
  double exponentEV = (m_AfferentBaroreceptor_Hz - m_BaroreceptorFrequencyBaseline_Hz) / kEV_Hz;
  m_Vagal_Hz = (baseEV_Hz + infEV_Hz * std::exp(exponentEV)) / (1.0 + std::exp(exponentEV));
  m_Vagal_Hz += (wEV_AC * m_ChemoreceptorFiringRate_Hz - wEV_AP * m_AfferentPulmonary_Hz + wEV_AT * m_AfferentThermal_Hz - thetaEV_Hz);
  //Could look in to adding respiratory sinus arrhythmia later
  m_data.GetDataTrack().Probe("VagalSignal", m_Vagal_Hz);
  //Tranlate vagal signal to non-dimensionalized neurotransmitter concentrations
  double tauAchV_s = 1.32;
  double dAchV = (1.0 / tauAchV_s) * (-m_AchVagal + m_Vagal_Hz);
  m_AchVagal += dAchV * m_dt_s;
  m_data.GetDataTrack().Probe("VagalACH", m_AchVagal);
}

void Nervous::EfferentResponse()
{
  //Modulation parameters for Ach/Nor concentrations at effector sites
  double kAchV = 1.0, kNorSN = 1.0, kNorV = 0.1, kNorSV = 1.0;
  double wAM_Ach = 0.166;
  double tauAM_s = 1.32;
  double cAchV = kAchV * m_AchVagal;
  double cNorSN = kNorSN * m_AchSympatheticNode + kNorV * m_AchVagal;
  double cNorSV = kNorSV * m_AchSympatheticPeripheral;
  double dAdrenalMedulla = (1.0 / tauAM_s) * (-m_AdrenalMedulla + wAM_Ach * (m_SympathethicNode_Hz + 5.0 * m_SympatheticPeripheral_Hz)); //5.0 because this value applies to vascular, brain, muscle, splanchnic, extrasplanchnic locations
  m_AdrenalMedulla += dAdrenalMedulla * m_dt_s;
  double cEpiAM = 0.8 * m_AdrenalMedulla; //epineprine from adrenal medulla
  double cNorAM = 0.2 * m_AdrenalMedulla; //norepinephrine from adrenal medulla
  m_data.GetDataTrack().Probe("AM_Epi", cEpiAM);
  m_data.GetDataTrack().Probe("AM_Nor", cNorAM);
  m_data.GetDataTrack().Probe("SympatheticNode_Nor", cNorSN);
  m_data.GetDataTrack().Probe("SympatheticPeripheral_Nor", cNorSV);

  //--------------------------Effector responses----------------------------------
  //Heart Rate
  double tauHrAchV_s = 1.5, tauHrNorSN_s = 2.0, tauHrNorAM_s = 2.0;
  double kSN_V = 0.25;  //1.0;
  double minNorSN = 2.66, minNorAM = 1.0; //2.66;
  double gainHrAchV = 0.09, gainHrNorSN = -0.23, gainHrNorAM = -0.13;
  double fRSA_Hz = 0.15, gainRSA = 0.01;
  double patientAge = m_data.GetPatient().GetAge().GetValue(TimeUnit::yr);
  double heartPeriod0 = 1.0 / (1.97 - 9.5e-3 * patientAge);
  double hrAchV = m_HeartRateEffectors[0];
  double hrNorSN = m_HeartRateEffectors[1];
  double hrNorAM = m_HeartRateEffectors[2];

  double dHrAchV = (1.0 / tauHrAchV_s) * (-hrAchV + gainHrAchV * cAchV);
  double dHrNorSN = 0.0;
  if (cNorSN > (minNorSN + kSN_V * cAchV)) {
    dHrNorSN = (1.0 / tauHrNorSN_s) * (-hrNorSN + gainHrNorSN * std::log(cNorSN - kSN_V * cAchV - minNorSN + 1.0));
    //dHrNorSN = (1.0 / tauHrNorSN_s) * (-hrNorSN + gainHrNorSN * std::log(cNorSN - minNorSN + 1.0));
  }
  double dHrNorAM = 0.0;
  if (cNorAM > minNorAM) {
	dHrNorAM = (1.0 / tauHrNorAM_s) * (-hrNorAM + gainHrNorAM * std::log(cNorAM - minNorAM + 1.0));
  }
  m_HeartRateEffectors[0] = hrAchV + dHrAchV * m_dt_s;
  m_HeartRateEffectors[1] = hrNorSN + dHrNorSN * m_dt_s;
  m_HeartRateEffectors[2] = hrNorAM + dHrNorAM * m_dt_s;

  double nextHeartPeriod = heartPeriod0 + GeneralMath::VectorSum(m_HeartRateEffectors);
  double nextHeartRate = 60.0 / nextHeartPeriod;
  m_data.GetDataTrack().Probe("HeartRateEffectors", m_HeartRateEffectors);
  m_data.GetDataTrack().Probe("HeartRate_Next", nextHeartRate);
  

  //Heart elastance
  double tauElAchV_s = 1.5, tauElNorSN_s = 8.0, tauElNorAM_s = 8.0, tauElEpiAM_s = 8.0;
  double minEpiAM = 2.66;
  //Baseline gains match lit model output but BG does not change elastance as much source model, so we need to tone these down
  double gainElAchV = -0.09 / 20.0, gainElNorSN = 0.475/ 20.0, gainElNorAM = 0.475/20.0, gainElEpiAM = 0.475/20.0;
  double elAchV = m_HeartElastanceEffectors[0];
  double elNorSN = m_HeartElastanceEffectors[1];
  double elNorAM = m_HeartElastanceEffectors[2];
  double elEpiAM = m_HeartElastanceEffectors[3];

  double dElAchV = (1.0 / tauElAchV_s) * (-elAchV + gainElAchV * cAchV);
  double dElNorSN = 0.0;
  double dElNorAM = 0.0;
  double dElEpiAM = 0.0;
  if (cNorSN > minNorSN) {
    dElNorSN = (1.0 / tauElNorSN_s) * (-elNorSN + gainElNorSN * std::log(cNorSN - minNorSN + 1.0));
  }
  if (cNorAM > minNorAM) {
    dElNorAM = (1.0 / tauElNorAM_s) * (-elNorAM + gainElNorAM * std::log(cNorAM - minNorAM + 1.0));
  }
  if (cEpiAM > minEpiAM) {
    dElEpiAM = (1.0 / tauElEpiAM_s) * (-elEpiAM + gainElEpiAM * std::log(cEpiAM - minEpiAM + 1.0));
  }
  m_HeartElastanceEffectors[0] = elAchV + dElAchV * m_dt_s;
  m_HeartElastanceEffectors[1] = elNorSN + dElNorSN * m_dt_s;
  m_HeartElastanceEffectors[2] = elNorAM + dElNorAM * m_dt_s;
  m_HeartElastanceEffectors[3] = elEpiAM + dElEpiAM * m_dt_s;

  double nextHeartElastance = 2.392 + GeneralMath::VectorSum(m_HeartElastanceEffectors);		//2.392 = baseline value

  m_data.GetDataTrack().Probe("HeartElastanceEffectors", m_HeartElastanceEffectors);
  m_data.GetDataTrack().Probe("HeartElastance_Next", nextHeartElastance);

  //Resistance
  double tauRNorSV_s = 6.0, tauREpiAM_s = 6.0;
  double gainRNorSV = 1.94, gainREpiAM = 1.94;
  double minNorSV = 2.66;
  double rNorSV = m_ResistanceEffectors[0];
  double rEpiAM = m_ResistanceEffectors[1];
  double dRNorSN = 0.0;
  double dREpiAM = 0.0;
  if (cNorSN > minNorSN) {
    dRNorSN = (1.0 / tauRNorSV_s) * (-rNorSV + gainRNorSV * std::log(cNorSV - minNorSV + 1.0));
  }
  if (cEpiAM > minEpiAM) {
    dREpiAM = (1.0 / tauREpiAM_s) * (-rEpiAM + gainREpiAM * std::log(cEpiAM - minEpiAM + 1.0));
  }
  m_ResistanceEffectors[0] = rNorSV + dRNorSN * m_dt_s;
  m_ResistanceEffectors[1] = rEpiAM + dREpiAM * m_dt_s;

  double offset1 = 3.476;
  double offset2 = 0.79;
  double peripheralResistanceMod = GeneralMath::VectorSum(m_ResistanceEffectors)-offset1;
  double brainResistanceMod = m_ResistanceEffectors[0] - m_ResistanceEffectors[1]-offset2;

  
  double muscleBase = 2.106;
  double splanchnicBase = 2.49;
  double extraSplanchnicBase = 1.655;
  double brainBase = 6.57;
  double ventricleBase = 2.392;

  m_data.GetDataTrack().Probe("ResistanceEffectors", m_ResistanceEffectors);
  m_data.GetDataTrack().Probe("ResistanceFractionChange_muscle", peripheralResistanceMod / muscleBase);
  m_data.GetDataTrack().Probe("ResistanceFractionChange_splanchnic", peripheralResistanceMod / splanchnicBase);
  m_data.GetDataTrack().Probe("ResistanceFractionChange_extraSplanchnic", peripheralResistanceMod / extraSplanchnicBase);
  m_data.GetDataTrack().Probe("ResistanceFractionChange_brain", brainResistanceMod / brainBase);
  m_data.GetDataTrack().Probe("ResistanceFractionChange_ventricle", brainResistanceMod / ventricleBase);

  if (m_FeedbackActive) {
    GetBaroreceptorHeartRateScale().SetValue(nextHeartRate / m_data.GetPatient().GetHeartRateBaseline(FrequencyUnit::Per_min));
    GetBaroreceptorHeartElastanceScale().SetValue(nextHeartElastance / 2.42);
	GetBaroreceptorResistanceScale().SetValue(1.0 + peripheralResistanceMod / splanchnicBase);
  }

  
}
}