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
  m_CentralVentilationDelta_L_Per_min = 0.0;
  m_ChemoreceptorFiringRate_Hz = 3.65;
  m_ChemoreceptorFiringRateSetPoint_Hz = m_ChemoreceptorFiringRate_Hz;
  m_PeripheralBloodGasInteractionBaseline_Hz = 0.0;
  m_PeripheralVentilationDelta_L_Per_min = 0.0;

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

  GetResistanceScaleCerebral().SetValue(1.0);
  GetResistanceScaleExtrasplanchnic().SetValue(1.0);
  GetResistanceScaleMuscle().SetValue(1.0);
  GetResistanceScaleSplanchnic().SetValue(1.0);
  GetResistanceScaleVentricle().SetValue(1.0);
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
  m_CentralVentilationDelta_L_Per_min = in.ChemoreceptorCentralVentilationDelta_L_Per_min();
  m_ChemoreceptorFiringRate_Hz = in.ChemoreceptorFiringRate_Hz();
  m_ChemoreceptorFiringRateSetPoint_Hz = in.ChemoreceptorFiringRateSetPoint_Hz();
  m_PeripheralBloodGasInteractionBaseline_Hz = in.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz();
  m_PeripheralVentilationDelta_L_Per_min = in.ChemoreceptorPeripheralVentilationDelta_L_Per_min();

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
  data.ChemoreceptorCentralVentilationDelta_L_Per_min(m_CentralVentilationDelta_L_Per_min);
  data.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz(m_PeripheralBloodGasInteractionBaseline_Hz);
  data.ChemoreceptorFiringRate_Hz(m_ChemoreceptorFiringRate_Hz);
  data.ChemoreceptorFiringRateSetPoint_Hz(m_ChemoreceptorFiringRateSetPoint_Hz);
  data.ChemoreceptorPeripheralVentilationDelta_L_Per_min(m_PeripheralVentilationDelta_L_Per_min);

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

  //Test Values
  m_HeartRateModSympathetic = 0.0;
  m_HeartRateModVagal = 0.0;
  m_ElastanceMod = 0.0;
  m_ComplianceMod = 0.0;
  m_ResistanceMod = 0.0;
  m_AfferentBaroreceptor = 25.0;
  m_BaroreceptorBaseline = 0.38;
  m_AfferentPulmonary = 5.75;
  m_FilteredPulmonaryVenousPressure_Hz = 6.0; //May need to tune depending on BG baseline for LeftAtrium
  m_AfferentAtrial = 9.0; //Same as above
  m_SympatheticHeartSignal = 0.175;
  m_SympatheticPeripheralSignal = 0.175;
  m_VagalSignal = 0.80;
  m_SympatheticHeartSignal_Baseline = m_SympatheticHeartSignal;
  m_SympatheticPeripheralSignal_Baseline = m_SympatheticPeripheralSignal;
  m_VagalSignal_Baseline = m_VagalSignal;
  m_AfferentStrain = 0.04;
  m_PressureHalfMax = m_data.GetPatient().GetSystolicArterialPressureBaseline(PressureUnit::mmHg);
  m_HeartElastanceEffector = 0.0;
  m_CerebralBloodFlowBaseline_mL_Per_s = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralCapillariesToCerebralVeins1)->GetFlow(VolumePerTimeUnit::mL_Per_s);
  m_CerebralBloodFlowFilter = m_CerebralBloodFlowBaseline_mL_Per_s;
  m_CerebralAutoEffect = 0.0;
  m_CerebralCO2Baseline_mmHg = m_data.GetCompartments().GetExtracellularFluid(*m_data.GetCompartments().GetTissueCompartment(BGE::TissueCompartment::Brain)).GetSubstanceQuantity(m_data.GetSubstances().GetCO2())->GetPartialPressure(PressureUnit::mmHg);
}

void Nervous::AtSteadyState()
{
  if (m_data.GetState() == EngineState::AtInitialStableState) {
    m_FeedbackActive = true;
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

  m_SympatheticHeartSignal_Baseline = m_SympatheticHeartSignal;
  m_SympatheticPeripheralSignal_Baseline = m_SympatheticPeripheralSignal;
  m_VagalSignal_Baseline = m_VagalSignal;
 

  m_CerebralBloodFlowBaseline_mL_Per_s = m_CerebralBloodFlowFilter;
  m_CerebralCO2Baseline_mmHg = m_data.GetCompartments().GetExtracellularFluid(*m_data.GetCompartments().GetTissueCompartment(BGE::TissueCompartment::Brain)).GetSubstanceQuantity(m_data.GetSubstances().GetCO2())->GetPartialPressure(PressureUnit::mmHg);
  m_CerebralAutoEffect = 0.0;
  m_PressureHalfMax = m_data.GetCardiovascular().GetSystolicArterialPressure(PressureUnit::mmHg);
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
  CentralSignalProcess();
  EfferentResponse();
  CerebralAutoregulation();

  //BaroreceptorFeedback();
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
  m_data.GetDataTrack().Probe("BaroreceptorComponents", GetBaroreceptorFrequencyComponents(FrequencyUnit::Hz));

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
  //Randall Model
  double A = 5.0;
  double voigtK = 0.1;
  double voigtTau = 0.9;
  double slope = 0.04;

  double systolicPressure_mmHg = m_data.GetCardiovascular().GetSystolicArterialPressure(PressureUnit::mmHg);
  double strainExp = std::exp(-slope * (systolicPressure_mmHg - m_PressureHalfMax));
  double wallStrain = 1.0 - std::sqrt((1.0 + strainExp)/(A+strainExp));
  m_data.GetDataTrack().Probe("WallStrain", wallStrain);

  double dStrain = (1.0 / voigtTau) * (-m_AfferentStrain + voigtK * wallStrain);
  m_AfferentStrain += (dStrain * m_dt_s);
  double strainSignal = wallStrain - m_AfferentStrain;

  //Convert strain signal to be on same basis as Ursino model--this puts deviations in wall strain on same basis as other signals
  double fBMax = 47.78;
  double fBMin = 2.52;
  double kB = 0.075;
  double baroExponent = std::exp((strainSignal - m_BaroreceptorBaseline) / kB);

  m_AfferentBaroreceptor = (fBMin + fBMax * baroExponent) / (1.0 + baroExponent);


  //Update setpoint
  if (m_data.GetState() > EngineState::SecondaryStabilization) {
    double kBaro = 7.0e-5;
    double dSetpointAdjust = kBaro * (systolicPressure_mmHg - m_PressureHalfMax);
    m_PressureHalfMax += (dSetpointAdjust * m_dt_s);
  }
  m_data.GetDataTrack().Probe("Randall_Setpoint", m_PressureHalfMax);


  if (!m_FeedbackActive) {
    return;
  }

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
  //GetBaroreceptorComplianceScale().SetValue(normalizedCompliance);
#ifdef VERBOSE
  m_data.GetDataTrack().Probe("normalizedHeartRate", normalizedHeartRate);
  m_data.GetDataTrack().Probe("normalizedHeartElastance", normalizedHeartElastance);
  m_data.GetDataTrack().Probe("normalizedResistance", normalizedResistance);
  m_data.GetDataTrack().Probe("normalizedCompliance", normalizedCompliance);
  m_data.GetDataTrack().Probe("meanArterialPressureSetPoint_mmHg", meanArterialPressureSetPoint_mmHg);
#endif
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

  //GetChemoreceptorHeartRateScale().SetValue(maxHeartRateDelta * modifier);

  // Calculate the normalized change in heart elastance
  double normalizedHeartElastance = 1.0;
  /// \todo Compute and apply chemoreceptor-mediated contractility changes
  //GetChemoreceptorHeartElastanceScale().SetValue(normalizedHeartElastance);
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
  m_data.GetDataTrack().Probe("Afferent_Baroreceptor", m_AfferentBaroreceptor);
  //Generate afferent chemoreceptor signal (*ac = afferent chemoreceptor)
  ChemoreceptorFeedback();
  double afferentChemoreceptor_Hz = m_ChemoreceptorFiringRate_Hz;
  m_data.GetDataTrack().Probe("Afferent_Chemoreceptor", afferentChemoreceptor_Hz);

  //Generate afferent lung stretch receptor signal (*ap = afferent pulmonary)
  double tauAP_s = 2.0;
  double gainAP_Hz_Per_L = 11.25;
  double tidalVolume_L = m_data.GetRespiratory().GetTidalVolume(VolumeUnit::L);
  double dFrequencyAP_Hz = (1.0 / tauAP_s) * (-m_AfferentPulmonary + tidalVolume_L * gainAP_Hz_Per_L);
  if (m_FeedbackActive) {
    //Tidal volume jumps around so much at beginning
    m_AfferentPulmonary += dFrequencyAP_Hz * m_dt_s;
  }

  m_data.GetDataTrack().Probe("Afferent_PulmonaryStretch", m_AfferentPulmonary);

  //Afferent atrial stretch receptors (*aa = afferent atrial)
  double tauAA_s = 6.37;
  double maxAA_Hz = 18.0;
  double kAA_mmHg = 3.429;
  double pulmonarySetpoint_mmHg = 6.0;
  double pulmonaryVenousPressure_mmHg = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetNode(BGE::CardiovascularNode::LeftAtrium1)->GetPressure(PressureUnit::mmHg);
  double dFilteredPulmonaryVenousPressure = (1.0 / tauAA_s) * (pulmonaryVenousPressure_mmHg - m_FilteredPulmonaryVenousPressure_Hz);
  m_FilteredPulmonaryVenousPressure_Hz += dFilteredPulmonaryVenousPressure * m_dt_s;
  double exponentAA = (m_FilteredPulmonaryVenousPressure_Hz - pulmonarySetpoint_mmHg) / kAA_mmHg;
  m_AfferentAtrial = (maxAA_Hz * std::exp(exponentAA)) / (1.0 + std::exp(exponentAA));

  m_data.GetDataTrack().Probe("AtrialVenousPressure", pulmonaryVenousPressure_mmHg);
  m_data.GetDataTrack().Probe("AtrialVenousPressure_Filtered", m_FilteredPulmonaryVenousPressure_Hz);
  m_data.GetDataTrack().Probe("Afferent_Atrial", m_AfferentAtrial);

  //Afferent thermal receptors (*AT)--may do these later, but for now leave at baseline so as to not disrupt model
  double AfferentThermal_Hz = 5.0;
  m_AfferentThermal_Hz = AfferentThermal_Hz;
}

void Nervous::CentralSignalProcess()
{
  //Randall model
  double kPB = 5.0;
  double kS = 5.0;
  double tauPB = 1.8;
  double tauS = 10.0;
  double qPB = 10.0;
  double qS = 10.0;
  double sPB = 0.54;
  double sS = 0.05;

  //double gPB = 1.0 / (1.0 + std::exp(-qPB * (m_AfferentBaroreceptor - sPB)));
  //double gS = 1.0 / (1.0 + std::exp(qS * (m_AfferentBaroreceptor - sS)));
  //double dTPB = (-m_VagalSignal + kPB * gPB) / (tauPB);
  //double dTS = (-m_SympatheticPeripheralSignal + kS * gS) / (tauS);
  //m_VagalSignal += (dTPB * m_dt_s);
  //m_SympatheticPeripheralSignal += (dTS * m_dt_s);

  //m_data.GetDataTrack().Probe("Randall_TPar", m_VagalSignal);
  //m_data.GetDataTrack().Probe("Randall_TSym", m_SympatheticPeripheralSignal);

  //Consolidate test
  double fEs0 = 16.11;
  double fEsMin = 2.66;
  double fEsInf = 2.1;
  double fEsMax = 60.0;
  double kES = 0.0675;

  //Update hypoxia thresholds for sympathetic signals
  double xMinSN = -49.38 - 18.75, xMaxSN = 3.59 - 18.75, xMinSP = 7.33 - 15.26, xMaxSP = 13.32 - 15.26;
  double o2SetSN = 45.0, o2SetSP = 30.0, kO2SN = 6.0, kO2SP = 2.0, tauISC = 30.0;
  double arterialO2 = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg);
  double expActivationSN = std::exp((arterialO2 - o2SetSN) / kO2SN);
  double expActivationSP = std::exp((arterialO2 - o2SetSP) / kO2SP);
  double hypoxiaInputSN = (xMinSN + xMaxSN * expActivationSN) / (1.0 + expActivationSN);
  double hypoxiaInputSP = (xMinSP + xMaxSP * expActivationSP) / (1.0 + expActivationSP);
  //Sympathetic signal to SA/AV nodes
  double wSN_AB = -1.0, wSN_AC = 1.0, wSN_AT = -1.0; //Weights from arterial baroreceptors, chemoreceptors, pulmonary stretch, atrial stretch, thermal feedback
  double thetaSN = 3.59;
  double exponentSN = kES * (wSN_AB * m_AfferentBaroreceptor + wSN_AC * m_ChemoreceptorFiringRate_Hz - thetaSN);
  double testSympatheticNode = std::exp(exponentSN);
  m_SympatheticHeartSignal = testSympatheticNode;
  ULIM(testSympatheticNode, fEsMax);
  //Sympathetic signal to peripheral vascular beds
  double wSV_AB = -1.13, wSV_AC = 1.716, wSV_AP = -0.34, wSV_AA = -1.0, wSV_AT = 1.0;
  double thetaSV = 0.0;
  double exponentSV = kES * (wSV_AB * m_AfferentBaroreceptor + wSV_AC * m_ChemoreceptorFiringRate_Hz + wSV_AP * m_AfferentPulmonary - thetaSV);
  double testSympatheticPeripheral = std::exp(exponentSV);
  m_SympatheticPeripheralSignal = testSympatheticPeripheral;
  ULIM(testSympatheticPeripheral, fEsMax);

  m_data.GetDataTrack().Probe("Ursino_SympatheticNode", testSympatheticNode);
  m_data.GetDataTrack().Probe("Ursino_SympatheticPeripheral", testSympatheticPeripheral);
  m_data.GetDataTrack().Probe("SympatheticNodeExponent", exponentSN);
  m_data.GetDataTrack().Probe("SympatheticPeripheralExponent", exponentSV);

  //Parasympathetic / vagal signal path (*EV)
  double fEvInf = 6.3;
  double fEv0 = 3.2;
  double kEV_Hz = 7.06;
  double wEV_AC = 0.2, wEV_AP = 0.103;
  
  double thetaEV_Hz = 0.0; //Derivation: -0.68 (Ursino) + .105 * 5 (thermal addition) + .105 * 3.5 (pulmonary offset for BG)
  double exponentEV = (m_AfferentBaroreceptor - 25.0) / kEV_Hz;
  double testVagal = 1.2 * std::exp(exponentEV) / (1.0 + std::exp(exponentEV));
  testVagal += (wEV_AC * m_ChemoreceptorFiringRate_Hz - wEV_AP * m_AfferentPulmonary - thetaEV_Hz);
  m_VagalSignal = testVagal;
  m_data.GetDataTrack().Probe("Ursino_Parasympathetic", testVagal);

}

void Nervous::EfferentResponse()
{
  //Heart Rate
  double hPB = 0.5;
  double hS = 0.3;
  double baselineSignal = -hPB * m_VagalSignal_Baseline + hS * m_SympatheticHeartSignal_Baseline;
  double HR0 = m_data.GetPatient().GetHeartRateBaseline(FrequencyUnit::Per_min) / (1.0 + baselineSignal);
  double HRNext = HR0 * (1.0 - hPB * m_VagalSignal + hS * m_SympatheticHeartSignal);
  m_data.GetDataTrack().Probe("Randall_HR", HRNext);

  //Heart elastance
  double leftElastance0 = 2.41;
  double rightElastance0 = 0.474;
  double leftGain = 0.45;
  double rightGain = 0.282;
  double leftTau = 2.0;
  double rightTau = 1.5;

  double dLeftElastance = (1.0 / leftTau) * (-m_HeartElastanceEffector + leftGain * m_SympatheticHeartSignal);
  double dRightElastance = (1.0 / rightTau) * (-m_ElastanceMod + rightGain * m_SympatheticHeartSignal);
  m_HeartElastanceEffector += (dLeftElastance * m_dt_s);
  m_ElastanceMod += (dRightElastance * m_dt_s);

  m_data.GetDataTrack().Probe("LeftHeartElastance_Mod", m_HeartElastanceEffector);
  m_data.GetDataTrack().Probe("LeftHeartElastance_Next", m_HeartElastanceEffector + leftElastance0);
  m_data.GetDataTrack().Probe("RightHeartElastance_Mod", m_ElastanceMod);
  m_data.GetDataTrack().Probe("RightHeartElastance_Next", m_ElastanceMod + rightElastance0);

  //Resistance
  double tauR_s = 3.0;
  double gainR = 0.6;
  double baseR_mmHg_s_Per_mL = 1.0 - gainR * m_SympatheticPeripheralSignal_Baseline;
  
  double dR = (1.0 / tauR_s) * (-m_ResistanceMod + gainR * m_SympatheticPeripheralSignal);
  m_ResistanceMod += (dR * m_dt_s);
  m_data.GetDataTrack().Probe("Resistance_Mod", m_ResistanceMod);
  m_data.GetDataTrack().Probe("Resistance_Next", baseR_mmHg_s_Per_mL + m_ResistanceMod);

  //Venous Compliance
  double baselineVolume = 3200.0;
  double vol0 = 3248.1;
  double gainVolume = -275.0;
  double tauVolume = 10.0;

  double dVolume = (1.0 / tauVolume) * (-m_ComplianceMod + gainVolume * m_SympatheticPeripheralSignal);
  m_ComplianceMod += (dVolume * m_dt_s);

  m_data.GetDataTrack().Probe("Volume_Mod", m_ComplianceMod);
  m_data.GetDataTrack().Probe("Compliance_Fraction", (vol0 + m_ComplianceMod) / baselineVolume);

  if (m_FeedbackActive) {
    GetBaroreceptorHeartRateScale().SetValue(HRNext / m_data.GetPatient().GetHeartRateBaseline(FrequencyUnit::Per_min));
    GetBaroreceptorHeartElastanceScale().SetValue((m_HeartElastanceEffector + leftElastance0) / 2.49);
    GetResistanceScaleExtrasplanchnic().SetValue(baseR_mmHg_s_Per_mL + m_ResistanceMod);
    GetResistanceScaleMuscle().SetValue(baseR_mmHg_s_Per_mL + m_ResistanceMod);
    GetResistanceScaleSplanchnic().SetValue(baseR_mmHg_s_Per_mL + m_ResistanceMod);
    GetBaroreceptorComplianceScale().SetValue((vol0 + m_ComplianceMod) / baselineVolume);
  }

}

void Nervous::CerebralAutoregulation()
{
  double tauAuto = 40.0;
  double complianceGainHigh = 2.87, complianceGainLow = 0.11;
  double Kr = 4.96e4;
  double complianceMid = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralArteries1ToSpinalFluid)->GetComplianceBaseline(FlowComplianceUnit::mL_Per_mmHg);
  double complianceGain = 0.0;
  double complianceSlope = 0.0;
  double filterConstant = 0.5;

  double cranialCO2_mmHg = m_data.GetCompartments().GetExtracellularFluid(*m_data.GetCompartments().GetTissueCompartment(BGE::TissueCompartment::Brain)).GetSubstanceQuantity(m_data.GetSubstances().GetCO2())->GetPartialPressure(PressureUnit::mmHg);
  double cbfBaseDelta = 1.8 / (1.0 + std::exp(-0.06 * (cranialCO2_mmHg - 57.0))) - 0.6;
  cbfBaseDelta = 0.0;
  double cbfBase = m_CerebralBloodFlowBaseline_mL_Per_s * (1.0 + cbfBaseDelta);
  double kAuto = 18.0 / (1.0 + std::exp(2.0 * (cranialCO2_mmHg - m_CerebralCO2Baseline_mmHg) / m_CerebralCO2Baseline_mmHg));
  double cerebralBloodFlow_mL_Per_s = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralCapillariesToCerebralVeins1)->GetFlow(VolumePerTimeUnit::mL_Per_s);
  double dCerebralBloodFlow = filterConstant * (-m_CerebralBloodFlowFilter + cerebralBloodFlow_mL_Per_s);
  m_CerebralBloodFlowFilter += (dCerebralBloodFlow * m_dt_s);

  double dAuto = (1.0 / tauAuto) * (-m_CerebralAutoEffect + kAuto * (m_CerebralBloodFlowFilter - 12.0) / 12.0);
  m_CerebralAutoEffect += (dAuto * m_dt_s);

  if (m_CerebralAutoEffect < 0.0) {
    complianceGain = complianceGainHigh;
    complianceSlope = complianceGainHigh / 4.0;
  } else {
    complianceGain = complianceGainLow;
    complianceSlope = complianceGainLow / 4.0;
  }

  double nextCompliance = ((complianceMid - 0.5 * complianceGain) + (complianceMid + 0.5 * complianceGain) * std::exp(-m_CerebralAutoEffect / complianceSlope)) / (1.0 + std::exp(-m_CerebralAutoEffect / complianceSlope));
  double cerebralArteryVolume_mL = m_data.GetCompartments().GetLiquidCompartment(BGE::VascularCompartment::CerebralArteries)->GetVolume(VolumeUnit::mL);
  double cerebralArteryVolumeBase = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetNode(BGE::CerebralNode::CerebralArteries1)->GetVolumeBaseline(VolumeUnit::mL);
  double cerebralResistanceBase = m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralArteries2ToCapillaries)->GetResistanceBaseline(FlowResistanceUnit::mmHg_s_Per_mL);
  Kr = cerebralResistanceBase / (std::pow(complianceMid / 12.5, 2.0));
  double nextResistance = Kr * std::pow(complianceMid / cerebralArteryVolume_mL, 2.0);

  if (m_FeedbackActive) {
    m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralArteries1ToSpinalFluid)->GetNextCompliance().SetValue(nextCompliance, FlowComplianceUnit::mL_Per_mmHg);
    m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath(BGE::CerebralPath::CerebralArteries2ToCapillaries)->GetNextResistance().SetValue(nextResistance, FlowResistanceUnit::mmHg_s_Per_mL);
  }

  m_data.GetDataTrack().Probe("Cerebral_CO2", cranialCO2_mmHg);
  m_data.GetDataTrack().Probe("Cerebral_AutoEffect", m_CerebralAutoEffect);
  m_data.GetDataTrack().Probe("Cerebral_Compliance", nextCompliance);
  m_data.GetDataTrack().Probe("Cerebral_Resistance", nextResistance);
  m_data.GetDataTrack().Probe("Cerebral_Volume", cerebralArteryVolume_mL);
  m_data.GetDataTrack().Probe("Cerebral_FilterFlow", m_CerebralBloodFlowFilter);
}
}