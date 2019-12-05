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
  m_blockActive = false;

  m_BaroreceptorFatigueScale = 0.0;
  m_ChemoreceptorFiringRate_Hz = 2.0;
  m_ChemoreceptorFiringRateSetPoint_Hz = m_ChemoreceptorFiringRate_Hz;
  m_CentralFrequencyDelta_Per_min = 0.0;
  m_CentralPressureDelta_cmH2O = 0.0;
  m_PeripheralFrequencyDelta_Per_min = 0.0;
  m_PeripheralPressureDelta_cmH2O = 0.0;
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
  m_ChemoreceptorFiringRate_Hz = in.ChemoreceptorFiringRate_Hz();
  m_ChemoreceptorFiringRateSetPoint_Hz = in.ChemoreceptorFiringRateSetPoint_Hz();
  m_PeripheralBloodGasInteractionBaseline_Hz = in.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz();
  m_CentralFrequencyDelta_Per_min = in.CentralFrequencyDelta_Per_min();
  m_CentralPressureDelta_cmH2O = in.CentralPressureDelta_cmH2O(); //Meant to be cmH2O, correct later
  m_PeripheralFrequencyDelta_Per_min = in.PeripheralFrequencyDelta_Per_min();
  m_PeripheralPressureDelta_cmH2O = in.PeripheralPressureDelta_cmH2O(); //Meant to be cmH2O, correct later

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
  data.ChemoreceptorPeripheralBloodGasInteractionBaseline_Hz(m_PeripheralBloodGasInteractionBaseline_Hz);
  data.ChemoreceptorFiringRate_Hz(m_ChemoreceptorFiringRate_Hz);
  data.ChemoreceptorFiringRateSetPoint_Hz(m_ChemoreceptorFiringRateSetPoint_Hz);
  data.CentralFrequencyDelta_Per_min(m_CentralFrequencyDelta_Per_min);
  data.CentralPressureDelta_cmH2O(m_CentralPressureDelta_cmH2O);
  data.PeripheralFrequencyDelta_Per_min(m_PeripheralFrequencyDelta_Per_min);
  data.PeripheralPressureDelta_cmH2O(m_PeripheralPressureDelta_cmH2O);
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
  m_PreviousCO2Concentration = 0.45;
  m_DynamicCO2Signal = 0.0;
}

void Nervous::AtSteadyState()
{
  if (m_data.GetState() == EngineState::AtInitialStableState)
    m_FeedbackActive = true;

  // The set-points (Baselines) get reset at the end of each stabilization period.
  m_ArterialOxygenSetPoint_mmHg = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg);
  m_ArterialCarbonDioxideSetPoint_mmHg = m_data.GetBloodChemistry().GetArterialCarbonDioxidePressure(PressureUnit::mmHg);

  //Central and peripheral ventilation changes are set to 0 because patient baseline ventilation is updated to include
  //their contributions at steady state.

  //The chemoreceptor firing rate and its setpoint are reset so that central and peripheral derivatives will evaluate to 0
  //the first time step after stabilization (and will stay that way, assuming no other perturbations to blood gas levels)
  m_ChemoreceptorFiringRateSetPoint_Hz = m_ChemoreceptorFiringRate_Hz;

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
  BaroreceptorFeedback();
  ChemoreceptorFeedback();
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
  if (!m_FeedbackActive) {
    return;
  }

  //First calculate the sympathetic and parasympathetic firing rates:
  double nu = m_data.GetConfiguration().GetResponseSlope();
  double meanArterialPressure_mmHg = m_data.GetCardiovascular().GetMeanArterialPressure(PressureUnit::mmHg);
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
  normalizedHeartRate += deltaNormalizedHeartRate;
  GetBaroreceptorHeartRateScale().SetValue(normalizedHeartRate);

  //Calculate the normalized change in heart elastance
  double normalizedHeartElastance = GetBaroreceptorHeartElastanceScale().GetValue();
  double tauElastance_s = m_data.GetConfiguration().GetHeartElastanceDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedHeartElastance = (1.0 / tauElastance_s) * (-normalizedHeartElastance + m_normalizedAlphaElastance * sympatheticFraction + m_normalizedGammaElastance) * m_dt_s;
  normalizedHeartElastance += deltaNormalizedHeartElastance;
  GetBaroreceptorHeartElastanceScale().SetValue(normalizedHeartElastance);


  //Calculate the normalized change in flow resistance for any cardiovascular resistor
  double normalizedResistance = GetBaroreceptorResistanceScale().GetValue();
  double tauResistance_s = m_data.GetConfiguration().GetSystemicResistanceDistributedTimeDelay(TimeUnit::s);
  double deltaNormalizedResistance = (1.0 / tauResistance_s) * (-normalizedResistance + (m_normalizedAlphaResistance - m_BaroreceptorFatigueScale) * sympatheticFraction + m_normalizedGammaResistance) * m_dt_s;
  normalizedResistance += (deltaNormalizedResistance);
  GetBaroreceptorResistanceScale().SetValue(normalizedResistance);

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
    traumaPain *= 20.0;   //25% TBSA burn will give pain scale = 5, 40% TBSA will give pain scale = 8.0
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
  double drugCNSModifier = m_data.GetDrugs().GetCentralNervousResponse().GetValue();
  //-------Respiratory Feedback:  This is active throughtout the simulation (including stabilization)------------------------------

  //Chemoreceptor time constants (all in s):  _p/c_P/F --> peripheral(p)/central(c) chemoreceptor effect on respiratory driver pressure(P)/frequency(F)
  double tau_p_P = 100.0, tau_p_F = 110.0, tau_c_P = 180.0, tau_c_F = 180.;
  //Chemoreceptor gains -- Tuned to data from Reynold, 1972 and Reynolds, 1973 (cited in Cheng, 2016)
  double gain_p_P = 0.8, gain_p_F = 0.50, gain_c_P = 0.75, gain_c_F = 0.70;

  //Inputs
  double arterialO2Pressure_mmHg = m_data.GetBloodChemistry().GetArterialOxygenPressure(PressureUnit::mmHg);
  double arterialCO2Pressure_mmHg = m_data.GetBloodChemistry().GetArterialCarbonDioxidePressure(PressureUnit::mmHg);

  //Magosso and Ursino cite findings that central chemoreceptors are less sensitive at sub-normal levels of CO2 than to super-normal levels
  if (arterialCO2Pressure_mmHg < m_ArterialCarbonDioxideSetPoint_mmHg && m_data.GetState() == EngineState::Active) {
    gain_c_P *= 0.067;
    gain_c_F *= 0.067;
  }

  //Dissociation curves--Derive apparent O2 and CO2 saturations from arterial pressures (sat calculator not good enough to use here)
  double a1 = 0.3836, alpha1 = 3.198e-2, K1 = 15.0, beta1 = 8.275e-3, C1 = 9.0, Z = 2.277e-2, Ca_O2max = 0.2;		//O2 params
  double a2 = 1.819, alpha2 = 5.591e-2, K2 = 194.4, beta2 = 3.255e-2, C2 = 86.11;									//CO2 params
  double bohrHaldaneO2 = arterialO2Pressure_mmHg * (1.0 + beta1 * arterialCO2Pressure_mmHg) / (K1 * (1.0 + alpha1 * arterialCO2Pressure_mmHg));
  double bohrHaldaneCO2 = arterialCO2Pressure_mmHg * (1.0 + beta2 * arterialO2Pressure_mmHg) / (K2 * (1.0 + alpha2 * arterialO2Pressure_mmHg));
  double o2Concentration = Z * C1 * std::pow(bohrHaldaneO2, 1.0 / a1) / (1.0 + std::pow(bohrHaldaneO2, 1.0 / a1));
  double co2Concentration = Z * C2 * std::pow(bohrHaldaneCO2, 1.0 / a2) / (1.0 + std::pow(bohrHaldaneCO2, 1.0 / a2));
  double o2Saturation = o2Concentration / Ca_O2max;
  //Peripheral (carotid) chemoreceptors -- afferent signal generation
  double A = 600.0, B = 10.18, KO2 = 200, KCO2 = 1.0, Ct = 0.36, Kstat = 20.0, Kdyn = 45.0;
  double tauPH = 3.5, tauZH = 600.0, tauPL = 3.5, dynWeight = 0.25;
  double xO2 = A * (1.0 - o2Saturation) + B;
  double psiO2 = KO2 * (1.0 - std::exp(-xO2 / KO2));
  double psiCO2 = KCO2 * (co2Concentration - Ct);
  double psiCombined = psiO2 * psiCO2;
  double staticSignal = 0.0;
  if (co2Concentration > Ct) {
    staticSignal = Kstat * (1.0 - exp(-psiCombined / Kstat));
  }
  double dCO2Concentration = dynWeight * (co2Concentration - m_PreviousCO2Concentration);
  m_PreviousCO2Concentration += dCO2Concentration * m_dt_s;
  double dDynamicCO2Signal = (1.0 / tauPH) * (tauZH * dCO2Concentration - m_DynamicCO2Signal);
  m_DynamicCO2Signal += dDynamicCO2Signal * m_dt_s;
  double afferentFiringInput = staticSignal + Kdyn * (1.0 - std::exp(-m_DynamicCO2Signal / Kdyn)); 
  double dAfferentFiring = (1.0 / tauPL) * (afferentFiringInput - m_ChemoreceptorFiringRate_Hz);

  //Update peripheral and central inputs to respiratory driver frequency and pressure
  double centralInput = arterialCO2Pressure_mmHg - m_ArterialCarbonDioxideSetPoint_mmHg;
  double peripheralInput = m_ChemoreceptorFiringRate_Hz - m_ChemoreceptorFiringRateSetPoint_Hz;
  double dFrequency_c = (-m_CentralFrequencyDelta_Per_min + gain_c_F * centralInput) / tau_c_F * m_dt_s;
  double dPressure_c = (-m_CentralPressureDelta_cmH2O + gain_c_P * centralInput) / tau_c_P * m_dt_s;
  double dFrequency_p = (-m_PeripheralFrequencyDelta_Per_min + gain_p_F * peripheralInput) / tau_p_F * m_dt_s;
  double dPressure_p = (-m_PeripheralPressureDelta_cmH2O + gain_p_P * peripheralInput) / tau_p_P * m_dt_s;

  m_data.GetDataTrack().Probe("O2Sat_Calc", o2Saturation);
  m_data.GetDataTrack().Probe("Chemo_StaticSignal", staticSignal);
  m_data.GetDataTrack().Probe("Chemo_DynamicSignal", m_DynamicCO2Signal);
  m_data.GetDataTrack().Probe("CO2_Concentration", m_PreviousCO2Concentration);
  m_data.GetDataTrack().Probe("Chem_PeripheralSignal", m_ChemoreceptorFiringRate_Hz);
  m_data.GetDataTrack().Probe("CentralFrequencyDelta", m_CentralFrequencyDelta_Per_min);
  m_data.GetDataTrack().Probe("CentralPressureDelta", m_CentralPressureDelta_cmH2O);
  m_data.GetDataTrack().Probe("PeripheralFrequencyDelta", m_PeripheralFrequencyDelta_Per_min);
  m_data.GetDataTrack().Probe("PeripheralPressureDelta", m_PeripheralPressureDelta_cmH2O);

  //Update Respiratory metrics
  double baselineRespirationRate_Per_min = m_data.GetPatient().GetRespirationRateBaseline(FrequencyUnit::Per_min);
  double baselineDrivePressure_cmH2O = m_Patient->GetRespiratoryDriverAmplitudeBaseline(PressureUnit::cmH2O);
  double nextRespirationRate_Per_min;
  double nextDrivePressure_cmH2O;

  if (m_data.GetState() < EngineState::AtSecondaryStableState) {
    nextDrivePressure_cmH2O = baselineDrivePressure_cmH2O - dPressure_c - dPressure_p;
    m_data.GetRespiratory().GetRespirationDriverPressure().SetValue(nextDrivePressure_cmH2O, PressureUnit::cmH2O);
    m_Patient->GetRespiratoryDriverAmplitudeBaseline().SetValue(nextDrivePressure_cmH2O, PressureUnit::cmH2O);
  } else {
    m_ChemoreceptorFiringRate_Hz += dAfferentFiring * m_dt_s;
    m_ChemoreceptorFiringRate_Hz = std::max(0.0, m_ChemoreceptorFiringRate_Hz);
    m_CentralFrequencyDelta_Per_min += dFrequency_c;
    m_CentralPressureDelta_cmH2O += dPressure_c;
    m_PeripheralFrequencyDelta_Per_min += dFrequency_p;
    m_PeripheralPressureDelta_cmH2O += dPressure_p;

    nextRespirationRate_Per_min = baselineRespirationRate_Per_min + m_CentralFrequencyDelta_Per_min + m_PeripheralFrequencyDelta_Per_min;
    nextDrivePressure_cmH2O = baselineDrivePressure_cmH2O - m_CentralPressureDelta_cmH2O - m_PeripheralPressureDelta_cmH2O;

    //This tuning of CNS modifier was done so that standard Fentanyl and Morphine runs approximately match main engine baselines.
    if (drugCNSModifier > ZERO_APPROX) {
      for (auto drug : m_data.GetSubstances().GetSubstances()) {
        if (drug->GetClassification() == CDM::enumSubstanceClass::Opioid) {
          nextRespirationRate_Per_min *= (1.0 - drugCNSModifier / 5.0);
          nextDrivePressure_cmH2O *= (1.0 - drugCNSModifier / 2.0);
          break; //Don't evaluate more than once if for some reason we give someone morphine and fentanyl
        }
      }
    }

    m_data.GetRespiratory().GetRespirationDriverFrequency().SetValue(nextRespirationRate_Per_min, FrequencyUnit::Per_min);
    m_data.GetRespiratory().GetRespirationDriverPressure().SetValue(nextDrivePressure_cmH2O, PressureUnit::cmH2O);
  }

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
}