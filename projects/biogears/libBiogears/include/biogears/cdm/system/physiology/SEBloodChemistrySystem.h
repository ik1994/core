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

#pragma once
#include <biogears/exports.h>
#include <biogears/cdm/patient/actions/SEInflammationState.h>
#include <biogears/cdm/system/SESystem.h>
#include <biogears/schema/biogears/BioGearsPhysiology.hxx>


namespace biogears {

class SEScalarMass;
class MassUnit;
class SEScalarMassPerVolume;
class MassPerVolumeUnit;
class SEScalarHeatCapacitancePerMass;
class HeatCapacitancePerMassUnit;
class SEScalarFraction;
class SEScalarVolume;
class VolumeUnit;
class SEScalarAmountPerVolume;
class AmountPerVolumeUnit;
class SEScalarPressure;
class PressureUnit;
class SESepsisState;
class SEInflammationState;


/** @copydoc Physiology_BloodChemistrySystemData
  @nosubgrouping */
class BIOGEARS_API SEBloodChemistrySystem : public SESystem {
public:
  SEBloodChemistrySystem(Logger* logger);
  ~SEBloodChemistrySystem() override;

  static size_t TypeHash() { return reinterpret_cast<size_t>(&TypeHash); }  //! Hopefully this returns a unique ID for every type
  static constexpr char const * const  TypeTag() { return "SEBloodChemistrySystem"; }
  const char* classname() const override { return TypeTag(); }
  size_t hash_code() const override { return TypeHash(); }

  void Clear() override; /**< @copydoc DOXY_CDM_CLEAR */

  const SEScalar* GetScalar(const char* name) override;
  const SEScalar* GetScalar(const std::string& name) override; /**< @copydoc DOXY_CDM_GET_SCALAR */

  /**  @name Serialization */ //@{
  bool Load(const CDM::BloodChemistrySystemData& in); /**< @copydoc DOXY_CDM_LOAD */
  CDM::BloodChemistrySystemData* Unload() const override; /**< @copydoc DOXY_CDM_UNLOAD */

  Tree<const char*> GetPhysiologyRequestGraph() const override;
protected:
  void Unload(CDM::BloodChemistrySystemData& data) const; /**< @copydoc DOXY_CDM_UNLOAD_TO */ //@}
public:
  /**  @name BloodDensity */ //@{ @copybrief Physiology_BloodChemistrySystemData_BloodDensity
  bool HasBloodDensity() const; /**< @copydoc DOXY_CDM_HAS */
  SEScalarMassPerVolume& GetBloodDensity(); /**< @copydoc DOXY_CDM_GET */
  double GetBloodDensity(const MassPerVolumeUnit& unit) const; /**< @copydoc DOXY_CDM_GET_VALUE */ //@}

 
  /** @name ArterialBloodPH
  *  @brief @copybrief Physiology_BloodChemistrySystemData_ArterialBloodPH
  *  @{*/
  /// %Test if member has been allocated
  bool HasArterialBloodPH() const;
  /// Get member class, allocate if nullptr
  SEScalar& GetArterialBloodPH();
  double GetArterialBloodPH() const;
  //@}

  /** @name ArterialBloodPHBaseline
  *  @brief @copybrief Physiology_BloodChemistrySystemData_ArterialBloodPHBaseline
  *  @{*/
  /// %Test if member has been allocated
  bool HasArterialBloodPHBaseline() const;
  /// Get member class, allocate if nullptr
  SEScalar& GetArterialBloodPHBaseline();
  double GetArterialBloodPHBaseline() const;
  //@}

  /** @name VenousBloodPH
  *  @brief @copybrief Physiology_BloodChemistrySystemData_VenousBloodPH
  *  @{*/
  /// %Test if member has been allocated
  bool HasVenousBloodPH() const;
  /// Get member class, allocate if nullptr
  SEScalar& GetVenousBloodPH();
  double GetVenousBloodPH() const;
  //@}

  /** @name BloodSpecificHeat
  *  @brief @copybrief Physiology_BloodChemistrySystemData_BloodSpecificHeat
  *  @{*/
  /// %Test if member has been allocated
  bool HasBloodSpecificHeat() const;
  /// Get member class, allocate if nullptr
  SEScalarHeatCapacitancePerMass& GetBloodSpecificHeat();
  double GetBloodSpecificHeat(const HeatCapacitancePerMassUnit& unit) const;
  //@}

  /** @name BloodUreaNitrogenConcentration
  *  @brief @copybrief Physiology_BloodChemistrySystemData_BloodUreaNitrogenConcentration
  *  @{*/
  /// %Test if member has been allocated
  bool HasBloodUreaNitrogenConcentration() const;
  /// Get member class, allocate if nullptr
  SEScalarMassPerVolume& GetBloodUreaNitrogenConcentration();
  double GetBloodUreaNitrogenConcentration(const MassPerVolumeUnit& unit) const;
  //@}

  /** @name CarbonDioxideSaturation
  *  @brief @copybrief Physiology_BloodChemistrySystemData_CarbonDioxideSaturation
  *  @{*/
  /// %Test if member has been allocated
  bool HasCarbonDioxideSaturation() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetCarbonDioxideSaturation();
  double GetCarbonDioxideSaturation() const;
  //@}

  /** @name CarbonMonoxideSaturation
  *  @brief @copybrief Physiology_BloodChemistrySystemData_CarbonMonoxideSaturation
  *  @{*/
  /// %Test if member has been allocated
  bool HasCarbonMonoxideSaturation() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetCarbonMonoxideSaturation();
  double GetCarbonMonoxideSaturation() const;
  //@}

  /** @name Hematocrit
  *  @brief @copybrief Physiology_BloodChemistrySystemData_Hematocrit
  *  @{*/
  /// %Test if member has been allocated
  bool HasHematocrit() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetHematocrit();
  double GetHematocrit() const;
  //@}

  /** @name HemoglobinContent
  *  @brief @copybrief Physiology_BloodChemistrySystemData_HemoglobinContent
  *  @{*/
  /// %Test if member has been allocated
  bool HasHemoglobinContent() const;
  /// Get member class, allocate if nullptr
  SEScalarMass& GetHemoglobinContent();
  double GetHemoglobinContent(const MassUnit& unit) const;
  //@}

  /** @name OxygenSaturation
   *  @brief @copybrief Physiology_BloodChemistrySystemData_OxygenSaturation
   *  @{*/
  /// %Test if member has been allocated
  bool HasOxygenSaturation() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetOxygenSaturation();
  double GetOxygenSaturation() const;
  //@}

  /** @name OxygenSaturation
   *  @brief @copybrief Physiology_BloodChemistrySystemData_OxygenVenousSaturation
   *  @{*/
  /// %Test if member has been allocated
  bool HasOxygenVenousSaturation() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetOxygenVenousSaturation();
  double GetOxygenVenousSaturation() const;
  //@}

  /** @name Phosphate
  *  @brief @copybrief Physiology_BloodChemistrySystemData_Phosphate
  *  @{*/
  /// %Test if member has been allocated
  bool HasPhosphate() const;
  /// Get member class, allocate if nullptr
  SEScalarAmountPerVolume& GetPhosphate();
  double GetPhosphate(const AmountPerVolumeUnit& unit) const;
  //@}

  /** @name PlasmaVolume
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PlasmaVolume
  *  @{*/
  /// %Test if member has been allocated
  bool HasPlasmaVolume() const;
  /// Get member class, allocate if nullptr
  SEScalarVolume& GetPlasmaVolume();
  double GetPlasmaVolume(const VolumeUnit& unit) const;
  //@}

  /** @name PulseOximetry
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PulseOximetry
  *  @{*/
  /// %Test if member has been allocated
  bool HasPulseOximetry() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetPulseOximetry();
  double GetPulseOximetry() const;
  //@}

  /** @name RedBloodCellAcetylcholinesteraseConcentration
  *  @brief @copybrief Physiology_BloodChemistrySystemData_RedBloodCellAcetylcholinesteraseConcentration
  *  @{*/
  /// %Test if member has been allocated
  bool HasRedBloodCellAcetylcholinesterase() const;
  /// Get member class, allocate if nullptr
  SEScalarAmountPerVolume& GetRedBloodCellAcetylcholinesterase();
  double GetRedBloodCellAcetylcholinesterase(const AmountPerVolumeUnit& unit) const;
  //@}

  /** @name RedBloodCellCount
   *  @brief @copybrief Physiology_BloodChemistrySystemData_RedBloodCellCount
   *  @{*/
  /// %Test if member has been allocated
  bool HasRedBloodCellCount() const;
  /// Get member class, allocate if nullptr
  SEScalarAmountPerVolume& GetRedBloodCellCount();
  double GetRedBloodCellCount(const AmountPerVolumeUnit& unit) const;
  //@}

  /** @name StrongIonDifference
  *  @brief @copybrief Physiology_BloodChemistrySystemData_StrongIonDifference
  *  @{*/
  /// %Test if member has been allocated
  bool HasStrongIonDifference() const;
  /// Get member class, allocate if nullptr
  SEScalarAmountPerVolume& GetStrongIonDifference();
  double GetStrongIonDifference(const AmountPerVolumeUnit& unit) const;
  //@}

  /** @name ShuntFraction
   *  @brief @copybrief Physiology_BloodChemistrySystemData_ShuntFraction
   *  @{*/
  /// %Test if member has been allocated
  bool HasShuntFraction() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetShuntFraction();
  double GetShuntFraction() const;
  //@}

  bool HasTotalBilirubin() const;
  SEScalarMassPerVolume& GetTotalBilirubin();
  double GetTotalBilirubin(const MassPerVolumeUnit& unit) const;

  /** @name TotalProteinConcentration
  *  @brief @copybrief Physiology_BloodChemistrySystemData_TotalProteinConcentration
  *  @{*/
  /// %Test if member has been allocated
  bool HasTotalProteinConcentration() const;
  /// Get member class, allocate if nullptr
  SEScalarMassPerVolume& GetTotalProteinConcentration();
  double GetTotalProteinConcentration(const MassPerVolumeUnit& unit) const;
  //@}

  /** @name VolumeFractionNeutralPhospholipidInPlasma
  *  @brief @copybrief Physiology_BloodChemistrySystemData_VolumeFractionNeutralPhospholipidInPlasma
  *  @{*/
  /// %Test if member has been allocated
  bool HasVolumeFractionNeutralPhospholipidInPlasma() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetVolumeFractionNeutralPhospholipidInPlasma();
  double GetVolumeFractionNeutralPhospholipidInPlasma() const;
  //@}

  /** @name VolumeFractionNeutralLipidInPlasma
  *  @brief @copybrief Physiology_BloodChemistrySystemData_VolumeFractionNeutralLipidInPlasma
  *  @{*/
  /// %Test if member has been allocated
  bool HasVolumeFractionNeutralLipidInPlasma() const;
  /// Get member class, allocate if nullptr
  SEScalarFraction& GetVolumeFractionNeutralLipidInPlasma();
  double GetVolumeFractionNeutralLipidInPlasma() const;
  //@}

  /** @name ArterialCarbonDioxidePressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_ArterialCarbonDioxidePressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasArterialCarbonDioxidePressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetArterialCarbonDioxidePressure();
  double GetArterialCarbonDioxidePressure(const PressureUnit& unit) const;
  //@}

  /** @name ArterialOxygenPressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_ArterialOxygenPressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasArterialOxygenPressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetArterialOxygenPressure();
  double GetArterialOxygenPressure(const PressureUnit& unit) const;
  //@}

  /** @name PulmonaryArterialCarbonDioxidePressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PulmonaryArterialCarbonDioxidePressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasPulmonaryArterialCarbonDioxidePressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetPulmonaryArterialCarbonDioxidePressure();
  double GetPulmonaryArterialCarbonDioxidePressure(const PressureUnit& unit) const;
  //@}

  /** @name PulmonaryArterialOxygenPressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PulmonaryArterialOxygenPressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasPulmonaryArterialOxygenPressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetPulmonaryArterialOxygenPressure();
  double GetPulmonaryArterialOxygenPressure(const PressureUnit& unit) const;
  //@}

  /** @name PulmonaryVenousCarbonDioxidePressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PulmonaryVenousCarbonDioxidePressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasPulmonaryVenousCarbonDioxidePressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetPulmonaryVenousCarbonDioxidePressure();
  double GetPulmonaryVenousCarbonDioxidePressure(const PressureUnit& unit) const;
  //@}

  /** @name PulmonaryVenousOxygenPressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_PulmonaryVenousOxygenPressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasPulmonaryVenousOxygenPressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetPulmonaryVenousOxygenPressure();
  double GetPulmonaryVenousOxygenPressure(const PressureUnit& unit) const;
  //@}

  /** @name VenousCarbonDioxidePressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_VenousCarbonDioxidePressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasVenousCarbonDioxidePressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetVenousCarbonDioxidePressure();
  double GetVenousCarbonDioxidePressure(const PressureUnit& unit) const;
  //@}
  /** @name VenousOxygenPressure
  *  @brief @copybrief Physiology_BloodChemistrySystemData_VenousOxygenPressure
  *  @{*/
  /// %Test if member has been allocated
  bool HasVenousOxygenPressure() const;
  /// Get member class, allocate if nullptr
  SEScalarPressure& GetVenousOxygenPressure();
  double GetVenousOxygenPressure(const PressureUnit& unit) const;
  //@}
  virtual bool HasAcuteInflammatoryResponse() const;
  virtual SEInflammationState& GetAcuteInflammatoryResponse();
protected:
  SEScalarMassPerVolume* m_BloodDensity;
  SEScalar* m_ArterialBloodPH;
  SEScalar* m_ArterialBloodPHBaseline;
  SEScalar* m_VenousBloodPH;
  SEScalarHeatCapacitancePerMass* m_BloodSpecificHeat;
  SEScalarMassPerVolume* m_BloodUreaNitrogenConcentration;
  SEScalarFraction* m_CarbonDioxideSaturation;
  SEScalarFraction* m_CarbonMonoxideSaturation;
  SEScalarFraction* m_Hematocrit;
  SEScalarMass* m_HemoglobinContent;
  SEScalarFraction* m_OxygenSaturation;
  SEScalarFraction* m_OxygenVenousSaturation;
  SEScalarAmountPerVolume* m_Phosphate;
  SEScalarVolume* m_PlasmaVolume;
  SEScalarFraction* m_PulseOximetry;
  SEScalarAmountPerVolume* m_RedBloodCellAcetylcholinesterase;
  SEScalarAmountPerVolume* m_RedBloodCellCount;
  SEScalarFraction* m_ShuntFraction;
  SEScalarAmountPerVolume* m_StrongIonDifference;
  SEScalarMassPerVolume* m_TotalBilirubin;
  SEScalarMassPerVolume* m_TotalProteinConcentration;
  SEScalarFraction* m_VolumeFractionNeutralPhospholipidInPlasma;
  SEScalarFraction* m_VolumeFractionNeutralLipidInPlasma;

  SEScalarPressure* m_ArterialCarbonDioxidePressure;
  SEScalarPressure* m_ArterialOxygenPressure;
  SEScalarPressure* m_PulmonaryArterialCarbonDioxidePressure;
  SEScalarPressure* m_PulmonaryArterialOxygenPressure;
  SEScalarPressure* m_PulmonaryVenousCarbonDioxidePressure;
  SEScalarPressure* m_PulmonaryVenousOxygenPressure;
  SEScalarPressure* m_VenousCarbonDioxidePressure;
  SEScalarPressure* m_VenousOxygenPressure;
  SEInflammationState* m_AcuteInflammatoryResponse;
};
}

