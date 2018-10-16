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
#include <biogears/cdm/CommonDataModel.h>
#include <biogears/cdm/engine/PhysiologyEngineConfiguration.h>
#include <biogears/exports.h>
#include <biogears/schema/biogears/OverrideConfig.hxx>
#include <biogears/schema/cdm/Properties.hxx>
#include <biogears/cdm/properties/SEScalarTemperature.h>
#include <biogears/cdm/properties/SEScalarPressure.h>
#include <biogears/cdm/properties/SEScalarFrequency.h>

CDM_BIND_DECL(OverrideConfigData)

namespace biogears {
class BIOGEARS_API OverrideConfig : public Loggable {
public:
  OverrideConfig();
  virtual ~OverrideConfig();

  virtual void Clear(); // Deletes all members

  // virtual void Merge(const OverrideConfig& from);

  virtual bool LoadOverride(const std::string& file);
  virtual bool Load(const CDM::OverrideConfigData& in);
  virtual CDM::OverrideConfigData* Unload() const;
  // bool ReadOverrideParmeters(const std::string& overrideParameterFile);
  //virtual bool HasOverride() const;
  //virtual SEScalarPressure& GetOverride();
  //virtual double GetOverride(const PressureUnit& unit) const;

  // Cardiovascular Parameters
  virtual bool HasMeanArterialPressureOverride() const;
  virtual SEScalarPressure& GetMeanArterialPressureOverride();
  virtual double GetMeanArterialPressureOverride(const PressureUnit& unit) const;
  virtual bool HasHeartRateOverride() const;
  virtual SEScalarFrequency& GetHeartRateOverride();
  virtual double GetHeartRateOverride(const FrequencyUnit& unit) const;

  virtual bool HasEnableCardiovascularOverride() const
  {
    return m_overrideMode != (CDM::enumOnOff::Off);
  }
  virtual bool IsCardiovascularOverrideEnabled() const
  {
    return m_overrideMode == CDM::enumOnOff::On;
  }
  virtual void EnableCardiovascularOverride(CDM::enumOnOff::value s) { m_overrideMode = s; }

  //Energy Parameters
  virtual bool HasCoreTemperatureOverride() const;
  virtual SEScalarTemperature& GetCoreTemperatureOverride();
  virtual double GetCoreTemperatureOverride(const TemperatureUnit& unit) const;
  virtual bool HasSkinTemperatureOverride() const;
  virtual SEScalarTemperature& GetSkinTemperatureOverride();
  virtual double GetSkinTemperatureOverride(const TemperatureUnit& unit) const;

  virtual bool HasEnableEnergyOverride() const
  {
    return m_overrideMode != (CDM::enumOnOff::Off);
  }
  virtual bool IsEnergyOverrideEnabled() const
  {
    return m_overrideMode == CDM::enumOnOff::On;
  }
  virtual void EnableEnergyOverride(CDM::enumOnOff::value s) { m_overrideMode = s; }

protected:
  virtual void Unload(CDM::OverrideConfigData& data) const;
  bool ReadOverrideParameters(const std::string& overrideParameterFile);

  CDM::enumOnOff::value m_overrideMode;
  SEScalarPressure* m_MeanArterialPressureOverride;
  SEScalarTemperature* m_CoreTemperatureOverride;
  SEScalarTemperature* m_SkinTemperatureOverride;
  SEScalarFrequency* m_HeartRateOverride;

};
}
