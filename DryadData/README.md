# Commodity canola and seed canola visitation and plant data
---

Three csv files containing:

- Insect visitation at plots of hybrid seed canola flowers, as well as pollen deposition (seedVisitation.csv)
- Seed yield data from commodity (commodityPlants.csv) and hybrid seed (seedPlants.csv) plants

## Description of the data and file structure

Column descriptions for seedVisitation.csv

- Field: unique ID for field
- Distance: distance (m) from honey bee (A. mellifera) hives at edge of planted crop
- minDist: distance (m) from nearest alfalfa leafcutting bee (M. rotundata) shelter
- Bay: was observation from hermaphroditic ("male": M) or male-sterile ("female": F) bay?
- EdgeCent: was observation from the edge of the female bay or center?
- Date: sampling date in YYYY-MM-DD format
- StartTime: start of observation in HH:MM format (24 hour time)
- TotalTime: time taken for observations (minutes)
- FlDens: flower count per m2
- lbee: observed M. rotundata visits. Individual bees were not counted, only total visits
- hbee: observed A. mellifera visits. Individual bees were not counted, only total visits
- otherBee: observed visits from other bee species
- hFly: observed visits from hoverfly (Diptera, family Syrphidae)
- AirTemp: air temperature (C) at the beginning of the survey
- WindSp: wind speed (km/hr) at the beginning of the survey
- RH: relative humidity (%)
- Year: year of sampling
- PlDens: plant stem count per m2 (number of individual plants)
- Pollen1 - Pollen5: Brassica pollen counts on each stigma sampled from that plot

Column descriptions for commodityPlants.csv

- Year: year of sampling
- Field: unique ID for field
- Distance: distance (m) from honey bee (A. mellifera) hives at edge of planted crop
- Plant: unique ID for plant
- VegMass: total dried vegetative mass (no seeds) of plant
- SeedMass: total dried seed mass
- Branch: number of primary branches coming off of the main stem
- Pods: number of pods per plant
- Missing: number of pedicels without a corresponding pod (i.e. aborted flowers)
- SeedCount: total number of seeds per plant - counted using a seed counting machine
- PodCount1 - PodCount5: Number of seeds per pod for 5 randomly sampled pods
- PodMass1 - PodMass5: Weight of seeds for 5 randomly sampled pods (using the same pods as above)

Column descriptions for seedPlants.csv

- Year: year of sampling
- Field: unique ID for field
- Distance: distance (m) from honey bee (A. mellifera) hives at edge of planted crop
- EdgeCent: was plant from the edge of the female bay or center?
- Plant: unique ID for plant
- VegMass: total dried vegetative mass (no seeds) of plant
- SeedMass: total dried seed mass
- Branch: number of primary branches coming off of the main stem
- Pods: number of pods per plant
- Missing: number of pedicels without a corresponding pod (i.e. aborted flowers)
- PodCount1 - PodCount5: Number of seeds per pod for 5 randomly sampled pods
- PodMass1 - PodMass5: Weight of seeds for 5 randomly sampled pods (using the same pods as above)

## Sharing/Access information

NA

## Code/Software

NA