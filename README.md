# Data from: Taxon- and size-specific variation in coral growth between tropical and subtropical-temperate reefs of the western Pacific

---

This README.md was generated on 2026-09-24 by Yuting Vicky Lin ([vicky.linyuting@gmail.com](mailto:vicky.linyuting@gmail.com)).

1. **Author Information**

* First & corresponding author
  * Name: Dr. Yuting Vicky Lin
  * Institution: Sesoko Marine Research Station, Tropical Biosphere Research Center, University of the Ryukyus, Okinawa 905-0227, Japan
* Second author
  * Name: Dr. Munasik Munasik
  * Institution: Department of Marine Science, Faculty of Fisheries and Marine Science, Diponegoro University, Semarang, Indonesia
* Third author
  * Name: Dr. Kakaskasen Andreas Roeroe
  * Institution: Department of Marine Science, Faculty of Fisheries and Marine Sciences, Sam Ratulangi University, Manado, Indonesia
* Fourth author
  * Name: Dr. Seiji Arakaki
  * Institution: Amakusa Marine Biological Laboratory, Kyushu University, Kumamoto, Japan
* Fifth author
  * Name: Dr. Takuma Mezaki
  * Institution: Kuroshio Biological Research Foundation, Kochi, Japan
* Sixth author
  * Name: Dr. Takashi Kawai
  * Institution: Tokyo Kyuei Co., Ltd., Saitama, Japan
* Seventh author
  * Name: Dr. Jean J.B. Tanangonan
  * Institution: Department of Environmental Management, Faculty of Agriculture, Kindai University, Nara, Japan
* Eighth author
  * Name: Dr. Diah Permata Wijayanti
  * Institution: Department of Marine Science, Faculty of Fisheries and Marine Science, Diponegoro University, Semarang, Indonesia
* Nineth & corresponding author
  * Name: Dr. Yoko Nozawa
  * Institution: Sesoko Marine Research Station, Tropical Biosphere Research Center, University of the Ryukyus, Okinawa 905-0227, Japan; Biodiversity Research Center, Academia Sinica, Taipei, Taiwan
  * Email: [nozaway@cs.u-ryukyu.ac.jp](mailto:nozaway@cs.u-ryukyu.ac.jp) 

       2.**Date of data collection**: 2012-2014

       3.**Geographic location of data collection**: Western Pacific

**&#xA0;     &#xA0;**&#x34;.**Funding sources that supported the collection of the data**:

The study was funded by the Thematic Research Grant of Academia Sinica (23-2g) and an Internal Research Grant of Biodiversity Research Center, Academia Sinica to Y.N. Y.V.L. received a fellowship from Postdoctoral Research Abroad Program (113-2917-I-564-034), funded by National Science and Technology Council, Taiwan.

**&#xA0;     &#xA0;**&#x35;.**Recommended citation for this dataset**

Lin YV, Munasik M, Roeroe K, Arakaki S, Mezaki T, Kawai T, Tanangonan J, Wijayanti DP, Nozawa Y. Taxon- and size-specific variation in coral growth between tropical and subtropical-temperate reefs of the western Pacific. Accepted by the *Journal of Biogeography* at 22 September 2026

<br />

## DESCRIPTION OF THE DATA AND FILE STRUCTURE

### DATA & FILE OVERVIEW

1. **Description of dataset**

   These data were used to analyze the growth pattern for five coral genera in the tropical and subtropical-temperate regions in the western Pacific.
2. **File List**
   * File 1: env.csv
     * File 1 description: 16 environmental parameters across 14 sites in the western Pacific
   * File 2: coral_growth.csv
     * File 2 description: Country, location, latitude, longitude, latitudinal region and other sampling information for 1258 coral colonies sampled in the western Pacific

### METHODOLOGICAL INFORMATION

A detailed description of data acquisition and processing can be found in the XXXX

#### DATA-SPECIFIC INFORMATION

##### **env.csv**

1. Number of variables/columns: 21
2. Number of cases/rows: 14
3. Missing data codes

   None
4. Variable List
   * Column A - No
   * Column B - Site
   * Column C - Lat
   * Column D - Lon
   * Column E – Lat_GP
   * Column F – SST_mean
   * Column G – SST_sd
   * Column H – SST_max
   * Column I – SST_min
   * Column J – PAR_mean
   * Column K – PAR_sd
   * Column L – PAR_max
   * Column M – PAR_min
   * Column N – chla_mean
   * Column O – chla_sd
   * Column P – chla_max
   * Column Q – chla_min
   * Column R – kd490_mean
   * Column S – kd490_sd
   * Column T – kd490_max
   * Column U – kd490_min
5. Abbreviations used
   * for variables:
     * No: Number
     * Lat: Latitude
     * Lon: Longitude
     * Lat_GP: Latitudinal region
     * SST: Sea surface temperature
     * PAR: Photosynthetically active radiation
     * chla: Chlorophyll a
     * kd490: Diffuse attenuation coefficient at 490 nm
     * sd: standard deviation

##### **coral_growth.csv**

1. Number of variables/columns: 18
2. Number of cases/rows: 1258
3. Missing data codes

   None
4. Variable List
   * Column A - Country
   * Column B - Location
   * Column C - Site
   * Column D - Latitude
   * Column E – Longitude
   * Column F – Lat_reg
   * Column G – Quadrat_number
   * Column H – CoralID
   * Column I – Genus
   * Column J – Morpho
   * Column K – Sampling_1
   * Column L – Area_1
   * Column M – Sampling_2
   * Column N – Area_2
   * Column O – Sampling_3
   * Column P – Area_3
   * Column Q – Interval_days_12
   * Column R – Interval_days_23
5. Abbreviations used
   * for environmental factors:
     * Lat_reg: Latitudinal region
     * Sampling_1: First sampling date
     * Area_1: Size at the first sampling (cm^2)
     * Sampling_2: Second sampling date
     * Area_2: Size at the second sampling (cm^2)
     * Sampling_3: Third sampling date
     * Area_3: Size at the third sampling (cm^2)
     * Interval_days_12: Number of days between first and second samplings
     * Interval_days_23: Number of days between second and third samplings

## Usage notes

R is required to open R_analysis.R and function.R; the script was created using version 4.4.0. function.R contains all  necessary functions for conducting analyses, therefore, run function.R before running the analysis in R_analysis.R. Microsoft Excel can be used to view env.csv and coral_growth.csv.
