#! /usr/bin/env pyhton3
# ONLY TO BE CALLED FROM SNAKEMAKE!!!
# Run pip3 install -r requirements.txt first
from pandas import concat, read_excel

# Load data
bird_names = read_excel(snakemake.input[0], sheetname="Bird species")
mammal_names = read_excel(snakemake.input[0], sheetname="Mammal species")
# Mammals don't have english names, so create column and populate with empty
# strings
mammal_names['eng_name'] = ''

# Change the column names to match
bird_names.columns = ['sci_name', 'eng_name']
mammal_names.columns = ['sci_name', 'eng_name']

# Fill in the English names for mammals
mammal_names.loc[mammal_names.sci_name=='Alces alces', 'eng_name'] = 'Elk'
mammal_names.loc[mammal_names.sci_name=='Ammotragus lervia', 'eng_name'] = 'Barbary sheep'
mammal_names.loc[mammal_names.sci_name=='Bison bonasus', 'eng_name'] = 'European bison'
mammal_names.loc[mammal_names.sci_name=='Canis lupus', 'eng_name'] = 'Grey wolf'
mammal_names.loc[mammal_names.sci_name=='Capra aegagrus', 'eng_name'] = 'Wild goat'
mammal_names.loc[mammal_names.sci_name=='Capra ibex', 'eng_name'] = 'Alpine ibex'
mammal_names.loc[mammal_names.sci_name=='Capra pyrenaica', 'eng_name'] = 'Iberian ibex'
mammal_names.loc[mammal_names.sci_name=='Capreolus capreolus', 'eng_name'] = 'European roe deer'
mammal_names.loc[mammal_names.sci_name=='Castor fiber', 'eng_name'] = 'Eurasian beaver'
mammal_names.loc[mammal_names.sci_name=='Cervus elaphus', 'eng_name'] = 'Red deer'
mammal_names.loc[mammal_names.sci_name=='Dama dama', 'eng_name'] = 'Fallow deer'
mammal_names.loc[mammal_names.sci_name=='Erignathus barbatus', 'eng_name'] = 'Bearded seal'
mammal_names.loc[mammal_names.sci_name=='Gulo gulo', 'eng_name'] = 'Wolverine'
mammal_names.loc[mammal_names.sci_name=='Halichoerus grypus', 'eng_name'] = 'Grey seal'
mammal_names.loc[mammal_names.sci_name=='Hystrix cristata', 'eng_name'] = 'Crested porcupine'
mammal_names.loc[mammal_names.sci_name=='Lynx lynx', 'eng_name'] = 'Eurasian lynx'
mammal_names.loc[mammal_names.sci_name=='Lynx pardinus', 'eng_name'] = 'Iberian lynx'
mammal_names.loc[mammal_names.sci_name=='Macaca sylvanus', 'eng_name'] = 'Barbary macaque'
mammal_names.loc[mammal_names.sci_name=='Meles meles', 'eng_name'] = 'European badger'
mammal_names.loc[mammal_names.sci_name=='Monachus monachus', 'eng_name'] = 'Mediterranean monk seal'
mammal_names.loc[mammal_names.sci_name=='Odobenus rosmarus', 'eng_name'] = 'Walrus'
mammal_names.loc[mammal_names.sci_name=='Phoca groenlandica', 'eng_name'] = 'Harp seal'
mammal_names.loc[mammal_names.sci_name=='Phoca hispida', 'eng_name'] = 'Ringed seal'
mammal_names.loc[mammal_names.sci_name=='Phoca vitulina', 'eng_name'] = 'Harbor seal'
mammal_names.loc[mammal_names.sci_name=='Rangifer tarandus', 'eng_name'] = 'Reindeer'
mammal_names.loc[mammal_names.sci_name=='Rupicapra pyrenaica', 'eng_name'] = 'Pyrenean chamois'
mammal_names.loc[mammal_names.sci_name=='Rupicapra rupicapra', 'eng_name'] = 'Alpine chamois'
mammal_names.loc[mammal_names.sci_name=='Sus scrofa', 'eng_name'] = 'Wild boar'
mammal_names.loc[mammal_names.sci_name=='Ursus arctos', 'eng_name'] = 'Brown bear'
mammal_names.loc[mammal_names.sci_name=='Ursus maritimus', 'eng_name'] = 'Polar bear'

# Add columns identifying the taxon
bird_names['taxon'] = 'birds'
mammal_names['taxon'] = 'mammals'
all_names = concat([bird_names, mammal_names])

# Write output
all_names.to_csv(snakemake.output[0])
